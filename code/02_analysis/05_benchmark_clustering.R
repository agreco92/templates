set.seed(0110)
options(future.globals.maxSize = 1000 * 1024^2)

library(tidyverse, quietly = T, verbose = F, warn.conflicts = F)

library(Seurat, quietly = T, verbose = F, warn.conflicts = F)
library(bluster, quietly = T, verbose = F, warn.conflicts = F)
library(scran, quietly = T, verbose = F, warn.conflicts = F)
library(igraph, quietly = T, verbose = F, warn.conflicts = F)

library(ggplot2, quietly = T, verbose = F, warn.conflicts = F)
library(ggrastr, quietly = T, verbose = F, warn.conflicts = F)
library(ggbeeswarm, quietly = T, verbose = F, warn.conflicts = F)
library(patchwork, quietly = T, verbose = F, warn.conflicts = F)
library(pheatmap, quietly = T, verbose = F, warn.conflicts = F)

library(optparse, quietly = T, verbose = F, warn.conflicts = F)
library(yaml, quietly = T, verbose = F, warn.conflicts = F)
library(glue, quietly = T, verbose = F, warn.conflicts = F)


# command line arguments -------------------------------------------------

is_interactive <- interactive()

if (is_interactive) {
  args <- list(
    dir = getwd(),
    output = NULL,
    config_name = "05_benchmark_clustering",
    module_name = "02_analysis"
  )
} else {
  option_list <- list(
    make_option(
      c("-d", "--dir"),
      type = "character",
      default = getwd(),
      help = "working directory (all paths in script are relative to this)",
      metavar = "character"
    ),
    make_option(
      c("-o", "--output"),
      type = "character",
      default = NULL,
      help = "output directory (overrides config, for Nextflow compatibility)",
      metavar = "character"
    ),
    make_option(
      c("-c", "--config_name"),
      type = "character",
      default = "05_benchmark_clustering",
      help = "name of config block to use from YAML",
      metavar = "character"
    ),
    make_option(
      c("-m", "--module_name"),
      type = "character",
      default = "02_analysis",
      help = "name of YAML config file (without .yaml extension)",
      metavar = "character"
    )
  )
  args_parser <- OptionParser(option_list = option_list)
  args <- parse_args(args_parser)
}

# helpers live in templates2; set TEMPLATES2_HOME in your .Renviron (see README)
if (!nzchar(Sys.getenv("TEMPLATES2_HOME"))) {
  stop("TEMPLATES2_HOME is not set. Add it to your .Renviron (see README).")
}
.tmpl_code <- file.path(Sys.getenv("TEMPLATES2_HOME"), "code")
source(file.path(.tmpl_code, "utils", "functions.R"))
source(file.path(.tmpl_code, "utils", "graphics_setup.R"))

# load config files ------------------------------------------------------

config <- load_config(
  project_dir = args$dir,
  module_name = args$module_name,
  script_name = args$config_name
)

cat("Configuration loaded:\n")
cat(yaml::as.yaml(config))

## resolve output directories ---------------------------------------------

if (!is.null(args$output)) {
  config$dir_output <- args$output
  config$dir_figures <- args$output
}

config$dir_output <- config$dir_output %||% "."
config$dir_figures <- config$dir_figures %||% "."

dir.create(config$dir_output, recursive = TRUE, showWarnings = FALSE)
dir.create(config$dir_figures, recursive = TRUE, showWarnings = FALSE)

message("Output directory: ", config$dir_output)
message("Figures directory: ", config$dir_figures)

# load data --------------------------------------------------------------

seu <- readRDS(config$input_path)

## resolve parameters ------------------------------------------------------

reduction <- config$reduction %||% "pca"
n_dims <- config$dims %||% 30
umap_name <- config$umap_reduction %||% "umap"
purity_k <- config$purity_k %||% 50

if (!reduction %in% names(seu@reductions)) {
  stop(
    "Reduction '",
    reduction,
    "' not found. Available: ",
    paste(names(seu@reductions), collapse = ", ")
  )
}

emb <- Embeddings(seu, reduction)
n_dims <- min(n_dims, ncol(emb))
emb <- emb[, 1:n_dims]

message(sprintf(
  "Using reduction '%s' with %d dims (%d cells)",
  reduction,
  n_dims,
  nrow(emb)
))

## resolve clustering columns ----------------------------------------------

cluster_columns <- config$cluster_columns
if (is.character(cluster_columns)) {
  cluster_columns <- as.list(cluster_columns)
}
cluster_columns <- unlist(cluster_columns)

missing_cols <- setdiff(cluster_columns, colnames(seu@meta.data))
if (length(missing_cols) > 0) {
  stop(
    "Clustering columns not found in metadata: ",
    paste(missing_cols, collapse = ", ")
  )
}

message("Evaluating clusterings: ", paste(cluster_columns, collapse = ", "))

# per-cluster quality metrics --------------------------------------------

metrics_list <- list()

for (clust_col in cluster_columns) {
  message("\n--- Metrics for: ", clust_col, " ---")
  clusters <- factor(seu@meta.data[[clust_col]])

  ## silhouette width -------------------------------------------------------

  message("  Approximate silhouette widths...")
  sil <- approxSilhouette(emb, clusters = clusters)

  ## neighbor purity --------------------------------------------------------

  message("  Neighbor purity (k=", purity_k, ")...")
  pur <- neighborPurity(emb, clusters = clusters, k = purity_k)

  ## cluster RMSD -----------------------------------------------------------

  message("  Cluster RMSD...")
  rmsd <- clusterRMSD(emb, clusters = clusters)

  # per-cell metrics
  cell_metrics <- data.frame(
    cell = colnames(seu),
    clustering = clust_col,
    cluster = as.character(clusters),
    silhouette_width = sil$width,
    nearest_other = as.character(sil$other),
    purity = pur$purity,
    max_neighbor = as.character(pur$maximum)
  )

  # per-cluster summary
  cluster_summary <- data.frame(
    clustering = clust_col,
    cluster = levels(clusters),
    n_cells = as.integer(table(clusters)),
    pct_cells = as.numeric(prop.table(table(clusters))) * 100,
    median_silhouette = as.numeric(tapply(sil$width, clusters, median)),
    mean_silhouette = as.numeric(tapply(sil$width, clusters, mean)),
    pct_negative_sil = as.numeric(tapply(sil$width < 0, clusters, mean)) * 100,
    median_purity = as.numeric(tapply(pur$purity, clusters, median)),
    rmsd = rmsd
  )

  metrics_list[[clust_col]] <- list(
    cell_metrics = cell_metrics,
    cluster_summary = cluster_summary,
    sil = sil,
    pur = pur,
    rmsd = rmsd
  )
}

# graph modularity -------------------------------------------------------

mod_list <- list()

if (isTRUE(config$modularity %||% TRUE)) {
  for (clust_col in cluster_columns) {
    clusters <- factor(seu@meta.data[[clust_col]])
    g <- NULL

    if (isTRUE(config$rebuild_graph)) {
      # build fresh SNN from reduction
      snn_k <- config$snn_k %||% 20
      snn_type <- config$snn_type %||% "jaccard"
      message("  Building SNN graph (k=", snn_k, ", type=", snn_type, ")...")
      g <- makeSNNGraph(emb, k = snn_k, type = snn_type)
    } else {
      # try to use existing graph from Seurat object
      graph_name <- config$graph_name
      if (is.null(graph_name)) {
        snn_names <- grep(
          "snn",
          names(seu@graphs),
          value = TRUE,
          ignore.case = TRUE
        )
        if (length(snn_names) > 0) graph_name <- snn_names[1]
      }

      if (!is.null(graph_name) && graph_name %in% names(seu@graphs)) {
        message("  Using stored graph: ", graph_name)
        adj <- seu@graphs[[graph_name]]
        g <- graph_from_adjacency_matrix(
          as(adj, "dgCMatrix"),
          weighted = TRUE,
          mode = "undirected"
        )
      } else {
        # fall back to building from reduction
        snn_k <- config$snn_k %||% 20
        snn_type <- config$snn_type %||% "jaccard"
        message("  No SNN graph found — building (k=", snn_k, ")...")
        g <- makeSNNGraph(emb, k = snn_k, type = snn_type)
      }
    }

    message("  Pairwise modularity for: ", clust_col)
    mod <- pairwiseModularity(g, clusters = clusters, as.ratio = TRUE)
    mod_list[[clust_col]] <- mod
  }
}

# cross-clustering comparison --------------------------------------------

ari_mat <- NULL
jacc_list <- list()

if (length(cluster_columns) > 1) {
  message("\n--- Comparing clusterings ---")

  ## adjusted rand index ----------------------------------------------------

  n_clust <- length(cluster_columns)
  ari_mat <- matrix(
    NA,
    n_clust,
    n_clust,
    dimnames = list(cluster_columns, cluster_columns)
  )

  for (i in seq_along(cluster_columns)) {
    for (j in seq_along(cluster_columns)) {
      ari_mat[i, j] <- pairwiseRand(
        seu@meta.data[[cluster_columns[i]]],
        seu@meta.data[[cluster_columns[j]]],
        mode = "index"
      )
    }
  }

  ## jaccard similarity -----------------------------------------------------

  for (i in seq_len(n_clust - 1)) {
    for (j in (i + 1):n_clust) {
      pair_name <- paste(cluster_columns[i], "vs", cluster_columns[j])
      jacc_list[[pair_name]] <- linkClustersMatrix(
        seu@meta.data[[cluster_columns[i]]],
        seu@meta.data[[cluster_columns[j]]]
      )
    }
  }
}

# parameter sweep (optional) ---------------------------------------------

sweep_results <- NULL

if (isTRUE(config$sweep$enabled)) {
  message("\n--- Parameter sweep ---")
  k_values <- config$sweep$k_values %||% c(5, 10, 15, 20, 25, 30, 35, 40)
  sweep_funs <- config$sweep$cluster_fun %||% c("louvain", "walktrap")

  if (is.character(k_values)) {
    k_values <- as.integer(k_values)
  }
  sweep_funs <- unlist(sweep_funs)

  sweep_results <- clusterSweep(
    emb,
    NNGraphParam(),
    k = as.integer(k_values),
    cluster.fun = sweep_funs
  )

  # compute metrics for each parameter combination
  sweep_metrics <- data.frame()
  for (i in seq_len(ncol(sweep_results$clusters))) {
    clust_i <- factor(sweep_results$clusters[, i])
    sil_i <- approxSilhouette(emb, clusters = clust_i)

    sweep_metrics <- rbind(
      sweep_metrics,
      data.frame(
        k = sweep_results$parameters$k[i],
        cluster_fun = sweep_results$parameters$cluster.fun[i],
        n_clusters = length(unique(clust_i)),
        mean_silhouette = mean(sil_i$width),
        median_silhouette = median(sil_i$width),
        pct_negative_sil = mean(sil_i$width < 0) * 100
      )
    )
  }
  sweep_results$metrics <- sweep_metrics
}

# summary table -----------------------------------------------------------

summary_all <- bind_rows(lapply(metrics_list, `[[`, "cluster_summary"))

write.csv(
  summary_all,
  file.path(config$dir_output, "clustering_benchmark.csv"),
  row.names = FALSE
)

cat("\n=== Clustering Benchmark Summary ===\n")
print(summary_all, row.names = FALSE)

# plots -------------------------------------------------------------------

use_subfolders <- length(cluster_columns) > 1

# adaptive point size based on dataset size
n_cells <- ncol(seu)
pt_cluster <- 3 * max(0.5, min(5, 10000 / n_cells))
pt_feature <- pt_cluster * 1.5

# helper: sort cluster labels numerically if possible, else alphabetically
sort_clusters <- function(x) {
  u <- unique(as.character(x))
  u_num <- suppressWarnings(as.numeric(u))
  if (all(!is.na(u_num))) {
    lvls <- as.character(sort(u_num))
  } else {
    lvls <- sort(u)
  }
  factor(x, levels = lvls)
}

for (clust_col in cluster_columns) {
  ml <- metrics_list[[clust_col]]
  cm <- ml$cell_metrics
  cs <- ml$cluster_summary
  safe_name <- gsub("[^a-zA-Z0-9_]", "_", clust_col)

  # consistent cluster ordering
  cm$cluster <- sort_clusters(cm$cluster)
  cs$cluster <- sort_clusters(cs$cluster)

  # subfolder per clustering when comparing multiple resolutions
  if (use_subfolders) {
    sub_figdir <- file.path(config$dir_figures, safe_name)
    sub_outdir <- file.path(config$dir_output, safe_name)
  } else {
    sub_figdir <- config$dir_figures
    sub_outdir <- config$dir_output
  }
  dir.create(sub_figdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(sub_outdir, recursive = TRUE, showWarnings = FALSE)

  # save per-clustering summary
  write.csv(
    cs,
    file.path(sub_outdir, "cluster_summary.csv"),
    row.names = FALSE
  )

  ## individual plots (build objects) ----------------------------------------

  p_sil <- ggplot(cm, aes(x = cluster, y = silhouette_width)) +
    geom_quasirandom(
      aes(color = nearest_other),
      size = 0.3, alpha = 0.6, groupOnX = TRUE
    ) +
    geom_boxplot(
      width = 0.15, outlier.shape = NA, alpha = 0.4, fill = "grey90"
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "red",
      linewidth = 0.3
    ) +
    coord_flip() +
    labs(
      x = "Cluster", y = "Silhouette width",
      title = "Silhouette width", color = "Nearest\nother"
    )

  p_pur <- ggplot(cm, aes(x = cluster, y = purity, fill = cluster)) +
    geom_boxplot(outlier.size = 0.3, alpha = 0.7) +
    coord_flip() +
    labs(x = "Cluster", y = "Neighbor purity", title = "Neighbor purity") +
    theme(legend.position = "none")

  p_rmsd <- ggplot(cs, aes(x = cluster, y = rmsd, fill = cluster)) +
    geom_col(alpha = 0.8) +
    coord_flip() +
    labs(x = "Cluster", y = "RMSD", title = "Cluster RMSD") +
    theme(legend.position = "none")

  ## modularity as ggplot (for patchwork) -----------------------------------

  p_mod <- NULL
  if (clust_col %in% names(mod_list)) {
    mod_mat <- mod_list[[clust_col]]
    mod_df <- as.data.frame(as.table(mod_mat))
    colnames(mod_df) <- c("from", "to", "ratio")
    mod_df$from <- sort_clusters(mod_df$from)
    mod_df$to <- sort_clusters(mod_df$to)
    mod_df$label <- sprintf("%.1f", mod_df$ratio)

    p_mod <- ggplot(mod_df, aes(x = from, y = to, fill = log10(ratio + 1))) +
      geom_tile() +
      geom_text(aes(label = label), size = 2.5) +
      scale_fill_gradient(
        low = "white",
        high = "#d73027",
        name = "log10(ratio+1)"
      ) +
      coord_equal() +
      labs(x = NULL, y = NULL, title = "Pairwise modularity") +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
  }

  ## combined diagnostic panel (10 x 5 in — full slide) ---------------------

  if (!is.null(p_mod)) {
    p_diagnostic <- (p_sil | p_rmsd) /
      (p_mod | p_pur) +
      plot_annotation(title = clust_col)
  } else {
    p_diagnostic <- (p_sil | p_rmsd) /
      p_pur +
      plot_annotation(title = clust_col)
  }

  ggsave(
    file.path(sub_figdir, "diagnostic_panel.pdf"),
    p_diagnostic,
    width = 10,
    height = 7.5
  )

  ## UMAP by cluster (5 x 4 in — half slide) --------------------------------

  if (umap_name %in% names(seu@reductions)) {
    p_umap_cluster <- DimPlot(
      seu,
      group.by = clust_col,
      reduction = umap_name,
      label = TRUE,
      pt.size = pt_cluster,
      raster = TRUE,
      raster.dpi = c(800, 800)
    ) +
      ggtitle(clust_col) +
      theme(legend.position = "right")

    ggsave(
      file.path(sub_figdir, "umap_clusters.pdf"),
      p_umap_cluster,
      width = 5,
      height = 4
    )

    ## UMAP diagnostics (10 x 3.5 in — full width) --------------------------

    seu$.bench_sil <- ml$sil$width
    seu$.bench_pur <- ml$pur$purity

    p_umap_sil <- FeaturePlot(
      seu,
      features = ".bench_sil",
      reduction = umap_name,
      pt.size = pt_feature,
      raster = TRUE,
      raster.dpi = c(800, 800)
    ) +
      scale_color_gradient2(
        low = "#d73027",
        mid = "#ffffbf",
        high = "#1a9850",
        midpoint = 0,
        name = "Silhouette"
      ) +
      ggtitle("Silhouette width")

    p_umap_pur <- FeaturePlot(
      seu,
      features = ".bench_pur",
      reduction = umap_name,
      pt.size = pt_feature,
      raster = TRUE,
      raster.dpi = c(800, 800)
    ) +
      scale_color_viridis_c(option = "viridis", name = "Purity") +
      ggtitle("Neighbor purity")

    p_umap_diag <- p_umap_cluster | p_umap_sil | p_umap_pur

    ggsave(
      file.path(sub_figdir, "umap_diagnostics.pdf"),
      p_umap_diag,
      width = 10,
      height = 3.5
    )

    seu$.bench_sil <- NULL
    seu$.bench_pur <- NULL
  }

  ## modularity heatmap (standalone — half slide) ---------------------------

  if (clust_col %in% names(mod_list)) {
    mod_mat <- mod_list[[clust_col]]

    pdf(
      file.path(sub_figdir, "modularity_heatmap.pdf"),
      width = 5,
      height = 4.5
    )
    pheatmap(
      log10(mod_mat + 1),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      display_numbers = TRUE,
      number_format = "%.2f",
      color = colorRampPalette(c("white", "#fee08b", "#d73027"))(100),
      main = paste0("Pairwise modularity (log10) \u2014 ", clust_col)
    )
    dev.off()
  }
}

## comparison heatmaps (top-level, across clusterings) ---------------------

if (!is.null(ari_mat)) {
  pdf(
    file.path(config$dir_figures, "comparison_ari.pdf"),
    width = 5,
    height = 4.5
  )
  pheatmap(
    ari_mat,
    display_numbers = TRUE,
    number_format = "%.3f",
    color = colorRampPalette(c("white", "#2166ac"))(100),
    main = "Adjusted Rand Index"
  )
  dev.off()
}

if (length(jacc_list) > 0) {
  for (pair_name in names(jacc_list)) {
    safe_pair <- gsub("[^a-zA-Z0-9_]", "_", pair_name)
    pdf(
      file.path(
        config$dir_figures,
        paste0("comparison_jaccard_", safe_pair, ".pdf")
      ),
      width = 5,
      height = 4.5
    )
    pheatmap(
      jacc_list[[pair_name]],
      display_numbers = TRUE,
      number_format = "%.2f",
      color = colorRampPalette(c("white", "#1a9850"))(100),
      main = paste0("Jaccard similarity \u2014 ", pair_name)
    )
    dev.off()
  }
}

## sweep plots -------------------------------------------------------------

if (!is.null(sweep_results)) {
  sm <- sweep_results$metrics

  p_ncl <- ggplot(sm, aes(x = k, y = n_clusters, color = cluster_fun)) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1.5) +
    labs(x = "k (neighbors)", y = "Number of clusters", color = "Algorithm") +
    ggtitle("Parameter sweep \u2014 cluster count")

  p_sil_sweep <- ggplot(
    sm,
    aes(x = k, y = mean_silhouette, color = cluster_fun)
  ) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1.5) +
    labs(
      x = "k (neighbors)",
      y = "Mean silhouette width",
      color = "Algorithm"
    ) +
    ggtitle("Parameter sweep \u2014 silhouette")

  p_neg <- ggplot(sm, aes(x = k, y = pct_negative_sil, color = cluster_fun)) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1.5) +
    labs(
      x = "k (neighbors)",
      y = "% negative silhouette",
      color = "Algorithm"
    ) +
    ggtitle("Parameter sweep \u2014 misassigned cells")

  p_sweep <- p_ncl / p_sil_sweep / p_neg + plot_layout(guides = "collect")

  ggsave(
    file.path(config$dir_figures, "sweep_metrics.pdf"),
    p_sweep,
    width = 5,
    height = 7
  )

  write.csv(
    sm,
    file.path(config$dir_output, "sweep_metrics.csv"),
    row.names = FALSE
  )
}

# save outputs ------------------------------------------------------------

yaml::write_yaml(config, file.path(config$dir_output, "config.yaml"))
save.image(file.path(config$dir_output, "workspace.rda"))

message("\n=== Clustering benchmark complete! ===")
message("Outputs: ", config$dir_output)
message("Figures: ", config$dir_figures)
