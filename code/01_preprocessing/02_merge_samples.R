set.seed(0110)
options(future.globals.maxSize = 1000 * 1024^2)
library(reticulate, quietly = T, verbose = F, warn.conflicts = F)
# python interpreter comes from RETICULATE_PYTHON (see .Renviron.example)

library(tidyverse, quietly = T, verbose = F, warn.conflicts = F)

library(Seurat, quietly = T, verbose = F, warn.conflicts = F)
library(ggplot2, quietly = T, verbose = F, warn.conflicts = F)
library(ggthemes, quietly = T, verbose = F, warn.conflicts = F)
library(patchwork, quietly = T, verbose = F, warn.conflicts = F)

library(optparse, quietly = T, verbose = F, warn.conflicts = F)
library(jsonlite, quietly = T, verbose = F, warn.conflicts = F)
library(yaml, quietly = T, verbose = F, warn.conflicts = F)
library(glue, quietly = T, verbose = F, warn.conflicts = F)

# command line arguments -------------------------------------------------

is_interactive <- interactive()

if (is_interactive) {
  args <- list(
    dir = getwd(),
    output = NULL,
    config_name = "02_merge_samples",
    module_name = "01_preprocessing"
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
      default = "02_merge_samples",
      help = "name of config block to use from YAML",
      metavar = "character"
    ),
    make_option(
      c("-m", "--module_name"),
      type = "character",
      default = "01_preprocessing",
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

# load config files  ----------------------------------------------------

config <- load_config(
  project_dir = args$dir,
  module_name = args$module_name,
  script_name = args$config_name
)

cat("Configuration loaded:\n")
cat(yaml::as.yaml(config))

## resolve output directories ---------------------------------------------

# CLI --output takes priority (for Nextflow compatibility)
if (!is.null(args$output)) {
  config$dir_output <- args$output
  config$dir_figures <- args$output
}

# Default to current directory if not set anywhere
config$dir_output <- config$dir_output %||% "."
config$dir_figures <- config$dir_figures %||% "."

dir.create(config$dir_output, recursive = TRUE, showWarnings = FALSE)
dir.create(config$dir_figures, recursive = TRUE, showWarnings = FALSE)

message("Output directory: ", config$dir_output)
message("Figures directory: ", config$dir_figures)

# load samples -----------------------------------------------------------

samples <- unlist(lapply(config$sample_dir, function(d) {
  dir(d, recursive = TRUE, pattern = ".rds$", full.names = TRUE)
}))

# Exclude samples by regex pattern if specified
if (!is.null(config$sample_exclude_regex)) {
  n_before <- length(samples)
  # need sample names to exclude based on sample regex
  sample_names <- unlist(lapply(config$sample_dir, function(d) {
    dir(d)
  }))

  samples <- samples[!grepl(config$sample_exclude_regex, sample_names)]
  cat(sprintf(
    "Sample exclusion: %d -> %d samples (excluded %d matching '%s')\n",
    n_before,
    length(samples),
    n_before - length(samples),
    config$sample_exclude_regex
  ))
}

seu_list <- lapply(samples, function(s) {
  seu <- readRDS(s)
  seu@active.assay <- "RNA"
  # seu@assays[["SCT"]] <- NULL
  return(seu)
})

seu <- Reduce(merge, seu_list)


## join layers (if needed) -----------------------------------------------

# Join layers only if there are multiple layers per sample
# (i.e., if layers haven't been automatically joined during merge)
if (length(Layers(seu[["RNA"]])) > 2) {
  message("Joining ", length(Layers(seu[["RNA"]])), " layers...")
  seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
} else {
  message(
    "Layers already joined (",
    length(Layers(seu[["RNA"]])),
    " layers present)"
  )
}

# Convert to V3 assay if needed
if (inherits(seu[["RNA"]], "Assay5")) {
  message("Converting RNA assay from V5 to V3...")
  seu <- scCustomize::Convert_Assay(seu, assay = "RNA", convert_to = "V3")
} else {
  message("RNA assay is already V3 format")
}

# remove contaminants ----------------------------------------------------

cat("Computing contamination scores and filtering...\n")
if ("contamination_markers" %in% names(config)) {
  # quick normalization only for this task (will be repeated later)
  seu <- NormalizeData(seu)
  # Loop over contamination marker sets
  contamination_types <- names(config[["contamination_markers"]])
  seu$exclude_contamination <- FALSE

  for (contamination_type in contamination_types) {
    # Get markers and threshold for this contamination type
    markers <- config[["contamination_markers"]][[contamination_type]]
    threshold <- config[["contamination_thresholds"]][[contamination_type]]

    cat(sprintf("  Processing %s contamination...\n", contamination_type))

    # Score the module
    seu <- AddModuleScore(
      seu,
      features = list(markers)
    )

    # Rename from Cluster1 to <contamination_type>_score
    score_col <- paste0(contamination_type, "_score")
    seu[[score_col]] <- seu$Cluster1
    seu$Cluster1 <- NULL

    # Flag cells exceeding threshold
    exclude_col <- paste0("exclude_", contamination_type)
    seu[[exclude_col]] <- seu[[score_col]] > threshold

    # Update overall exclusion flag
    seu$exclude_contamination <- seu$exclude_contamination | seu[[exclude_col]]

    # Report stats
    n_excluded <- sum(seu[[exclude_col]])
    cat(sprintf(
      "    Flagged %d cells (%.2f%%) with score > %.2f\n",
      n_excluded,
      100 * n_excluded / ncol(seu),
      threshold
    ))

    # Diagnostic density plot of contamination score with threshold line
    df_score <- data.frame(
      score = as.numeric(seu@meta.data[, score_col, drop = T])
    )
    p_score <- ggplot(df_score, aes(x = score)) +
      geom_density(fill = "steelblue", alpha = 0.4) +
      coord_transform(y = "log1p") +
      geom_vline(
        xintercept = threshold,
        color = "red",
        linetype = "dashed",
        size = 1
      ) +
      labs(
        title = paste0("Contamination score: ", contamination_type),
        x = "Score",
        y = "Density",
        subtitle = paste0("Threshold = ", threshold)
      ) +
      theme_minimal()

    ggsave(
      p_score,
      filename = paste0(
        config$dir_figures,
        "/contamination_",
        contamination_type,
        "_score_density.pdf"
      ),
      width = 6,
      height = 4
    )
  }

  # Report total exclusions and filter
  n_total_excluded <- sum(seu$exclude_contamination)
  cat(sprintf(
    "  Total cells to exclude: %d (%.2f%%)\n",
    n_total_excluded,
    100 * n_total_excluded / ncol(seu)
  ))

  seu <- seu[, !seu$exclude_contamination]
  cat(sprintf("  Cells remaining after filtering: %d\n\n", ncol(seu)))
}


# preprocessing (normalization to PCA) -----------------------------------
seu <- pool_normalize(seu, do_clustering = T, batch_col = config$batch_column)
seu <- FindVariableFeatures(seu, verbose = FALSE)

# Optional: drop named genes from the variable features before PCA/integration,
# so they cannot shape the embedding, while keeping them in the data layer for
# markers/DGE. `variable_features_exclude` is either an inline list of gene names
# or a path to a file with one gene per line. NULL (default) = no change.
if (!is.null(config$variable_features_exclude)) {
  vf_exclude <- config$variable_features_exclude
  if (length(vf_exclude) == 1 && is.character(vf_exclude) && file.exists(vf_exclude)) {
    vf_exclude <- readLines(vf_exclude)
  }
  vf_exclude <- unique(unlist(vf_exclude))
  vf <- VariableFeatures(seu)
  hit <- intersect(vf, vf_exclude)
  VariableFeatures(seu) <- setdiff(vf, vf_exclude)
  cat(sprintf(
    "Variable-feature exclusion: %d -> %d HVGs (removed %d of %d requested genes present in the HVG set)\n",
    length(vf), length(VariableFeatures(seu)), length(hit), length(vf_exclude)
  ))
}

seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 30, verbose = FALSE)

# Default linear DR is PCA
dr_linear <- "pca" # will change to rpca or cca if integration is performed

# integration (optional) -------------------------------------------------

if (
  "integration" %in%
    names(config) &&
    isTRUE(config[["integration"]][["enabled"]])
) {
  cat("Splitting layers for integration...\n")
  seu[["RNA"]] <- split(
    seu[["RNA"]],
    f = seu@meta.data[, config$batch_column, drop = TRUE]
  )
  cat("  Number of layers:", length(Layers(seu[["RNA"]])), "\n")
  cat("Starting integration...\n")

  if (config[["integration"]][["method"]] == "rpca") {
    cat("  Performing RPCA integration...\n")
    seu <- IntegrateLayers(
      object = seu,
      method = RPCAIntegration,
      orig.reduction = "pca",
      new.reduction = "integrated.rpca",
      dims = 1:30,
      verbose = TRUE
    )
    dr_linear <- "integrated.rpca"
  } else if (config[["integration"]][["method"]] == "cca") {
    cat("  Performing CCA integration...\n")
    seu <- IntegrateLayers(
      object = seu,
      method = CCAIntegration,
      orig.reduction = "pca",
      new.reduction = "integrated.cca",
      dims = 1:30,
      verbose = TRUE
    )
    dr_linear <- "integrated.cca"
  } else {
    warning(
      "Unknown integration method: ",
      config[["integration"]][["method"]],
      ". Using PCA."
    )
  }

  cat("  Integration complete. Using reduction:", dr_linear, "\n")

  # Re-join layers after integration
  seu <- JoinLayers(seu)

  # Convert to V3 assay if needed (after integration)
  if (inherits(seu[["RNA"]], "Assay5")) {
    message("Converting RNA assay from V5 to V3 after integration...")
    seu <- scCustomize::Convert_Assay(seu, assay = "RNA", convert_to = "V3")
  } else {
    message("RNA assay is already V3 format")
  }
}

# neighbors, UMAP, clustering -------------------------------------------

cat("Computing neighbors and UMAP using:", dr_linear, "\n")
seu <- FindNeighbors(seu, reduction = dr_linear, dims = 1:30)
seu <- RunUMAP(seu, reduction = dr_linear, dims = 1:30)

seu <- FindClusters(
  seu,
  resolution = as.numeric(config$cluster_res),
  verbose = FALSE,
  algorithm = 4,
  cluster.name = paste0("merged_leiden_res_", config$cluster_res)
)

## markers ----------------------------------------------------------------

markers <- FindAllMarkers(
  seu,
  only.pos = TRUE,
  group.by = paste0("merged_leiden_res_", config$cluster_res)
)
markers |>
  group_by(cluster) |>
  dplyr::filter(avg_log2FC > 1)

markers |>
  group_by(cluster) |>
  dplyr::filter(avg_log2FC > 1, p_val_adj < 0.01) |>
  slice_head(n = 20) |>
  ungroup() -> top20
write_csv(
  top20 |> select(gene, cluster, p_val_adj, avg_log2FC),
  file = paste0(config$dir_output, "/cluster_markers.csv")
)

markers |>
  group_by(cluster) |>
  dplyr::filter(avg_log2FC > 1) |>
  slice_head(n = 5) |>
  ungroup() -> top5


# assign condition  ------------------------------------------------------------

if ("condition_regex" %in% names(config)) {
  # Initialize with default condition

  seu$condition <- config[["condition_regex"]][["default"]]

  # Get the column to match against

  match_col <- config[["condition_regex"]][["column"]]

  # Iterate over condition_regex$dict: keys are condition names, values are regex patterns

  for (condition_name in names(config[["condition_regex"]][["dict"]])) {
    pattern <- config[["condition_regex"]][["dict"]][[condition_name]]
    matches <- grepl(pattern, seu@meta.data[, match_col, drop = TRUE])
    seu@meta.data[matches, "condition"] <- condition_name
  }
}

# plots -----------------------------------------------------------------

p1 <- DimPlot(seu, group.by = config$sample_column, split.by = "condition") +
  scale_color_tableau()
p1

ggsave(
  p1,
  filename = paste0(config$dir_figures, "/umap_sample.pdf"),
  width = 10,
  height = 9
)

p2 <- DimPlot(
  seu,
  group.by = "seurat_clusters",
  split.by = config$sample_column,
  ncol = 2
) +
  scale_color_tableau("Tableau 20")
p2
ggsave(
  p2,
  filename = paste0(config$dir_figures, "/umap_sample_clusters.pdf"),
  width = 15,
  height = 18
)

p3 <- DimPlot(
  seu,
  group.by = paste0("merged_leiden_res_", config$cluster_res),
  split.by = "condition",
  label = T,
  ncol = 2
) +
  scale_color_tableau("Tableau 20")

ggsave(
  p3,
  filename = paste0(config$dir_figures, "/umap_condition_clusters.pdf"),
  width = 13,
  height = 9
)

p4 <- DoHeatmap(seu, features = top5$gene) + NoLegend()
ggsave(
  p4,
  filename = paste0(config$dir_figures, "/merged_marker_heatmap.pdf"),
  width = 8,
  height = 9
)

p5 <- DotPlot(
  seu,
  features = unique(top5$gene),
  group.by = paste0("merged_leiden_res_", config$cluster_res)
) +
  coord_flip()
ggsave(
  p5,
  filename = paste0(config$dir_figures, "/merged_marker_dotplot.pdf"),
  width = 8,
  height = 13
)

# save object ------------------------------------------------------------

saveRDS(object = seu, file = paste0(config$dir_output, "/merged.rds"))

writeLines(
  yaml::as.yaml(config),
  con = file.path(config$dir_output, "config_used.yaml")
)
