set.seed(0110)
options(future.globals.maxSize = 1000 * 1024^2)
library(reticulate, quietly = T, verbose = F, warn.conflicts = F)
# python interpreter comes from RETICULATE_PYTHON (see .Renviron.example)
suppressPackageStartupMessages({
  library(tidyverse, quietly = T, verbose = F, warn.conflicts = F)

  library(Seurat, quietly = T, verbose = F, warn.conflicts = F)
  library(DropletUtils, quietly = T, verbose = F, warn.conflicts = F)
  library(scDblFinder, quietly = T, verbose = F, warn.conflicts = F)
  library(scuttle, quietly = T, verbose = F, warn.conflicts = F)
  library(scran, quietly = T, verbose = F, warn.conflicts = F)
  library(ggplot2, quietly = T, verbose = F, warn.conflicts = F)
  library(ggthemes, quietly = T, verbose = F, warn.conflicts = F)
  library(patchwork, quietly = T, verbose = F, warn.conflicts = F)

  library(optparse, quietly = T, verbose = F, warn.conflicts = F)
  library(jsonlite, quietly = T, verbose = F, warn.conflicts = F)
  library(yaml, quietly = T, verbose = F, warn.conflicts = F)
  library(glue, quietly = T, verbose = F, warn.conflicts = F)
})

# command line arguments -------------------------------------------------

is_interactive <- interactive()

if (is_interactive) {
  args <- list(
    dir = getwd(),
    output = NULL,
    sample_id = "IMQ_Rep1",
    config_name = "01_qc_normalize_dr",
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
      c("-s", "--sample_id"),
      type = "character",
      default = NULL,
      help = "sample name",
      metavar = "character"
    ),
    make_option(
      c("-c", "--config_name"),
      type = "character",
      default = "01_qc_normalize_dr",
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
  script_name = args$config_name,
  sample_id = args$sample_id
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

# resolve input file -----------------------------------------------------

sample_path <- file.path(config$sample_dir, config$sample_id)

if (is.null(config$input_filename) || config$input_filename == "auto") {
  rds_files <- list.files(sample_path, pattern = "\\.rds$", full.names = FALSE)
  if (length(rds_files) == 0) {
    stop("No .rds files found in ", sample_path)
  } else if (length(rds_files) > 1) {
    stop(
      "Multiple .rds files found in ",
      sample_path,
      ": ",
      paste(rds_files, collapse = ", "),
      ". Set input_filename in config."
    )
  }
  input_file <- file.path(sample_path, rds_files[1])
  message("Auto-detected input file: ", input_file)
} else {
  input_file <- file.path(sample_path, config$input_filename)
}

# QC --------------------------------------------------------------------

seu <- readRDS(input_file)

# drop HTO assay
seu@active.assay <- "RNA"
seu@assays[["HTO"]] <- NULL

#
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
seu[["percent.hb"]] <- PercentageFeatureSet(seu, pattern = "^Hb[^(p)]")
seu[["nCount_RNA_log10"]] <- log10(seu$nCount_RNA)

# add metadata from config -----------------------------------------------
if (!is.null(config$metadata)) {
  message("Adding metadata from config...")
  for (meta_key in names(config$metadata)) {
    meta_value <- config$metadata[[meta_key]]
    message(sprintf("  - Adding %s: %s", meta_key, meta_value))
    seu[[meta_key]] <- meta_value
  }
}

# add metadata from csv --------------------------------------------------
# config$metadata_csv:
#   path: path to CSV file (relative to project dir or absolute)
#   sample_column: column in CSV used to match sample_id (default: "sample_id")
#   columns: list of columns to add (null = all columns except sample_column)
if (!is.null(config$metadata_csv)) {
  csv_path <- config$metadata_csv$path
  if (!startsWith(csv_path, "/")) {
    csv_path <- file.path(args$dir, csv_path)
  }
  message("Loading metadata from CSV: ", csv_path)
  meta_csv <- read.csv(csv_path, stringsAsFactors = FALSE)

  sample_col <- config$metadata_csv$sample_column %||% "sample_id"
  if (!sample_col %in% colnames(meta_csv)) {
    stop("Column '", sample_col, "' not found in metadata CSV.")
  }

  row_idx <- which(meta_csv[[sample_col]] == config$sample_id)
  if (length(row_idx) == 0) {
    stop("sample_id '", config$sample_id, "' not found in column '", sample_col, "' of metadata CSV.")
  }
  if (length(row_idx) > 1) {
    warning("Multiple rows match sample_id '", config$sample_id, "' — using first.")
    row_idx <- row_idx[1]
  }

  cols_to_add <- config$metadata_csv$columns %||% setdiff(colnames(meta_csv), sample_col)
  for (col in cols_to_add) {
    if (!col %in% colnames(meta_csv)) {
      warning("Column '", col, "' not found in CSV, skipping.")
      next
    }
    val <- as.character(meta_csv[row_idx, col])
    message(sprintf("  - Adding %s: %s (from CSV)", col, val))
    seu[[col]] <- val
  }
}


qc_plot_x <- ifelse(!is.null(config$split_by), config$split_by, "sample_id")
p1 <- seu@meta.data |>
  ggplot(aes(x = !!sym(qc_plot_x), y = percent.mt)) +
  ggbeeswarm::geom_quasirandom(
    aes(col = !!sym(qc_plot_x)),
    stroke = 0
  ) +
  geom_hline(yintercept = config[["qc_thresholds"]][["max_mito_pct"]]) +
  scale_y_sqrt()


p2 <- seu@meta.data |>
  ggplot(aes(x = !!sym(qc_plot_x), y = percent.hb)) +
  ggbeeswarm::geom_quasirandom(
    aes(col = !!sym(qc_plot_x)),
    stroke = 0
  ) +
  geom_hline(yintercept = config[["qc_thresholds"]][["max_hb_pct"]]) +
  scale_y_sqrt()

p3 <- seu@meta.data |>
  ggplot(aes(x = !!sym(qc_plot_x), y = nCount_RNA_log10)) +
  ggbeeswarm::geom_quasirandom(
    aes(col = !!sym(qc_plot_x)),
    stroke = 0
  ) +
  geom_hline(yintercept = config[["qc_thresholds"]][["min_counts_log10"]]) +
  scale_y_sqrt()

p_qc <- p1 + p2 + p3 + plot_layout(ncol = 3, guides = "collect")
ggsave(
  paste0(config$dir_figures, "/", config$sample_id, "_qc.pdf"),
  p_qc,
  width = 10,
  height = 4
)
exclude <- (seu$nCount_RNA_log10 <
  config[["qc_thresholds"]][["min_counts_log10"]]) +
  (seu$percent.mt > config[["qc_thresholds"]][["max_mito_pct"]]) +
  (seu$percent.hb > config[["qc_thresholds"]][["max_hb_pct"]])
table(exclude)

seu <- seu[, exclude == 0]

# further exclusion of cells based on metadata listed in exclude list
if (!is.null(config$exclude)) {
  for (i in 1:length(config$exclude)) {
    seu <- seu[, seu[[names(config$exclude)[i]]] != config$exclude[[i]]]
  }
}

# doublet removal --------------------------------------------------------
if (config$remove_doublets) {
  tmp <- scDblFinder(seu@assays$RNA$counts)
  seu@meta.data[, colnames(colData(tmp))] <- colData(tmp) |> as.data.frame()
}

# normalization ------------------------------------------------------

if (!is.null(config$ambient_corrected_counts)) {
  noamb_path <- paste0(
    config$sample_dir,
    "/",
    config$sample_id,
    "/",
    config$ambient_corrected_counts
  )
  noamb <- readRDS(noamb_path)
  seu@assays[["RNA"]] <- CreateAssayObject(counts = noamb[, Cells(seu)])
}

seu <- pool_normalize(
  seurat_obj = seu,
  do_clustering = TRUE,
  cluster_col = NULL,
  batch_col = NULL
)

# clustering and dr ------------------------------------------------------
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:20, verbose = FALSE, densmap = T)

seu <- FindNeighbors(seu, dims = 1:20, verbose = FALSE)
seu <- FindClusters(
  seu,
  resolution = config$cluster_res,
  verbose = FALSE
)
p4 <- DimPlot(seu, label = TRUE)

p4

ggsave(
  paste0(config$dir_figures, "/", config$sample_id, "_umap_clusters.pdf"),
  p4,
  width = 6,
  height = 5
)

markers <- FindAllMarkers(seu, only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

p5 <- DoHeatmap(seu, features = top5$gene) + NoLegend()

ggsave(
  paste0(config$dir_figures, "/", config$sample_id, "_marker_heatmap.pdf"),
  p5,
  width = 7,
  height = 8
)

p_cluster_qc <- seu@meta.data |>
  ggplot(aes(x = seurat_clusters, y = nCount_RNA_log10)) +
  ggbeeswarm::geom_quasirandom(aes(col = seurat_clusters)) +
  scale_color_tableau("Tableau 20")

ggsave(
  paste0(config$dir_figures, "/", config$sample_id, "_cluster_qc.pdf"),
  p_cluster_qc,
  width = 6,
  height = 4
)
if (config$remove_doublets) {
  p_cluster_doublet <- seu@meta.data |>
    ggplot(aes(x = seurat_clusters, y = scDblFinder.score)) +
    ggbeeswarm::geom_quasirandom(aes(col = scDblFinder.class)) +
    scale_color_calc()

  ggsave(
    paste0(config$dir_figures, "/", config$sample_id, "_cluster_doublets.pdf"),
    p_cluster_doublet,
    width = 6,
    height = 4
  )
}

# cluster vs split_by heatmap --------------------------------------------

if (!is.null(config$split_by) && config$split_by %in% colnames(seu@meta.data)) {
  # Create contingency table
  cluster_split_table <- table(
    seu$seurat_clusters,
    seu@meta.data[[config$split_by]]
  )

  # Convert to proportions (by cluster)
  cluster_split_prop <- prop.table(cluster_split_table, margin = 1)

  # Convert to long format for ggplot
  cluster_split_df <- as.data.frame(cluster_split_prop)
  colnames(cluster_split_df) <- c("cluster", "split_var", "proportion")

  p_cluster_split <- ggplot(
    cluster_split_df,
    aes(x = split_var, y = cluster, fill = proportion)
  ) +
    geom_tile() +
    geom_text(
      aes(label = sprintf("%.2f", proportion)),
      color = "white",
      size = 3
    ) +
    scale_fill_viridis_c(option = "magma") +
    labs(
      x = config$split_by,
      y = "Cluster",
      fill = "Proportion",
      title = paste0("Cluster composition by ", config$split_by)
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )

  ggsave(
    paste0(
      config$dir_figures,
      "/",
      config$sample_id,
      "_cluster_split_heatmap.pdf"
    ),
    p_cluster_split,
    width = 8,
    height = 6
  )
}


# remove doublets --------------------------------------------------------
if (config$remove_doublets) {
  seu <- seu[, seu$scDblFinder.class != "doublet"]
}

# save object ------------------------------------------------------------
saveRDS(
  seu,
  paste0(config$dir_output, "/qc_normalized.rds")
)
