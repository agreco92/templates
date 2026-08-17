set.seed(0110)
library(reticulate, quietly = T, verbose = F)
# python interpreter comes from RETICULATE_PYTHON (see .Renviron.example)
suppressPackageStartupMessages({
  library(tidyverse, quietly = T, verbose = F)

  library(Seurat, quietly = T, verbose = F)
  library(DropletUtils, quietly = T, verbose = F)
  library(SoupX, quietly = T, verbose = F)
  library(ggplot2, quietly = T, verbose = F)
  library(ggthemes, quietly = T, verbose = F)
  library(patchwork, quietly = T, verbose = F)

  library(optparse, quietly = T, verbose = F)
  library(jsonlite, quietly = T, verbose = F)
  library(yaml, quietly = T, verbose = F)
  library(glue, quietly = T, verbose = F)
})

# command line arguments -------------------------------------------------

is_interactive <- interactive()

if (is_interactive) {
  args <- list(
    dir = getwd(),
    output = NULL,
    sample_id = "I1",
    config_name = "00_cellcalling_demux",
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
      help = "sample name (REQUIRED)",
      metavar = "character"
    ),
    make_option(
      c("-c", "--config_name"),
      type = "character",
      default = "00_cellcalling_demux",
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
  renv::load(args$dir)

  # Check that sample_id is provided
  if (is.null(args$sample_id)) {
    stop("Error: sample_id (-s/--sample_id) is required!")
  }
}

# Load utility functions
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

# Build input path for this sample ----------------------------------------

input_folder <- file.path(
  config$sample_dir,
  config$sample_id,
  config$sample_subpath
)

cat("Processing sample:", config$sample_id, "\n")
cat("Input folder (raw matrix):", input_folder, "\n")
cat("Output folder:", config$dir_output, "\n")
cat("Figures folder:", config$dir_figures, "\n\n")

# Check if input folder exists
if (!dir.exists(input_folder)) {
  stop("Error: Input folder does not exist: ", input_folder)
}

# read data --------------------------------------------------------------

cat("Reading raw data...\n")

# Determine input type (default to 10x for backwards compatibility)
input_type <- config[["input"]][["type"]] %||% "10x"
cat("Input type:", input_type, "\n")

if (input_type == "10x") {
  # Standard 10x Genomics format using Read10X
  cat("Using Read10X for 10x Genomics format...\n")
  data_10x_raw <- Read10X(data.dir = input_folder)

  # Check if data contains multiple modalities
  if (is.list(data_10x_raw)) {
    # Support both Antibody Capture (HTO) and Multiplexing Capture (CMO)
    if ("Antibody Capture" %in% names(data_10x_raw)) {
      hto_counts_raw <- data_10x_raw[["Antibody Capture"]]
    } else if ("Multiplexing Capture" %in% names(data_10x_raw)) {
      hto_counts_raw <- data_10x_raw[["Multiplexing Capture"]]
    } else {
      hto_counts_raw <- NULL
    }
    gex_counts_raw <- data_10x_raw[["Gene Expression"]]
  } else {
    gex_counts_raw <- data_10x_raw
    hto_counts_raw <- NULL
  }
} else {
  # Generic matrix format using ReadMtx
  cat("Using ReadMtx for generic matrix format...\n")

  # Get file names from config, with sensible defaults

  mtx_file <- config[["input"]][["mtx"]] %||% "matrix.mtx.gz"
  barcodes_file <- config[["input"]][["barcodes"]] %||% "barcodes.tsv.gz"
  features_file <- config[["input"]][["features"]] %||% "features.tsv.gz"
  feature_col <- config[["input"]][["feature_col"]] %||% 2

  cat("  Matrix file:", mtx_file, "\n")
  cat("  Barcodes file:", barcodes_file, "\n")
  cat("  Features file:", features_file, "\n")
  cat("  Feature column:", feature_col, "\n")

  gex_counts_raw <- ReadMtx(
    mtx = file.path(input_folder, mtx_file),
    cells = file.path(input_folder, barcodes_file),
    features = file.path(input_folder, features_file),
    feature.column = feature_col
  )

  # Generic format typically doesn't have HTO
  hto_counts_raw <- NULL
}

cat("Raw Gene Expression features:", nrow(gex_counts_raw), "\n")
cat("Raw Gene Expression barcodes:", ncol(gex_counts_raw), "\n")

# cell calling  -------------------------------------------------------------

cat("Running emptyDrops for cell calling...\n")
e.out <- emptyDrops(gex_counts_raw, test.ambient = TRUE)
fdr_threshold <- config$dropletutils_fdr

# Identify cells passing FDR threshold
is_cell <- which(e.out$FDR <= fdr_threshold)

gex_low_umi_q01 <- log10(quantile(
  e.out$Total[e.out$FDR < fdr_threshold],
  probs = .01,
  na.rm = TRUE
))
bc_rank <- barcodeRanks(gex_counts_raw)

cat(
  "Cells called (FDR <",
  fdr_threshold,
  "):",
  length(is_cell),
  "\n\n"
)

# Filter to called cells
gex_counts <- gex_counts_raw[, is_cell]
if (!is.null(hto_counts_raw)) {
  hto_counts <- hto_counts_raw[, is_cell]
}

cat("Filtered Gene Expression barcodes:", ncol(gex_counts), "\n")

# SoupX ambient RNA correction --------------------------------------------

if (isTRUE(config$ambient_correction)) {
  cat("\n==================================================\n")
  cat("Running SoupX for ambient RNA correction...\n")
  cat("==================================================\n\n")
  soupX_cluster_res <- config$soupx_cluster_res %||% 0.1
  # Create SoupChannel object
  # tod = table of droplets (raw data)
  # toc = table of counts (filtered cells)
  cat("Creating SoupChannel object...\n")
  sc <- SoupChannel(tod = gex_counts_raw, toc = gex_counts)

  cat("SoupChannel created with:\n")
  cat("  Genes:", nrow(sc$toc), "\n")
  cat("  Cells:", ncol(sc$toc), "\n")

  # Perform basic clustering for SoupX
  # SoupX works better with clustering information
  cat("Performing quick clustering for SoupX...\n")
  seu_temp <- CreateSeuratObject(
    counts = gex_counts,
    project = config$sample_id
  )
  seu_temp <- NormalizeData(seu_temp, verbose = FALSE)
  seu_temp <- FindVariableFeatures(seu_temp, verbose = FALSE)
  seu_temp <- ScaleData(seu_temp, verbose = FALSE)
  seu_temp <- RunPCA(seu_temp, verbose = FALSE)
  seu_temp <- FindNeighbors(seu_temp, dims = 1:20, verbose = FALSE)
  seu_temp <- FindClusters(seu_temp, resolution = 0.1, verbose = FALSE)

  # Add clustering information to SoupChannel
  clusters <- setNames(seu_temp$seurat_clusters, colnames(seu_temp))
  sc <- setClusters(sc, clusters)

  cat("Clusters assigned:", length(unique(clusters)), "clusters\n\n")

  # Estimate contamination fraction automatically
  cat("Estimating contamination fraction with autoEstCont...\n")
  pdf(
    file.path(config$dir_figures, "soupX_autoEstCont_diagnostic.pdf"),
    width = 8,
    height = 6
  )
  sc <- autoEstCont(sc)
  dev.off()

  cat(
    "Estimated contamination fraction (rho):",
    round(sc$metaData$rho[1], 4),
    "\n\n"
  )

  # Adjust counts to remove ambient RNA
  cat("Adjusting counts to remove ambient RNA...\n")
  gex_counts_corrected <- adjustCounts(sc, roundToInt = TRUE)

  cat("Ambient RNA correction complete!\n")
  cat("  Original total counts:", sum(gex_counts), "\n")
  cat("  Corrected total counts:", sum(gex_counts_corrected), "\n")
  cat("  Counts removed:", sum(gex_counts) - sum(gex_counts_corrected), "\n")
  cat(
    "  Percent removed:",
    round(
      100 * (sum(gex_counts) - sum(gex_counts_corrected)) / sum(gex_counts),
      2
    ),
    "%\n\n"
  )

  # Save the corrected matrix
  cat("Saving ambient-corrected count matrix...\n")
  saveRDS(
    gex_counts_corrected,
    file.path(config$dir_output, "ambient_corrected_counts.rds")
  )

  # Save SoupChannel object for reference
  saveRDS(
    sc,
    file.path(config$dir_output, "soupX_channel.rds")
  )

  cat("✓ Ambient-corrected matrix saved\n\n")
} else {
  cat("\n==================================================\n")
  cat("Skipping SoupX ambient RNA correction (disabled in config)\n")
  cat("==================================================\n\n")
  sc <- NULL
  gex_counts_corrected <- NULL
}

## create Seurat object with original counts ------------------------------

cat("Creating Seurat object with original (non-corrected) counts...\n")
seu <- CreateSeuratObject(
  counts = gex_counts,
  assay = "RNA",
  project = config$sample_id
)

# Add cellcalling stats to seurat object
e.out_merge <- e.out[, c("LogProb", "FDR")] |> as.data.frame()
seu@meta.data[, paste0("du_", colnames(e.out_merge))] <- e.out_merge[
  rownames(seu@meta.data),
]

# Add SoupX info to metadata (if ambient correction was run)
if (isTRUE(config$ambient_correction) && !is.null(sc)) {
  seu$soupX_rho <- sc$metaData[colnames(seu), "rho"]
  seu$soupX_cluster <- sc$metaData[colnames(seu), "clusters"]
}

# Add sample ID to metadata
seu$sample_id <- config$sample_id

if (isTRUE(config$hto_demultiplex)) {
  if (is.null(hto_counts_raw)) {
    stop(
      "Error: HTO demultiplexing requested but no Antibody Capture data found"
    )
  }

  cat("HTO features:", nrow(hto_counts), "\n")
  cat("HTO barcodes:", ncol(hto_counts), "\n\n")
  # Add HTO assay
  seu[["HTO"]] <- CreateAssayObject(
    counts = hto_counts + 1
  )

  seu@active.assay <- "HTO"

  # Normalize HTO data using CLR (centered log-ratio) transformation
  cat("Normalizing HTO data...\n")
  seu <- NormalizeData(
    seu,
    assay = "HTO",
    normalization.method = "CLR"
  )

  # Check if hto_counts are identically 0, skip demultiplexing if so
  max_hto <- max(hto_counts)
  cat("Maximum HTO count:", max_hto, "\n")

  if (max_hto > 0) {
    cat("Running HTODemux...\n")

    seu <- HTODemux(
      seu,
      assay = "HTO",
      positive.quantile = 0.99
    )

    cat("Demultiplexing results:\n")
    print(table(seu$HTO_classification.global))
    print(table(seu$HTO_classification))
    cat("\n")

    ## hto umap ---------------------------------------------------------------

    cat("Computing HTO UMAP...\n")
    seu <- ScaleData(seu)
    seu <- RunUMAP(
      seu,
      features = rownames(seu),
      densmap = TRUE,
      reduction.name = "hto_umap"
    )

    cat("✓ Demultiplexing complete\n\n")
  } else {
    cat("WARNING: All HTO counts are zero - skipping demultiplexing\n\n")
  }
}


# save object ------------------------------------------------------------

cat("Saving results...\n")
saveRDS(
  seu,
  file.path(config$dir_output, "cell_call_demux.rds")
)

saveRDS(
  e.out,
  file.path(config$dir_output, "dropletutils.rds")
)

# Save summary statistics
summary_stats <- data.frame(
  sample_id = config$sample_id,
  total_barcodes = nrow(e.out),
  cells_called = sum(e.out$FDR <= fdr_threshold, na.rm = TRUE),
  fdr_threshold = fdr_threshold,
  final_cells = ncol(seu),
  ambient_correction = isTRUE(config$ambient_correction),
  total_counts_original = sum(gex_counts)
)

if (isTRUE(config$ambient_correction)) {
  summary_stats$mean_contamination_fraction <- mean(sc$metaData$rho)
  summary_stats$total_counts_corrected <- sum(gex_counts_corrected)
  summary_stats$counts_removed <- sum(gex_counts) - sum(gex_counts_corrected)
  summary_stats$percent_removed <- 100 *
    (sum(gex_counts) - sum(gex_counts_corrected)) /
    sum(gex_counts)
}

write.csv(
  summary_stats,
  file.path(config$dir_output, "demux_soupX_summary.csv"),
  row.names = FALSE
)

cat("✓ All data saved\n\n")

# plots -------------------------------------------------------------------

cat("Generating diagnostic plots...\n")

## cell calling plots -----------------------------------------------------

p_bcrank <- bc_rank |>
  as.data.frame() |>
  filter(!duplicated(bc_rank$rank)) |>
  ggplot(aes(x = log10(rank), y = log10(total))) +
  ggrastr::geom_point_rast(stroke = 0) +
  geom_hline(yintercept = gex_low_umi_q01, col = "firebrick") +
  geom_hline(yintercept = log10(metadata(bc_rank)$knee), col = "dodgerblue") +
  annotate(
    "text",
    x = 0,
    y = log10(metadata(bc_rank)$knee) + 0.1,
    label = "kneepoint",
    col = "dodgerblue",
    hjust = 0
  ) +
  annotate(
    "text",
    x = 5,
    y = gex_low_umi_q01 + 0.1,
    label = "q01 dropletutils",
    col = "firebrick",
    hjust = 0
  ) +
  ggtitle(paste0(config$sample_id, " - Barcode Rank Plot"))

ggsave(
  file.path(config$dir_figures, "barcode_rank.pdf"),
  p_bcrank,
  width = 5,
  height = 5
)

p_fdr_hist <- e.out |>
  as.data.frame() |>
  filter(!is.na(FDR), Total > 300) |>
  mutate(
    min_pos = suppressWarnings(min(FDR[FDR > 0], na.rm = TRUE)),
    min_pos = ifelse(is.finite(min_pos), min_pos, 1e-300),
    FDR_adj = if_else(FDR > 0, FDR, min_pos / 10)
  ) |>
  ggplot(aes(x = FDR_adj)) +
  geom_histogram(
    bins = 30,
    color = "grey20",
    fill = "grey70",
    linewidth = 0.2
  ) +
  geom_vline(xintercept = fdr_threshold, linetype = "dashed") +
  scale_x_log10() +
  scale_y_continuous(trans = "log10") +
  labs(
    title = paste0(config$sample_id, " - FDR Histogram"),
    subtitle = "Droplets with >300 UMIs",
    x = "FDR (log10)",
    y = "Count (log10)"
  )

ggsave(
  file.path(config$dir_figures, "fdr_histogram.pdf"),
  p_fdr_hist,
  width = 5,
  height = 5
)

## HTO plots --------------------------------------------------------------

if (isTRUE(config$hto_demultiplex) && exists("hto_counts")) {
  # HTO library size plot
  p_hto_libsize <- seu@assays$HTO@counts |>
    t() |>
    as.data.frame() |>
    rownames_to_column("cell") |>
    pivot_longer(cols = !cell, names_to = "antibody") |>
    ggplot(aes(x = antibody, y = log10(value))) +
    geom_violin(scale = "width") +
    ggtitle(paste0(config$sample_id, " - HTO Library Size"))

  ggsave(
    file.path(config$dir_figures, "hto_libsize.pdf"),
    p_hto_libsize,
    width = 6,
    height = 5
  )

  if (
    max(hto_counts) > 0 && "HTO_classification" %in% colnames(seu@meta.data)
  ) {
    # Ridge plot
    n_htos <- min(4, nrow(seu[["HTO"]]))
    p_hto_ridge <- RidgePlot(
      seu,
      assay = "HTO",
      features = rownames(seu[["HTO"]])[1:n_htos],
      ncol = 2
    ) +
      ggtitle(paste0(config$sample_id, " - HTO Ridge Plot"))

    ggsave_rds(
      file.path(config$dir_figures, "hto_ridge.pdf"),
      p_hto_ridge,
      width = 9,
      height = 9
    )

    # HTO statistics barplot
    p_hto_stats <- seu@meta.data |>
      select(HTO_classification, HTO_classification.global) |>
      mutate(
        HTO = if_else(
          HTO_classification.global != "Doublet",
          HTO_classification,
          "doublet"
        )
      ) |>
      group_by(HTO) |>
      summarise(n = n()) |>
      ungroup() |>
      ggplot(aes(x = HTO, y = n)) +
      geom_col() +
      geom_text(aes(label = n), vjust = -0.5) +
      ggtitle(paste0(config$sample_id, " - HTO Classification Counts"))

    ggsave_rds(
      file.path(config$dir_figures, "hto_stats.pdf"),
      p_hto_stats,
      width = 9,
      height = 9
    )

    # HTO UMAP plots
    p_hto_umap1 <- DimPlot(
      seu,
      reduction = "hto_umap",
      group.by = "HTO_classification.global"
    ) +
      scale_color_canva() +
      ggtitle(paste0(config$sample_id, " - HTO Global Classification"))

    p_hto_umap2 <- DimPlot(
      seu,
      reduction = "hto_umap",
      group.by = "HTO_classification"
    ) +
      scale_color_tableau("Tableau 20") +
      ggtitle(paste0(config$sample_id, " - HTO Classification"))

    p_hto_umap <- p_hto_umap1 + p_hto_umap2
    ggsave_rds(
      file.path(config$dir_figures, "hto_umap.pdf"),
      p_hto_umap,
      width = 10,
      height = 6
    )
  }
}

cat("✓ All plots saved\n\n")

# done --------------------------------------------------------------------

cat("==================================================\n")
cat("Processing complete for sample:", config$sample_id, "\n")
cat("==================================================\n")
cat("Output files saved to:", config$dir_output, "\n")
cat("Figures saved to:", config$dir_figures, "\n")
cat("Final cell count:", ncol(seu), "\n")
if (isTRUE(config$ambient_correction)) {
  cat("Mean contamination fraction:", round(mean(sc$metaData$rho), 4), "\n")
}
cat("==================================================\n")
