set.seed(0110)
options(expressions = 500000)
library(reticulate, quietly = T, verbose = F, warn.conflicts = F)
# python interpreter comes from RETICULATE_PYTHON (see .Renviron.example)

library(tidyverse, quietly = T, verbose = F, warn.conflicts = F)

library(Seurat, quietly = T, verbose = F, warn.conflicts = F)
library(DropletUtils, quietly = T, verbose = F, warn.conflicts = F)
library(ggplot2, quietly = T, verbose = F, warn.conflicts = F)
library(ggrastr, quietly = T, verbose = F, warn.conflicts = F)
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
    config_name = "00_map_to_reference",
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
      default = "00_map_to_reference",
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

# load data --------------------------------------------------------------

reference <- readRDS(config$reference_path)
reference@active.assay <- "RNA"

query <- readRDS(config$query_path)
query@active.assay <- config$assay_use %||% "RNA"

# get column names from config
sample_column <- config$sample_column %||% "sample_id" # remove
batch_column <- config$batch_column # can be NULL

# remove unused assays from query
assay_use <- config$assay_use %||% "RNA"
other_assays <- names(query@assays)[names(query@assays) != assay_use]
for (assay in other_assays) {
  query[[assay]] <- NULL
}

# mapping ----------------------------------------------------------------

# if batch column is specified, split layers
if (!is.null(batch_column)) {
  query[[assay_use]] <- split(
    query[[assay_use]],
    f = query@meta.data[, batch_column, drop = TRUE]
  )
}

# Anchor finding settings
anchor_reduction <- config$ref_anchors$dr
anchor_dims <- config$ref_anchors$n_dims %||% 20

anchors <- FindTransferAnchors(
  reference = reference,
  query = query,
  dims = 1:anchor_dims,
  reference.reduction = anchor_reduction
)

# UMAP configuration for reference
# Three modes supported (in priority order):
#   1. precomputed: use existing UMAP reduction in reference (!needs model)
#   2. precomputed_nn: compute UMAP from pre-computed NN graph (e.g., WNN)
#   3. dr: compute UMAP from a dimensionality reduction

if ("precomputed" %in% names(config$ref_umap)) {
  # Mode 1: Use existing UMAP
  ref_umap_name <- config$ref_umap[["precomputed"]]
  if (!ref_umap_name %in% names(reference@reductions)) {
    stop("Precomputed UMAP '", ref_umap_name, "' not found in reference")
  }
  message("Using precomputed UMAP: ", ref_umap_name)
} else if ("precomputed_nn" %in% names(config$ref_umap)) {
  # Mode 2: Compute UMAP from NN graph
  ref_umap_name <- config$ref_umap[["name"]] %||% "umap"
  if ("dr" %in% names(config$ref_umap)) {
    warning("Both 'precomputed_nn' and 'dr' specified; using 'precomputed_nn'")
  }
  message("Computing UMAP from NN graph: ", config$ref_umap[["precomputed_nn"]])
  reference <- RunUMAP(
    reference,
    nn.name = config$ref_umap[["precomputed_nn"]],
    reduction.name = ref_umap_name,
    reduction.key = paste0(gsub("\\.", "", ref_umap_name), "_"),
    return.model = TRUE
  )
} else if ("dr" %in% names(config$ref_umap)) {
  # Mode 3: Compute UMAP from dimensionality reduction
  ref_umap_name <- config$ref_umap[["name"]] %||% "umap"
  umap_dims <- config$ref_umap[["n_dims"]] %||% anchor_dims
  message(
    "Computing UMAP from reduction: ",
    config$ref_umap[["dr"]],
    " with dims 1:",
    umap_dims
  )
  reference <- RunUMAP(
    reference,
    dims = 1:umap_dims,
    reduction = config$ref_umap[["dr"]],
    reduction.name = ref_umap_name,
    reduction.key = paste0(gsub("\\.", "", ref_umap_name), "_"),
    return.model = TRUE,
    densmap = isTRUE(config$ref_umap[["densmap"]])
  )
} else {
  stop(
    "No UMAP method specified. Set one of: ref_umap$precomputed, ref_umap$precomputed_nn, or ref_umap$dr"
  )
}

# build refdata list from config
refdata <- list()
for (name in names(config$refdata)) {
  refdata[[name]] <- config$refdata[[name]]
}

# map query to reference
query <- MapQuery(
  anchorset = anchors,
  query = query,
  reference = reference,
  refdata = refdata,
  reference.reduction = anchor_reduction,
  reduction.model = ref_umap_name
)

## extract prediction score delta for each mapping -----------------------

for (mapping_name in names(config$refdata)) {
  # skip non-categorical mappings (e.g., ADT)
  prediction_assay <- paste0("prediction.score.", mapping_name)
  if (!prediction_assay %in% names(query@assays)) {
    next
  }

  # Get the delta between highest and second highest score for each cell
  prediction_score_delta <- apply(
    query@assays[[prediction_assay]]@data,
    2,
    function(x) {
      sorted_scores <- sort(x, decreasing = TRUE)
      sorted_scores[1] - sorted_scores[2]
    }
  )

  # Rename columns: predicted.<mapping_name> -> <prefix>_<mapping_name>
  prefix <- config$output_prefix %||% "ref"
  old_pattern <- paste0("^predicted.", mapping_name)
  new_pattern <- paste0(prefix, "_", mapping_name)

  cols_rename <- grep(old_pattern, colnames(query@meta.data))
  colnames(query@meta.data)[cols_rename] <- gsub(
    old_pattern,
    new_pattern,
    colnames(query@meta.data)[cols_rename]
  )

  # Add score delta
  query@meta.data[, paste0(
    new_pattern,
    ".score_delta"
  )] <- prediction_score_delta
}

# save outputs -----------------------------------------------------------

# Save reductions
saveRDS(
  query@reductions$ref.umap,
  file = paste0(config$dir_output, "/ref_mapped_umap.rds")
)
saveRDS(
  query@reductions$ref.pca,
  file = paste0(config$dir_output, "/ref_mapped_pca.rds")
)

# Save predicted assays (e.g., ADT) if present
for (mapping_name in names(config$refdata)) {
  predicted_assay <- paste0("predicted_", mapping_name)
  if (predicted_assay %in% names(query@assays)) {
    saveRDS(
      query@assays[[predicted_assay]],
      file = paste0(config$dir_output, "/ref_mapped_", mapping_name, ".rds")
    )
  }
}

# Save metadata
saveRDS(
  query@meta.data,
  file = paste0(config$dir_output, "/ref_mapped_metadata.rds")
)

# plots ------------------------------------------------------------------

prefix <- config$output_prefix %||% "ref"

## UMAP plots ------------------------------------------------------------

# Reference vs Query UMAP
for (mapping_name in names(config$refdata)) {
  ref_col <- config$refdata[[mapping_name]]
  query_col <- paste0(prefix, "_", mapping_name)

  # skip if column doesn't exist
  if (!query_col %in% colnames(query@meta.data)) {
    next
  }
  if (!ref_col %in% colnames(reference@meta.data)) {
    next
  }

  # Create consistent color scale
  cell_types <- unique(c(
    reference@meta.data[[ref_col]],
    query@meta.data[[query_col]]
  ))
  n_colors <- min(length(cell_types), 20)
  colscale <- tableau_color_pal("Tableau 20")(n_colors)
  names(colscale) <- cell_types[1:n_colors]

  p_ref <- DimPlot(
    reference,
    reduction = ref_umap_name,
    raster = TRUE,
    pt.size = 1,
    raster.dpi = c(600, 600),
    group.by = ref_col,
    label = TRUE,
    label.size = 6 / .pt
  ) +
    theme_get() +
    scale_color_manual(values = colscale, na.value = "grey80") +
    ggtitle("Reference")

  p_query <- DimPlot(
    query,
    reduction = "ref.umap",
    pt.size = 1.25,
    raster = TRUE,
    raster.dpi = c(600, 600),
    group.by = query_col,
    label = FALSE
  ) +
    theme_get() +
    theme(legend.position = "none") +
    scale_color_manual(values = colscale, na.value = "grey80") +
    ggtitle("Query (mapped)")

  p_umap <- (p_ref | p_query) + patchwork::plot_layout(guides = "collect")

  ggsave(
    paste0(config$dir_figures, "/umap_", mapping_name, "_projection.pdf"),
    p_umap,
    width = 18,
    height = 8,
    units = "cm"
  )
}

# save config used -------------------------------------------------------

writeLines(
  yaml::as.yaml(config),
  con = file.path(config$dir_output, "config_used.yaml")
)

message("Reference mapping complete!")
