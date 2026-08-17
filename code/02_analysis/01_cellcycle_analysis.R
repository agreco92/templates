set.seed(0110)
options(future.globals.maxSize = 1000 * 1024^2)
options(expressions = 500000)
library(reticulate, quietly = T, verbose = F, warn.conflicts = F)
# python interpreter comes from RETICULATE_PYTHON (see .Renviron.example)

library(tidyverse, quietly = T, verbose = F, warn.conflicts = F)

library(Seurat, quietly = T, verbose = F, warn.conflicts = F)
library(tricycle, quietly = T, verbose = F, warn.conflicts = F)

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
    config_name = "03_cellcycle_analysis",
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
      default = "01_cellcycle_analysis",
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
  script_name = args$config_name,
  sample_id = NULL
)

cat("Configuration loaded:\n")
cat(yaml::as.yaml(config))


## resolve output directories --------------------------------------------

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

seu <- readRDS(config$input_path)

# Get column names from config (NULL if not specified)
group_column <- config$group_column
sample_column <- config$sample_column
condition_column <- config$condition_column
control_condition <- config$control_condition

# Tricycle thresholds for G1 classification
g1_lower <- config[["tricycle_thresholds"]][["g1_lower"]] %||% 1.5
g1_upper <- config[["tricycle_thresholds"]][["g1_upper"]] %||% 5.0

# Minimum cells per group for plotting
min_cells_trend <- config$min_cells_trend %||% 500


# load cell cycle genes --------------------------------------------------

# Support custom paths or default location
g2m_genes_path <- config$g2m_genes_path %||%
  paste0(args$dir, "/data/other/g2m.genes")
s_genes_path <- config$s_genes_path %||%
  paste0(args$dir, "/data/other/s.genes")

g2m_genes <- readRDS(g2m_genes_path)
s_genes <- readRDS(s_genes_path)


# score cell cycle -------------------------------------------------------

seu <- AddModuleScore(seu, list(G2M = g2m_genes, S = s_genes))
seu@meta.data <- seu@meta.data |>
  dplyr::rename(G2M_score = Cluster1, S_score = Cluster2)

seu$tricycle_position <- tricycle::estimate_cycle_position(
  seu@assays$RNA$data,
  gname.type = "SYMBOL"
)

# Classify cells as cycling or G1
seu$tricycle_cycling <- "cycling"
seu$tricycle_cycling[seu$tricycle_position < g1_lower] <- "G1"
seu$tricycle_cycling[seu$tricycle_position > g1_upper] <- "G1"


# save annotated object --------------------------------------------------

saveRDS(seu, file = paste0(config$dir_output, "/cc_annotated.rds"))


# plots ------------------------------------------------------------------

## cell cycle trend plot -------------------------------------------------

# Filter to control condition if specified
if (!is.null(control_condition) && !is.null(condition_column)) {
  trend_data <- seu@meta.data |>
    filter(.data[[condition_column]] == control_condition)
} else {
  trend_data <- seu@meta.data
}

# Build faceting formula based on available columns
facet_vars <- c(sample_column, group_column)
facet_vars <- facet_vars[!sapply(facet_vars, is.null)]

# Build grouping columns (excluding NULLs)
grouping_cols <- c(group_column, sample_column)
grouping_cols <- grouping_cols[!sapply(grouping_cols, is.null)]

# Apply min_cells filter only if group_column is specified
if (!is.null(group_column)) {
  trend_data <- trend_data |>
    add_count(.data[[group_column]]) |>
    filter(n >= min_cells_trend)
}

p_cellcycle_trends <- trend_data |>
  pivot_longer(
    cols = c(S_score, G2M_score),
    names_to = "Phase",
    values_to = "score"
  ) |>
  group_by(Phase, across(all_of(grouping_cols))) |>
  mutate(bin = cut_interval(tricycle_position, n = config$nbins_trend)) |>
  group_by(Phase, across(all_of(grouping_cols)), bin) |>
  summarise(
    position = median(tricycle_position),
    score = median(score),
    n = n(),
    .groups = "drop"
  ) |>
  ggplot(aes(x = position, y = score)) +
  geom_point(aes(col = Phase, size = n), alpha = .5, stroke = 0) +
  geom_smooth(aes(col = Phase), se = F) +
  geom_vline(xintercept = g1_lower, col = "grey70") +
  geom_vline(xintercept = g1_upper, col = "grey70") +
  {
    if (length(facet_vars) > 0) facet_wrap(facet_vars) else NULL
  } +
  scale_color_wsj() +
  labs(
    title = "Cell cycle phase scores across tricycle position",
    x = "Tricycle position",
    y = "Phase score"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  paste0(config$dir_figures, "/cell_cycle_trend.pdf"),
  p_cellcycle_trends,
  width = 12,
  height = 8
)

cat("\nAnalysis complete. Results saved to:", config$dir_output, "\n")
