set.seed(0110)
options(future.globals.maxSize = 1000 * 1024^2)
options(expressions = 500000)

library(tidyverse, quietly = T, verbose = F, warn.conflicts = F)

library(Seurat, quietly = T, verbose = F, warn.conflicts = F)
library(SingleCellExperiment, quietly = T, verbose = F, warn.conflicts = F)
library(scran, quietly = T, verbose = F, warn.conflicts = F)
library(edgeR, quietly = T, verbose = F, warn.conflicts = F)

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
    config_name = "02_celltype_DEG_edgeR",
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
      default = "02_celltype_DEG_edgeR",
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


# convert to SingleCellExperiment ----------------------------------------

counts_matrix <- GetAssayData(seu, slot = "counts", assay = "RNA")
metadata <- seu@meta.data

sce <- SingleCellExperiment(
  assays = list(counts = counts_matrix),
  colData = metadata
)


# create pseudobulk samples ----------------------------------------------

# Aggregate counts across cells with same combination of celltype and sample
summed <- aggregateAcrossCells(
  sce,
  id = colData(sce)[, c(config$celltype_column, config$sample_column)]
)


# pseudobulk DE analysis loop --------------------------------------------

# Get thresholds with defaults
fdr_threshold <- config$de_thresholds$fdr %||% 0.05
logfc_threshold <- config$de_thresholds$logfc %||% 1
top_genes_n <- config$volcano$top_genes_label %||% 10
min_cells <- config$min_cells %||% 10

volcano_plots <- list()
ma_plots <- list()
de_results_list <- list()

for (celltype in unique(colData(summed)[[config$celltype_column]])) {
  cat("Processing celltype:", celltype, "\n")

  # Subset to current celltype
  current <- summed[, colData(summed)[[config$celltype_column]] == celltype]

  # Filter out samples with too few cells
  discarded <- current$ncells < min_cells
  if (all(discarded)) {
    cat("  Skipping - insufficient cells in all samples\n")
    next
  }

  current <- current[, !discarded]

  # Check if we have samples from both conditions
  conditions_present <- unique(colData(current)[[config$condition_column]])
  if (length(conditions_present) < 2) {
    cat("  Skipping - only one condition present\n")
    next
  }

  # Create DGEList object
  y <- DGEList(counts(current), samples = colData(current))

  # Filter lowly expressed genes
  keep <- filterByExpr(y, group = colData(current)[[config$condition_column]])
  y <- y[keep, ]

  if (sum(keep) < 10) {
    cat("  Skipping - too few genes after filtering\n")
    next
  }

  # Calculate normalization factors
  y <- calcNormFactors(y)

  # Set control as reference level
  colData(current)[[config$condition_column]] <- factor(
    colData(current)[[config$condition_column]],
    levels = c(config$control_condition, config$case_condition)
  )

  # Create design matrix
  design <- model.matrix(
    as.formula(paste0("~", config$condition_column)),
    data = colData(current)
  )

  # Estimate dispersions
  y <- estimateDisp(y, design)

  # Fit GLM
  fit <- glmQLFit(y, design, robust = TRUE)

  # Test for DE
  res <- glmQLFTest(fit, coef = ncol(design))

  # Extract results
  tt <- topTags(res, n = Inf)$table
  tt$gene <- rownames(tt)

  # Save results
  write_csv(tt, paste0(config$dir_output, "/", celltype, ".csv"))

  de_results_list[[celltype]] <- tt

  ## volcano plot --------------------------------------------------------

  volcano_data <- tt |>
    mutate(
      neg_log10_p = -log10(FDR),
      significant = FDR < fdr_threshold & abs(logFC) > logfc_threshold,
      direction = case_when(
        FDR < fdr_threshold & logFC > logfc_threshold ~ "Up",
        FDR < fdr_threshold & logFC < -logfc_threshold ~ "Down",
        TRUE ~ "NS"
      )
    )

  notable_genes <- volcano_data |>
    filter(direction != "NS") |>
    slice_min(FDR, n = top_genes_n) |>
    pull(gene)

  p_volcano <- ggplot(
    volcano_data,
    aes(x = logFC, y = neg_log10_p, color = direction)
  ) +
    geom_point(alpha = 0.75, size = 1) +
    ggrepel::geom_text_repel(
      data = volcano_data |> filter(gene %in% notable_genes),
      aes(label = gene),
      size = 4,
      max.overlaps = 20
    ) +
    scale_color_manual(
      values = c("Up" = "red", "Down" = "blue", "NS" = "grey"),
      name = "Expression"
    ) +
    geom_hline(
      yintercept = -log10(fdr_threshold),
      linetype = "dashed",
      color = "black"
    ) +
    geom_vline(
      xintercept = c(-logfc_threshold, logfc_threshold),
      linetype = "dashed",
      color = "black"
    ) +
    labs(
      title = paste0("Pseudobulk Volcano: ", celltype),
      subtitle = paste0("n_samples = ", ncol(current)),
      x = "Log2 Fold Change",
      y = "-log10(FDR)"
    ) +
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )

  volcano_plots[[celltype]] <- p_volcano

  ggsave(
    paste0(config$dir_figures, "/", gsub(" ", "_", celltype), "_volcano.pdf"),
    p_volcano,
    width = 8,
    height = 6
  )

  ## MA plot -------------------------------------------------------------

  p_ma <- plotMD(y, column = 1, main = paste0(celltype, ": MA plot"))
  ma_plots[[celltype]] <- recordPlot()

  cat("  Completed - found", sum(volcano_data$direction != "NS"), "DEGs\n")
}


# summary statistics -----------------------------------------------------

summary_df <- map_dfr(names(de_results_list), function(ct) {
  res <- de_results_list[[ct]]
  tibble(
    celltype = ct,
    n_genes_tested = nrow(res),
    n_sig = sum(res$FDR < fdr_threshold, na.rm = TRUE),
    n_up = sum(
      res$FDR < fdr_threshold & res$logFC > logfc_threshold,
      na.rm = TRUE
    ),
    n_down = sum(
      res$FDR < fdr_threshold & res$logFC < -logfc_threshold,
      na.rm = TRUE
    )
  )
})

write_csv(summary_df, paste0(config$dir_output, "/pseudobulk_summary.csv"))
print(summary_df)


# DEG summary bubble plot ------------------------------------------------

# Calculate cell counts and abundance changes per celltype
celltype_stats <- metadata |>

  group_by(
    celltype = .data[[config$celltype_column]],
    condition = .data[[config$condition_column]]
  ) |>
  summarise(n_cells = n(), .groups = "drop") |>
  pivot_wider(
    names_from = condition,
    values_from = n_cells,
    values_fill = 0
  )

# Calculate percentage change in abundance (case vs control)
celltype_stats <- celltype_stats |>
  mutate(
    control_cells = .data[[config$control_condition]],
    case_cells = .data[[config$case_condition]],
    abundance_logfc = log2(case_cells / control_cells)
  )

# Merge with DEG summary
bubble_data <- summary_df |>
  left_join(celltype_stats, by = "celltype") |>
  filter(!is.na(control_cells))

p_bubble <- ggplot(
  bubble_data,
  aes(
    x = control_cells,
    y = n_sig,
    color = abundance_logfc,
    size = abs(abundance_logfc)
  )
) +
  geom_point(alpha = 0.8) +
  ggrepel::geom_text_repel(
    aes(label = celltype),
    size = 3,
    max.overlaps = 20,
    box.padding = 0.5
  ) +
  scale_color_gradient2(
    low = "blue",
    mid = "grey50",
    high = "red",
    midpoint = 0,
    name = "change\nin abundance (log2)"
  ) +
  scale_size_continuous(
    name = "|% change|",
    range = c(2, 8)
  ) +
  scale_x_log10() +
  labs(
    title = "DEG counts per cell type",
    subtitle = paste0(
      "Color: abundance change (",
      config$case_condition,
      " vs ",
      config$control_condition,
      ")"
    ),
    x = paste0("# cells in ", config$control_condition, " (log10)"),
    y = "# significant DEGs"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )
p_bubble
ggsave(
  paste0(config$dir_figures, "/DEG_summary_bubble.pdf"),
  p_bubble,
  width = 6,
  height = 5
)


# DEG up/down bar plot ---------------------------------------------------

bar_data <- summary_df |>
  select(celltype, n_up, n_down) |>
  pivot_longer(
    cols = c(n_up, n_down),
    names_to = "direction",
    values_to = "n_genes"
  ) |>
  mutate(
    direction = factor(
      direction,
      levels = c("n_up", "n_down"),
      labels = c("Upregulated", "Downregulated")
    )
  )

p_bar <- ggplot(
  bar_data,
  aes(x = n_genes, y = reorder(celltype, n_genes), fill = direction)
) +
  geom_col(position = "dodge", alpha = 0.8) +
  scale_fill_manual(
    values = c("Upregulated" = "red", "Downregulated" = "blue"),
    name = "Direction"
  ) +
  labs(
    title = "Differentially expressed genes per cell type",
    subtitle = paste0(config$case_condition, " vs ", config$control_condition),
    x = "# DEGs",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "top"
  )

p_bar

ggsave(
  paste0(config$dir_figures, "/DEG_summary_bar.pdf"),
  p_bar,
  width = 6,
  height = 4
)

# save results -----------------------------------------------------------

saveRDS(
  de_results_list,
  paste0(config$dir_output, "/pseudobulk_de_results.rds")
)

cat("\nAnalysis complete. Results saved to:", config$dir_output, "\n")
