set.seed(0110)
library(Seurat)
library(tidyverse)
library(nichenetr)
library(cowplot)
library(ggpubr)

library(optparse)
library(jsonlite)
library(yaml)
library(glue)

# command line arguments -------------------------------------------------

is_interactive <- interactive()

if (is_interactive) {
  args <- list(
    dir = getwd(),
    output = NULL,
    sample_id = NULL,
    config_name = "00_nichenet_agnostic",
    module_name = "03_signalling"
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
      help = "sample name (optional, for receiver-specific runs)",
      metavar = "character"
    ),
    make_option(
      c("-c", "--config_name"),
      type = "character",
      default = "00_nichenet_agnostic",
      help = "name of config block to use from YAML",
      metavar = "character"
    ),
    make_option(
      c("-m", "--module_name"),
      type = "character",
      default = "03_signalling",
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

# Re-glue output paths with receiver variable now available
receiver <- config$receiver
config$dir_output <- glue(config$dir_output %||% ".")
config$dir_figures <- glue(config$dir_figures %||% ".")

dir.create(config$dir_output, recursive = TRUE, showWarnings = FALSE)
dir.create(config$dir_figures, recursive = TRUE, showWarnings = FALSE)


# read input -------------------------------------------------------------

seu <- readRDS(config$seu_path)

lr_network <- readRDS(config$lr_network_path) |>
  distinct(from, to)

ligand_target_matrix <- readRDS(config$ligand_target_matrix_path)

weighted_networks <- readRDS(config$weighted_networks_path)


# preprocess -------------------------------------------------------------

celltype_col <- config$celltype_col
pct_expressed_thr <- config$pct_expressed_thr %||% 0.1

Idents(seu) <- seu@meta.data[, celltype_col] |> as.factor()

expressed_genes_receiver <- get_expressed_genes(
  receiver,
  seu,
  pct = pct_expressed_thr
)

# all available receptors
all_receptors <- unique(lr_network$to)
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)

# extract potential ligands using lr net and expressed receptors
potential_ligands <- lr_network |>
  filter(to %in% expressed_receptors) |>
  pull(from) |>
  unique()

cat("Potential ligands:", length(potential_ligands), "\n")


# Define gene set of interest --------------------------------------------

seu_receiver <- subset(seu, idents = receiver)

DE_table_receiver <- read.csv(config$de_receiver_path)

de_fdr_threshold <- config$de_fdr_threshold %||% 0.05
de_logfc_threshold <- config$de_logfc_threshold %||% 0.5

geneset_oi <- DE_table_receiver |>
  filter(
    FDR <= de_fdr_threshold & (abs(logFC) >= de_logfc_threshold),
    gene %in% rownames(ligand_target_matrix)
  ) |>
  pull(gene)


# define background genes ------------------------------------------------

background_expressed_genes <- expressed_genes_receiver[
  expressed_genes_receiver %in% rownames(ligand_target_matrix)
]

cat("Background expressed genes:", length(background_expressed_genes), "\n")
cat("Gene set of interest:", length(geneset_oi), "\n")


# activity analysis ------------------------------------------------------

n_top_ligands <- config$n_top_ligands %||% 30
n_target_genes <- config$n_target_genes %||% 50
regulatory_potential_cutoff <- config$regulatory_potential_cutoff %||% 0.25

ligand_activities <- predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = background_expressed_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands
)

ligand_activities <- ligand_activities |>
  arrange(-aupr_corrected) |>
  mutate(rank = rank(desc(aupr_corrected)))

# histogram of ligand activities
p_hist_lig_activity <- ggplot(ligand_activities, aes(x = aupr_corrected)) +
  geom_histogram(color = "black", fill = "darkorange", bins = 50) +
  geom_vline(
    aes(
      xintercept = min(
        ligand_activities |>
          top_n(n_top_ligands, aupr_corrected) |>
          pull(aupr_corrected)
      )
    ),
    color = "red",
    linetype = "dashed",
    linewidth = 1
  ) +
  labs(x = "ligand activity (AUPR)", y = "# ligands") +
  coord_transform(y = "log1p") +
  theme_classic()

top_ligands <- ligand_activities |>
  top_n(n_top_ligands, aupr_corrected) |>
  arrange(-aupr_corrected) |>
  pull(test_ligand) |>
  unique()

# ligand activity heatmap
vis_ligand_aupr <- ligand_activities |>
  filter(test_ligand %in% top_ligands) |>
  column_to_rownames("test_ligand") |>
  select(aupr_corrected) |>
  arrange(aupr_corrected) |>
  as.matrix(ncol = 1)

p_ligand_aupr <- make_heatmap_ggplot(
  vis_ligand_aupr,
  "Prioritized ligands",
  "Ligand activity",
  legend_title = "AUPR",
  color = "darkorange"
) +
  theme(axis.text.x.top = element_blank())


## targets --------------------------------------------------------------
# active target genes: target genes in the gset_of_interest with highest regulatory potential

active_ligand_target_links_df <- top_ligands |>
  lapply(
    get_weighted_ligand_target_links,
    geneset = geneset_oi,
    ligand_target_matrix = ligand_target_matrix,
    n = n_target_genes
  ) |>
  bind_rows() |>
  drop_na()

cat("Active ligand-target links:", nrow(active_ligand_target_links_df), "\n")

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = regulatory_potential_cutoff
)

order_ligands <- intersect(
  top_ligands,
  colnames(active_ligand_target_links)
) |>
  rev()

order_targets <- active_ligand_target_links_df$target |>
  unique() |>
  intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets, order_ligands])

p_ligand_target <- make_heatmap_ggplot(
  vis_ligand_target,
  "Prioritized ligands",
  "Predicted target genes",
  color = "purple",
  legend_title = "Regulatory potential"
) +
  scale_fill_gradient2(low = "whitesmoke", high = "purple")


## receptors of top ligands ------------------------------------------------

ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  top_ligands,
  expressed_receptors,
  lr_network,
  weighted_networks$lr_sig
)

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  top_ligands,
  order_hclust = "both"
)

p_ligand_receptor <- make_heatmap_ggplot(
  t(vis_ligand_receptor_network),
  y_name = "Ligands",
  x_name = "Receptors",
  color = "mediumvioletred",
  legend_title = "Prior interaction potential"
)


## combined plot -----------------------------------------------------------

figures_without_legend <- plot_grid(
  p_ligand_aupr + theme(legend.position = "none"),
  p_ligand_target +
    theme(
      legend.position = "none",
      axis.title.y = element_blank()
    ),
  p_ligand_receptor +
    theme(
      legend.position = "none",
      axis.title.y = element_blank()
    ),
  align = "hv",
  nrow = 1,
  rel_widths = c(1, 3, 2)
)

legends <- plot_grid(
  as_ggplot(get_legend(p_ligand_aupr)),
  as_ggplot(get_legend(p_ligand_target)),
  as_ggplot(get_legend(p_ligand_receptor)),
  nrow = 1,
  align = "h"
)

combined_plot <- plot_grid(
  figures_without_legend,
  legends,
  rel_heights = c(10, 3),
  nrow = 2,
  align = "hv"
)


# save outputs -----------------------------------------------------------

# Save individual plots
ggsave(
  file.path(config$dir_figures, "ligand_activity_histogram.pdf"),
  p_hist_lig_activity,
  width = 8,
  height = 6
)

ggsave(
  file.path(config$dir_figures, "ligand_activity_heatmap.pdf"),
  p_ligand_aupr,
  width = 6,
  height = 10
)

ggsave(
  file.path(config$dir_figures, "ligand_target_heatmap.pdf"),
  p_ligand_target,
  width = 14,
  height = 10
)

ggsave(
  file.path(config$dir_figures, "ligand_receptor_heatmap.pdf"),
  p_ligand_receptor,
  width = 10,
  height = 8
)

ggsave(
  file.path(config$dir_figures, "combined_nichenet_summary.pdf"),
  combined_plot,
  width = 18,
  height = 12
)

# Save data outputs
saveRDS(
  ligand_activities,
  file.path(config$dir_output, "ligand_activities.rds")
)

saveRDS(
  active_ligand_target_links_df,
  file.path(config$dir_output, "active_ligand_target_links.rds")
)

saveRDS(
  ligand_receptor_links_df,
  file.path(config$dir_output, "ligand_receptor_links.rds")
)

cat("\nNicheNet analysis complete!\n")
cat("Output saved to:", config$dir_output, "\n")
cat("Figures saved to:", config$dir_figures, "\n")
