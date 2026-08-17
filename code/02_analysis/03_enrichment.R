set.seed(0110)

suppressPackageStartupMessages({
  library(tidyverse, quietly = T, verbose = F)

  library(msigdbr, quietly = T, verbose = F)
  library(clusterProfiler, quietly = T, verbose = F)
  library(enrichplot, quietly = T, verbose = F)

  library(ggplot2, quietly = T, verbose = F)
  library(cowplot, quietly = T, verbose = F)

  library(optparse, quietly = T, verbose = F)
  library(yaml, quietly = T, verbose = F)
  library(glue, quietly = T, verbose = F)
})

# command line arguments -------------------------------------------------

is_interactive <- interactive()

if (is_interactive) {
  args <- list(
    dir = getwd(),
    output = NULL,
    sample_id = NULL,
    config_name = "03_enrichment",
    module_name = "02_analysis"
  )
} else {
  option_list <- list(
    make_option(
      c("-d", "--dir"),
      type = "character",
      default = getwd(),
      help = "working directory (all paths are relative to this)",
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
      help = "sample/analysis name (optional, for sample-specific config overrides)",
      metavar = "character"
    ),
    make_option(
      c("-c", "--config_name"),
      type = "character",
      default = "03_enrichment",
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

# load config files ------------------------------------------------------

config <- load_config(
  project_dir = args$dir,
  module_name = args$module_name,
  script_name = args$config_name,
  sample_id = args$sample_id
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

# helper functions -------------------------------------------------------

#' Safe dotplot
#'
#' Returns a dotplot or an empty placeholder if enrichment result is NULL.
safe_dotplot <- function(enrich_obj, ...) {
  if (is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0) {
    ggplot() +
      annotate(
        "text",
        x = 0.5,
        y = 0.5,
        label = "No enriched terms",
        size = 5
      ) +
      theme_void() +
      theme(plot.background = element_rect(fill = "white", color = NA))
  } else {
    dotplot(enrich_obj, ...)
  }
}

#' Safe treeplot with termsim computation
#'
#' Computes pairwise term similarity and generates a treeplot.
#' Returns an empty placeholder plot if there are no enriched terms.
safe_treeplot <- function(gsea_obj, ...) {
  tryCatch(
    {
      obj <- pairwise_termsim(gsea_obj)
      treeplot(obj, ...)
    },
    error = function(e) {
      ggplot() +
        annotate(
          "text",
          x = 0.5,
          y = 0.5,
          label = "No enriched terms",
          size = 5
        ) +
        theme_void() +
        theme(plot.background = element_rect(fill = "white", color = NA))
    }
  )
}

# load genesets ----------------------------------------------------------

message("Loading gene sets from MSigDB...")

# Canonical pathways (M2 collection)
canonical_sets <- msigdbr::msigdbr(
  db_species = config$msigdb$db_species,
  species = config$msigdb$species,
  collection = config$msigdb$canonical$collection,
  subcollection = config$msigdb$canonical$subcollection
) |>
  filter(gs_subcollection != "CGP") |>
  mutate(gs_name = str_remove(gs_name, "^[^_]+_") |> str_to_lower())

# GO terms (M5 collection)
go_sets <- msigdbr::msigdbr(
  db_species = config$msigdb$db_species,
  species = config$msigdb$species,
  collection = config$msigdb$go$collection,
  subcollection = config$msigdb$go$subcollection
) |>
  mutate(gs_name = str_remove(gs_name, "^GO[A-Z]*_") |> str_to_lower())

message(
  "  Canonical pathways: ",
  n_distinct(canonical_sets$gs_name),
  " gene sets"
)
message("  GO terms: ", n_distinct(go_sets$gs_name), " gene sets")

# load DEG files ---------------------------------------------------------

# Get all CSV files from the specified directory
de_dir <- config$de_input_dir
if (is.null(de_dir)) {
  stop("No DE input directory specified in config$de_input_dir")
}

de_files <- list.files(de_dir, pattern = "\\.csv$", full.names = TRUE)

# Apply exclusion pattern if specified
if (!is.null(config$file_exclude)) {
  de_files <- de_files[!grepl(config$file_exclude, basename(de_files))]
}

if (length(de_files) == 0) {
  stop("No CSV files found in ", de_dir)
}

message("\nFound ", length(de_files), " DE file(s) to analyze")

# extract column names from config
col_gene <- config$de_columns$gene
col_logfc <- config$de_columns$logfc
col_pval <- config$de_columns$pval_adj %||% config$de_columns$pval

# extract thresholds
thresh_logfc <- config$de_thresholds$logfc %||% 0.5
thresh_pval <- config$de_thresholds$pval_adj %||% 0.05

# main loop --------------------------------------------------------------

for (de_file in de_files) {
  analysis_name <- tools::file_path_sans_ext(basename(de_file))

  message("\n", strrep("=", 60))
  message("Analyzing: ", analysis_name)
  message(strrep("=", 60))

  # Create output directories for this analysis
  dir_out_analysis <- file.path(config$dir_output, analysis_name)
  dir_fig_canonical <- file.path(config$dir_figures, analysis_name, "canonical")
  dir_fig_go <- file.path(config$dir_figures, analysis_name, "go")

  dir.create(dir_out_analysis, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir_fig_canonical, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir_fig_go, recursive = TRUE, showWarnings = FALSE)

  # load DEG data --------------------------------------------------------

  message("  Loading differential expression results...")

  gs <- read.csv(de_file) |>
    dplyr::rename(
      gene = !!sym(col_gene),
      logFC = !!sym(col_logfc),
      pval = !!sym(col_pval)
    ) |>
    dplyr::arrange(desc(logFC)) |>
    select(gene, logFC, pval)

  # Create ranked gene list for GSEA
  gene_ranked <- gs$logFC |> setNames(gs$gene)

  # Get significant up/down genes for ORA
  genes_up <- gs |>
    filter(logFC > thresh_logfc, pval < thresh_pval) |>
    pull(gene)
  genes_dn <- gs |>
    filter(logFC < (-thresh_logfc), pval < thresh_pval) |>
    pull(gene)

  message("    Total genes: ", nrow(gs))
  message(
    "    Upregulated genes (logFC > ",
    thresh_logfc,
    ", p < ",
    thresh_pval,
    "): ",
    length(genes_up)
  )
  message(
    "    Downregulated genes (logFC < -",
    thresh_logfc,
    ", p < ",
    thresh_pval,
    "): ",
    length(genes_dn)
  )

  # enrichment analysis: canonical ---------------------------------------

  message("  Running ORA on canonical pathways...")

  enrich_up_canon <- enricher(
    genes_up,
    universe = gs$gene |> unique(),
    TERM2GENE = canonical_sets[, c("gs_name", "gene_symbol")]
  )

  enrich_dn_canon <- enricher(
    genes_dn,
    universe = gs$gene |> unique(),
    TERM2GENE = canonical_sets[, c("gs_name", "gene_symbol")]
  )

  # enrichment analysis: GO ----------------------------------------------

  message("  Running ORA on GO terms...")

  enrich_up_go <- enricher(
    genes_up,
    universe = gs$gene |> unique(),
    TERM2GENE = go_sets[, c("gs_name", "gene_symbol")]
  )

  enrich_dn_go <- enricher(
    genes_dn,
    universe = gs$gene |> unique(),
    TERM2GENE = go_sets[, c("gs_name", "gene_symbol")]
  )

  # GSEA -----------------------------------------------------------------

  message("  Running GSEA...")

  gsea_canon <- GSEA(
    gene_ranked,
    TERM2GENE = canonical_sets[, c("gs_name", "gene_symbol")],
    verbose = FALSE,
    pvalueCutoff = config$gsea$pvalue_cutoff %||% 0.05
  )

  gsea_go <- GSEA(
    gene_ranked,
    TERM2GENE = go_sets[, c("gs_name", "gene_symbol")],
    verbose = FALSE,
    pvalueCutoff = config$gsea$pvalue_cutoff %||% 0.05
  )

  # save results ---------------------------------------------------------

  message("  Saving enrichment results...")

  # ORA results
  write.csv(
    as.data.frame(enrich_up_canon),
    file.path(dir_out_analysis, "ora_up_canonical.csv"),
    row.names = FALSE
  )
  write.csv(
    as.data.frame(enrich_dn_canon),
    file.path(dir_out_analysis, "ora_dn_canonical.csv"),
    row.names = FALSE
  )
  write.csv(
    as.data.frame(enrich_up_go),
    file.path(dir_out_analysis, "ora_up_go.csv"),
    row.names = FALSE
  )
  write.csv(
    as.data.frame(enrich_dn_go),
    file.path(dir_out_analysis, "ora_dn_go.csv"),
    row.names = FALSE
  )

  # GSEA results
  write.csv(
    as.data.frame(gsea_canon),
    file.path(dir_out_analysis, "gsea_canonical.csv"),
    row.names = FALSE
  )
  write.csv(
    as.data.frame(gsea_go),
    file.path(dir_out_analysis, "gsea_go.csv"),
    row.names = FALSE
  )

  # visualization: canonical ---------------------------------------------

  message("  Generating canonical pathway plots...")

  # ORA dotplots
  p_ora_up_canon <- safe_dotplot(enrich_up_canon, showCategory = 15) +
    labs(title = paste0(analysis_name, " - ORA Upregulated"))

  ggsave(
    file.path(dir_fig_canonical, "ora_up_dotplot.pdf"),
    p_ora_up_canon,
    width = 8,
    height = 6
  )

  p_ora_dn_canon <- safe_dotplot(enrich_dn_canon, showCategory = 15) +
    labs(title = paste0(analysis_name, " - ORA Downregulated"))

  ggsave(
    file.path(dir_fig_canonical, "ora_dn_dotplot.pdf"),
    p_ora_dn_canon,
    width = 8,
    height = 6
  )

  # GSEA treeplots (up and down separately)
  if (nrow(as.data.frame(gsea_canon)) > 0) {
    gsea_canon_up <- gsea_canon |> filter(NES > 0)
    gsea_canon_dn <- gsea_canon |> filter(NES < 0)

    p_gsea_tree_up_canon <- safe_treeplot(gsea_canon_up, nWords = 0) +
      labs(title = paste0(analysis_name, " - GSEA Up (Canonical)")) +
      theme(legend.position = "bottom")

    ggsave(
      file.path(dir_fig_canonical, "gsea_treeplot_up.pdf"),
      p_gsea_tree_up_canon,
      width = 10,
      height = 8
    )

    p_gsea_tree_dn_canon <- safe_treeplot(gsea_canon_dn, nWords = 0) +
      labs(title = paste0(analysis_name, " - GSEA Down (Canonical)")) +
      theme(legend.position = "bottom")

    ggsave(
      file.path(dir_fig_canonical, "gsea_treeplot_dn.pdf"),
      p_gsea_tree_dn_canon,
      width = 10,
      height = 8
    )
  }

  # visualization: GO ----------------------------------------------------

  message("  Generating GO term plots...")

  # ORA dotplots
  p_ora_up_go <- safe_dotplot(enrich_up_go, showCategory = 15) +
    labs(title = paste0(analysis_name, " - ORA Upregulated (GO)"))

  ggsave(
    file.path(dir_fig_go, "ora_up_dotplot.pdf"),
    p_ora_up_go,
    width = 8,
    height = 6
  )

  p_ora_dn_go <- safe_dotplot(enrich_dn_go, showCategory = 15) +
    labs(title = paste0(analysis_name, " - ORA Downregulated (GO)"))

  ggsave(
    file.path(dir_fig_go, "ora_dn_dotplot.pdf"),
    p_ora_dn_go,
    width = 8,
    height = 6
  )

  # GSEA treeplots (up and down separately)
  if (nrow(as.data.frame(gsea_go)) > 0) {
    gsea_go_up <- gsea_go |> filter(NES > 0)
    gsea_go_dn <- gsea_go |> filter(NES < 0)

    p_gsea_tree_up_go <- safe_treeplot(gsea_go_up, nWords = 0) +
      labs(title = paste0(analysis_name, " - GSEA Up (GO)")) +
      theme(legend.position = "bottom")

    ggsave(
      file.path(dir_fig_go, "gsea_treeplot_up.pdf"),
      p_gsea_tree_up_go,
      width = 10,
      height = 8
    )

    p_gsea_tree_dn_go <- safe_treeplot(gsea_go_dn, nWords = 0) +
      labs(title = paste0(analysis_name, " - GSEA Down (GO)")) +
      theme(legend.position = "bottom")

    ggsave(
      file.path(dir_fig_go, "gsea_treeplot_dn.pdf"),
      p_gsea_tree_dn_go,
      width = 10,
      height = 8
    )
  }

  message("  Done with ", analysis_name, "!")
}

message("\n", strrep("=", 60))
message("Enrichment analysis complete!")
message(strrep("=", 60))

# save workspace ---------------------------------------------------------

workspace_file <- file.path(config$dir_output, "enrichment_workspace.RData")
save.image(file = workspace_file)
message("Workspace saved to: ", workspace_file)
