set.seed(0110)
library(reticulate, quietly = T, verbose = F)
# python interpreter comes from RETICULATE_PYTHON (see .Renviron.example)

seu <- readRDS(
  "./data/03_cdaa_response_analysis/01_cellcycle_analysis/cc_annotated.rds"
)
dir_output <- ("./data/03_cdaa_response_analysis/01_cellcycle_analysis")

scCustomize::as.anndata(
  seu,
  file_path = dir_output,
  file_name = "cc_annotated.h5ad"
)
