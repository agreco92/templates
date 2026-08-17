# Attach reference mapping results to Seurat object
# This snippet loads mapping metadata and dimensional reductions
# and attaches them to an existing Seurat object

# Define paths to mapping results
lca_mapping_dir <- "./data/02_differential_exp/00_map_to_LCA/"
lca_mapping_df <- paste0(lca_mapping_dir, "ref_mapped_metadata.rds")
lca_mapping_umap <- paste0(lca_mapping_dir, "ref_mapped_umap.rds")
lca_mapping_pca <- paste0(lca_mapping_dir, "ref_mapped_pca.rds")

# Load Seurat object
seu <- readRDS(seu_path)

# Load and attach mapping metadata
lca_mapping_result <- readRDS(lca_mapping_df)
seu@meta.data <- lca_mapping_result

# Load and attach dimensional reductions
lca_umap <- readRDS(lca_mapping_umap)
lca_pca <- readRDS(lca_mapping_pca)

seu@reductions[["lca_umap"]] <- lca_umap
seu@reductions[["lca_pca"]] <- lca_pca
