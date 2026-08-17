glue_recursive <- function(x, env) {
  if (is.list(x)) {
    # Recursively process each element if it's a list
    return(lapply(x, glue_recursive, env))
  } else if (is.character(x)) {
    # Only apply glue() to single strings, not vectors of strings
    if (length(x) == 1) {
      # Try-catch to handle potential glue errors
      tryCatch(
        {
          return(glue::glue(x, .envir = env))
        },
        error = function(e) {
          # If glue fails, return the original string
          message("Warning: Could not glue string: ", x)
          return(x)
        }
      )
    } else if (length(x) > 1) {
      # For character vectors, process each element separately
      return(sapply(x, function(s) glue_recursive(s, env)))
    } else {
      # Empty character vector
      return(x)
    }
  } else {
    # Return non-character values as is (e.g., booleans, numbers)
    return(x)
  }
}

load_config_old <- function(
  project_dir,
  module_name,
  script_name,
  sample_id = NULL
) {
  # Load global config
  config <- jsonlite::read_json(file.path(project_dir, "params.json"))
  config$projectDir <- project_dir

  # Load module-specific YAML config
  yaml_path <- file.path(project_dir, "config", paste0(module_name, ".yaml"))

  if (!file.exists(yaml_path)) {
    warning("Module config file not found: ", yaml_path)
    return(config)
  }

  yaml_config <- yaml::read_yaml(yaml_path)

  # Load script-specific config
  if (script_name %in% names(yaml_config)) {
    script_config <- yaml_config[[script_name]] %>%
      glue_recursive(env = list2env(config))

    config[names(script_config)] <- script_config
  }

  # Load sample-specific config overrides if sample_id provided
  if (!is.null(sample_id) && sample_id %in% names(yaml_config)) {
    if (script_name %in% names(yaml_config[[sample_id]])) {
      sample_config <- yaml_config[[sample_id]][[script_name]] %>%
        glue_recursive(env = list2env(config))

      config[names(sample_config)] <- sample_config
    }
  }

  return(config)
}

load_config <- function(
  project_dir,
  module_name,
  script_name,
  sample_id = NULL
) {
  # Load global config (obsolete, will remove in favour of inheriting project_dir from script args (see 01_normalize_dr))
  # config <- jsonlite::read_json(file.path(project_dir, "params.json"))
  config <- c()
  config$projectDir <- project_dir
  config$sample_id <- sample_id

  # Load module-specific YAML config
  yaml_path <- file.path(project_dir, "config", paste0(module_name, ".yaml"))

  if (!file.exists(yaml_path)) {
    warning("Module config file not found: ", yaml_path)
    return(config)
  }

  yaml_config <- yaml::read_yaml(yaml_path)

  # Load script-specific config
  if (script_name %in% names(yaml_config)) {
    script_config <- yaml_config[[script_name]] %>%
      glue_recursive(env = list2env(config))

    config[names(script_config)] <- script_config
  }

  # Load sample-specific config overrides if sample_id provided
  if (!is.null(sample_id) && sample_id %in% names(yaml_config)) {
    if (script_name %in% names(yaml_config[[sample_id]])) {
      sample_config <- yaml_config[[sample_id]][[script_name]] %>%
        glue_recursive(env = list2env(config))

      config[names(sample_config)] <- sample_config
    }
  }

  return(config)
}
# Save plot as both image and RDS
ggsave_rds <- function(filename, plot = last_plot(), ...) {
  # Save image
  ggsave(filename, plot, ...)

  # Save RDS
  rds_filename <- sub("\\.[^.]+$", ".rds", filename)
  saveRDS(plot, rds_filename)

  invisible(plot)
}

pool_normalize <- function(
  seurat_obj,
  do_clustering = TRUE,
  cluster_col = NULL,
  batch_col = NULL,
  assay = "RNA",
  min.mean = 0.1,
  sf_exclude = NULL
) {
  # Load required libraries
  require(scuttle)
  require(scran)
  require(batchelor)
  require(SingleCellExperiment)

  # ========== INPUT VALIDATION ==========
  # Issue warning if both do_clustering and cluster_col are specified
  if (do_clustering && !is.null(cluster_col)) {
    warning(
      "Both do_clustering=TRUE and cluster_col are specified. ",
      "Using provided cluster_col and setting do_clustering=FALSE."
    )
    do_clustering <- FALSE
  }

  # Validate cluster_col exists if provided
  if (!is.null(cluster_col)) {
    if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
      stop(sprintf("cluster_col '%s' not found in metadata", cluster_col))
    }
  }

  # Validate batch_col exists if provided
  if (!is.null(batch_col)) {
    if (!batch_col %in% colnames(seurat_obj@meta.data)) {
      stop(sprintf("batch_col '%s' not found in metadata", batch_col))
    }
  }

  # ========== EXTRACT COUNTS ==========
  counts <- seurat_obj@assays[[assay]]$counts
  sce <- SingleCellExperiment(assays = list(counts = counts))
  # colnames(sce) <- colnames(seurat_obj)

  # ========== GENES FOR SIZE-FACTOR ESTIMATION ==========
  # sf_exclude (a regex, e.g. "^mt-|^Rp[sl]") drops genes from quickCluster and
  # computePooledFactors ONLY. logNormCounts still normalises ALL genes using the
  # resulting per-cell size factors, and the counts matrix is left untouched
  # (mito/ribo QC metrics stay computable). NULL keeps subset.row = NULL, the
  # documented default, so behaviour is unchanged when the argument is not used.
  sf_subset <- NULL
  if (!is.null(sf_exclude)) {
    drop <- grepl(sf_exclude, rownames(sce))
    if (sum(!drop) < 2) {
      stop(sprintf(
        "sf_exclude '%s' leaves only %d genes for size-factor estimation",
        sf_exclude,
        sum(!drop)
      ))
    }
    sf_subset <- rownames(sce)[!drop]
    message(sprintf(
      "Size-factor estimation excludes %d/%d genes (pattern: %s)",
      sum(drop),
      nrow(sce),
      sf_exclude
    ))
  }

  # ========== DETERMINE CLUSTERS ==========
  clusters <- NULL
  use_multibatchnorm <- FALSE

  # batch_col takes priority over cluster_col
  if (!is.null(batch_col)) {
    # Issue warning if cluster_col is also specified
    if (!is.null(cluster_col)) {
      warning(
        "Both batch_col and cluster_col specified. ",
        "batch_col takes priority; cluster_col will be ignored."
      )
    }

    if (do_clustering) {
      # Batch-aware clustering
      message("Using quickCluster with batch blocking by '", batch_col, "'...")
      batch_ids <- seurat_obj@meta.data[[batch_col]]
      clusters <- scran::quickCluster(
        sce,
        min.mean = min.mean,
        block = batch_ids,
        subset.row = sf_subset
      )
    } else {
      # Special case: batch without clustering -> use multiBatchNorm
      message(
        "Using multiBatchNorm (batches as clusters, no additional clustering)"
      )
      use_multibatchnorm <- TRUE
      batch_ids <- seurat_obj@meta.data[[batch_col]]
    }
  } else if (!is.null(cluster_col)) {
    # Use provided clusters (only if batch_col not specified)
    message("Using provided clusters from '", cluster_col, "'")
    clusters <- seurat_obj@meta.data[[cluster_col]]
  } else if (do_clustering) {
    # Regular clustering without batch
    message("Using quickCluster...")
    clusters <- scran::quickCluster(sce, min.mean = min.mean, subset.row = sf_subset)
  } else {
    # No clustering at all - single cluster
    message("No clustering: treating all cells as a single cluster")
    clusters <- rep(1, ncol(sce))
  }

  # ========== COMPUTE SIZE FACTORS AND NORMALIZE ==========
  if (use_multibatchnorm) {
    # Use multiBatchNorm for the special case
    sce <- batchelor::multiBatchNorm(
      sce,
      batch = batch_ids,
      subset.row = sf_subset,
      normalize.all = TRUE
    )
    size_factors <- BiocGenerics::sizeFactors(sce)
    logcounts_mat <- SummarizedExperiment::assay(sce, "logcounts")
  } else {
    # Standard workflow: compute pooled size factors
    message("Computing pooled size factors...")
    sce <- scuttle::computePooledFactors(
      sce,
      clusters = clusters,
      min.mean = min.mean,
      subset.row = sf_subset
    )

    # Log-normalize
    sce <- scuttle::logNormCounts(sce)
    size_factors <- BiocGenerics::sizeFactors(sce)
    logcounts_mat <- SummarizedExperiment::assay(sce, "logcounts")
  }

  # ========== UPDATE SEURAT OBJECT ==========
  # Store normalized data
  seurat_obj@assays[[assay]]$data <- logcounts_mat

  # Store size factors with proper names
  names(size_factors) <- colnames(seurat_obj)
  seurat_obj$size_factors <- size_factors

  # Store normalization clusters in metadata
  if (use_multibatchnorm) {
    # For multiBatchNorm, clusters are the batch IDs
    seurat_obj$pool_norm_clusters <- as.character(batch_ids)
  } else if (!is.null(clusters)) {
    # For all other methods, store the clusters used
    seurat_obj$pool_norm_clusters <- as.character(clusters)
  } else {
    # Should not happen, but set to "1" if no clusters somehow
    seurat_obj$pool_norm_clusters <- "1"
  }

  message("Normalization complete!")
  return(seurat_obj)
}

subset_using_list <- function(
  seurat_obj,
  keep_list = NULL,
  remove_list = NULL
) {
  # Start with all cells selected
  cells_to_keep <- rep(TRUE, ncol(seurat_obj))

  # Process keep_list: Keep cells matching ALL column conditions (AND across columns)
  # Within each column, match ANY value (OR within column using %in%)
  if (!is.null(keep_list)) {
    for (col_name in names(keep_list)) {
      values <- keep_list[[col_name]]
      # AND across columns: cells must match at least one value per column
      cells_to_keep <- cells_to_keep &
        (seurat_obj@meta.data[[col_name]] %in% values)
    }
  }

  # Process remove_list: Exclude cells matching ANY column condition
  if (!is.null(remove_list)) {
    cells_to_remove <- rep(FALSE, ncol(seurat_obj))

    for (col_name in names(remove_list)) {
      values <- remove_list[[col_name]]
      # OR logic: remove cells matching any of the conditions
      cells_to_remove <- cells_to_remove |
        (seurat_obj@meta.data[[col_name]] %in% values)
    }

    # Exclude cells marked for removal
    cells_to_keep <- cells_to_keep & !cells_to_remove
  }

  # Slice the Seurat object using logical indexing
  return(seurat_obj[, cells_to_keep])
}
