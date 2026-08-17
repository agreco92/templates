#!/usr/bin/env python3
"""
sciRED factor analysis — template module.

Runs sciRED (Single-Cell Interpretable REsidual Decomposition) to identify
factor-covariate associations in scRNA-seq data.

Pipeline:
  raw counts → Poisson GLM (lib size) → Pearson residuals →
  StandardScaler → PCA → varimax rotation → FCA + FIS

Input:  AnnData (.h5ad) with raw counts in adata.X.
Output: FCA heatmap, FIS table, rotated scores/loadings, top genes.

Usage (standalone):
    python 04_sciRED_factor_analysis.py \
        --dir /path/to/project \
        --config_name 04_sciRED_factor_analysis \
        --module_name 02_analysis

Usage (with shell wrapper):
    See config_examples/02_analysis.yaml for YAML config structure.
"""

import argparse
import os
import sys
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc
import yaml

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.pipeline import Pipeline

from sciRED import ensembleFCA as efca
from sciRED import glm
from sciRED import rotations as rot
from sciRED import metrics as met
from sciRED.utils import preprocess as proc
from sciRED.utils import visualize as vis

# ---------------------------------------------------------------------------
# CLI arguments (mirrors R template pattern)
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description="sciRED factor analysis (YAML-configured module)"
)
parser.add_argument(
    "-d", "--dir", type=str, default=os.getcwd(),
    help="Project root directory (all paths relative to this)"
)
parser.add_argument(
    "-o", "--output", type=str, default=None,
    help="Output directory (overrides config, for Nextflow compatibility)"
)
parser.add_argument(
    "-s", "--sample_id", type=str, default=None,
    help="Sample ID (for sample-level config overrides)"
)
parser.add_argument(
    "-c", "--config_name", type=str, default="04_sciRED_factor_analysis",
    help="Config block name in YAML"
)
parser.add_argument(
    "-m", "--module_name", type=str, default="02_analysis",
    help="YAML config file name (without .yaml)"
)
args = parser.parse_args()

# ---------------------------------------------------------------------------
# Load configuration
# ---------------------------------------------------------------------------
sys.path.insert(0, os.path.join(args.dir, "code", "utils"))
try:
    from load_config import load_config
except ImportError:
    # Fallback: look in templates directory
    templates_utils = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "utils"
    )
    sys.path.insert(0, templates_utils)
    from load_config import load_config

config = load_config(
    project_dir=args.dir,
    module_name=args.module_name,
    script_name=args.config_name,
    sample_id=args.sample_id,
)

print("Configuration loaded:")
print(yaml.dump(config, default_flow_style=False))

# ---------------------------------------------------------------------------
# Resolve output directories (CLI --output takes priority)
# ---------------------------------------------------------------------------
if args.output is not None:
    config["dir_output"] = args.output
    config["dir_figures"] = args.output

config.setdefault("dir_output", ".")
config.setdefault("dir_figures", ".")

os.makedirs(config["dir_output"], exist_ok=True)
os.makedirs(config["dir_figures"], exist_ok=True)

print(f"Output directory: {config['dir_output']}")
print(f"Figures directory: {config['dir_figures']}")

# ---------------------------------------------------------------------------
# Resolve parameters with defaults
# ---------------------------------------------------------------------------
input_path = config.get("input_path")
if input_path is None:
    raise ValueError("config must specify 'input_path' (path to .h5ad with raw counts)")

n_components = config.get("n_components", 30)
n_genes = config.get("n_genes", 2000)
covariates = config.get("covariates", ["age", "sex", "genotype", "batch"])
technical_covariates = config.get("technical_covariates", None)
lib_size_col = config.get("lib_size_col", "nCount_RNA")
seed = config.get("seed", 10)
n_top_genes = config.get("n_top_genes", 30)
fca_scale = config.get("fca_scale", "standard")
fca_mean = config.get("fca_mean", "arithmatic")
fca_time_eff = config.get("fca_time_eff", True)

np.random.seed(seed)

# ---------------------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------------------
# NORMALIZATION NOTE (from sciRED paper, Methods p.8-9):
# sciRED expects RAW COUNTS in adata.X — NOT log-normalized data.
# The Poisson GLM with library size as covariate IS the normalization.
# Pearson residuals are then z-scored per gene (StandardScaler) before PCA.
print(f"Loading {input_path} ...")
adata = sc.read_h5ad(input_path)
print(f"  {adata.n_obs} cells x {adata.n_vars} genes")
print(f"  Metadata columns: {list(adata.obs.columns)}")

# Sanity check: warn if data looks log-normalized
max_val = adata.X.max() if not hasattr(adata.X, "toarray") else adata.X.toarray().max()
if max_val < 20:
    print(f"  WARNING: max(adata.X) = {max_val:.1f} — looks log-normalized.")
    print("  sciRED expects raw counts. Export RNA$counts from Seurat, not RNA$data.")

# Validate covariates
for cov in covariates:
    if cov not in adata.obs.columns:
        raise ValueError(
            f"Covariate '{cov}' not in adata.obs. Available: {list(adata.obs.columns)}"
        )

# ---------------------------------------------------------------------------
# 2. Subset to top variable genes (by variance, on raw counts)
# ---------------------------------------------------------------------------
print(f"Subsetting to top {n_genes} variable genes ...")
adata_sub, gene_idx = proc.get_sub_data(adata, num_genes=n_genes)
y, genes, num_cells, num_genes_actual = proc.get_data_array(adata_sub)
print(f"  Data matrix: {num_cells} cells x {num_genes_actual} genes")

# ---------------------------------------------------------------------------
# 3. Poisson GLM — regress out library size (+ optional technical covariates)
# ---------------------------------------------------------------------------
print("Fitting Poisson GLM ...")

x_lib = proc.get_library_design_mat(adata_sub, lib_size=lib_size_col)

if technical_covariates:
    tech_mats = []
    for tc in technical_covariates:
        print(f"  Adding technical covariate: {tc}")
        x_tc = proc.get_design_mat(metadata_col=tc, data=adata_sub)
        tech_mats.append(x_tc)
    x = np.column_stack([x_lib] + tech_mats)
else:
    x = x_lib

print(f"  Design matrix shape: {x.shape}")
glm_fit = glm.poissonGLM(y, x)
resid_pearson = glm_fit["resid_pearson"]
y_resid = resid_pearson.T  # (num_cells, num_genes)
print("  Pearson residuals computed.")

# ---------------------------------------------------------------------------
# 4. PCA + varimax rotation
# ---------------------------------------------------------------------------
print(f"Running PCA ({n_components} components) + varimax rotation ...")

pipeline = Pipeline([
    ("scaling", StandardScaler()),
    ("pca", PCA(n_components=n_components))
])
pca_scores = pipeline.fit_transform(y_resid)
pca_model = pipeline.named_steps["pca"]
pca_loadings = pca_model.components_

# Scree plot
var_explained = pca_model.explained_variance_ratio_
fig, ax = plt.subplots(figsize=(8, 4))
ax.bar(range(1, n_components + 1), var_explained)
ax.set_xlabel("Factor")
ax.set_ylabel("Variance explained")
ax.set_title("Scree plot (pre-rotation)")
fig.tight_layout()
fig.savefig(os.path.join(config["dir_figures"], "scree_plot.pdf"), dpi=150)
plt.close(fig)
print(f"  Total variance explained: {var_explained.sum():.3f}")

# Varimax rotation
rotation_results = rot.varimax(pca_loadings.T)
rotated_loadings = rotation_results["rotloading"]
rotation_matrix = rotation_results["rotmat"]
rotated_scores = rot.get_rotated_scores(pca_scores, rotation_matrix)
print("  Varimax rotation complete.")

# ---------------------------------------------------------------------------
# 5. Factor-Covariate Association Table (FCAT)
# ---------------------------------------------------------------------------
print("Computing Factor-Covariate Association (FCA) table ...")

fcat_dict = {}
for cov in covariates:
    print(f"  Processing covariate: {cov}")
    cov_vec = adata_sub.obs[cov].astype(str).reset_index(drop=True)
    fcat_cov = efca.FCAT(
        covariate_vec=cov_vec,
        factor_scores=rotated_scores,
        scale=fca_scale,
        mean=fca_mean,
        time_eff=fca_time_eff,
    )
    fcat_dict[cov] = fcat_cov

fcat_all = pd.concat(fcat_dict.values(), axis=0)
fcat_all.to_csv(os.path.join(config["dir_output"], "FCA_table.csv"))
print(f"  FCA table saved: {fcat_all.shape}")

# ---------------------------------------------------------------------------
# 6. Plot FCA heatmap
# ---------------------------------------------------------------------------
print("Plotting FCA heatmaps ...")

vis.plot_FCAT(
    fcat_all,
    title="Factor-Covariate Association",
    color="RdBu_r",
    x_axis_fontsize=20, y_axis_fontsize=16, title_fontsize=18,
    x_axis_tick_fontsize=14, y_axis_tick_fontsize=14, legend_fontsize=14,
    save=True,
    save_path=os.path.join(config["dir_figures"], "FCA_heatmap.pdf"),
)
plt.close("all")

# Otsu threshold for significant associations
fcat_values = fcat_all.values.flatten()
threshold = efca.get_otsu_threshold(fcat_values)
print(f"  Otsu threshold: {threshold:.3f}")

with open(os.path.join(config["dir_output"], "FCA_threshold.txt"), "w") as f:
    f.write(f"Otsu threshold: {threshold:.4f}\n")
    for cov in covariates:
        fcat_cov = fcat_dict[cov]
        matched = (fcat_cov.values > threshold).any(axis=1)
        matched_levels = fcat_cov.index[matched].tolist()
        f.write(f"\n{cov} — significant levels:\n")
        for lev in matched_levels:
            max_factor = fcat_cov.columns[fcat_cov.loc[lev].argmax()]
            max_score = fcat_cov.loc[lev].max()
            f.write(f"  {lev}: {max_factor} (score={max_score:.3f})\n")

# Per-covariate heatmaps
for cov in covariates:
    vis.plot_FCAT(
        fcat_dict[cov],
        title=f"FCA — {cov}",
        color="RdBu_r",
        x_axis_fontsize=18, y_axis_fontsize=14, title_fontsize=16,
        x_axis_tick_fontsize=12, y_axis_tick_fontsize=12, legend_fontsize=12,
        save=True,
        save_path=os.path.join(config["dir_figures"], f"FCA_heatmap_{cov}.pdf"),
    )
    plt.close("all")

# ---------------------------------------------------------------------------
# 7. Factor Interpretability Scores (FIS)
# ---------------------------------------------------------------------------
print("Computing Factor Interpretability Scores (FIS) ...")

all_metrics = {}
all_metrics["Bimodality"] = met.bimodality_index(rotated_scores)
all_metrics["Effect size"] = met.factor_variance(rotated_scores)

for cov in covariates:
    cov_vec = adata_sub.obs[cov].astype(str).reset_index(drop=True)
    all_metrics[f"Specificity ({cov})"] = met.simpson_diversity_index(fcat_dict[cov])
    all_metrics[f"Homogeneity ({cov})"] = met.average_scaled_var(
        rotated_scores, cov_vec, mean_type="arithmetic"
    )

fis_df = pd.DataFrame(all_metrics)
fis_df.index = [f"F{i+1}" for i in range(n_components)]

fis_scaled = met.get_scaled_metrics(fis_df)
fis_scaled = pd.DataFrame(fis_scaled, index=fis_df.index, columns=fis_df.columns)
fis_scaled.to_csv(os.path.join(config["dir_output"], "FIS_table.csv"))

vis.plot_FIST(
    fis_scaled.T,
    title="Factor Interpretability Scores",
    xticks_fontsize=14, yticks_fontsize=12, legend_fontsize=12,
    save=True,
    save_path=os.path.join(config["dir_figures"], "FIS_heatmap.pdf"),
)
plt.close("all")
print("  FIS table saved.")

# ---------------------------------------------------------------------------
# 8. Save rotated scores and loadings
# ---------------------------------------------------------------------------
print("Saving factor scores and loadings ...")

factor_names = [f"F{i+1}" for i in range(n_components)]

scores_df = pd.DataFrame(rotated_scores, index=adata_sub.obs_names, columns=factor_names)
scores_df.to_csv(
    os.path.join(config["dir_output"], "rotated_factor_scores.csv.gz"),
    compression="gzip",
)

loadings_df = pd.DataFrame(rotated_loadings, index=genes, columns=factor_names)
loadings_df.to_csv(
    os.path.join(config["dir_output"], "rotated_factor_loadings.csv.gz"),
    compression="gzip",
)

# ---------------------------------------------------------------------------
# 9. Top genes per factor
# ---------------------------------------------------------------------------
print("Extracting top loaded genes per factor ...")

with open(os.path.join(config["dir_output"], "top_genes_per_factor.txt"), "w") as f:
    for i in range(n_components):
        fname = factor_names[i]
        lvec = loadings_df[fname]
        top_pos = lvec.nlargest(n_top_genes)
        top_neg = lvec.nsmallest(n_top_genes)
        f.write(f"\n{'='*60}\n{fname}\n{'='*60}\n")
        f.write("Top positive loadings:\n")
        for gene, val in zip(top_pos.index, top_pos.values):
            f.write(f"  {gene:20s}  {val:+.4f}\n")
        f.write("\nTop negative loadings:\n")
        for gene, val in zip(top_neg.index, top_neg.values):
            f.write(f"  {gene:20s}  {val:+.4f}\n")

# ---------------------------------------------------------------------------
# 10. Scaled variance tables per covariate
# ---------------------------------------------------------------------------
print("Computing scaled variance tables ...")

for cov in covariates:
    cov_vec = adata_sub.obs[cov].astype(str).reset_index(drop=True)
    sv_table = met.scaled_var_table(rotated_scores, cov_vec)
    sv_table.to_csv(os.path.join(config["dir_output"], f"scaled_variance_{cov}.csv"))

    vis.plot_relativeVar(
        sv_table,
        title=f"Relative variance — {cov}",
        x_axis_fontsize=18, y_axis_fontsize=14, title_fontsize=16,
        x_axis_tick_fontsize=12, y_axis_tick_fontsize=12, legend_fontsize=12,
        save=True,
        save_path=os.path.join(config["dir_figures"], f"relative_variance_{cov}.pdf"),
    )
    plt.close("all")

# ---------------------------------------------------------------------------
# Save config
# ---------------------------------------------------------------------------
with open(os.path.join(config["dir_output"], "config.yaml"), "w") as f:
    yaml.dump(config, f, default_flow_style=False)

print(f"\nDone. All outputs saved to: {config['dir_output']}")
