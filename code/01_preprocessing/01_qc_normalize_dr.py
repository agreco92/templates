# %% imports
import argparse
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import warnings
import yaml


# command line arguments 
parser = argparse.ArgumentParser(
    description="QC, normalization and dimensionality reduction"
)

parser.add_argument(
    "-d", "--dir",
    dest="project_dir",
    default=str(Path.cwd()),
    help="Project root containing config/"
)
parser.add_argument(
    "-o", "--output",
    default=None,
    help="Override output and figure directories"
)
parser.add_argument(
    "-s", "--sample_id",
    required=True,
    help="Sample identifier"
)
parser.add_argument(
    "-c", "--config_name",
    default="01_qc_normalize_dr"
)
parser.add_argument(
    "-m", "--module_name",
    default="01_preprocessing"
)

# %% interactive overrides
# Only active in an IPython/Positron kernel; a shell run keeps its real CLI args.
if "ipykernel" in sys.modules:
    sys.argv = [
        "01_qc_normalize_dr.py",
        "--dir", os.getcwd(),
        "--module_name", "01_preprocessing_python",
        "--sample_id", "cdaa_lsec",
    ]

args = parser.parse_args()
project_dir = Path(args.project_dir).resolve()

templates_home = os.environ.get("TEMPLATES2_HOME")


# %% load config files 
if not templates_home:
    raise RuntimeError(
        "TEMPLATES2_HOME is not set. Export it in your shell; "
        "Python does not read .Renviron."
    )

utils_dir = Path(templates_home) / "code" / "utils"
sys.path.insert(0, str(utils_dir))

from load_config import load_config

cfg = load_config(
    project_dir=str(project_dir),
    module_name=args.module_name,
    script_name=args.config_name,
    sample_id=args.sample_id,
)


print("Configuration loaded:")
print(yaml.safe_dump(cfg, sort_keys=False))

## resolve output directories

if args.output is not None:
    cfg["dir_output"] = args.output
    cfg["dir_figures"] = args.output

output_dir = Path(cfg.get("dir_output") or ".")
figures_dir = Path(cfg.get("dir_figures") or ".")
sc.settings.figdir = figures_dir

output_dir.mkdir(parents=True, exist_ok=True)
figures_dir.mkdir(parents=True, exist_ok=True)

# resolve input file

sample_dir = Path(cfg["sample_dir"])

if not sample_dir.is_absolute():
    sample_dir = project_dir / sample_dir

sample_path = sample_dir / cfg["sample_id"]
input_filename = cfg.get("input_filename")

if input_filename is None or input_filename == "auto":
    input_files = list(sample_path.glob("*.h5ad"))

    if len(input_files) == 0:
        raise FileNotFoundError(f"No .h5ad files found in {sample_path}")

    if len(input_files) > 1:
        raise RuntimeError(
            f"Multiple .h5ad files found in {sample_path}; "
            "set input_filename in the config."
        )

    input_file = input_files[0]
else:
    input_file = sample_path / input_filename

# %% QC
adata = sc.read_h5ad(input_file)

adata.layers["counts"] = adata.X.copy()
# remove empty cells
keep = np.asarray(adata.layers["counts"].sum(axis=1)).ravel() > 0

print(f"Removing {(~keep).sum()} zero-count cells")
adata = adata[keep].copy()

# TODO: drop HTO assay
# add metadata from config
if cfg.get("metadata"):
    #   message("Adding metadata from config...")
    for (key,value) in cfg["metadata"].items():
        adata.obs[key] = value

# add metadata from csv
if cfg.get("metadata_csv"):
    metadata_cfg = cfg["metadata_csv"]
    csv_path = Path(metadata_cfg["path"])

    if not csv_path.is_absolute():
        csv_path = project_dir / csv_path
    
    meta_csv = pd.read_csv(csv_path)    

    sample_col = metadata_cfg.get("sample_column")

    if sample_col is None:
        sample_col = "sample_id"

    if sample_col not in meta_csv.columns:
        raise ValueError(
            f"Column {sample_col!r} not found in metadata CSV."
        )

    matching_rows = meta_csv.loc[
        meta_csv[sample_col] == cfg["sample_id"]
    ]

    if matching_rows.empty:
        raise ValueError(
            f"sample_id {cfg['sample_id']!r} not found in "
            f"column {sample_col!r} of metadata CSV."
        )

    if len(matching_rows) > 1:
        warnings.warn(
            f"Multiple rows match sample_id {cfg['sample_id']!r}; "
            "using the first."
        )
    sample_row = matching_rows.iloc[0]
    
    cols_to_add = metadata_cfg.get("columns")
    if cols_to_add is None:
        cols_to_add = [
            col for col in meta_csv.columns
            if col!= sample_col
        ]

    for col in cols_to_add:
        if col not in meta_csv.columns:
            warnings.warn(f"Column {col!r} not found in CSV; skipping.")
            continue

        value = sample_row[col]
        print(f"  - Adding {col}: {value} (from CSV)")
        adata.obs[col] = value

## %% doublet removal
if cfg.get("remove_doublets") == True:
    sc.pp.scrublet(adata)
    adata.obs["predicted_doublet"] = adata.obs["predicted_doublet"].eq(True)

# TODO: add diagnostic plot

# %% libsize, mito, ribo
# mitochondrial genes, "MT-" for human, "Mt-" for mouse
mito_pattern = cfg.get("qc_patterns", {}).get("mito") or "^mt-"
adata.var["mt"] = adata.var["gene_symbol"].str.match(mito_pattern)
# ribosomal genes
ribo_pattern = cfg.get("qc_patterns",{}).get("ribo") or r"^RP[SL]\d+[AXY]?L?\d*$|^RPSA$|^RPLP\d$"
adata.var["ribo"] = adata.var["gene_symbol"].str.match(ribo_pattern, case=False)
# hemoglobin genes
hb_pattern = cfg.get("qc_patterns",{}).get("hb") or "^Hb[^(p)]"
adata.var["hb"] = adata.var["gene_symbol"].str.contains(hb_pattern, case=False)
adata.var["Malat1"] = adata.var["gene_symbol"].str.fullmatch("^Malat1", case=False)

sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb","Malat1"], 
    inplace=True, log1p=False)

adata.obs["log10_total_counts"] = np.log10(adata.obs["total_counts"])
adata.obs["nuclearity"] = np.log2((adata.obs["total_counts_Malat1"]+1)/(adata.obs["total_counts_ribo"]+1))

sc.pl.violin(
    adata,
    ["n_genes_by_counts", "log10_total_counts", "pct_counts_mt","nuclearity"],
    jitter=0.4,
    multi_panel=True,
    show=False,
)
plt.savefig(figures_dir / f"{cfg['sample_id']}_qc.pdf", bbox_inches="tight")
plt.close()
# remove low QC and doublets
qc_exclusion_df = pd.DataFrame({
    "low_counts": adata.obs["log10_total_counts"] < cfg["qc_thresholds"]["min_counts_log10"],
    "low_nfeatures": adata.obs["n_genes"] < cfg["qc_thresholds"]["min_features"],
    # "doublet": adata.obs["predicted_doublet"].eq(True),
    "high_mito": adata.obs["pct_counts_mt"] > cfg["qc_thresholds"]["max_mito_pct"],
    "high_nuclearity": adata.obs["nuclearity"] > cfg["qc_thresholds"]["nuclearity_above"],
    "low_nuclearity": adata.obs["nuclearity"] < cfg["qc_thresholds"]["nuclearity_below"]
})

if cfg.get("remove_doublets") == True:
    qc_exclusion_df["doublet"] = adata.obs["predicted_doublet"].eq(True)

# quick UMAP to check clustering vs QC metric
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
sc.pp.pca(adata, n_comps=20)
sc.pp.neighbors(adata)
sc.tl.umap(adata, key_added = "qc_umap")

qc_cols = ["log10_total_counts","doublet_score","pct_counts_mt",
 "nuclearity","cellbender_background_fraction", "pct_counts_ribo"]
qc_cols = [c for c in qc_cols if c in adata.obs.columns]

sc.pl.embedding(adata, color=qc_cols,
    ncols= 2, cmap = "rainbow", basis="qc_umap", show=False)
plt.savefig(figures_dir / f"{cfg['sample_id']}_qc_umap.pdf", bbox_inches="tight")
plt.close()

adata.obs["exclude_reason"] = qc_exclusion_df.apply(lambda row: ",".join(row.index[row]), axis=1)

sc.pl.embedding(adata, color="exclude_reason", basis="qc_umap", show=False)
plt.savefig(figures_dir / f"{cfg['sample_id']}_qc_umap_exclusion.pdf", bbox_inches="tight")
plt.close()

adata.obs["keep"] = ~qc_exclusion_df.any(axis=1)
adata = adata[adata.obs["keep"]].copy()

if cfg.get("exclude"):
    for key,value in cfg.get("exclude").items():
        adata = adata[~adata.obs[key].isin(np.atleast_1d(value))].copy()


# TODO: how to handle ambient-corrected counts?
# if cfg.get("ambient_corrected_counts"):
#     noamb_path = sample_path / cfg.get("ambient_corrected_counts"))
#     # noamb_path = cfg.get("sample_dir") + cfg.get("sample_id") + cfg.get("ambient_corrected_counts")
#     if noamb_path.endswith(".h5"):
#         # cellbender case
#         noamb = sc.read_10x_h5(noamb_path,gex_only=True)
#         adata = adata 
#         adata.layers["X_cellbender"] = noamb.X
# note: this does not work yet. cellbender does cellcalling too, so we need to subset the object prior to QC
# %% normalization

adata.X = sc.pp.normalize_total(
        adata,
        layer="counts",
        exclude_highly_expressed=True,
        inplace=False)["X"]
adata.uns.pop("log1p", None)
# # clustering and dr
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
sc.pp.pca(adata)
# sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pl.umap(adata, color = qc_cols, cmap = "rainbow", ncols=2, show=False)
plt.savefig(figures_dir / f"{cfg['sample_id']}_umap_qc.pdf", bbox_inches="tight")
plt.close()

# cluster vs split_by heatmap

# save object 
adata.write_h5ad(filename=output_dir / "qc_normalized.h5ad")