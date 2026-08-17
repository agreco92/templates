#!/bin/bash
# Proof of concept: fan a per-sample template out over every sample in a
# project. Real multi-sample orchestration is a job for the Nextflow port;
# this is just the shell-loop version of the same idea. Copy and edit the
# four variables below for your own template/stage.
set -euo pipefail
: "${TEMPLATES2_HOME:?set TEMPLATES2_HOME (see README)}"

PROJECT_DIR="/path/to/project"
SAMPLE_DIR="$PROJECT_DIR/data/some/previous/stage"
TEMPLATE="$TEMPLATES2_HOME/code/01_preprocessing/00_cellcalling_demux_soupX.R"
CONFIG_NAME="00_cellcalling_demux"
MODULE_NAME="01_preprocessing"

for sample in $(find "$SAMPLE_DIR" -maxdepth 1 -mindepth 1 -type d -printf '%f\n' | sort); do
  echo "Processing $sample..."
  Rscript "$TEMPLATE" \
    --sample_id "$sample" \
    --dir "$PROJECT_DIR" \
    --config_name "$CONFIG_NAME" \
    --module_name "$MODULE_NAME"
done
