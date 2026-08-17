# templates2: scRNA-seq template library

templates2 is the self-contained home of the scRNA-seq processing template
library. The analysis scripts live under `code/`.

Goal: evolve this from a hand-driven template + YAML system into something that
can scaffold configs and run pipelines. A Nextflow port is the next step.

## The template contract (how the system works)

Every template is a standalone R or Python script with one calling convention:

```
Rscript <template>.R \
  --dir <projectDir> \           # project root; all relative paths resolve here
  --module_name <module> \       # which <module>.yaml file under <projectDir>/config/
  --config_name <block> \        # which top-level block inside that YAML
  --sample_id <sample> \         # optional; selects per-sample override block
  --output <dir>                 # optional; overrides dir_output/dir_figures (Nextflow-style)
```

`load_config(project_dir, module_name, script_name, sample_id)` (in
`code/utils/functions.R`, with a Python port in `code/utils/load_config.py`):

1. reads `<projectDir>/config/<module_name>.yaml`
2. pulls the `<config_name>` block
3. applies **glue interpolation**: `{projectDir}`, `{sample_id}` in any string
4. overlays `<sample_id>.<config_name>` if `sample_id` matches a top-level key
   (this is how per-sample QC thresholds are tuned, see `config_examples/liu_skin.yaml`)

Returns a flat config object. Scripts read `config$key` (R) / `config["key"]`.

## Path / data-flow conventions (the implicit DAG)

There is no explicit pipeline definition. Stages chain by **path convention**:
a stage writes to `{projectDir}/data/<module>/<step>/{sample_id}/` and the next
stage reads from that same directory. Per-sample stages fan out, one run per
sample (see `examples/run_over_samples.sh`); aggregate stages (merge, DEG,
enrichment) run once.

`registry/templates.yaml` is a plain catalog of the stages: script location,
R/Python, per-sample vs aggregate, which environment it needs, and a pointer to
the documented config block. It is an index only, parameter docs stay in
`config_examples/`.

## Gotchas

- **R partial matching**: `config$remove` would match `config$remove_doublets`.
  Config-driven `$` access is risky; prefer `config[["key"]]`.
- **utils path**: scripts source `code/utils/` from `$TEMPLATES2_HOME`, which
  the user sets in `.Renviron`. There is no detection fallback by design, an
  unset variable is a hard error. `--dir` locates the project's config and data
  only, never the helpers.
- **Python for reticulate**: resolved from `$RETICULATE_PYTHON`. Scripts no
  longer name a conda env or a micromamba binary. Only 4 of the 10 stages touch
  Python at all, and only for `FindClusters(algorithm = 4)` (leidenalg) and
  `RunUMAP(densmap = TRUE)` (umap-learn). The registry's `conda_env` field
  records which, `null` means pure R.
- **Environment files**: `envs/reticulate.yml` and `envs/scired.yml` are
  hand-written minimal specs. Do not regenerate them with `micromamba env
  export` from a working env; the shared envs carry unrelated packages.
- **Nextflow readiness**: `--output` and the per-stage `conda_env` exist for the
  port. The remaining gap is that stage *inputs* are still derived by path
  convention rather than taken as explicit staged paths.

## Layout

- `code/`: all template scripts, helpers (`utils/`), snippets.
- `config_examples/`: documented config blocks; the source of truth for parameters.
- `registry/templates.yaml`: catalog/index of the stages.
- `envs/`: conda environment files for the stages that need Python.
- `examples/`: minimal usage examples.
- `DESCRIPTION`: R package dependencies, readable by renv.
- `.Renviron.example`: the environment variables templates2 expects.

## When adding a template

Read the script under `code/<module>/`, extract its CLI defaults
(`config_name`/`module_name`), every `config$...`/`config[[...]]` access, the
input path it reads, and the output(s) it writes. Add a catalog entry to
`registry/templates.yaml` (including `conda_env`) and a documented block to
`config_examples/`.
