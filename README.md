# templates2

> Work in progress. Paths, conventions and the config contract can still change on `main`. Pin a tagged release if a project depends on specific behavior.

## Human readme

As a staff bioinformatician I kept meeting datasets that needed the same preprocessing steps, from cell calling and QC through clustering, differential expression, integration and reference mapping.

So I am collecting the R and Python scripts I have tried and tested, the ones that run almost identically across projects, and isolating their parameters in separate `.yaml` files.

As the project has grown it has become a bit more general but less immediate. In its current version a couple of environment variables need to be set by the user (I suspect me only, for the time being :) ) and the Python environments have been added in an `envs/` folder.

The next obvious steps are to turn these into proper Nextflow modules and to containerise the environments with Docker/Singularity.

------------------------------------------------------------------------

## Technical readme (thanks Claude!)

### Pinning a version

Tags mark points where the contract is stable enough to build on.

``` bash
git clone --branch v0.1.0 <repo-url>       # fresh clone at a tag
git checkout v0.1.0                        # or, in an existing clone
```

On GitHub, tags also appear under Releases as downloadable archives. `git tag -n99` lists what changed at each tag.

`v0.1.0` is the first public release. `TEMPLATES2_HOME` is required from it onward.

### How the template system works (30-second version)

Every template is a standalone R/Python script with one calling convention:

``` bash
Rscript <template>.R \
  --dir <projectDir> \        # project root; all relative paths resolve here
  --module_name <module> \    # which <projectDir>/config/<module>.yaml file
  --config_name <block> \     # which block inside that YAML
  --sample_id <sample>        # optional; selects a per-sample override block
```

`load_config()` reads `<projectDir>/config/<module>.yaml`, pulls the named block, interpolates `{projectDir}` / `{sample_id}`, and overlays a per-sample override block if one matches. The script reads everything from `config`.

A generic script plus a project YAML gives one concrete analysis. To reuse a template on a new dataset, write a new config block instead of editing the script.

The full contract is in [CLAUDE.md](CLAUDE.md).

### What's in here

| Path | What it is |
|---------------------------|---------------------------------------------|
| `code/` | All template scripts and helpers (`utils/`, snippets). |
| `examples/` | Minimal usage examples, such as fanning a template out over samples. |
| `config_examples/` | Documented config blocks, the source of truth for each template's parameters. |
| `envs/` | Conda environment files for the stages that need Python. |
| `registry/templates.yaml` | Catalog of every stage: script path, R/Python, per-sample vs run-once, environment, pointer to its config block. |
| `CLAUDE.md` | The contract: CLI args, `load_config` behavior, path conventions. |

### Running a template in a project

Single run: see above. To run a per-sample stage over many samples, the current pattern is a shell loop with one `Rscript` call per sample. See [examples/run_over_samples.sh](examples/run_over_samples.sh) for a minimal version. A runner that reads the catalog and fans out automatically is a future option.

### Setting `TEMPLATES2_HOME`

Templates load their shared helpers (`functions.R`, `graphics_setup.R`) from templates2 itself rather than from a per-project copy, so no symlink is needed. R has to be told where templates2 lives. The helpers should eventually end up in a package, but for now set `TEMPLATES2_HOME` in your `.Renviron`:

``` bash
echo 'TEMPLATES2_HOME=/path/to/templates2' >> ~/.Renviron
```

Restart R afterwards. `.Renviron.example` at the repo root lists everything templates2 expects to find in the environment.

One gotcha: R reads only the first `.Renviron` it finds, checking the working directory before your home directory, so a project-local `.Renviron` shadows `~/.Renviron` rather than adding to it.

`--dir` is used only to find the project's config and data, not the helpers.

### Next steps

- A Nextflow port, turning the catalog into a real pipeline. `--output` and the per-stage `conda_env` field already exist to support it.
- Docker/Singularity images per stage, keyed off `conda_env`.
- scverse templates, Python/scanpy-native alongside the current Seurat-centric ones.
- A public template-authoring skill, so a raw script can be conformed to the config contract without hand-editing.

###