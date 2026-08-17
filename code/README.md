# code/

All template scripts and their helpers, organized by module:

```
code/
  01_preprocessing/   cell calling, QC/normalize/DR, merge
  02_analysis/        mapping, cell cycle, DEG, enrichment, sciRED, benchmarking
  03_signalling/      NicheNet
  snippets/           small ad-hoc scripts (no config contract)
  utils/              functions.R (load_config etc.), graphics, load_config.py
```

Each script (except snippets) follows the YAML-config contract, see
`../CLAUDE.md`. After adding a script here, give it an entry in
`../registry/templates.yaml`. See `../examples/run_over_samples.sh` for how
a per-sample stage is fanned out over a project's samples.
