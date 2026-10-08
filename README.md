# ANCOM-BC benchmarks

Benchmarks for the R and Python implementations of ANCOM-BC and related
differential-abundance methods.

## Data and output layout

The datasets are published separately on Zenodo. Download the files listed in
[`data/manifest.yaml`](data/manifest.yaml) into `data/`. The scripts still
accept the legacy layout where data files live in the repository root, but the
`data/` layout is preferred.

All generated checkpoints, per-run results, summaries, and plots are written
to `results/`. Paths are resolved relative to the repository rather than the
shell's current working directory, so scripts can be launched from elsewhere.
CSV checkpoints are replaced atomically after each run.

## Entry points

```text
python benchmark.py          # simulated-data Python benchmark
Rscript benchmark.r          # simulated-data R benchmark
python benchmark_real.py     # single-cell Python benchmark
Rscript benchmark_real.r     # single-cell R benchmark
python benchmark_emp.py      # EMP500 Python benchmark
Rscript run_emp500_benchmark.R  # EMP500 R-method benchmark
```

The EMP500 R runner reads `feature_emp500_subset.csv` and
`metadata_emp500_subset.csv` directly; no intermediate RDS input is required.
