# ANCOM-BC benchmarks

Benchmarks for the R and Python implementations of ANCOM-BC and related
differential-abundance methods.

## Data and output layout

The datasets are published separately on
[Zenodo](https://doi.org/10.5281/zenodo.23224668). `fetch_data.py` reads
[`data/manifest.yaml`](data/manifest.yaml), streams each requested file into
`data/`, verifies its SHA-256 checksum, and only then moves it into place. It
uses only the Python standard library.

List the available dataset groups:

```bash
python fetch_data.py --list
```

Download everything, or only the datasets needed for a benchmark:

```bash
python fetch_data.py
python fetch_data.py synthetic
python fetch_data.py emp500 sc
```

Files whose checksums already match are skipped. Use `--force` to download
them again, or `--output-dir PATH` to use a different destination. Run
`python fetch_data.py --help` for all options.

## Data verification

After downloading, use `verify_data.py` to check that every selected file is
present and matches the SHA-256 checksum recorded in the manifest. Verification
does not access the network.

```bash
python verify_data.py # verify every dataset or selected datasets
python verify_data.py synthetic
python verify_data.py emp500 sc
python verify_data.py --data-dir PATH # verify a custom download location
```

The verifier reports valid, missing, mismatched, and unreadable files. It exits
with status `0` when every selected file is valid and status `1` if verification
fails. Use `--quiet` to hide successful files.

Dataset groups correspond to the benchmark entry points as follows:

| Dataset | Used by |
| --- | --- |
| `synthetic` | `benchmark.py`, `benchmark.r` |
| `emp500` | `benchmark_emp.py`, `run_emp500_benchmark.R` |
| `sc` | `benchmark_real.py`, `benchmark_real.r` |
| `pseq` | exploratory analysis in `benchmark_res.ipynb` |

The benchmark scripts still accept the legacy layout where data files live in
the repository root, but the `data/` layout is preferred.

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
python verify_data.py        # verify downloaded dataset checksums
```


