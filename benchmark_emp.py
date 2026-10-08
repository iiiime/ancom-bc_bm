#!/usr/bin/env python3
"""Benchmark Python differential-abundance methods on the EMP500 dataset.

The available method runners are:

  1. ANCOM-BC and ANCOM-BC2 (scikit-bio)
  2. PyDESeq2 (pydeseq2 0.5.4)
  3. edgePython (edgepython 0.2.6)
  4. pyedger (pyedger 0.1.0)

The default configuration runs ANCOM-BC and ANCOM-BC2 for three replicates.
The model formula is ``~ empo_3 + env_biome``.
"""

import argparse
import json
import resource
import subprocess
import sys
import time
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from benchmark_io import (
    atomic_pickle_dump,
    atomic_write_csv,
    input_path,
    load_pickle,
    output_path,
)

warnings.filterwarnings("ignore")


METHODS = ["ancombc", "ancombc2"]
METHOD_LABELS = {
    "ancombc": "ANCOM-BC (Python)",
    "ancombc2": "ANCOM-BC2 (Python)",
    "pydeseq2": "PyDESeq2",
    "edgepython": "edgePython",
    "pyedger": "pyedger",
}
N_REPLICATES = 3
TIMEOUTS = {
    "ancombc": 600,
    "ancombc2": 1200,
    "pydeseq2": 10800,
    "edgepython": 1800,
    "pyedger": 1800,
}
RESULTS_CSV = output_path("emp500_benchmark_py_ancombc.csv")
PREPARED_DATA = output_path("cache/emp500_prepared.pkl")


def fix_collinearity(meta_df):
    """
    Iteratively merge collinear env_biome levels into 'Other' until the
    design matrix ~empo_3 + env_biome is full rank.
    """
    from patsy import dmatrix

    for iteration in range(50):
        formula = "~ empo_3 + env_biome"
        try:
            design = dmatrix(formula, meta_df)
        except Exception:
            # patsy may fail on perfectly collinear data
            design = None

        if design is not None:
            mat = np.asarray(design, dtype=float)
            rank = np.linalg.matrix_rank(mat)
            n_cols = mat.shape[1]
            if rank == n_cols:
                print(f"  Collinearity fixed after {iteration} merges. "
                      f"Design: {n_cols} cols, rank {rank}")
                return meta_df

        # Find the env_biome level to merge
        eb = meta_df["env_biome"]
        levels = eb.cat.categories.tolist() if hasattr(eb, "cat") else sorted(eb.unique())
        merged = False
        for level in levels:
            if level == "Other":
                continue
            trial = meta_df.copy()
            # Add "Other" as a category before assigning (if not already present)
            if "Other" not in trial["env_biome"].cat.categories:
                trial["env_biome"] = trial["env_biome"].cat.add_categories(["Other"])
            trial.loc[trial["env_biome"] == level, "env_biome"] = "Other"
            trial["env_biome"] = trial["env_biome"].cat.remove_unused_categories()
            try:
                d = dmatrix("~ empo_3 + env_biome", trial)
                m = np.asarray(d, dtype=float)
                if np.linalg.matrix_rank(m) > rank:
                    meta_df = trial
                    print(f"  Merged env_biome '{level}' -> 'Other' "
                          f"(rank {rank} -> {np.linalg.matrix_rank(m)})")
                    merged = True
                    break
            except Exception:
                continue

        if not merged:
            # Force merge the smallest non-Other level
            counts_by_level = meta_df["env_biome"].value_counts()
            for level, cnt in counts_by_level.items():
                if level == "Other":
                    continue
                if "Other" not in meta_df["env_biome"].cat.categories:
                    meta_df["env_biome"] = meta_df["env_biome"].cat.add_categories(["Other"])
                meta_df.loc[meta_df["env_biome"] == level, "env_biome"] = "Other"
                meta_df["env_biome"] = meta_df["env_biome"].cat.remove_unused_categories()
                print(f"  Force-merged '{level}' ({cnt} samples) -> 'Other'")
                merged = True
                break

        if not merged:
            print("  WARNING: Could not fix collinearity completely")
            return meta_df

    return meta_df


def prepare_data():
    """Clean EMP500 CSV inputs and cache the prepared DataFrames."""

    counts_csv = input_path("feature_emp500_subset.csv")
    meta_csv = input_path("metadata_emp500_subset.csv")

    counts_df = pd.read_csv(counts_csv, index_col=0)
    meta_df = pd.read_csv(meta_csv, index_col=0)
    meta_df = pd.DataFrame({
        "empo_3": pd.Categorical(meta_df["empo_3"]),
        "env_biome": pd.Categorical(meta_df["env_biome"]),
    }, index=meta_df.index)

    keep = meta_df["empo_3"].notna() & meta_df["env_biome"].notna()
    if keep.sum() < len(keep):
        print(f"  Removing {(~keep).sum()} samples with missing metadata")
        meta_df = meta_df[keep]
        counts_df = counts_df.loc[keep]

    common = counts_df.index.intersection(meta_df.index)
    counts_df = counts_df.loc[common]
    meta_df = meta_df.loc[common]

    nonzero = counts_df.sum(axis=0) > 0
    n_removed = (~nonzero).sum()
    if n_removed > 0:
        print(f"  Removing {n_removed} all-zero taxa")
        counts_df = counts_df.loc[:, nonzero]

    print(f"  Final: {counts_df.shape[0]} samples x {counts_df.shape[1]} taxa")
    print(f"  empo_3: {meta_df['empo_3'].nunique()} levels")
    print(f"  env_biome: {meta_df['env_biome'].nunique()} levels (before collinearity fix)")
    print(f"  Sparsity: {(counts_df == 0).sum().sum() / counts_df.size * 100:.1f}%")

    meta_df = fix_collinearity(meta_df)
    print(f"  env_biome: {meta_df['env_biome'].nunique()} levels (after collinearity fix)")

    atomic_pickle_dump({"counts": counts_df, "meta": meta_df}, PREPARED_DATA)
    print(f"  Saved prepared data: {PREPARED_DATA}")

    return counts_df, meta_df


def load_prepared_data():
    """Load prepared data from pickle."""
    data = load_pickle(PREPARED_DATA)
    return data["counts"], data["meta"]


def bonferroni_combine(pvals):
    """Bonferroni combine: min(p) * n, capped at 1.0."""
    pvals = np.asarray(pvals, dtype=float)
    pvals = pvals[~np.isnan(pvals)]
    if len(pvals) == 0:
        return np.nan
    return min(1.0, pvals.min() * len(pvals))


def bh_adjust(pvals):
    """Benjamini-Hochberg FDR adjustment."""
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    if n == 0:
        return pvals
    order = np.argsort(pvals)
    ranked = pvals[order]
    adjusted = ranked * n / (np.arange(n) + 1)
    for i in range(n - 2, -1, -1):
        adjusted[i] = min(adjusted[i], adjusted[i + 1])
    adjusted = np.clip(adjusted, 0, 1)
    result = np.empty(n)
    result[order] = adjusted
    return result


def run_ancombc(counts_df, meta_df):
    """ANCOM-BC via scikit-bio 0.7.3."""
    from skbio.stats.composition import ancombc

    table = counts_df + 1

    result = ancombc(table=table, metadata=meta_df, formula="empo_3 + env_biome").result

    covariates = result.index.get_level_values("Covariate")
    empo3_mask = covariates.str.startswith("empo_3")
    empo3_res = result[empo3_mask]

    taxa = empo3_res.index.get_level_values("FeatureID").unique()

    combined_p = {}
    for taxon in taxa:
        taxon_pvals = empo3_res.loc[taxon, "pvalue"].values
        combined_p[taxon] = bonferroni_combine(taxon_pvals)

    pvals = np.array([combined_p.get(t, np.nan) for t in taxa])
    qvals = bh_adjust(pvals)

    return pd.DataFrame({
        "taxon": taxa,
        "p_value": pvals,
        "q_value": qvals,
    })


def sensitivity(main, results):
    result = main.copy()
    signif = pd.concat([x for x in results], axis=1)
    result['Pass'] = signif.eq(signif.iloc[:, 0], axis=0).all(axis=1)
    result['Robust'] = result['Signif'] & result['Pass']
    return result


def run_ancombc2(counts_df, meta_df):
    from skbio.stats.composition._ancombc2 import ancombc2

    table = counts_df + 1
    fits = [ancombc2(table, meta_df, formula="empo_3+env_biome", pseudocount=p).res['Signif'] for p in (0.1, 0.5, 1)]
    result = ancombc2(table=table, metadata=meta_df, formula="empo_3 + env_biome").result
    sensitivity(result, fits)

    covariates = result.index.get_level_values("Covariate")
    empo3_mask = covariates.str.startswith("empo_3")
    empo3_res = result[empo3_mask]

    taxa = empo3_res.index.get_level_values("FeatureID").unique()

    combined_p = {}
    for taxon in taxa:
        taxon_pvals = empo3_res.loc[taxon, "pvalue"].values
        combined_p[taxon] = bonferroni_combine(taxon_pvals)

    pvals = np.array([combined_p.get(t, np.nan) for t in taxa])
    qvals = bh_adjust(pvals)

    return pd.DataFrame({
        "taxon": taxa,
        "p_value": pvals,
        "q_value": qvals,
    })


def run_pydeseq2(counts_df, meta_df):
    """PyDESeq2 0.5.4 — Wald tests for each empo_3 level vs reference."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats

    inference = DefaultInference(n_cpus=1)

    dds = DeseqDataSet(
        counts=counts_df,
        metadata=meta_df,
        design="~empo_3 + env_biome",
        inference=inference,
        quiet=True,
    )
    dds.deseq2()

    empo3_levels = sorted(meta_df["empo_3"].unique().tolist())
    ref_level = empo3_levels[0]  # first alphabetically
    test_levels = empo3_levels[1:]

    all_pvals = {}
    for level in test_levels:
        try:
            stat = DeseqStats(
                dds,
                contrast=["empo_3", level, ref_level],
                inference=inference,
                quiet=True,
            )
            stat.summary()
            for taxon in stat.results_df.index:
                if taxon not in all_pvals:
                    all_pvals[taxon] = []
                all_pvals[taxon].append(stat.results_df.loc[taxon, "pvalue"])
        except Exception as e:
            print(f"    PyDESeq2 contrast {level} vs {ref_level} failed: {e}")
            continue

    taxa = list(all_pvals.keys())
    pvals = np.array([bonferroni_combine(all_pvals[t]) for t in taxa])
    qvals = bh_adjust(pvals)

    return pd.DataFrame({
        "taxon": taxa,
        "p_value": pvals,
        "q_value": qvals,
    })


def run_edgepython(counts_df, meta_df):
    """edgePython 0.2.6 — multi-coefficient LRT with fixed dispersion.

    Note: estimate_disp hangs on 97% sparse EMP500 data (same issue as R edgeR's
    estimateDisp with locfit). Using fixed dispersion=0.1, consistent with the
    R edgeR benchmark workaround.
    """
    import edgepython
    from patsy import dmatrix

    dge = edgepython.make_dgelist(
        counts=counts_df.T.values.astype(float),
        group=meta_df["empo_3"].values,
    )
    dge = edgepython.calc_norm_factors(dge)

    design = dmatrix("~empo_3 + env_biome", meta_df)

    dispersion = 0.1
    print(f"    edgePython: using fixed dispersion=0.1 (estimate_disp hangs on sparse data)")

    fit = edgepython.glm_fit(dge, design=design, dispersion=dispersion)

    coef_names = list(design.design_info.column_names)
    empo3_indices = [i for i, name in enumerate(coef_names) if name.startswith("empo_3")]

    if not empo3_indices:
        raise ValueError("No empo_3 coefficients found in design matrix")

    lrt = edgepython.glm_lrt(fit, coef=empo3_indices)
    tt = edgepython.top_tags(lrt, n=np.inf)

    if isinstance(tt, dict):
        table = tt.get("table", tt.get("table", pd.DataFrame()))
    elif isinstance(tt, pd.DataFrame):
        table = tt
    else:
        table = pd.DataFrame(tt)

    pvals = table["PValue"].values
    qvals = table.get("FDR", bh_adjust(pvals)).values

    return pd.DataFrame({
        "taxon": table.index.astype(str),
        "p_value": pvals,
        "q_value": qvals,
    })


def run_pyedger(counts_df, meta_df):
    """pyedger 0.1.0 — joint LRT via manual null-model fit with fixed dispersion.

    pyedger's glmLRT only supports single-coefficient tests (no list of coefs),
    so running 13 individual LRTs would require 13 null-model refits. Instead,
    we perform a single joint LRT by fitting the null model (all empo_3
    coefficients dropped) once and computing LR = null_deviance - full_deviance
    with df = n_empo3_coefs. This is equivalent to edgeR's multi-coefficient
    glmLRT and consistent with the edgepython runner.

    Note: estimateDisp hangs on 97% sparse EMP500 data (same issue as R edgeR's
    estimateDisp with locfit). Using fixed dispersion=0.1, consistent with the
    R edgeR benchmark workaround.
    """
    import pyedger
    from patsy import dmatrix
    from scipy import stats

    dge = pyedger.DGEList(
        counts=counts_df.T.values.astype(float),
        group=meta_df["empo_3"].values,
    )
    dge = pyedger.calcNormFactors(dge)

    design = dmatrix("~empo_3 + env_biome", meta_df)
    design_mat = np.asarray(design, dtype=float)
    coef_names = list(design.design_info.column_names)

    empo3_indices = [i for i, name in enumerate(coef_names) if name.startswith("empo_3")]
    if not empo3_indices:
        raise ValueError("No empo_3 coefficients found in design matrix")

    print(f"    pyedger: using fixed dispersion=0.1, joint LRT on {len(empo3_indices)} empo_3 coefs")

    fit_full = pyedger.glmFit(dge, design=design, dispersion=0.1, prior_count=0)

    null_cols = [i for i in range(design_mat.shape[1]) if i not in empo3_indices]
    design_null = pd.DataFrame(design_mat[:, null_cols],
                               columns=[coef_names[i] for i in null_cols])
    fit_null = pyedger.glmFit(dge, design=design_null, dispersion=0.1, prior_count=0)

    lr = fit_null.deviance - fit_full.deviance
    df_test = len(empo3_indices)
    pvals = stats.chi2.sf(lr, df_test)

    qvals = bh_adjust(pvals)
    taxa = [str(i) for i in range(len(pvals))]

    return pd.DataFrame({
        "taxon": taxa,
        "p_value": pvals,
        "q_value": qvals,
    })


METHOD_RUNNERS = {
    "ancombc": run_ancombc,
    "ancombc2": run_ancombc2,
    "pydeseq2": run_pydeseq2,
    "edgepython": run_edgepython,
    "pyedger": run_pyedger,
}


def run_single_method(method, replicate):
    """Run a single method in the current process and output JSON results."""
    print(f"Running {method} (rep {replicate})...")

    counts_df, meta_df = load_prepared_data()
    runner = METHOD_RUNNERS[method]

    start_time = time.perf_counter()
    status = "SUCCESS"
    error_msg = ""
    n_sig = 0
    n_total = counts_df.shape[1]

    try:
        results = runner(counts_df, meta_df)
        runtime = time.perf_counter() - start_time
        n_sig = int((results["q_value"] < 0.05).sum())

        result_file = output_path(
            f"emp500_{method}_rep{replicate}_results.csv.gz"
        )
        atomic_write_csv(results, result_file, index=False)
        print(f"  Results saved: {result_file} ({n_sig} sig / {n_total} taxa)")
    except Exception as e:
        runtime = time.perf_counter() - start_time
        status = "ERROR"
        error_msg = str(e)[:500]
        print(f"  FAILED: {e}")
        import traceback

        traceback.print_exc()

    peak_mem_raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == "darwin":
        peak_mem_mb = peak_mem_raw / (1024.0 * 1024.0)
    else:
        peak_mem_mb = peak_mem_raw / 1024.0

    output = {
        "method": method,
        "replicate": replicate,
        "runtime_s": round(runtime, 2),
        "peak_memory_MB": round(peak_mem_mb, 1),
        "status": status,
        "error_message": error_msg,
        "n_sig": n_sig,
        "n_total": n_total,
        "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
    }

    print(json.dumps(output))
    return output



def load_existing_results():
    """Load existing results CSV for resume capability."""
    if RESULTS_CSV.exists():
        return pd.read_csv(RESULTS_CSV)
    return pd.DataFrame(columns=[
        "method", "replicate", "runtime_s", "peak_memory_MB",
        "status", "error_message", "n_sig", "n_total", "timestamp"
    ])


def save_results(df):
    atomic_write_csv(df, RESULTS_CSV, index=False)


def run_full_benchmark():
    if PREPARED_DATA.exists():
        print(f"=== Using cached prepared data: {PREPARED_DATA} ===")
    else:
        print("=== Preparing EMP500 Data ===")
        prepare_data()

    existing = load_existing_results()
    completed = set()
    if len(existing) > 0:
        for _, row in existing.iterrows():
            if row["status"] == "SUCCESS":
                completed.add((row["method"], row["replicate"]))
        print(f"\n=== Resuming: {len(completed)} runs already completed ===")

    all_results = []
    if len(existing) > 0:
        all_results.extend(existing.to_dict("records"))

    total = len(METHODS) * N_REPLICATES
    run_idx = 0

    for method in METHODS:
        for rep in range(1, N_REPLICATES + 1):
            run_idx += 1
            if (method, rep) in completed:
                print(f"\n[{run_idx}/{total}] {method} rep{rep} — SKIPPED (already completed)")
                continue

            print(f"\n[{run_idx}/{total}] {method} rep{rep}")
            timeout = TIMEOUTS.get(method, 3600)

            cmd = [
                sys.executable, str(Path(__file__).resolve()),
                "--run-method", method,
                "--replicate", str(rep),
            ]

            try:
                proc = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=timeout,
                )

                output_lines = proc.stdout.strip().split("\n")
                result_json = None
                for line in reversed(output_lines):
                    line = line.strip()
                    if line.startswith("{") and line.endswith("}"):
                        try:
                            result_json = json.loads(line)
                            break
                        except json.JSONDecodeError:
                            continue

                if result_json:
                    all_results.append(result_json)
                    df = pd.DataFrame(all_results)
                    save_results(df)
                    print(f"  -> Status: {result_json['status']}, "
                          f"Runtime: {result_json['runtime_s']}s, "
                          f"Memory: {result_json['peak_memory_MB']} MB, "
                          f"Sig: {result_json['n_sig']}")
                else:
                    error_msg = proc.stderr[:500] if proc.stderr else "Unknown error"
                    result = {
                        "method": method,
                        "replicate": rep,
                        "runtime_s": timeout,
                        "peak_memory_MB": 0,
                        "status": "ERROR",
                        "error_message": error_msg,
                        "n_sig": 0,
                        "n_total": 8351,
                        "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    }
                    all_results.append(result)
                    df = pd.DataFrame(all_results)
                    save_results(df)
                    print(f"  -> FAILED: {error_msg[:200]}")

            except subprocess.TimeoutExpired:
                result = {
                    "method": method,
                    "replicate": rep,
                    "runtime_s": timeout,
                    "peak_memory_MB": 0,
                    "status": "TIMEOUT",
                    "error_message": f"Exceeded {timeout}s timeout",
                    "n_sig": 0,
                    "n_total": 8351,
                    "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                }
                all_results.append(result)
                df = pd.DataFrame(all_results)
                save_results(df)
                print(f"  -> TIMEOUT after {timeout}s")

    print("\n=== Generating Analysis and Plots ===")
    df = pd.DataFrame(all_results)
    generate_analysis(df)



def generate_analysis(df):
    """Generate summary table and plots."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    summary_rows = []
    for method in METHODS:
        mdf = df[df["method"] == method]
        if len(mdf) == 0:
            continue
        success_df = mdf[mdf["status"] == "SUCCESS"]
        status = "SUCCESS" if len(success_df) > 0 else mdf["status"].iloc[0]

        summary_rows.append({
            "method": method,
            "method_label": METHOD_LABELS[method],
            "n_runs": len(mdf),
            "n_success": len(success_df),
            "mean_runtime_s": success_df["runtime_s"].mean() if len(success_df) > 0 else mdf["runtime_s"].mean(),
            "sd_runtime_s": success_df["runtime_s"].std() if len(success_df) > 1 else 0,
            "mean_memory_MB": success_df["peak_memory_MB"].mean() if len(success_df) > 0 else mdf["peak_memory_MB"].mean(),
            "sd_memory_MB": success_df["peak_memory_MB"].std() if len(success_df) > 1 else 0,
            "mean_n_sig": success_df["n_sig"].mean() if len(success_df) > 0 else None,
            "sd_n_sig": success_df["n_sig"].std() if len(success_df) > 1 else 0,
            "status": status,
        })

    summary_df = pd.DataFrame(summary_rows)
    summary_csv = output_path("emp500_python_benchmark_summary.csv")
    atomic_write_csv(summary_df, summary_csv, index=False)
    print(f"Summary saved: {summary_csv}")


    print("\n=== Summary ===")
    print(f"{'Method':<22} {'Runtime':>12} {'Memory':>12} {'Sig taxa':>10} {'Status':>10}")
    print("-" * 70)
    for _, row in summary_df.iterrows():
        rt = row["mean_runtime_s"]
        rt_str = f"{rt:.0f}s" if rt < 60 else f"{rt/60:.1f}min" if rt < 3600 else f"{rt/3600:.1f}h"
        mem = row["mean_memory_MB"]
        mem_str = f"{mem/1024:.1f}GB" if mem > 1024 else f"{mem:.0f}MB"
        sig = f"{row['mean_n_sig']:.0f}" if pd.notna(row["mean_n_sig"]) else "N/A"
        print(f"{row['method_label']:<22} {rt_str:>12} {mem_str:>12} {sig:>10} {row['status']:>10}")

    success_summary = summary_df[summary_df["status"] == "SUCCESS"].sort_values("mean_runtime_s")
    if len(success_summary) > 0:
        fig, ax = plt.subplots(figsize=(10, 5))
        labels = [METHOD_LABELS[m] for m in success_summary["method"]]
        runtimes = success_summary["mean_runtime_s"].values
        errors = summary_df.set_index("method").loc[success_summary["method"], "sd_runtime_s"].fillna(0).values

        bars = ax.barh(range(len(labels)), runtimes, xerr=errors,
                       edgecolor="black", linewidth=0.5, height=0.5)
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=11)
        ax.set_xlabel("Runtime (seconds, log scale)", fontsize=12)
        ax.set_title("EMP500 Python Benchmark: Runtime", fontsize=14)
        ax.set_xscale("log")

        for i, val in enumerate(runtimes):
            label = f"{val:.0f}s" if val < 60 else f"{val/60:.1f}min" if val < 3600 else f"{val/3600:.1f}h"
            ax.text(val * 1.15, i, label, va="center", fontsize=9)

        ax.invert_yaxis()
        plt.tight_layout()
        plt.savefig(output_path("emp500_python_runtime_barplot.png"), dpi=150)
        plt.close()
        print("Runtime barplot saved")


    if len(success_summary) > 0:
        mem_summary = summary_df[summary_df["status"] == "SUCCESS"].sort_values("mean_memory_MB")
        fig, ax = plt.subplots(figsize=(10, 5))
        labels = [METHOD_LABELS[m] for m in mem_summary["method"]]
        memories = mem_summary["mean_memory_MB"].values
        errors = summary_df.set_index("method").loc[mem_summary["method"], "sd_memory_MB"].fillna(0).values

        bars = ax.barh(range(len(labels)), memories, xerr=errors,
                       edgecolor="black", linewidth=0.5, height=0.5)
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=11)
        ax.set_xlabel("Peak Memory (MB, log scale)", fontsize=12)
        ax.set_title("EMP500 Python Benchmark: Peak Memory", fontsize=14)
        ax.set_xscale("log")

        for i, val in enumerate(memories):
            label = f"{val/1024:.1f}GB" if val > 1024 else f"{val:.0f}MB"
            ax.text(val * 1.15, i, label, va="center", fontsize=9)

        ax.invert_yaxis()
        plt.tight_layout()
        plt.savefig(output_path("emp500_python_memory_barplot.png"), dpi=150)
        plt.close()
        print("Memory barplot saved")

    if len(success_summary) > 0:
        fig, ax = plt.subplots(figsize=(8, 6))
        for _, row in success_summary.iterrows():
            m = row["method"]
            ax.scatter(row["mean_runtime_s"], row["mean_memory_MB"], s=150, 
                       )
            ax.annotate(METHOD_LABELS[m], (row["mean_runtime_s"], row["mean_memory_MB"]),
                        xytext=(row["mean_runtime_s"], row["mean_memory_MB"]),
                        fontsize=10, ha="left", va="bottom")

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Runtime (seconds, log scale)", fontsize=12)
        ax.set_ylabel("Peak Memory (MB, log scale)", fontsize=12)
        ax.set_title("EMP500 Python Benchmark: Runtime vs Memory", fontsize=14)
        ax.grid(True, alpha=0.3, which="both")
        plt.tight_layout()
        plt.savefig(output_path("emp500_python_runtime_memory_scatter.png"), dpi=150)
        plt.close()
        print("Scatter plot saved")


    sig_summary = summary_df[summary_df["mean_n_sig"].notna()].sort_values("mean_n_sig", ascending=False)
    if len(sig_summary) > 0:
        fig, ax = plt.subplots(figsize=(10, 5))
        labels = [METHOD_LABELS[m] for m in sig_summary["method"]]
        sig_vals = sig_summary["mean_n_sig"].values

        bars = ax.barh(range(len(labels)), sig_vals,
                       edgecolor="black", linewidth=0.5, height=0.5)
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=11)
        ax.set_xlabel("Number of Significant Taxa (q < 0.05)", fontsize=12)
        ax.set_title("EMP500 Python Benchmark: Significant Taxa", fontsize=14)

        for i, val in enumerate(sig_vals):
            ax.text(val + 50, i, f"{int(val)}", va="center", fontsize=9)

        ax.invert_yaxis()
        plt.tight_layout()
        plt.savefig(output_path("emp500_python_sig_taxa_barplot.png"), dpi=150)
        plt.close()
        print("Sig taxa barplot saved")


def generate_r_vs_python_comparison(py_df, r_csv):
    """Generate R vs Python comparison plot."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt


    r_df = pd.read_csv(r_csv)
    r_df = r_df[r_df["status"] == "SUCCESS"]

    r_to_py = {
        "ancombc": "ancombc",
        "deseq2": "pydeseq2",
        "edger": "edgepython",
    }

    pairs = []
    for r_method, py_method in r_to_py.items():
        r_data = r_df[r_df["method"] == r_method]
        py_data = py_df[(py_df["method"] == py_method) & (py_df["status"] == "SUCCESS")]
        if len(r_data) > 0 and len(py_data) > 0:
            pairs.append({
                "r_method": r_method,
                "py_method": py_method,
                "r_runtime": r_data["runtime_s"].mean(),
                "py_runtime": py_data["runtime_s"].mean(),
                "r_memory": r_data["peak_memory_MB"].mean(),
                "py_memory": py_data["peak_memory_MB"].mean(),
                "r_sig": r_data["n_sig"].mean() if "n_sig" in r_data.columns else None,
                "py_sig": py_data["n_sig"].mean(),
            })

    if not pairs:
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Runtime comparison
    ax1 = axes[0]
    x = np.arange(len(pairs))
    width = 0.35
    r_runtimes = [p["r_runtime"] for p in pairs]
    py_runtimes = [p["py_runtime"] for p in pairs]
    ax1.bar(x - width/2, r_runtimes, width, label="R", color="#0279EE", edgecolor="black", linewidth=0.5)
    ax1.bar(x + width/2, py_runtimes, width, label="Python", color="#E9ED4C", edgecolor="black", linewidth=0.5)
    ax1.set_ylabel("Runtime (seconds, log scale)", fontsize=11)
    ax1.set_title("Runtime: R vs Python", fontsize=13)
    ax1.set_yscale("log")
    ax1.set_xticks(x)
    ax1.set_xticklabels([f"{p['r_method']} →\n{p['py_method']}" for p in pairs], fontsize=9)
    ax1.legend(fontsize=10)

    for i, (r_rt, py_rt) in enumerate(zip(r_runtimes, py_runtimes)):
        speedup = r_rt / py_rt if py_rt > 0 else 0
        ax1.text(i, max(r_rt, py_rt) * 1.5, f"{speedup:.1f}x",
                ha="center", fontsize=10, fontweight="bold", color="#75A025")

    # Memory comparison
    ax2 = axes[1]
    r_memories = [p["r_memory"] for p in pairs]
    py_memories = [p["py_memory"] for p in pairs]

    ax2.bar(x - width/2, r_memories, width, label="R", color="#0279EE", edgecolor="black", linewidth=0.5)
    ax2.bar(x + width/2, py_memories, width, label="Python", color="#E9ED4C", edgecolor="black", linewidth=0.5)
    ax2.set_ylabel("Peak Memory (MB, log scale)", fontsize=11)
    ax2.set_title("Memory: R vs Python", fontsize=13)
    ax2.set_yscale("log")
    ax2.set_xticks(x)
    ax2.set_xticklabels([f"{p['r_method']} →\n{p['py_method']}" for p in pairs], fontsize=9)
    ax2.legend(fontsize=10)

    plt.suptitle("EMP500: R vs Python DA Method Comparison", fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(
        output_path("emp500_python_r_vs_comparison.png"),
        dpi=150,
        bbox_inches="tight",
    )
    plt.close()
    print("R vs Python comparison plot saved")


def main():
    parser = argparse.ArgumentParser(description="EMP500 Python DA Method Benchmark")
    parser.add_argument("--run-method", type=str, help="Run a single method")
    parser.add_argument("--replicate", type=int, default=1, help="Replicate number")
    parser.add_argument("--analyze", action="store_true", help="Only generate plots from existing results")
    args = parser.parse_args()

    if args.run_method:
        run_single_method(args.run_method, args.replicate)
    elif args.analyze:
        df = load_existing_results()
        if len(df) == 0:
            print("No results found. Run the benchmark first.")
            return
        generate_analysis(df)
    else:
        run_full_benchmark()


if __name__ == "__main__":
    main()
