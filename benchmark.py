#!/usr/bin/env python3
"""Benchmark Python ancombc and ancombc2.

vary samples, vary features, vary covariates, full 10k×10k
2 implementations: ancombc, ancombc2
3 repeats per scenario (after 1 warm-up)
Memory: tracemalloc (paper method) + psutil RSS (OS-level)
"""

import os
# Single-thread enforcement
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import sys
import time
import json
import threading
import tracemalloc
import importlib.util
import numpy as np
import pandas as pd
import psutil


from skbio.stats.composition._ancombc2 import ancombc, ancombc2

COUNTS_FILE = "sim_counts_10k_20_meta.csv.gz"
META_FILE = "sim_metadata_10k_20_meta.csv"
OUTPUT = "bench_py_results_pseudo_sens.csv"
N_REPEATS = 3
WARMUP = 1


scenarios_vary_n = pd.DataFrame({
    "N": [10, 100, 1000, 2000, 4000, 6000, 8000, 10000],
    "P": [1000] * 8,
    "Regime": ["Vary_Samples"] * 8,
    "N_Cov": [2] * 8,
})

scenarios_vary_p = pd.DataFrame({
    "N": [1000] * 8,
    "P": [10, 100, 1000, 2000, 4000, 6000, 8000, 10000],
    "Regime": ["Vary_Features"] * 8,
    "N_Cov": [2] * 8,
})

scenarios_vary_cov = pd.DataFrame({
    "N": [1000] * 5,
    "P": [1000] * 5,
    "Regime": ["Vary_Covariates"] * 5,
    "N_Cov": [2, 4, 6, 8, 10],
})

scenarios_full = pd.DataFrame({
    "N": [10000],
    "P": [10000],
    "Regime": ["Full_Matrix"],
    "N_Cov": [10],
})

scenarios = pd.concat([
    scenarios_vary_n,
    scenarios_vary_p,
    scenarios_vary_cov,
    scenarios_full,
]).drop_duplicates(subset=["N", "P", "N_Cov"]).reset_index(drop=True)


def build_formula(n_cov):
    """Build formula string for n_cov covariates (2, 4, 6, 8, 10)."""
    parts = []
    for j in range(1, n_cov // 2 + 1):
        parts.append(f"cont_cov_{j}")
        parts.append(f"cat_cov_{j}")
    return " + ".join(parts)


def sensitivity(main, results):
    result = main.copy()
    signif = pd.concat([x for x in results], axis=1)
    result['Pass'] = signif.eq(signif.iloc[:, 0], axis=0).all(axis=1)
    result['Robust'] = result['Signif'] & result['Pass']
    return result


class RSSSampler:
    """Sample process RSS in a background thread."""

    def __init__(self, interval_ms=50):
        self.interval = interval_ms / 1000.0
        self._thread = None
        self._stop_event = threading.Event()
        self._samples = []
        self._baseline = 0

    def start(self):
        proc = psutil.Process()
        self._baseline = proc.memory_info().rss
        self._samples = [self._baseline]
        self._stop_event.clear()
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()

    def _run(self):
        proc = psutil.Process()
        while not self._stop_event.is_set():
            try:
                self._samples.append(proc.memory_info().rss)
            except Exception:
                pass
            self._stop_event.wait(self.interval)

    def stop(self):
        self._stop_event.set()
        if self._thread:
            self._thread.join(timeout=2)
        proc = psutil.Process()
        try:
            self._samples.append(proc.memory_info().rss)
        except Exception:
            pass
        peak = max(self._samples) if self._samples else self._baseline
        return {
            "baseline_mb": self._baseline / (1024 * 1024),
            "peak_mb": peak / (1024 * 1024),
            "delta_mb": (peak - self._baseline) / (1024 * 1024),
        }



from threadpoolctl import threadpool_limits


def run_one(table, metadata, formula, impl):
    """Run a single benchmark measurement."""
    tracemalloc.start()
    rss_sampler = RSSSampler(interval_ms=50)
    rss_sampler.start()

    start_time = time.time()
    success = False
    error_msg = None

    try:
        with threadpool_limits(limits=1):
            if impl == "ancombc":
                res = ancombc(
                    table + 1,
                    metadata,
                    formula=formula,
                    p_adjust="holm",
                    alpha=0.05,
                ).result
            elif impl == "ancombc2":
                fits = [ancombc2(table+1, metadata, formula=formula, pseudocount=p).res['Signif'] for p in (0.1, 0.5, 1)]
                res = ancombc2(
                    table+1,
                    metadata,
                    formula=formula,
                    p_adjust="holm",
                    alpha=0.05,
                    pseudocount=1,
                ).result
                res_sens = sensitivity(res, fits)
            else:
                raise ValueError(f"Unknown impl: {impl}")
        success = True
    except Exception as e:
        error_msg = str(e)[:200]

    end_time = time.time()
    rss_info = rss_sampler.stop()
    _, peak_mem = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    elapsed = end_time - start_time
    tracemalloc_mb = peak_mem / (1024 * 1024)

    return {
        "result": res,
        "success": success,
        "time_sec": elapsed,
        "tracemalloc_mb": tracemalloc_mb,
        "rss_baseline_mb": rss_info["baseline_mb"],
        "rss_peak_mb": rss_info["peak_mb"],
        "rss_delta_mb": rss_info["delta_mb"],
        "error": error_msg,
    }


def main():
    print("Loading data...")
    counts_df = pd.read_csv(COUNTS_FILE, index_col=0)
    meta_df = pd.read_csv(META_FILE)
    meta_df.index = "S" + (meta_df.index + 1).astype(str)
    meta_df = meta_df.loc[counts_df.index]
    print(f"  Data loaded: {counts_df.shape[0]} samples × {counts_df.shape[1]} features")

    results_log = []
    # Load existing results if checkpoint exists (for resume)
    if os.path.exists(OUTPUT):
        existing = pd.read_csv(OUTPUT)
        results_log = existing.to_dict("records")
        done_keys = set(
            (r["Regime"], r["Samples"], r["Features"], r["Impl"], r["Run_ID"])
            for r in results_log
            if r.get("Status") == "Success"
        )
        print(f"  Resuming: {len(results_log)} existing records, {len(done_keys)} completed runs")
    else:
        done_keys = set()

    total_scenarios = len(scenarios)
    for idx, row in scenarios.iterrows():
        curr_n = int(row["N"])
        curr_p = int(row["P"])
        regime = row["Regime"]
        n_cov = int(row["N_Cov"])
        formula = build_formula(n_cov)

        print(f"\n{'='*60}")
        print(f"[{idx+1}/{total_scenarios}] {regime} | {curr_n}×{curr_p} | cov={n_cov}")
        print(f"  Formula: {formula}")

        # Subset data
        sub_counts = counts_df.iloc[:curr_n, :curr_p]
        sub_meta = meta_df.iloc[:curr_n]

        for impl in ["ancombc", "ancombc2"]:
            print(f"\n  --- {impl} ---")

            # Warm-up
            print(f"  warm-up... ", end="", flush=True)
            w = run_one(sub_counts, sub_meta, formula, impl)
            if w["success"]:
                w["result"].to_csv(f"~/proj/ancombc/ancom-bc_bm-main/benchmark_sim_res/impl_{impl}_s_{curr_n}_f_{curr_p}_cov_{n_cov}.csv")
                print(f"OK ({w['time_sec']:.2f}s)")
            else:
                print(f"FAILED: {w['error']}")
                continue

            # Measured runs
            for i in range(1, N_REPEATS + 1):
                key = (regime, curr_n, curr_p, impl, i)
                if key in done_keys:
                    print(f"  run {i}/{N_REPEATS} (cached)")
                    continue

                print(f"  run {i}/{N_REPEATS}... ", end="", flush=True)
                r = run_one(sub_counts, sub_meta, formula, impl)

                status = "Success" if r["success"] else "Failed"
                if r["success"]:
                    print(
                        f"{r['time_sec']:.2f}s | "
                        f"tracemalloc: {r['tracemalloc_mb']:.1f}MB | "
                        f"RSS delta: {r['rss_delta_mb']:.1f}MB"
                    )
                else:
                    print(f"FAILED: {r['error']}")

                entry = {
                    "Regime": regime,
                    "Samples": curr_n,
                    "Features": curr_p,
                    "N_Cov": n_cov,
                    "Impl": f"Py_{impl}",
                    "Run_ID": i,
                    "Time_Sec": r["time_sec"],
                    "Mem_Trace_MB": r["tracemalloc_mb"],
                    "RSS_Baseline_MB": r["rss_baseline_mb"],
                    "RSS_Peak_MB": r["rss_peak_mb"],
                    "RSS_Delta_MB": r["rss_delta_mb"],
                    "Status": status,
                    "Error": r["error"] if not r["success"] else "",
                    "Timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
                }
                results_log.append(entry)

                # Checkpoint after each run
                pd.DataFrame(results_log).to_csv(OUTPUT, index=False)

        # Force garbage collection between scenarios
        import gc
        gc.collect()

    # Final save
    pd.DataFrame(results_log).to_csv(OUTPUT, index=False)
    print(f"\n{'='*60}")
    print(f"Benchmark complete. Results saved to {OUTPUT}")

    # Summary
    df = pd.DataFrame(results_log)
    if not df.empty:
        success_df = df[df["Status"] == "Success"]
        summary = success_df.groupby(["Regime", "Samples", "Features", "Impl"]).agg(
            Avg_Time_Sec=("Time_Sec", "mean"),
            Avg_Mem_Trace_MB=("Mem_Trace_MB", "mean"),
            Avg_RSS_Delta_MB=("RSS_Delta_MB", "mean"),
            N_Runs=("Status", "count"),
        ).reset_index()
        print("\nSummary:")
        print(summary.to_string())


if __name__ == "__main__":
    main()
