import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1" 
os.environ["OPENBLAS_NUM_THREADS"] = "1" 
os.environ["VECLIB_MAXIMUM_THREADS"] = "1" 
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import pandas as pd
from skbio.stats.composition import ancombc
import time
import tracemalloc
from threadpoolctl import threadpool_limits
from benchmark_io import atomic_write_csv, input_path, output_path

# setup
COUNTS_FILE = "feature_sc.csv"
META_FILE = "metadata_sc.csv"
OUTPUT = output_path("benchmark_sc_py.csv")
N_REPEATS = 3



def run_benchmark_iteration(table, metadata, formula):
    tracemalloc.start()

    start_time = time.perf_counter()
    success = False
    error_msg = None
    res = None
    
    try:
        with threadpool_limits(limits=1):
            res = ancombc(
                table + 1, 
                metadata, 
                formula=formula, 
                p_adjust='holm', 
                alpha=0.05
            )
        success = True
    except Exception as e:
        error_msg = str(e)
    
    end_time = time.perf_counter()
    _, peak_mem = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    elapsed_sec = end_time - start_time
    peak_mem_mb = peak_mem / (1024 * 1024)

    return success, elapsed_sec, peak_mem_mb, error_msg, res


def main():
    print("loading data")
    counts_df = pd.read_csv(input_path(COUNTS_FILE), index_col=0)
    meta_df = pd.read_csv(input_path(META_FILE), index_col=0)

    print("\n=================================================")
    formula_str = "cond"
    results_log = []

    print("   > Warm-up... ", end="", flush=True)
    warmup_success, _, _, warmup_error, _ = run_benchmark_iteration(
        counts_df, meta_df, formula_str
    )
    if warmup_success:
        print("OK")
    else:
        print(f"FAILED. Error: {warmup_error}")

    for i in range(1, N_REPEATS + 1):
        print(f"   > Run {i}/{N_REPEATS}... ", end="", flush=True)
        success, duration, mem_mb, error, _ = run_benchmark_iteration(
            counts_df, meta_df, formula_str
        )
        status = "Success" if success else "Failed"
        if success:
            print(f"{duration:.2f}s | Peak RAM: {mem_mb:.1f} MB")
        else:
            print(f"FAILED. Error: {error}")

        results_log.append({
            "Dataset": "single cell",
            "Run_ID": i,
            "Time_Sec": duration,
            "Memory_Peak_MB": mem_mb,
            "Status": status,
            "Error": error or "",
            "Timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        })
        atomic_write_csv(pd.DataFrame(results_log), OUTPUT, index=False)

    print("\nbenchmark completed")
    results_df = pd.DataFrame(results_log)
    success_df = results_df[results_df["Status"] == "Success"]
    if not success_df.empty:
        summary = success_df.groupby(["Dataset"]).agg(
            Avg_Time_Sec=("Time_Sec", "mean"),
            Avg_Mem_MB=("Memory_Peak_MB", "mean"),
            N_Runs=("Status", "count"),
        ).reset_index()
        print("\nSummary Table:")
        print(summary.to_string())
    print(f"\nFull metrics saved to: {OUTPUT}")


if __name__ == "__main__":
    main()
