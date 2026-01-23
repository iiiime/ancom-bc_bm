import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1" 
os.environ["OPENBLAS_NUM_THREADS"] = "1" 
os.environ["VECLIB_MAXIMUM_THREADS"] = "1" 
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import numpy as np
import pandas as pd
from skbio.stats.composition import ancombc
import time
import tracemalloc
import sys
from threadpoolctl import threadpool_limits

# setup
COUNTS_FILE = "feature_sc.csv"
META_FILE = "metadata_sc.csv"
OUTPUT = "benchmark_sc_py.csv"
#MAX_CPU = os.cpu_count()
N_REPEATS = 3



print("loading data")
if not os.path.exists(COUNTS_FILE):
    sys.exit(f"Error: {COUNTS_FILE} not found.")


counts_df = pd.read_csv(COUNTS_FILE, index_col=0)
meta_df = pd.read_csv(META_FILE, index_col=0)
#meta_df.index = 'S' + (meta_df.index + 1).astype(str)
#meta_df = meta_df.loc[counts_df.index]


def run_benchmark_iteration(table, metadata, formula):
    tracemalloc.start()
    
    start_time = time.time()
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
    
    end_time = time.time()
    _, peak_mem = tracemalloc.get_traced_memory() # TODO: use memory profiler
    # use whichtime to track the results entire script
    # $(which time) python benchmark_cores.py 
    # OMP_NUM_THREADS=4 python script.py for simple
    tracemalloc.stop()
    elapsed_sec = end_time - start_time # or use timeit
    peak_mem_mb = peak_mem / (1024 * 1024) # convert to MB
    
    return success, elapsed_sec, peak_mem_mb, error_msg, res


print(f"\n benchmark number of cores")
results_log = []


print(f"\n=================================================")
formula_str = "cond"


# warm-up run
print("   > Warm-up... ", end="", flush=True)
w_success, _, _, w_err, _ = run_benchmark_iteration(counts_df, meta_df, formula_str)

# measured run
for i in range(1, N_REPEATS + 1):
    print(f"   > Run {i}/{N_REPEATS}... ", end="", flush=True)
    
    success, duration, mem_mb, err, out_res = run_benchmark_iteration(counts_df, meta_df, formula_str)
    
    status = "Success" if success else "Failed"
    
    if success:
        print(f"{duration:.2f}s | Peak RAM: {mem_mb:.1f} MB")
    else:
        print(f"FAILED. Error: {err}")
        
    results_log.append({
        'Dataset': "single cell",
        'Run_ID': i,
        'Time_Sec': duration,
        'Memory_Peak_MB': mem_mb,
        'Status': status,
        'Timestamp': time.strftime("%Y-%m-%d %H:%M:%S")
    })
        
pd.DataFrame(results_log).to_csv(OUTPUT, index=False)
print("\nbenchmark completed")

results_df = pd.DataFrame(results_log)
if not results_df.empty:
    summary = results_df[results_df['Status'] == 'Success'].groupby(['Dataset']).agg(
        Avg_Time_Sec=('Time_Sec', 'mean'),
        Avg_Mem_MB=('Memory_Peak_MB', 'mean'),
        Success_Rate=('Status', lambda x: (x == 'Success').sum() / N_REPEATS)
    ).reset_index()
    
    print("\nSummary Table:")
    print(summary.to_string())
    print(f"\nFull metrics saved to: {OUTPUT}")
