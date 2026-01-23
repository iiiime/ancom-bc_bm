import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1" 
os.environ["OPENBLAS_NUM_THREADS"] = "1" 
os.environ["VECLIB_MAXIMUM_THREADS"] = "1" 
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import pandas as pd
import numpy as np
from skbio.stats.composition import ancombc
import time
import tracemalloc
from threadpoolctl import threadpool_limits

# setup
COUNTS_FILE = "sim_counts_10k_20_meta.csv"
META_FILE = "sim_metadata_10k_20_meta.csv"
OUTPUT = "benchmark_metrics_python.csv"
# MAX_CPU = os.cpu_count()
N_REPEATS = 3

scenarios_vary_n = pd.DataFrame({
    'N': [10, 100, 1000, 2000, 4000, 6000, 8000, 10000],
    'P': [1000] * 8,
    'Regime': ['Vary_Samples'] * 8
})

scenarios_vary_p = pd.DataFrame({
    'N': [1000] * 8,
    'P': [10, 100, 1000, 2000, 4000, 6000, 8000, 10000],
    'Regime': ['Vary_Features'] * 8
})

scenarios_full_mat = pd.DataFrame({
    'N': [10000],
    'P': [10000],
    'Regime': ['Full_matrix']
    })
scenarios = pd.concat([scenarios_vary_n, scenarios_vary_p, scenarios_full_mat]).drop_duplicates().reset_index(drop=True)


print("loading data")
counts_df = pd.read_csv(COUNTS_FILE, index_col=0)
meta_df = pd.read_csv(META_FILE)
meta_df.index = 'S' + (meta_df.index + 1).astype(str)
meta_df = meta_df.loc[counts_df.index]


def run_benchmark(table, metadata, formula):
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
    _, peak_mem = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    elapsed_sec = end_time - start_time
    peak_mem_mb = peak_mem / (1024 * 1024) # convert to MB
    
    return success, elapsed_sec, peak_mem_mb, error_msg


print(f"\nstarting benchmarking")
results_log = []

for idx, row in scenarios.iterrows():
    curr_n = int(row['N'])
    curr_p = int(row['P'])
    regime = row['Regime']

    print(f"\n=================================================")
    print(f"{regime} | {curr_n} samples x {curr_p} features")
    
    # subset data
    sub_counts = counts_df.iloc[:curr_n, :curr_p]
    sub_meta = meta_df.iloc[:curr_n]
    
    # warm-up run
    print("warm-up run........... ")
    w_success, _, _, w_err = run_benchmark(sub_counts, sub_meta, "cont_cov_1 + cat_cov_1")
    if w_success:
        print("OK.")
    else:
        print(f"FAILED ({w_err})")
        continue

    # measured runs
    for i in range(1, N_REPEATS + 1):
        print(f"measured run {i}/{N_REPEATS}... ")
        
        success, duration, mem_mb, err = run_benchmark(sub_counts, sub_meta, "cont_cov_1 + cat_cov_1")
        
        status = "Success" if success else "Failed"
        
        if success:
            print(f"{duration:.2f}s | Peak RAM: {mem_mb:.1f} MB")
        else:
            print(f"FAILED. Error: {err}")
            
        # intermediate results
        results_log.append({
            'Regime': regime,
            'Samples': curr_n,
            'Features': curr_p,
            'Run_ID': i,
            'Time_Sec': duration,
            'Memory_Peak_MB': mem_mb,
            'Status': status,
            'Timestamp': time.strftime("%Y-%m-%d %H:%M:%S")
        })
        

pd.DataFrame(results_log).to_csv(OUTPUT, index=False)
print("\nbenchmark dompleted")

results_df = pd.DataFrame(results_log)
if not results_df.empty:
    summary = results_df[results_df['Status'] == 'Success'].groupby(['Regime', 'Samples', 'Features']).agg(
        Avg_Time_Sec=('Time_Sec', 'mean'),
        Avg_Mem_MB=('Memory_Peak_MB', 'mean'),
        Success_Rate=('Status', lambda x: (x == 'Success').sum() / N_REPEATS)
    ).reset_index()
    
    print("\nsummary table:")
    print(summary.to_string())
    print(f"\nFull metrics saved to: {OUTPUT}")
