library(ANCOMBC)
library(TreeSummarizedExperiment)
library(tidyverse)


feature <- "sim_counts_10k_20_meta.csv"
meta <- "sim_metadata_10k_20_meta.csv"
OUTPUT <- "benchmark_metrics_ancombc.csv"
N_REPEATS <- 3
CORES     <- 1

scenarios_vary_n <- data.frame(
  N = c(10, 100, 1000, 2000, 4000, 6000, 8000, 10000),
  P = 1000,
  Regime = "Vary_Samples"
)

scenarios_vary_p <- data.frame(
  N = 1000,
  P = c(10, 100, 1000, 2000, 4000, 6000, 8000, 10000), 
  Regime = "Vary_Features"
)

scenarios_full_mat <- data.frame(N=10000, P=10000, Regime="full_matrix")

scenarios <- rbind(scenarios_vary_n, scenarios_vary_p, scenarios_full_mat)
scenarios <- scenarios[!duplicated(scenarios[, c("N", "P")]), ] # remove ducplicated 1000 * 1000

cat("loading dataset\n")
if(!file.exists(feature)) stop("file not exist")
counts_df <- read.csv(feature, row.names = 1)
counts_df <- as.matrix(counts_df)

meta_df <- read.csv(meta)
rownames(meta_df) = paste0("S", seq_len(10000))
master_meta <- meta_df

counts_df <- t(counts_df) 

total_p <- nrow(counts_df)
total_n <- ncol(counts_df)

cat(sprintf("   - Master Data Loaded: %d Features x %d Samples\n", total_p, total_n))



run_ancombc <- function(tse) {
  ancombc(
    data = tse,
    assay_name = "counts",
    formula = "cont_cov_1 + cat_cov_1",
    p_adj_method = "holm",
    prv_cut = 0.0,
    lib_cut = 0,    
    group = "cat_cov_1",
    struc_zero = FALSE,
    neg_lb = FALSE,
    alpha = 0.05,
    n_cl = CORES,
    verbose = FALSE
  )
}

cat("\nbenchmark number of samples and features\n")
results_log <- data.frame()

for (idx in 1:nrow(scenarios)) {
  
  curr_n <- scenarios$N[idx]
  curr_p <- scenarios$P[idx]
  regime <- scenarios$Regime[idx]
  
  cat(sprintf("\n=================================================\n"))
  cat(sprintf("%d Samples * %d Features\n", regime, curr_n, curr_p))
  
  sub_counts <- counts_df[1:curr_p, 1:curr_n, drop=FALSE]
  sub_meta   <- master_meta[1:curr_n, , drop=FALSE]
  
  tse_sub <- TreeSummarizedExperiment(
    assays = list(counts = sub_counts),
    colData = sub_meta
  )
  
  # warm-up runs
  cat("Warm-up run......... ")
  tryCatch({
    invisible(
run_ancombc(tse_sub))
    cat("OK.\n")
  }, error = function(e) {
    cat("FAILED.\n")
    cat("     Error:", e$message, "\n")
  })
  
  # measured runs
  for (i in 1:N_REPEATS) {
    cat(sprintf("Run %d/%d... ", i, N_REPEATS))
    
    gc(reset = TRUE)
    initial_mem <- sum(gc()[, 2])
    
    start_time <- proc.time()
    
    status <- "Success"
    error_msg <- NA
    out <- NULL
    
    tryCatch({
      out <- 
run_ancombc(tse_sub)
    }, error = function(e) {
      status <<- "Failed"
      error_msg <<- e$message
    })
    
    end_time <- proc.time()
    elapsed_sec <- (end_time - start_time)["elapsed"]
    
    final_gc <- gc()
    peak_mem_mb <- sum(final_gc[, 6])
    mem_diff_mb <- peak_mem_mb - initial_mem
    
    if (status == "Success") {
      cat(sprintf("%.2fs | RAM: %.1f MB\n", elapsed_sec, mem_diff_mb))
    } else {
      cat(sprintf("FAILED (%s)\n", error_msg))
    }
    
    # Log Entry
    entry <- data.frame(
      Regime = regime,
      Samples = curr_n,
      Features = curr_p,
      Run_ID = i,
      Time_Sec = as.numeric(elapsed_sec),
      Memory_Used_MB = mem_diff_mb,
      Status = status,
      Timestamp = as.character(Sys.time())
    )
    
    results_log <- rbind(results_log, entry)
    
  }
}

write.csv(results_log, OUTPUT, row.names = FALSE)
cat("\n benchmark completed\n")

if (nrow(results_log) > 0) {
  summary_table <- results_log %>%
    filter(Status == "Success") %>%
    group_by(Regime, Samples, Features) %>%
    summarise(
      Avg_Time_Sec = mean(Time_Sec),
      Avg_Mem_MB = mean(Memory_Used_MB),
      Success_Rate = n() / N_REPEATS
    ) %>%
    arrange(Regime, Samples, Features)
  
  print(summary_table)
  cat(sprintf("\nresults saved to: %s\n", OUTPUT))
}