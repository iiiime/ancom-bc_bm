library(ANCOMBC)
library(TreeSummarizedExperiment)
library(tidyverse)
library(parallel)

# setup
abn <- "sim_counts_10k.csv"
meta <- "sim_metadata_10k.csv"
OUT <- "benchmark_n_cores.csv"
N_REPEATS <- 3
N_FEATURES <- 1000
N_SAMPLES <- 1000

cat(sprintf("\n=================================================\n"))
cat("loading dataset...\n")
if(!file.exists(abn)) stop("Count file not found.")
counts_df <- read.csv(abn, row.names = 1) # sample * feature
counts_df <- as.matrix(counts_df)
counts_df <- t(counts_df) 

meta_df <- read.csv(meta)
rownames(meta_df) = paste0("S", seq_len(10000))

# subset data
tse <- TreeSummarizedExperiment(
  assays = list(counts = counts_df[1:N_FEATURES, 1:N_SAMPLES, drop=FALSE]),
  colData = meta_df[1:N_SAMPLES, , drop=FALSE]
)

run_ancombc <- function(formula) {
  ancombc(
    data = tse,
    formula = formula,
    p_adj_method = "holm",
    prv_cut = 0.0,
    lib_cut = 0,    
    group = "cat_cov_1",
    struc_zero = FALSE,
    neg_lb = FALSE,
    alpha = 0.05,
    n_cl = n_cores,
    verbose = FALSE
  )
}

results_log <- data.frame()

# test on different number of cores
for (idx in 1:5) {
  indices <- 1:idx
  formula = paste0("cont_cov_", indices, " + cat_cov_", indices)
  formula = paste('~', paste(formula, collapse = " + "))
  
  cat(sprintf("number of cov: %d \n ", idx*2))
  # warm-up run
  cat("warm-up run \n")
  tryCatch({
    invisible(run_ancombc(formula))
  }, error = function(e) {
    cat("FAILED.\n")
    cat("     Error:", e$message, "\n")
  })
  
  # measured runs
  for (i in 1:N_REPEATS) {
    cat(sprintf("Run %d/%d", i, N_REPEATS))
    
    gc(reset = TRUE)
    initial_mem <- sum(gc()[, 2])
    start_time <- proc.time()
    
    status <- "Success"
    error_msg <- NA
    out <- NULL
    
    tryCatch({
      out <- run_ancombc(formula)
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
      Cores = idx,
      Run_ID = i,
      Time_Sec = as.numeric(elapsed_sec),
      Memory_Used_MB = mem_diff_mb,
      Status = status,
      Timestamp = as.character(Sys.time())
    )
    
    results_log <- rbind(results_log, entry)
    write.csv(results_log, OUT, row.names = FALSE)
  }
}

cat("\n benchmark completed \n")

if (nrow(results_log) > 0) {
  summary_table <- results_log %>%
    filter(Status == "Success") %>%
    group_by(Cores) %>%
    summarise(
      Avg_Time_Sec = mean(Time_Sec),
      Avg_Mem_MB = mean(Memory_Used_MB)
    ) %>%
    arrange(Cores)
  
  print(summary_table)
  cat(sprintf("\n results saved to: %s\n", OUT))
}