#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = FALSE))
} else {
  normalizePath(".", mustWork = FALSE)
}
source(file.path(script_dir, "benchmark_io.R"))

library(ANCOMBC)
library(tidyverse)

# setup
abn <- input_path("feature_sc.csv")
meta <- input_path("metadata_sc.csv")
OUT <- output_path("benchmark_real_sc.csv")
N_REPEATS <- 3

cat(sprintf("\n=================================================\n"))
cat("loading dataset...\n")

counts_df <- read.csv(abn, row.names = 1) # sample * feature
#counts_df <- as.matrix(counts_df)
counts_df <- t(counts_df) 

meta_df <- read.csv(meta, row.names = 1)
formula <- "cond"

run_ancombc <- function(counts_df, meta_df, formula) {
  ancombc(
    data = counts_df,
    meta_data = meta_df,
    formula = formula,
    p_adj_method = "holm",
    prv_cut = 0.0,
    lib_cut = 0,
    struc_zero = FALSE,
    neg_lb = TRUE,
    tol = 1e-5,
    max_iter = 100,
    conserve = FALSE,
    alpha = 0.05,
    n_cl = 1,
    global = FALSE,
    verbose = FALSE
  )
}

results_log <- data.frame()

# test on real data

# warm-up run
cat("warm-up run \n")
tryCatch({
  invisible(run_ancombc(counts_df, meta_df, formula))
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
    out <- run_ancombc(counts_df, meta_df, formula)
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
    Dataset = abn,
    Run_ID = i,
    Time_Sec = as.numeric(elapsed_sec),
    Memory_Used_MB = mem_diff_mb,
    Status = status,
    Error = if (status == "Failed") as.character(error_msg) else "",
    Timestamp = as.character(Sys.time())
  )
  
  results_log <- rbind(results_log, entry)
  atomic_write_csv(results_log, OUT)
}


cat("\n benchmark completed \n")

if (nrow(results_log) > 0) {
  summary_table <- results_log %>%
    filter(Status == "Success") %>%
    group_by(Dataset) %>%
    summarise(
      Avg_Time_Sec = mean(Time_Sec),
      Avg_Mem_MB = mean(Memory_Used_MB)
    ) %>%
    arrange(Dataset)
  
  print(summary_table)
  cat(sprintf("\n results saved to: %s\n", OUT))
}
