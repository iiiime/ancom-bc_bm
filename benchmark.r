#!/usr/bin/env Rscript
# Benchmark R ancombc() and ancombc2()
# 22 scenarios: vary samples, vary features, vary covariates, full 10k×10k*10cov
# 2 implementations: ancombc, ancombc2
# 3 repeats per scenario (after 1 warm-up)
# Memory: gc() + OS-level RSS
# Single-thread: n_cl=1

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = FALSE))
} else {
  normalizePath(".", mustWork = FALSE)
}
source(file.path(script_dir, "benchmark_io.R"))

library(ANCOMBC)
library(tidyverse)

COUNTS_FILE <- input_path("sim_counts_10k_20_meta.csv.gz")
META_FILE <- input_path("sim_metadata_10k_20_meta.csv")
OUTPUT <- output_path("bench_r_results_pseudo_sens.csv")
N_REPEATS <- 3


scenarios_vary_n <- data.frame(
  N = c(10, 100, 1000, 2000, 4000, 6000, 8000, 10000),
  P = 1000,
  Regime = "Vary_Samples",
  N_Cov = 2
)

scenarios_vary_p <- data.frame(
  N = 1000,
  P = c(10, 100, 1000, 2000, 4000, 6000, 8000, 10000),
  Regime = "Vary_Features",
  N_Cov = 2
)

scenarios_vary_cov <- data.frame(
  N = 1000,
  P = 1000,
  Regime = "Vary_Covariates",
  N_Cov = c(2, 4, 6, 8, 10)
)

scenarios_full <- data.frame(
  N = 10000,
  P = 10000,
  Regime = "Full_Matrix",
  N_Cov = 10
)

scenarios <- rbind(scenarios_vary_n, scenarios_vary_p, scenarios_vary_cov, scenarios_full)
scenarios <- scenarios[!duplicated(scenarios[, c("N", "P", "N_Cov")]), ]


build_formula <- function(n_cov) {
  parts <- c()
  for (j in 1:(n_cov %/% 2)) {
    parts <- c(parts, paste0("cont_cov_", j), paste0("cat_cov_", j))
  }
  paste(parts, collapse = " + ")
}


read_proc_status <- function() {
  if (!file.exists("/proc/self/status")) {
    return(list(VmRSS_MB = NA_real_))
  }
  lines <- readLines("/proc/self/status", warn = FALSE)
  vmrss_line <- grep("^VmRSS:", lines, value = TRUE)

  parse_kb <- function(line) {
    val <- gsub("[^0-9]", "", line)
    as.numeric(val) / 1024  # KB -> MB
  }

  vmrss_mb <- if (length(vmrss_line) > 0) parse_kb(vmrss_line[1]) else NA

  list(VmRSS_MB = vmrss_mb)
}


run_ancombc_r <- function(counts_mat, meta_df, formula_str) {
  res <- ancombc(
    data = counts_mat,
    meta_data = meta_df,
    formula = formula_str,
    p_adj_method = "holm",
    prv_cut = 0.0,
    lib_cut = 0,
    group = "cat_cov_1",
    struc_zero = FALSE,
    neg_lb = FALSE,
    alpha = 0.05,
    n_cl = 1,
    verbose = FALSE
  )
  return(res$res)
}


run_ancombc2_r <- function(counts_mat, meta_df, formula_str) {
  meta_df$cat_cov_1 <- factor(meta_df$cat_cov_1)
  res <- ancombc2(
    data = counts_mat,
    meta_data = meta_df,
    fix_formula = formula_str,
    p_adj_method = "holm",
    prv_cut = 0.0,
    lib_cut = 0,
    group = "cat_cov_1",
    struc_zero = FALSE,
    neg_lb = FALSE,
    alpha = 0.05,
    pseudo = 1,
    pseudo_sens = TRUE,
    n_cl = 1,
    verbose = FALSE
  )
  return(res$res)
}


run_one <- function(counts_mat, meta_df, formula_str, impl) {
  baseline_status <- read_proc_status()

  gc(reset = TRUE)
  initial_gc <- sum(gc()[, 2])

  start_time <- proc.time()

  status <- "Success"
  error_msg <- NA
  res <- NULL

  tryCatch({
    if (impl == "ancombc") {
      res <- run_ancombc_r(counts_mat, meta_df, formula_str)
    } else if (impl == "ancombc2") {
      res <- run_ancombc2_r(counts_mat, meta_df, formula_str)
    }
  }, error = function(e) {
    status <<- "Failed"
    error_msg <<- substr(as.character(e$message), 1, 200)
    res <<- NULL
  })

  end_time <- proc.time()
  elapsed_sec <- as.numeric((end_time - start_time)["elapsed"])

  final_gc <- gc()
  gc_peak_mb <- sum(final_gc[, 6])
  gc_delta_mb <- gc_peak_mb - initial_gc

  # RSS
  final_status <- read_proc_status()
  rss_baseline_mb <- baseline_status$VmRSS_MB
  rss_after_mb <- final_status$VmRSS_MB
  rss_delta_mb <- rss_after_mb - rss_baseline_mb

  list(
    res = res,
    success = (status == "Success"),
    time_sec = elapsed_sec,
    gc_peak_mb = gc_peak_mb,
    gc_delta_mb = gc_delta_mb,
    rss_baseline_mb = rss_baseline_mb,
    rss_peak_mb = rss_after_mb,
    rss_delta_mb = rss_delta_mb,
    error = error_msg
  )
}


# Main
cat("Loading data...\n")
counts_df <- read.csv(COUNTS_FILE, row.names = 1)
counts_df <- as.matrix(counts_df)
counts_df <- t(counts_df)  # features × samples (taxa_are_rows = TRUE)

meta_df <- read.csv(META_FILE)
rownames(meta_df) <- paste0("S", seq_len(nrow(meta_df)))

total_p <- nrow(counts_df)
total_n <- ncol(counts_df)
cat(sprintf("  Data loaded: %d features × %d samples\n", total_p, total_n))

# Load existing results
results_log <- data.frame()
done_keys <- list()
if (file.exists(OUTPUT)) {
  existing <- read.csv(OUTPUT, stringsAsFactors = FALSE)
  results_log <- existing
  success_existing <- existing[existing$Status == "Success", ]
  for (i in seq_len(nrow(success_existing))) {
    r <- success_existing[i, ]
    done_keys[[paste(r$Regime, r$Samples, r$Features, r$Impl, r$Run_ID, sep = "|")]] <- TRUE
  }
  cat(sprintf("  Resuming: %d existing records\n", nrow(results_log)))
}

total_scenarios <- nrow(scenarios)

for (idx in seq_len(nrow(scenarios))) {
  curr_n <- scenarios$N[idx]
  curr_p <- scenarios$P[idx]
  regime <- scenarios$Regime[idx]
  n_cov <- scenarios$N_Cov[idx]
  formula_str <- build_formula(n_cov)

  cat(sprintf("\n%s\n", paste(rep("=", 60), collapse = "")))
  cat(sprintf("[%d/%d] %s | %dx%d | cov=%d\n", idx, total_scenarios, regime, curr_n, curr_p, n_cov))
  cat(sprintf("  Formula: %s\n", formula_str))

  sub_counts <- counts_df[1:curr_p, 1:curr_n, drop = FALSE]
  sub_meta <- meta_df[1:curr_n, , drop = FALSE]

  for (impl in c("ancombc", "ancombc2")) {
    cat(sprintf("\n  --- %s ---\n", impl))

    cat("  warm-up... ")
    w <- run_one(sub_counts, sub_meta, formula_str, impl)
    if (w$success) {
      cat(sprintf("OK (%.2fs)\n", w$time_sec))
      fname <- output_path(file.path(
        "simulation",
        sprintf("impl_%s_s_%d_f_%d_cov_%d.csv", impl, curr_n, curr_p, n_cov)
      ))
      atomic_write_csv(w$res, fname, row.names = TRUE)
    } else {
      cat(sprintf("FAILED: %s\n", w$error))
      next
    }

    # Measured runs
    for (i in 1:N_REPEATS) {
      key <- paste(regime, curr_n, curr_p, paste0("R_", impl), i, sep = "|")
      if (!is.null(done_keys[[key]])) {
        cat(sprintf("  run %d/%d (cached)\n", i, N_REPEATS))
        next
      }

      cat(sprintf("  run %d/%d... ", i, N_REPEATS))
      r <- run_one(sub_counts, sub_meta, formula_str, impl)

      if (r$success) {
        cat(sprintf(
          "%.2fs | gc delta: %.1fMB | RSS delta: %.1fMB\n",
          r$time_sec, r$gc_delta_mb, r$rss_delta_mb
        ))
      } else {
        cat(sprintf("FAILED: %s\n", r$error))
      }

      entry <- data.frame(
        Regime = regime,
        Samples = curr_n,
        Features = curr_p,
        N_Cov = n_cov,
        Impl = paste0("R_", impl),
        Run_ID = i,
        Time_Sec = r$time_sec,
        Mem_Trace_MB = r$gc_delta_mb,
        RSS_Baseline_MB = r$rss_baseline_mb,
        RSS_Peak_MB = r$rss_peak_mb,
        RSS_Delta_MB = r$rss_delta_mb,
        Status = if (r$success) "Success" else "Failed",
        Error = if (!r$success) as.character(r$error) else "",
        Timestamp = as.character(Sys.time()),
        stringsAsFactors = FALSE
      )

      results_log <- rbind(results_log, entry)
      atomic_write_csv(results_log, OUTPUT)
    }
  }

  gc()
}

atomic_write_csv(results_log, OUTPUT)
cat(sprintf("\n%s\n", paste(rep("=", 60), collapse = "")))
cat(sprintf("Benchmark complete. Results saved to %s\n", OUTPUT))

if (nrow(results_log) > 0) {
  success_df <- results_log[results_log$Status == "Success", ]
  if (nrow(success_df) > 0) {
    summary_table <- success_df %>%
      group_by(Regime, Samples, Features, Impl) %>%
      summarise(
        Avg_Time_Sec = mean(Time_Sec),
        Avg_Mem_Trace_MB = mean(Mem_Trace_MB),
        Avg_RSS_Delta_MB = mean(RSS_Delta_MB),
        N_Runs = n(),
        .groups = "drop"
      ) %>%
      arrange(Regime, Samples, Features, Impl)

    cat("\nSummary:\n")
    print(summary_table)
  }
}
