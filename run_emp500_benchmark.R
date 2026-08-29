#!/usr/bin/env Rscript
# run_emp500_benchmark.R
# Main orchestration: runs all 15 DA methods on EMP500 data (3 replicates each).
# 14 R methods via Rscript + 1 Python method via python3.
# Each method run is a separate subprocess via /usr/bin/time -v
# to measure wall-clock runtime and peak memory (RSS).

# --- Configuration ---
data_dir <- "../data"
output_dir <- "../results"
wrapper_path <- "method_wrappers_emp500.R"
single_method_script <- "../scripts/run_single_method_emp500.R"
python_script <- "../scripts/run_ancombc_python.py"
results_csv <- file.path(output_dir, "emp500_benchmark_results.csv")
default_timeout <- 7200  # 120 minutes per method
# Per-method timeout overrides (in seconds)
method_timeouts <- list(
  orm = 18000,       # 5 hours (fits 8351 separate models)
  maaslin2 = 18000,  # 3 hours (slow with plotting)
  maaslin3 = 7200,   # 2 hours
  corncob = 7200,    # 2 hours
  ancombc = 7200,    # 2 hours
  ancombc2 = 7200,   # 2 hours
  deseq2 = 7200      # 2 hours
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


all_methods <- c("limma_voom", "zicoseq", "edger",
                 "metagenomeseq", "aldex2", "ancombc2", "ancombc",
                 "deseq2", "maaslin3", "corncob", "orm")
replicates <- 1:3

# --- Helper: parse /usr/bin/time -v output ---
parse_time_output <- function(output_lines) {
  runtime_s <- NA
  peak_mem_kb <- NA
  for (line in output_lines) {
    if (grepl("Elapsed.*wall clock.*time", line, ignore.case = TRUE)) {
      time_str <- sub(".*:\\s*", "", line)
      time_str <- trimws(time_str)
      parts <- strsplit(time_str, ":")[[1]]
      if (length(parts) == 2) {
        runtime_s <- as.numeric(parts[1]) * 60 + as.numeric(parts[2])
      } else if (length(parts) == 3) {
        runtime_s <- as.numeric(parts[1]) * 3600 + as.numeric(parts[2]) * 60 + as.numeric(parts[3])
      }
    }
    if (grepl("Maximum resident set size", line, ignore.case = TRUE)) {
      mem_str <- sub(".*:\\s*", "", line)
      peak_mem_kb <- as.numeric(trimws(mem_str))
    }
  }
  list(runtime_s = runtime_s, peak_mem_kb = peak_mem_kb)
}

# --- Helper: run a single R method with timeout ---
run_single_r <- function(method, replicate) {
  cat(sprintf("\n[%s] Running %s (rep %d)...\n",
              format(Sys.time(), "%H:%M:%S"), method, replicate))

  timeout_seconds <- if (!is.null(method_timeouts[[method]])) method_timeouts[[method]] else default_timeout

  rscript_args <- c(
    single_method_script,
    "--method", method,
    "--replicate", as.character(replicate),
    "--data-dir", data_dir,
    "--output-dir", output_dir,
    "--wrapper-path", wrapper_path
  )

  start_time <- Sys.time()
  output_lines <- system2("/usr/bin/time",
                          args = c("-v", "Rscript", rscript_args),
                          stdout = TRUE, stderr = TRUE,
                          timeout = timeout_seconds)
  exit_code <- attr(output_lines, "status")
  if (is.null(exit_code)) exit_code <- 0
  if (is.null(output_lines)) output_lines <- character(0)
  end_time <- Sys.time()
  wall_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  parsed <- parse_time_output(output_lines)

  if (exit_code == 124) {
    status <- "TIMEOUT"
  } else if (exit_code == 2) {
    status <- "INSTALL_FAILED"
  } else if (exit_code != 0) {
    status <- "ERROR"
  } else {
    status <- "SUCCESS"
  }

  error_msg <- ""
  if (status != "SUCCESS") {
    error_lines <- grep("^\\[run_single_method\\] (ERROR|WARNING)", output_lines, value = TRUE)
    if (length(error_lines) > 0) {
      error_msg <- paste(error_lines, collapse = "; ")
    } else {
      error_msg <- paste(tail(output_lines, 5), collapse = "; ")
    }
  }

  runtime_s <- if (!is.na(parsed$runtime_s)) parsed$runtime_s else wall_time
  peak_mem_mb <- if (!is.na(parsed$peak_mem_kb)) parsed$peak_mem_kb / 1024 else NA

  cat(sprintf("  -> Status: %s, Runtime: %.1fs, Peak Memory: %.1f MB\n",
              status, runtime_s, peak_mem_mb))

  data.frame(
    method = method, replicate = replicate,
    runtime_s = runtime_s, peak_memory_MB = peak_mem_mb,
    status = status, error_message = substr(error_msg, 1, 500),
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE
  )
}

# --- Helper: run a single Python method with timeout ---
run_single_python <- function(method, replicate) {
  cat(sprintf("\n[%s] Running %s (rep %d)...\n",
              format(Sys.time(), "%H:%M:%S"), method, replicate))

  timeout_seconds <- if (!is.null(method_timeouts[[method]])) method_timeouts[[method]] else default_timeout

  py_args <- c(
    python_script,
    "--replicate", as.character(replicate),
    "--data-dir", data_dir,
    "--output-dir", output_dir
  )

  start_time <- Sys.time()
  output_lines <- system2("/usr/bin/time",
                          args = c("-v", "python3", py_args),
                          stdout = TRUE, stderr = TRUE,
                          timeout = timeout_seconds)
  exit_code <- attr(output_lines, "status")
  if (is.null(exit_code)) exit_code <- 0
  if (is.null(output_lines)) output_lines <- character(0)
  end_time <- Sys.time()
  wall_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  parsed <- parse_time_output(output_lines)

  if (exit_code == 124) {
    status <- "TIMEOUT"
  } else if (exit_code != 0) {
    status <- "ERROR"
  } else {
    status <- "SUCCESS"
  }

  error_msg <- ""
  if (status != "SUCCESS") {
    error_lines <- grep("^\\[run_ancombc_python\\] (ERROR|WARNING)", output_lines, value = TRUE)
    if (length(error_lines) > 0) {
      error_msg <- paste(error_lines, collapse = "; ")
    } else {
      error_msg <- paste(tail(output_lines, 10), collapse = "; ")
    }
  }

  runtime_s <- if (!is.na(parsed$runtime_s)) parsed$runtime_s else wall_time
  peak_mem_mb <- if (!is.na(parsed$peak_mem_kb)) parsed$peak_mem_kb / 1024 else NA

  cat(sprintf("  -> Status: %s, Runtime: %.1fs, Peak Memory: %.1f MB\n",
              status, runtime_s, peak_mem_mb))

  data.frame(
    method = method, replicate = replicate,
    runtime_s = runtime_s, peak_memory_MB = peak_mem_mb,
    status = status, error_message = substr(error_msg, 1, 500),
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE
  )
}

# --- Main benchmark loop ---
cat("=== EMP500 DA Method Benchmark ===\n")
cat(sprintf("R methods: %d (%s)\n", length(all_methods), paste(all_methods, collapse = ", ")))

cat(sprintf("Replicates: %d\n", length(replicates)))
cat(sprintf("Total runs: %d\n", length(all_methods) * length(replicates)))
cat(sprintf("Default timeout: %d seconds per run\n", default_timeout))
cat(sprintf("Output: %s\n\n", results_csv))

# Initialize results CSV
if (!file.exists(results_csv)) {
  write.csv(data.frame(method = character(), replicate = integer(),
                        runtime_s = numeric(), peak_memory_MB = numeric(),
                        status = character(), error_message = character(),
                        timestamp = character(),
                        stringsAsFactors = FALSE),
            results_csv, row.names = FALSE)
}

# Read existing results to skip completed runs
existing <- read.csv(results_csv, stringsAsFactors = FALSE)
completed <- paste(existing$method, existing$replicate, sep = ":")

total_runs <- length(all_methods) * length(replicates)
run_count <- 0

for (rep in replicates) {
  for (meth in all_methods) {
    run_count <- run_count + 1
    run_key <- paste(meth, rep, sep = ":")

    if (run_key %in% completed) {
      cat(sprintf("\n[%d/%d] SKIP (already done): %s rep %d\n",
                  run_count, total_runs, meth, rep))
      next
    }

    cat(sprintf("\n[%d/%d] ", run_count, total_runs))

    result_row <- run_single_r(meth, rep)

    write.table(result_row, results_csv, append = TRUE,
                sep = ",", row.names = FALSE, col.names = FALSE)
    completed <- c(completed, run_key)

    gc()
    Sys.sleep(2)
  }
}

cat("\n\n=== Benchmark Complete ===\n")
cat(sprintf("Results saved to: %s\n", results_csv))

final_results <- read.csv(results_csv, stringsAsFactors = FALSE)
cat(sprintf("Total runs: %d\n", nrow(final_results)))
cat(sprintf("Successful: %d\n", sum(final_results$status == "SUCCESS")))
cat(sprintf("Errors: %d\n", sum(final_results$status == "ERROR")))
cat(sprintf("Timeouts: %d\n", sum(final_results$status == "TIMEOUT")))
cat(sprintf("Install failures: %d\n", sum(final_results$status == "INSTALL_FAILED")))
