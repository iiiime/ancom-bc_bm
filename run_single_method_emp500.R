#!/usr/bin/env Rscript
# run_single_method_emp500.R
# Run a single DA method on EMP500 data.
# Called as subprocess via /usr/bin/time -v for memory measurement.
#
# Usage:
#   Rscript run_single_method_emp500.R --method aldex2 --replicate 1 \
#     --data-dir /workspace/benchmark/data --output-dir /workspace/benchmark/results

# --- Parse command-line arguments ---
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  opt <- list(
    method = "aldex2",
    replicate = 1L,
    data_dir = "/workspace/benchmark/data",
    output_dir = "/workspace/benchmark/results",
    wrapper_path = "method_wrappers_emp500.R"
  )
  i <- 1
  while (i <= length(args)) {
    if (args[i] %in% c("--method", "--data-dir", "--output-dir", "--wrapper-path")) {
      key <- sub("^--", "", args[i])
      key <- gsub("-", "_", key)
      opt[[key]] <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--replicate") {
      opt$replicate <- as.integer(args[i + 1])
      i <- i + 2
    } else {
      i <- i + 1
    }
  }
  opt
}

opt <- parse_args()

cat(sprintf("[run_single_method] method=%s replicate=%d\n", opt$method, opt$replicate))
cat(sprintf("[run_single_method] data-dir=%s output-dir=%s\n", opt$data_dir, opt$output_dir))

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load data ---
data_file <- file.path(opt$data_dir, "emp500_sim.rds")
if (!file.exists(data_file)) {
  stop(sprintf("Data file not found: %s", data_file))
}
sim_data <- readRDS(data_file)
counts <- sim_data$counts       # taxa x samples
metadata <- sim_data$metadata   # samples x vars

cat(sprintf("[run_single_method] Data: %d taxa x %d samples\n", nrow(counts), ncol(counts)))

# --- Filter all-zero taxa only (no prevalence filter per user choice) ---
n_before <- nrow(counts)
counts <- counts[rowSums(counts) > 0, , drop = FALSE]  # remove all-zero taxa only
n_after <- nrow(counts)
if (n_before != n_after) {
  cat(sprintf("[run_single_method] Removed %d all-zero taxa (%d -> %d)\n",
              n_before - n_after, n_before, n_after))
}

# --- Source method wrappers ---
wrapper_file <- opt$wrapper_path
if (!file.exists(wrapper_file)) {
  stop(sprintf("method_wrappers_emp500.R not found at: %s", wrapper_file))
}
source(wrapper_file)

# --- Preprocess: remove zero-count samples, group rare factor levels ---
pp <- preprocess_emp500(counts, metadata, min_level_size = 5)
counts <- pp$counts
metadata <- pp$metadata
cat(sprintf("[run_single_method] After preprocessing: %d taxa x %d samples, empo_3=%d levels, env_biome=%d levels\n",
            nrow(counts), ncol(counts), nlevels(metadata$empo_3), nlevels(metadata$env_biome)))

# --- Get method function ---
method_funcs <- get_method_functions()
method_pkgs <- get_method_packages()

if (!(opt$method %in% names(method_funcs))) {
  stop(sprintf("Unknown method: %s. Available: %s", opt$method,
               paste(names(method_funcs), collapse = ", ")))
}

# --- Check if required package is installed ---
required_pkg <- method_pkgs[[opt$method]]
if (!requireNamespace(required_pkg, quietly = TRUE)) {
  cat(sprintf("[run_single_method] ERROR: Package '%s' not installed for method '%s'\n",
              required_pkg, opt$method))
  quit(status = 2)
}

# --- Run the method ---
cat(sprintf("[run_single_method] Running %s...\n", opt$method))

set.seed(1000 + opt$replicate)
set_single_threaded()

result <- tryCatch(
  withCallingHandlers({
    method_func <- method_funcs[[opt$method]]
    res <- method_func(counts, metadata)

    out_file <- file.path(opt$output_dir,
                          sprintf("emp500_%s_rep%d_results.rds",
                                  opt$method, opt$replicate))
    saveRDS(res, out_file)
    cat(sprintf("[run_single_method] Results saved to %s\n", out_file))
    cat(sprintf("[run_single_method] Results: %d taxa, %d significant (q<0.05)\n",
                nrow(res), sum(res$q_value < 0.05, na.rm = TRUE)))
    "SUCCESS"
  }, warning = function(w) {
    cat(sprintf("[run_single_method] WARNING: %s\n", conditionMessage(w)))
    invokeRestart("muffleWarning")
  }),
  error = function(e) {
    cat(sprintf("[run_single_method] ERROR: %s\n", conditionMessage(e)))
    "ERROR"
  }
)

cat(sprintf("[run_single_method] Status: %s\n", result))
quit(status = if (result == "SUCCESS") 0 else 1)
