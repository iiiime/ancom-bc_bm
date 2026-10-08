# Shared, repository-relative I/O helpers for the benchmark scripts.

command_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", command_args, value = TRUE)
PROJECT_ROOT <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = FALSE))
} else {
  normalizePath(".", mustWork = FALSE)
}
DATA_DIR <- file.path(PROJECT_ROOT, "data")
RESULTS_DIR <- file.path(PROJECT_ROOT, "results")

input_path <- function(filename, data_dir = DATA_DIR) {
  preferred <- file.path(data_dir, filename)
  legacy <- file.path(PROJECT_ROOT, filename)
  if (file.exists(preferred)) return(preferred)
  if (file.exists(legacy)) return(legacy)
  stop(sprintf(
    "Dataset not found: %s. Download it using data/manifest.yaml.",
    preferred
  ), call. = FALSE)
}

output_path <- function(filename, output_dir = RESULTS_DIR) {
  path <- file.path(output_dir, filename)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path
}

atomic_write_csv <- function(data, path, row.names = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(".", basename(path), "."),
                        tmpdir = dirname(path), fileext = ".csv")
  on.exit(unlink(temporary), add = TRUE)
  write.csv(data, temporary, row.names = row.names)
  if (!file.rename(temporary, path)) {
    stop(sprintf("Could not replace output file: %s", path), call. = FALSE)
  }
  invisible(path)
}
