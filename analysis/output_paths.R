# Set ALFALFA_RESULTS_DIR to compare regenerated outputs in a separate directory.
# Input matrices always come from data/processed/ in the repository.
results_dir <- Sys.getenv("ALFALFA_RESULTS_DIR", unset = here::here("results"))
result_path <- function(...) {
  path <- file.path(results_dir, ...)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path
}
