# ==============================================================================
# 08b_fix_fastlmm_failures.R
# Re-run fastlmm_score on any (cell, rep) where the recorded fastlmm has
# elapsed < 0.5s AND K_hat == 0 (indicative of a silent crash inside the
# worker — most likely OOM at n=3000 with too many parallel workers, or
# a rrBLUP requireNamespace failure).
#
# This is idempotent: only failed entries are re-run.  Successful entries
# are preserved.  Use a small `--cores` value to keep memory headroom.
#
# Usage:
#   Rscript sim/bench_full/08b_fix_fastlmm_failures.R --cores 8
# ==============================================================================

source("sim/bench_full/00_config.R")
source("sim/bench_full/08_add_fastlmm.R")   # imports seed_for_cell, context_for_cell, augment_one

FAILURE_ELAPSED_S <- 0.5   # below this, treat as a silent crash

# -- Identify failures ---------------------------------------------------------
list_failures <- function() {
  dirs <- c("01" = "01_scaling_n", "02" = "02_scaling_rho",
              "03" = "03_arch",     "04" = "04_scaling_m")
  rows <- list()
  for (dir_id in names(dirs)) {
    rdir <- file.path("results/bench_full", dirs[[dir_id]])
    if (!dir.exists(rdir)) next
    files <- list.files(rdir, pattern = "_b\\d+\\.rds$", full.names = TRUE)
    for (f in files) {
      p <- tryCatch(readRDS(f), error = function(e) NULL)
      if (is.null(p$results$fastlmm)) next
      e <- p$results$fastlmm$elapsed
      k <- p$results$fastlmm$K_hat
      if (!is.null(e) && e < FAILURE_ELAPSED_S && k == 0) {
        rows[[length(rows) + 1L]] <- list(file = f, dir_id = dir_id,
                                             elapsed = e, K_hat = k)
      }
    }
  }
  rows
}

# -- Re-run on one failure -----------------------------------------------------
rerun_one <- function(item) {
  # Strip existing fastlmm entry, then re-augment
  p <- readRDS(item$file)
  p$metrics <- p$metrics[p$metrics$method != "fastlmm", , drop = FALSE]
  p$results$fastlmm <- NULL
  saveRDS(p, item$file)
  augment_one(item$file, item$dir_id)
}


if (sys.nframe() == 0L) {
  args <- parse_args()
  n_cores <- if (!is.null(args$cores)) as.integer(args$cores) else 8L
  cat(sprintf("=== 08b_fix_fastlmm_failures: cores = %d ===\n", n_cores))

  Sys.setenv(OMP_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
              OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

  # CRITICAL: pre-load rrBLUP into the parent process BEFORE forking.
  # Without this, each mclapply worker calls requireNamespace("rrBLUP")
  # concurrently on its first invocation; under heavy load this can race
  # and return FALSE silently, causing fastlmm_score to bail with K_hat=0
  # in milliseconds.  Pre-loading guarantees every fork inherits an already
  # initialised rrBLUP namespace.
  cat("Pre-loading rrBLUP in parent process before fork...\n")
  if (!requireNamespace("rrBLUP", quietly = TRUE))
    stop("rrBLUP is not installed — install it before running this script.")
  loadNamespace("rrBLUP")  # forces full load, not just namespace registration
  cat("  rrBLUP version:", as.character(packageVersion("rrBLUP")), "\n")

  failures <- list_failures()
  cat(sprintf("Found %d fastlmm failures (elapsed < %.1fs and K_hat == 0)\n",
              length(failures), FAILURE_ELAPSED_S))
  if (length(failures) == 0L) {
    cat("Nothing to do.\n"); quit("no", status = 0L)
  }
  # Show breakdown
  dir_counts <- table(sapply(failures, function(x) x$dir_id))
  for (d in names(dir_counts))
    cat(sprintf("  dir %s : %d failures\n", d, dir_counts[[d]]))

  cat("Re-running...\n")
  results <- parallel::mclapply(failures, function(it) {
    tryCatch(rerun_one(it),
              error = function(e) {
                message("ERR ", basename(it$file), ": ",
                         conditionMessage(e))
                "error"
              })
  }, mc.cores = n_cores, mc.preschedule = FALSE)
  tbl <- table(unlist(results))
  cat(sprintf("Result: %s\n",
              paste(names(tbl), tbl, sep = "=", collapse = ", ")))
}
