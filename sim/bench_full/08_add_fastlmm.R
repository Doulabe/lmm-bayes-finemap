# ==============================================================================
# 08_add_fastlmm.R
# Augment existing per-(cell, rep) RDS files in results/bench_full/0[1-4]_*/
# with FaST-LMM-equivalent results (single-pass LMM score test + Bonferroni).
#
# Rationale: the original B=100 sweep used REML_stepwise as the LMM-frequentist
# comparator, but greedy stepwise on block-AR(1) with rho=0.95 selects all
# proximal SNPs of the top hit before reaching the next causal — yielding F1=0
# uniformly. fastlmm_score() is the standard Lippert et al. 2011 approach:
# a single null LMM fit on the full X-derived GRM, then per-SNP whitened score
# test with Bonferroni at alpha=0.05.
#
# This script regenerates the data with the SAME seed used by the original
# script (01-04), runs ONLY fastlmm_score, and appends:
#   - one new row to payload$metrics with method = "fastlmm"
#   - one new entry to payload$results$fastlmm
# Other methods are untouched — no need to re-run B=100.
#
# Usage:
#   Rscript sim/bench_full/08_add_fastlmm.R                       # all 4 dirs, default cores
#   Rscript sim/bench_full/08_add_fastlmm.R --dir 01              # only dir 01
#   Rscript sim/bench_full/08_add_fastlmm.R --cores 30            # 30-worker mclapply
#   Rscript sim/bench_full/08_add_fastlmm.R --dir 02 --cores 30   # combine
# ==============================================================================

source("sim/bench_full/00_config.R")

# -- Seed reconstruction (must match 01-04 scripts exactly) --------------------
seed_for_cell <- function(dir_id, cell_tag, rep_id) {
  base <- 20260425L + 1000L * (rep_id - 1L)
  if (dir_id == "01") {
    # n%04d_sg%.1f  -> seed = base + n + round(1000 * sg)
    n  <- as.integer(sub("^n(\\d{4})_sg.*$", "\\1", cell_tag))
    sg <- as.numeric(sub("^n\\d{4}_sg([0-9.]+)$", "\\1", cell_tag))
    return(base + n + as.integer(round(1000 * sg)))
  } else if (dir_id == "02") {
    # rho%.2f_sg%.1f -> seed = base + round(1000 * rho) + round(100 * sg)
    rho <- as.numeric(sub("^rho([0-9.]+)_sg.*$", "\\1", cell_tag))
    sg  <- as.numeric(sub("^rho[0-9.]+_sg([0-9.]+)$", "\\1", cell_tag))
    return(base + as.integer(round(1000 * rho)) + as.integer(round(100 * sg)))
  } else if (dir_id == "03") {
    # sig{weak,medium,strong}_sg%.1f -> base + switch + round(100 * sg)
    signal <- sub("^sig([a-z]+)_sg.*$", "\\1", cell_tag)
    sg     <- as.numeric(sub("^sig[a-z]+_sg([0-9.]+)$", "\\1", cell_tag))
    return(base + switch(signal, weak = 1L, medium = 2L, strong = 3L) +
             as.integer(round(100 * sg)))
  } else if (dir_id == "04") {
    # m%05d_sg%.1f -> base + m + round(100 * sg)
    m  <- as.integer(sub("^m(\\d{5})_sg.*$", "\\1", cell_tag))
    sg <- as.numeric(sub("^m\\d{5}_sg([0-9.]+)$", "\\1", cell_tag))
    return(base + m + as.integer(round(100 * sg)))
  }
  stop("Unknown dir_id: ", dir_id)
}

# -- Cell context from tag + dir -----------------------------------------------
context_for_cell <- function(dir_id, cell_tag) {
  if (dir_id == "01") {
    n  <- as.integer(sub("^n(\\d{4})_sg.*$", "\\1", cell_tag))
    sg <- as.numeric(sub("^n\\d{4}_sg([0-9.]+)$", "\\1", cell_tag))
    return(list(n = n, m = ANCHOR$m, rho = ANCHOR$rho,
                  K_true = ANCHOR$K_true, beta_true = ANCHOR$beta_true,
                  sigma_g2 = sg))
  } else if (dir_id == "02") {
    rho <- as.numeric(sub("^rho([0-9.]+)_sg.*$", "\\1", cell_tag))
    sg  <- as.numeric(sub("^rho[0-9.]+_sg([0-9.]+)$", "\\1", cell_tag))
    return(list(n = ANCHOR$n, m = ANCHOR$m, rho = rho,
                  K_true = ANCHOR$K_true, beta_true = ANCHOR$beta_true,
                  sigma_g2 = sg))
  } else if (dir_id == "03") {
    signal <- sub("^sig([a-z]+)_sg.*$", "\\1", cell_tag)
    sg     <- as.numeric(sub("^sig[a-z]+_sg([0-9.]+)$", "\\1", cell_tag))
    beta_true <- switch(signal,
                          weak    = c(0.4, 0.2, 0.2, 0.1, 0.1),
                          medium  = c(0.8, 0.4, 0.4, 0.2, 0.2),
                          strong  = c(1.6, 0.8, 0.8, 0.4, 0.4))
    return(list(n = ANCHOR$n, m = ANCHOR$m, rho = ANCHOR$rho,
                  K_true = length(beta_true), beta_true = beta_true,
                  sigma_g2 = sg))
  } else if (dir_id == "04") {
    m  <- as.integer(sub("^m(\\d{5})_sg.*$", "\\1", cell_tag))
    sg <- as.numeric(sub("^m\\d{5}_sg([0-9.]+)$", "\\1", cell_tag))
    return(list(n = ANCHOR$n, m = m, rho = ANCHOR$rho,
                  K_true = ANCHOR$K_true, beta_true = ANCHOR$beta_true,
                  sigma_g2 = sg))
  }
  stop("Unknown dir_id: ", dir_id)
}

# -- Process one RDS file ------------------------------------------------------
augment_one <- function(rds_path, dir_id) {
  payload <- readRDS(rds_path)
  if ("fastlmm" %in% payload$metrics$method) {
    return(invisible("cached"))
  }
  fn <- basename(rds_path)
  # Parse cell_tag and rep
  m <- regmatches(fn, regexec("^(.+)_b(\\d+)\\.rds$", fn))[[1]]
  if (length(m) != 3L) {
    warning("Cannot parse filename: ", fn); return(invisible("skip"))
  }
  cell_tag <- m[2]; rep_id <- as.integer(m[3])

  seed <- seed_for_cell(dir_id, cell_tag, rep_id)
  ctx  <- context_for_cell(dir_id, cell_tag)

  d <- gen_dataset(n = ctx$n, m = ctx$m, rho = ctx$rho,
                     K_true = ctx$K_true, beta_true = ctx$beta_true,
                     sigma_g2 = ctx$sigma_g2,
                     block_size = ANCHOR$block_size, seed = seed)
  if (!identical(as.numeric(d$truth), as.numeric(payload$truth))) {
    warning(sprintf("Seed mismatch in %s: truth differs", fn))
    return(invisible("mismatch"))
  }
  t0 <- Sys.time()
  res <- tryCatch(fastlmm_score(d$y, d$X, alpha = 0.05),
                    error = function(e) list(indices = integer(0),
                                              K_hat = 0L, scores = numeric(0)))
  elapsed <- as.numeric(Sys.time() - t0, units = "secs")
  tp <- length(intersect(res$indices, d$truth))
  fp <- res$K_hat - tp
  rec  <- tp / ctx$K_true
  prec <- if (res$K_hat > 0) tp / res$K_hat else 0
  f1   <- if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0

  # Build new metrics row matching existing column structure
  new_row <- payload$metrics[1, ]                # template
  new_row$method <- "fastlmm"
  new_row$K_hat  <- res$K_hat
  new_row$tp     <- tp
  new_row$fp     <- fp
  new_row$recall <- rec
  new_row$precision <- prec
  new_row$f1     <- f1
  new_row$elapsed <- elapsed
  # n, m, rho, sigma_g2, rep are preserved from template

  payload$metrics <- rbind(payload$metrics, new_row)
  payload$results$fastlmm <- list(method = "fastlmm",
                                      indices = res$indices,
                                      K_hat = res$K_hat,
                                      scores = res$scores,
                                      elapsed = elapsed)
  saveRDS(payload, rds_path)
  invisible("done")
}


# -- Main ----------------------------------------------------------------------
if (sys.nframe() == 0L) {
  # Force single-threaded BLAS to avoid Accelerate thread contention under mclapply
  Sys.setenv(OMP_NUM_THREADS = "1",
              OPENBLAS_NUM_THREADS = "1",
              MKL_NUM_THREADS = "1",
              VECLIB_MAXIMUM_THREADS = "1")

  # CRITICAL: pre-load rrBLUP into the parent process BEFORE forking.
  # Otherwise each mclapply worker calls requireNamespace("rrBLUP")
  # concurrently on its first invocation; under heavy load this can race
  # and silently return FALSE, making fastlmm_score bail with K_hat=0
  # in ~1ms.  Pre-loading guarantees every fork inherits an already
  # initialised rrBLUP namespace.
  if (!requireNamespace("rrBLUP", quietly = TRUE))
    stop("rrBLUP is not installed — install it before running this script.")
  loadNamespace("rrBLUP")
  cat(sprintf("Pre-loaded rrBLUP %s in parent process.\n",
              as.character(packageVersion("rrBLUP"))))

  args <- parse_args()
  dirs <- if (!is.null(args$dir)) sprintf("%02d", as.integer(args$dir)) else
            c("01", "02", "03", "04")
  dir_names <- c("01" = "01_scaling_n", "02" = "02_scaling_rho",
                  "03" = "03_arch",      "04" = "04_scaling_m")

  # CLI-configurable worker count.  Default: leave 2 cores free for the OS.
  n_cores <- if (!is.null(args$cores)) as.integer(args$cores) else
               max(1L, parallel::detectCores() - 2L)
  cat(sprintf("=== Using mc.cores = %d (detected %d, --cores override = %s) ===\n",
              n_cores, parallel::detectCores(),
              if (!is.null(args$cores)) args$cores else "none"))

  for (dir_id in dirs) {
    dn <- dir_names[[dir_id]]
    rdir <- file.path("results/bench_full", dn)
    if (!dir.exists(rdir)) {
      cat(sprintf("  [skip] %s not found\n", rdir)); next
    }
    files <- list.files(rdir, pattern = "_b\\d+\\.rds$", full.names = TRUE)
    cat(sprintf("=== %s : %d files ===\n", dn, length(files)))
    out <- parallel::mclapply(files, function(f)
      tryCatch(augment_one(f, dir_id),
                error = function(e) {
                  message("ERR ", basename(f), ": ", conditionMessage(e))
                  "error"
                }),
      mc.cores = n_cores, mc.preschedule = FALSE)
    tbl <- table(unlist(out))
    cat(sprintf("  %s\n", paste(names(tbl), tbl, sep = "=",
                                  collapse = ", ")))
  }
  cat("\nDone. Re-run sim/bench_full/99_aggregate.R to refresh summaries.\n")
}
