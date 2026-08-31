# ==============================================================================
# 30_exact_rerun.R
#
# Recompute the primary CBF-LMM rows of the controlled benchmark with the
# EXACT candidate-specific implementation (R/CBF_LMM_exact.R):
#   * model kernel K_jk (candidate excluded from its own background),
#     evaluated by rank-one downdates of the pool eigendecomposition;
#   * numerical marginalisation of delta (15 nodes, half-Cauchy prior);
#   * TRUE profile-ML eBIC stopping on the full accepted model, first
#     candidate included (the empty set is a possible return).
#
# The datasets are reconstructed with the ORIGINAL per-axis seed formulas, so
# every comparator row (SuSiE, BSLMM, BayesR, FaST-LMM, ...) stored in the
# original payloads remains valid: this script clones each payload, replaces
# the MS_L_eBIC entry by the exact fit (the shared-kernel fit is kept under
# MS_L_eBIC_sharedGk), recomputes the metrics table, and writes to
# results/bench_full_exact/<axis>/<same filename>. Originals are not touched.
#
# Checkpointed per replicate (skip if output exists) — safe to kill/resume.
#
# Usage:
#   Rscript sim/bench_full/30_exact_rerun.R [--cores 5]
# ==============================================================================

Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressPackageStartupMessages(library(parallel))
source("sim/bench_full/00_config.R")
source("R/CBF_LMM_exact.R")

IN_ROOT  <- "results/bench_full"
OUT_ROOT <- "results/bench_full_exact"

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(flag, default) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) as.integer(args[i + 1L]) else default
}
N_CORES <- arg_val("--cores", 5L)

# ── work-list: every payload of the four controlled axes ─────────────────────
AXES <- c("01_scaling_n", "02_scaling_rho", "03_arch", "04_scaling_m")
work <- do.call(rbind, lapply(AXES, function(ax) {
  fs <- list.files(file.path(IN_ROOT, ax), pattern = "_b[0-9]+\\.rds$")
  if (!length(fs)) return(NULL)
  data.frame(axis = ax, file = fs, stringsAsFactors = FALSE)
}))

parse_cell <- function(axis, file) {
  rep_id <- as.integer(sub(".*_b([0-9]+)\\.rds$", "\\1", file))
  sg     <- as.numeric(sub(".*_sg([0-9.]+)_b[0-9]+\\.rds$", "\\1", file))
  out <- list(rep = rep_id, sg = sg,
              n = ANCHOR$n, m = ANCHOR$m, rho = ANCHOR$rho,
              beta = ANCHOR$beta_true)
  if (axis == "01_scaling_n") {
    out$n <- as.integer(sub("^n([0-9]+)_.*", "\\1", file))
    out$seed <- 20260425L + 1000L * (rep_id - 1L) + out$n +
                as.integer(round(1000 * sg))
  } else if (axis == "02_scaling_rho") {
    out$rho <- as.numeric(sub("^rho([0-9.]+)_sg.*", "\\1", file))
    out$seed <- 20260425L + 1000L * (rep_id - 1L) +
                as.integer(round(1000 * out$rho)) +
                as.integer(round(100 * sg))
  } else if (axis == "03_arch") {
    sig <- sub("^sig([a-z]+)_sg.*", "\\1", file)
    out$beta <- switch(sig,
                       weak   = c(0.4, 0.2, 0.2, 0.1, 0.1),
                       medium = c(0.8, 0.4, 0.4, 0.2, 0.2),
                       strong = c(1.6, 0.8, 0.8, 0.4, 0.4))
    out$seed <- 20260425L + 1000L * (rep_id - 1L) +
                switch(sig, weak = 1L, medium = 2L, strong = 3L) +
                as.integer(round(100 * sg))
  } else {
    out$m <- as.integer(sub("^m([0-9]+)_.*", "\\1", file))
    out$seed <- 20260425L + 1000L * (rep_id - 1L) + out$m +
                as.integer(round(100 * sg))
  }
  out
}

one_file <- function(axis, file) {
  out_dir <- file.path(OUT_ROOT, axis)
  out_fn  <- file.path(out_dir, file)
  if (file.exists(out_fn) && file.size(out_fn) > 0) return("cached")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  pc <- parse_cell(axis, file)
  d  <- gen_dataset(n = pc$n, m = pc$m, rho = pc$rho,
                    K_true = length(pc$beta), beta_true = pc$beta,
                    sigma_g2 = pc$sg, block_size = ANCHOR$block_size,
                    seed = pc$seed)

  t0  <- Sys.time()
  fit <- CBF_LMM_stepwise_exact(d$y, d$X, tau2 = ANCHOR$tau2,
                                  K_max = ANCHOR$K_max,
                                  delta_eval = "marginal",
                                  n_nodes = ANCHOR$N_delta)
  el  <- as.numeric(Sys.time() - t0, units = "secs")

  # readRDS occasionally fails transiently under concurrent forked reads
  # ("error reading from connection"); retry with backoff.
  orig <- NULL
  for (att in 1:6) {
    orig <- tryCatch(readRDS(file.path(IN_ROOT, axis, file)),
                     error = function(e) NULL)
    if (!is.null(orig)) break
    Sys.sleep(2 * att)
  }
  if (is.null(orig)) stop("unreadable after retries: ", file)
  stopifnot(identical(sort(orig$truth), sort(d$truth)))   # seed sanity check
  res  <- orig$results
  res$MS_L_eBIC_sharedGk <- res$MS_L_eBIC
  res$MS_L_eBIC <- list(method = "MS_L_eBIC", indices = fit$indices,
                        K_hat = fit$K_hat, scores = fit$log_BF_at_step,
                        elapsed = el)
  metrics <- compute_metrics(res, orig$truth)
  extra   <- setdiff(names(orig$metrics), names(metrics))
  for (cn in extra) metrics[[cn]] <- orig$metrics[[cn]][1L]
  saveRDS(list(metrics = metrics, truth = orig$truth, results = res), out_fn)
  "done"
}

# ── ordering: anchor-sized cells first, n=3000 last ──────────────────────────
prio <- function(axis, file) {
  n_val <- if (axis == "01_scaling_n")
             as.integer(sub("^n([0-9]+)_.*", "\\1", file)) else ANCHOR$n
  m_val <- if (axis == "04_scaling_m")
             as.integer(sub("^m([0-9]+)_.*", "\\1", file)) else ANCHOR$m
  n_val * 10 + (m_val > 5000)          # small n first, big m later, n=3000 last
}
work <- work[order(mapply(prio, work$axis, work$file)), ]

message(sprintf("exact rerun: %d replicates on %d cores", nrow(work), N_CORES))
invisible(mclapply(seq_len(nrow(work)), function(i) {
  w <- work[i, ]
  r <- tryCatch(one_file(w$axis, w$file),
                error = function(e) sprintf("FAIL %s/%s: %s",
                                            w$axis, w$file,
                                            conditionMessage(e)))
  if (startsWith(r, "FAIL")) message(r)
  r
}, mc.cores = N_CORES, mc.preschedule = FALSE))
message("done.")
