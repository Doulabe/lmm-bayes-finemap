# ==============================================================================
# 40_score_ablation_exact.R
# Bayes-factor ablation: Conditional Score-LMM = the EXACT CBF-LMM machinery
# (candidate-specific K_jk, rank-one evaluation, profile-ML eBIC stopping,
# empty set allowed) with the ranking replaced by the frequentist Q_jk^2 at
# the per-candidate REML plug-in node instead of the conditional Bayes factor.
# Answers: how much of the gain comes from the Bayes factor itself, rather
# than from conditional LMM ranking plus eBIC stopping?
#
# Cells (same seeds as the corresponding axis scripts):
#   anchor  (n=1000, m=5000, rho=0.95)  x sg in {0, 0.5}  x 100 reps
#   rho098  (n=1000, m=5000, rho=0.98)  x sg in {0, 0.5}  x 100 reps
# The CBF arm is NOT recomputed: it is read from the exact-campaign payloads
# (results/bench_full_exact/01_scaling_n and 02_scaling_rho).
#
# Output: results/bench_full_exact/40_score_ablation/ + summary CSV
# Usage:  Rscript sim/bench_full/40_score_ablation_exact.R [--cores 5]
# ==============================================================================

Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressPackageStartupMessages(library(parallel))
source("sim/bench_full/00_config.R")
source("R/CBF_LMM_exact.R")

OUT_DIR <- "results/bench_full_exact/40_score_ablation"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(flag, default) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) as.integer(args[i + 1L]) else default
}
N_CORES <- arg_val("--cores", 5L)

one_rep <- function(cell, sg, b) {
  ck <- file.path(OUT_DIR, sprintf("%s_sg%.1f_b%03d.rds", cell, sg, b))
  if (file.exists(ck) && file.size(ck) > 0) return("cached")
  if (cell == "anchor") {
    rho  <- 0.95
    seed <- 20260425L + 1000L * (b - 1L) + 1000L + as.integer(round(1000 * sg))
  } else {
    rho  <- 0.98
    seed <- 20260425L + 1000L * (b - 1L) +
            as.integer(round(1000 * rho)) + as.integer(round(100 * sg))
  }
  d <- gen_dataset(n = 1000L, m = ANCHOR$m, rho = rho,
                   K_true = ANCHOR$K_true, beta_true = ANCHOR$beta_true,
                   sigma_g2 = sg, block_size = ANCHOR$block_size, seed = seed)
  fit <- CBF_LMM_stepwise_exact(d$y, d$X, tau2 = ANCHOR$tau2,
                                  K_max = ANCHOR$K_max,
                                  n_nodes = ANCHOR$N_delta,
                                  rank_by = "score")
  tp <- length(intersect(fit$indices, d$truth))
  saveRDS(list(cell = cell, sg = sg, b = b, indices = fit$indices,
               K_hat = fit$K_hat, tp = tp,
               recall = tp / length(d$truth),
               precision = if (fit$K_hat > 0) tp / fit$K_hat else 0), ck)
  "done"
}

grid <- expand.grid(cell = c("anchor", "rho098"), sg = c(0, 0.5),
                    b = 1:100, stringsAsFactors = FALSE)
message(sprintf("Score-LMM exact ablation: %d replicates on %d cores",
                nrow(grid), N_CORES))
invisible(mclapply(seq_len(nrow(grid)), function(i)
  tryCatch(one_rep(grid$cell[i], grid$sg[i], grid$b[i]),
           error = function(e) message(sprintf("FAIL %s sg%.1f b%d: %s",
             grid$cell[i], grid$sg[i], grid$b[i], conditionMessage(e)))),
  mc.cores = N_CORES, mc.preschedule = FALSE))

# ── summary: Score arm vs stored CBF-exact arm on the same datasets ──────────
read_cbf <- function(axis, pat) {
  fs <- list.files(file.path("results/bench_full_exact", axis),
                   pattern = pat, full.names = TRUE)
  do.call(rbind, lapply(fs, function(f) {
    m <- readRDS(f)$metrics
    m <- m[m$method == "MS_L_eBIC",
           c("K_hat", "recall", "precision", "f1")]
    m$sg <- readRDS(f)$metrics$sigma_g2[1]; m
  }))
}
sc_fs <- list.files(OUT_DIR, pattern = "^(anchor|rho098).*rds$",
                    full.names = TRUE)
sc <- do.call(rbind, lapply(sc_fs, function(f) {
  r <- readRDS(f)
  f1 <- if (r$precision + r$recall > 0)
          2 * r$precision * r$recall / (r$precision + r$recall) else 0
  data.frame(cell = r$cell, sg = r$sg, K_hat = r$K_hat, recall = r$recall,
             precision = r$precision, f1 = f1)
}))
cat("\n════ Score-LMM (exact) ════\n")
print(aggregate(cbind(K_hat, recall, precision, f1) ~ cell + sg, sc, mean),
      digits = 3, row.names = FALSE)
cbf_anchor <- read_cbf("01_scaling_n", "^n1000_sg[0-9.]+_b[0-9]+[.]rds$")
cbf_rho    <- read_cbf("02_scaling_rho", "^rho0[.]98_sg[0-9.]+_b[0-9]+[.]rds$")
cbf_anchor$cell <- "anchor"; cbf_rho$cell <- "rho098"
cb <- rbind(cbf_anchor, cbf_rho)
cat("\n════ CBF-LMM (exact, campagne) ════\n")
print(aggregate(cbind(K_hat, recall, precision, f1) ~ cell + sg, cb, mean),
      digits = 3, row.names = FALSE)
write.csv(aggregate(cbind(K_hat, recall, precision, f1) ~ cell + sg, sc, mean),
          file.path(OUT_DIR, "score_ablation_summary.csv"), row.names = FALSE)
message("done.")
