# ==============================================================================
# 31_exact_reml_anchor.R
# delta-sensitivity arm for the exact implementation (Table tab_delta):
# on the anchor cells (n=1000, m=5000, sg in {0, 0.5}, 100 reps, same seeds
# as 01_scaling_n), run CBF_LMM_stepwise_exact with delta_eval = "reml"
# (per-candidate grid-profile delta_hat_jk under H0_jk), same ML-eBIC
# stopping. The marginal-exact fit is already in results/bench_full_exact/.
#
# Output: results/bench_full_exact/31_reml_anchor/n1000_sg{sg}_b{b}.rds
# Usage:  Rscript sim/bench_full/31_exact_reml_anchor.R [--cores 5]
# ==============================================================================

Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressPackageStartupMessages(library(parallel))
source("sim/bench_full/00_config.R")
source("R/CBF_LMM_exact.R")

OUT_DIR <- "results/bench_full_exact/31_reml_anchor"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(flag, default) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) as.integer(args[i + 1L]) else default
}
N_CORES <- arg_val("--cores", 5L)

one_rep <- function(sg, b) {
  ck <- file.path(OUT_DIR, sprintf("n1000_sg%.1f_b%02d.rds", sg, b))
  if (file.exists(ck) && file.size(ck) > 0) return("cached")
  seed <- 20260425L + 1000L * (b - 1L) + 1000L + as.integer(round(1000 * sg))
  d <- gen_dataset(n = 1000L, m = ANCHOR$m, rho = ANCHOR$rho,
                   K_true = ANCHOR$K_true, beta_true = ANCHOR$beta_true,
                   sigma_g2 = sg, block_size = ANCHOR$block_size, seed = seed)
  t0  <- Sys.time()
  fit <- CBF_LMM_stepwise_exact(d$y, d$X, tau2 = ANCHOR$tau2,
                                  K_max = ANCHOR$K_max,
                                  delta_eval = "reml",
                                  n_nodes = ANCHOR$N_delta)
  el <- as.numeric(Sys.time() - t0, units = "secs")
  tp <- length(intersect(fit$indices, d$truth))
  saveRDS(list(sg = sg, b = b, indices = fit$indices, K_hat = fit$K_hat,
               tp = tp, recall = tp / length(d$truth),
               precision = if (fit$K_hat > 0) tp / fit$K_hat else 0,
               elapsed = el), ck)
  "done"
}

grid <- expand.grid(sg = c(0, 0.5), b = seq_len(100L))
message(sprintf("REML-exact anchor: %d replicates on %d cores",
                nrow(grid), N_CORES))
invisible(mclapply(seq_len(nrow(grid)), function(i)
  tryCatch(one_rep(grid$sg[i], grid$b[i]),
           error = function(e) message(sprintf("FAIL sg%.1f b%d: %s",
                                               grid$sg[i], grid$b[i],
                                               conditionMessage(e)))),
  mc.cores = N_CORES, mc.preschedule = FALSE))

fs <- list.files(OUT_DIR, pattern = "rds$", full.names = TRUE)
res <- do.call(rbind, lapply(fs, function(f) {
  r <- readRDS(f)
  f1 <- if (r$precision + r$recall > 0)
          2 * r$precision * r$recall / (r$precision + r$recall) else 0
  data.frame(sg = r$sg, K_hat = r$K_hat, recall = r$recall,
             precision = r$precision, f1 = f1)
}))
print(aggregate(cbind(K_hat, recall, precision, f1) ~ sg, res, mean), digits = 3)
message("done.")
