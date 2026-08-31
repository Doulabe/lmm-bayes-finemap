# ==============================================================================
# 39_exact_tau2_gamma_anchor.R
# Exact-implementation sensitivity sweeps on the anchor configuration
# (n=1000, m=5000, rho=0.95, medium architecture; same seeds as axis 01):
#   * slab variance tau2 in {0.01, 0.10, 0.25}  (0.04 = main campaign, reused)
#   * eBIC penalty gamma in {0.5, 1.5}          (1    = main campaign, reused)
# 100 replicates per (setting, regime). Checkpointed, resumable.
# Output: results/bench_full_exact/39_tau2_gamma/ + summary CSV
# Usage:  Rscript sim/bench_full/39_exact_tau2_gamma_anchor.R [--cores 5]
# ==============================================================================

Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressPackageStartupMessages(library(parallel))
source("sim/bench_full/00_config.R")
source("R/CBF_LMM_exact.R")

OUT_DIR <- "results/bench_full_exact/39_tau2_gamma"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(flag, default) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) as.integer(args[i + 1L]) else default
}
N_CORES <- arg_val("--cores", 5L)

one_rep <- function(kind, val, sg, b) {
  ck <- file.path(OUT_DIR, sprintf("%s%.2f_sg%.1f_b%03d.rds", kind, val, sg, b))
  if (file.exists(ck) && file.size(ck) > 0) return("cached")
  seed <- 20260425L + 1000L * (b - 1L) + 1000L + as.integer(round(1000 * sg))
  d <- gen_dataset(n = 1000L, m = ANCHOR$m, rho = ANCHOR$rho,
                   K_true = ANCHOR$K_true, beta_true = ANCHOR$beta_true,
                   sigma_g2 = sg, block_size = ANCHOR$block_size, seed = seed)
  fit <- if (kind == "tau")
    CBF_LMM_stepwise_exact(d$y, d$X, tau2 = val, K_max = ANCHOR$K_max,
                             n_nodes = ANCHOR$N_delta)
  else
    CBF_LMM_stepwise_exact(d$y, d$X, tau2 = ANCHOR$tau2, gamma = val,
                             K_max = ANCHOR$K_max, n_nodes = ANCHOR$N_delta)
  tp <- length(intersect(fit$indices, d$truth))
  saveRDS(list(kind = kind, val = val, sg = sg, b = b,
               K_hat = fit$K_hat, tp = tp,
               recall = tp / length(d$truth),
               precision = if (fit$K_hat > 0) tp / fit$K_hat else 0), ck)
  "done"
}

grid <- rbind(expand.grid(kind = "tau",   val = c(0.01, 0.10, 0.25),
                          sg = c(0, 0.5), b = 1:100,
                          stringsAsFactors = FALSE),
              expand.grid(kind = "gamma", val = c(0.5, 1.5),
                          sg = c(0, 0.5), b = 1:100,
                          stringsAsFactors = FALSE))
message(sprintf("exact tau2/gamma anchor: %d replicates on %d cores",
                nrow(grid), N_CORES))
invisible(mclapply(seq_len(nrow(grid)), function(i)
  tryCatch(one_rep(grid$kind[i], grid$val[i], grid$sg[i], grid$b[i]),
           error = function(e) message(sprintf("FAIL %s%.2f sg%.1f b%d: %s",
             grid$kind[i], grid$val[i], grid$sg[i], grid$b[i],
             conditionMessage(e)))),
  mc.cores = N_CORES, mc.preschedule = FALSE))

fs <- list.files(OUT_DIR, pattern = "rds$", full.names = TRUE)
res <- do.call(rbind, lapply(fs, function(f) {
  r <- readRDS(f)
  f1 <- if (r$precision + r$recall > 0)
          2 * r$precision * r$recall / (r$precision + r$recall) else 0
  data.frame(kind = r$kind, val = r$val, sg = r$sg, K_hat = r$K_hat,
             recall = r$recall, precision = r$precision, f1 = f1)
}))
summ <- aggregate(cbind(K_hat, recall, precision, f1) ~ kind + val + sg,
                  res, mean)
print(summ, digits = 3, row.names = FALSE)
write.csv(summ, file.path(OUT_DIR, "tau2_gamma_summary.csv"),
          row.names = FALSE)
message("done.")
