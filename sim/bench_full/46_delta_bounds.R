# ==============================================================================
# 46_delta_bounds.R
# Sensitivity of the marginalized Bayes factor to the delta truncation bounds:
# primary grid (15 nodes on [1e-4, 1e3]) versus a wider, denser grid
# (19 nodes on [1e-5, 1e4], comparable node spacing). Anchor configuration,
# both polygenic regimes, campaign seed formula. Decisions compared per
# replicate.
# Output: results/bench_full_exact/46_delta_bounds/
# Usage:  Rscript sim/bench_full/46_delta_bounds.R [--cores 2] [--B 100]
# ==============================================================================

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressPackageStartupMessages(library(parallel))
source("sim/bench_full/00_config.R")
source("R/CBF_LMM_exact.R")

OUT_DIR <- "results/bench_full_exact/46_delta_bounds"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(flag, default) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) as.integer(args[i + 1L]) else default
}
N_CORES <- arg_val("--cores", 2L)
B       <- arg_val("--B", 100L)

REGIMES <- c(0, 0.5)
jobs <- expand.grid(sg = REGIMES, b = seq_len(B))

one_job <- function(i) {
  sg <- jobs$sg[i]; b <- jobs$b[i]
  ckpt <- file.path(OUT_DIR, sprintf("sg%.1f_b%03d.rds", sg, b))
  if (file.exists(ckpt) && file.size(ckpt) > 0) return("cached")
  pc <- ANCHOR; pc$sigma_g2 <- sg
  seed <- 20260101L + 100L * b + as.integer(round(10 * sg))
  d <- gen_dataset(n = pc$n, m = pc$m, rho = pc$rho,
                   K_true = pc$K_true, beta_true = pc$beta,
                   sigma_g2 = sg, seed = seed)
  fit0 <- CBF_LMM_stepwise_exact(d$y, d$X, tau2 = pc$tau2, K_max = pc$K_max,
                                   n_nodes = 15L,
                                   delta_lo = 1e-4, delta_hi = 1e3)
  fitW <- CBF_LMM_stepwise_exact(d$y, d$X, tau2 = pc$tau2, K_max = pc$K_max,
                                   n_nodes = 19L,
                                   delta_lo = 1e-5, delta_hi = 1e4)
  out <- data.frame(sigma_g2 = sg, rep = b,
                    K_primary = length(fit0$indices),
                    K_wide    = length(fitW$indices),
                    identical = identical(sort(fit0$indices),
                                          sort(fitW$indices)))
  saveRDS(list(row = out, sel_primary = fit0$indices,
               sel_wide = fitW$indices), ckpt)
  "done"
}

message(sprintf("46_delta_bounds: %d jobs on %d cores", nrow(jobs), N_CORES))
invisible(mclapply(seq_len(nrow(jobs)), function(i)
  tryCatch(one_job(i), error = function(e)
    message("FAIL job ", i, ": ", conditionMessage(e))),
  mc.cores = N_CORES, mc.preschedule = FALSE))

fs <- list.files(OUT_DIR, pattern = "^sg.*rds$", full.names = TRUE)
res <- do.call(rbind, lapply(fs, function(f) readRDS(f)$row))
agg <- aggregate(cbind(identical, K_primary, K_wide) ~ sigma_g2, res, mean)
print(agg, digits = 4, row.names = FALSE)
write.csv(agg, file.path(OUT_DIR, "delta_bounds_summary.csv"),
          row.names = FALSE)
message("done 46.")
