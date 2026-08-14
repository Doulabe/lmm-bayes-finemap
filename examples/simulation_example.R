# ==============================================================================
# examples/simulation_example.R
#
# Single anchor-cell simulation matching the benchmark configuration
# (n = 1000, m = 5000, K_true = 5, block-AR(1) rho = 0.95) with optional
# polygenic background.  Compares MS_L_eBIC, JS_L_eBIC, and JointPosterior
# stopping variants.
#
# Run from the repository root:
#   Rscript examples/simulation_example.R
# ==============================================================================

source("R/LMM_core.R")
source("R/LMM_stepwise_fast.R")
source("sim/bench_full/00_config.R")    # provides gen_dataset() and ANCHOR

# Choose the polygenic regime: 0 for homogeneous noise, 0.5 for polygenic
SIGMA_G2 <- 0.5

set.seed(20260425L)
d <- gen_dataset(n = ANCHOR$n, m = ANCHOR$m, rho = ANCHOR$rho,
                   K_true = ANCHOR$K_true, beta_true = ANCHOR$beta_true,
                   sigma_g2 = SIGMA_G2, block_size = ANCHOR$block_size,
                   seed = 20260425L)
cat(sprintf("Generated: n=%d, m=%d, K_true=%d, sigma_g^2=%.1f\n",
              nrow(d$X), ncol(d$X), length(d$truth), SIGMA_G2))
cat("True causal indices:", d$truth, "\n\n")

#  Marginal variant with eBIC stopping
cat("=== MS_L_eBIC ===\n")
t0 <- Sys.time()
res <- MS_L_LMM_stepwise_fast(d$y, d$X, tau2 = 0.04,
                                  K_max = 10L, criterion = "eBIC",
                                  theta = 0.99, n_nodes = 15L)
elapsed <- as.numeric(Sys.time() - t0, units = "secs")
tp <- length(intersect(res$indices, d$truth))
cat(sprintf("  K_hat = %d, true positives = %d / %d  (recall = %.2f)\n",
              res$K_hat, tp, length(d$truth), tp / length(d$truth)))
cat(sprintf("  elapsed: %.1f s\n", elapsed))
cat("  q_at_step:", round(res$q_at_step, 3), "\n\n")

#  Joint variant with eBIC
cat("=== JS_L_eBIC ===\n")
t0 <- Sys.time()
res <- JS_L_LMM_stepwise_fast(d$y, d$X, tau2 = 0.04,
                                  K_max = 10L, criterion = "eBIC",
                                  theta = 0.99, n_nodes = 15L)
elapsed <- as.numeric(Sys.time() - t0, units = "secs")
tp <- length(intersect(res$indices, d$truth))
cat(sprintf("  K_hat = %d, true positives = %d / %d  (recall = %.2f)\n",
              res$K_hat, tp, length(d$truth), tp / length(d$truth)))
cat(sprintf("  elapsed: %.1f s\n", elapsed))
cat("  q_at_step:", round(res$q_at_step, 3), "\n\n")

#  Joint variant with JointPosterior (more recall-oriented)
cat("=== JS_L_JP99 ===\n")
t0 <- Sys.time()
res <- JS_L_LMM_stepwise_fast(d$y, d$X, tau2 = 0.04,
                                  K_max = 10L, criterion = "JointPosterior",
                                  theta = 0.99, n_nodes = 15L)
elapsed <- as.numeric(Sys.time() - t0, units = "secs")
tp <- length(intersect(res$indices, d$truth))
cat(sprintf("  K_hat = %d, true positives = %d / %d  (recall = %.2f)\n",
              res$K_hat, tp, length(d$truth), tp / length(d$truth)))
cat(sprintf("  elapsed: %.1f s\n", elapsed))
cat("  q_at_step:", round(res$q_at_step, 3), "\n\n")

cat("Done. Expect eBIC variants to be parsimonious (K_hat ~ 3-5);\n")
cat("the JP99 variant typically selects a few additional secondary hits.\n")
