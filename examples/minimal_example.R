# ==============================================================================
# minimal_example.R — self-contained CBF-LMM demonstration on simulated data
#
# Simulates a block-AR(1) region with 3 causal variants and a polygenic
# background, runs the exact CBF-LMM procedure (candidate-specific kernel
# K_jk + profile-ML eBIC stopping), and reports the selected variants and
# the per-step wall-clock time at two sample sizes to illustrate the O(n^3)
# scaling of the per-step eigendecomposition.
#
# Usage (from the repository root):
#   Rscript examples/minimal_example.R
# Runs in a few minutes on a laptop; no external data required.
# ==============================================================================

source("R/CBF_LMM_exact.R")

simulate_region <- function(n, m, rho = 0.95, block = 10L,
                            truth = c(1L, 11L, 21L),
                            beta = c(0.8, 0.5, 0.4),
                            sigma_g2 = 0.5, seed = 1L) {
  set.seed(seed)
  Z <- matrix(rnorm(n * m), n, m)
  S <- outer(seq_len(block), seq_len(block), function(a, b) rho^abs(a - b))
  L <- chol(S + 1e-8 * diag(block))
  for (bk in seq_len(m %/% block)) {
    cols <- ((bk - 1L) * block + 1L):(bk * block)
    Z[, cols] <- Z[, cols] %*% L
  }
  X <- scale(Z)
  nc <- setdiff(seq_len(m), truth)
  u <- as.numeric(X[, nc] %*% rnorm(length(nc), 0, sqrt(sigma_g2 / length(nc))))
  y <- as.numeric(X[, truth] %*% beta) + u + rnorm(n)
  list(X = X, y = y, truth = truth)
}

run_one <- function(n, m) {
  d <- simulate_region(n, m)
  t0 <- Sys.time()
  fit <- CBF_LMM_stepwise_exact(d$y, d$X, tau2 = 0.04, K_max = 10L,
                                  n_nodes = 15L)
  elapsed <- as.numeric(Sys.time() - t0, units = "secs")
  cat(sprintf("n = %4d, m = %4d | selected: %-18s | truth: %s | %5.1f s (%.2f s/step)\n",
              n, m, paste(sort(fit$indices), collapse = ","),
              paste(d$truth, collapse = ","),
              elapsed, elapsed / max(1L, length(fit$indices) + 1L)))
  invisible(elapsed)
}

cat("CBF-LMM minimal example (block-AR(1), K* = 3, polygenic background)\n\n")
t1 <- run_one(n =  500L, m = 2000L)
t2 <- run_one(n = 1000L, m = 2000L)
cat(sprintf("\nPer-step cost is dominated by one n x n eigendecomposition, O(n^3):\nobserved ratio %.1f for n doubling (theoretical 8, lower in practice\nbecause the candidate-level O(N_delta n^2 c_k) work scales as n^2).\n",
            t2 / t1))
