# ==============================================================================
# examples/quickstart.R
#
# Minimal 30-line demo of the Bayesian stepwise marginal-likelihood framework.
# Simulates a small dataset (n = 500, m = 300, K_true = 3) and runs both the
# marginal and joint variants with eBIC stopping.
#
# Run from the repository root:
#   Rscript examples/quickstart.R
# ==============================================================================

source("R/LMM_core.R")
source("R/LMM_stepwise_fast.R")

set.seed(42)
n <- 500L; m <- 300L
truth <- c(50L, 120L, 200L)
beta_true <- c(0.8, 0.4, 0.4)

# Block-AR(1) correlated X
block_size <- 10L
X <- matrix(rnorm(n * m), n, m)
ar1_block <- function(p, r) outer(seq_len(p), seq_len(p),
                                       function(a, b) r^abs(a - b))
L <- chol(ar1_block(block_size, 0.95) + 1e-8 * diag(block_size))
for (b in seq_len(m %/% block_size)) {
  cols <- ((b - 1) * block_size + 1):(b * block_size)
  X[, cols] <- X[, cols] %*% L
}
X <- scale(X)

# Outcome with three causal effects plus iid noise
y <- as.numeric(X[, truth] %*% beta_true) + rnorm(n)

cat("=== Marginal variant (MS_L_eBIC) ===\n")
res_ms <- MS_L_LMM_stepwise_fast(y, X, tau2 = 0.04, K_max = 10L,
                                      criterion = "eBIC", theta = 0.99,
                                      n_nodes = 15L)
cat("Selected indices:", res_ms$indices, "\n")
cat("Truth:           ", truth, "\n")
cat("K_hat:", res_ms$K_hat, " q_at_step:",
     round(res_ms$q_at_step, 3), "\n\n")

cat("=== Joint variant (JS_L_eBIC) ===\n")
res_js <- JS_L_LMM_stepwise_fast(y, X, tau2 = 0.04, K_max = 10L,
                                      criterion = "eBIC", theta = 0.99,
                                      n_nodes = 15L)
cat("Selected indices:", res_js$indices, "\n")
cat("K_hat:", res_js$K_hat, "\n")

cat("\nDone. Under the rho = 0.95 within-block LD of this demo, a selected index\nmay be an immediate LD neighbour of the true causal (e.g. 119 for 120):\nthis is expected fine-mapping behaviour at n = 500, not an error.\n")
