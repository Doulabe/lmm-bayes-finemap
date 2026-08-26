# ==============================================================================
# examples/quickstart.R
#
# Minimal demo of CBF-LMM (conditional Bayes factor + marginalized delta +
# eBIC stopping).  Simulates a small dataset (n = 500, m = 300, K_true = 3),
# runs the primary procedure, then two optional runs: the plug-in REML
# sensitivity evaluation and the exploratory joint Schur-complement score.
#
# Run from the repository root:
#   Rscript examples/quickstart.R
# ==============================================================================

source("R/LMM_core.R")
source("R/LMM_stepwise_fast.R")
source("R/LMM_reml_bf.R")
source("R/CBF_LMM.R")

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

cat("=== CBF-LMM (primary: marginalized delta + eBIC) ===\n")
res <- CBF_LMM_stepwise(y, X, tau2 = 0.04, K_max = 10L, n_nodes = 15L)
cat("Selected indices:", res$indices, "\n")
cat("Truth:           ", truth, "\n")
cat("K_hat:", res$K_hat, " q_at_step:",
     round(res$q_at_step, 3), "\n\n")

cat("=== Sensitivity: plug-in REML evaluation of delta ===\n")
res_reml <- CBF_LMM_stepwise(y, X, tau2 = 0.04, K_max = 10L,
                             delta_eval = "reml")
cat("Selected indices:", res_reml$indices, "\n")
cat("K_hat:", res_reml$K_hat, "\n\n")

cat("=== Exploratory: joint Schur-complement score (future work) ===\n")
res_js <- JS_L_LMM_stepwise_fast(y, X, tau2 = 0.04, K_max = 10L,
                                      criterion = "eBIC", theta = 0.99,
                                      n_nodes = 15L)
cat("Selected indices:", res_js$indices, "\n")
cat("K_hat:", res_js$K_hat, "\n")

cat("\nDone. Under the rho = 0.95 within-block LD of this demo, a selected index\nmay be an immediate LD neighbour of the true causal (e.g. 119 for 120):\nthis is expected fine-mapping behaviour at n = 500, not an error.\n")
