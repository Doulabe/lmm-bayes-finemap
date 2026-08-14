# ==============================================================================
# 22_gamma_sensitivity.R
#
# Sensitivity of the eBIC stopping rule to the penalty parameter gamma
# (Supplementary Note S7.4). The criterion is
#     eBIC(k) = -2 ell_k + k (log n + 2 gamma log m),
# with the default gamma = 1 (Chen & Chen, 2008). This script maps the effect of
# gamma in {0, 0.5, 1} on K_hat, F1, precision and recall.
#
# Setting: anchor (n=1000, m=1000, rho=0.95, K_true=5, beta=(0.8,0.4,0.4,0.2,
# 0.2)) crossed with sigma_g2 in {0, 0.5}. B=12 replicates per cell. tau^2=0.04
# fixed; delta estimated by profile-likelihood REML (single eigendecomposition
# per step).
#
# Output: results/bench_full/22_gamma_sensitivity/gamma_summary.rds (+ table)
#
# Usage:
#   Rscript sim/bench_full/22_gamma_sensitivity.R
# ==============================================================================

Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressPackageStartupMessages(library(parallel))
source("sim/bench_full/00_config.R")

OUT_DIR <- "results/bench_full/22_gamma_sensitivity"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

B        <- 12L
N_CORES  <- max(1L, min(4L, detectCores() - 1L))
GAMMA_GRID <- c(0.0, 0.5, 1.0)
SG_GRID    <- c(0.0, 0.5)

metrics <- function(sel, truth, K_star = 5L) {
  tp <- length(intersect(sel, truth)); K <- length(sel)
  rec  <- tp / K_star
  prec <- if (K > 0) tp / K else 0
  f1   <- if (rec + prec > 0) 2 * rec * prec / (rec + prec) else 0
  c(K_hat = K, recall = rec, precision = prec, f1 = f1)
}

reml_delta_profile <- function(Uty, UtF, s, n, p) {
  negll <- function(log_delta) {
    d <- 1 + exp(log_delta) * s; Dinv <- 1 / d
    FtDF <- crossprod(UtF * Dinv, UtF); FtDy <- crossprod(UtF * Dinv, Uty)
    yDy  <- sum(Uty^2 * Dinv)
    ch <- tryCatch(chol(FtDF), error = function(e) NULL); if (is.null(ch)) return(1e10)
    sol <- backsolve(ch, forwardsolve(t(ch), FtDy))
    rss <- as.numeric(yDy - crossprod(FtDy, sol)); if (!is.finite(rss) || rss <= 0) return(1e10)
    0.5 * (sum(log(d)) + 2 * sum(log(diag(ch))) + (n - p) * log(rss))
  }
  o <- tryCatch(optimize(negll, c(log(1e-4), log(1e4))), error = function(e) NULL)
  if (is.null(o)) 0.04 else max(exp(o$minimum), 1e-6)
}

# Plug-in REML stepwise with eBIC penalty parametrised by gamma.
reml_stepwise_gamma <- function(y, X, gamma, K_max = 8L, tau2 = 0.04) {
  n <- nrow(X); m <- ncol(X)
  selected <- integer(0); remaining <- seq_len(m); F_k <- matrix(1, n, 1L)
  ebic_path <- c(Inf, rep(NA, K_max))
  for (k_step in seq_len(K_max)) {
    if (length(remaining) == 0L) break
    rem <- X[, remaining, drop = FALSE]; K <- tcrossprod(rem) / max(ncol(rem) - 1L, 1L)
    ee <- eigen(K, symmetric = TRUE); s <- pmax(ee$values, 0); U <- ee$vectors
    Uty0 <- crossprod(U, y); UtF0 <- crossprod(U, F_k); UtX0 <- crossprod(U, rem)
    delta_hat <- reml_delta_profile(Uty0, UtF0, s, n, ncol(F_k))
    w <- 1 / sqrt(1 + delta_hat * s)
    Uty <- Uty0 * w; UtFk <- UtF0 * w; UtX <- UtX0 * w
    A <- crossprod(UtFk); Ainv <- solve(A); beta <- Ainv %*% crossprod(UtFk, Uty)
    ty <- Uty - UtFk %*% beta; projX <- UtFk %*% (Ainv %*% crossprod(UtFk, UtX)); tXj <- UtX - projX
    D_jk <- colSums(tXj^2); u_j <- as.numeric(crossprod(tXj, ty))
    rss0 <- as.numeric(crossprod(ty)); df <- n - ncol(F_k); Q2_j <- u_j^2 / (D_jk * rss0 / df)
    v_j <- tau2 * D_jk; lb <- -0.5 * log1p(v_j) - 0.5 * df * log1p(-v_j * Q2_j / df / (1 + v_j))
    j_best <- which.max(lb); if (!is.finite(lb[j_best])) break
    sel_new <- c(selected, remaining[j_best]); F_new <- cbind(F_k, X[, remaining[j_best], drop = FALSE])
    proj2 <- solve(crossprod(F_new), crossprod(F_new, y)); rss1 <- as.numeric(crossprod(y - F_new %*% proj2))
    ll_new <- -0.5 * (n - ncol(F_new)) * log(rss1)
    ebic_new <- -2 * ll_new + length(sel_new) * (log(n) + 2 * gamma * log(m))
    if (k_step >= 2L && ebic_new > ebic_path[k_step]) break
    selected <- sel_new; F_k <- F_new; ebic_path[k_step + 1L] <- ebic_new; remaining <- remaining[-j_best]
  }
  selected
}

run_cell <- function(sg, g, m = 1000L, B = 12L) {
  res <- mclapply(seq_len(B), function(b) {
    tryCatch({
      sd <- 20260951L + b + round(10 * sg)
      d <- gen_dataset(n = 1000L, m = m, rho = 0.95, K_true = 5L,
                       beta_true = c(0.8, 0.4, 0.4, 0.2, 0.2),
                       sigma_g2 = sg, block_size = 10L, seed = sd)
      metrics(reml_stepwise_gamma(d$y, d$X, gamma = g, K_max = 8L), d$truth)
    }, error = function(e) NULL)
  }, mc.cores = N_CORES)
  res <- res[!vapply(res, is.null, logical(1))]
  c(colMeans(do.call(rbind, res), na.rm = TRUE), n_ok = length(res))
}

GA <- list()
for (sg in SG_GRID) for (g in GAMMA_GRID) {
  GA[[sprintf("sg%.1f_g%.1f", sg, g)]] <- run_cell(sg, g)
  saveRDS(GA, file.path(OUT_DIR, "gamma_summary.rds"))
}

cat("=== eBIC gamma sensitivity (REML-profile stepwise, anchor m=1000) ===\n\n")
cat(sprintf("  %-12s %8s %8s %8s %8s %6s\n", "cell", "K_hat", "F1", "prec", "recall", "n_ok"))
for (nm in names(GA)) {
  M <- GA[[nm]]
  cat(sprintf("  %-12s %8.2f %8.3f %8.3f %8.3f %6d\n",
              nm, M["K_hat"], M["f1"], M["precision"], M["recall"], as.integer(M["n_ok"])))
}
cat(sprintf("\nSaved: %s\n", file.path(OUT_DIR, "gamma_summary.rds")))
