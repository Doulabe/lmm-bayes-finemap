# ==============================================================================
# 21_eb_tau2_sensitivity.R
#
# Empirical-Bayes estimation of the slab variance tau^2 (Supplementary Note S5).
#
# The conditional Bayes factor at a fixed variance ratio delta is a smooth scalar
# function of tau^2, so an empirical-Bayes slab is available in closed form: at
# each selection step tau^2 is set to the maximiser of log BF for the leading
# candidate by a single 1-D Brent search (EB-tau2). This script compares EB-tau2
# with the fixed default tau^2 = 0.04 and records the distribution of the
# estimated tau^2 along the greedy path.
#
# Setting: anchor (n=1000, rho in {0.95, 0.98}, K_true=5, beta=(0.8,0.4,0.4,0.2,
# 0.2)) crossed with sigma_g2 in {0, 0.5}; m=1000 keeps the per-step
# eigendecomposition affordable for the sweep (the conclusion is m-robust).
# B=10 replicates per cell. Variance ratio delta is estimated by a profile-
# likelihood REML using a single eigendecomposition per step.
#
# Output: results/bench_full/21_eb_tau2/eb_tau2_summary.rds (+ printed table)
#
# Usage:
#   Rscript sim/bench_full/21_eb_tau2_sensitivity.R
# ==============================================================================

Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressPackageStartupMessages(library(parallel))
source("sim/bench_full/00_config.R")

OUT_DIR <- "results/bench_full/21_eb_tau2"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

B          <- 10L
N_CORES    <- max(1L, min(4L, detectCores() - 1L))
M_SWEEP    <- 1000L

metrics <- function(sel, truth, K_star = 5L) {
  tp <- length(intersect(sel, truth)); K <- length(sel)
  rec  <- tp / K_star
  prec <- if (K > 0) tp / K else 0
  f1   <- if (rec + prec > 0) 2 * rec * prec / (rec + prec) else 0
  c(K_hat = K, recall = rec, precision = prec, f1 = f1)
}

# Fast REML estimate of delta from a single eigendecomposition (profile REML).
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

# Stepwise scan; tau2_mode in {"fixed","eb"}. Returns selected indices and the
# per-step EB-tau2 path.
reml_stepwise <- function(y, X, K_max = 8L, tau2_mode = "fixed", tau2_fixed = 0.04) {
  n <- nrow(X); m <- ncol(X)
  selected <- integer(0); remaining <- seq_len(m); F_k <- matrix(1, n, 1L)
  eb_tau2_path <- numeric(0); ebic_path <- c(Inf, rep(NA, K_max))
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
    logBF_tau <- function(tau2, j) { v <- tau2 * D_jk[j]; -0.5 * log1p(v) - 0.5 * df * log1p(-v * Q2_j[j] / df / (1 + v)) }
    if (tau2_mode == "eb") {
      v0 <- 0.04 * D_jk; lb0 <- -0.5 * log1p(v0) - 0.5 * df * log1p(-v0 * Q2_j / df / (1 + v0)); j_best <- which.max(lb0)
      opt <- tryCatch(optimize(function(t) logBF_tau(t, j_best), c(1e-3, 1.0), maximum = TRUE),
                      error = function(e) list(maximum = 0.04))
      tau2_use <- opt$maximum; eb_tau2_path <- c(eb_tau2_path, tau2_use)
    } else {
      tau2_use <- tau2_fixed; v0 <- tau2_use * D_jk
      lb0 <- -0.5 * log1p(v0) - 0.5 * df * log1p(-v0 * Q2_j / df / (1 + v0)); j_best <- which.max(lb0)
    }
    if (!is.finite(lb0[j_best])) break
    sel_new <- c(selected, remaining[j_best]); F_new <- cbind(F_k, X[, remaining[j_best], drop = FALSE])
    proj2 <- solve(crossprod(F_new), crossprod(F_new, y)); rss1 <- as.numeric(crossprod(y - F_new %*% proj2))
    ll_new <- -0.5 * (n - ncol(F_new)) * log(rss1); ebic_new <- -2 * ll_new + length(sel_new) * (log(n) + 2 * log(m))
    if (k_step >= 2L && ebic_new > ebic_path[k_step]) break
    selected <- sel_new; F_k <- F_new; ebic_path[k_step + 1L] <- ebic_new; remaining <- remaining[-j_best]
  }
  list(indices = selected, eb_tau2 = eb_tau2_path)
}

run_cell <- function(n, rho, sg, m = M_SWEEP, B = 10L, seed0 = 20260901L) {
  res <- mclapply(seq_len(B), function(b) {
    tryCatch({
      sd <- seed0 + b + as.integer(n) + round(100 * rho) + round(10 * sg)
      d <- gen_dataset(n = n, m = m, rho = rho, K_true = 5L,
                       beta_true = c(0.8, 0.4, 0.4, 0.2, 0.2),
                       sigma_g2 = sg, block_size = 10L, seed = sd)
      eb <- reml_stepwise(d$y, d$X, K_max = 8L, tau2_mode = "eb")
      fx <- reml_stepwise(d$y, d$X, K_max = 8L, tau2_mode = "fixed", tau2_fixed = 0.04)
      list(eb_m = metrics(eb$indices, d$truth), fx_m = metrics(fx$indices, d$truth),
           eb_tau2 = eb$eb_tau2)
    }, error = function(e) NULL)
  }, mc.cores = N_CORES)
  res <- res[!vapply(res, is.null, logical(1))]
  eb_M <- do.call(rbind, lapply(res, `[[`, "eb_m"))
  fx_M <- do.call(rbind, lapply(res, `[[`, "fx_m"))
  all_tau2 <- unlist(lapply(res, `[[`, "eb_tau2"))
  list(eb = colMeans(eb_M, na.rm = TRUE), fx = colMeans(fx_M, na.rm = TRUE),
       n_ok = length(res), tau2_q = quantile(all_tau2, c(0.1, 0.5, 0.9), na.rm = TRUE))
}

CELLS <- list(c(1000, 0.95, 0.0, "anchor sg=0"),
              c(1000, 0.95, 0.5, "anchor sg=0.5"),
              c(1000, 0.98, 0.5, "high LD sg=0.5"))
EB <- list()
for (cfg in CELLS) {
  EB[[cfg[4]]] <- run_cell(as.integer(cfg[1]), as.numeric(cfg[2]), as.numeric(cfg[3]))
  saveRDS(EB, file.path(OUT_DIR, "eb_tau2_summary.rds"))
}

cat("=== EB-tau2 vs fixed tau2=0.04 (REML-profile stepwise, eBIC) ===\n\n")
for (nm in names(EB)) {
  r <- EB[[nm]]
  cat(sprintf("--- %s (n_ok=%d) ---\n", nm, r$n_ok))
  cat(sprintf("  EB-tau2 : F1=%.3f  K_hat=%.2f  prec=%.3f  rec=%.3f\n",
              r$eb["f1"], r$eb["K_hat"], r$eb["precision"], r$eb["recall"]))
  cat(sprintf("  fix 0.04: F1=%.3f  K_hat=%.2f  prec=%.3f  rec=%.3f\n",
              r$fx["f1"], r$fx["K_hat"], r$fx["precision"], r$fx["recall"]))
  cat(sprintf("  EB-tau2 est: q10=%.3f  median=%.3f  q90=%.3f\n\n",
              r$tau2_q[1], r$tau2_q[2], r$tau2_q[3]))
}
cat(sprintf("Saved: %s\n", file.path(OUT_DIR, "eb_tau2_summary.rds")))
