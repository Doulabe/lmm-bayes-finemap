# ==============================================================================
# 24_shared_kernel_validation.R
#
# Validation of the shared-kernel eigendecomposition approximation
# (Supplementary Note S3). The derivation uses a candidate-specific leave-one-out
# random-effect kernel R_jk (background = all markers except the candidate j),
# but the implementation amortises a single shared step-level kernel across all
# candidates, since dropping one column perturbs the n x n sample covariance
# eigenstructure by O(1/m). This script quantifies the resulting bias on a small
# design by comparing, per candidate at step 1:
#   * exact : log BF using the leave-candidate-out kernel K_{-j} (one eigen per j)
#   * shared: log BF using the single shared kernel K (one eigen for all j)
# and reporting max |delta log BF|, the Spearman rank correlation of the two
# log-BF vectors, and the Jaccard overlap of the top-K* selected sets.
#
# Output: results/bench_full/24_shared_kernel/shared_kernel_validation.rds (+ table)
#
# Usage:
#   Rscript sim/bench_full/24_shared_kernel_validation.R
# ==============================================================================

Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressPackageStartupMessages(library(parallel))
source("sim/bench_full/00_config.R")

OUT_DIR <- "results/bench_full/24_shared_kernel"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

TAU2 <- 0.04
B    <- 5L
N_CORES <- max(1L, min(4L, detectCores() - 1L))

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

# Step-1 conditional log BF for a set of candidate columns Xc, given a kernel K
# and conditioning matrix F_k. Vectorised over the columns of Xc.
logbf_given_kernel <- function(y, Xc, F_k, K, tau2 = TAU2) {
  n <- nrow(Xc)
  ee <- eigen(K, symmetric = TRUE); s <- pmax(ee$values, 0); U <- ee$vectors
  Uty0 <- crossprod(U, y); UtF0 <- crossprod(U, F_k); UtX0 <- crossprod(U, Xc)
  delta <- reml_delta_profile(Uty0, UtF0, s, n, ncol(F_k))
  w <- 1 / sqrt(1 + delta * s)
  Uty <- Uty0 * w; UtF <- UtF0 * w; UtX <- UtX0 * w
  A <- crossprod(UtF); Ainv <- solve(A); beta <- Ainv %*% crossprod(UtF, Uty)
  ty <- Uty - UtF %*% beta
  projX <- UtF %*% (Ainv %*% crossprod(UtF, UtX)); tX <- UtX - projX
  D <- colSums(tX^2); u <- as.numeric(crossprod(tX, ty))
  rss0 <- as.numeric(crossprod(ty)); df <- n - ncol(F_k); Q2 <- u^2 / (D * rss0 / df)
  v <- tau2 * D
  as.numeric(-0.5 * log1p(v) - 0.5 * df * log1p(-v * Q2 / df / (1 + v)))
}

run_rep <- function(n, m, rho, sg, seed) {
  d <- gen_dataset(n = n, m = m, rho = rho, K_true = 4L,
                   beta_true = c(0.8, 0.4, 0.4, 0.2),
                   sigma_g2 = sg, block_size = 10L, seed = seed)
  X <- d$X; y <- d$y; F_k <- matrix(1, n, 1L)

  # Shared kernel: one eigendecomposition from all markers, used for all candidates.
  K_shared <- tcrossprod(X) / max(m - 1L, 1L)
  lbf_shared <- logbf_given_kernel(y, X, F_k, K_shared, TAU2)

  # Exact kernel: leave-candidate-out background, one eigendecomposition per candidate.
  lbf_exact <- numeric(m)
  for (j in seq_len(m)) {
    Kj <- tcrossprod(X[, -j, drop = FALSE]) / max(m - 2L, 1L)
    lbf_exact[j] <- logbf_given_kernel(y, X[, j, drop = FALSE], F_k, Kj, TAU2)
  }

  ok <- is.finite(lbf_shared) & is.finite(lbf_exact)
  Kstar <- 4L
  topk_shared <- order(lbf_shared, decreasing = TRUE)[seq_len(Kstar)]
  topk_exact  <- order(lbf_exact,  decreasing = TRUE)[seq_len(Kstar)]
  jac <- length(intersect(topk_shared, topk_exact)) /
         length(union(topk_shared, topk_exact))
  c(max_abs_dlogbf = max(abs(lbf_shared[ok] - lbf_exact[ok])),
    spearman       = suppressWarnings(cor(lbf_shared[ok], lbf_exact[ok], method = "spearman")),
    jaccard_topK   = jac)
}

CELLS <- list(
  list(tag = "anchor, m=200",         n = 300L, m = 200L,  rho = 0.95, sg = 0.0),
  list(tag = "anchor, m=1000",        n = 300L, m = 1000L, rho = 0.95, sg = 0.0),
  list(tag = "polygenic, m=200",      n = 300L, m = 200L,  rho = 0.95, sg = 0.5),
  list(tag = "high LD, m=200",        n = 300L, m = 200L,  rho = 0.98, sg = 0.0)
)

RES <- list()
for (cc in CELLS) {
  reps <- mclapply(seq_len(B), function(b)
    tryCatch(run_rep(cc$n, cc$m, cc$rho, cc$sg, seed = 20260990L + b + round(10 * cc$sg)),
             error = function(e) NULL), mc.cores = N_CORES)
  reps <- reps[!vapply(reps, is.null, logical(1))]
  M <- do.call(rbind, reps)
  RES[[cc$tag]] <- colMeans(M, na.rm = TRUE)
  saveRDS(RES, file.path(OUT_DIR, "shared_kernel_validation.rds"))
}

cat("=== Shared-kernel vs exact leave-candidate-out kernel (step 1) ===\n")
cat(sprintf("  %-26s %14s %12s %12s\n", "cell", "max|dlogBF|", "Spearman", "Jaccard@K*"))
for (nm in names(RES)) {
  r <- RES[[nm]]
  cat(sprintf("  %-26s %14.4f %12.4f %12.3f\n",
              nm, r["max_abs_dlogbf"], r["spearman"], r["jaccard_topK"]))
}
cat(sprintf("\nSaved: %s\n", file.path(OUT_DIR, "shared_kernel_validation.rds")))

## ---- Step-2 validation (high LD): condition on the top-1 step-1 SNP -------
## The more conditioned regime, where a correlated neighbour of the candidate
## may already sit in F_k. Exact kernel excludes the selected SNP and the
## candidate; shared kernel excludes only the selected SNP.
run_rep_step2 <- function(seed, n = 300L, m = 200L, rho = 0.98, sg = 0.0) {
  d <- gen_dataset(n = n, m = m, rho = rho, K_true = 4L,
                   beta_true = c(0.8, 0.4, 0.4, 0.2),
                   sigma_g2 = sg, block_size = 10L, seed = seed)
  X <- d$X; y <- d$y
  lbf1 <- logbf_given_kernel(y, X, matrix(1, n, 1L), tcrossprod(X) / max(m - 1L, 1L))
  sel  <- which.max(lbf1)
  Fk   <- cbind(1, X[, sel]); rem <- setdiff(seq_len(m), sel)
  lbf_sh <- logbf_given_kernel(y, X[, rem, drop = FALSE], Fk,
                               tcrossprod(X[, -sel, drop = FALSE]) / max(m - 2L, 1L))
  lbf_ex <- numeric(length(rem))
  for (i in seq_along(rem)) {
    j <- rem[i]
    Kj <- tcrossprod(X[, -c(sel, j), drop = FALSE]) / max(m - 3L, 1L)
    lbf_ex[i] <- logbf_given_kernel(y, X[, j, drop = FALSE], Fk, Kj)
  }
  ok <- is.finite(lbf_sh) & is.finite(lbf_ex); Kstar <- 3L
  ts <- rem[order(lbf_sh, decreasing = TRUE)[seq_len(Kstar)]]
  te <- rem[order(lbf_ex, decreasing = TRUE)[seq_len(Kstar)]]
  c(max_abs_dlogbf = max(abs(lbf_sh[ok] - lbf_ex[ok])),
    spearman = suppressWarnings(cor(lbf_sh[ok], lbf_ex[ok], method = "spearman")),
    jaccard_topK = length(intersect(ts, te)) / length(union(ts, te)))
}
reps2 <- mclapply(seq_len(B), function(b)
  tryCatch(run_rep_step2(20260995L + b), error = function(e) NULL), mc.cores = N_CORES)
reps2 <- reps2[!vapply(reps2, is.null, logical(1))]
M2 <- colMeans(do.call(rbind, reps2), na.rm = TRUE)
saveRDS(M2, file.path(OUT_DIR, "shared_kernel_step2_highLD.rds"))
cat(sprintf("\n=== Step 2, high LD (rho=0.98), 1 SNP conditioned ===\n  max|dlogBF|=%.2f  Spearman=%.4f  Jaccard@K*=%.3f\n",
            M2["max_abs_dlogbf"], M2["spearman"], M2["jaccard_topK"]))
