# ==============================================================================
# 28_path_kernel_validation.R
#
# Path-level validation of the shared-kernel approximation (extends Note S3).
# Script 24 compares shared vs leave-candidate-out kernels within a single
# step (ranking agreement). Because the selection is recursive
# (S_0 -> S_1 -> ...), a small step-1 difference could in principle change
# K_2 and propagate. This script runs the FULL stepwise procedure twice on
# the same data:
#   * shared: one step-level kernel K_k = X_Ck X_Ck' / c_k (the implementation)
#   * loo   : exact leave-candidate-out kernel K_jk^(-j) for every candidate
# with identical conditioning, REML-profile delta, tau2, and the incremental
# eBIC rule (accept iff 2 logBF > log n + 2 log m; first candidate always
# accepted, K_max = 10), and reports per replicate:
#   Jaccard(S_shared, S_loo), K_shared, K_loo, identical (0/1).
#
# Designs (small, where the exact computation is feasible):
#   cell A: n=300, m=100, K*=4, B=50   (~40 s/rep on 1 core)
#   cell B: n=500, m=200, K*=4, B=12   (~5 min/rep on 1 core)
# both crossed with sigma_g2 in {0, 0.5}. Roughly ~1 h total on 4 cores.
#
# Output: results/bench_full/28_path_kernel/ (per-replicate .rds checkpoints,
#         resumable) + path_kernel_summary.csv
#
# Usage:
#   Rscript sim/bench_full/28_path_kernel_validation.R [--cores 4]
# ==============================================================================

Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressPackageStartupMessages(library(parallel))
source("sim/bench_full/00_config.R")

OUT_DIR <- "results/bench_full/28_path_kernel"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

TAU2  <- 0.04
K_MAX <- 10L

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(flag, default) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) as.integer(args[i + 1L]) else default
}
N_CORES <- arg_val("--cores", max(1L, min(4L, detectCores() - 1L)))

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

# Conditional log BF for candidate columns Xc given kernel K and block F_k
# (same evaluation as script 24; vectorised over the columns of Xc).
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

# Full stepwise path under either kernel mode, with the incremental eBIC rule.
stepwise_path <- function(y, X, mode = c("shared", "loo")) {
  mode <- match.arg(mode)
  n <- nrow(X); m <- ncol(X)
  pen <- log(n) + 2 * log(m)          # incremental eBIC charge per variant
  selected <- integer(0)
  repeat {
    k_step <- length(selected) + 1L
    if (k_step > K_MAX) break
    C_k <- setdiff(seq_len(m), selected)
    if (!length(C_k)) break
    F_k <- cbind(matrix(1, n, 1L), X[, selected, drop = FALSE])
    if (mode == "shared") {
      K <- tcrossprod(X[, C_k, drop = FALSE]) / length(C_k)
      lbf <- logbf_given_kernel(y, X[, C_k, drop = FALSE], F_k, K)
    } else {
      lbf <- vapply(seq_along(C_k), function(i) {
        B_jk <- C_k[-i]
        Kj <- tcrossprod(X[, B_jk, drop = FALSE]) / length(B_jk)
        logbf_given_kernel(y, X[, C_k[i], drop = FALSE], F_k, Kj)
      }, numeric(1L))
    }
    j_best <- which.max(lbf)
    if (!is.finite(lbf[j_best])) break
    # incremental eBIC: first candidate always accepted (as implemented)
    if (k_step >= 2L && 2 * lbf[j_best] <= pen) break
    selected <- c(selected, C_k[j_best])
  }
  selected
}

one_rep <- function(cell, sg, b) {
  ck <- file.path(OUT_DIR, sprintf("%s_sg%.1f_b%03d.rds", cell$id, sg, b))
  if (file.exists(ck) && file.size(ck) > 0) return(invisible("cached"))
  d <- gen_dataset(n = cell$n, m = cell$m, rho = 0.95, K_true = 4L,
                   beta_true = c(0.8, 0.4, 0.4, 0.2), sigma_g2 = sg,
                   block_size = 10L,
                   seed = 20270280L + 10000L * match(cell$id, c("A", "B")) +
                          1000L * as.integer(round(10 * sg)) + b)
  S_sh  <- stepwise_path(d$y, d$X, "shared")
  S_loo <- stepwise_path(d$y, d$X, "loo")
  jac <- if (length(union(S_sh, S_loo)) == 0L) 1 else
           length(intersect(S_sh, S_loo)) / length(union(S_sh, S_loo))
  saveRDS(list(cell = cell$id, sg = sg, b = b,
               S_shared = S_sh, S_loo = S_loo,
               K_shared = length(S_sh), K_loo = length(S_loo),
               jaccard = jac,
               identical = as.integer(setequal(S_sh, S_loo))), ck)
  invisible("done")
}

CELLS <- list(list(id = "A", n = 300L, m = 100L, B = 50L),
              list(id = "B", n = 500L, m = 200L, B = 12L))

grid <- do.call(rbind, lapply(CELLS, function(cl)
  expand.grid(cell = cl$id, sg = c(0, 0.5), b = seq_len(cl$B),
              stringsAsFactors = FALSE)))

message(sprintf("path-kernel validation: %d replicates on %d cores",
                nrow(grid), N_CORES))
invisible(mclapply(seq_len(nrow(grid)), function(i) {
  row  <- grid[i, ]
  cell <- CELLS[[match(row$cell, c("A", "B"))]]
  tryCatch(one_rep(cell, row$sg, row$b),
           error = function(e) message(sprintf("rep %s/%s/%d failed: %s",
                                               row$cell, row$sg, row$b,
                                               conditionMessage(e))))
}, mc.cores = N_CORES, mc.preschedule = FALSE))

# ---- summary ----------------------------------------------------------------
files <- list.files(OUT_DIR, pattern = "^[AB]_sg.*\\.rds$", full.names = TRUE)
res <- do.call(rbind, lapply(files, function(f) {
  r <- readRDS(f)
  # prefix: is the shorter path the exact ordered prefix of the longer one?
  a <- r$S_shared; b <- r$S_loo; n1 <- min(length(a), length(b))
  data.frame(cell = r$cell, sg = r$sg, K_shared = r$K_shared,
             K_loo = r$K_loo, jaccard = r$jaccard, identical = r$identical,
             prefix = as.integer(n1 == 0L ||
                                 identical(a[seq_len(n1)], b[seq_len(n1)])))
}))
summ <- aggregate(cbind(jaccard, identical, prefix, K_shared, K_loo) ~
                    cell + sg, data = res, FUN = mean)
names(summ)[3:7] <- c("mean_jaccard", "prop_identical", "prop_prefix",
                      "mean_K_shared", "mean_K_loo")
write.csv(summ, file.path(OUT_DIR, "path_kernel_summary.csv"),
          row.names = FALSE)
print(summ)
