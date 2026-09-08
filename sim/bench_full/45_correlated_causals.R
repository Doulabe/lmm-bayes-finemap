# ==============================================================================
# 45_correlated_causals.R
# Same-block correlated-causal stress test (reviewer request): K*=3 causal
# variants inside ONE AR(1) block (positions 1, 3, 5 of block 1), so pairwise
# causal-causal |r| = rho^2 and rho^4. Anchor-like design otherwise.
# Exact CBF engine + external comparators; exact-identity metrics.
# Output: results/bench_full_exact/45_correlated_causals/
# Usage:  Rscript sim/bench_full/45_correlated_causals.R [--cores 5] [--B 100]
# ==============================================================================

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressPackageStartupMessages(library(parallel))
source("sim/bench_full/00_config.R")
source("R/CBF_LMM_exact.R")
for (p in c("BGLR", "hibayes", "rrBLUP", "susieR")) loadNamespace(p)

OUT_DIR <- "results/bench_full_exact/45_correlated_causals"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(flag, default) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) as.integer(args[i + 1L]) else default
}
N_CORES <- arg_val("--cores", 5L)
B       <- arg_val("--B", 100L)

N <- 1000L; M <- 5000L
TRUTH <- c(1L, 3L, 5L)                    # same block, pairwise r = rho^2, rho^4
BETA  <- c(0.8, 0.4, 0.4)                 # anchor-scale effects
RHOS  <- c(0.95, 0.98)
REGIMES <- c(0, 0.5)

gen_correlated <- function(n, m, rho, truth, beta, sigma_g2, seed,
                           block_size = 10L) {
  set.seed(seed)
  Z <- matrix(rnorm(n * m), n, m)
  ar1 <- function(p, r) outer(seq_len(p), seq_len(p), function(a, b) r^abs(a - b))
  L <- chol(ar1(block_size, rho) + 1e-8 * diag(block_size))
  for (bk in seq_len(m %/% block_size)) {
    cols <- ((bk - 1) * block_size + 1):(bk * block_size)
    Z[, cols] <- Z[, cols] %*% L
  }
  X <- scale(Z)
  u <- if (sigma_g2 > 0) {
    nc <- setdiff(seq_len(m), truth)
    as.numeric(X[, nc] %*% rnorm(length(nc), 0, sqrt(sigma_g2 / length(nc))))
  } else rep(0, n)
  y <- as.numeric(X[, truth] %*% beta) + u + rnorm(n)
  list(X = X, y = y, truth = truth, sigma_g2 = sigma_g2)
}

exact_metrics <- function(sel, truth) {
  K_hat <- length(sel); tp <- length(intersect(sel, truth))
  rec  <- tp / length(truth)
  prec <- if (K_hat > 0) tp / K_hat else 0
  f1   <- if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0
  data.frame(K_hat = K_hat, tp = tp, precision = prec, recall = rec, f1 = f1)
}

jobs <- expand.grid(rho = RHOS, sg = REGIMES, b = seq_len(B))

one_job <- function(i) {
  rho <- jobs$rho[i]; sg <- jobs$sg[i]; b <- jobs$b[i]
  ckpt <- file.path(OUT_DIR, sprintf("rho%.2f_sg%.1f_b%03d.rds", rho, sg, b))
  if (file.exists(ckpt) && file.size(ckpt) > 0) return("cached")
  seed <- 20260907L + 10000L * b + as.integer(round(100 * rho)) +
          as.integer(round(10 * sg))
  d <- gen_correlated(N, M, rho, TRUTH, BETA, sg, seed)

  sels <- list()
  fit <- CBF_LMM_stepwise_exact(d$y, d$X, tau2 = ANCHOR$tau2,
                                  K_max = ANCHOR$K_max,
                                  n_nodes = ANCHOR$N_delta)
  sels$MS_L_eBIC <- fit$indices
  cmp <- run_methods(d, c("SuSiE", "BayesR", "fastlmm", "BSLMM"),
                     K_max = ANCHOR$K_max, theta = ANCHOR$theta,
                     N_delta = ANCHOR$N_delta, tau2 = ANCHOR$tau2)
  for (mn in names(cmp)) sels[[mn]] <- cmp[[mn]]$indices

  rows <- do.call(rbind, lapply(names(sels), function(mn) {
    r <- exact_metrics(sels[[mn]], d$truth); r$method <- mn; r }))
  rows$rho <- rho; rows$sigma_g2 <- sg; rows$rep <- b
  saveRDS(list(rows = rows, indices = sels), ckpt)
  "done"
}

message(sprintf("45_correlated_causals: %d jobs on %d cores",
                nrow(jobs), N_CORES))
invisible(mclapply(seq_len(nrow(jobs)), function(i)
  tryCatch(one_job(i), error = function(e)
    message("FAIL job ", i, ": ", conditionMessage(e))),
  mc.cores = N_CORES, mc.preschedule = FALSE))

fs <- list.files(OUT_DIR, pattern = "^rho.*rds$", full.names = TRUE)
res <- do.call(rbind, lapply(fs, function(f) readRDS(f)$rows))
agg <- aggregate(cbind(K_hat, tp, precision, recall, f1) ~ method + rho + sigma_g2,
                 res, mean)
print(agg[order(agg$rho, agg$sigma_g2, -agg$f1), ], digits = 3, row.names = FALSE)
write.csv(agg, file.path(OUT_DIR, "correlated_causals_summary.csv"),
          row.names = FALSE)
message("done 45.")
