# ==============================================================================
# 43_pip_thresholds.R
# Operating-point sensitivity for the PIP-thresholded comparators (reviewer
# request): SuSiE and BayesR selections at PIP > {0.5, 0.9, 0.99} on the
# anchor configuration, both polygenic regimes, campaign seed formula.
# Exact-identity metrics.
# Output: results/bench_full_exact/43_pip_thresholds/
# Usage:  Rscript sim/bench_full/43_pip_thresholds.R [--cores 5] [--B 100]
# ==============================================================================

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressPackageStartupMessages(library(parallel))
source("sim/bench_full/00_config.R")
for (p in c("hibayes", "susieR")) loadNamespace(p)

OUT_DIR <- "results/bench_full_exact/43_pip_thresholds"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(flag, default) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) as.integer(args[i + 1L]) else default
}
N_CORES <- arg_val("--cores", 5L)
B       <- arg_val("--B", 100L)
REGIMES <- c(0, 0.5)
THRESH  <- c(0.5, 0.9, 0.99)

metrics_at <- function(pip, thr, truth) {
  sel <- which(pip > thr)
  tp <- length(intersect(sel, truth)); K <- length(sel)
  pr <- if (K > 0) tp / K else 0
  rc <- tp / length(truth)
  f1 <- if (pr + rc > 0) 2 * pr * rc / (pr + rc) else 0
  data.frame(threshold = thr, K_hat = K, precision = pr, recall = rc, f1 = f1)
}

jobs <- expand.grid(sg = REGIMES, b = seq_len(B))

one_job <- function(i) {
  sg <- jobs$sg[i]; b <- jobs$b[i]
  ckpt <- file.path(OUT_DIR, sprintf("sg%.1f_b%03d.rds", sg, b))
  if (file.exists(ckpt) && file.size(ckpt) > 0) return("cached")
  pc <- ANCHOR
  seed <- 20260101L + 100L * b + as.integer(round(10 * sg))
  d <- gen_dataset(n = pc$n, m = pc$m, rho = pc$rho,
                   K_true = pc$K_true, beta_true = pc$beta,
                   sigma_g2 = sg, seed = seed)

  su <- tryCatch(susieR::susie(d$X, d$y, L = pc$K_max, verbose = FALSE),
                 error = function(e) NULL)
  pip_su <- if (!is.null(su)) as.numeric(su$pip) else rep(0, ncol(d$X))

  Xm <- d$X; rownames(Xm) <- paste0("S", seq_len(nrow(d$X)))
  df_y <- data.frame(id = rownames(Xm), y = d$y, stringsAsFactors = FALSE)
  br <- tryCatch(suppressMessages(suppressWarnings(
          hibayes::ibrm(formula = y ~ 1, data = df_y, M = Xm, M.id = df_y$id,
                        method = "BayesR",
                        Pi = c(0.95, 0.02, 0.02, 0.01),
                        fold = c(0, 1e-4, 1e-3, 1e-2),
                        niter = 2000L, nburn = 400L,
                        threads = 1L, verbose = FALSE))),
        error = function(e) NULL)
  pip_br <- if (!is.null(br) && !is.null(br$pip)) as.numeric(br$pip)
            else rep(0, ncol(d$X))

  rows <- rbind(
    do.call(rbind, lapply(THRESH, function(th) {
      r <- metrics_at(pip_su, th, d$truth); r$method <- "SuSiE"; r })),
    do.call(rbind, lapply(THRESH, function(th) {
      r <- metrics_at(pip_br, th, d$truth); r$method <- "BayesR"; r })))
  rows$sigma_g2 <- sg; rows$rep <- b
  saveRDS(rows, ckpt)
  "done"
}

message(sprintf("43_pip_thresholds: %d jobs on %d cores", nrow(jobs), N_CORES))
invisible(mclapply(seq_len(nrow(jobs)), function(i)
  tryCatch(one_job(i), error = function(e)
    message("FAIL job ", i, ": ", conditionMessage(e))),
  mc.cores = N_CORES, mc.preschedule = FALSE))

fs <- list.files(OUT_DIR, pattern = "^sg.*rds$", full.names = TRUE)
res <- do.call(rbind, lapply(fs, readRDS))
agg <- aggregate(cbind(K_hat, precision, recall, f1) ~ method + threshold + sigma_g2,
                 res, mean)
print(agg[order(agg$method, agg$sigma_g2, agg$threshold), ],
      digits = 3, row.names = FALSE)
write.csv(agg, file.path(OUT_DIR, "pip_threshold_summary.csv"),
          row.names = FALSE)
message("done 43.")
