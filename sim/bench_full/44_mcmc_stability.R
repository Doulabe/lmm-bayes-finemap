# ==============================================================================
# 44_mcmc_stability.R
# MCMC convergence control for the two sampling comparators (reviewer request):
# BayesB (BGLR) and BayesR (hibayes) at the benchmark settings
# (2,000 iterations / 400 burn-in) versus a long run (10,000 / 2,000),
# paired on the same anchor datasets. Reports per-replicate agreement of the
# PIP vector (correlation), the Top-K* set, and the thresholded selection.
# Output: results/bench_full_exact/44_mcmc_stability/
# Usage:  Rscript sim/bench_full/44_mcmc_stability.R [--cores 5] [--B 50]
# ==============================================================================

# --- reviewer-reserve guard: create sim/bench_full/SKIP_44 to defer this run ---
if (file.exists("sim/bench_full/SKIP_44")) {
  message("44_mcmc_stability deferred (reviewer reserve); remove sim/bench_full/SKIP_44 to run.")
  quit(save = "no")
}

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressPackageStartupMessages(library(parallel))
source("sim/bench_full/00_config.R")
for (p in c("BGLR", "hibayes")) loadNamespace(p)

OUT_DIR <- "results/bench_full_exact/44_mcmc_stability"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(flag, default) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) as.integer(args[i + 1L]) else default
}
N_CORES <- arg_val("--cores", 5L)
B       <- arg_val("--B", 50L)
REGIMES <- c(0, 0.5)

bayesb_pip <- function(y, X, nIter, burnIn, seed) {
  set.seed(seed)
  fit <- suppressMessages(suppressWarnings(
    BGLR::BGLR(y = y, ETA = list(list(X = X, model = "BayesB")),
               nIter = nIter, burnIn = burnIn, verbose = FALSE)))
  as.numeric(fit$ETA[[1]]$d)
}
bayesr_pip <- function(y, X, niter, nburn, seed) {
  set.seed(seed)
  Xm <- X; rownames(Xm) <- paste0("S", seq_len(nrow(X)))
  df_y <- data.frame(id = rownames(Xm), y = y, stringsAsFactors = FALSE)
  fit <- suppressMessages(suppressWarnings(
    hibayes::ibrm(formula = y ~ 1, data = df_y, M = Xm, M.id = df_y$id,
                  method = "BayesR",
                  Pi = c(0.95, 0.02, 0.02, 0.01),
                  fold = c(0, 1e-4, 1e-3, 1e-2),
                  niter = niter, nburn = nburn,
                  threads = 1L, verbose = FALSE)))
  as.numeric(fit$pip)
}

cmp_row <- function(pip_s, pip_l, truth, thr, method, sg, b) {
  Ks <- length(truth)
  top_s <- order(-pip_s)[seq_len(Ks)]; top_l <- order(-pip_l)[seq_len(Ks)]
  sel_s <- which(pip_s > thr);         sel_l <- which(pip_l > thr)
  f1 <- function(sel) {
    tp <- length(intersect(sel, truth)); K <- length(sel)
    pr <- if (K > 0) tp / K else 0; rc <- tp / Ks
    if (pr + rc > 0) 2 * pr * rc / (pr + rc) else 0
  }
  data.frame(method = method, sigma_g2 = sg, rep = b,
             pip_cor = suppressWarnings(cor(pip_s, pip_l)),
             pip_spear = suppressWarnings(cor(pip_s, pip_l, method = "spearman")),
             topK_overlap = length(intersect(top_s, top_l)) / Ks,
             K_short = length(sel_s), K_long = length(sel_l),
             sel_jaccard = if (length(union(sel_s, sel_l)) > 0)
               length(intersect(sel_s, sel_l)) / length(union(sel_s, sel_l)) else 1,
             f1_short = f1(sel_s), f1_long = f1(sel_l))
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
  pb_s <- bayesb_pip(d$y, d$X, 2000L, 400L,  seed + 11L)
  pb_l <- bayesb_pip(d$y, d$X, 10000L, 2000L, seed + 13L)
  pr_s <- bayesr_pip(d$y, d$X, 2000L, 400L,  seed + 17L)
  pr_l <- bayesr_pip(d$y, d$X, 10000L, 2000L, seed + 19L)
  rows <- rbind(cmp_row(pb_s, pb_l, d$truth, 0.5,  "BayesB", sg, b),
                cmp_row(pr_s, pr_l, d$truth, 0.99, "BayesR", sg, b))
  saveRDS(rows, ckpt)
  "done"
}

message(sprintf("44_mcmc_stability: %d jobs on %d cores", nrow(jobs), N_CORES))
invisible(mclapply(seq_len(nrow(jobs)), function(i)
  tryCatch(one_job(i), error = function(e)
    message("FAIL job ", i, ": ", conditionMessage(e))),
  mc.cores = N_CORES, mc.preschedule = FALSE))

fs <- list.files(OUT_DIR, pattern = "^sg.*rds$", full.names = TRUE)
res <- do.call(rbind, lapply(fs, readRDS))
agg <- aggregate(cbind(pip_cor, pip_spear, topK_overlap, K_short, K_long,
                       sel_jaccard, f1_short, f1_long) ~ method + sigma_g2,
                 res, mean)
print(agg, digits = 3, row.names = FALSE)
write.csv(agg, file.path(OUT_DIR, "mcmc_stability_summary.csv"),
          row.names = FALSE)
message("done 44.")
