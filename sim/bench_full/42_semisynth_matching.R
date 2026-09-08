# ==============================================================================
# 42_semisynth_matching.R
# Rerun the full semi-synthetic 1000G benchmark (all methods) on the canonical
# exact-era datasets (seed formulas of 33_semisynth_exact.R) and score with
# BOTH metrics:
#   - one-to-one maximum-matching LD-aware metrics (new primary)
#   - legacy bidirectional tagging metrics (continuity check)
# Selected indices are stored per replicate this time.
# Integrity check: the legacy-metric CBF rows must reproduce the current
# Table 4 values exactly (same datasets, same engine, same metric).
# Output: results/bench_full_exact/42_semisynth_matching/
# Usage:  Rscript sim/bench_full/42_semisynth_matching.R [--cores 5] [--B 200]
# ==============================================================================

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressPackageStartupMessages(library(parallel))
source("sim/bench_full/00_config.R")
source("sim/bench_full/matching_lib.R")
source("R/CBF_LMM_exact.R")
for (p in c("BGLR", "hibayes", "rrBLUP", "susieR")) loadNamespace(p)

LOCI <- list(
  chr1 = "/Users/kossi/Desktop/Dossier_these/Redaction/Projet4_Binary_Individual/results_realdata/polyfun_1000g/locus_1000g_polyfun.rds",
  chr6 = "/Users/kossi/Desktop/Dossier_these/Redaction/Projet1_donnees_brutes/locus2_chr6_1000g_polyfun.rds")
BETA <- c(0.7, 0.5, 0.4); K_TRUE <- 3L; TAU <- 0.5
REGIMES <- c(0, 0.5)
OUT_DIR <- "results/bench_full_exact/42_semisynth_matching"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(flag, default) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) as.integer(args[i + 1L]) else default
}
N_CORES <- arg_val("--cores", 5L)
B       <- arg_val("--B", 200L)

# --- canonical dataset construction (33_semisynth_exact.R) ---
pick_causals <- function(R, K, seed, r_thr = 0.3) {
  set.seed(seed); m <- nrow(R); sel <- integer(0)
  while (length(sel) < K) {
    cand <- sample.int(m, 1L)
    if (!length(sel) || all(abs(R[cand, sel]) < r_thr)) sel <- c(sel, cand)
  }
  sel
}
make_y <- function(X, causal, beta, sigma_g2, seed) {
  set.seed(seed + 7L); n <- nrow(X)
  u <- if (sigma_g2 > 0) {
    nc <- setdiff(seq_len(ncol(X)), causal)
    Z <- X[, nc, drop = FALSE]
    as.numeric(Z %*% rnorm(ncol(Z), 0, sqrt(sigma_g2 / ncol(Z))))
  } else rep(0, n)
  as.numeric(X[, causal] %*% beta) + u + rnorm(n)
}

PAN <- lapply(names(LOCI), function(nm) {
  x <- readRDS(LOCI[[nm]]); X <- scale(x$G); X[!is.finite(X)] <- 0
  list(name = nm, X = X, R = x$R)
})
names(PAN) <- names(LOCI)

jobs <- expand.grid(locus = names(LOCI), sg = REGIMES, b = seq_len(B),
                    stringsAsFactors = FALSE)

one_job <- function(i) {
  locus <- jobs$locus[i]; sg <- jobs$sg[i]; b <- jobs$b[i]
  ckpt <- file.path(OUT_DIR, sprintf("%s_sg%.1f_b%02d.rds", locus, sg, b))
  if (file.exists(ckpt) && file.size(ckpt) > 0) return("cached")
  seed <- 20260620L + 1000L * b + as.integer(round(100 * sg)) +
          (if (locus == "chr6") 500000L else 0L)
  P <- PAN[[locus]]
  causal <- pick_causals(P$R, K_TRUE, seed)
  y <- make_y(P$X, causal, BETA, sg, seed)
  d <- list(X = P$X, y = y, truth = causal, sigma_g2 = sg)

  sels <- list()
  fit <- CBF_LMM_stepwise_exact(y, P$X, tau2 = ANCHOR$tau2,
                                  K_max = ANCHOR$K_max,
                                  n_nodes = ANCHOR$N_delta)
  sels$MS_L_eBIC <- fit$indices
  cmp <- run_methods(d, c("SuSiE", "BayesR", "fastlmm", "BSLMM"),
                     K_max = ANCHOR$K_max, theta = ANCHOR$theta,
                     N_delta = ANCHOR$N_delta, tau2 = ANCHOR$tau2)
  for (mn in names(cmp)) sels[[mn]] <- cmp[[mn]]$indices

  rows <- do.call(rbind, lapply(names(sels), function(mn) {
    mm <- matching_metrics(sels[[mn]], causal, P$R, TAU)
    tg <- tagging_metrics(sels[[mn]],  causal, P$R, TAU)
    data.frame(method = mn, K_hat = mm$K_hat, TP = mm$TP,
               precision = mm$precision, recall = mm$recall, f1 = mm$f1,
               precision_tag = tg$precision, recall_tag = tg$recall,
               f1_tag = tg$f1)
  }))
  rows$locus <- locus; rows$sigma_g2 <- sg; rows$rep <- b
  saveRDS(list(rows = rows, indices = sels, causal = causal), ckpt)
  "done"
}

message(sprintf("42_semisynth_matching: %d jobs on %d cores",
                nrow(jobs), N_CORES))
invisible(mclapply(seq_len(nrow(jobs)), function(i)
  tryCatch(one_job(i), error = function(e)
    message("FAIL job ", i, ": ", conditionMessage(e))),
  mc.cores = N_CORES, mc.preschedule = FALSE))

fs <- list.files(OUT_DIR, pattern = "^chr.*rds$", full.names = TRUE)
res <- do.call(rbind, lapply(fs, function(f) readRDS(f)$rows))
agg <- aggregate(cbind(K_hat, precision, recall, f1,
                       precision_tag, recall_tag, f1_tag) ~ method + sigma_g2,
                 res, mean)
print(agg[order(agg$sigma_g2, -agg$f1), ], digits = 3, row.names = FALSE)
write.csv(agg, file.path(OUT_DIR, "matching_summary.csv"), row.names = FALSE)
aggl <- aggregate(cbind(K_hat, precision, recall, f1) ~ method + sigma_g2 + locus,
                  res, mean)
write.csv(aggl, file.path(OUT_DIR, "matching_bylocus.csv"), row.names = FALSE)
message("done 42.")
