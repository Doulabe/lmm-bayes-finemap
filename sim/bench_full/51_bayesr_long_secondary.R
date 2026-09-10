# ==============================================================================
# 51_bayesr_long_secondary.R
# BayesR at the long budget (10,000/2,000) for the secondary benchmarks:
#   C. semisynthetic matching benchmark  (42, 800 replicates, in-place)
#   D. correlated-causal stress test     (45, 400 replicates, in-place)
#   E. PIP-threshold sensitivity         (43, 200 replicates, in-place)
# Same in-place atomic-update pattern as 50; short-chain results preserved
# under *_2k names. Idempotent.
# Usage: Rscript sim/bench_full/51_bayesr_long_secondary.R [--cores 5]
# ==============================================================================

Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressPackageStartupMessages(library(parallel))
source("sim/bench_full/00_config.R")
source("sim/bench_full/matching_lib.R")
loadNamespace("hibayes")

NITER <- 10000L; NBURN <- 2000L
args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(flag, default) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) as.integer(args[i + 1L]) else default
}
N_CORES <- arg_val("--cores", 5L)
LIMIT   <- arg_val("--limit", 0L)

# retry-hardened readRDS: rides out transient iCloud/FS read failures
read_retry <- function(fn, attempts = 10L) {
  for (att in seq_len(attempts)) {
    p <- tryCatch(readRDS(fn), error = function(e) NULL)
    if (!is.null(p)) return(p)
    Sys.sleep(min(2 * att, 20))
  }
  stop("unreadable after ", attempts, " attempts: ", fn)
}

bayesr_long <- function(y, X, seed) {
  set.seed(seed)
  Xm <- X; rownames(Xm) <- paste0("S", seq_len(nrow(X)))
  df_y <- data.frame(id = rownames(Xm), y = y, stringsAsFactors = FALSE)
  fit <- suppressMessages(suppressWarnings(
    hibayes::ibrm(formula = y ~ 1, data = df_y, M = Xm, M.id = df_y$id,
                  method = "BayesR",
                  Pi = c(0.95, 0.02, 0.02, 0.01),
                  fold = c(0, 1e-4, 1e-3, 1e-2),
                  niter = NITER, nburn = NBURN,
                  threads = 1L, verbose = FALSE)))
  as.numeric(fit$pip)
}
atomic_save <- function(obj, fn) {
  tmp <- paste0(fn, ".tmp"); saveRDS(obj, tmp); file.rename(tmp, fn)
}

## ── C. semisynthetic (42) ────────────────────────────────────────────────────
LOCI <- list(
  chr1 = "/Users/kossi/Desktop/Dossier_these/Redaction/Projet4_Binary_Individual/results_realdata/polyfun_1000g/locus_1000g_polyfun.rds",
  chr6 = "/Users/kossi/Desktop/Dossier_these/Redaction/Projet1_donnees_brutes/locus2_chr6_1000g_polyfun.rds")
BETA <- c(0.7, 0.5, 0.4); TAU <- 0.5
DIR42 <- "results/bench_full_exact/42_semisynth_matching"
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
  x <- read_retry(LOCI[[nm]]); X <- scale(x$G); X[!is.finite(X)] <- 0
  list(name = nm, X = X, R = x$R)
})
names(PAN) <- names(LOCI)

fs42 <- list.files(DIR42, pattern = "^chr.*_b[0-9]+\\.rds$")
one_42 <- function(f) {
  fn <- file.path(DIR42, f)
  p <- read_retry(fn)
  if ("BayesR_2k" %in% names(p$indices)) return("cached")
  locus <- sub("_sg.*$", "", f)
  sg    <- as.numeric(sub("^chr[16]_sg([0-9.]+)_b.*$", "\\1", f))
  b     <- as.integer(sub(".*_b([0-9]+)\\.rds$", "\\1", f))
  seed  <- 20260620L + 1000L * b + as.integer(round(100 * sg)) +
           (if (locus == "chr6") 500000L else 0L)
  P <- PAN[[locus]]
  causal <- pick_causals(P$R, 3L, seed)
  stopifnot(identical(causal, p$causal))
  y <- make_y(P$X, causal, BETA, sg, seed)
  pip <- bayesr_long(y, P$X, seed + 23L)
  sel <- which(pip > 0.99)
  p$indices$BayesR_2k <- p$indices$BayesR
  p$indices$BayesR <- sel
  mm <- matching_metrics(sel, causal, P$R, TAU)
  tg <- tagging_metrics(sel,  causal, P$R, TAU)
  new_row <- data.frame(method = "BayesR", K_hat = mm$K_hat, TP = mm$TP,
                        precision = mm$precision, recall = mm$recall,
                        f1 = mm$f1, precision_tag = tg$precision,
                        recall_tag = tg$recall, f1_tag = tg$f1,
                        locus = locus, sigma_g2 = sg, rep = b)
  old <- p$rows[p$rows$method == "BayesR", ]
  old$method <- "BayesR_2k"
  p$rows <- rbind(p$rows[p$rows$method != "BayesR", ],
                  new_row[, names(p$rows)], old)
  atomic_save(p, fn)
  "done"
}

## ── D. correlated causals (45) ───────────────────────────────────────────────
DIR45 <- "results/bench_full_exact/45_correlated_causals"
TRUTH45 <- c(1L, 3L, 5L); BETA45 <- c(0.8, 0.4, 0.4)
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
  list(X = X, y = y, truth = truth)
}
exact_metrics45 <- function(sel, truth) {
  K_hat <- length(sel); tp <- length(intersect(sel, truth))
  rec  <- tp / length(truth)
  prec <- if (K_hat > 0) tp / K_hat else 0
  f1   <- if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0
  data.frame(K_hat = K_hat, tp = tp, precision = prec, recall = rec, f1 = f1)
}
fs45 <- list.files(DIR45, pattern = "^rho.*_b[0-9]+\\.rds$")
one_45 <- function(f) {
  fn <- file.path(DIR45, f)
  p <- read_retry(fn)
  if ("BayesR_2k" %in% names(p$indices)) return("cached")
  rho <- as.numeric(sub("^rho([0-9.]+)_sg.*$", "\\1", f))
  sg  <- as.numeric(sub("^rho[0-9.]+_sg([0-9.]+)_b.*$", "\\1", f))
  b   <- as.integer(sub(".*_b([0-9]+)\\.rds$", "\\1", f))
  seed <- 20260907L + 10000L * b + as.integer(round(100 * rho)) +
          as.integer(round(10 * sg))
  d <- gen_correlated(1000L, 5000L, rho, TRUTH45, BETA45, sg, seed)
  pip <- bayesr_long(d$y, d$X, seed + 23L)
  sel <- which(pip > 0.99)
  p$indices$BayesR_2k <- p$indices$BayesR
  p$indices$BayesR <- sel
  r <- exact_metrics45(sel, TRUTH45); r$method <- "BayesR"
  r$rho <- rho; r$sigma_g2 <- sg; r$rep <- b
  old <- p$rows[p$rows$method == "BayesR", ]
  old$method <- "BayesR_2k"
  p$rows <- rbind(p$rows[p$rows$method != "BayesR", ],
                  r[, names(p$rows)], old)
  atomic_save(p, fn)
  "done"
}

## ── E. PIP thresholds (43) ───────────────────────────────────────────────────
DIR43 <- "results/bench_full_exact/43_pip_thresholds"
THRESH <- c(0.5, 0.9, 0.99)
metrics_at <- function(pip, thr, truth) {
  sel <- which(pip > thr)
  tp <- length(intersect(sel, truth)); K <- length(sel)
  pr <- if (K > 0) tp / K else 0
  rc <- tp / length(truth)
  f1 <- if (pr + rc > 0) 2 * pr * rc / (pr + rc) else 0
  data.frame(threshold = thr, K_hat = K, precision = pr, recall = rc, f1 = f1)
}
fs43 <- list.files(DIR43, pattern = "^sg.*_b[0-9]+\\.rds$")
one_43 <- function(f) {
  fn <- file.path(DIR43, f)
  res <- read_retry(fn)
  if ("BayesR_2k" %in% res$method) return("cached")
  sg <- as.numeric(sub("^sg([0-9.]+)_b.*$", "\\1", f))
  b  <- as.integer(sub(".*_b([0-9]+)\\.rds$", "\\1", f))
  pc <- ANCHOR
  seed <- 20260101L + 100L * b + as.integer(round(10 * sg))
  d <- gen_dataset(n = pc$n, m = pc$m, rho = pc$rho,
                   K_true = pc$K_true, beta_true = pc$beta,
                   sigma_g2 = sg, seed = seed)
  pip <- bayesr_long(d$y, d$X, seed + 23L)
  new_rows <- do.call(rbind, lapply(THRESH, function(th) {
    r <- metrics_at(pip, th, d$truth); r$method <- "BayesR"; r }))
  new_rows$sigma_g2 <- sg; new_rows$rep <- b
  old <- res[res$method == "BayesR", ]
  old$method <- "BayesR_2k"
  res <- rbind(res[res$method != "BayesR", ],
               new_rows[, names(res)], old)
  atomic_save(res, fn)
  "done"
}

## ── run all three, checkpointed ─────────────────────────────────────────────
run_block <- function(fs, fun, label) {
  if (LIMIT > 0L) fs <- head(fs, LIMIT)
  message(sprintf("51/%s: %d payloads on %d cores", label, length(fs), N_CORES))
  invisible(mclapply(fs, function(f)
    tryCatch(fun(f), error = function(e)
      message("FAIL ", label, "/", f, ": ", conditionMessage(e))),
    mc.cores = N_CORES, mc.preschedule = FALSE))
}
run_block(fs42, one_42, "semisynth")
run_block(fs45, one_45, "correlated")
run_block(fs43, one_43, "thresholds")
message("done 51.")
