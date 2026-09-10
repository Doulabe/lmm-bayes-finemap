# ==============================================================================
# 50_bayesr_long_controlled.R
# Rerun BayesR ONLY at the long MCMC budget (10,000 iterations / 2,000 burn-in)
# on the controlled benchmark (four axes, 2,200 replicates), following the
# chain-length control of Appendix D. Datasets are reconstructed with the
# ORIGINAL per-axis seed formulas (30_exact_rerun.R); all other method rows
# are untouched. Payloads in results/bench_full_exact/<axis>/ are updated IN
# PLACE atomically (tmp + rename): the short-chain row is preserved as
# BayesR_2k, the BayesR row is replaced by the long-chain fit, and the
# metrics table is recomputed. Idempotent: payloads already carrying
# BayesR_2k are skipped. A second pass patches the Top-K* payloads
# (results/bench_full_exact/topk) from the stored long-chain rankings, with a
# truth-identity check and a direct recomputation fallback.
# Usage: Rscript sim/bench_full/50_bayesr_long_controlled.R [--cores 5]
# ==============================================================================

Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressPackageStartupMessages(library(parallel))
source("sim/bench_full/00_config.R")
loadNamespace("hibayes")

EX <- "results/bench_full_exact"
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
  NULL
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

# ── per-axis seed formulas (identical to 30_exact_rerun.R) ───────────────────
parse_cell <- function(axis, file) {
  rep_id <- as.integer(sub(".*_b([0-9]+)\\.rds$", "\\1", file))
  sg     <- as.numeric(sub(".*_sg([0-9.]+)_b[0-9]+\\.rds$", "\\1", file))
  out <- list(rep = rep_id, sg = sg,
              n = ANCHOR$n, m = ANCHOR$m, rho = ANCHOR$rho,
              beta = ANCHOR$beta_true)
  if (axis == "01_scaling_n") {
    out$n <- as.integer(sub("^n([0-9]+)_.*", "\\1", file))
    out$seed <- 20260425L + 1000L * (rep_id - 1L) + out$n +
                as.integer(round(1000 * sg))
  } else if (axis == "02_scaling_rho") {
    out$rho <- as.numeric(sub("^rho([0-9.]+)_sg.*", "\\1", file))
    out$seed <- 20260425L + 1000L * (rep_id - 1L) +
                as.integer(round(1000 * out$rho)) +
                as.integer(round(100 * sg))
  } else if (axis == "03_arch") {
    sig <- sub("^sig([a-z]+)_sg.*", "\\1", file)
    out$beta <- switch(sig,
                       weak   = c(0.4, 0.2, 0.2, 0.1, 0.1),
                       medium = c(0.8, 0.4, 0.4, 0.2, 0.2),
                       strong = c(1.6, 0.8, 0.8, 0.4, 0.4))
    out$seed <- 20260425L + 1000L * (rep_id - 1L) +
                switch(sig, weak = 1L, medium = 2L, strong = 3L) +
                as.integer(round(100 * sg))
  } else {
    out$m <- as.integer(sub("^m([0-9]+)_.*", "\\1", file))
    out$seed <- 20260425L + 1000L * (rep_id - 1L) + out$m +
                as.integer(round(100 * sg))
  }
  out
}

AXES <- c("01_scaling_n", "02_scaling_rho", "03_arch", "04_scaling_m")
work <- do.call(rbind, lapply(AXES, function(ax) {
  fs <- list.files(file.path(EX, ax), pattern = "_b[0-9]+\\.rds$")
  if (!length(fs)) return(NULL)
  data.frame(axis = ax, file = fs, stringsAsFactors = FALSE)
}))

one_file <- function(i) {
  axis <- work$axis[i]; file <- work$file[i]
  fn <- file.path(EX, axis, file)
  p <- read_retry(fn)
  if (is.null(p)) stop("unreadable: ", fn)
  if ("BayesR_2k" %in% names(p$results)) return("cached")

  pc <- parse_cell(axis, file)
  d  <- gen_dataset(n = pc$n, m = pc$m, rho = pc$rho,
                    K_true = length(pc$beta), beta_true = pc$beta,
                    sigma_g2 = pc$sg, block_size = ANCHOR$block_size,
                    seed = pc$seed)
  stopifnot(identical(sort(p$truth), sort(d$truth)))

  t0 <- Sys.time()
  pip <- bayesr_long(d$y, d$X, seed = pc$seed + 23L)
  el  <- as.numeric(Sys.time() - t0, units = "secs")
  sel <- which(pip > 0.99)
  Ks  <- length(p$truth)

  p$results$BayesR_2k <- p$results$BayesR
  p$results$BayesR <- list(method = "BayesR", indices = sel,
                           K_hat = length(sel), scores = pip[sel],
                           elapsed = el,
                           top_idx = order(-pip)[seq_len(Ks)],
                           niter = NITER, nburn = NBURN)
  metrics <- compute_metrics(p$results, p$truth)
  extra   <- setdiff(names(p$metrics), names(metrics))
  for (cn in extra) metrics[[cn]] <- p$metrics[[cn]][1L]
  p$metrics <- metrics

  tmp <- paste0(fn, ".tmp")
  saveRDS(p, tmp); file.rename(tmp, fn)
  "done"
}

if (LIMIT > 0L) work <- head(work, LIMIT)
message(sprintf("50_bayesr_long_controlled: %d payloads on %d cores",
                nrow(work), N_CORES))
invisible(mclapply(seq_len(nrow(work)), function(i)
  tryCatch(one_file(i), error = function(e)
    message("FAIL ", work$axis[i], "/", work$file[i], ": ",
            conditionMessage(e))),
  mc.cores = N_CORES, mc.preschedule = FALSE))

## ── pass 2: patch the Top-K* payloads from the stored long-chain rankings ──
TOPK_MAP <- list(
  anchor    = c("01_scaling_n",  "n1000"),
  n3000     = c("01_scaling_n",  "n3000"),
  n500      = c("01_scaling_n",  "n0500"),
  rho098    = c("02_scaling_rho","rho0.98"),
  rho080    = c("02_scaling_rho","rho0.80"),
  sigweak   = c("03_arch",       "sigweak"),
  sigstrong = c("03_arch",       "sigstrong"),
  m10000    = c("04_scaling_m",  "m10000"))

tk_fs <- list.files(file.path(EX, "topk"), pattern = "_b[0-9]+\\.rds$",
                    full.names = TRUE)
patched <- 0L; fallback <- 0L; unread <- 0L
for (tf in tk_fs) {
  tp <- read_retry(tf)
  if (is.null(tp)) { unread <- unread + 1L
    message("SKIP unreadable topk: ", tf); next }
  if (!is.null(tp$metrics$BayesR$top_K_recall_2k)) next
  tag <- sub("_sg[0-9.]+$", "", tp$cell_tag)
  mm  <- TOPK_MAP[[tag]]
  src <- file.path(EX, mm[1],
                   sprintf("%s_sg%.1f_b%02d.rds", mm[2], tp$sigma_g2, tp$rep))
  ok <- FALSE
  if (file.exists(src)) {
    ap <- read_retry(src)
    if (!is.null(ap) && identical(sort(ap$truth), sort(tp$truth)) &&
        !is.null(ap$results$BayesR$top_idx)) {
      top_idx <- ap$results$BayesR$top_idx
      ok <- TRUE
    }
  }
  if (!ok) { fallback <- fallback + 1L; next }
  Ks <- length(tp$truth)
  tp$metrics$BayesR$top_K_recall_2k <- tp$metrics$BayesR$top_K_recall
  tp$metrics$BayesR$top_K_recall <-
    length(intersect(top_idx, tp$truth)) / Ks
  tmp <- paste0(tf, ".tmp"); saveRDS(tp, tmp); file.rename(tmp, tf)
  patched <- patched + 1L
}
message(sprintf(
  "topk patch: %d patched, %d skipped (axis payload not ready), %d unreadable",
  patched, fallback, unread))

## quick aggregate check
M <- do.call(rbind, lapply(AXES, function(ax) {
  fs <- list.files(file.path(EX, ax), full.names = TRUE,
                   pattern = "_b[0-9]+\\.rds$")
  do.call(rbind, lapply(fs, function(f) {
    p <- read_retry(f, attempts = 3L)
    if (is.null(p)) return(NULL)
    m <- p$metrics
    m[m$method %in% c("BayesR", "BayesR_2k"), c("method", "K_hat", "f1")]
  }))
}))
print(aggregate(cbind(K_hat, f1) ~ method, M, mean), digits = 3)
message("done 50.")
