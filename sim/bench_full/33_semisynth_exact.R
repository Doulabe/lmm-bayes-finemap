# ==============================================================================
# 33_semisynth_exact.R
# Recompute the CBF-LMM rows of the semi-synthetic 1000G benchmark with the
# exact candidate-specific engine (K_jk + profile-ML eBIC), on the same
# datasets (same seed formula as 25_semisynth_1000g.R). The original
# checkpoint rows (all methods) are cloned; the MS_L_eBIC row is replaced and
# the old one kept as MS_L_eBIC_sharedGk.
# Output: results/bench_full_exact/25_semisynth_1000g/<same files> + summary
# Usage:  Rscript sim/bench_full/33_semisynth_exact.R [--cores 5]
# ==============================================================================

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressPackageStartupMessages(library(parallel))
source("sim/bench_full/00_config.R")
source("R/CBF_LMM_exact.R")

DATA_DIR <- "data/1000g"
LOCI <- list(chr1 = file.path(DATA_DIR, "locus_chr1_1000g.rds"),
             chr6 = file.path(DATA_DIR, "locus_chr6_1000g.rds"))
if (!is.na(i <- match("--locus_chr1", args <- commandArgs(TRUE))) && i < length(args))
  LOCI$chr1 <- args[i + 1L]
if (!is.na(i <- match("--locus_chr6", args)) && i < length(args))
  LOCI$chr6 <- args[i + 1L]
missing <- LOCI[!file.exists(unlist(LOCI))]
if (length(missing))
  stop("1000G locus panel(s) not found:\n  ",
       paste(unlist(missing), collapse = "\n  "),
       "\nBuild them first (see README) or pass --locus_chr1 / --locus_chr6.")
BETA   <- c(0.7, 0.5, 0.4); K_TRUE <- 3L; TAU <- 0.5
IN_DIR  <- "results/bench_full/25_semisynth_1000g"
OUT_DIR <- "results/bench_full_exact/25_semisynth_1000g"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(flag, default) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) as.integer(args[i + 1L]) else default
}
N_CORES <- arg_val("--cores", 5L)

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
ld_row <- function(sel, causal, R, method) {
  K_hat <- length(sel)
  if (K_hat == 0L)
    return(data.frame(K_hat = 0L, recall = 0, precision = 0, f1 = 0,
                      method = method))
  rec  <- mean(vapply(causal, function(cj) any(abs(R[sel, cj]) >= TAU), logical(1)))
  prec <- mean(vapply(sel,   function(s)  any(abs(R[s, causal]) >= TAU), logical(1)))
  f1   <- if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0
  data.frame(K_hat = K_hat, recall = rec, precision = prec, f1 = f1,
             method = method)
}

PAN <- lapply(names(LOCI), function(nm) {
  x <- readRDS(LOCI[[nm]]); X <- scale(x$G); X[!is.finite(X)] <- 0
  list(name = nm, X = X, R = x$R)
})
names(PAN) <- names(LOCI)

files <- list.files(IN_DIR, pattern = "^chr[16]_sg[0-9.]+_b[0-9]+\\.rds$")
one_file <- function(f) {
  out_fn <- file.path(OUT_DIR, f)
  if (file.exists(out_fn) && file.size(out_fn) > 0) return("cached")
  locus <- sub("_sg.*$", "", f)
  sg    <- as.numeric(sub("^chr[16]_sg([0-9.]+)_b.*$", "\\1", f))
  b     <- as.integer(sub(".*_b([0-9]+)\\.rds$", "\\1", f))
  seed  <- 20260620L + 1000L * b + as.integer(round(100 * sg)) +
           (if (locus == "chr6") 500000L else 0L)
  P <- PAN[[locus]]
  causal <- pick_causals(P$R, K_TRUE, seed)
  y <- make_y(P$X, causal, BETA, sg, seed)

  fit <- CBF_LMM_stepwise_exact(y, P$X, tau2 = ANCHOR$tau2,
                                  K_max = ANCHOR$K_max,
                                  n_nodes = ANCHOR$N_delta)

  rows <- readRDS(file.path(IN_DIR, f))
  old <- rows[rows$method == "MS_L_eBIC", ]
  old$method <- "MS_L_eBIC_sharedGk"
  new <- ld_row(fit$indices, causal, P$R, "MS_L_eBIC")
  new$locus <- locus; new$sigma_g2 <- sg; new$rep <- b
  rows <- rbind(rows[rows$method != "MS_L_eBIC", ],
                new[, names(rows)], old)
  saveRDS(rows, out_fn)
  "done"
}

message(sprintf("semisynth exact: %d replicates on %d cores",
                length(files), N_CORES))
invisible(mclapply(files, function(f)
  tryCatch(one_file(f), error = function(e)
    message("FAIL ", f, ": ", conditionMessage(e))),
  mc.cores = N_CORES, mc.preschedule = FALSE))

fs <- list.files(OUT_DIR, pattern = "^chr.*rds$", full.names = TRUE)
res <- do.call(rbind, lapply(fs, readRDS))
agg <- aggregate(cbind(K_hat, recall, precision, f1) ~ method + sigma_g2,
                 res, mean)
print(agg[order(agg$sigma_g2, -agg$f1), ], digits = 3, row.names = FALSE)
write.csv(agg, file.path(OUT_DIR, "semisynth_exact_summary.csv"),
          row.names = FALSE)
message("done.")
