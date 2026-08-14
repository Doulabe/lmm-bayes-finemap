# ==============================================================================
# 25_semisynth_1000g.R
# Semi-synthetic benchmark on REAL 1000G EUR genotype panels (real LD),
# addressing the limitation that the main simulation uses block-AR(1) LD only.
#
# Design (validated for n=503 real panels):
#   - 2 loci: chr1 (m=1493), chr6 (m=2143); 1000G phase3 EUR, n=503.
#   - K*=3 causals with low mutual LD (|r|<0.3), beta=(0.7,0.5,0.4)
#     (recalibrated from the block-AR anchor so that all signals are
#      detectable at n=503).
#   - sigma_g2 in {0, 0.5}; B replicates.
#   - LD-AWARE recovery (field standard for real-LD fine-mapping): a causal is
#     recovered if some selected SNP has r^2 >= 0.25 (|r| >= 0.5) with it; a
#     selection is a true positive if it tags some causal at r^2 >= 0.25.
#     Exact-index recovery is impossible here (causals have perfect-LD twins).
#   - 8-method set via run_methods() (MBF evaluation, same as the main figures).
#
# Reuses run_methods() from 00_config.R unchanged (X-agnostic).
#
# Usage (from project root):
#   Rscript sim/bench_full/25_semisynth_1000g.R --B 50 --cores 8
# ==============================================================================

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
suppressPackageStartupMessages({ library(parallel) })
source("sim/bench_full/00_config.R")
for (p in c("BGLR","hibayes","rrBLUP")) loadNamespace(p)   # pre-load before fork

# 1000G locus panels. Each .rds is a list with G (n x m, 0/1/2 genotypes),
# R (m x m LD matrix), pos, snpvar.  See README, "Semi-synthetic 1000G
# benchmark", for how to build them.  Override with --locus_chr1/--locus_chr6.
DATA_DIR <- "data/1000g"
LOCI <- list(chr1 = file.path(DATA_DIR, "locus_chr1_1000g.rds"),
             chr6 = file.path(DATA_DIR, "locus_chr6_1000g.rds"))

BETA    <- c(0.7, 0.5, 0.4)            # K*=3, recalibrated for n=503
K_TRUE  <- length(BETA)
REGIMES <- c(0, 0.5)
METHODS <- c("JS_L","MS_L","SuSiE","BayesR","fastlmm","BSLMM")   # -> 8 labels
TAU     <- 0.5                         # LD-aware threshold (|r| >= 0.5, r^2 >= 0.25)
OUTDIR  <- "results/bench_full/25_semisynth_1000g"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

pick_causals <- function(R, K, seed, r_thr = 0.3) {
  set.seed(seed); m <- nrow(R); sel <- integer(0)
  for (j in sample(m)) {
    if (length(sel) == 0L || all(abs(R[j, sel]) < r_thr)) sel <- c(sel, j)
    if (length(sel) == K) break
  }
  sort(sel)
}

make_y <- function(X, causal, beta, sigma_g2, seed) {
  set.seed(seed + 7L); n <- nrow(X)
  if (sigma_g2 > 0) {
    nc  <- setdiff(seq_len(ncol(X)), causal)
    Rbg <- tcrossprod(X[, nc]) / length(nc)
    ee  <- eigen(Rbg, symmetric = TRUE)
    half <- ee$vectors %*% diag(sqrt(pmax(sigma_g2 * ee$values, 0))) %*% t(ee$vectors)
    u <- as.numeric(half %*% rnorm(n))
  } else u <- rep(0, n)
  as.numeric(X[, causal] %*% beta) + u + rnorm(n)
}

# LD-aware recall / precision / F1 (a selection tags a causal if |r| >= TAU)
ld_metrics <- function(sel, causal, R) {
  K_hat <- length(sel)
  if (K_hat == 0L) return(data.frame(K_hat = 0L, recall = 0, precision = 0, f1 = 0))
  rec  <- mean(vapply(causal, function(cj) any(abs(R[sel, cj]) >= TAU), logical(1)))
  prec <- mean(vapply(sel,   function(s)  any(abs(R[s, causal]) >= TAU), logical(1)))
  f1   <- if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0
  data.frame(K_hat = K_hat, recall = rec, precision = prec, f1 = f1)
}

one_rep <- function(locus_name, X, R, sigma_g2, b) {
  ckpt <- file.path(OUTDIR, sprintf("%s_sg%.1f_b%02d.rds", locus_name, sigma_g2, b))
  if (file.exists(ckpt) && file.size(ckpt) > 0) return(invisible("cached"))
  seed <- 20260620L + 1000L * b + as.integer(round(100 * sigma_g2)) +
          (if (locus_name == "chr6") 500000L else 0L)
  causal <- pick_causals(R, K_TRUE, seed)
  y <- make_y(X, causal, BETA, sigma_g2, seed)
  d <- list(X = X, y = y, truth = causal, sigma_g2 = sigma_g2)
  res <- tryCatch(run_methods(d, METHODS, K_max = ANCHOR$K_max,
                              theta = ANCHOR$theta, N_delta = ANCHOR$N_delta,
                              tau2 = ANCHOR$tau2),
                  error = function(e) { message("ERR ", locus_name, " sg", sigma_g2,
                                                " b", b, ": ", conditionMessage(e)); NULL })
  if (is.null(res)) return(invisible("error"))
  rows <- do.call(rbind, lapply(names(res), function(mn) {
    m <- ld_metrics(res[[mn]]$indices, causal, R); m$method <- mn; m }))
  rows$locus <- locus_name; rows$sigma_g2 <- sigma_g2; rows$rep <- b
  saveRDS(rows, ckpt)
  invisible("done")
}

args <- parse_args()
B    <- if (!is.null(args$B)) as.integer(args$B) else 200L
if (!is.null(args$locus_chr1)) LOCI$chr1 <- args$locus_chr1
if (!is.null(args$locus_chr6)) LOCI$chr6 <- args$locus_chr6
missing <- LOCI[!file.exists(unlist(LOCI))]
if (length(missing))
  stop("1000G locus panel(s) not found:\n  ", paste(unlist(missing), collapse = "\n  "),
       "\nBuild them first (see README) or pass --locus_chr1 / --locus_chr6.")
NC   <- if (!is.null(args$cores)) as.integer(args$cores) else max(1L, detectCores() - 2L)

PAN <- lapply(names(LOCI), function(nm) {
  x <- readRDS(LOCI[[nm]]); X <- scale(x$G); X[!is.finite(X)] <- 0
  list(name = nm, X = X, R = x$R)
}); names(PAN) <- names(LOCI)

cells <- expand.grid(locus = names(LOCI), sigma_g2 = REGIMES, b = seq_len(B),
                     stringsAsFactors = FALSE)
cat(sprintf("Full run: %d loci x %d regimes x B=%d = %d cells, mc.cores=%d\n",
            length(LOCI), length(REGIMES), B, nrow(cells), NC))
t0 <- Sys.time()
out <- mclapply(seq_len(nrow(cells)), function(i) {
  cl <- cells[i, ]; p <- PAN[[cl$locus]]
  tryCatch(one_rep(cl$locus, p$X, p$R, cl$sigma_g2, cl$b),
           error = function(e) { message("CELL ERR ", i, ": ", conditionMessage(e)); "error" })
}, mc.cores = NC, mc.preschedule = FALSE)
cat(sprintf("Done %d cells in %.1f min: %s\n", nrow(cells),
            as.numeric(Sys.time() - t0, units = "mins"),
            paste(names(table(unlist(out))), table(unlist(out)), sep = "=", collapse = ", ")))

files <- list.files(OUTDIR, pattern = "_b[0-9]+\\.rds$", full.names = TRUE)
all <- do.call(rbind, lapply(files, readRDS))
agg <- aggregate(cbind(f1, recall, precision, K_hat) ~ method + sigma_g2, all, mean)
write.csv(all, file.path(OUTDIR, "semisynth_1000g_raw.csv"), row.names = FALSE)
write.csv(agg, file.path(OUTDIR, "semisynth_1000g_summary.csv"), row.names = FALSE)
cat("\n=== Aggregate (LD-aware, pooled over 2 loci x", B, "reps) ===\n")
print(format(agg[order(agg$sigma_g2, -agg$f1), ], digits = 3), row.names = FALSE)
