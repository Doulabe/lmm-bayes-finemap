# ==============================================================================
# 52_bayesr_long_realdata.R
# BayesR at the long budget (10,000/2,000) for the two real-data analyses:
#   F1. GEUVADIS cis-eQTL loci (8 genes)  — in-place payload update
#   F2. heterogeneous-stock mouse BMI      — standalone result file
# Short-chain GEUVADIS rows preserved as BayesR_2k. Idempotent.
# Usage: Rscript sim/bench_full/52_bayesr_long_realdata.R
# ==============================================================================

Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
source("sim/bench_full/00_config.R")
loadNamespace("hibayes")

NITER <- 10000L; NBURN <- 2000L
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

## ── F1. GEUVADIS ─────────────────────────────────────────────────────────────
D       <- "data/geuvadis"
LOC_DIR <- file.path(D, "loci")
OUT_DIR <- "results/bench_full_exact/26_geuvadis"
SAMP    <- file.path(D, "geuvadis_eur_overlap.txt")
EXPR_GZ <- file.path(D, "GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz")

GENES <- data.frame(
  gene   = c("ENSG00000164308.12", "ENSG00000197728.5", "ENSG00000166750.4",
             "ENSG00000203875.4",  "ENSG00000198468.2", "ENSG00000124587.9",
             "ENSG00000230658.1",  "ENSG00000174652.12"),
  lead_pos = c(96252589L, 56401085L, 33571546L, 86387888L,
               213049214L, 42944850L, 23143113L, 9544276L),
  stringsAsFactors = FALSE)

hdr <- strsplit(readLines(gzfile(EXPR_GZ), n = 1L), "\t")[[1]]
expr_rows <- list()
con <- gzfile(EXPR_GZ); open(con); readLines(con, n = 1L)
while (length(l <- readLines(con, n = 4000L))) {
  for (g in GENES$gene) {
    hit <- grep(paste0("^", g, "\t"), l)
    if (length(hit)) expr_rows[[g]] <- strsplit(l[hit[1]], "\t")[[1]]
  }
  if (length(expr_rows) == nrow(GENES)) break
}
close(con)
stopifnot(length(expr_rows) == nrow(GENES))

gt2d <- function(v) { d <- integer(length(v))
  d[v %in% c("0|1", "1|0", "0/1", "1/0")] <- 1L
  d[v %in% c("1|1", "1/1")] <- 2L; d }
ord <- readLines(SAMP)

for (i in seq_len(nrow(GENES))) { tryCatch({
  g   <- GENES[i, ]; tag <- sub("\\..*$", "", g$gene)
  ck  <- file.path(OUT_DIR, paste0(tag, ".rds"))
  p   <- read_retry(ck)
  if ("BayesR_2k" %in% names(p$indices)) { cat("[cached]", tag, "\n"); next }

  sp  <- strsplit(readLines(file.path(LOC_DIR, paste0(tag, ".gt"))),
                  "\t", fixed = TRUE)
  pos <- as.integer(vapply(sp, `[`, "", 1L))
  G   <- vapply(sp, function(r) gt2d(r[-1L]), integer(length(ord)))
  maf <- colMeans(G) / 2; maf <- pmin(maf, 1 - maf)
  keep <- maf >= 0.05 & apply(G, 2, sd) > 0
  G <- G[, keep, drop = FALSE]; pos <- pos[keep]
  j_lead <- match(g$lead_pos, pos)
  stopifnot(identical(p$j_lead, j_lead))

  y_row <- expr_rows[[g$gene]]
  expr  <- as.numeric(y_row[-(1:4)]); names(expr) <- hdr[-(1:4)]
  y <- as.numeric(scale(as.numeric(expr[ord])))
  X <- scale(G); X[!is.finite(X)] <- 0

  t0  <- Sys.time()
  pip <- bayesr_long(y, X, seed = 20260908L + i)
  el  <- as.numeric(Sys.time() - t0, units = "secs")
  sel <- which(pip > 0.99)
  g_lead <- if (!is.na(j_lead)) G[, j_lead] else NULL
  r2b <- if (length(sel) > 0 && !is.null(g_lead))
           max(vapply(sel, function(s) cor(G[, s], g_lead)^2, numeric(1)))
         else NA_real_

  summ <- p$summary
  old_row <- summ[summ$method == "BayesR", ]
  old_row$method <- "BayesR_2k"
  summ[summ$method == "BayesR",
       c("K_hat", "r2_lead", "elapsed")] <- list(length(sel), r2b, el)
  summ <- rbind(summ, old_row)
  p$summary <- summ
  p$indices$BayesR_2k <- p$indices$BayesR
  p$indices$BayesR <- sel
  tmp <- paste0(ck, ".tmp"); saveRDS(p, tmp); file.rename(tmp, ck)
  cat(sprintf("[done] %s BayesR long: K=%d r2_lead=%s (%.0fs)\n",
              tag, length(sel),
              ifelse(is.na(r2b), "NA", sprintf("%.3f", r2b)), el))
}, error = function(e) message("FAIL gene ", i, ": ", conditionMessage(e))) }
fs <- list.files(OUT_DIR, pattern = "^ENSG.*rds$", full.names = TRUE)
summ <- do.call(rbind, lapply(fs, function(f) readRDS(f)$summary))
write.csv(summ, file.path(OUT_DIR, "geuvadis_summary.csv"), row.names = FALSE)

## ── F2. mouse BMI (autosomes) ────────────────────────────────────────────────
ckm <- "results/bench_full_exact/12_mouse_autosomes/bayesr_long.rds"
if (!file.exists(ckm)) {
  library(BGLR); data(mice)
  autosome_idx <- which(mice.map$chr != "X")
  X_aut <- scale(mice.X[, autosome_idx]); X_aut[!is.finite(X_aut)] <- 0
  y <- as.numeric(scale(mice.pheno$Obesity.BMI))
  t0 <- Sys.time()
  pip <- bayesr_long(y, X_aut, seed = 20260908L)
  el  <- as.numeric(Sys.time() - t0, units = "secs")
  sel <- which(pip > 0.99)
  saveRDS(list(indices = sel, K_hat = length(sel),
               pip_top = sort(pip, decreasing = TRUE)[1:10],
               elapsed = el, niter = NITER, nburn = NBURN), ckm)
  cat(sprintf("[done] mouse BayesR long: K=%d (%.0fs); top PIPs: %s\n",
              length(sel), el,
              paste(sprintf("%.2f", sort(pip, decreasing = TRUE)[1:5]),
                    collapse = " ")))
} else cat("[cached] mouse\n")
message("done 52.")
