# ==============================================================================
# 34_geuvadis_exact.R
# Recompute the CBF-LMM rows of the GEUVADIS cis-eQTL analysis with the exact
# candidate-specific engine (K_jk + profile-ML eBIC). Uses the cached .gt
# locus files and the local expression matrix; comparator rows are cloned from
# the original per-gene RDS (results/bench_full/26_geuvadis/), with the old
# shared-kernel row kept as MS_L_eBIC_sharedGk.
# Output: results/bench_full_exact/26_geuvadis/<gene>.rds + geuvadis_summary.csv
# Usage:  Rscript sim/bench_full/34_geuvadis_exact.R
# ==============================================================================

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
source("sim/bench_full/00_config.R")
source("R/CBF_LMM_exact.R")

D       <- "data/geuvadis"
LOC_DIR <- file.path(D, "loci")
IN_DIR  <- "results/bench_full/26_geuvadis"
OUT_DIR <- "results/bench_full_exact/26_geuvadis"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
SAMP    <- file.path(D, "geuvadis_eur_overlap.txt")
EXPR_GZ <- file.path(D, "GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz")

GENES <- data.frame(
  gene   = c("ENSG00000164308.12", "ENSG00000197728.5", "ENSG00000166750.4",
             "ENSG00000203875.4",  "ENSG00000198468.2", "ENSG00000124587.9",
             "ENSG00000230658.1",  "ENSG00000174652.12"),
  symbol = c("ERAP2", "RPS26", "SLFN5", "SNHG5", "FLVCR1-AS1(1q32)",
             "PEX6-region", "TRA2A-AS(7p15)", "ZNF266"),
  lead_rs  = c("rs2910686", "rs10876864", "rs11080327", "rs1059307",
               "rs12123978", "rs6907751", "rs10233039", "rs10420709"),
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

ord <- readLines(SAMP)   # sample order of the cached .gt files (n = 358)

for (i in seq_len(nrow(GENES))) {
  g   <- GENES[i, ]; tag <- sub("\\..*$", "", g$gene)
  ck  <- file.path(OUT_DIR, paste0(tag, ".rds"))
  if (file.exists(ck)) { cat("[cached]", g$symbol, "\n"); next }
  gt  <- file.path(LOC_DIR, paste0(tag, ".gt"))
  sp  <- strsplit(readLines(gt), "\t", fixed = TRUE)
  pos <- as.integer(vapply(sp, `[`, "", 1L))
  G   <- vapply(sp, function(r) gt2d(r[-1L]), integer(length(ord)))
  stopifnot(nrow(G) == length(ord))
  maf <- colMeans(G) / 2; maf <- pmin(maf, 1 - maf)
  keep <- maf >= 0.05 & apply(G, 2, sd) > 0
  G <- G[, keep, drop = FALSE]; pos <- pos[keep]
  j_lead <- match(g$lead_pos, pos)

  y_row <- expr_rows[[g$gene]]
  expr  <- as.numeric(y_row[-(1:4)]); names(expr) <- hdr[-(1:4)]
  y <- as.numeric(scale(as.numeric(expr[ord])))
  X <- scale(G); X[!is.finite(X)] <- 0

  orig <- readRDS(file.path(IN_DIR, paste0(tag, ".rds")))
  stopifnot(identical(orig$j_lead, j_lead), nrow(G) == orig$summary$n[1])

  t0  <- Sys.time()
  fit <- CBF_LMM_stepwise_exact(y, X, tau2 = ANCHOR$tau2,
                                  K_max = ANCHOR$K_max,
                                  n_nodes = ANCHOR$N_delta)
  el  <- as.numeric(Sys.time() - t0, units = "secs")

  g_lead <- if (!is.na(j_lead)) G[, j_lead] else NULL
  r2b <- if (fit$K_hat > 0 && !is.null(g_lead))
           max(vapply(fit$indices, function(s) cor(G[, s], g_lead)^2,
                      numeric(1)))
         else NA_real_

  summ <- orig$summary
  old_row <- summ[summ$method == "MS_L_eBIC", ]
  old_row$method <- "MS_L_eBIC_sharedGk"
  summ[summ$method == "MS_L_eBIC",
       c("K_hat", "r2_lead", "elapsed")] <- list(fit$K_hat, r2b, el)
  summ <- rbind(summ, old_row)

  idx <- orig$indices
  idx$MS_L_eBIC_sharedGk <- idx$MS_L_eBIC
  idx$MS_L_eBIC <- fit$indices

  saveRDS(list(summary = summ, indices = idx, pos = pos, j_lead = j_lead), ck)
  cat(sprintf("[done] %-16s exact: K=%d r2_lead=%.3f (%.1fs) | shared: K=%d r2=%.3f\n",
              g$symbol, fit$K_hat, r2b, el,
              old_row$K_hat, old_row$r2_lead))
}

fs <- list.files(OUT_DIR, pattern = "^ENSG.*rds$", full.names = TRUE)
summ <- do.call(rbind, lapply(fs, function(f) readRDS(f)$summary))
write.csv(summ, file.path(OUT_DIR, "geuvadis_summary.csv"), row.names = FALSE)
cat("\n=== GEUVADIS exact vs shared ===\n")
print(summ[summ$method %in% c("MS_L_eBIC", "MS_L_eBIC_sharedGk"),
           c("gene", "method", "K_hat", "r2_lead")], row.names = FALSE)
