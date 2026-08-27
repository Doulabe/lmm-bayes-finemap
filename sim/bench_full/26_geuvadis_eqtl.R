# ==============================================================================
# 26_geuvadis_eqtl.R
# Human cis-eQTL fine-mapping illustration on GEUVADIS (public data).
#
# For each of 8 genes with strong, replicated cis-eQTLs (official EUR373
# FDR5 best-association list), fine-map the cis-window (TSS +/- 500 kb,
# extended to cover the official lead) in the n = 358 individuals present in
# both GEUVADIS and the 1000G phase-3 EUR panel:
#   y = PEER-normalized gene expression (GD462 resk10), standardized
#   X = phase-3 genotype dosages (0/1/2), biallelic SNPs, MAF >= 0.05
#   ground truth = the official GEUVADIS EUR373 best cis-eQTL (by position)
# Methods (paper settings): CBF-LMM (= MS_L_eBIC, marginalized delta + eBIC),
# SuSiE (PIP > 0.99), BSLMM (PIP > 0.5), BayesR (PIP > 0.99),
# FaST-LMM (Bonferroni 0.05).
#
# Inputs (see README / data/geuvadis):
#   data/geuvadis/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz
#   data/geuvadis/EUR373.gene.cis.FDR5.best.rs137.txt.gz
#   data/geuvadis/geuvadis_eur_overlap.txt      (358 sample IDs)
# Genotypes are streamed remotely from the 1000G phase-3 FTP via bcftools
# and cached per locus under data/geuvadis/loci/.
#
# Output: results/bench_full/26_geuvadis/<gene>.rds + geuvadis_summary.csv
# Usage:  Rscript sim/bench_full/26_geuvadis_eqtl.R
# ==============================================================================

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
source("sim/bench_full/00_config.R")

D       <- "data/geuvadis"
LOC_DIR <- file.path(D, "loci");  dir.create(LOC_DIR, showWarnings = FALSE)
OUT_DIR <- "results/bench_full/26_geuvadis"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
SAMP    <- file.path(D, "geuvadis_eur_overlap.txt")
EXPR_GZ <- file.path(D, "GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz")
URLB    <- "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
W       <- 500000L
METHODS <- c("MS_L_eBIC", "SuSiE", "BayesR", "fastlmm", "BSLMM")

# 8 strongest rs-lead autosomal eQTLs from EUR373.gene.cis.FDR5.best.rs137
GENES <- data.frame(
  gene   = c("ENSG00000164308.12", "ENSG00000197728.5", "ENSG00000166750.4",
             "ENSG00000203875.4",  "ENSG00000198468.2", "ENSG00000124587.9",
             "ENSG00000230658.1",  "ENSG00000174652.12"),
  symbol = c("ERAP2", "RPS26", "SLFN5", "SNHG5", "FLVCR1-AS1(1q32)",
             "PEX6-region", "TRA2A-AS(7p15)", "ZNF266"),
  chr      = c("5", "12", "17", "6", "1", "6", "7", "19"),
  lead_rs  = c("rs2910686", "rs10876864", "rs11080327", "rs1059307",
               "rs12123978", "rs6907751", "rs10233039", "rs10420709"),
  lead_pos = c(96252589L, 56401085L, 33571546L, 86387888L,
               213049214L, 42944850L, 23143113L, 9544276L),
  stringsAsFactors = FALSE)

## ---- expression: load once, all 8 rows --------------------------------------
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
tss <- vapply(expr_rows, function(r) as.integer(r[4]), integer(1))[GENES$gene]
cat("Expression rows loaded for", length(expr_rows), "genes.\n")

gt2d <- function(v) { d <- integer(length(v))
  d[v %in% c("0|1", "1|0", "0/1", "1/0")] <- 1L
  d[v %in% c("1|1", "1/1")] <- 2L; d }

stream_locus <- function(chr, lo, hi, tag) {
  gt <- file.path(LOC_DIR, paste0(tag, ".gt"))
  url <- sprintf(URLB, chr)
  for (try in 1:3) {
    if (file.exists(gt) && file.size(gt) > 0) break
    cmd <- sprintf(paste0(
      "bcftools view -r %s:%d-%d -S %s --force-samples -m2 -M2 -v snps %s 2>/dev/null",
      " | bcftools query -f '%%POS[\\t%%GT]\\n' > %s 2>/dev/null"),
      chr, lo, hi, shQuote(SAMP), shQuote(url), shQuote(gt))
    system(cmd)
    if (!file.exists(gt) || file.size(gt) == 0) Sys.sleep(10)
  }
  if (!file.exists(gt) || file.size(gt) == 0) return(NULL)
  gt
}

run_gene <- function(i) {
  g <- GENES[i, ]; tag <- sub("\\..*$", "", g$gene)
  ck <- file.path(OUT_DIR, paste0(tag, ".rds"))
  if (file.exists(ck)) { cat("[cached]", g$symbol, "\n"); return(readRDS(ck)) }
  t_all <- Sys.time()
  lo <- min(tss[g$gene] - W, g$lead_pos - 100000L)
  hi <- max(tss[g$gene] + W, g$lead_pos + 100000L)
  gt <- stream_locus(g$chr, lo, hi, tag)
  if (is.null(gt)) { cat("[STREAM FAIL]", g$symbol, "\n"); return(NULL) }

  ord <- system(sprintf("bcftools query -l -S %s --force-samples %s 2>/dev/null",
                        shQuote(SAMP), shQuote(sprintf(URLB, g$chr))), intern = TRUE)
  sp  <- strsplit(readLines(gt), "\t", fixed = TRUE)
  pos <- as.integer(vapply(sp, `[`, "", 1L))
  G   <- vapply(sp, function(r) gt2d(r[-1L]), integer(length(ord)))
  maf <- colMeans(G) / 2; maf <- pmin(maf, 1 - maf)
  keep <- maf >= 0.05 & matrixStats::colSds(G) > 0
  if (!requireNamespace("matrixStats", quietly = TRUE))
    keep <- maf >= 0.05 & apply(G, 2, sd) > 0
  G <- G[, keep, drop = FALSE]; pos <- pos[keep]
  j_lead <- match(g$lead_pos, pos)
  y_row <- expr_rows[[g$gene]]
  expr  <- as.numeric(y_row[-(1:4)]); names(expr) <- hdr[-(1:4)]
  y <- as.numeric(scale(as.numeric(expr[ord])))
  X <- scale(G); X[!is.finite(X)] <- 0

  d <- list(X = X, y = y, truth = if (is.na(j_lead)) integer(0) else j_lead,
            sigma_g2 = NA)
  res <- run_methods(d, METHODS, K_max = ANCHOR$K_max, theta = ANCHOR$theta,
                     N_delta = ANCHOR$N_delta, tau2 = ANCHOR$tau2)

  g_lead <- if (!is.na(j_lead)) G[, j_lead] else NULL
  rows <- do.call(rbind, lapply(names(res), function(mn) {
    r <- res[[mn]]
    r2b <- if (r$K_hat > 0 && !is.null(g_lead))
             max(vapply(r$indices, function(s) cor(G[, s], g_lead)^2, numeric(1)))
           else NA_real_
    data.frame(gene = g$symbol, lead = g$lead_rs, method = mn,
               K_hat = r$K_hat, r2_lead = r2b,
               n = nrow(G), m = ncol(G), elapsed = r$elapsed)
  }))
  out <- list(summary = rows, indices = lapply(res, `[[`, "indices"),
              pos = pos, j_lead = j_lead)
  saveRDS(out, ck)
  cat(sprintf("[done] %-16s m=%d, lead present=%s, %.1f min\n", g$symbol,
              ncol(G), !is.na(j_lead),
              as.numeric(Sys.time() - t_all, units = "mins")))
  out
}

ALL <- lapply(seq_len(nrow(GENES)), function(i)
  tryCatch(run_gene(i), error = function(e) {
    cat("[ERROR]", GENES$symbol[i], ":", conditionMessage(e), "\n"); NULL }))
ALL <- ALL[!vapply(ALL, is.null, logical(1))]
summ <- do.call(rbind, lapply(ALL, `[[`, "summary"))
write.csv(summ, file.path(OUT_DIR, "geuvadis_summary.csv"), row.names = FALSE)
cat("\n=== SUMMARY (r2 of best selection to official lead) ===\n")
print(summ[order(summ$gene, summ$method),
           c("gene", "method", "K_hat", "r2_lead", "m")], row.names = FALSE)
