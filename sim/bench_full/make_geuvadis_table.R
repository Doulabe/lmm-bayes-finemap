# ==============================================================================
# make_geuvadis_table.R
# LaTeX table for the GEUVADIS cis-eQTL illustration (26_geuvadis_eqtl.R),
# annotated with the identity of each method's best selected variant.
#
# Cell format: rsID of the best selection (highest r^2 to the official lead);
# bold when r^2 = 1.00; otherwise the r^2 follows in parentheses; "+k" marks
# k additional selected variants; --- marks an empty selection.
#
# rsIDs come from data/geuvadis/pos2rs.csv (position -> dbSNP rsID on GRCh37,
# SNV-disambiguated via the Ensembl GRCh37 REST API; official lead rsIDs are
# reproduced exactly by this annotation).
#
# Input : results/bench_full/26_geuvadis/{<gene>.rds, geuvadis_summary.csv}
#         data/geuvadis/{pos2rs.csv, loci/<gene>.gt}
# Output: <out_dir>/tab_geuvadis_eqtl.tex
# Usage : Rscript sim/bench_full/make_geuvadis_table.R [--out_dir results/bench_full/tables_cbf]
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
i <- match("--out_dir", args)
OUT <- if (!is.na(i)) args[i + 1L] else "results/bench_full/tables_cbf"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

RES <- "results/bench_full/26_geuvadis"
D   <- "data/geuvadis"
p2r_file <- if (file.exists(file.path(D, "pos2rs.csv"))) file.path(D, "pos2rs.csv") else "sim/bench_full/geuvadis_pos2rs.csv"
p2r <- read.csv(p2r_file, header = FALSE,
                col.names = c("chr", "pos", "rs"), stringsAsFactors = FALSE)
s   <- read.csv(file.path(RES, "geuvadis_summary.csv"), stringsAsFactors = FALSE)

GENES <- data.frame(
  tag    = c("ENSG00000164308", "ENSG00000197728", "ENSG00000166750",
             "ENSG00000203875", "ENSG00000198468", "ENSG00000124587",
             "ENSG00000230658", "ENSG00000174652"),
  gene   = c("ERAP2", "RPS26", "SLFN5", "SNHG5", "FLVCR1-AS1(1q32)",
             "PEX6-region", "TRA2A-AS(7p15)", "ZNF266"),
  lead_pos = c(96252589L, 56401085L, 33571546L, 86387888L,
               213049214L, 42944850L, 23143113L, 9544276L),
  stringsAsFactors = FALSE)

LAB  <- c(MS_L_eBIC = "\\textsc{CBF-LMM}", SuSiE = "SuSiE", BSLMM = "BSLMM",
          BayesR = "BayesR", fastlmm = "FaST-LMM")
MORD <- names(LAB)
gt2d <- function(v) { d <- integer(length(v))
  d[v %in% c("0|1", "1|0", "0/1", "1/0")] <- 1L
  d[v %in% c("1|1", "1/1")] <- 2L; d }

cells <- list()
for (gi in seq_len(nrow(GENES))) {
  g <- GENES[gi, ]
  r <- readRDS(file.path(RES, paste0(g$tag, ".rds")))
  sp  <- strsplit(readLines(file.path(D, "loci", paste0(g$tag, ".gt"))),
                  "\t", fixed = TRUE)
  pos_all <- as.integer(vapply(sp, `[`, "", 1L))
  getG <- function(p) gt2d(sp[[match(p, pos_all)]][-1L])
  g_lead <- getG(g$lead_pos)
  for (m in MORD) {
    idx <- r$indices[[m]]
    key <- paste(g$gene, m, sep = "|")
    if (length(idx) == 0L) { cells[[key]] <- "---"; next }
    r2 <- vapply(idx, function(sx) cor(getG(r$pos[sx]), g_lead)^2, numeric(1))
    b  <- idx[which.max(r2)]; r2b <- max(r2)
    rs <- p2r$rs[match(r$pos[b], p2r$pos)]
    if (is.na(rs)) rs <- sprintf("chr?:%d", r$pos[b])
    extra <- if (length(idx) > 1L) sprintf("\\,{+%d}", length(idx) - 1L) else ""
    cells[[key]] <- if (r2b >= 0.995) sprintf("\\textbf{%s}%s", rs, extra)
                    else sprintf("%s%s (%.2f)", rs, extra, r2b)
  }
}

con <- file(file.path(OUT, "tab_geuvadis_eqtl.tex"), "w")
wl <- function(...) writeLines(paste0(...), con)
wl("\\begin{tabular}{l l ", paste(rep("l", length(MORD)), collapse = " "), "}")
wl("\\toprule")
wl(" & & \\multicolumn{", length(MORD),
   "}{c}{Best selected variant (bold: tags the official lead at $r^2 = 1.00$)} \\\\")
wl("\\cmidrule(lr){3-", 2 + length(MORD), "}")
wl("Gene & Lead eQTL & ", paste(LAB[MORD], collapse = " & "), " \\\\")
wl("\\midrule")
for (gi in seq_len(nrow(GENES))) {
  g <- GENES[gi, ]
  r0 <- s[s$gene == g$gene, ][1, ]
  wl(sprintf("%s & %s & %s \\\\", sub("\\(.*$", "", g$gene), r0$lead,
             paste(vapply(MORD, function(m) cells[[paste(g$gene, m, sep = "|")]],
                          ""), collapse = " & ")))
}
wl("\\bottomrule"); wl("\\end{tabular}")
close(con)
cat("Wrote", file.path(OUT, "tab_geuvadis_eqtl.tex"), "\n")
