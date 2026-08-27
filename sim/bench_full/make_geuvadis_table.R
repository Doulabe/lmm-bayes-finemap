# ==============================================================================
# make_geuvadis_table.R
# LaTeX table for the GEUVADIS cis-eQTL illustration (26_geuvadis_eqtl.R).
# Rows = genes; per method: K_hat and r^2 of the best selection to the
# official EUR373 lead eQTL.
# Output: <out_dir>/tab_geuvadis_eqtl.tex
# Usage:  Rscript sim/bench_full/make_geuvadis_table.R [--out_dir results/bench_full/tables_cbf]
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
i <- match("--out_dir", args)
OUT <- if (!is.na(i)) args[i + 1L] else "results/bench_full/tables_cbf"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

s <- read.csv("results/bench_full/26_geuvadis/geuvadis_summary.csv",
              stringsAsFactors = FALSE)
LAB <- c(MS_L_eBIC = "\\textsc{CBF-LMM}", SuSiE = "SuSiE", BSLMM = "BSLMM",
         BayesR = "BayesR", fastlmm = "FaST-LMM")
MORD <- names(LAB)
GORD <- c("ERAP2", "RPS26", "SLFN5", "SNHG5", "FLVCR1-AS1(1q32)",
          "PEX6-region", "TRA2A-AS(7p15)", "ZNF266")
GORD <- GORD[GORD %in% s$gene]

cell <- function(g, m) {
  r <- s[s$gene == g & s$method == m, ]
  if (nrow(r) == 0) return("---")
  if (r$K_hat == 0) return("$0$ (---)")
  sprintf("$%d$ ($%.2f$)", r$K_hat, r$r2_lead)
}
con <- file(file.path(OUT, "tab_geuvadis_eqtl.tex"), "w")
wl <- function(...) writeLines(paste0(...), con)
wl("\\begin{tabular}{l l r ", paste(rep("c", length(MORD)), collapse = " "), "}")
wl("\\toprule")
wl(" & & & \\multicolumn{", length(MORD),
   "}{c}{$\\hat K$ ($r^2$ of best selection to the official lead)} \\\\")
wl("\\cmidrule(lr){4-", 3 + length(MORD), "}")
wl("Gene & Lead eQTL & $m$ & ", paste(LAB[MORD], collapse = " & "), " \\\\")
wl("\\midrule")
for (g in GORD) {
  r0 <- s[s$gene == g, ][1, ]
  gname <- sub("\\(.*$", "", g)   # short display name
  wl(sprintf("%s & %s & %d & %s \\\\", gname, r0$lead, r0$m,
             paste(vapply(MORD, function(m) cell(g, m), ""), collapse = " & ")))
}
wl("\\bottomrule"); wl("\\end{tabular}")
close(con)
cat("Wrote", file.path(OUT, "tab_geuvadis_eqtl.tex"), "\n")
