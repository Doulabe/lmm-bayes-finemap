# ==============================================================================
# make_threshold_strat_tables.R
# Axis-stratified Top-K* recall tables for the threshold-independent
# comparison (Supplementary Note S6):
#   tab_threshold_indep_by_rho.tex   rho in {0.80, 0.95, 0.98}
#   tab_threshold_indep_by_n.tex     n   in {500, 1000, 3000}
#   tab_threshold_indep_by_arch.tex  weak / medium / strong
#   tab_threshold_indep_by_m.tex     m   in {5000, 10000}
# The anchor cell (n=1000, rho=0.95, medium, m=5000) supplies the middle
# level of every axis.  Entries are mean Top-K* recall +/- Monte Carlo SE
# over the B=100 replicates of each cell.
#
# Input : results/bench_full/19_threshold_full_grid/threshold_full_grid_long.csv
# Output: <out_dir>/tab_threshold_indep_by_{rho,n,arch,m}.tex
#
# Usage (from project root):
#   Rscript sim/bench_full/make_threshold_strat_tables.R [--out_dir results/bench_full/tables]
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
i <- match("--out_dir", args)
OUT <- if (!is.na(i)) args[i + 1L] else "results/bench_full/tables"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

d <- read.csv("results/bench_full/19_threshold_full_grid/threshold_full_grid_long.csv",
              stringsAsFactors = FALSE)
d$base_cell <- sub("_sg[0-9.]+$", "", d$cell)

LAB <- c(JS_L_eBIC = "\\textsc{JS-LMM}-eBIC", MS_L_eBIC = "\\textsc{MS-LMM}-eBIC",
         BSLMM = "BSLMM", BayesR = "BayesR", SuSiE = "SuSiE",
         fastlmm = "FaST-LMM")
ORD <- names(LAB)

cell_stat <- function(base_cell, sg, method) {
  v <- d$topK[d$base_cell == base_cell & d$sigma_g2 == sg & d$method == method]
  if (!length(v)) return("---")
  sprintf("$%.3f \\pm %.3f$", mean(v), sd(v) / sqrt(length(v)))
}

make_table <- function(file, header, levels) {
  # levels: named list  display_label -> base_cell tag
  con <- file(file.path(OUT, file), "w")
  wl <- function(...) writeLines(paste0(...), con)
  wl("\\begin{tabular}{l l cc}"); wl("\\toprule")
  wl(header, " & Method & $\\sigma_g^2=0$ & $\\sigma_g^2=0.5$ \\\\")
  wl("\\midrule")
  first <- TRUE
  for (lv in names(levels)) {
    if (!first) wl("\\midrule")
    first <- FALSE
    for (m in ORD)
      wl(sprintf("%s & %s & %s & %s \\\\",
                 lv, LAB[m],
                 cell_stat(levels[[lv]], 0,   m),
                 cell_stat(levels[[lv]], 0.5, m)))
  }
  wl("\\bottomrule"); wl("\\end{tabular}")
  close(con)
  cat("Wrote", file.path(OUT, file), "\n")
}

make_table("tab_threshold_indep_by_rho.tex", "$\\rho$",
           list("0.80" = "rho080", "0.95" = "anchor", "0.98" = "rho098"))
make_table("tab_threshold_indep_by_n.tex", "$n$",
           list("500" = "n500", "1000" = "anchor", "3000" = "n3000"))
make_table("tab_threshold_indep_by_arch.tex", "Signal",
           list("weak" = "sigweak", "medium" = "anchor", "strong" = "sigstrong"))
make_table("tab_threshold_indep_by_m.tex", "$m$",
           list("5{,}000" = "anchor", "10{,}000" = "m10000"))
