# ==============================================================================
# make_semisynth_tables.R
# LaTeX tables for the semi-synthetic 1000G benchmark (25_semisynth_1000g.R):
#   - tab_semisynth_1000g.tex          pooled 8-method x 2-regime table
#                                      (main text, "second design" section)
#   - tab_semisynth_1000g_bylocus.tex  per-locus mean F1 (supplement)
#
# Input : results/bench_full/25_semisynth_1000g/semisynth_1000g_raw.csv
# Output: <out_dir>/tab_semisynth_1000g{,_bylocus}.tex
#
# Usage (from project root):
#   Rscript sim/bench_full/make_semisynth_tables.R [--out_dir results/bench_full/tables]
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
i <- match("--out_dir", args)
OUT <- if (!is.na(i)) args[i + 1L] else "results/bench_full/tables"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

raw <- read.csv("results/bench_full/25_semisynth_1000g/semisynth_1000g_raw.csv")

LAB <- c(JS_L_eBIC = "\\textsc{JS-LMM}-eBIC", MS_L_eBIC = "\\textsc{MS-LMM}-eBIC",
         JS_L_JP99 = "\\textsc{JS-LMM}-JP99", MS_L_JP99 = "\\textsc{MS-LMM}-JP99",
         BSLMM = "BSLMM", BayesR = "BayesR", SuSiE = "SuSiE", fastlmm = "FaST-LMM")
ORD <- names(LAB)
se <- function(v) sd(v) / sqrt(length(v))

## ---- pooled table (main text) ----------------------------------------------
panel <- function(sg) {
  d <- raw[raw$sigma_g2 == sg, ]
  vapply(ORD, function(m) {
    s <- d[d$method == m, ]
    sprintf("%s & $%.3f \\pm %.3f$ & $%.3f$ & $%.3f$ & $%.1f$ \\\\",
            LAB[m], mean(s$f1), se(s$f1), mean(s$recall),
            mean(s$precision), mean(s$K_hat))
  }, character(1))
}
con <- file(file.path(OUT, "tab_semisynth_1000g.tex"), "w")
wl <- function(...) writeLines(paste0(...), con)
wl("\\begin{tabular}{l rrrr}"); wl("\\toprule")
wl("\\multicolumn{5}{c}{$\\sigma_g^2 = 0$ (homogeneous noise)}\\\\"); wl("\\midrule")
wl("Method & $F_1 \\pm$ SE & Recall & Precision & $\\widehat K$ \\\\"); wl("\\midrule")
writeLines(panel(0), con); wl("\\midrule")
wl("\\multicolumn{5}{c}{$\\sigma_g^2 = 0.5$ (polygenic background)}\\\\"); wl("\\midrule")
wl("Method & $F_1 \\pm$ SE & Recall & Precision & $\\widehat K$ \\\\"); wl("\\midrule")
writeLines(panel(0.5), con); wl("\\bottomrule"); wl("\\end{tabular}")
close(con)

## ---- per-locus F1 table (supplement) ---------------------------------------
f1m <- function(loc, sg) {
  d <- raw[raw$locus == loc & raw$sigma_g2 == sg, ]
  setNames(vapply(ORD, function(m) mean(d$f1[d$method == m]), numeric(1)), ORD)
}
c1a <- f1m("chr1", 0); c1b <- f1m("chr1", 0.5)
c6a <- f1m("chr6", 0); c6b <- f1m("chr6", 0.5)
con <- file(file.path(OUT, "tab_semisynth_1000g_bylocus.tex"), "w")
wl <- function(...) writeLines(paste0(...), con)
wl("\\begin{tabular}{l rr rr}"); wl("\\toprule")
wl(" & \\multicolumn{2}{c}{chr1 ($m=1{,}493$)} & \\multicolumn{2}{c}{chr6 ($m=2{,}143$)} \\\\")
wl("\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}")
wl("Method & $\\sigma_g^2{=}0$ & $\\sigma_g^2{=}0.5$ & $\\sigma_g^2{=}0$ & $\\sigma_g^2{=}0.5$ \\\\")
wl("\\midrule")
for (m in ORD)
  wl(sprintf("%s & $%.3f$ & $%.3f$ & $%.3f$ & $%.3f$ \\\\",
             LAB[m], c1a[m], c1b[m], c6a[m], c6b[m]))
wl("\\bottomrule"); wl("\\end{tabular}")
close(con)

cat("Wrote", file.path(OUT, "tab_semisynth_1000g.tex"), "and _bylocus.tex\n")
