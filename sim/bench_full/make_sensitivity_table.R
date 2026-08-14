# ============================================================
# make_sensitivity_table.R
# Convert results/bench_full/07_prior_sensitivity_summary.csv to
# theory/tab_supp_sensitivity_ig.tex (LaTeX subtable for Appendix C.2).
# ============================================================

suppressPackageStartupMessages(library(dplyr))

IN  <- "results/bench_full/07_prior_sensitivity_summary.csv"
dir.create("theory", showWarnings = FALSE)
OUT <- "theory/tab_supp_sensitivity_ig.tex"

if (!file.exists(IN)) {
  message("(skipped: ", IN, " not found)"); quit("no", status = 0L)
}
df <- read.csv(IN)
df$method <- factor(df$method, levels = c("MS_L_eBIC", "JS_L_eBIC"),
                     labels = c("\\textsc{MS-LMM}", "\\textsc{JS-LMM}"))
df$cell <- factor(df$cell, levels = unique(df$cell))

fmt2 <- function(x, d = 3) {
  if (is.na(x)) "---" else formatC(x, digits = d, format = "f")
}

rows <- c()
prev_cell <- ""
for (i in seq_len(nrow(df))) {
  r <- df[i, ]
  cell_lbl <- if (as.character(r$cell) != prev_cell) {
    prev_cell <- as.character(r$cell)
    sprintf("\\multirow{2}{*}{%s}", prev_cell)
  } else ""
  rows <- c(rows, sprintf(
    "%s & %s & %s & %s & %s & %s & %s \\\\",
    cell_lbl,
    as.character(r$method),
    fmt2(r$F1_hc),
    fmt2(r$F1_ig),
    sprintf("%+.3f", r$delta_F1),
    fmt2(r$jaccard, 2),
    fmt2(r$max_dlogBF, 2)))
  if (i %% 2L == 0L && i < nrow(df)) rows <- c(rows, "\\midrule")
}

body <- c(
  "\\begin{tabular}{l l rr r rr}",
  "\\toprule",
  "Cell & Method & $F_1^{\\mathrm{HC}}$ & $F_1^{\\mathrm{IG}}$ & ",
  "  $\\Delta F_1$ & Jaccard & $\\max |\\Delta\\log\\mathrm{BF}|$ \\\\",
  "\\midrule",
  rows,
  "\\bottomrule",
  "\\end{tabular}"
)

writeLines(body, OUT)
message("Wrote ", OUT)

# Summary stats for paper text
cat("\n--- Sensitivity summary (over ", nrow(df), " entries) ---\n",
    "Mean |Delta F1|       = ", formatC(mean(abs(df$delta_F1)), 4, format="f"), "\n",
    "Max  |Delta F1|       = ", formatC(max(abs(df$delta_F1)),  4, format="f"), "\n",
    "Mean Jaccard          = ", formatC(mean(df$jaccard),       3, format="f"), "\n",
    "Min  Jaccard          = ", formatC(min(df$jaccard),        3, format="f"), "\n",
    "Mean max|dlogBF|      = ", formatC(mean(df$max_dlogBF, na.rm=TRUE), 3, format="f"), "\n",
    sep="")
