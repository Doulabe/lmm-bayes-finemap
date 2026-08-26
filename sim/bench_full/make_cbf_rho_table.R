x <- read.csv("results/bench_full/02_scaling_rho_summary.csv")
keep <- c(MS_L_eBIC="\\textsc{CBF-LMM}", BSLMM="BSLMM", BayesR="BayesR",
          SuSiE="SuSiE", fastlmm="FaST-LMM")
dir.create("results/bench_full/tables_cbf", recursive = TRUE, showWarnings = FALSE)
con <- file("results/bench_full/tables_cbf/tab_supp_by_rho_primary.tex", "w")
wl <- function(...) writeLines(paste0(...), con)
wl("\\begin{tabular}{lrrrr}"); wl("\\toprule")
wl("Method & $\\overline{F_1}$ & Recall & Precision & $\\overline{\\hat K}$ \\\\")
wl("\\midrule")
first <- TRUE
for (r in c(0.80, 0.95, 0.98)) for (sg in c(0, 0.5)) {
  if (!first) wl("\\midrule"); first <- FALSE
  wl(sprintf("\\multicolumn{5}{c}{\\textbf{$\\rho=%.2f$, $\\sigma_g^2=%s$}} \\\\", r, sg))
  for (m in names(keep)) {
    s <- x[x$rho==r & x$sigma_g2==sg & x$method==m, ]
    wl(sprintf("%s & %.2f & %.2f & %.2f & %.2f \\\\",
               keep[m], s$mean_f1, s$mean_recall, s$mean_precision, s$mean_K_hat))
  }
}
wl("\\bottomrule"); wl("\\end{tabular}"); close(con)
cat("wrote tab_supp_by_rho_primary.tex\n")
