# ==============================================================================
# 48_make_reviewer_tables.R
# Manuscript tables for the correlated-causal stress test (45) and the
# PIP-threshold sensitivity (43).
# The summaries are rebuilt here from the per-replicate checkpoints (and the
# summary CSVs refreshed), so in-place payload updates such as the BayesR
# long-budget rerun (51) can never leave these tables stale.
# ==============================================================================
D45 <- "results/bench_full_exact/45_correlated_causals"
res45 <- do.call(rbind, lapply(
  list.files(D45, pattern = "^rho.*rds$", full.names = TRUE),
  function(f) readRDS(f)$rows))
IN45 <- aggregate(cbind(K_hat, tp, precision, recall, f1) ~
                    method + rho + sigma_g2, res45, mean)
write.csv(IN45, file.path(D45, "correlated_causals_summary.csv"),
          row.names = FALSE)

D43 <- "results/bench_full_exact/43_pip_thresholds"
res43 <- do.call(rbind, lapply(
  list.files(D43, pattern = "^sg.*rds$", full.names = TRUE), readRDS))
IN43 <- aggregate(cbind(K_hat, precision, recall, f1) ~
                    method + threshold + sigma_g2, res43, mean)
write.csv(IN43, file.path(D43, "pip_threshold_summary.csv"),
          row.names = FALSE)
OUT  <- "CBF_LMM_restructured/tables"
DISP <- c(MS_L_eBIC = "\\textsc{CBF-LMM}", BSLMM = "BayesB", BayesR = "BayesR",
          SuSiE = "SuSiE", fastlmm = "LMM scan")
ROSTER <- c("MS_L_eBIC", "BSLMM", "BayesR", "SuSiE", "fastlmm")

lines <- c("\\begin{tabular}{l rrrrr}", "\\toprule",
  "Method & $\\widehat K$ & TP & Precision & Recall & $F_1$ \\\\")
for (rho in c(0.95, 0.98)) for (sg in c(0, 0.5)) {
  lines <- c(lines, "\\midrule",
    sprintf("\\multicolumn{6}{c}{$\\rho=%.2f$, $\\sigma_g^2=%s$}\\\\", rho, sg),
    "\\midrule")
  for (mm in ROSTER) {
    s <- IN45[IN45$method == mm & abs(IN45$rho - rho) < 1e-9 &
              abs(IN45$sigma_g2 - sg) < 1e-9, ]
    lines <- c(lines,
      sprintf("%s & $%.1f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.3f$ \\\\",
              DISP[mm], s$K_hat, s$tp, s$precision, s$recall, s$f1))
  }
}
writeLines(c(lines, "\\bottomrule", "\\end{tabular}"),
           file.path(OUT, "tab_correlated_causals.tex"))

l2 <- c("\\begin{tabular}{l c rrrr}", "\\toprule",
  "Method & PIP threshold & $\\widehat K$ & Precision & Recall & $F_1$ \\\\")
for (sg in c(0, 0.5)) {
  l2 <- c(l2, "\\midrule",
    sprintf("\\multicolumn{6}{c}{$\\sigma_g^2=%s$}\\\\", sg), "\\midrule")
  for (mm in c("SuSiE", "BayesR")) for (th in c(0.5, 0.9, 0.99)) {
    s <- IN43[IN43$method == mm & abs(IN43$threshold - th) < 1e-9 &
              abs(IN43$sigma_g2 - sg) < 1e-9, ]
    l2 <- c(l2,
      sprintf("%s & $%.2f$ & $%.1f$ & $%.2f$ & $%.2f$ & $%.3f$ \\\\",
              DISP[mm], th, s$K_hat, s$precision, s$recall, s$f1))
  }
}
writeLines(c(l2, "\\bottomrule", "\\end{tabular}"),
           file.path(OUT, "tab_pip_thresholds.tex"))
cat("tables ok\n")
