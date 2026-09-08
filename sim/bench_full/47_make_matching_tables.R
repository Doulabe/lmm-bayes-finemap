# ==============================================================================
# 47_make_matching_tables.R
# Build the manuscript semisynthetic tables from the one-to-one matching rerun
# (42_semisynth_matching): pooled Table (tab_primary_semisynth.tex) and the
# by-locus table (tab_semisynth_bylocus_primary.tex). Deterministic bootstrap
# SEs (38_final_pooling convention).
# ==============================================================================
R_BOOT <- 2000L
IN  <- "results/bench_full_exact/42_semisynth_matching"
OUT <- "/Users/kossi/Desktop/Dossier_these/Redaction/LMM_Gaussian_Bayesian/CBF_LMM_restructured/tables"

fs <- list.files(IN, pattern = "^chr.*rds$", full.names = TRUE)
M <- do.call(rbind, lapply(fs, function(f) readRDS(f)$rows))
stopifnot(nrow(M) == 800 * 5)

boot_se <- function(x) {
  x <- x[is.finite(x)]
  set.seed(as.integer(sum(abs(x) * 1e4) %% 2147483647))
  sd(replicate(R_BOOT, mean(x[sample.int(length(x), replace = TRUE)])))
}
ROSTER  <- c("MS_L_eBIC","BSLMM","BayesR","SuSiE","fastlmm")
DISPLAY <- c(MS_L_eBIC="\\textsc{CBF-LMM}", BSLMM="BayesB", BayesR="BayesR",
             SuSiE="SuSiE", fastlmm="LMM scan")

row_line <- function(mm, sg) {
  s <- M[M$method == mm & M$sigma_g2 == sg, ]
  sprintf("%s & $%.3f\\pm%.3f$ & $%.3f$ & $%.3f$ & $%.1f$ \\\\",
          DISPLAY[mm], mean(s$f1), boot_se(s$f1),
          mean(s$recall), mean(s$precision), mean(s$K_hat))
}
block <- function(sg, titre) c(
  sprintf("\\multicolumn{5}{c}{$\\sigma_g^2=%s$ (%s)}\\\\", sg, titre),
  "\\midrule",
  "Method & $F_1\\pm$ SE & Recall & Precision & $\\widehat K$ \\\\",
  "\\midrule",
  vapply(ROSTER, row_line, character(1), sg = as.numeric(sg)))
tab <- c("\\begin{tabular}{l rrrr}", "\\toprule",
         block("0", "homogeneous noise"), "\\midrule",
         block("0.5", "polygenic background"), "\\bottomrule",
         "\\end{tabular}")
writeLines(tab, file.path(OUT, "tab_primary_semisynth.tex"))

cell <- function(mm, locus, sg) {
  s <- M[M$method == mm & M$sigma_g2 == sg & M$locus == locus, ]
  sprintf("$%.3f$", mean(s$f1))
}
byl <- c("\\begin{tabular}{l rr rr}", "\\toprule",
  " & \\multicolumn{2}{c}{chr1 ($m=1{,}493$)} & \\multicolumn{2}{c}{chr6 ($m=2{,}143$)} \\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
  "Method & $\\sigma_g^2{=}0$ & $\\sigma_g^2{=}0.5$ & $\\sigma_g^2{=}0$ & $\\sigma_g^2{=}0.5$ \\\\",
  "\\midrule",
  vapply(ROSTER, function(mm) sprintf("%s & %s & %s & %s & %s \\\\",
    DISPLAY[mm], cell(mm,"chr1",0), cell(mm,"chr1",0.5),
    cell(mm,"chr6",0), cell(mm,"chr6",0.5)), character(1)),
  "\\bottomrule", "\\end{tabular}")
writeLines(byl, file.path(OUT, "tab_semisynth_bylocus_primary.tex"))
cat("tables ecrites\n")
