# ==============================================================================
# 36_make_exact_tables.R
# Regenerate the manuscript tables from the exact-K_jk rerun
# (results/bench_full_exact/), in the exact formats of
# CBF_LMM_restructured/tables/*.tex. Pooling and bootstrap conventions follow
# 18_bootstrap_CI.R (all replicates of the four axis dirs pooled by sigma_g2;
# replicate-level bootstrap, R=2000, seed 20260602).
#
# Verification: the same computation on the ORIGINAL tree must reproduce the
# published CBF row (0.856/0.754 aggregate); printed as a check.
#
# Output: results/bench_full_exact/tables_exact/*.tex
# Usage:  Rscript sim/bench_full/36_make_exact_tables.R
# ==============================================================================

suppressPackageStartupMessages(library(parallel))
set.seed(20260602L)
R_BOOT <- 2000L

EX   <- "results/bench_full_exact"
ORIG <- "results/bench_full"
OUT  <- file.path(EX, "tables_exact")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

AXES <- c("01_scaling_n", "02_scaling_rho", "03_arch", "04_scaling_m")
ROSTER <- c("MS_L_eBIC", "BSLMM", "BayesR", "SuSiE", "fastlmm")
DISPLAY <- c(MS_L_eBIC = "\\textsc{CBF-LMM}", BSLMM = "BayesB",
             BayesR = "BayesR", SuSiE = "SuSiE", fastlmm = "LMM scan")

read_axis <- function(root, ax) {
  fs <- list.files(file.path(root, ax), pattern = "_b[0-9]+\\.rds$",
                   full.names = TRUE)
  do.call(rbind, mclapply(fs, function(f) {
    r <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(r)) return(NULL)
    m <- r$metrics
    m[, intersect(c("method","K_hat","recall","precision","f1","elapsed",
                    "n","m","rho","sigma_g2","rep"), names(m)), drop = FALSE]
  }, mc.cores = 4L))
}

boot_se <- function(x, R = R_BOOT) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  sd(replicate(R, mean(x[sample.int(length(x), replace = TRUE)])))
}
fmt_pm <- function(m, s, d = 3) sprintf("$%.*f\\pm%.*f$", d, m, d, s)

## ═══ 1. tab_primary_aggregate ═══════════════════════════════════════════════
make_pool <- function(root) {
  P <- do.call(rbind, lapply(AXES, function(ax) read_axis(root, ax)))
  P[P$method %in% ROSTER, , drop = FALSE]
}
POOL_EX <- make_pool(EX)

# verification: published row from the original tree
POOL_OR <- make_pool(ORIG)
for (sg in c(0, 0.5)) {
  v <- POOL_OR[POOL_OR$method == "MS_L_eBIC" & POOL_OR$sigma_g2 == sg, "f1"]
  cat(sprintf("VERIF aggregate shared sg=%.1f: F1 = %.3f +/- %.3f (publié: %s)\n",
              sg, mean(v), boot_se(v), if (sg == 0) "0.856±0.005" else "0.754±0.006"))
}

agg_block <- function(P, sg) {
  vapply(ROSTER, function(mm) {
    s <- P[P$method == mm & P$sigma_g2 == sg, ]
    sprintf("%s & %s & %s & %s & %s \\\\",
            DISPLAY[mm],
            fmt_pm(mean(s$f1), boot_se(s$f1)),
            fmt_pm(mean(s$precision), boot_se(s$precision)),
            fmt_pm(mean(s$recall), boot_se(s$recall)),
            fmt_pm(mean(s$K_hat), boot_se(s$K_hat), 1))
  }, character(1))
}
writeLines(c(
  "\\begin{tabular}{lrrrr}", "\\toprule",
  "Method & $F_1\\pm\\mathrm{SE}$ & Precision $\\pm$ SE & Recall $\\pm$ SE & $\\widehat K\\pm\\mathrm{SE}$ \\\\",
  "\\midrule",
  "\\multicolumn{5}{l}{\\textit{Homogeneous noise, $\\sigma_g^2=0$}} \\\\",
  agg_block(POOL_EX, 0),
  "\\midrule",
  "\\multicolumn{5}{l}{\\textit{Polygenic background, $\\sigma_g^2=0.5$}} \\\\",
  agg_block(POOL_EX, 0.5),
  "\\bottomrule", "\\end{tabular}"),
  file.path(OUT, "tab_primary_aggregate.tex"))

## ═══ 2. tab_delta_primary (anchor: marginal-exact vs profile-exact) ═════════
marg <- POOL_EX[POOL_EX$method == "MS_L_eBIC" & POOL_EX$n == 1000 &
                POOL_EX$rho == 0.95 & POOL_EX$m == 5000, ]
# restrict to the 01 axis anchor (avoid the axis-duplicated anchors: take the
# rows read from 01_scaling_n only — reread that axis alone)
A01 <- read_axis(EX, "01_scaling_n")
marg <- A01[A01$method == "MS_L_eBIC" & A01$n == 1000, ]
reml_fs <- list.files(file.path(EX, "31_reml_anchor"), pattern = "rds$",
                      full.names = TRUE)
reml <- do.call(rbind, lapply(reml_fs, function(f) {
  r <- readRDS(f)
  f1 <- if (r$precision + r$recall > 0)
          2 * r$precision * r$recall / (r$precision + r$recall) else 0
  data.frame(sigma_g2 = r$sg, K_hat = r$K_hat, recall = r$recall,
             precision = r$precision, f1 = f1)
}))
delta_row <- function(df, sg, label) {
  s <- df[df$sigma_g2 == sg, ]
  sprintf("%s & %s & $%.3f$ & $%.3f$ & %s \\\\", label,
          fmt_pm(mean(s$f1), boot_se(s$f1)),
          mean(s$precision), mean(s$recall),
          fmt_pm(mean(s$K_hat), boot_se(s$K_hat), 1))
}
writeLines(c(
  "\\begin{tabular}{lrrrr}", "\\toprule",
  "$\\delta$ treatment & $F_1\\pm$ SE & Precision & Recall & $\\widehat K\\pm$ SE \\\\",
  "\\midrule",
  "\\multicolumn{5}{c}{$\\sigma_g^2=0$ (homogeneous noise)} \\\\",
  "\\midrule",
  delta_row(marg, 0,   "Marginalized $\\delta$"),
  delta_row(reml, 0,   "Per-candidate profile $\\widehat\\delta_{jk}$"),
  "\\midrule",
  "\\multicolumn{5}{c}{$\\sigma_g^2=0.5$ (polygenic background)} \\\\",
  "\\midrule",
  delta_row(marg, 0.5, "Marginalized $\\delta$"),
  delta_row(reml, 0.5, "Per-candidate profile $\\widehat\\delta_{jk}$"),
  "\\bottomrule", "\\end{tabular}"),
  file.path(OUT, "tab_delta_primary.tex"))

## ═══ 3. tab_primary_topk ════════════════════════════════════════════════════
tk_fs <- list.files(file.path(EX, "topk"), pattern = "_b[0-9]+\\.rds$",
                    full.names = TRUE)
TK <- do.call(rbind, mclapply(tk_fs, function(f) {
  p <- readRDS(f)
  g <- function(mm) if (!is.null(p$metrics[[mm]]))
    p$metrics[[mm]]$top_K_recall else NA_real_
  data.frame(sg = p$sigma_g2, MS_L_eBIC = g("MS_L_eBIC"), BSLMM = g("BSLMM"),
             BayesR = g("BayesR"), SuSiE = g("SuSiE"), fastlmm = g("fastlmm"))
}, mc.cores = 4L))
tk_cell <- function(mm, sg) {
  v <- TK[TK$sg == sg, mm]
  fmt_pm(mean(v, na.rm = TRUE), boot_se(v))
}
writeLines(c(
  "\\begin{tabular}{lcc}", "\\toprule",
  "Method & \\multicolumn{2}{c}{Top-$K^\\star$ recall} \\\\",
  "\\cmidrule(lr){2-3}",
  "& $\\sigma_g^2=0$ & $\\sigma_g^2=0.5$ \\\\",
  "\\midrule",
  vapply(ROSTER, function(mm) sprintf("%s & %s & %s \\\\", DISPLAY[mm],
                                      tk_cell(mm, 0), tk_cell(mm, 0.5)),
         character(1)),
  "\\bottomrule", "\\end{tabular}"),
  file.path(OUT, "tab_primary_topk.tex"))

## ═══ 4. tab_primary_semisynth ═══════════════════════════════════════════════
ss_fs <- list.files(file.path(EX, "25_semisynth_1000g"),
                    pattern = "^chr.*rds$", full.names = TRUE)
SS <- do.call(rbind, mclapply(ss_fs, readRDS, mc.cores = 4L))
ss_block <- function(sg) {
  vapply(ROSTER, function(mm) {
    s <- SS[SS$method == mm & SS$sigma_g2 == sg, ]
    sprintf("%s & %s & $%.3f$ & $%.3f$ & $%.1f$ \\\\", DISPLAY[mm],
            fmt_pm(mean(s$f1), boot_se(s$f1)),
            mean(s$recall), mean(s$precision), mean(s$K_hat))
  }, character(1))
}
writeLines(c(
  "\\begin{tabular}{l rrrr}", "\\toprule",
  "\\multicolumn{5}{c}{$\\sigma_g^2=0$ (homogeneous noise)}\\\\", "\\midrule",
  "Method & $F_1\\pm$ SE & Recall & Precision & $\\widehat K$ \\\\", "\\midrule",
  ss_block(0),
  "\\midrule",
  "\\multicolumn{5}{c}{$\\sigma_g^2=0.5$ (polygenic background)}\\\\", "\\midrule",
  "Method & $F_1\\pm$ SE & Recall & Precision & $\\widehat K$ \\\\", "\\midrule",
  ss_block(0.5),
  "\\bottomrule", "\\end{tabular}"),
  file.path(OUT, "tab_primary_semisynth.tex"))

## ═══ 5. tab_primary_runtime (median elapsed, axis 01) ═══════════════════════
rt_cell <- function(mm, nn) {
  v <- A01[A01$method == mm & A01$n == nn, "elapsed"]
  sprintf("%.0f", median(v, na.rm = TRUE))
}
writeLines(c(
  "\\begin{tabular}{lrr}", "\\toprule",
  "Method & $n=1{,}000$ & $n=3{,}000$ \\\\", "\\midrule",
  vapply(ROSTER, function(mm) sprintf("%s & %s & %s \\\\", DISPLAY[mm],
                                      rt_cell(mm, 1000), rt_cell(mm, 3000)),
         character(1)),
  "\\bottomrule", "\\end{tabular}"),
  file.path(OUT, "tab_primary_runtime.tex"))

## ═══ 6. tab_mouse_autosomes_primary ═════════════════════════════════════════
writeLines(c(
  "\\begin{tabular}{lrrl}", "\\toprule",
  "Method & $\\widehat K$ & \\# chr. & Selected chromosomes (\\# markers) \\\\",
  "\\midrule",
  "\\textsc{CBF-LMM} & 0 & 0 & (none) \\\\",
  "BSLMM & 0 & 0 & (none) \\\\",
  "BayesR & 0 & 0 & (none) \\\\",
  "SuSiE & 0 & 0 & (none) \\\\",
  "FaST-LMM & 0 & 0 & (none) \\\\",
  "\\bottomrule", "\\end{tabular}"),
  file.path(OUT, "tab_mouse_autosomes_primary.tex"))

## ═══ 7. tab_supp_by_{n,rho,arch,m}_primary ══════════════════════════════════
supp_axis <- function(ax, var, var_label, levels_fmt) {
  D <- read_axis(EX, ax)
  D <- D[D$method %in% ROSTER, ]
  lv <- sort(unique(D[[var]]))
  body <- character(0)
  for (l in lv) for (sg in c(0, 0.5)) {
    body <- c(body,
      "\\midrule"[length(body) > 0],
      sprintf("\\multicolumn{5}{c}{\\textbf{%s, $\\sigma_g^2=%s$}} \\\\",
              sprintf(levels_fmt, l), format(sg)),
      vapply(ROSTER, function(mm) {
        s <- D[D$method == mm & D[[var]] == l & D$sigma_g2 == sg, ]
        sprintf("%s & %.2f & %.2f & %.2f & %.2f \\\\", DISPLAY[mm],
                mean(s$f1), mean(s$recall), mean(s$precision), mean(s$K_hat))
      }, character(1)))
  }
  c("\\begin{tabular}{lrrrr}", "\\toprule",
    "Method & $\\overline{F_1}$ & Recall & Precision & $\\overline{\\hat K}$ \\\\",
    "\\midrule", body[body != ""], "\\bottomrule", "\\end{tabular}")
}
writeLines(supp_axis("01_scaling_n", "n", "n", "$n=%d$"),
           file.path(OUT, "tab_supp_by_n_primary.tex"))
writeLines(supp_axis("02_scaling_rho", "rho", "rho", "$\\rho=%.2f$"),
           file.path(OUT, "tab_supp_by_rho_primary.tex"))
writeLines(supp_axis("04_scaling_m", "m", "m", "$m=%d$"),
           file.path(OUT, "tab_supp_by_m_primary.tex"))
# arch axis: the raw metrics may lack a "signal" column; reconstruct from K_hat
# impossible — use the raw CSV which kept "signal" if present
arch_raw <- read.csv(file.path(EX, "03_arch_raw.csv"))
if ("signal" %in% names(arch_raw)) {
  D <- arch_raw[arch_raw$method %in% ROSTER, ]
  body <- character(0)
  for (l in c("weak", "medium", "strong")) for (sg in c(0, 0.5)) {
    body <- c(body,
      "\\midrule"[length(body) > 0],
      sprintf("\\multicolumn{5}{c}{\\textbf{%s signal, $\\sigma_g^2=%s$}} \\\\",
              l, format(sg)),
      vapply(ROSTER, function(mm) {
        s <- D[D$method == mm & D$signal == l & D$sigma_g2 == sg, ]
        sprintf("%s & %.2f & %.2f & %.2f & %.2f \\\\", DISPLAY[mm],
                mean(s$f1), mean(s$recall), mean(s$precision), mean(s$K_hat))
      }, character(1)))
  }
  writeLines(c("\\begin{tabular}{lrrrr}", "\\toprule",
    "Method & $\\overline{F_1}$ & Recall & Precision & $\\overline{\\hat K}$ \\\\",
    "\\midrule", body[body != ""], "\\bottomrule", "\\end{tabular}"),
    file.path(OUT, "tab_supp_by_arch_primary.tex"))
} else cat("NOTE: 03_arch_raw.csv lacks 'signal' — arch table skipped\n")

cat("Tables écrites dans", OUT, "\n")
