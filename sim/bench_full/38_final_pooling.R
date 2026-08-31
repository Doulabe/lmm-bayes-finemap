# ==============================================================================
# 38_final_pooling.R
# THE pooling rule of the paper, applied uniformly to every pooled statistic:
#   the 16 unique design cells, anchor counted once (its sample-size-axis
#   instance), all replicates weighted equally (800 per polygenic regime).
# Cells per regime: n500, n1000(anchor), n3000 (axis 01); rho080, rho098
# (axis 02); sigweak, sigstrong (axis 03); m10000 (axis 04).
# Computes, on this common set:
#   * pooled aggregate for the roster + exact vs shared + varbvs reference
#   * pooled Top-K* recall (same 16 cells) for exact, shared, comparators
# and writes tab_primary_aggregate.tex (paper format), the S3 ablation table,
# and the S6 varbvs reference table.
# Usage: Rscript sim/bench_full/38_final_pooling.R
# ==============================================================================

suppressPackageStartupMessages(library(parallel))
set.seed(20260602L)
R_BOOT <- 2000L
EX  <- "results/bench_full_exact"
OUT <- file.path(EX, "tables_exact")

KEEP <- list("01_scaling_n"  = c("n0500", "n1000", "n3000"),
             "02_scaling_rho" = c("rho0.80", "rho0.98"),
             "03_arch"        = c("sigweak", "sigstrong"),
             "04_scaling_m"   = c("m10000"))

read_cells <- function() {
  do.call(rbind, unlist(recursive = FALSE, lapply(names(KEEP), function(ax) {
    lapply(KEEP[[ax]], function(cell) {
      fs <- list.files(file.path(EX, ax),
                       pattern = paste0("^", gsub("\\.", "[.]", cell),
                                        "_sg[0-9.]+_b[0-9]+[.]rds$"),
                       full.names = TRUE)
      do.call(rbind, mclapply(fs, function(f) {
        m <- readRDS(f)$metrics
        m$cell <- cell
        m[, intersect(c("method","K_hat","recall","precision","f1",
                        "sigma_g2","cell"), names(m))]
      }, mc.cores = 4L))
    })
  })))
}
M <- read_cells()
cat("réplicats poolés par régime:",
    nrow(M[M$method == "MS_L_eBIC" & M$sigma_g2 == 0, ]), "\n")

boot_se <- function(x) {
  x <- x[is.finite(x)]
  # deterministic per data vector: identical statistics get identical SEs
  set.seed(as.integer(sum(abs(x) * 1e4) %% 2147483647))
  sd(replicate(R_BOOT, mean(x[sample.int(length(x), replace = TRUE)])))
}
fmt_pm <- function(m, s, d = 3) sprintf("$%.*f\\pm%.*f$", d, m, d, s)
fmt_K  <- function(m, s) sprintf("$%.1f\\pm%.2f$", m, s)
stat_row <- function(mm, sg, cols = c("f1","precision","recall","K_hat")) {
  s <- M[M$method == mm & M$sigma_g2 == sg, ]
  sapply(cols, function(cn) c(mean(s[[cn]]), boot_se(s[[cn]])))
}

ROSTER <- c("MS_L_eBIC","BSLMM","BayesR","SuSiE","fastlmm")
DISPLAY <- c(MS_L_eBIC="\\textsc{CBF-LMM}", BSLMM="BayesB", BayesR="BayesR",
             SuSiE="SuSiE", fastlmm="LMM scan",
             MS_L_eBIC_sharedGk="Pool-kernel variant", varbvs="varbvs")

## ── 1. tab_primary_aggregate (paper rule) ───────────────────────────────────
blk <- function(sg) vapply(ROSTER, function(mm) {
  r <- stat_row(mm, sg)
  sprintf("%s & %s & %s & %s & %s \\\\", DISPLAY[mm],
          fmt_pm(r[1,"f1"], r[2,"f1"]), fmt_pm(r[1,"precision"], r[2,"precision"]),
          fmt_pm(r[1,"recall"], r[2,"recall"]), fmt_K(r[1,"K_hat"], r[2,"K_hat"]))
}, character(1))
writeLines(c("\\begin{tabular}{lrrrr}", "\\toprule",
  "Method & $F_1\\pm\\mathrm{SE}$ & Precision $\\pm$ SE & Recall $\\pm$ SE & $\\widehat K\\pm\\mathrm{SE}$ \\\\",
  "\\midrule",
  "\\multicolumn{5}{l}{\\textit{Homogeneous noise, $\\sigma_g^2=0$}} \\\\",
  blk(0), "\\midrule",
  "\\multicolumn{5}{l}{\\textit{Polygenic background, $\\sigma_g^2=0.5$}} \\\\",
  blk(0.5), "\\bottomrule", "\\end{tabular}"),
  file.path(OUT, "tab_primary_aggregate.tex"))

## pooled headline numbers for the text
for (mm in c("MS_L_eBIC", "MS_L_eBIC_sharedGk", "BSLMM", "varbvs"))
  for (sg in c(0, 0.5)) {
    r <- stat_row(mm, sg)
    cat(sprintf("%-20s sg=%.1f  F1=%.3f±%.3f  prec=%.3f rec=%.3f K=%.2f\n",
                mm, sg, r[1,"f1"], r[2,"f1"], r[1,"precision"],
                r[1,"recall"], r[1,"K_hat"]))
  }

## ── 2. Top-K on the SAME 16 cells ───────────────────────────────────────────
tk_fs <- list.files(file.path(EX, "topk"), pattern = "_b[0-9]+[.]rds$",
                    full.names = TRUE)
TK <- do.call(rbind, mclapply(tk_fs, function(f) {
  p <- readRDS(f)
  g <- function(mm) if (!is.null(p$metrics[[mm]]))
    p$metrics[[mm]]$top_K_recall else NA_real_
  data.frame(cell = p$cell_tag, sg = p$sigma_g2,
             MS_L_eBIC = g("MS_L_eBIC"),
             MS_L_eBIC_sharedGk = g("MS_L_eBIC_sharedGk"),
             BSLMM = g("BSLMM"), BayesR = g("BayesR"),
             SuSiE = g("SuSiE"), fastlmm = g("fastlmm"))
}, mc.cores = 4L))
# topk cell tags: anchor, n3000, rho098, n500, rho080, sigweak, sigstrong,
# m10000 — exactly the 16 unique cells, anchor once. All kept.
tk_cell <- function(mm, sg) {
  v <- TK[TK$sg == sg, mm]
  c(mean(v, na.rm = TRUE), boot_se(v))
}
writeLines(c("\\begin{tabular}{lcc}", "\\toprule",
  "Method & \\multicolumn{2}{c}{Top-$K^\\star$ recall} \\\\",
  "\\cmidrule(lr){2-3}", "& $\\sigma_g^2=0$ & $\\sigma_g^2=0.5$ \\\\",
  "\\midrule",
  vapply(ROSTER, function(mm) {
    a <- tk_cell(mm, 0); b <- tk_cell(mm, 0.5)
    sprintf("%s & %s & %s \\\\", DISPLAY[mm],
            fmt_pm(a[1], a[2]), fmt_pm(b[1], b[2]))
  }, character(1)),
  "\\bottomrule", "\\end{tabular}"),
  file.path(OUT, "tab_primary_topk.tex"))
for (mm in c("MS_L_eBIC", "MS_L_eBIC_sharedGk"))
  cat(sprintf("topk %-20s sg0=%.3f sg0.5=%.3f\n", mm,
              tk_cell(mm, 0)[1], tk_cell(mm, 0.5)[1]))

## ── 3. S3 ablation table: exact vs pool-kernel variant ──────────────────────
ab_row <- function(mm, sg) {
  r <- stat_row(mm, sg)
  tk <- tk_cell(mm, sg)
  sprintf("%s & %s & $%.3f$ & $%.3f$ & %s & %s \\\\",
          if (mm == "MS_L_eBIC") "Candidate-specific $K_{jk}$ (primary)"
          else "Pool kernel $G_k$ for all candidates",
          fmt_pm(r[1,"f1"], r[2,"f1"]), r[1,"precision"], r[1,"recall"],
          fmt_K(r[1,"K_hat"], r[2,"K_hat"]), fmt_pm(tk[1], tk[2]))
}
writeLines(c("\\begin{tabular}{lrrrrr}", "\\toprule",
  "Kernel & $F_1\\pm$ SE & Precision & Recall & $\\widehat K\\pm$ SE & Top-$K^\\star\\pm$ SE \\\\",
  "\\midrule",
  "\\multicolumn{6}{c}{$\\sigma_g^2=0$ (homogeneous noise)} \\\\", "\\midrule",
  ab_row("MS_L_eBIC", 0), ab_row("MS_L_eBIC_sharedGk", 0),
  "\\midrule",
  "\\multicolumn{6}{c}{$\\sigma_g^2=0.5$ (polygenic background)} \\\\", "\\midrule",
  ab_row("MS_L_eBIC", 0.5), ab_row("MS_L_eBIC_sharedGk", 0.5),
  "\\bottomrule", "\\end{tabular}"),
  file.path(OUT, "tab_kernel_ablation.tex"))

## ── 4. S6 varbvs reference table ────────────────────────────────────────────
vb_row <- function(sg) {
  r <- stat_row("varbvs", sg)
  c <- stat_row("MS_L_eBIC", sg)
  sprintf("$%s$ & %s & $%.3f$ & $%.3f$ & %s & %s \\\\", format(sg),
          fmt_pm(r[1,"f1"], r[2,"f1"]), r[1,"precision"], r[1,"recall"],
          fmt_K(r[1,"K_hat"], r[2,"K_hat"]),
          fmt_pm(c[1,"f1"], c[2,"f1"]))
}
writeLines(c("\\begin{tabular}{lrrrrr}", "\\toprule",
  "$\\sigma_g^2$ & $F_1\\pm$ SE & Precision & Recall & $\\widehat K\\pm$ SE & \\textsc{CBF-LMM} $F_1$ \\\\",
  "\\midrule", vb_row(0), vb_row(0.5),
  "\\bottomrule", "\\end{tabular}"),
  file.path(OUT, "tab_varbvs_reference.tex"))

cat("tables écrites dans", OUT, "\n")
