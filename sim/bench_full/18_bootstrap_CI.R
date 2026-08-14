# ==============================================================================
#  18_bootstrap_CI.R
#  Paired-bootstrap confidence intervals on benchmark metrics from the
#  existing per-(cell, replicate) RDS files.  No new compute required.
#
#  Outputs (under results/bench_full/):
#    bootstrap_CI_per_cell.csv   per (cell, method) with mean + SE + CI95
#    bootstrap_CI_aggregate.csv  per (sigma_g2, method) pooled over cells
#
#  And under theory/overleaf_compact/tables/:
#    tab_aggregate_main_CI.tex   updated main aggregate with CIs
#
#  Run from project root:
#    Rscript sim/bench_full/18_bootstrap_CI.R [--R 2000]
# ==============================================================================

suppressPackageStartupMessages({
  library(parallel)
})
args <- commandArgs(trailingOnly = TRUE)
.arg <- function(flag, def) {
  i <- which(args == flag); if (length(i)) args[i + 1] else def
}
R_BOOT  <- as.integer(.arg("--R", 2000L))
N_CORES <- as.integer(.arg("--cores", 4L))
set.seed(20260602L)

# --- Read all per-(cell, rep) RDS from the 4 MBF axis directories AND, if
#     present, the 4 REML-BF mirror directories produced by 0Xb_*_remlbf.R.
axis_dirs <- c("01_scaling_n", "02_scaling_rho", "03_arch", "04_scaling_m")
remlbf_dirs <- c("01b_scaling_n_remlbf", "02b_scaling_rho_remlbf",
                 "03b_arch_remlbf", "04b_scaling_m_remlbf")
# Drop REML-BF dirs that do not yet exist (they will be empty until 0Xb runs)
remlbf_dirs <- Filter(function(d) {
  dir.exists(file.path("results/bench_full", d)) &&
    length(list.files(file.path("results/bench_full", d),
                      pattern = "[.]rds$")) > 0
}, remlbf_dirs)
axis_dirs <- c(axis_dirs, remlbf_dirs)
cat(sprintf("Reading axis dirs: %s\n", paste(axis_dirs, collapse = ", ")))
# Common column set; extra columns (e.g. 03_arch's "signal") are dropped
.COMMON_COLS <- c("method", "K_hat", "tp", "fp", "recall", "precision",
                  "f1", "elapsed", "n", "m", "rho", "sigma_g2", "rep")
read_one <- function(f, axis_tag) {
  r <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(r) || is.null(r$metrics)) return(NULL)
  cell_tag <- sub("_b[0-9]+[.]rds$", "", basename(f))
  out <- r$metrics
  out <- out[, intersect(.COMMON_COLS, names(out)), drop = FALSE]
  out$axis <- axis_tag
  out$cell <- cell_tag
  out
}
metrics_all <- list()
for (d in axis_dirs) {
  path <- file.path("results/bench_full", d)
  files <- list.files(path, pattern = "[.]rds$", full.names = TRUE)
  files <- grep(" 2[.]rds$", files, invert = TRUE, value = TRUE)
  cat(sprintf("Reading %s: %d files (cores=%d)\n", d, length(files), N_CORES))
  flush.console()
  t0 <- Sys.time()
  rows <- mclapply(files, read_one, axis_tag = d, mc.cores = N_CORES)
  rows <- rows[!sapply(rows, is.null)]
  metrics_all[[d]] <- do.call(rbind, rows)
  cat(sprintf("  -> %d rows in %.1f s\n", nrow(metrics_all[[d]]),
              as.numeric(Sys.time() - t0, units = "secs")))
  flush.console()
}
M <- do.call(rbind, metrics_all)
cat(sprintf("Total rows: %d\n", nrow(M)))

# Drop the broken REML_stepwise rows for the headline table (its K_hat = 0
# everywhere indicates a runner bug; will be re-run separately for the
# ablation table in Block 3).
M <- M[M$method != "REML_stepwise", , drop = FALSE]

# Compute FDR row-wise (= 1 - precision when K_hat > 0)
M$fdr <- ifelse(M$K_hat > 0, 1 - M$precision, NA_real_)

# --- Paired bootstrap per (cell, method) --------------------------------------
metrics_to_boot <- c("f1", "recall", "precision", "fdr", "K_hat", "elapsed")
boot_one <- function(x, R = R_BOOT) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(c(mean = mean(x), se = NA, lo = NA, hi = NA))
  bs <- replicate(R, mean(x[sample.int(length(x), replace = TRUE)]))
  c(mean = mean(x),
    se   = sd(bs),
    lo   = unname(quantile(bs, 0.025)),
    hi   = unname(quantile(bs, 0.975)))
}

per_cell_rows <- list()
for (cell_tag in unique(M$cell)) {
  M_c <- M[M$cell == cell_tag, , drop = FALSE]
  for (mth in unique(M_c$method)) {
    sub <- M_c[M_c$method == mth, , drop = FALSE]
    for (met in metrics_to_boot) {
      b <- boot_one(sub[[met]])
      per_cell_rows[[length(per_cell_rows) + 1]] <- data.frame(
        cell = cell_tag, method = mth, metric = met,
        n_rep = nrow(sub),
        n = sub$n[1], m = sub$m[1], rho = sub$rho[1], sigma_g2 = sub$sigma_g2[1],
        mean = unname(b["mean"]),
        se   = unname(b["se"]),
        lo   = unname(b["lo"]),
        hi   = unname(b["hi"]),
        stringsAsFactors = FALSE
      )
    }
  }
}
per_cell <- do.call(rbind, per_cell_rows)
out_dir <- "results/bench_full"
write.csv(per_cell, file.path(out_dir, "bootstrap_CI_per_cell.csv"),
          row.names = FALSE)
cat(sprintf("Wrote %s (%d rows)\n",
            file.path(out_dir, "bootstrap_CI_per_cell.csv"), nrow(per_cell)))

# --- Aggregate by sigma_g2 (pooled over cells, each replicate weighted equally)
# Paired bootstrap at the replicate level pooled across cells within sigma_g2.
agg_rows <- list()
for (sg2 in unique(M$sigma_g2)) {
  M_s <- M[M$sigma_g2 == sg2, , drop = FALSE]
  for (mth in unique(M_s$method)) {
    sub <- M_s[M_s$method == mth, , drop = FALSE]
    for (met in metrics_to_boot) {
      b <- boot_one(sub[[met]])
      agg_rows[[length(agg_rows) + 1]] <- data.frame(
        sigma_g2 = sg2, method = mth, metric = met,
        n_rep = nrow(sub),
        mean = unname(b["mean"]), se = unname(b["se"]),
        lo = unname(b["lo"]), hi = unname(b["hi"]),
        stringsAsFactors = FALSE
      )
    }
  }
}
agg <- do.call(rbind, agg_rows)
write.csv(agg, file.path(out_dir, "bootstrap_CI_aggregate.csv"),
          row.names = FALSE)
cat(sprintf("Wrote %s (%d rows)\n",
            file.path(out_dir, "bootstrap_CI_aggregate.csv"), nrow(agg)))

# --- LaTeX table: main aggregate with 95% CI ----------------------------------
# Methods in display order
method_order <- c("JS_L_eBIC", "MS_L_eBIC", "REML_BF_eBIC", "JS_L_JP99", "MS_L_JP99",
                  "BSLMM", "BayesR", "SuSiE", "fastlmm")
display_name <- c(JS_L_eBIC = "\\textsc{JS\\_L\\_eBIC}",
                  MS_L_eBIC = "\\textsc{MS\\_L\\_eBIC}",
                  REML_BF_eBIC = "\\textsc{JS\\_L\\_REML-BF\\_eBIC}",
                  JS_L_JP99 = "\\textsc{JS\\_L\\_JP99}",
                  MS_L_JP99 = "\\textsc{MS\\_L\\_JP99}",
                  BSLMM     = "BSLMM",
                  BayesR    = "BayesR",
                  SuSiE     = "SuSiE",
                  fastlmm   = "FaST-LMM")

fmt <- function(m, se, dig = 2) {
  if (!is.finite(m)) return("---")
  sprintf("$%.*f \\pm %.*f$", dig, m, dig, se)
}

build_block <- function(sg2_val) {
  out <- character(length(method_order))
  for (i in seq_along(method_order)) {
    mth <- method_order[i]
    row <- function(met) agg[agg$method == mth & agg$sigma_g2 == sg2_val
                              & agg$metric == met, , drop = FALSE]
    f1   <- row("f1");   khat <- row("K_hat")
    rec  <- row("recall"); pre <- row("precision")
    out[i] <- sprintf("%s & %s & %s & %s & %s",
      display_name[mth],
      if (nrow(f1))   fmt(f1$mean,   f1$se, 3)   else "---",
      if (nrow(pre))  fmt(pre$mean,  pre$se, 3)  else "---",
      if (nrow(rec))  fmt(rec$mean,  rec$se, 3)  else "---",
      if (nrow(khat)) fmt(khat$mean, khat$se, 1) else "---")
  }
  out
}

tex_path <- "theory/overleaf_compact/tables/tab_aggregate_main_CI.tex"
dir.create(dirname(tex_path), recursive = TRUE, showWarnings = FALSE)
sink(tex_path)
cat("% Auto-generated by sim/bench_full/18_bootstrap_CI.R\n")
cat("% Paired bootstrap (R = ", R_BOOT, ") on per-replicate metrics,\n", sep = "")
cat("% pooled across the 18 axis cells within each sigma_g^2 regime.\n")
cat("\\begin{tabular}{lcccc}\n\\toprule\n")
cat("Method & $\\overline{F_1} \\pm$ SE & $\\overline{\\mathrm{Prec}} \\pm$ SE & ",
    "$\\overline{\\mathrm{Rec}} \\pm$ SE & $\\overline{\\hat K} \\pm$ SE \\\\\n", sep="")
cat("\\midrule\n\\multicolumn{5}{l}{\\emph{Homogeneous noise, $\\sigma_g^2 = 0$}}\\\\\n")
cat(paste(build_block(0.0), collapse = " \\\\\n"))
cat(" \\\\\n\\midrule\n")
cat("\\multicolumn{5}{l}{\\emph{Polygenic background, $\\sigma_g^2 = 0.5$}}\\\\\n")
cat(paste(build_block(0.5), collapse = " \\\\\n"))
cat(" \\\\\n\\bottomrule\n\\end{tabular}\n")
sink()
cat(sprintf("Wrote %s\n", tex_path))

cat("\nDone.\n")
