# ==============================================================================
#  19_threshold_full_grid.R
#  Extends the threshold-independent comparison of 10_threshold_independent.R
#  from the original 6-cell focused subset to the full 18-cell axis grid.
#  We reuse the per-SNP score extractors from 10_threshold_independent.R
#  and add the 12 missing cells (n=500 x 2, rho=0.80 x 2, weak signal x 2,
#  strong signal x 2, m=10000 x 2 = 10 new cells; 6 already done in
#  results/bench_full/10_threshold_indep/).
#
#  Output:
#    results/bench_full/19_threshold_full_grid/<cell>_b<rep>.rds
#    theory/overleaf_compact/tables/tab_threshold_indep_full.tex
#
#  Run from project root:
#    Rscript sim/bench_full/19_threshold_full_grid.R --cores 8 --B 100
# ==============================================================================
suppressPackageStartupMessages({
  library(parallel)
})

# --- BLAS thread pinning (CRITICAL for mclapply efficiency) ------------------
# Without this, each mclapply worker spawns its own BLAS thread pool, causing
# severe contention on multi-core runs (~10x slowdown on 6 cores observed).
Sys.setenv(OMP_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1",
           BLIS_NUM_THREADS = "1")
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
  cat("BLAS thread pinning applied via RhpcBLASctl.\n")
} else {
  message("RhpcBLASctl not installed -- relying on env vars only.")
}

# Source the existing per-SNP score extractors and aggregator from
# 10_threshold_independent.R (defines topK_*, score_pip_*, compute_metrics).
# We invoke it with sys.source so its main {if (sys.nframe()==0L)} block
# does not auto-run.
local({
  e <- new.env(parent = .GlobalEnv)
  sys.source("sim/bench_full/10_threshold_independent.R", envir = e,
             toplevel.env = e)
  # Export functions of interest to global
  for (nm in c("topK_MSL", "topK_JSL", "topK_SuSiE",
               "topK_BSLMM", "topK_fastlmm",
               "score_pip_SuSiE", "score_pip_BSLMM",
               "score_fastlmm",
               "compute_topK_metrics", "compute_PR_AUC", "compute_brier",
               "run_cell_threshold_indep")) {
    if (exists(nm, envir = e, inherits = FALSE))
      assign(nm, get(nm, envir = e), envir = .GlobalEnv)
  }
})
source("sim/bench_full/00_config.R")

args <- commandArgs(trailingOnly = TRUE)
.arg <- function(flag, def) {
  i <- which(args == flag); if (length(i)) args[i + 1] else def
}
N_CORES <- as.integer(.arg("--cores", 4L))
B       <- as.integer(.arg("--B",     100L))

OUT_DIR <- "results/bench_full/19_threshold_full_grid"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --- Build the 10 NEW cells (the 6 already-done are reused as-is) ----------
build_extra_cells <- function() {
  cells <- list()
  # Axis 1 (n): n=500 only (anchor=n=1000 and n=3000 already in 10_threshold_indep/)
  for (sg in c(0, 0.5))
    cells[[length(cells)+1L]] <- list(
      tag = sprintf("n500_sg%.1f", sg),
      n = 500L, m = 5000L, rho = 0.95,
      K_true = 5L, beta_true = c(0.8, 0.4, 0.4, 0.2, 0.2),
      sigma_g2 = sg, block_size = 10L,
      seed_offset = 10L)
  # Axis 2 (rho): rho=0.80 only (anchor=0.95 and rho=0.98 already done)
  for (sg in c(0, 0.5))
    cells[[length(cells)+1L]] <- list(
      tag = sprintf("rho080_sg%.1f", sg),
      n = 1000L, m = 5000L, rho = 0.80,
      K_true = 5L, beta_true = c(0.8, 0.4, 0.4, 0.2, 0.2),
      sigma_g2 = sg, block_size = 10L,
      seed_offset = 20L)
  # Axis 3 (signal): weak (0.4) and strong (1.6), anchor=med (0.8) already done
  for (sg in c(0, 0.5)) {
    cells[[length(cells)+1L]] <- list(
      tag = sprintf("sigweak_sg%.1f", sg),
      n = 1000L, m = 5000L, rho = 0.95,
      K_true = 5L, beta_true = c(0.4, 0.2, 0.2, 0.1, 0.1),
      sigma_g2 = sg, block_size = 10L,
      seed_offset = 30L)
    cells[[length(cells)+1L]] <- list(
      tag = sprintf("sigstrong_sg%.1f", sg),
      n = 1000L, m = 5000L, rho = 0.95,
      K_true = 5L, beta_true = c(1.6, 0.8, 0.8, 0.4, 0.4),
      sigma_g2 = sg, block_size = 10L,
      seed_offset = 40L)
  }
  # Axis 4 (m): m=10000 only (anchor=m=5000 already done)
  for (sg in c(0, 0.5))
    cells[[length(cells)+1L]] <- list(
      tag = sprintf("m10000_sg%.1f", sg),
      n = 1000L, m = 10000L, rho = 0.95,
      K_true = 5L, beta_true = c(0.8, 0.4, 0.4, 0.2, 0.2),
      sigma_g2 = sg, block_size = 10L,
      seed_offset = 50L)
  cells
}

cells <- build_extra_cells()
cat(sprintf("=== 19_threshold_full_grid: %d new cells x B=%d ===\n",
            length(cells), B))

if (!exists("run_cell_threshold_indep", envir = .GlobalEnv, inherits = FALSE)) {
  stop("run_cell_threshold_indep() not exported from 10_threshold_independent.R\n",
       "  Check that the local() sys.source() picked up the function list.")
}

# --- Main loop: per cell, per replicate ------------------------------------
jobs <- list()
for (i in seq_along(cells)) {
  for (b in seq_len(B)) jobs[[length(jobs)+1L]] <- list(cell = cells[[i]], rep = b)
}
cat(sprintf("Total jobs: %d (cells=%d x B=%d, cores=%d)\n",
            length(jobs), length(cells), B, N_CORES))

t0 <- Sys.time()
mclapply(jobs, function(job) {
  out_fn <- file.path(OUT_DIR,
                      sprintf("%s_b%03d.rds", job$cell$tag, job$rep))
  if (file.exists(out_fn)) return(invisible(NULL))
  cell <- job$cell
  seed <- 20260425L + 1000L * (job$rep - 1L) + as.integer(cell$n) +
          round(1000 * cell$sigma_g2) + cell$seed_offset
  # Redirect output to OUR dir, not the original 10_threshold_indep dir,
  # by setting OUT_DIR in the local env (the function reads it from there).
  old_OUT_DIR <- get("OUT_DIR", envir = .GlobalEnv)
  assign("OUT_DIR", OUT_DIR, envir = .GlobalEnv)
  tryCatch(run_cell_threshold_indep(
             n = cell$n, m = cell$m, rho = cell$rho,
             K_true = cell$K_true, beta_true = cell$beta_true,
             sigma_g2 = cell$sigma_g2, block_size = cell$block_size,
             seed = seed, cell_tag = cell$tag, rep_id = job$rep),
           error = function(e) {
             message(sprintf("Err %s b%d: %s", cell$tag,
                             job$rep, conditionMessage(e)))
           })
  assign("OUT_DIR", old_OUT_DIR, envir = .GlobalEnv)
  NULL
}, mc.cores = N_CORES)
cat(sprintf("Completed in %.1f min\n",
            as.numeric(Sys.time() - t0, units = "mins")))

# --- Aggregate the 6 existing + 10 new cells into one full-grid table ------
collect <- function(dirs) {
  files <- unlist(lapply(dirs, function(d)
    list.files(d, "[.]rds$", full.names = TRUE)))
  files <- grep(" 2[.]rds$", files, invert = TRUE, value = TRUE)
  rows <- mclapply(files, function(f) {
    r <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(r) || is.null(r$metrics)) return(NULL)
    # Flatten metrics into a long data frame
    do.call(rbind, lapply(names(r$metrics), function(mth) {
      data.frame(cell = r$cell_tag, rep = r$rep, n = r$n, m = r$m,
                 rho = r$rho, sigma_g2 = r$sigma_g2, method = mth,
                 topK = r$metrics[[mth]]$top_K_recall %||% r$metrics[[mth]]$topK %||% NA,
                 pr_auc = r$metrics[[mth]]$pr_auc %||% NA,
                 brier  = r$metrics[[mth]]$brier  %||% NA,
                 stringsAsFactors = FALSE)
    }))
  }, mc.cores = N_CORES)
  rows <- rows[!sapply(rows, is.null)]
  do.call(rbind, rows)
}
`%||%` <- function(x, y) if (is.null(x)) y else x

M <- collect(c("results/bench_full/10_threshold_indep", OUT_DIR))
cat(sprintf("Pool of %d (cell, rep, method) rows from both directories\n",
            nrow(M)))
write.csv(M, file.path(OUT_DIR, "threshold_full_grid_long.csv"),
          row.names = FALSE)

# Aggregate per (sigma_g2, method) and per (cell, method) for the table
boot_one <- function(x, R = 2000L) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(c(NA, NA, NA, NA))
  bs <- replicate(R, mean(x[sample.int(length(x), replace = TRUE)]))
  c(mean(x), sd(bs), quantile(bs, 0.025, names = FALSE),
    quantile(bs, 0.975, names = FALSE))
}
# Fairness fix: Top-K* recall is available for every method on all 16 cells.
# PR-AUC and Brier require the proposed framework's per-SNP iterative q^cond,
# which is only exposed on the 6-cell diagnostic subset (anchor, n3000, rho098).
# To avoid an apples-to-oranges pooled comparison (proposed PR-AUC over 6 cells
# vs competitors over 16), we restrict PR-AUC and Brier to the cells where the
# proposed framework has a non-NA value, for ALL methods.
ref_method <- "JS_L_eBIC"
matched_cells <- unique(M$cell[M$method == ref_method & !is.na(M$pr_auc)])
cat(sprintf("Matched cells for PR-AUC/Brier (proposed has values): %s\n",
            paste(matched_cells, collapse = ", ")))

agg <- list()
for (sg2 in c(0, 0.5)) {
  for (mth in unique(M$method)) {
    for (met in c("topK", "pr_auc", "brier")) {
      if (met == "topK") {
        sub <- M[M$sigma_g2 == sg2 & M$method == mth, , drop = FALSE]
      } else {
        # restrict to matched cells for fair comparison
        sub <- M[M$sigma_g2 == sg2 & M$method == mth &
                 M$cell %in% matched_cells, , drop = FALSE]
      }
      b <- boot_one(sub[[met]])
      agg[[length(agg)+1L]] <- data.frame(
        sigma_g2 = sg2, method = mth, metric = met,
        n_rep = nrow(sub),
        mean = b[1], se = b[2], lo = b[3], hi = b[4])
    }
  }
}
agg <- do.call(rbind, agg)
write.csv(agg, file.path(OUT_DIR, "threshold_full_grid_aggregate.csv"),
          row.names = FALSE)

# --- Build LaTeX table -------------------------------------------------------
tex_path <- "theory/overleaf_compact/tables/tab_threshold_indep_full.tex"
display <- c("MS_L_eBIC"   = "\\textsc{MS\\_L\\_eBIC}",
             "JS_L_eBIC"   = "\\textsc{JS\\_L\\_eBIC}",
             "BSLMM"       = "BSLMM",
             "BayesR"      = "BayesR",
             "SuSiE"       = "SuSiE",
             "fastlmm"     = "FaST-LMM")
order <- c("JS_L_eBIC", "MS_L_eBIC", "BSLMM", "BayesR", "SuSiE", "fastlmm")
fmt <- function(m, se, d = 3) {
  if (is.na(m)) "---" else sprintf("$%.*f \\pm %.*f$", d, m, d, se)
}

sink(tex_path)
cat("% Auto-generated by sim/bench_full/19_threshold_full_grid.R\n")
cat("% Top-K* recall only (selected-set metric, all 16 cells, all methods).\n")
cat("\\begin{tabular}{lcc}\n\\toprule\n")
cat("Method & \\multicolumn{2}{c}{Top-$K^\\star$ recall} \\\\\n")
cat("\\cmidrule(lr){2-3}\n")
cat(" & $\\sigma_g^2{=}0$ & $\\sigma_g^2{=}0.5$ \\\\\n")
cat("\\midrule\n")
for (mth in order) {
  row_vals <- character(2)
  for (j in seq_along(c(0, 0.5))) {
    sg2 <- c(0, 0.5)[j]
    r <- agg[agg$method == mth & agg$sigma_g2 == sg2 & agg$metric == "topK", ,
            drop = FALSE]
    row_vals[j] <- if (nrow(r)) fmt(r$mean, r$se) else "---"
  }
  cat(sprintf("%s & %s & %s \\\\\n",
              display[mth] %||% mth, row_vals[1], row_vals[2]))
}
cat("\\bottomrule\n\\end{tabular}\n")
sink()
cat(sprintf("Wrote %s\n", tex_path))
cat("Done.\n")
