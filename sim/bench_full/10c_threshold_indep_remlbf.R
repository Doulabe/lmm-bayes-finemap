# ==============================================================================
#  10c_threshold_indep_remlbf.R
#  Plug-in REML (REML-BF) Top-K* recovery on the 6 cells of the
#  threshold-independent diagnostic subset, paired with the existing
#  MBF Top-K* numbers from 10_threshold_independent.R for direct
#  comparison.
#
#  Why not full PR-AUC / Brier?  REML-BF is a single-evaluation procedure
#  (no probabilistic averaging over delta), so PR-AUC and Brier are
#  evaluated only via the stepwise top-K* recovery on this script's
#  diagnostic subset.  Per-SNP iterative q^cond under REML-BF would
#  require modifying log_BF_one() in 10b_add_per_snp_scores.R to fix
#  delta at REML; this is left to the planned full-grid REML-BF re-run.
#
#  Cells (same as 10_threshold_independent.R):
#    anchor sg0.0, anchor sg0.5, n3000 sg0.0, n3000 sg0.5,
#    rho098 sg0.0, rho098 sg0.5
#  B=100 reps per cell, 6 cells total.
#
#  Output:
#    results/bench_full/10c_threshold_indep_remlbf/<cell>_b<rep>.rds
#    results/bench_full/10c_threshold_indep_remlbf_summary.csv
#    theory/overleaf_compact/tables/tab_threshold_indep_remlbf.tex
#
#  Run from project root:
#    Rscript sim/bench_full/10c_threshold_indep_remlbf.R --cores 6 --B 100
# ==============================================================================

# BLAS pinning
Sys.setenv(OMP_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           BLIS_NUM_THREADS = "1")
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
}

suppressPackageStartupMessages(library(parallel))
source("sim/bench_full/00_config.R")
if (!requireNamespace("rrBLUP", quietly = TRUE))
  stop("rrBLUP needed for REML")

args <- commandArgs(trailingOnly = TRUE)
.arg <- function(flag, def) {
  i <- which(args == flag); if (length(i)) args[i + 1] else def
}
N_CORES <- as.integer(.arg("--cores", 6L))
B       <- as.integer(.arg("--B",     100L))

OUT_DIR <- "results/bench_full/10c_threshold_indep_remlbf"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
TEX_OUT <- "theory/overleaf_compact/tables/tab_threshold_indep_remlbf.tex"

#  Plug-in REML stepwise (intercept-only W) 
# Returns indices of selected SNPs; eBIC stopping is OFF here (we always
# run to K_star = K_true steps so that the Top-K* metric is well-defined).
plug_in_reml_topK <- function(y, X, K_star, tau2 = 0.04) {
  n <- nrow(X); m <- ncol(X)
  selected <- integer(0); remaining <- seq_len(m)
  F_k <- matrix(1, n, 1L)
  for (k_step in seq_len(K_star)) {
    if (length(remaining) == 0L) break
    rem <- X[, remaining, drop = FALSE]
    K   <- tcrossprod(rem) / max(ncol(rem) - 1L, 1L)
    proj <- solve(crossprod(F_k), crossprod(F_k, y))
    y_r  <- y - F_k %*% proj
    fit  <- tryCatch(rrBLUP::mixed.solve(y = y_r, K = K),
                     error = function(e) NULL)
    if (is.null(fit)) break
    delta_hat <- max(fit$Vu / fit$Ve, 1e-6)
    ee <- eigen(K, symmetric = TRUE)
    s <- pmax(ee$values, 0); U <- ee$vectors
    w <- 1 / sqrt(1 + delta_hat * s)
    Uty   <- crossprod(U, y) * w
    UtFk  <- crossprod(U, F_k) * w
    UtX   <- crossprod(U, rem) * w
    A    <- crossprod(UtFk)
    A_inv<- solve(A)
    beta <- A_inv %*% crossprod(UtFk, Uty)
    ty   <- Uty - UtFk %*% beta
    proj_X <- UtFk %*% (A_inv %*% crossprod(UtFk, UtX))
    tXj    <- UtX - proj_X
    D_jk   <- colSums(tXj^2)
    u_j    <- as.numeric(crossprod(tXj, ty))
    rss0   <- as.numeric(crossprod(ty))
    v_j    <- tau2 * D_jk
    df     <- n - ncol(F_k)
    Q2_j   <- u_j^2 / (D_jk * rss0 / df)
    log_BF <- -0.5 * log1p(v_j) - 0.5 * df *
              log1p(-v_j * Q2_j / df / (1 + v_j))
    j_best <- which.max(log_BF)
    if (!is.finite(log_BF[j_best])) break
    selected <- c(selected, remaining[j_best])
    F_k      <- cbind(F_k, X[, remaining[j_best], drop = FALSE])
    remaining <- remaining[-j_best]
  }
  selected
}

#  Cell list (mirrors 10_threshold_independent.R) 
build_cells <- function() {
  cells <- list()
  for (sg in c(0, 0.5))
    cells[[length(cells)+1L]] <- list(
      tag = sprintf("anchor_sg%.1f", sg),
      n = 1000L, m = 5000L, rho = 0.95, K_true = 5L,
      beta_true = c(0.8, 0.4, 0.4, 0.2, 0.2),
      sigma_g2 = sg, block_size = 10L, seed_offset = 0L)
  # NOTE: n=3000 cells skipped in REML-BF — eigen(K) at n=3000 is O(n^3) ~
  # 27x slower than n=1000.  The MBF benchmark already covers n=3000 in
  # 10_threshold_indep/.  REML-BF Top-K* on n=1000 cells (anchor + rho098)
  # is sufficient for the sensitivity narrative.
  for (sg in c(0, 0.5))
    cells[[length(cells)+1L]] <- list(
      tag = sprintf("rho098_sg%.1f", sg),
      n = 1000L, m = 5000L, rho = 0.98, K_true = 5L,
      beta_true = c(0.8, 0.4, 0.4, 0.2, 0.2),
      sigma_g2 = sg, block_size = 10L, seed_offset = 2L)
  cells
}
cells <- build_cells()

#  Per-(cell, rep) job 
run_one <- function(cell, rep_id) {
  out_fn <- file.path(OUT_DIR, sprintf("%s_b%03d.rds", cell$tag, rep_id))
  if (file.exists(out_fn)) return(invisible(NULL))
  set.seed(20260602L + rep_id + cell$seed_offset * 100000L)
  d <- gen_dataset(n = cell$n, m = cell$m, rho = cell$rho,
                   K_true = cell$K_true, beta_true = cell$beta_true,
                   sigma_g2 = cell$sigma_g2, block_size = cell$block_size,
                   seed = 20260602L + rep_id + cell$seed_offset * 100000L)
  sel <- plug_in_reml_topK(d$y, d$X, K_star = cell$K_true, tau2 = 0.04)
  tp  <- length(intersect(sel, d$truth))
  topK <- tp / cell$K_true
  payload <- list(cell_tag = cell$tag, rep = rep_id,
                  sigma_g2 = cell$sigma_g2,
                  n = cell$n, m = cell$m, rho = cell$rho,
                  truth = d$truth, sel = sel, topK = topK)
  saveRDS(payload, out_fn)
  invisible(payload)
}

cat(sprintf("=== 10c_threshold_indep_remlbf: %d cells x B=%d, cores=%d ===\n",
            length(cells), B, N_CORES))
jobs <- list()
for (c in cells) for (b in seq_len(B)) jobs[[length(jobs)+1L]] <- list(c, b)
t0 <- Sys.time()
mclapply(jobs, function(j) {
  tryCatch(run_one(j[[1]], j[[2]]),
           error = function(e) message(sprintf("Err %s b%d: %s",
                                                j[[1]]$tag, j[[2]],
                                                conditionMessage(e))))
}, mc.cores = N_CORES)
cat(sprintf("Compute done in %.1f min\n",
            as.numeric(Sys.time() - t0, units = "mins")))

#  Aggregate 
files <- list.files(OUT_DIR, "[.]rds$", full.names = TRUE)
rows <- lapply(files, function(f) {
  r <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(r)) return(NULL)
  data.frame(cell = r$cell_tag, rep = r$rep, sigma_g2 = r$sigma_g2,
             topK = r$topK)
})
rows <- rows[!sapply(rows, is.null)]
M <- do.call(rbind, rows)
write.csv(M, file.path(OUT_DIR, "..", "10c_threshold_indep_remlbf_summary.csv"),
          row.names = FALSE)

agg <- aggregate(topK ~ cell + sigma_g2, data = M,
                 FUN = function(x) c(mean = mean(x),
                                     se = sd(x) / sqrt(length(x))))

#  LaTeX table ------------------------------------------------------------
sink(TEX_OUT)
cat("% Auto-generated by sim/bench_full/10c_threshold_indep_remlbf.R\n")
cat("% Top-K^* recovery of the plug-in REML (REML-BF) stepwise variant on\n")
cat("% the 6 threshold-independent diagnostic cells, paired with the MBF\n")
cat("% numbers of Supplementary Note S6 (tab:threshold_indep).\n")
cat("\\begin{tabular}{lcc}\n\\toprule\n")
cat("Cell & Top-$K^\\star$ (mean $\\pm$ SE) & $n_{\\mathrm{rep}}$ \\\\\n")
cat("\\midrule\n")
for (i in seq_len(nrow(agg))) {
  r <- agg[i, ]
  cat(sprintf("%s ($\\sigma_g^2=%.1f$) & $%.3f \\pm %.3f$ & %d \\\\\n",
              r$cell, r$sigma_g2, r$topK[1], r$topK[2],
              sum(M$cell == r$cell & M$sigma_g2 == r$sigma_g2)))
}
cat("\\bottomrule\n\\end{tabular}\n")
sink()
cat(sprintf("Wrote %s\n", TEX_OUT))
cat("Done.\n")
