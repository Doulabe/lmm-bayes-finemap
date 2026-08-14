# ==============================================================================
#  04b_scaling_m_remlbf.R
#  Plug-in REML (REML-BF) mirror of 04_scaling_m.R
#  m in {5000, 10000} x sigma_g2 in {0, 0.5} -> 4 cells x B reps
# ==============================================================================
Sys.setenv(OMP_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           BLIS_NUM_THREADS = "1")
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
  cat("BLAS thread pinning applied.\n")
}

source("sim/bench_full/00_config.R")
source("sim/bench_full/_remlbf_utils.R")
if (!requireNamespace("rrBLUP", quietly = TRUE))
  stop("rrBLUP needed")

args <- commandArgs(trailingOnly = TRUE)
.arg <- function(flag, def) { i <- which(args == flag); if (length(i)) args[i+1] else def }
N_CORES <- as.integer(.arg("--cores", 6L))
B       <- as.integer(.arg("--B",     100L))

OUT_DIR <- "results/bench_full/04b_scaling_m_remlbf"
cells <- list()
for (mval in c(5000L, 10000L))
  for (sg in c(0.0, 0.5))
    cells[[length(cells)+1L]] <- list(
      tag = sprintf("m%05d_sg%.1f", mval, sg),
      n = ANCHOR$n, m = mval, rho = ANCHOR$rho,
      K_true = ANCHOR$K_true, beta_true = ANCHOR$beta_true,
      sigma_g2 = sg, block_size = ANCHOR$block_size,
      seed_offset = mval)

cat(sprintf("=== 04b_scaling_m_remlbf: m axis under REML-BF + eBIC ===\n"))
run_axis_remlbf(cells, B = B, out_dir = OUT_DIR, n_cores = N_CORES,
                K_max = ANCHOR$K_max, tau2 = ANCHOR$tau2,
                anchor_seed = 20260425L)
cat("Done.\n")
