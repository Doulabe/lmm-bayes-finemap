# ==============================================================================
#  01b_scaling_n_remlbf.R
#  Plug-in REML (REML-BF) mirror of 01_scaling_n.R
#  n in {500, 1000, 3000} x sigma_g2 in {0, 0.5} -> 6 cells x B reps
#
#  Run: Rscript sim/bench_full/01b_scaling_n_remlbf.R --cores 6 --B 100
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

OUT_DIR <- "results/bench_full/01b_scaling_n_remlbf"
cells <- list()
# NOTE: n=3000 cell skipped in REML-BF mirror -- eigen(K) at n=3000 is
# O(n^3) ~27x more expensive than at n=1000 and would dominate wall time.
# The MBF benchmark in 01_scaling_n/ already covers n=3000; here we only
# need the REML-BF complement at n in {500, 1000}.
for (nval in c(500L, 1000L))
  for (sg in c(0.0, 0.5))
    cells[[length(cells)+1L]] <- list(
      tag = sprintf("n%04d_sg%.1f", nval, sg),
      n = nval, m = ANCHOR$m, rho = ANCHOR$rho,
      K_true = ANCHOR$K_true, beta_true = ANCHOR$beta_true,
      sigma_g2 = sg, block_size = ANCHOR$block_size,
      seed_offset = 0L)

cat(sprintf("=== 01b_scaling_n_remlbf: n axis under REML-BF + eBIC ===\n"))
run_axis_remlbf(cells, B = B, out_dir = OUT_DIR, n_cores = N_CORES,
                K_max = ANCHOR$K_max, tau2 = ANCHOR$tau2,
                anchor_seed = 20260425L)
cat("Done.\n")
