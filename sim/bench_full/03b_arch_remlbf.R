# ==============================================================================
#  03b_arch_remlbf.R
#  Plug-in REML (REML-BF) mirror of 03_arch.R
#  signal architecture in {weak, med, strong} x sigma_g2 in {0, 0.5}
#  -> 6 cells x B reps. Each architecture preserves the 1+2+2 shape.
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

OUT_DIR <- "results/bench_full/03b_arch_remlbf"
# Signal architectures (matching 03_arch.R)
betas <- list(
  weak   = c(0.4, 0.2, 0.2, 0.1, 0.1),
  med    = c(0.8, 0.4, 0.4, 0.2, 0.2),     # = ANCHOR
  strong = c(1.6, 0.8, 0.8, 0.4, 0.4)
)
cells <- list()
for (sig_name in names(betas))
  for (sg in c(0.0, 0.5))
    cells[[length(cells)+1L]] <- list(
      tag = sprintf("sig%s_sg%.1f", sig_name, sg),
      n = ANCHOR$n, m = ANCHOR$m, rho = ANCHOR$rho,
      K_true = ANCHOR$K_true, beta_true = betas[[sig_name]],
      sigma_g2 = sg, block_size = ANCHOR$block_size,
      seed_offset = match(sig_name, names(betas)) * 1000L)

cat(sprintf("=== 03b_arch_remlbf: signal architecture axis under REML-BF + eBIC ===\n"))
run_axis_remlbf(cells, B = B, out_dir = OUT_DIR, n_cores = N_CORES,
                K_max = ANCHOR$K_max, tau2 = ANCHOR$tau2,
                anchor_seed = 20260425L)
cat("Done.\n")
