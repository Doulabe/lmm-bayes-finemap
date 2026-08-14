# __________________________________________________________________________
# 01_scaling_n.R
# Scaling in n at fixed (m=5000, rho=0.95, K_true=5, mixed beta).
# n in {500, 1000, 3000} x sigma_g2 in {0, 0.5}  -> 6 cells x B reps.
#
# CLI:
#   Rscript 01_scaling_n.R                    # run all 6 cells x B reps
#   Rscript 01_scaling_n.R --n 1000 --sg 0.5  # one cell, all reps
#   Rscript 01_scaling_n.R --n 1000 --sg 0.5 --rep 3  # single (cell, rep)
#   Rscript 01_scaling_n.R --include-bslmm    # add BSLMM to method set
#   Rscript 01_scaling_n.R --B 5              # override B
# ___________________________________________________________________________

source("sim/bench_full/00_config.R")

OUT_DIR <- "results/bench_full/01_scaling_n"

run_cell <- function(n, sg, rep_id, B, methods) {
  cell_tag <- sprintf("n%04d_sg%.1f", n, sg)
  if (already_done(OUT_DIR, cell_tag, rep_id)) {
    cat(sprintf("  [skip] %s b%02d already done\n", cell_tag, rep_id))
    return(invisible(NULL))
  }
  seed <- 20260425L + 1000L * (rep_id - 1L) + as.integer(n) +
            round(1000 * sg)
  d <- gen_dataset(n = n, m = ANCHOR$m, rho = ANCHOR$rho,
                    K_true = ANCHOR$K_true, beta_true = ANCHOR$beta_true,
                    sigma_g2 = sg, block_size = ANCHOR$block_size, seed = seed)
  res <- run_methods(d, methods, K_max = ANCHOR$K_max,
                     theta = ANCHOR$theta, N_delta = ANCHOR$N_delta,
                     tau2 = ANCHOR$tau2)
  metrics <- compute_metrics(res, d$truth)
  metrics$n <- n; metrics$sigma_g2 <- sg; metrics$rep <- rep_id
  metrics$rho <- ANCHOR$rho; metrics$m <- ANCHOR$m
  payload <- list(metrics = metrics, truth = d$truth, results = res)
  save_cell_rep(OUT_DIR, cell_tag, rep_id, payload)
  cat(sprintf("  [done] %s b%02d: F1 = %s\n", cell_tag, rep_id,
              paste(sprintf("%s=%.2f", metrics$method, metrics$f1),
                    collapse = ", ")))
  invisible(NULL)
}


# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------
if (sys.nframe() == 0L) {
  args <- parse_args()
  B <- if (!is.null(args$B)) as.integer(args$B) else ANCHOR$B
  include_bslmm <- isTRUE(args[["include-bslmm"]])
  methods <- default_methods(include_bslmm)

  n_grid  <- c(500L, 1000L, 3000L)
  sg_grid <- c(0.0, 0.5)
  if (!is.null(args$n))  n_grid  <- as.integer(args$n)
  if (!is.null(args$sg)) sg_grid <- as.numeric(args$sg)
  rep_grid <- if (!is.null(args$rep)) as.integer(args$rep) else seq_len(B)

  cat(sprintf("=== 01_scaling_n: %d cells x %d reps, methods = %s ===\n",
              length(n_grid) * length(sg_grid), length(rep_grid),
              paste(methods, collapse = ",")))
  for (n in n_grid) for (sg in sg_grid) for (rep_id in rep_grid) {
    run_cell(n, sg, rep_id, B = B, methods = methods)
  }
  cat("\nDone.\n")
}
