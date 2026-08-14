# ==============================================================================
# 03_arch.R
# Signal-strength scaling around the mixed beta architecture (Q1=(c) anchor).
# Three signal levels (weak, medium, strong) all preserve the
# 1-strong + 2-moderate + 2-weak ratio.
#
# Levels:
#   weak:    beta = (0.4, 0.2, 0.2, 0.1, 0.1)
#   medium:  beta = (0.8, 0.4, 0.4, 0.2, 0.2)   (= ANCHOR)
#   strong:  beta = (1.6, 0.8, 0.8, 0.4, 0.4)
#
# 3 signal levels x sigma_g2 in {0, 0.5} -> 6 cells x B reps.
#
# CLI:
#   Rscript 03_arch.R --signal medium --sg 0 --rep 3
# ==============================================================================

source("sim/bench_full/00_config.R")

OUT_DIR <- "results/bench_full/03_arch"

beta_for_signal <- function(level) {
  switch(level,
         weak    = c(0.4, 0.2, 0.2, 0.1, 0.1),
         medium  = c(0.8, 0.4, 0.4, 0.2, 0.2),
         strong  = c(1.6, 0.8, 0.8, 0.4, 0.4),
         stop("Unknown signal level: ", level))
}


run_cell <- function(signal, sg, rep_id, B, methods) {
  cell_tag <- sprintf("sig%s_sg%.1f", signal, sg)
  if (already_done(OUT_DIR, cell_tag, rep_id)) {
    cat(sprintf("  [skip] %s b%02d already done\n", cell_tag, rep_id))
    return(invisible(NULL))
  }
  beta_true <- beta_for_signal(signal)
  seed <- 20260425L + 1000L * (rep_id - 1L) +
            switch(signal, weak = 1L, medium = 2L, strong = 3L) +
            as.integer(round(100 * sg))
  d <- gen_dataset(n = ANCHOR$n, m = ANCHOR$m, rho = ANCHOR$rho,
                    K_true = length(beta_true), beta_true = beta_true,
                    sigma_g2 = sg, block_size = ANCHOR$block_size, seed = seed)
  res <- run_methods(d, methods, K_max = ANCHOR$K_max,
                     theta = ANCHOR$theta, N_delta = ANCHOR$N_delta,
                     tau2 = ANCHOR$tau2)
  metrics <- compute_metrics(res, d$truth)
  metrics$n <- ANCHOR$n; metrics$m <- ANCHOR$m; metrics$rho <- ANCHOR$rho
  metrics$signal <- signal; metrics$sigma_g2 <- sg; metrics$rep <- rep_id
  payload <- list(metrics = metrics, truth = d$truth, results = res,
                    beta_true = beta_true)
  save_cell_rep(OUT_DIR, cell_tag, rep_id, payload)
  cat(sprintf("  [done] %s b%02d: F1 = %s\n", cell_tag, rep_id,
              paste(sprintf("%s=%.2f", metrics$method, metrics$f1),
                    collapse = ", ")))
  invisible(NULL)
}

if (sys.nframe() == 0L) {
  args <- parse_args()
  B <- if (!is.null(args$B)) as.integer(args$B) else ANCHOR$B
  include_bslmm <- isTRUE(args[["include-bslmm"]])
  methods <- default_methods(include_bslmm)

  signal_grid <- c("weak", "medium", "strong")
  sg_grid     <- c(0.0, 0.5)
  if (!is.null(args$signal)) signal_grid <- as.character(args$signal)
  if (!is.null(args$sg))     sg_grid     <- as.numeric(args$sg)
  rep_grid <- if (!is.null(args$rep)) as.integer(args$rep) else seq_len(B)

  cat(sprintf("=== 03_arch: %d cells x %d reps ===\n",
              length(signal_grid) * length(sg_grid), length(rep_grid)))
  for (signal in signal_grid) for (sg in sg_grid) for (rep_id in rep_grid) {
    run_cell(signal, sg, rep_id, B = B, methods = methods)
  }
  cat("\nDone.\n")
}
