# ==============================================================================
# 04_scaling_m.R
# Scaling in m at fixed (n=1000, rho=0.95, K_true=5, mixed beta).
# m in {5000, 10000} x sigma_g2 in {0, 0.5}  -> 4 cells x B reps.
#
# CLI:
#   Rscript 04_scaling_m.R --m 10000 --sg 0.5 --rep 3
# ==============================================================================

source("sim/bench_full/00_config.R")

OUT_DIR <- "results/bench_full/04_scaling_m"

run_cell <- function(m, sg, rep_id, B, methods) {
  cell_tag <- sprintf("m%05d_sg%.1f", m, sg)
  if (already_done(OUT_DIR, cell_tag, rep_id)) {
    cat(sprintf("  [skip] %s b%02d already done\n", cell_tag, rep_id))
    return(invisible(NULL))
  }
  seed <- 20260425L + 1000L * (rep_id - 1L) +
            as.integer(m) + as.integer(round(100 * sg))
  d <- gen_dataset(n = ANCHOR$n, m = m, rho = ANCHOR$rho,
                    K_true = ANCHOR$K_true, beta_true = ANCHOR$beta_true,
                    sigma_g2 = sg, block_size = ANCHOR$block_size, seed = seed)
  res <- run_methods(d, methods, K_max = ANCHOR$K_max,
                     theta = ANCHOR$theta, N_delta = ANCHOR$N_delta,
                     tau2 = ANCHOR$tau2)
  metrics <- compute_metrics(res, d$truth)
  metrics$n <- ANCHOR$n; metrics$m <- m; metrics$rho <- ANCHOR$rho
  metrics$sigma_g2 <- sg; metrics$rep <- rep_id
  payload <- list(metrics = metrics, truth = d$truth, results = res)
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

  m_grid  <- c(5000L, 10000L)
  sg_grid <- c(0.0, 0.5)
  if (!is.null(args$m))  m_grid  <- as.integer(args$m)
  if (!is.null(args$sg)) sg_grid <- as.numeric(args$sg)
  rep_grid <- if (!is.null(args$rep)) as.integer(args$rep) else seq_len(B)

  cat(sprintf("=== 04_scaling_m: %d cells x %d reps ===\n",
              length(m_grid) * length(sg_grid), length(rep_grid)))
  for (m in m_grid) for (sg in sg_grid) for (rep_id in rep_grid) {
    run_cell(m, sg, rep_id, B = B, methods = methods)
  }
  cat("\nDone.\n")
}
