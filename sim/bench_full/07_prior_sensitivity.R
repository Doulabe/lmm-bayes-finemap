# ==============================================================================
# 07_prior_sensitivity.R
# Sensitivity analysis: half-Cauchy vs inverse-gamma prior on delta.
#
# For a representative subset of benchmark cells (anchor + neighbours), regenerate
# the data with the same seed and run MS_L_eBIC + JS_L_eBIC under TWO priors:
#   (a) half_cauchy on sqrt(delta) with scale = 1   (default)
#   (b) inv_gamma(shape=1, scale=1) on delta        (alternative)
# Then compare (i) log BF at each step, (ii) selected sets, (iii) F1.
#
# Output: results/bench_full/07_prior_sensitivity_summary.csv
# ==============================================================================

# Try multiple paths for 00_config.R so the script works from project
# root or from sim/bench_full/.
.find_config <- function() {
  for (p in c("sim/bench_full/00_config.R", "00_config.R",
              "bench_full/00_config.R", "../bench_full/00_config.R")) {
    if (file.exists(p)) { source(p); return(p) }
  }
  stop("00_config.R not found from CWD = ", getwd())
}
BENCH_DIR <- dirname(normalizePath(.find_config()))

OUT_DIR <- "results/bench_full"
OUT_CSV <- file.path(OUT_DIR, "07_prior_sensitivity_summary.csv")


run_sensitivity_cell <- function(n, m, rho, K_true, beta_true,
                                    sigma_g2, block_size, seed,
                                    cell_label) {
  d <- gen_dataset(n = n, m = m, rho = rho, K_true = K_true,
                     beta_true = beta_true, sigma_g2 = sigma_g2,
                     block_size = block_size, seed = seed)
  res_list <- list()
  for (prior in c("half_cauchy", "inv_gamma")) {
    for (method in c("MS_L_eBIC", "JS_L_eBIC")) {
      fn <- if (startsWith(method, "JS_L")) JS_L_LMM_stepwise_fast
            else MS_L_LMM_stepwise_fast
      t0 <- Sys.time()
      res <- fn(d$y, d$X, tau2 = ANCHOR$tau2, K_max = ANCHOR$K_max,
                  criterion = "eBIC",
                  theta = ANCHOR$theta, n_nodes = ANCHOR$N_delta,
                  prior = prior,
                  prior_params = if (prior == "half_cauchy")
                                    list(scale = 1.0)
                                  else list(shape = 1, scale = 1))
      elapsed <- as.numeric(Sys.time() - t0, units = "secs")
      tp <- length(intersect(res$indices, d$truth))
      fp <- res$K_hat - tp
      rec <- tp / K_true
      prec <- if (res$K_hat > 0) tp / res$K_hat else 0
      f1 <- if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0
      res_list[[paste(method, prior, sep="|")]] <- list(
        method = method, prior = prior, indices = res$indices,
        K_hat = res$K_hat, F1 = f1, recall = rec, precision = prec,
        log_BF_path = res$log_BF_at_step, elapsed = elapsed)
    }
  }
  # Compare half_cauchy vs inv_gamma per method
  rows <- list()
  for (method in c("MS_L_eBIC", "JS_L_eBIC")) {
    hc <- res_list[[paste0(method, "|half_cauchy")]]
    ig <- res_list[[paste0(method, "|inv_gamma")]]
    common <- length(intersect(hc$indices, ig$indices))
    union  <- length(union(hc$indices, ig$indices))
    jaccard <- if (union > 0) common / union else 1
    # log-BF path comparison (first min(K_hc, K_ig) steps)
    K_common <- min(hc$K_hat, ig$K_hat)
    max_dBF <- if (K_common > 0)
      max(abs(hc$log_BF_path[seq_len(K_common)] -
                ig$log_BF_path[seq_len(K_common)])) else NA
    rows[[length(rows) + 1L]] <- data.frame(
      cell = cell_label, method = method,
      K_hc = hc$K_hat, K_ig = ig$K_hat,
      F1_hc = hc$F1,   F1_ig = ig$F1,
      delta_F1 = ig$F1 - hc$F1,
      jaccard = jaccard,
      n_common = common, n_union = union,
      max_dlogBF = max_dBF,
      time_hc = hc$elapsed, time_ig = ig$elapsed,
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}


if (sys.nframe() == 0L) {
  saved_wd <- getwd()
  setwd(file.path(BENCH_DIR, "..", ".."))   # project root
  on.exit(setwd(saved_wd), add = TRUE)

  # Representative cells: anchor + neighbours (3 along each of 4 axes × 2 sg)
  cells <- list()
  for (n_val in c(500L, 1000L, 3000L))
    for (sg in c(0, 0.5))
      cells[[length(cells)+1L]] <- list(
        label = sprintf("n=%d,sg=%.1f", n_val, sg),
        n = n_val, m = ANCHOR$m, rho = ANCHOR$rho, K_true = ANCHOR$K_true,
        beta_true = ANCHOR$beta_true, sigma_g2 = sg,
        block_size = ANCHOR$block_size, seed = 20260425L + n_val + 100*sg)
  for (rho_val in c(0.80, 0.98))
    for (sg in c(0, 0.5))
      cells[[length(cells)+1L]] <- list(
        label = sprintf("rho=%.2f,sg=%.1f", rho_val, sg),
        n = ANCHOR$n, m = ANCHOR$m, rho = rho_val, K_true = ANCHOR$K_true,
        beta_true = ANCHOR$beta_true, sigma_g2 = sg,
        block_size = ANCHOR$block_size,
        seed = 20260425L + as.integer(round(1000*rho_val)) + 100*sg)
  # Signal-architecture cells skipped: under sigma_g^2 = 0.5 the eBIC stopping
  # rule is computationally pathological for ALL three signal levels (weak,
  # medium, strong) — best_log_BF stays positive for many spurious SNPs
  # at every step, so K saturates at K_max and MC integration becomes
  # the bottleneck. The 10 cells across the (n, ρ) grid are sufficient for
  # the sensitivity claim.

  cat(sprintf("=== Prior sensitivity: %d cells (half_cauchy vs inv_gamma) ===\n",
              length(cells)))

  # Checkpoint per cell: each worker writes its rows to a separate RDS
  ck_dir <- file.path("results/bench_full/07_prior_sensitivity_ck")
  dir.create(ck_dir, showWarnings = FALSE, recursive = TRUE)

  out_list <- parallel::mclapply(cells, function(cc) {
    ck_fn <- file.path(ck_dir, paste0("cell_", gsub("[^A-Za-z0-9]", "_",
                                                       cc$label), ".rds"))
    if (file.exists(ck_fn) && file.size(ck_fn) > 0) {
      message(sprintf("  %s  (cached)", cc$label))
      return(readRDS(ck_fn))
    }
    message(sprintf("  %s ...", cc$label))
    args <- cc; args$cell_label <- cc$label; args$label <- NULL
    res <- do.call(run_sensitivity_cell, args)
    saveRDS(res, ck_fn)
    message(sprintf("  %s done", cc$label))
    res
  }, mc.cores = min(8L, parallel::detectCores() - 2L),
     mc.preschedule = FALSE)
  out <- do.call(rbind, out_list)

  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
  write.csv(out, OUT_CSV, row.names = FALSE)
  cat(sprintf("\nWrote %s (%d rows).\n", OUT_CSV, nrow(out)))

  # Summary statistics
  cat("\n=== Summary ===\n")
  cat(sprintf("Mean |delta F1|:           %.4f\n", mean(abs(out$delta_F1))))
  cat(sprintf("Max  |delta F1|:           %.4f\n", max(abs(out$delta_F1))))
  cat(sprintf("Mean Jaccard(selected):    %.3f\n", mean(out$jaccard)))
  cat(sprintf("Min  Jaccard:              %.3f\n", min(out$jaccard)))
  cat(sprintf("Mean max(|delta log BF|):  %.3f\n",
              mean(out$max_dlogBF, na.rm = TRUE)))
  cat(sprintf("Cells with |delta logBF|<0.05: %d / %d\n",
              sum(out$max_dlogBF < 0.05, na.rm = TRUE), nrow(out)))
  cat(sprintf("Cells with |delta logBF|<0.20: %d / %d\n",
              sum(out$max_dlogBF < 0.20, na.rm = TRUE), nrow(out)))
}
