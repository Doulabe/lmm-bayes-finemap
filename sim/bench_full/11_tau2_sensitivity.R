# ==============================================================================
# 11_tau2_sensitivity.R
#
# Sensitivity of the proposed framework (MS_L_eBIC, JS_L_eBIC) to the slab
# variance tau^2 in the point-normal g-prior on the candidate effect.
#
# Setting: anchor cell (n=1000, m=5000, rho=0.95, K_true=5, beta=(0.8,0.4,0.4,
# 0.2,0.2)), crossed with sigma_g2 in {0, 0.5} and tau^2 in {0.01, 0.04,
# 0.10, 0.25}. B=20 replicates per cell.
#
# Output: results/bench_full/11_tau2_sensitivity/<cell>_b<rep>.rds
#         results/bench_full/11_tau2_sensitivity_summary.csv (aggregated)
#
# Usage:
#   Rscript sim/bench_full/11_tau2_sensitivity.R --cores 16
# ==============================================================================

source("sim/bench_full/00_config.R")

OUT_DIR <- "results/bench_full/11_tau2_sensitivity"
SUMMARY_CSV <- "results/bench_full/11_tau2_sensitivity_summary.csv"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

TAU2_GRID  <- c(0.01, 0.04, 0.10, 0.25)
SG_GRID    <- c(0, 0.5)
B          <- 20L


run_one <- function(tau2, sg, rep_id) {
  cell_tag <- sprintf("tau%.2f_sg%.1f", tau2, sg)
  fn <- file.path(OUT_DIR, sprintf("%s_b%02d.rds", cell_tag, rep_id))
  if (file.exists(fn) && file.size(fn) > 0) return(invisible("cached"))

  seed <- 20260425L + 1000L * (rep_id - 1L) +
            as.integer(round(1000 * tau2)) + as.integer(round(100 * sg))
  d <- gen_dataset(n = ANCHOR$n, m = ANCHOR$m, rho = ANCHOR$rho,
                     K_true = ANCHOR$K_true, beta_true = ANCHOR$beta_true,
                     sigma_g2 = sg, block_size = ANCHOR$block_size,
                     seed = seed)

  res_list <- list()
  for (method in c("MS_L_eBIC", "JS_L_eBIC")) {
    fn_run <- if (startsWith(method, "JS_L")) JS_L_LMM_stepwise_fast
              else MS_L_LMM_stepwise_fast
    t0 <- Sys.time()
    res <- fn_run(d$y, d$X, tau2 = tau2, K_max = ANCHOR$K_max,
                    criterion = "eBIC",
                    theta = ANCHOR$theta, n_nodes = ANCHOR$N_delta)
    elapsed <- as.numeric(Sys.time() - t0, units = "secs")
    tp <- length(intersect(res$indices, d$truth))
    fp <- res$K_hat - tp
    rec <- tp / ANCHOR$K_true
    prec <- if (res$K_hat > 0) tp / res$K_hat else 0
    f1 <- if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0
    res_list[[method]] <- list(
      method = method, tau2 = tau2, sigma_g2 = sg, rep = rep_id,
      indices = res$indices, K_hat = res$K_hat, tp = tp, fp = fp,
      recall = rec, precision = prec, f1 = f1, elapsed = elapsed)
  }
  payload <- list(cell_tag = cell_tag, rep = rep_id,
                   tau2 = tau2, sigma_g2 = sg, truth = d$truth,
                   results = res_list)
  saveRDS(payload, fn)
  invisible("done")
}


aggregate_results <- function() {
  files <- list.files(OUT_DIR, pattern = "\\.rds$", full.names = TRUE)
  if (length(files) == 0L) return(invisible(NULL))
  rows <- list()
  for (f in files) {
    p <- readRDS(f)
    for (m_name in names(p$results)) {
      r <- p$results[[m_name]]
      rows[[length(rows)+1L]] <- data.frame(
        tau2 = r$tau2, sigma_g2 = r$sigma_g2, method = m_name, rep = r$rep,
        f1 = r$f1, recall = r$recall, precision = r$precision,
        K_hat = r$K_hat, elapsed = r$elapsed,
        stringsAsFactors = FALSE)
    }
  }
  df <- do.call(rbind, rows)
  suppressPackageStartupMessages(library(dplyr))
  s <- df %>% group_by(method, tau2, sigma_g2) %>%
    summarise(mean_f1 = mean(f1, na.rm = TRUE),
              sd_f1 = sd(f1, na.rm = TRUE),
              mean_recall = mean(recall, na.rm = TRUE),
              mean_precision = mean(precision, na.rm = TRUE),
              mean_K_hat = mean(K_hat, na.rm = TRUE),
              n_reps = n(),
              .groups = "drop") %>%
    arrange(sigma_g2, method, tau2)
  write.csv(s, SUMMARY_CSV, row.names = FALSE)
  cat(sprintf("Wrote %s\n", SUMMARY_CSV))
  print(as.data.frame(s))
  invisible(s)
}


if (sys.nframe() == 0L) {
  Sys.setenv(OMP_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
              OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
  if (!requireNamespace("rrBLUP", quietly = TRUE))
    stop("rrBLUP needed")
  loadNamespace("rrBLUP")

  args <- parse_args()
  n_cores <- if (!is.null(args$cores)) as.integer(args$cores) else 8L
  cat(sprintf("=== tau^2 sensitivity: cores = %d, B = %d ===\n", n_cores, B))

  jobs <- list()
  for (tau2 in TAU2_GRID)
    for (sg in SG_GRID)
      for (rep_id in seq_len(B))
        jobs[[length(jobs)+1L]] <- list(tau2 = tau2, sg = sg, rep = rep_id)
  cat(sprintf("  %d jobs (4 tau^2 x 2 sg x %d reps)\n", length(jobs), B))

  parallel::mclapply(jobs, function(j) {
    tryCatch(run_one(j$tau2, j$sg, j$rep),
              error = function(e) {
                message("ERR tau=", j$tau2, " sg=", j$sg, " b", j$rep,
                          ": ", conditionMessage(e))
                "error"
              })
  }, mc.cores = n_cores, mc.preschedule = FALSE)

  cat("\n=== Aggregating ===\n")
  aggregate_results()
}
