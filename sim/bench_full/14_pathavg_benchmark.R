# ==============================================================================
# 14_pathavg_benchmark.R
#
# Benchmark of the path-averaged (PA) variant against standard stepwise.
# Tests whether maintaining a beam of B paths and computing marginal PIPs
# improves calibration (Brier score) without sacrificing parsimony.
#
# Setting: 5 focused cells × B_reps=50 replicates × {PA-MS_L, PA-JS_L,
# standard MS_L for reference} × beam width B = 10.
#
# Output: results/bench_full/14_pathavg/<cell>_b<rep>.rds
#         results/bench_full/14_pathavg_summary.csv
#
# Usage: Rscript sim/bench_full/14_pathavg_benchmark.R --cores 16
# ==============================================================================

source("sim/bench_full/00_config.R")
source("R/LMM_stepwise_pathavg.R")

OUT_DIR <- "results/bench_full/14_pathavg"
SUM_CSV <- "results/bench_full/14_pathavg_summary.csv"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

BEAM_WIDTH <- 10L
K_MAX      <- 8L     # K* = 5 + 3 buffer
TAU2       <- 0.04
N_NODES    <- 15L

# Phase 2 (focused validation): 2 anchor cells only.
# If PA shows the predicted Brier reduction (~50%), we extend to more cells.
build_cells <- function() {
  list(
    list(tag = "anchor_sg0",   n = 1000L, m = 5000L, rho = 0.95,
          K_true = 5L, beta_true = c(0.8,0.4,0.4,0.2,0.2),
          sigma_g2 = 0,   seed_offset = 0L),
    list(tag = "anchor_sg0.5", n = 1000L, m = 5000L, rho = 0.95,
          K_true = 5L, beta_true = c(0.8,0.4,0.4,0.2,0.2),
          sigma_g2 = 0.5, seed_offset = 100L)
  )
}

#  Metrics 
top_K_recall <- function(scores, truth, K) {
  if (all(is.na(scores))) return(NA_real_)
  ord <- order(-scores, na.last = TRUE)
  length(intersect(ord[seq_len(min(K, length(scores)))], truth)) / length(truth)
}
pr_auc <- function(scores, truth) {
  scores[is.na(scores)] <- -Inf
  ord <- order(-scores)
  is_truth <- seq_along(scores) %in% truth
  tp_cum <- cumsum(is_truth[ord]); fp_cum <- cumsum(!is_truth[ord])
  precision <- tp_cum / pmax(tp_cum + fp_cum, 1)
  recall <- tp_cum / length(truth)
  recall <- c(0, recall); precision <- c(1, precision)
  sum(diff(recall) * (precision[-1] + precision[-length(precision)]) / 2)
}
brier_score <- function(probs, truth) {
  probs <- pmin(pmax(probs, 0), 1); probs[is.na(probs)] <- 0
  is_truth <- as.numeric(seq_along(probs) %in% truth)
  mean((probs - is_truth)^2)
}

#  Run one 
run_one <- function(cell, rep_id) {
  fn <- file.path(OUT_DIR, sprintf("%s_b%02d.rds", cell$tag, rep_id))
  if (file.exists(fn) && file.size(fn) > 0) return(invisible("cached"))

  seed <- 20260425L + 1000L * (rep_id - 1L) + cell$seed_offset
  d <- gen_dataset(n = cell$n, m = cell$m, rho = cell$rho,
                     K_true = cell$K_true, beta_true = cell$beta_true,
                     sigma_g2 = cell$sigma_g2,
                     block_size = ANCHOR$block_size, seed = seed)

  results <- list()

  # Standard MS_L (single greedy path) as baseline for paired comparison
  t0 <- Sys.time()
  res_std_ms <- tryCatch(MS_L_LMM_stepwise_fast(d$y, d$X,
                                                   tau2 = TAU2,
                                                   K_max = K_MAX,
                                                   criterion = "eBIC",
                                                   theta = 0.99,
                                                   n_nodes = N_NODES),
                           error = function(e) NULL)
  t_std_ms <- as.numeric(Sys.time() - t0, units = "secs")
  if (!is.null(res_std_ms)) {
    sel <- res_std_ms$indices
    tp <- length(intersect(sel, d$truth))
    rec <- tp / cell$K_true
    prec <- if (length(sel) > 0L) tp / length(sel) else 0
    f1 <- if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0
    # For standard, q at selection step = score
    q_std <- numeric(ncol(d$X))
    q_std[res_std_ms$indices] <- res_std_ms$q_at_step
    results$std_MS_L <- list(
      indices = sel, K_hat = length(sel), tp = tp,
      recall = rec, precision = prec, f1 = f1,
      top_K_recall = top_K_recall(q_std, d$truth, cell$K_true),
      pr_auc = pr_auc(q_std, d$truth),
      brier = brier_score(q_std, d$truth),
      elapsed = t_std_ms
    )
  }

  # PA-MS_L
  t0 <- Sys.time()
  res_pa_ms <- tryCatch(PA_MS_L_LMM_stepwise_fast(d$y, d$X,
                                                    tau2 = TAU2,
                                                    B = BEAM_WIDTH,
                                                    K_max = K_MAX,
                                                    n_nodes = N_NODES),
                          error = function(e) NULL)
  t_pa_ms <- as.numeric(Sys.time() - t0, units = "secs")
  if (!is.null(res_pa_ms)) {
    q <- res_pa_ms$q_marginal
    sel <- res_pa_ms$indices
    tp <- length(intersect(sel, d$truth))
    rec <- tp / cell$K_true
    prec <- if (length(sel) > 0L) tp / length(sel) else 0
    f1 <- if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0
    results$PA_MS_L <- list(
      indices = sel, K_hat = length(sel), tp = tp,
      recall = rec, precision = prec, f1 = f1,
      top_K_recall = top_K_recall(q, d$truth, cell$K_true),
      pr_auc = pr_auc(q, d$truth),
      brier = brier_score(q, d$truth),
      best_path = res_pa_ms$best_path_indices,
      elapsed = t_pa_ms
    )
  }

  # PA-JS_L
  t0 <- Sys.time()
  res_pa_js <- tryCatch(PA_JS_L_LMM_stepwise_fast(d$y, d$X,
                                                    tau2 = TAU2,
                                                    B = BEAM_WIDTH,
                                                    K_max = K_MAX,
                                                    n_nodes = N_NODES),
                          error = function(e) NULL)
  t_pa_js <- as.numeric(Sys.time() - t0, units = "secs")
  if (!is.null(res_pa_js)) {
    q <- res_pa_js$q_marginal
    sel <- res_pa_js$indices
    tp <- length(intersect(sel, d$truth))
    rec <- tp / cell$K_true
    prec <- if (length(sel) > 0L) tp / length(sel) else 0
    f1 <- if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0
    results$PA_JS_L <- list(
      indices = sel, K_hat = length(sel), tp = tp,
      recall = rec, precision = prec, f1 = f1,
      top_K_recall = top_K_recall(q, d$truth, cell$K_true),
      pr_auc = pr_auc(q, d$truth),
      brier = brier_score(q, d$truth),
      best_path = res_pa_js$best_path_indices,
      elapsed = t_pa_js
    )
  }

  payload <- list(cell_tag = cell$tag, rep = rep_id,
                    sigma_g2 = cell$sigma_g2, K_true = cell$K_true,
                    truth = d$truth, results = results)
  saveRDS(payload, fn)
  invisible("done")
}

#  Aggregate 
aggregate_summary <- function() {
  files <- list.files(OUT_DIR, pattern = "\\.rds$", full.names = TRUE)
  if (length(files) == 0L) return(invisible(NULL))
  rows <- list()
  for (f in files) {
    p <- readRDS(f)
    for (m_name in names(p$results)) {
      r <- p$results[[m_name]]
      rows[[length(rows)+1L]] <- data.frame(
        cell = p$cell_tag, rep = p$rep, sigma_g2 = p$sigma_g2,
        method = m_name,
        K_hat = r$K_hat, f1 = r$f1, recall = r$recall, precision = r$precision,
        top_K = r$top_K_recall, pr_auc = r$pr_auc, brier = r$brier,
        elapsed = r$elapsed,
        stringsAsFactors = FALSE)
    }
  }
  df <- do.call(rbind, rows)
  suppressPackageStartupMessages(library(dplyr))
  s <- df %>% group_by(cell, sigma_g2, method) %>%
    summarise(mean_f1 = mean(f1), mean_recall = mean(recall),
              mean_precision = mean(precision),
              mean_top_K = mean(top_K), mean_pr_auc = mean(pr_auc),
              mean_brier = mean(brier),
              mean_K_hat = mean(K_hat),
              mean_elapsed = mean(elapsed),
              n_reps = n(), .groups = "drop") %>%
    arrange(cell, method)
  write.csv(s, SUM_CSV, row.names = FALSE)
  cat(sprintf("Wrote %s\n", SUM_CSV))
  print(as.data.frame(s))
  invisible(s)
}


if (sys.nframe() == 0L) {
  Sys.setenv(OMP_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
              OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
  args <- parse_args()
  n_cores <- if (!is.null(args$cores)) as.integer(args$cores) else 8L
  B <- if (!is.null(args$B)) as.integer(args$B) else 50L
  cells <- build_cells()
  cat(sprintf("=== PA benchmark: %d cells x B=%d, beam=%d ===\n",
              length(cells), B, BEAM_WIDTH))

  jobs <- list()
  for (cc in cells)
    for (b in seq_len(B))
      jobs[[length(jobs)+1L]] <- list(cell = cc, rep_id = b)
  cat(sprintf("  %d jobs total\n", length(jobs)))

  parallel::mclapply(jobs, function(j)
    tryCatch(run_one(j$cell, j$rep_id),
              error = function(e) message("ERR ", j$cell$tag, " b", j$rep_id,
                                              ": ", conditionMessage(e))),
    mc.cores = n_cores, mc.preschedule = FALSE)

  cat("\n=== Aggregating ===\n")
  aggregate_summary()
}
