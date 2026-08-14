# ==============================================================================
# 13_mixtau_benchmark.R
#
# Benchmark of the mixture-slab variant (variant 4 of the methodological
# review) against the standard single-tau^2 framework.
#
# Setting: 5 focused cells × B=50 replicates.  Each cell tests MS_L and
# JS_L with three slab configurations:
#   * tau^2 = 0.04 (default, single value)
#   * tau^2 in {0.01, 0.04, 0.10, 0.25} mixture, uniform weights
#   * tau^2 in {0.04, 0.10, 0.25} mixture, polygenicity-biased weights
#
# Output: results/bench_full/13_mixtau/<cell>_b<rep>.rds
#         results/bench_full/13_mixtau_summary.csv
#
# Usage: Rscript sim/bench_full/13_mixtau_benchmark.R --cores 16
# ==============================================================================

source("sim/bench_full/00_config.R")
source("R/LMM_stepwise_mixtau.R")

OUT_DIR <- "results/bench_full/13_mixtau"
SUM_CSV <- "results/bench_full/13_mixtau_summary.csv"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

#  Cell list: 5 focused cells covering the polygenic dichotomy 
build_cells <- function() {
  list(
    list(tag = "anchor_sg0",        n = 1000L, m = 5000L, rho = 0.95,
          K_true = 5L, beta_true = c(0.8,0.4,0.4,0.2,0.2),
          sigma_g2 = 0,   seed_offset = 0L),
    list(tag = "anchor_sg0.5",      n = 1000L, m = 5000L, rho = 0.95,
          K_true = 5L, beta_true = c(0.8,0.4,0.4,0.2,0.2),
          sigma_g2 = 0.5, seed_offset = 100L),
    list(tag = "weak_sg0.5",        n = 1000L, m = 5000L, rho = 0.95,
          K_true = 5L, beta_true = c(0.4,0.2,0.2,0.1,0.1),
          sigma_g2 = 0.5, seed_offset = 200L),
    list(tag = "n3000_sg0.5",       n = 3000L, m = 5000L, rho = 0.95,
          K_true = 5L, beta_true = c(0.8,0.4,0.4,0.2,0.2),
          sigma_g2 = 0.5, seed_offset = 300L),
    list(tag = "rho098_sg0",        n = 1000L, m = 5000L, rho = 0.98,
          K_true = 5L, beta_true = c(0.8,0.4,0.4,0.2,0.2),
          sigma_g2 = 0,   seed_offset = 400L)
  )
}

#  Tau configurations to compare 
TAU_CONFIGS <- list(
  std_004      = list(name = "tau2=0.04 (default)",
                       tau2_mix = c(0.04),     weights = c(1)),
  mix_uniform  = list(name = "Mix uniform {0.01,0.04,0.10,0.25}",
                       tau2_mix = c(0.01, 0.04, 0.10, 0.25),
                       weights = c(0.25, 0.25, 0.25, 0.25)),
  mix_polybias = list(name = "Mix polybias {0.04,0.10,0.25}",
                       tau2_mix = c(0.04, 0.10, 0.25),
                       weights = c(0.4, 0.3, 0.3))
)

#  Run one (cell, rep, tau_config, variant) 
run_one <- function(cell, rep_id) {
  tag <- cell$tag
  fn <- file.path(OUT_DIR, sprintf("%s_b%02d.rds", tag, rep_id))
  if (file.exists(fn) && file.size(fn) > 0) return(invisible("cached"))

  seed <- 20260425L + 1000L * (rep_id - 1L) + cell$seed_offset
  d <- gen_dataset(n = cell$n, m = cell$m, rho = cell$rho,
                     K_true = cell$K_true, beta_true = cell$beta_true,
                     sigma_g2 = cell$sigma_g2,
                     block_size = ANCHOR$block_size, seed = seed)

  results <- list()
  for (variant in c("MS_L", "JS_L")) {
    for (conf_name in names(TAU_CONFIGS)) {
      conf <- TAU_CONFIGS[[conf_name]]
      method_tag <- sprintf("%s_%s", variant, conf_name)
      fn_run <- if (variant == "JS_L") JS_L_LMM_stepwise_fast_mixtau
                else MS_L_LMM_stepwise_fast_mixtau
      t0 <- Sys.time()
      res <- tryCatch(fn_run(d$y, d$X,
                                tau2_mix = conf$tau2_mix,
                                weights = conf$weights,
                                K_max = ANCHOR$K_max,
                                criterion = "eBIC",
                                theta = ANCHOR$theta,
                                n_nodes = ANCHOR$N_delta),
                        error = function(e) NULL)
      elapsed <- as.numeric(Sys.time() - t0, units = "secs")
      if (is.null(res)) {
        results[[method_tag]] <- list(
          variant = variant, config = conf_name, indices = integer(0),
          K_hat = 0L, tp = 0L, fp = 0L,
          recall = 0, precision = 0, f1 = 0,
          elapsed = elapsed)
        next
      }
      tp <- length(intersect(res$indices, d$truth))
      fp <- res$K_hat - tp
      rec <- tp / cell$K_true
      prec <- if (res$K_hat > 0) tp / res$K_hat else 0
      f1 <- if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0
      results[[method_tag]] <- list(
        variant = variant, config = conf_name,
        indices = res$indices, K_hat = res$K_hat,
        tp = tp, fp = fp, recall = rec, precision = prec, f1 = f1,
        elapsed = elapsed)
    }
  }
  payload <- list(cell_tag = tag, rep = rep_id,
                    sigma_g2 = cell$sigma_g2, truth = d$truth,
                    results = results)
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
    for (mtag in names(p$results)) {
      r <- p$results[[mtag]]
      rows[[length(rows)+1L]] <- data.frame(
        cell = p$cell_tag, rep = p$rep, sigma_g2 = p$sigma_g2,
        variant = r$variant, config = r$config,
        method_tag = mtag,
        f1 = r$f1, recall = r$recall, precision = r$precision,
        K_hat = r$K_hat, elapsed = r$elapsed,
        stringsAsFactors = FALSE)
    }
  }
  df <- do.call(rbind, rows)
  suppressPackageStartupMessages(library(dplyr))
  s <- df %>% group_by(cell, sigma_g2, variant, config) %>%
    summarise(mean_f1 = mean(f1), sd_f1 = sd(f1),
              mean_recall = mean(recall), mean_precision = mean(precision),
              mean_K_hat = mean(K_hat), mean_elapsed = mean(elapsed),
              n_reps = n(), .groups = "drop") %>%
    arrange(cell, variant, config)
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
  cat(sprintf("=== mixtau benchmark: %d cells x B=%d x 6 method-configs ===\n",
              length(cells), B))

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
