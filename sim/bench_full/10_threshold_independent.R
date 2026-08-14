# ==============================================================================
# 10_threshold_independent.R
#
# Threshold-independent comparison of methods:
#   * Top-K* (K*=5) recall/precision: pick each method's top-5 by per-SNP score
#   * PR-AUC: precision-recall area under curve from full per-SNP ranking
#   * Brier score: calibration of posterior-probability-based methods
#
# Re-runs each method at the anchor cell + axis-variation cells, B=100 reps,
# extracting per-SNP scores for ALL m candidates (not just selected).
#
# Output: results/bench_full/10_threshold_indep/<cell>_b<rep>.rds
# Each RDS stores: per-SNP scores for {MS_L step-1 q^Bayes, JS_L step-1 q^Bayes,
#   SuSiE PIPs, BSLMM PIPs, BayesR PIPs, fastlmm |Z|}, plus truth indices.
#
# Usage:
#   Rscript sim/bench_full/10_threshold_independent.R --cores 16
#   Rscript sim/bench_full/10_threshold_independent.R --cell anchor --rep 1
# ==============================================================================

source("sim/bench_full/00_config.R")

OUT_DIR <- "results/bench_full/10_threshold_indep"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -- Per-SNP score extractors -------------------------------------------------
#
# Design note (important!): the proposed framework is intrinsically stepwise
# --- its advantage comes from absorbing previously selected SNPs into the
# LMM kernel R_jk at each step.  A step-1 single-pass scoring of all
# candidates incurs the same proximal contamination as FaST-LMM (the kernel
# contains the causal SNPs).  We therefore report Top-K* by running the
# stepwise algorithm with K_max = K_star and treating the K_star selected
# SNPs as the top-K_star ranking.  PR-AUC and Brier are evaluated only on
# methods that produce a full per-SNP PIP vector (SuSiE, BSLMM, BayesR).
# --------------------------------------------------------------------------

# Our methods: run stepwise to K_max = K_star, return the K_star selected
# indices as the "top-K_star ranking".
topK_MSL <- function(y, X, K_star, tau2 = 0.04, n_nodes = 15L) {
  res <- tryCatch(MS_L_LMM_stepwise_fast(y, X, tau2 = tau2,
                                            K_max = K_star,
                                            criterion = "JointPosterior",
                                            theta = 0,     # force K_star steps
                                            n_nodes = n_nodes),
                    error = function(e) NULL)
  if (is.null(res)) return(integer(0))
  res$indices
}

topK_JSL <- function(y, X, K_star, tau2 = 0.04, n_nodes = 15L) {
  res <- tryCatch(JS_L_LMM_stepwise_fast(y, X, tau2 = tau2,
                                            K_max = K_star,
                                            criterion = "JointPosterior",
                                            theta = 0,
                                            n_nodes = n_nodes),
                    error = function(e) NULL)
  if (is.null(res)) return(integer(0))
  res$indices
}


# SuSiE: extract per-SNP PIP (max over L credible sets)
score_susie <- function(y, X, L = 10) {
  fit <- tryCatch(susieR::susie(X, y, L = L, verbose = FALSE),
                    error = function(e) NULL)
  if (is.null(fit)) return(rep(NA_real_, ncol(X)))
  as.numeric(susieR::susie_get_pip(fit))
}

# BSLMM via BGLR: extract per-SNP probabilities
score_bslmm <- function(y, X, nIter = 3000L, burnIn = 500L) {
  if (!requireNamespace("BGLR", quietly = TRUE))
    return(rep(NA_real_, ncol(X)))
  fit <- tryCatch({
    suppressMessages(suppressWarnings({
      eta <- list(list(X = X, model = "BayesB"))
      BGLR::BGLR(y = y, ETA = eta, nIter = nIter, burnIn = burnIn,
                  verbose = FALSE)
    }))
  }, error = function(e) NULL)
  if (is.null(fit) || is.null(fit$ETA[[1]]$d)) return(rep(NA_real_, ncol(X)))
  as.numeric(fit$ETA[[1]]$d)
}

# BayesR (hibayes::ibrm method="BayesR"): use fit$pip directly
score_bayesR <- function(y, X, niter = 2000L, nburn = 400L) {
  if (!requireNamespace("hibayes", quietly = TRUE))
    return(rep(NA_real_, ncol(X)))
  n_obs <- nrow(X); Xm <- X
  rownames(Xm) <- paste0("S", seq_len(n_obs))
  df_y <- data.frame(id = rownames(Xm), y = y, stringsAsFactors = FALSE)
  fit <- tryCatch(suppressMessages(suppressWarnings(
            hibayes::ibrm(formula = y ~ 1, data = df_y, M = Xm, M.id = df_y$id,
                           method = "BayesR",
                           Pi = c(0.95, 0.02, 0.02, 0.01),
                           fold = c(0, 1e-4, 1e-3, 1e-2),
                           niter = niter, nburn = nburn,
                           threads = 1L, verbose = FALSE))),
          error = function(e) NULL)
  if (is.null(fit) || is.null(fit$pip)) return(rep(NA_real_, ncol(X)))
  as.numeric(fit$pip)
}

# Bayesian Lasso via BGLR (model = "BL"); score = posterior mean |beta_hat|
# (continuous double-exponential prior lead to no exact zero, no PIP).  Used as a
# ranking score for threshold-independent metrics; treated similarly to
# FaST-LMM's |Z| score: Brier is NA (no calibrated inclusion probability).
score_BL <- function(y, X, nIter = 2000L, burnIn = 400L) {
  if (!requireNamespace("BGLR", quietly = TRUE))
    return(rep(NA_real_, ncol(X)))
  fit <- tryCatch({
    suppressMessages(suppressWarnings({
      eta <- list(list(X = X, model = "BL"))
      BGLR::BGLR(y = y, ETA = eta, nIter = nIter, burnIn = burnIn,
                 verbose = FALSE)
    }))
  }, error = function(e) NULL)
  if (is.null(fit) || is.null(fit$ETA[[1]]$b)) return(rep(NA_real_, ncol(X)))
  abs(as.numeric(fit$ETA[[1]]$b))
}

# FaST-LMM: use |Z| as score
score_fastlmm <- function(y, X) {
  if (!requireNamespace("rrBLUP", quietly = TRUE))
    return(rep(NA_real_, ncol(X)))
  n <- nrow(X); m <- ncol(X)
  K <- tcrossprod(X) / max(m - 1L, 1L)
  fit <- tryCatch(rrBLUP::mixed.solve(y = y, K = K),
                    error = function(e) NULL)
  if (is.null(fit)) return(rep(NA_real_, m))
  delta_hat <- fit$Vu / fit$Ve
  ee <- eigen(K, symmetric = TRUE)
  d_inv_sqrt <- 1 / sqrt(1 + delta_hat * ee$values)
  UtX <- crossprod(ee$vectors, X)
  Uty <- crossprod(ee$vectors, y)
  Z <- as.numeric(crossprod(UtX * d_inv_sqrt, Uty * d_inv_sqrt) /
                     sqrt(colSums((UtX * d_inv_sqrt)^2) + 1e-300))
  abs(Z)
}

# -- Threshold-independent metrics --------------------------------------------

# Top-K precision/recall: take top K SNPs by score
topK_recall <- function(scores, truth, K) {
  if (all(is.na(scores))) return(NA_real_)
  ord <- order(-scores, na.last = TRUE)
  top <- ord[seq_len(min(K, length(scores)))]
  length(intersect(top, truth)) / length(truth)
}

# PR-AUC via trapezoidal integration on sorted scores
pr_auc <- function(scores, truth, n_total) {
  if (all(is.na(scores))) return(NA_real_)
  scores[is.na(scores)] <- -Inf
  ord <- order(-scores)
  is_truth <- seq_along(scores) %in% truth
  tp_cum <- cumsum(is_truth[ord])
  fp_cum <- cumsum(!is_truth[ord])
  precision <- tp_cum / pmax(tp_cum + fp_cum, 1)
  recall <- tp_cum / length(truth)
  # Trapezoidal integration
  recall  <- c(0, recall)
  precision <- c(1, precision)
  sum(diff(recall) * (precision[-1] + precision[-length(precision)]) / 2)
}

# Brier score: only meaningful for methods producing probabilities in [0,1]
brier_score <- function(probs, truth, n_total) {
  if (all(is.na(probs))) return(NA_real_)
  probs <- pmin(pmax(probs, 0), 1)        # clip to [0,1]
  is_truth <- as.numeric(seq_along(probs) %in% truth)
  mean((probs - is_truth)^2)
}

# Expected Calibration Error (ECE) with 10 bins
ece <- function(probs, truth, n_bins = 10) {
  if (all(is.na(probs))) return(NA_real_)
  probs <- pmin(pmax(probs, 0), 1)
  is_truth <- as.numeric(seq_along(probs) %in% truth)
  bins <- cut(probs, breaks = seq(0, 1, length.out = n_bins + 1L),
                include.lowest = TRUE)
  out <- 0
  for (b in levels(bins)) {
    in_bin <- which(bins == b)
    if (length(in_bin) == 0) next
    avg_p <- mean(probs[in_bin])
    avg_y <- mean(is_truth[in_bin])
    out <- out + (length(in_bin) / length(probs)) * abs(avg_p - avg_y)
  }
  out
}

# -- Run one cell -------------------------------------------------------------

run_cell_threshold_indep <- function(n, m, rho, K_true, beta_true, sigma_g2,
                                       block_size, seed, cell_tag, rep_id) {
  fn <- file.path(OUT_DIR, sprintf("%s_b%02d.rds", cell_tag, rep_id))
  if (file.exists(fn) && file.size(fn) > 0) {
    cat(sprintf("  [skip] %s b%02d already done\n", cell_tag, rep_id))
    return(invisible(NULL))
  }
  d <- gen_dataset(n = n, m = m, rho = rho, K_true = K_true,
                     beta_true = beta_true, sigma_g2 = sigma_g2,
                     block_size = block_size, seed = seed)
  truth <- d$truth

  K_star <- K_true
  cat(sprintf("  %s b%02d: scoring...\n", cell_tag, rep_id))

  # Per-SNP scores for PIP-based methods
  pip_susie   <- score_susie(d$y, d$X, L = 10L)
  pip_bslmm   <- score_bslmm(d$y, d$X, nIter = 3000L, burnIn = 500L)
  pip_bayesR  <- score_bayesR(d$y, d$X, niter = 2000L, nburn = 400L)
  beta_BL     <- score_BL(d$y, d$X, nIter = 2000L, burnIn = 400L)
  abs_z_fast  <- score_fastlmm(d$y, d$X)

  # Top-K* selected indices for stepwise methods (no scoring, just selection)
  top_MSL <- topK_MSL(d$y, d$X, K_star = K_star)
  top_JSL <- topK_JSL(d$y, d$X, K_star = K_star)

  out <- list(
    MS_L_eBIC = list(top_K_recall = length(intersect(top_MSL, truth)) / K_star,
                       pr_auc = NA_real_, brier = NA_real_, ece = NA_real_),
    JS_L_eBIC = list(top_K_recall = length(intersect(top_JSL, truth)) / K_star,
                       pr_auc = NA_real_, brier = NA_real_, ece = NA_real_),
    SuSiE = list(top_K_recall = topK_recall(pip_susie, truth, K_star),
                   pr_auc = pr_auc(pip_susie, truth, n_total = m),
                   brier = brier_score(pip_susie, truth, n_total = m),
                   ece = ece(pip_susie, truth, n_bins = 10L)),
    BSLMM = list(top_K_recall = topK_recall(pip_bslmm, truth, K_star),
                   pr_auc = pr_auc(pip_bslmm, truth, n_total = m),
                   brier = brier_score(pip_bslmm, truth, n_total = m),
                   ece = ece(pip_bslmm, truth, n_bins = 10L)),
    BayesR = list(top_K_recall = topK_recall(pip_bayesR, truth, K_star),
                    pr_auc = pr_auc(pip_bayesR, truth, n_total = m),
                    brier = brier_score(pip_bayesR, truth, n_total = m),
                    ece = ece(pip_bayesR, truth, n_bins = 10L)),
    BL = list(top_K_recall = topK_recall(beta_BL, truth, K_star),
                pr_auc = pr_auc(beta_BL, truth, n_total = m),
                brier = NA_real_, ece = NA_real_),
    fastlmm = list(top_K_recall = topK_recall(abs_z_fast, truth, K_star),
                     pr_auc = pr_auc(abs_z_fast, truth, n_total = m),
                     brier = NA_real_, ece = NA_real_)
  )

  payload <- list(cell_tag = cell_tag, rep = rep_id,
                    n = n, m = m, rho = rho, sigma_g2 = sigma_g2,
                    K_true = K_true, truth = truth,
                    metrics = out)
  saveRDS(payload, fn)
  cat(sprintf("  [done] %s b%02d\n", cell_tag, rep_id))
  invisible(NULL)
}

# -- Cell list: anchor + 4 axis variations × 2 sg = 10 cells × B=100 = 1000 ---

build_cells <- function() {
  cells <- list()
  # Anchor: n=1000, m=5000, rho=0.95, signal=medium
  for (sg in c(0, 0.5))
    cells[[length(cells)+1L]] <- list(
      tag = sprintf("anchor_sg%.1f", sg),
      n = 1000L, m = 5000L, rho = 0.95,
      K_true = 5L, beta_true = c(0.8, 0.4, 0.4, 0.2, 0.2),
      sigma_g2 = sg, block_size = 10L,
      seed_offset = 0L)
  # n=3000 variation
  for (sg in c(0, 0.5))
    cells[[length(cells)+1L]] <- list(
      tag = sprintf("n3000_sg%.1f", sg),
      n = 3000L, m = 5000L, rho = 0.95,
      K_true = 5L, beta_true = c(0.8, 0.4, 0.4, 0.2, 0.2),
      sigma_g2 = sg, block_size = 10L,
      seed_offset = 1L)
  # rho=0.98 variation
  for (sg in c(0, 0.5))
    cells[[length(cells)+1L]] <- list(
      tag = sprintf("rho098_sg%.1f", sg),
      n = 1000L, m = 5000L, rho = 0.98,
      K_true = 5L, beta_true = c(0.8, 0.4, 0.4, 0.2, 0.2),
      sigma_g2 = sg, block_size = 10L,
      seed_offset = 2L)
  cells
}


if (sys.nframe() == 0L) {
  Sys.setenv(OMP_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
              OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
  if (!requireNamespace("rrBLUP", quietly = TRUE))
    stop("rrBLUP needed")
  loadNamespace("rrBLUP")
  loadNamespace("susieR")
  loadNamespace("BGLR")

  args <- parse_args()
  n_cores <- if (!is.null(args$cores)) as.integer(args$cores) else 8L
  B <- if (!is.null(args$B)) as.integer(args$B) else 100L

  cells <- build_cells()
  cat(sprintf("=== Threshold-independent: %d cells x B=%d ===\n",
              length(cells), B))

  # Build (cell, rep) job list
  jobs <- list()
  for (cc in cells) {
    for (b in seq_len(B)) {
      seed <- 20260425L + 1000L * (b - 1L) + cc$seed_offset
      jobs[[length(jobs)+1L]] <- c(cc, list(rep = b, seed = seed))
    }
  }
  cat(sprintf("  %d jobs total\n", length(jobs)))

  parallel::mclapply(jobs, function(job) {
    tryCatch(run_cell_threshold_indep(
      n = job$n, m = job$m, rho = job$rho,
      K_true = job$K_true, beta_true = job$beta_true,
      sigma_g2 = job$sigma_g2, block_size = job$block_size,
      seed = job$seed, cell_tag = job$tag, rep_id = job$rep),
      error = function(e) {
        message("ERR ", job$tag, " b", job$rep, ": ", conditionMessage(e))
      })
  }, mc.cores = n_cores, mc.preschedule = FALSE)
  cat("\nDone.\n")
}
