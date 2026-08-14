# ==============================================================================
# 00_config.R
# Common configuration, data generators, method runners, and metric helpers
# for the full benchmark suite.
#
# All scripts (01-05) source this file. Anchor configuration:
#   - K_true = 5, beta = (0.8, 0.4, 0.4, 0.2, 0.2)  -- 1 strong + 2 mod + 2 weak
#   - n_anchor = 1000, m_anchor = 5000, rho_anchor = 0.95
#   - sigma_g2 in {0, 0.5}, two polygenic regimes
#   - K_max = 10, theta = 0.99 (JP99), N_delta = 15
#   - JS_L/MS_L are run under both eBIC AND JP99 stopping
#     -> 8 method labels in total
#   - B  replicates per cell
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(susieR)
})

# Sourced by the orchestrator scripts; tries multiple paths for portability.
# Supports CWD = project root, sim/, or sim/bench_full/.
.try_source <- function(paths) {
  for (p in paths) if (file.exists(p)) { source(p); return(invisible(p)) }
  stop("Could not locate any of: ", paste(paths, collapse = ", "),
       "\n  CWD = ", getwd())
}
.try_source(c("R/LMM_core.R",
              "../R/LMM_core.R",
              "../../R/LMM_core.R"))
.try_source(c("R/LMM_stepwise_fast.R",
              "../R/LMM_stepwise_fast.R",
              "../../R/LMM_stepwise_fast.R"))


# ------------------------------------------------------------------------------
# Anchor configuration
# ------------------------------------------------------------------------------
ANCHOR <- list(
  n         = 1000L,
  m         = 5000L,
  rho       = 0.95,
  K_true    = 5L,
  beta_true = c(0.8, 0.4, 0.4, 0.2, 0.2),    # 1 strong + 2 mod + 2 weak
  sigma_g2  = 0.0,
  block_size = 10L,
  K_max     = 10L,
  theta     = 0.99,
  N_delta   = 15L,
  tau2      = 0.04,
  B         = 10L
)


# ------------------------------------------------------------------------------
# Data generator: block-AR(1) X with optional polygenic background
# ------------------------------------------------------------------------------
gen_dataset <- function(n, m, rho, K_true, beta_true,
                        sigma_g2 = 0, block_size = 10L,
                        seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Z <- matrix(rnorm(n * m), n, m)
  n_blocks <- m %/% block_size
  ar1_block <- function(p, r) outer(seq_len(p), seq_len(p),
                                     function(a, b) r^abs(a - b))
  S <- ar1_block(block_size, rho)
  L <- chol(S + 1e-8 * diag(block_size))
  for (bk in seq_len(n_blocks)) {
    cols <- ((bk - 1) * block_size + 1):(bk * block_size)
    Z[, cols] <- Z[, cols] %*% L
  }
  X <- scale(Z)

  # Place K_true causals at distinct positions (one per block to avoid
  # within-block degeneracy under high rho)
  truth <- seq(1, by = block_size, length.out = K_true)

  # Polygenic background from non-causal SNPs, if requested
  if (sigma_g2 > 0) {
    non_causal <- setdiff(seq_len(m), truth)
    R <- tcrossprod(X[, non_causal]) / length(non_causal)
    ee <- eigen(R, symmetric = TRUE)
    half <- ee$vectors %*%
            diag(sqrt(pmax(sigma_g2 * ee$values, 0))) %*%
            t(ee$vectors)
    u <- as.numeric(half %*% rnorm(n))
  } else {
    u <- rep(0, n)
  }

  y <- as.numeric(X[, truth] %*% beta_true) + u + rnorm(n)
  list(X = X, y = y, truth = truth, sigma_g2 = sigma_g2)
}


# ------------------------------------------------------------------------------
# Method runners
# Each runner returns: list(method, indices, K_hat, scores, elapsed)
# `methods` is a character vector subset of:
#   c("JS_L_eBIC", "JS_L_JP99",
#     "MS_L_eBIC", "MS_L_JP99",
#     "SuSiE", "BSLMM", "BayesR", "BL",
#     "fastlmm", "REML_stepwise")
# Note: "BL" = Bayesian Lasso via BGLR::BGLR with model="BL"; uses |beta_hat|
# as score, top-K selection (no PIP since prior is continuous double-exponential).
# Backward-compatible aliases:
#   "JS_L" -> "JS_L_JP99"  (the legacy posterior-threshold variant, but at 0.99)
#   "MS_L" -> "MS_L_JP99"
# ------------------------------------------------------------------------------
run_methods <- function(d, methods,
                          K_max = 10L, theta = 0.99, N_delta = 15L,
                          tau2 = 0.04) {
  out <- list()

  # Resolve aliases
  if ("JS_L" %in% methods)
    methods <- c(methods, "JS_L_eBIC", "JS_L_JP99")
  if ("MS_L" %in% methods)
    methods <- c(methods, "MS_L_eBIC", "MS_L_JP99")
  methods <- unique(setdiff(methods, c("JS_L", "MS_L")))

  if ("JS_L_eBIC" %in% methods) {
    t0 <- Sys.time()
    res <- JS_L_LMM_stepwise_fast(d$y, d$X, tau2 = tau2, K_max = K_max,
                                    criterion = "eBIC",
                                    theta = theta, n_nodes = N_delta)
    out$JS_L_eBIC <- list(method = "JS_L_eBIC", indices = res$indices,
                       K_hat = res$K_hat, scores = res$q_at_step,
                       elapsed = as.numeric(Sys.time() - t0, units = "secs"))
  }
  if ("JS_L_JP99" %in% methods) {
    t0 <- Sys.time()
    res <- JS_L_LMM_stepwise_fast(d$y, d$X, tau2 = tau2, K_max = K_max,
                                    criterion = "JointPosterior",
                                    theta = theta, n_nodes = N_delta)
    out$JS_L_JP99 <- list(method = "JS_L_JP99", indices = res$indices,
                       K_hat = res$K_hat, scores = res$q_at_step,
                       elapsed = as.numeric(Sys.time() - t0, units = "secs"))
  }

  if ("MS_L_eBIC" %in% methods) {
    t0 <- Sys.time()
    res <- MS_L_LMM_stepwise_fast(d$y, d$X, tau2 = tau2, K_max = K_max,
                                    criterion = "eBIC",
                                    theta = theta, n_nodes = N_delta)
    out$MS_L_eBIC <- list(method = "MS_L_eBIC", indices = res$indices,
                       K_hat = res$K_hat, scores = res$q_at_step,
                       elapsed = as.numeric(Sys.time() - t0, units = "secs"))
  }
  if ("MS_L_JP99" %in% methods) {
    t0 <- Sys.time()
    res <- MS_L_LMM_stepwise_fast(d$y, d$X, tau2 = tau2, K_max = K_max,
                                    criterion = "JointPosterior",
                                    theta = theta, n_nodes = N_delta)
    out$MS_L_JP99 <- list(method = "MS_L_JP99", indices = res$indices,
                       K_hat = res$K_hat, scores = res$q_at_step,
                       elapsed = as.numeric(Sys.time() - t0, units = "secs"))
  }

  if ("SuSiE" %in% methods) {
    t0 <- Sys.time()
    fit <- tryCatch(susieR::susie(d$X, d$y, L = K_max, verbose = FALSE),
                      error = function(e) NULL)
    if (!is.null(fit)) {
      sel <- which(fit$pip > theta)
      out$SuSiE <- list(method = "SuSiE", indices = sel,
                          K_hat = length(sel), scores = fit$pip[sel],
                          elapsed = as.numeric(Sys.time() - t0, units = "secs"))
    } else {
      out$SuSiE <- list(method = "SuSiE", indices = integer(0),
                          K_hat = 0L, scores = numeric(0),
                          elapsed = as.numeric(Sys.time() - t0, units = "secs"))
    }
  }


  if ("BSLMM" %in% methods) {
    t0 <- Sys.time()
    if (requireNamespace("BGLR", quietly = TRUE)) {
      fit <- tryCatch({
        suppressMessages(suppressWarnings({
          eta <- list(list(X = d$X, model = "BayesB"))
          BGLR::BGLR(y = d$y, ETA = eta, nIter = 2000L, burnIn = 400L,
                     verbose = FALSE)
        }))
      }, error = function(e) NULL)
      if (!is.null(fit) && !is.null(fit$ETA[[1]]$d)) {
        pip <- fit$ETA[[1]]$d
        sel <- which(pip > 0.5)
        out$BSLMM <- list(method = "BSLMM", indices = sel,
                            K_hat = length(sel), scores = pip[sel],
                            elapsed = as.numeric(Sys.time() - t0, units = "secs"))
      } else {
        out$BSLMM <- list(method = "BSLMM", indices = integer(0),
                            K_hat = 0L, scores = numeric(0),
                            elapsed = as.numeric(Sys.time() - t0, units = "secs"))
      }
    } else {
      out$BSLMM <- list(method = "BSLMM", indices = integer(0), K_hat = 0L,
                          scores = numeric(0), elapsed = NA)
    }
  }

  if ("BayesR" %in% methods) {
    t0 <- Sys.time()
    if (requireNamespace("hibayes", quietly = TRUE)) {
      n_obs <- nrow(d$X)
      Xm <- d$X
      rownames(Xm) <- paste0("S", seq_len(n_obs))
      df_y <- data.frame(id = rownames(Xm), y = d$y, stringsAsFactors = FALSE)
      fit <- tryCatch(suppressMessages(suppressWarnings(
                hibayes::ibrm(formula = y ~ 1, data = df_y,
                               M = Xm, M.id = df_y$id,
                               method = "BayesR",
                               Pi = c(0.95, 0.02, 0.02, 0.01),
                               fold = c(0, 1e-4, 1e-3, 1e-2),
                               niter = 2000L, nburn = 400L,
                               threads = 1L, verbose = FALSE))),
              error = function(e) NULL)
      if (!is.null(fit) && !is.null(fit$pip)) {
        pip <- as.numeric(fit$pip)
        sel <- which(pip > theta)
        out$BayesR <- list(method = "BayesR", indices = sel,
                            K_hat = length(sel), scores = pip[sel],
                            elapsed = as.numeric(Sys.time() - t0, units = "secs"))
      } else {
        out$BayesR <- list(method = "BayesR", indices = integer(0),
                            K_hat = 0L, scores = numeric(0),
                            elapsed = as.numeric(Sys.time() - t0, units = "secs"))
      }
    } else {
      out$BayesR <- list(method = "BayesR", indices = integer(0), K_hat = 0L,
                          scores = numeric(0), elapsed = NA)
    }
  }

  if ("BL" %in% methods) {
    # Bayesian Lasso via BGLR::BGLR (model = "BL"); MCMC with double-exponential
    # prior on the regression coefficients.  Posterior mean |beta_hat| is used
    # as a continuous score (no binary inclusion probability available, unlike
    # spike-and-slab methods).  For threshold-based selection we pick the top-K
    # largest |beta_hat| where K is the design-known K_true (see prose).
    t0 <- Sys.time()
    if (requireNamespace("BGLR", quietly = TRUE)) {
      fit <- tryCatch({
        suppressMessages(suppressWarnings({
          eta_BL <- list(list(X = d$X, model = "BL"))
          BGLR::BGLR(y = d$y, ETA = eta_BL, nIter = 2000L, burnIn = 400L,
                     verbose = FALSE)
        }))
      }, error = function(e) NULL)
      if (!is.null(fit) && !is.null(fit$ETA[[1]]$b)) {
        beta_abs <- abs(as.numeric(fit$ETA[[1]]$b))
        # Top-K_true selection (K_true assumed known via d$truth length)
        K_sel <- length(d$truth)
        sel <- order(-beta_abs)[seq_len(K_sel)]
        out$BL <- list(method = "BL", indices = sel,
                       K_hat = K_sel, scores = beta_abs[sel],
                       elapsed = as.numeric(Sys.time() - t0, units = "secs"))
      } else {
        out$BL <- list(method = "BL", indices = integer(0),
                       K_hat = 0L, scores = numeric(0),
                       elapsed = as.numeric(Sys.time() - t0, units = "secs"))
      }
    } else {
      out$BL <- list(method = "BL", indices = integer(0), K_hat = 0L,
                     scores = numeric(0), elapsed = NA)
    }
  }

  if ("REML_stepwise" %in% methods) {
    t0 <- Sys.time()
    res <- tryCatch(reml_stepwise_score(d$y, d$X, K_max = K_max),
                      error = function(e) list(indices = integer(0),
                                                K_hat = 0L, scores = numeric(0)))
    out$REML_stepwise <- list(method = "REML_stepwise",
                                indices = res$indices, K_hat = res$K_hat,
                                scores = res$scores,
                                elapsed = as.numeric(Sys.time() - t0,
                                                       units = "secs"))
  }

  if ("fastlmm" %in% methods) {
    t0 <- Sys.time()
    res <- tryCatch(fastlmm_score(d$y, d$X, alpha = 0.05),
                      error = function(e) list(indices = integer(0),
                                                K_hat = 0L, scores = numeric(0)))
    out$fastlmm <- list(method = "fastlmm",
                          indices = res$indices, K_hat = res$K_hat,
                          scores = res$scores,
                          elapsed = as.numeric(Sys.time() - t0,
                                                 units = "secs"))
  }

  out
}


# ------------------------------------------------------------------------------
# REML stepwise score (LMM-aware frequentist stepwise baseline)
# Uses rrBLUP::mixed.solve for delta_hat then BIC stopping with |Q|^2.
# ------------------------------------------------------------------------------
reml_stepwise_score <- function(y, X, K_max = 10L) {
  if (!requireNamespace("rrBLUP", quietly = TRUE)) {
    return(list(indices = integer(0), K_hat = 0L, scores = numeric(0)))
  }
  n <- nrow(X); m <- ncol(X)
  selected <- integer(0); scores <- numeric(0)
  remaining <- seq_len(m)
  for (k_step in seq_len(K_max)) {
    rem <- X[, remaining, drop = FALSE]
    K <- tcrossprod(rem) / max(ncol(rem) - 1L, 1L)
    fit_null <- tryCatch(rrBLUP::mixed.solve(y = y, K = K),
                            error = function(e) NULL)
    if (is.null(fit_null)) break
    delta_hat <- fit_null$Vu / fit_null$Ve
    # Whitening
    ee <- eigen(K, symmetric = TRUE)
    s <- ee$values; U <- ee$vectors
    d_inv_sqrt <- 1 / sqrt(1 + delta_hat * s)
    UtX_rem <- crossprod(U, rem)
    Uty     <- crossprod(U, y)
    Z_scores <- as.numeric(crossprod(UtX_rem * d_inv_sqrt,
                                       Uty * d_inv_sqrt) /
                              sqrt(colSums((UtX_rem * d_inv_sqrt)^2) + 1e-300))
    j_best <- which.max(abs(Z_scores))
    if (Z_scores[j_best]^2 < log(n)) break  # BIC-style stop
    selected <- c(selected, remaining[j_best])
    scores   <- c(scores, abs(Z_scores[j_best]))
    remaining <- remaining[-j_best]
  }
  list(indices = selected, K_hat = length(selected), scores = scores)
}


# ------------------------------------------------------------------------------
# FaST-LMM-equivalent: single-pass LMM score test with Bonferroni threshold.
# This is the standard frequentist LMM-GWAS approach (Lippert et al. 2011):
#   1. Fit null LMM y = mu + u + e with u ~ N(0, sigma_g^2 K) on full X-derived GRM
#   2. For each SNP j, compute score statistic under H_0: beta_j = 0
#   3. Select SNPs with |Z| > Bonferroni critical value at alpha = 0.05/m
# Avoids the stepwise greedy selection bias of REML_stepwise.
# ------------------------------------------------------------------------------
fastlmm_score <- function(y, X, alpha = 0.05) {
  if (!requireNamespace("rrBLUP", quietly = TRUE)) {
    return(list(indices = integer(0), K_hat = 0L, scores = numeric(0)))
  }
  n <- nrow(X); m <- ncol(X)
  # Single REML fit on full X-derived GRM
  K <- tcrossprod(X) / max(m - 1L, 1L)
  fit <- tryCatch(rrBLUP::mixed.solve(y = y, K = K),
                    error = function(e) NULL)
  if (is.null(fit)) return(list(indices = integer(0), K_hat = 0L,
                                    scores = numeric(0)))
  delta_hat <- fit$Vu / fit$Ve

  # Whitening on full V_delta = I + delta * K
  ee <- eigen(K, symmetric = TRUE)
  d_inv_sqrt <- 1 / sqrt(1 + delta_hat * ee$values)
  UtX <- crossprod(ee$vectors, X)
  Uty <- crossprod(ee$vectors, y)
  Z   <- as.numeric(crossprod(UtX * d_inv_sqrt, Uty * d_inv_sqrt) /
                       sqrt(colSums((UtX * d_inv_sqrt)^2) + 1e-300))

  # Bonferroni-corrected critical value
  z_crit <- qnorm(1 - alpha / (2 * m))
  sel <- which(abs(Z) > z_crit)
  ord <- order(-abs(Z)[sel])
  list(indices = sel[ord], K_hat = length(sel),
        scores  = abs(Z)[sel[ord]],
        delta_hat = delta_hat, z_crit = z_crit)
}


# ------------------------------------------------------------------------------
# Metrics
# ------------------------------------------------------------------------------
compute_metrics <- function(method_results, truth) {
  K_true <- length(truth)
  rows <- list()
  for (m_name in names(method_results)) {
    r <- method_results[[m_name]]
    tp <- length(intersect(r$indices, truth))
    fp <- r$K_hat - tp
    rec  <- tp / K_true
    prec <- if (r$K_hat > 0) tp / r$K_hat else 0
    f1   <- if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0
    rows[[length(rows) + 1L]] <- data.frame(
      method = m_name, K_hat = r$K_hat, tp = tp, fp = fp,
      recall = rec, precision = prec, f1 = f1,
      elapsed = r$elapsed, stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}


# ------------------------------------------------------------------------------
# CLI argument parser (lightweight, no optparse dependency)
#   Supports: --key value  and  --flag (boolean, presence implies TRUE)
# ------------------------------------------------------------------------------
parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    a <- args[i]
    if (startsWith(a, "--")) {
      key <- sub("^--", "", a)
      if (i + 1L <= length(args) && !startsWith(args[i + 1L], "--")) {
        v <- args[i + 1L]
        # Try numeric coercion
        nv <- suppressWarnings(as.numeric(v))
        out[[key]] <- if (!is.na(nv)) nv else v
        i <- i + 2L
      } else {
        out[[key]] <- TRUE
        i <- i + 1L
      }
    } else {
      i <- i + 1L
    }
  }
  out
}


# Default method set: everything except BSLMM (slow, MCMC).
# varbvs is no longer included. BayesR (hibayes) is a first-class comparator.
default_methods <- function(include_bslmm = FALSE) {
  base <- c("JS_L", "MS_L", "SuSiE", "REML_stepwise", "fastlmm", "BayesR")
  if (include_bslmm) c(base, "BSLMM") else base
}


# Cell ID hash for filenames: short, unambiguous, sortable.
cell_id <- function(...) {
  args <- list(...)
  paste(unlist(lapply(seq_along(args), function(i)
    sprintf("%s%s", names(args)[i], format(args[[i]], scientific = FALSE)))),
    collapse = "_")
}


# Save one (cell, rep) to results/bench_full/<script>/<cell_tag>_b<rep>.rds
save_cell_rep <- function(out_dir, cell_tag, rep_id, payload) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  fn <- file.path(out_dir, sprintf("%s_b%02d.rds", cell_tag, rep_id))
  saveRDS(payload, fn)
  invisible(fn)
}


# Skip if cell+rep already computed (checkpointing).
already_done <- function(out_dir, cell_tag, rep_id) {
  fn <- file.path(out_dir, sprintf("%s_b%02d.rds", cell_tag, rep_id))
  file.exists(fn) && file.size(fn) > 0
}
