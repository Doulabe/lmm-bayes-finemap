# ==============================================================================
#  _remlbf_utils.R
#  Shared utility: plug-in REML stepwise (REML-BF) selection with eBIC stopping,
#  plus a metrics helper.  Sourced by 0Xb_*_remlbf.R, 15b, and 10c.
#
#  This file should NOT be invoked directly with Rscript; it defines functions
#  and exits.  No top-level side effects.
# ==============================================================================

# Defer rrBLUP check to caller; required for REML fits below.

# --- Plug-in REML stepwise selection with eBIC stopping ---------------------
# Generic conditioning matrix W (default intercept only); supports sex /
# batch / etc.  Mirrors the algebra of plug_in_reml_eBIC() in 16_ablation.R
# but with arbitrary W.
plug_in_reml_eBIC_W <- function(y, X, W = NULL, K_max = 10L, tau2 = 0.04) {
  n <- nrow(X); m <- ncol(X)
  if (is.null(W)) W <- matrix(1, n, 1L)
  if (!is.matrix(W)) W <- matrix(W, ncol = 1L)
  kW <- ncol(W)

  selected <- integer(0); scores <- numeric(0); log_BFs <- numeric(0)
  remaining <- seq_len(m)
  ebic_path <- numeric(K_max + 1L); ebic_path[1L] <- Inf

  F_k <- W
  for (k_step in seq_len(K_max)) {
    if (length(remaining) == 0L) break
    rem <- X[, remaining, drop = FALSE]
    K   <- tcrossprod(rem) / max(ncol(rem) - 1L, 1L)
    proj <- solve(crossprod(F_k), crossprod(F_k, y))
    y_r  <- y - F_k %*% proj
    fit  <- tryCatch(rrBLUP::mixed.solve(y = y_r, K = K),
                     error = function(e) NULL)
    if (is.null(fit)) break
    delta_hat <- max(fit$Vu / fit$Ve, 1e-6)
    ee <- eigen(K, symmetric = TRUE)
    s <- pmax(ee$values, 0); U <- ee$vectors
    w_white <- 1 / sqrt(1 + delta_hat * s)
    Uty   <- crossprod(U, y) * w_white
    UtFk  <- crossprod(U, F_k) * w_white
    UtX   <- crossprod(U, rem) * w_white
    A    <- crossprod(UtFk)
    A_inv<- solve(A)
    beta <- A_inv %*% crossprod(UtFk, Uty)
    ty   <- Uty - UtFk %*% beta
    proj_X <- UtFk %*% (A_inv %*% crossprod(UtFk, UtX))
    tXj    <- UtX - proj_X

    D_jk   <- colSums(tXj^2)
    u_j    <- as.numeric(crossprod(tXj, ty))
    rss0   <- as.numeric(crossprod(ty))
    v_j    <- tau2 * D_jk
    df     <- n - ncol(F_k)
    Q2_j   <- u_j^2 / (D_jk * rss0 / df)
    log_BF <- -0.5 * log1p(v_j) - 0.5 * df *
              log1p(-v_j * Q2_j / df / (1 + v_j))
    j_best <- which.max(log_BF)
    if (!is.finite(log_BF[j_best])) break

    sel_new <- c(selected, remaining[j_best])
    F_new   <- cbind(F_k, X[, remaining[j_best], drop = FALSE])
    proj2   <- solve(crossprod(F_new), crossprod(F_new, y))
    rss1    <- as.numeric(crossprod(y - F_new %*% proj2))
    ll_new  <- -0.5 * (n - ncol(F_new)) * log(rss1)
    ebic_new <- -2 * ll_new + length(sel_new) * (log(n) + 2 * log(m))
    if (k_step >= 2L && ebic_new > ebic_path[k_step]) break

    selected <- sel_new
    scores   <- c(scores, log_BF[j_best])
    log_BFs  <- c(log_BFs, log_BF[j_best])
    F_k      <- F_new
    remaining <- remaining[-j_best]
    ebic_path[k_step + 1L] <- ebic_new
  }
  list(method = "REML_BF_eBIC",
       indices = selected, K_hat = length(selected),
       scores = scores, log_BFs = log_BFs)
}

# --- Per-replicate REML-BF on a simulated dataset --------------------------
# Returns a metrics data frame compatible with the MBF benchmark (same
# columns as compute_metrics() in 00_config.R: method, K_hat, tp, fp,
# recall, precision, f1, elapsed).  The "method" label is REML_BF_eBIC.
remlbf_metrics_one_rep <- function(d, K_max = 10L, tau2 = 0.04) {
  truth <- d$truth
  t0  <- Sys.time()
  fit <- plug_in_reml_eBIC_W(d$y, d$X, K_max = K_max, tau2 = tau2)
  elapsed <- as.numeric(Sys.time() - t0, units = "secs")
  tp  <- length(intersect(fit$indices, truth))
  fp  <- length(setdiff(fit$indices, truth))
  recall    <- if (length(truth) > 0) tp / length(truth) else NA
  precision <- if (length(fit$indices) > 0) tp / length(fit$indices) else NA
  f1 <- if (is.na(recall) || is.na(precision) || (recall + precision) == 0) 0
        else 2 * recall * precision / (recall + precision)
  data.frame(method = "REML_BF_eBIC",
             K_hat = fit$K_hat, tp = tp, fp = fp,
             recall = recall, precision = precision, f1 = f1,
             elapsed = elapsed,
             stringsAsFactors = FALSE)
}

# --- Generic axis runner ----------------------------------------------------
# axis_cells: list of lists each with (tag, n, m, rho, K_true, beta_true,
#                                     sigma_g2, block_size, seed_offset)
# Saves per-(cell, rep) RDS in out_dir.
run_axis_remlbf <- function(axis_cells, B, out_dir, n_cores,
                             K_max = 10L, tau2 = 0.04, anchor_seed = 20260425L) {
  suppressPackageStartupMessages(library(parallel))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  jobs <- list()
  for (c in axis_cells) for (b in seq_len(B)) jobs[[length(jobs)+1L]] <- list(c, b)
  cat(sprintf(" %d cells x B=%d = %d jobs (cores=%d)\n",
              length(axis_cells), B, length(jobs), n_cores))
  t0 <- Sys.time()
  mclapply(jobs, function(j) {
    cell <- j[[1]]; rep_id <- j[[2]]
    out_fn <- file.path(out_dir, sprintf("%s_b%03d.rds", cell$tag, rep_id))
    if (file.exists(out_fn)) return(invisible(NULL))
    seed <- anchor_seed + 1000L * (rep_id - 1L) +
            as.integer(cell$n) + round(1000 * cell$sigma_g2) +
            cell$seed_offset
    d <- tryCatch(gen_dataset(n = cell$n, m = cell$m, rho = cell$rho,
                              K_true = cell$K_true, beta_true = cell$beta_true,
                              sigma_g2 = cell$sigma_g2,
                              block_size = cell$block_size, seed = seed),
                  error = function(e) NULL)
    if (is.null(d)) return(invisible(NULL))
    metrics <- tryCatch(remlbf_metrics_one_rep(d, K_max = K_max, tau2 = tau2),
                        error = function(e) {
                          message(sprintf("Err %s b%d: %s", cell$tag, rep_id,
                                          conditionMessage(e)))
                          NULL
                        })
    if (is.null(metrics)) return(invisible(NULL))
    metrics$n <- cell$n; metrics$m <- cell$m; metrics$rho <- cell$rho
    metrics$sigma_g2 <- cell$sigma_g2; metrics$rep <- rep_id
    payload <- list(metrics = metrics, truth = d$truth)
    saveRDS(payload, out_fn)
    invisible(NULL)
  }, mc.cores = n_cores)
  cat(sprintf(" done in %.1f min\n",
              as.numeric(Sys.time() - t0, units = "mins")))
}
