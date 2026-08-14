# ==============================================================================
#  LMM_reml_bf.R
#  Plug-in REML evaluation of the conditional Bayes factor (REML-BF).
#
#  plug_in_reml_eBIC_W() performs stepwise selection scoring each candidate by
#  the closed-form conditional Bayes factor evaluated at the per-step restricted
#  maximum-likelihood estimate of the variance-component ratio delta, with eBIC
#  stopping.  This is the "REML-BF" operating point: the balanced default
#  recommended for routine refinement.  The marginalized (MBF) operating point
#  is provided by {MS,JS}_L_LMM_stepwise_fast(..., delta_eval = "marginal").
#
#  This is a verbatim mirror of plug_in_reml_eBIC_W() in
#  .../_remlbf_utils.R, the function that produced the REML-BF
#  results reported in the manuscript; keep the two in sync.
#
#  Requires rrBLUP (for the per-step REML fit).
# ==============================================================================

plug_in_reml_eBIC_W <- function(y, X, W = NULL, K_max = 10L, tau2 = 0.04) {
  n <- nrow(X); m <- ncol(X)
  if (is.null(W)) W <- matrix(1, n, 1L)
  if (!is.matrix(W)) W <- matrix(W, ncol = 1L)

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
       scores = scores, log_BFs = log_BFs,
       ebic_path = ebic_path[seq_len(length(selected) + 1L)])
}
