# ==============================================================================
# CBF_LMM_exact.R
#
# Exact candidate-specific implementation of CBF-LMM.
#
# Statistical model (per candidate j at step k):
#   y = F_k alpha + x_j beta + u_jk + eps,   u_jk ~ N(0, sg2 * K_jk),
#   K_jk = X_{B_jk} X_{B_jk}' / m_k,   B_jk = C_k \ {j},   m_k = m - k.
# The candidate is EXCLUDED from its own polygenic background.
#
# Computation: one eigendecomposition of the raw pool cross-product
# X_{C_k} X_{C_k}' per step. With Rt = X_{C_k} X_{C_k}' / m_k,
#   Z_jk(delta) = I + delta * K_jk = (I + delta * Rt) - a x_j x_j',
#   a = delta / m_k,
# so every candidate-specific covariance is a rank-one downdate of the SAME
# matrix A = I + delta * Rt, diagonal in the pool eigenbasis. Sherman-Morrison
# and the matrix determinant lemma then give all GLS quantities
#   D_jk   = x_j' P_jk x_j,   u_jk = x_j' P_jk y,   RSS0_jk = y' P_jk y,
#   P_jk   = Z^-1 - Z^-1 F (F' Z^-1 F)^-1 F' Z^-1,
# in closed vectorised form over all candidates — the whitening Z^{-1/2} is
# never formed. Cost per step ~ the shared-kernel implementation.
#
# Stopping: TRUE extended BIC on the profile MAXIMUM likelihood of the full
# accepted model (Chen & Chen 2008, sparse approximation):
#   after proposing S = S_{k-1} U {l_k}, fit
#     y = 1 alpha + X_S beta_S + u_S + eps,  u_S ~ N(0, sg2 * K(S)),
#     K(S) = X_{S^c} X_{S^c}' / (m - |S|),
#   maximise the profile ML log-likelihood over delta, and accept iff
#     eBIC(S) = -2 l_ML(S) + |S| {log n + 2 gamma log m}
#   improves on eBIC(S_{k-1}). The FIRST candidate must also pass the rule
#   (the procedure may return the empty set under the global null).
#   The eigendecomposition needed for K(S) is the next step's pool eigen,
#   obtained by a rank-one downdate of the current raw cross-product.
#
# delta handling for the conditional BF:
#   "marginal": numerical marginalisation, candidate-specific null
#               pi_jk(delta | y, H0_jk) handled implicitly by integrating
#               numerator and denominator under the same prior grid;
#   "reml"    : per-candidate grid-profile delta_hat_jk maximising the
#               H0_jk restricted likelihood over the quadrature grid.
# ==============================================================================

.try_source_exact <- function(paths) {
  for (p in paths) if (file.exists(p)) { source(p); return(invisible(p)) }
}
.try_source_exact(c("R/LMM_core.R", "../R/LMM_core.R", "LMM_core.R"))

# ──────────────────────────────────────────────────────────────────────────────
# 1. Exact per-candidate H0/H1 restricted log-likelihood pieces at fixed delta
#    All inputs are in the POOL eigenbasis (U, s_raw of X_C X_C').
# ──────────────────────────────────────────────────────────────────────────────

#' @param Uty,UtF,UtX  eigen-coordinates of y (n), F_k (n x k), X_cand (n x m_c)
#' @param s_rank       eigenvalues of Rt = X_C X_C' / m_k
#' @param a            delta / m_k  (rank-one downdate weight)
#' @param delta        variance-ratio value
#' @return list of m_c-vectors lp0, lp1 (restricted log-lik pieces, constants
#'         dropped) plus D, u, RSS0 for reuse.
exact_pieces_at_delta <- function(Uty, UtF, UtX, s_rank, delta, a,
                                    tau2, n, k) {
  d_inv <- 1 / (1 + delta * s_rank)                    # n
  WU_y  <- d_inv * Uty
  q0    <- sum(Uty * WU_y)                              # y'A^-1 y
  p0    <- as.numeric(crossprod(UtF, WU_y))             # k
  M0    <- crossprod(UtF, d_inv * UtF)                  # k x k
  L0    <- chol(M0 + 1e-12 * diag(k))
  w0    <- backsolve(L0, forwardsolve(t(L0), p0))       # M0^-1 p0
  S00   <- sum(p0 * w0)
  logdetM0 <- 2 * sum(log(diag(L0)))
  logdetA  <- sum(log1p(delta * s_rank))

  alpha <- colSums(d_inv * UtX^2)                       # m_c : x'A^-1 x
  beta  <- as.numeric(crossprod(UtX, WU_y))             # m_c : x'A^-1 y
  FWX   <- crossprod(UtF, d_inv * UtX)                  # k x m_c
  H     <- backsolve(L0, forwardsolve(t(L0), FWX))      # k x m_c : M0^-1 f_j
  phi   <- colSums(FWX * H)                             # m_c : f'M0^-1 f
  psi   <- as.numeric(crossprod(FWX, w0))               # m_c : f'M0^-1 p0

  t_j   <- pmax(1 - a * alpha, 1e-12)                   # det lemma factor
  g_j   <- a / t_j
  denom <- 1 + g_j * phi

  # y' Z^-1 y and S0 = (F'Z^-1 y)'(F'Z^-1 F)^-1(F'Z^-1 y)
  yZy   <- q0 + g_j * beta^2
  pMp   <- S00 + 2 * g_j * beta * psi + g_j^2 * beta^2 * phi
  hp    <- psi + g_j * beta * phi
  S0    <- pMp - g_j * hp^2 / denom
  RSS0  <- pmax(yZy - S0, 1e-300)

  D_j   <- pmax(alpha / t_j - phi / (denom * t_j^2), 1e-300)
  u_j   <- (beta - hp / denom) / t_j

  logdetZ <- logdetA + log(t_j)
  logdetM <- logdetM0 + log(denom)

  v     <- tau2 * D_j
  RSS1  <- pmax(RSS0 - v * u_j^2 / (D_j * (1 + v)), 1e-300)

  base  <- -0.5 * logdetZ - 0.5 * logdetM
  lp0   <- base - 0.5 * (n - k) * log(RSS0)
  lp1   <- base - 0.5 * (n - k) * log(RSS1) - 0.5 * log1p(v)

  list(lp0 = lp0, lp1 = lp1, D = D_j, u = u_j, RSS0 = RSS0)
}

# ──────────────────────────────────────────────────────────────────────────────
# 2. Exact conditional BF for all candidates at a step
# ──────────────────────────────────────────────────────────────────────────────

exact_bf_step <- function(y, X_cand, F_k, U, s_raw, m_k,
                            tau2 = 0.04,
                            delta_eval = c("marginal", "reml"),
                            prior = "half_cauchy",
                            prior_params = list(scale = 1.0),
                            n_nodes = 15L,
                            delta_lo = 1e-4, delta_hi = 1e3) {
  delta_eval <- match.arg(delta_eval)
  n <- length(y); m_c <- ncol(X_cand); k <- ncol(F_k)
  s_rank <- pmax(s_raw, 0) / m_k

  Uty <- as.numeric(crossprod(U, y))
  UtF <- crossprod(U, F_k)
  UtX <- crossprod(U, X_cand)

  grid    <- make_delta_grid(n_nodes, delta_lo, delta_hi)
  dv      <- grid$delta
  prior_d <- lmm_prior_delta(dv, prior = prior, params = prior_params)
  log_q_w <- log(grid$w_GL) + log(dv) + log(pmax(prior_d, 1e-300))
  log_q_w <- log_q_w - max(log_q_w)

  lp0_mat <- matrix(0, m_c, n_nodes)
  lp1_mat <- matrix(0, m_c, n_nodes)
  q2_mat  <- matrix(0, m_c, n_nodes)
  for (r in seq_len(n_nodes)) {
    pc <- exact_pieces_at_delta(Uty, UtF, UtX, s_rank, dv[r],
                                  a = dv[r] / m_k, tau2 = tau2, n = n, k = k)
    lp0_mat[, r] <- pc$lp0
    lp1_mat[, r] <- pc$lp1
    q2_mat[, r]  <- pc$u^2 / (pc$D * pc$RSS0 / (n - k))
  }

  # per-candidate grid-profile delta_hat_jk under H0_jk (REML plug-in node)
  idx <- max.col(lp0_mat, ties.method = "first")
  Q2_reml <- q2_mat[cbind(seq_len(m_c), idx)]

  if (delta_eval == "marginal") {
    lse <- function(v_) { M <- max(v_); M + log(sum(exp(v_ - M))) }
    log_BF <- vapply(seq_len(m_c), function(j)
      lse(lp1_mat[j, ] + log_q_w) - lse(lp0_mat[j, ] + log_q_w), numeric(1L))
  } else {
    log_BF <- lp1_mat[cbind(seq_len(m_c), idx)] -
              lp0_mat[cbind(seq_len(m_c), idx)]
  }
  list(log_BF = log_BF, Q2_reml = Q2_reml,
       lp0 = lp0_mat, lp1 = lp1_mat, grid = dv)
}

# ──────────────────────────────────────────────────────────────────────────────
# 3. Profile-ML eBIC of the full accepted model
# ──────────────────────────────────────────────────────────────────────────────

#' -2 * profile ML log-likelihood (constants dropped) of
#'   y = F alpha + u + eps,  u ~ N(0, sg2 * K),  K eigen = (U_K, s_K)
#' maximised over delta on the log scale.
neg2_profile_ml <- function(y, Fmat, U_K, s_K,
                              log_delta_lo = log(1e-6),
                              log_delta_hi = log(1e6)) {
  n <- length(y); p <- ncol(Fmat)
  Uty <- as.numeric(crossprod(U_K, y))
  UtF <- crossprod(U_K, Fmat)
  s_K <- pmax(s_K, 0)
  f <- function(log_delta) {
    delta <- exp(log_delta)
    d_inv <- 1 / (1 + delta * s_K)
    M     <- crossprod(UtF, d_inv * UtF)
    pv    <- as.numeric(crossprod(UtF, d_inv * Uty))
    ch    <- tryCatch(chol(M + 1e-12 * diag(p)), error = function(e) NULL)
    if (is.null(ch)) return(1e10)
    sol   <- backsolve(ch, forwardsolve(t(ch), pv))
    rss   <- sum(Uty^2 * d_inv) - sum(pv * sol)
    if (!is.finite(rss) || rss <= 0) return(1e10)
    sum(log1p(delta * s_K)) + n * log(rss / n)
  }
  opt <- optimize(f, c(log_delta_lo, log_delta_hi))
  min(opt$objective, f(log_delta_lo))   # include the near-zero boundary
}

# ──────────────────────────────────────────────────────────────────────────────
# 4. Full exact stepwise procedure
# ──────────────────────────────────────────────────────────────────────────────

CBF_LMM_stepwise_exact <- function(y, X, a_n = NULL,
                                     tau2 = 0.04, K_max = 10L,
                                     delta_eval = c("marginal", "reml"),
                                     gamma = 1,
                                     n_nodes = 15L,
                                     prior = "half_cauchy",
                                     prior_params = list(scale = 1.0),
                                     delta_lo = 1e-4, delta_hi = 1e3,
                                     force_K = NULL,
                                     rank_by = c("bf", "score")) {
  # rank_by = "score": Score-LMM ablation — identical K_jk, delta handling and
  # profile-ML eBIC stopping, but candidates ranked by the frequentist
  # Q_jk^2 at the per-candidate REML plug-in node instead of the Bayes factor.
  rank_by <- match.arg(rank_by)
  # force_K: run exactly force_K greedy steps with NO stopping test (used for
  # threshold-independent top-K ranking); the eBIC machinery is skipped.
  delta_eval <- match.arg(delta_eval)
  n <- nrow(X); m <- ncol(X)
  if (is.null(a_n)) a_n <- rep(1, n)
  if (!is.matrix(a_n)) a_n <- matrix(a_n, ncol = 1L)
  if (!is.null(force_K)) K_max <- as.integer(force_K)
  pen_unit <- log(n) + 2 * gamma * log(m)

  selected  <- integer(0)
  remaining <- seq_len(m)
  log_BF_at_step <- numeric(0)
  ebic_path <- numeric(0)

  # step-1 pool: raw cross-product of all markers (also the S_0 eBIC kernel)
  raw  <- tcrossprod(X)
  eig  <- eigen(raw, symmetric = TRUE)
  U    <- eig$vectors; s_raw <- pmax(eig$values, 0)

  # eBIC(S_0): intercept-only fixed part, kernel = all markers / m
  ebic_prev <- if (is.null(force_K)) neg2_profile_ml(y, a_n, U, s_raw / m)
               else Inf
  ebic_path <- c(ebic_path, ebic_prev)

  for (k_step in seq_len(K_max)) {
    if (length(remaining) <= 1L) break
    m_k <- m - k_step
    F_k <- if (length(selected) == 0L) a_n
           else cbind(a_n, X[, selected, drop = FALSE])

    bf <- exact_bf_step(y, X[, remaining, drop = FALSE], F_k,
                          U, s_raw, m_k,
                          tau2 = tau2, delta_eval = delta_eval,
                          prior = prior, prior_params = prior_params,
                          n_nodes = n_nodes,
                          delta_lo = delta_lo, delta_hi = delta_hi)
    rank_score <- if (rank_by == "bf") bf$log_BF else bf$Q2_reml
    best_idx <- which.max(rank_score)
    best_j   <- remaining[best_idx]
    if (!is.finite(rank_score[best_idx])) break

    # next pool raw cross-product: rank-one downdate, then eigen (needed for
    # the eBIC kernel of the proposal, and for the next step if accepted)
    xj       <- X[, best_j]
    raw_next <- raw - tcrossprod(xj)
    eig_next <- eigen(raw_next, symmetric = TRUE)
    U_next   <- eig_next$vectors
    s_next   <- pmax(eig_next$values, 0)

    S_prop  <- c(selected, best_j)
    if (is.null(force_K)) {
      F_prop  <- cbind(a_n, X[, S_prop, drop = FALSE])
      ebic_new <- neg2_profile_ml(y, F_prop, U_next,
                                    s_next / (m - length(S_prop))) +
                  length(S_prop) * pen_unit
      # every candidate, including the first, must pass the stopping rule
      if (ebic_new >= ebic_prev) break
    } else {
      ebic_new <- NA_real_
    }

    selected       <- S_prop
    log_BF_at_step <- c(log_BF_at_step, bf$log_BF[best_idx])
    ebic_path      <- c(ebic_path, ebic_new)
    ebic_prev      <- ebic_new
    remaining      <- setdiff(remaining, best_j)
    raw <- raw_next; U <- U_next; s_raw <- s_next
  }

  list(indices = selected,
       K_hat = length(selected),
       log_BF_at_step = log_BF_at_step,
       ebic_path = ebic_path,
       variant = paste0("exact_", delta_eval))
}

# ──────────────────────────────────────────────────────────────────────────────
# 5. Brute-force reference (validation only): one eigendecomposition per
#    candidate-specific kernel K_jk, same likelihood pieces, same quadrature.
# ──────────────────────────────────────────────────────────────────────────────

exact_bf_step_bruteforce <- function(y, X_cand, F_k, X_pool, m_k,
                                       tau2 = 0.04,
                                       prior = "half_cauchy",
                                       prior_params = list(scale = 1.0),
                                       n_nodes = 15L,
                                       delta_lo = 1e-4, delta_hi = 1e3) {
  n <- length(y); m_c <- ncol(X_cand); k <- ncol(F_k)
  grid    <- make_delta_grid(n_nodes, delta_lo, delta_hi)
  dv      <- grid$delta
  prior_d <- lmm_prior_delta(dv, prior = prior, params = prior_params)
  log_q_w <- log(grid$w_GL) + log(dv) + log(pmax(prior_d, 1e-300))
  log_q_w <- log_q_w - max(log_q_w)
  lse <- function(v_) { M <- max(v_); M + log(sum(exp(v_ - M))) }

  log_BF <- numeric(m_c)
  for (j in seq_len(m_c)) {
    # K_jk from scratch
    keep <- setdiff(seq_len(ncol(X_pool)), j)
    Kj   <- tcrossprod(X_pool[, keep, drop = FALSE]) / m_k
    ee   <- eigen(Kj, symmetric = TRUE)
    s    <- pmax(ee$values, 0); Uj <- ee$vectors
    Uty  <- as.numeric(crossprod(Uj, y))
    UtF  <- crossprod(Uj, F_k)
    Utx  <- as.numeric(crossprod(Uj, X_cand[, j]))
    lp0v <- numeric(n_nodes); lp1v <- numeric(n_nodes)
    for (r in seq_len(n_nodes)) {
      d_r  <- 1 / sqrt(1 + dv[r] * s)
      yv   <- d_r * Uty; Fv <- d_r * UtF; xv <- d_r * Utx
      FtF  <- crossprod(Fv); ch <- chol(FtF + 1e-12 * diag(k))
      th_y <- backsolve(ch, forwardsolve(t(ch), crossprod(Fv, yv)))
      th_x <- backsolve(ch, forwardsolve(t(ch), crossprod(Fv, xv)))
      ry   <- yv - as.numeric(Fv %*% th_y)
      rx   <- xv - as.numeric(Fv %*% th_x)
      D    <- sum(rx^2); u <- sum(rx * ry)
      rss0 <- sum(ry^2); v <- tau2 * D
      rss1 <- max(rss0 - v * u^2 / (D * (1 + v)), 1e-300)
      base <- -0.5 * sum(log1p(dv[r] * s)) - 0.5 * 2 * sum(log(diag(ch)))
      lp0v[r] <- base - 0.5 * (n - k) * log(rss0)
      lp1v[r] <- base - 0.5 * (n - k) * log(rss1) - 0.5 * log1p(v)
    }
    log_BF[j] <- lse(lp1v + log_q_w) - lse(lp0v + log_q_w)
  }
  list(log_BF = log_BF)
}
