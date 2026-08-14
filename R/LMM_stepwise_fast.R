# ==============================================================================
# LMM_stepwise_fast.R
# Vectorized / batched implementation of the marginal Bayes factor
# computation across all m candidates simultaneously.
#
# Strategy:
#   - The eigendecomposition R = U S U^T is shared across all candidates
#     and across all delta-grid nodes.
#   - For each step, precompute U^T X_cand (n x m_cand), U^T F_k (n x k),
#     and U^T y (n) ONCE at O(n^2 m_cand). This is the dominant cost.
#   - For each delta in the grid, the action of Z_delta^{-1/2} on a column
#     reduces to scaling the eigen-coordinate vector by 1/sqrt(1 + delta * s).
#   - The projection P_F^* = I - F^*(F^*'F^*)^{-1} F^*' factors are computed
#     analytically in eigen coordinates without leaving the (n,k) space.
#   - All RSS_0(delta), RSS_1(delta), Q*(delta) for ALL candidates are
#     produced as length-m_cand vectors via column-wise arithmetic.
#
#
# Returns log-BF marginalized for each candidate.
# ==============================================================================

# Source dependencies (must have LMM_core.R loaded)
.try_source <- function(paths) {
  for (p in paths) if (file.exists(p)) { source(p); return(invisible(p)) }
}
.try_source(c("R/LMM_core.R", "../R/LMM_core.R", "LMM_core.R"))


# ------------------------------------------------------------------------------
# 1. Batched RSS quantities at fixed delta for all candidates simultaneously
# ------------------------------------------------------------------------------

#' Vectorised RSS_0(delta) and RSS_1(delta) for all m candidates simultaneously.
#'
#' Inputs are precomputed eigen-coordinates:
#'   UtX_cand : n x m_cand  (=  U^T X_cand,    one-time precompute)
#'   UtF_k    : n x k       (=  U^T F_k)
#'   Uty      : n           (=  U^T y)
#'   s        : eigenvalues of R
#' The whitening factor at delta is d_r = 1/sqrt(1 + delta * s).
#'
#' Returns vectors RSS_0, RSS_1, Q_star, D_jk, ... all of length m_cand.
rss_at_delta_k_batched <- function(UtX_cand, UtF_k, Uty,
                                     s, delta, tau2 = 1.0,
                                     n, k) {
  d_r <- 1 / sqrt(1 + delta * s)               # n-vector

  # Whitened eigen-coordinates (all in eigen frame, no need to map back to X-frame)
  Uty_s   <- d_r * Uty                          # n
  UtFk_s  <- d_r * UtF_k                        # n x k     (broadcasting)
  UtXc_s  <- d_r * UtX_cand                     # n x m_cand

  # In eigen coordinates the inner products are preserved (unitary):
  #   (F_k^*)' (F_k^*) = (U^T F_k_s)' (U^T F_k_s) = (UtF_k_s)' UtF_k_s
  # because U is orthonormal and we kept the same basis.
  FtF      <- crossprod(UtFk_s)                 # k x k
  L_FtF    <- chol(FtF + 1e-10 * diag(k))
  Fty_s    <- crossprod(UtFk_s, Uty_s)          # k
  FtXc     <- crossprod(UtFk_s, UtXc_s)         # k x m_cand

  theta_y  <- backsolve(L_FtF, forwardsolve(t(L_FtF), Fty_s))         # k
  theta_X  <- backsolve(L_FtF, forwardsolve(t(L_FtF), FtXc))          # k x m_cand

  # Projected residual r_perp_eigen = Uty_s - UtFk_s %*% theta_y    (n)
  res_s    <- Uty_s - as.numeric(UtFk_s %*% theta_y)                  # n
  # Projected candidate Xj_perp_eigen = UtXc_s - UtFk_s %*% theta_X  (n x m_cand)
  Xc_perp  <- UtXc_s - UtFk_s %*% theta_X                             # n x m_cand

  # All inner products via column-wise sums in eigen frame:
  D_jk     <- colSums(Xc_perp^2)                                       # m_cand
  rss0_all <- as.numeric(crossprod(res_s))                             # scalar (= y^T P y)
  Xjy_p    <- as.numeric(crossprod(Xc_perp, res_s))                    # m_cand

  v_vec    <- tau2 * D_jk
  rss1_all <- rss0_all - (v_vec * Xjy_p^2) / (D_jk * (1 + v_vec))

  sigma2_hat <- rss0_all / (n - k)
  Q_star     <- Xjy_p / sqrt(D_jk * sigma2_hat)
  log_det_FtF <- 2 * sum(log(diag(L_FtF)))
  log_det_Z   <- sum(log1p(delta * s))

  list(rss0 = rep(rss0_all, length(D_jk)),       # broadcast for shape compatibility
       rss1 = rss1_all,
       Q_star = Q_star,
       D_jk   = D_jk,
       v      = v_vec,
       sigma2_hat = sigma2_hat,
       log_det_Z = log_det_Z,
       log_det_FtF = log_det_FtF)
}


# ------------------------------------------------------------------------------
# 2. Batched marginalized Bayes factor for all candidates at step k
# ------------------------------------------------------------------------------

#' Marginalized Bayes factor for all m candidates simultaneously, at step k.
#'
#' @param y          n-vector
#' @param X_cand     n x m_cand candidate matrix
#' @param F_k        n x k nuisance design (intercept + selected SNPs)
#' @param R_eig      list(U, s) eigendecomposition of R
#' @param tau2       slab variance
#' @param prior, prior_params, n_nodes, delta_lo, delta_hi: as in lmm_bf_marginalized_k
#'
#' @return list(log_BF = m_cand-vector, log_BF_at_delta = matrix m_cand x n_nodes)
lmm_bf_marginalised_k_batched <- function(y, X_cand, F_k, R_eig,
                                            tau2 = 1.0,
                                            prior = "half_cauchy",
                                            prior_params = list(scale = 1.0),
                                            n_nodes = 30L,
                                            delta_lo = 1e-4,
                                            delta_hi = 1e3) {
  n <- length(y); m_cand <- ncol(X_cand); k <- ncol(F_k)
  U <- R_eig$U; s <- R_eig$s

  # -- ONE-TIME precomputes (the dominant cost) -------------------------------
  UtX_cand <- crossprod(U, X_cand)                # n x m_cand   
  UtF_k    <- crossprod(U, F_k)                   # n x k       
  Uty      <- as.numeric(crossprod(U, y))         # n             

  # -- Quadrature setup -------------------------------------------------------
  grid    <- make_delta_grid(n_nodes, delta_lo, delta_hi)
  dv      <- grid$delta
  prior_d <- lmm_prior_delta(dv, prior = prior, params = prior_params)
  log_q_w <- log(grid$w_GL) + log(dv) + log(pmax(prior_d, 1e-300))
  log_q_w <- log_q_w - max(log_q_w)

  # -- Sweep over delta-grid; arithmetic stays vectorised over m_cand -------
  log_p_h0_mat <- matrix(0, nrow = m_cand, ncol = n_nodes)
  log_p_h1_mat <- matrix(0, nrow = m_cand, ncol = n_nodes)
  for (r in seq_len(n_nodes)) {
    rss <- rss_at_delta_k_batched(UtX_cand, UtF_k, Uty,
                                    s, dv[r], tau2, n, k)
    base <- -0.5 * rss$log_det_Z - 0.5 * rss$log_det_FtF
    log_p_h0_mat[, r] <- base - 0.5 * (n - k) * log(rss$rss0)
    log_p_h1_mat[, r] <- base - 0.5 * (n - k) * log(rss$rss1) -
                          0.5 * log1p(rss$v)
  }

  # -- log-sum-exp per candidate ----------------------------------------------
  lse_row <- function(row_vec, weights) {
    M <- max(row_vec + weights)
    M + log(sum(exp(row_vec + weights - M)))
  }
  log_BF <- numeric(m_cand)
  for (j in seq_len(m_cand)) {
    log_BF[j] <- lse_row(log_p_h1_mat[j, ], log_q_w) -
                  lse_row(log_p_h0_mat[j, ], log_q_w)
  }

  list(log_BF = log_BF,
       grid = dv,
       log_BF_at_delta = log_p_h1_mat - log_p_h0_mat)
}


# ------------------------------------------------------------------------------
# 3. Stepwise wrapper using the batched marginal BF
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 3a. Batched JOINT BF at fixed delta for all candidates simultaneously
# ------------------------------------------------------------------------------

#' Batched joint Bayes factor (Proposition prop:joint_BF) for all candidates
#' at fixed delta. Same precomputation pattern as the marginal batched form,
#' plus a Schur-complement step done in vectorised form.
#'
#' Inputs (all in eigen frame, precomputed):
#'   UtX_cand : n x m_cand
#'   UtX_S    : n x k_sel    (selected SNPs only, NOT intercept; n x 0 at k=1)
#'   UtF_k    : n x k_full   (intercept + selected)
#'   Uty      : n
#'   s        : eigenvalues
#' EXACT joint Bayes factor at fixed delta, batched over candidates.
#' Replaces the Wakefield asymptotic approximation with the exact closed
#' form via Sherman-Morrison-Woodbury under the Gaussian likelihood with
#' Jeffreys prior on sigma^2 and joint k-dim Gaussian slab on lambda.
#'
#' Formula (per candidate j, k_sel = number of selected SNPs):
#'   A    = X_S^T Z^{-1} X_S    (k_sel x k_sel, shared)
#'   p    = X_S^T Z^{-1} y_c    (k_sel-vector,  shared)
#'   v    = (I + tau^2 A)^{-1} p
#'   b_j  = X_S^T Z^{-1} X_j    (k_sel-vector,  per-candidate)
#'   u_j  = (I + tau^2 A)^{-1} b_j
#'   c_j  = X_j^T Z^{-1} X_j     (scalar,        per-candidate)
#'   q_j  = X_j^T Z^{-1} y_c     (scalar,        per-candidate)
#'   S_j  = 1 + tau^2(c_j - tau^2 b_j^T u_j)
#'   RSS0_joint = y_c^T Z^{-1} y_c - tau^2 p^T v        (scalar, shared)
#'   RSS1_joint_j = RSS0_joint - tau^2 (q_j - tau^2 b_j^T v)^2 / S_j
#'   log_BF_j = -0.5 log(S_j)
#'              + 0.5 (n-1) (log RSS0_joint - log RSS1_joint_j)
#'
#' At k_sel = 0, A is empty and the formula reduces algebraically to
#' the marginal Bayes factor (MS_L) — by construction, since the joint
#' slab on a 1-dim coefficient is the same as the marginal slab.
joint_bf_exact_k_batched <- function(UtX_cand, UtX_S, Uty_c,
                                       s, delta, tau2 = 1.0,
                                       n, k_sel, RSS_0) {
  d_inv <- 1 / (1 + delta * s)               # n-vector

  # Per-candidate primary statistics in V^{-1} metric
  Vinv_Xc <- d_inv * UtX_cand                 # n x m_cand (eigen-frame scaled)
  c_vec   <- colSums(Vinv_Xc * UtX_cand)      # m_cand: X_j^T V^{-1} X_j
  q_vec   <- as.numeric(crossprod(UtX_cand, d_inv * Uty_c))  # m_cand

  if (k_sel > 0L) {
    # Shared-across-candidates blocks
    Vinv_XS  <- d_inv * UtX_S                                 # n x k_sel
    A        <- crossprod(UtX_S, Vinv_XS)                     # k_sel x k_sel
    p_vec    <- as.numeric(crossprod(UtX_S, d_inv * Uty_c))   # k_sel
    M_aug    <- diag(k_sel) + tau2 * A
    L_aug    <- chol(M_aug + 1e-10 * diag(k_sel))
    # v = (I + tau^2 A)^{-1} p
    v_vec    <- backsolve(L_aug,
                            forwardsolve(t(L_aug), p_vec))     # k_sel
    pT_v     <- as.numeric(crossprod(p_vec, v_vec))            # scalar

    # Per-candidate quantities
    b_mat    <- crossprod(UtX_S, Vinv_Xc)                      # k_sel x m_cand
    # u_j = (I + tau^2 A)^{-1} b_j   (column-wise solve)
    u_mat    <- backsolve(L_aug,
                            forwardsolve(t(L_aug), b_mat))      # k_sel x m_cand
    bT_u     <- colSums(b_mat * u_mat)                          # m_cand
    bT_v     <- as.numeric(crossprod(b_mat, v_vec))             # m_cand

    S_vec        <- 1 + tau2 * (c_vec - tau2 * bT_u)            # m_cand
    num_vec      <- q_vec - tau2 * bT_v                         # m_cand
    RSS0_joint   <- RSS_0 - tau2 * pT_v                         # shared
  } else {
    S_vec      <- 1 + tau2 * c_vec
    num_vec    <- q_vec
    RSS0_joint <- RSS_0
  }
  # RSS_1_joint per candidate
  RSS1_vec <- RSS0_joint - tau2 * num_vec^2 / S_vec
  RSS1_vec <- pmax(RSS1_vec, 1e-300)

  log_bf <- -0.5 * log(S_vec) +
              0.5 * (n - 1) * (log(RSS0_joint) - log(RSS1_vec))

  list(log_BF = log_bf,
       S = S_vec, c = c_vec, q = q_vec,
       RSS0_joint = RSS0_joint,
       RSS1 = RSS1_vec,
       log_det_Z = sum(log1p(delta * s)))
}


# ------------------------------------------------------------------------------
# 3a. Batched JOINT BF at fixed delta for all candidates simultaneously
#     (Wakefield asymptotic form — kept for reference / regression testing)
# ------------------------------------------------------------------------------

joint_bf_at_delta_k_batched <- function(UtX_cand, UtX_S, UtF_k, Uty,
                                          s, delta, tau2 = 1.0,
                                          n, k_full, k_sel) {
  # Apply Z_delta^{-1} (NOT invsqrt) in eigen frame: scale by 1/(1+delta*s)
  d_inv <- 1 / (1 + delta * s)                       # n-vector

  # All inner products via the V^{-1} metric: A' V^{-1} B = (UtA)' diag(d_inv) (UtB)
  # We use weighted dot products in eigen frame.

  # d_jj^V = X_j' V^{-1} X_j  (un-scaled)
  d_jj_vec <- colSums((d_inv * UtX_cand) * UtX_cand)        # m_cand

  # Schur s_jj^V = d_jj^V - c_j^V' (I_{k-1}^V)^{-1} c_j^V
  if (k_sel > 0L) {
    # c_j^V matrix: (k_sel x m_cand) = (UtX_S)' diag(d_inv) UtX_cand
    c_mat   <- crossprod(UtX_S, d_inv * UtX_cand)            # k_sel x m_cand
    # I_{k-1}^V = (UtX_S)' diag(d_inv) UtX_S
    I_km1   <- crossprod(UtX_S, d_inv * UtX_S)               # k_sel x k_sel
    L_I     <- chol(I_km1 + 1e-8 * diag(k_sel))
    Lt_inv_c <- backsolve(L_I, c_mat, transpose = TRUE)      # k_sel x m_cand
    s_jj_vec <- pmax(d_jj_vec - colSums(Lt_inv_c^2), 1e-10)
    norm_c2  <- colSums(c_mat^2)
    v_jj_vec <- tau2 * (norm_c2 + d_jj_vec^2) / s_jj_vec
  } else {
    s_jj_vec <- d_jj_vec
    v_jj_vec <- tau2 * d_jj_vec
  }

  # Profiled theta_hat under V^{-1}: solve (F' V^{-1} F) theta = F' V^{-1} y
  FtVF <- crossprod(UtF_k, d_inv * UtF_k)                    # k_full x k_full
  FtVy <- crossprod(UtF_k, d_inv * Uty)                      # k_full
  theta_hat <- solve(FtVF + 1e-10 * diag(k_full), FtVy)

  # Residual r = y - F_k theta_hat ; in eigen frame: Utr = Uty - UtF_k theta
  Utr        <- Uty - as.numeric(UtF_k %*% theta_hat)        # n
  # X_j' V^{-1} r  for all j: (UtX_cand)' diag(d_inv) Utr
  Xj_Vinv_r  <- as.numeric(crossprod(UtX_cand, d_inv * Utr)) # m_cand
  # rss_null  = r' V^{-1} r
  rss_null   <- as.numeric(crossprod(Utr, d_inv * Utr))
  sigma2_hat <- rss_null / (n - k_full)

  Z_jV   <- Xj_Vinv_r / sqrt(s_jj_vec * sigma2_hat)
  log_bf <- -0.5 * log1p(v_jj_vec) +
             (v_jj_vec * Z_jV^2) / (2 * (1 + v_jj_vec))

  list(log_BF = log_bf, Z_jV = Z_jV, v_jj = v_jj_vec,
       s_jj = s_jj_vec, d_jj = d_jj_vec,
       sigma2_hat = sigma2_hat,
       rss_null = rss_null,
       log_det_Z = sum(log1p(delta * s)))
}


#' Batched marginalized joint Bayes factor: sweep delta, weight by posterior
#' under H_0 (same as the per-candidate version but vectorized).
lmm_joint_bf_marginalised_k_batched <- function(y, X_cand, F_k, X_S, R_eig,
                                                  tau2 = 1.0,
                                                  prior = "half_cauchy",
                                                  prior_params = list(scale = 1.0),
                                                  n_nodes = 30L,
                                                  delta_lo = 1e-4,
                                                  delta_hi = 1e3) {
  n <- length(y); m_cand <- ncol(X_cand)
  k_full <- ncol(F_k); k_sel <- ncol(X_S)
  U <- R_eig$U; s <- R_eig$s

  # Precomputes (one-time)
  UtX_cand <- crossprod(U, X_cand)
  UtX_S    <- if (k_sel > 0L) crossprod(U, X_S)
              else matrix(0, nrow = n, ncol = 0L)
  Uta      <- as.numeric(crossprod(U, F_k[, 1L]))   # intercept basis
  Uty      <- as.numeric(crossprod(U, y))

  grid    <- make_delta_grid(n_nodes, delta_lo, delta_hi)
  dv      <- grid$delta
  prior_d <- lmm_prior_delta(dv, prior = prior, params = prior_params)
  log_q_w <- log(grid$w_GL) + log(dv) + log(pmax(prior_d, 1e-300))
  log_q_w <- log_q_w - max(log_q_w)

  # Standard log-sum-exp accumulators for unconditional BF
  log_BF_at  <- matrix(0, nrow = m_cand, ncol = n_nodes)
  log_p_h0   <- numeric(n_nodes)   # log p(y | delta, H_0 joint), shared across cand
  for (r in seq_len(n_nodes)) {
    delta <- dv[r]
    d_inv <- 1 / (1 + delta * s)
    # Profile intercept under V^{-1}
    aV    <- as.numeric(crossprod(Uta, d_inv * Uta))  # a^T V^{-1} a
    aVy   <- as.numeric(crossprod(Uta, d_inv * Uty))
    beta0 <- aVy / aV
    Uty_c <- Uty - beta0 * Uta                        # eigen-frame centred y
    RSS_0 <- as.numeric(crossprod(Uty_c, d_inv * Uty_c))

    fit <- joint_bf_exact_k_batched(UtX_cand, UtX_S, Uty_c,
                                       s, delta, tau2, n, k_sel, RSS_0)
    log_BF_at[, r] <- fit$log_BF
    # log p(y | delta, H_0 joint) for the marginalisation weights:
    # |Z_d|^{-1/2} * |I+tau^2 A|^{-1/2} * (a^T V^{-1} a)^{-1/2} * RSS0_joint^{-(n-1)/2}
    # log_det |I+tau^2 A| only depends on shared X_S, recompute here
    if (k_sel > 0L) {
      A_local      <- crossprod(UtX_S, d_inv * UtX_S)
      log_det_IA   <- determinant(diag(k_sel) + tau2 * A_local,
                                     logarithm = TRUE)$modulus
    } else log_det_IA <- 0
    log_p_h0[r] <- -0.5 * fit$log_det_Z -
                    0.5 * log_det_IA -
                    0.5 * log(aV) -
                    0.5 * (n - 1) * log(fit$RSS0_joint)
  }

  # Marginalized log-BF per candidate via posterior over delta.
  # The posterior over delta under H_0 joint:
  #   log_w_post[r] = log_q_w[r] + log_p_h0[r]
  # Marginalized log-BF: log E_{post(delta|y,H0)}[ exp(log_BF(delta)) ]
  # = lse(log_q_w + log_p_h0 + log_BF) - lse(log_q_w + log_p_h0)
  log_w_total <- log_q_w + log_p_h0
  log_w_total <- log_w_total - max(log_w_total)

  lse <- function(x) { M <- max(x); M + log(sum(exp(x - M))) }
  log_BF_marg <- numeric(m_cand)
  log_den <- lse(log_w_total)
  for (j in seq_len(m_cand)) {
    log_BF_marg[j] <- lse(log_w_total + log_BF_at[j, ]) - log_den
  }

  list(log_BF       = log_BF_marg,
       log_BF_at_delta = log_BF_at,
       grid         = dv)
}


# ------------------------------------------------------------------------------
# 4. Stepwise wrapper using the batched joint BF
# ------------------------------------------------------------------------------

JS_L_LMM_stepwise_fast <- function(y, X, a_n = NULL,
                                     p_prior = NULL, tau2 = 0.04,
                                     K_max = 15L,
                                     criterion = c("eBIC", "JointPosterior"),
                                     delta_eval = c("marginal", "reml"),
                                     theta = 0.95,
                                     n_nodes = 20L,
                                     prior = "half_cauchy",
                                     prior_params = list(scale = 1.0)) {
  criterion  <- match.arg(criterion)
  delta_eval <- match.arg(delta_eval)
  # REML-BF operating point: under the plug-in REML evaluation the conditional
  # projected scoring is used (the joint Schur-complement regularisation is a
  # feature of the marginalized MBF path); dispatches to the validated REML-BF.
  if (delta_eval == "reml") {
    fit <- plug_in_reml_eBIC_W(y, X, W = a_n, K_max = K_max, tau2 = tau2)
    return(list(indices = fit$indices,
                q_at_step = rep(NA_real_, fit$K_hat),
                log_BF_at_step = fit$log_BFs,
                ebic_path = fit$ebic_path,
                K_hat = fit$K_hat,
                variant = "reml_bf"))
  }
  n <- nrow(X); m <- ncol(X)
  # a_n may be (i) NULL -> intercept only,
  #            (ii) a numeric vector of length n -> intercept-like,
  #            (iii) a matrix n x k_W -> intercept plus external covariates
  #                  (sex, batch, family, etc.).
  if (is.null(a_n))           a_n <- rep(1, n)
  if (!is.matrix(a_n))        a_n <- matrix(a_n, ncol = 1L)
  if (is.null(p_prior)) p_prior <- rep(0.5, m)

  selected <- integer(0)
  q_at_step <- numeric(K_max)
  log_BF_at_step <- numeric(K_max)
  ebic_path <- numeric(K_max + 1L); ebic_path[1L] <- Inf
  cum_log_BF <- 0                              # cumulative model log-BF
  remaining <- seq_len(m)

  for (k_step in seq_len(K_max)) {
    if (length(remaining) == 0L) break

    F_k <- if (length(selected) == 0L) a_n
           else cbind(a_n, X[, selected, drop = FALSE])
    X_S <- if (length(selected) == 0L) matrix(0, nrow = n, ncol = 0L)
           else X[, selected, drop = FALSE]
    rem_mat <- X[, remaining, drop = FALSE]
    R_step  <- tcrossprod(rem_mat) / max(ncol(rem_mat) - 1L, 1L)
    R_eig   <- eigendecompose_R(R_step)

    bf_all <- lmm_joint_bf_marginalised_k_batched(y, rem_mat, F_k, X_S, R_eig,
                                                    tau2 = tau2,
                                                    prior = prior,
                                                    prior_params = prior_params,
                                                    n_nodes = n_nodes)
    log_BF_vec <- bf_all$log_BF
    p_vec      <- p_prior[remaining]
    q_vec      <- (1 - p_vec) * exp(log_BF_vec) /
                  ((1 - p_vec) * exp(log_BF_vec) + p_vec + 1e-300)

    best_idx <- which.max(log_BF_vec)
    best_j   <- remaining[best_idx]
    best_log_BF <- log_BF_vec[best_idx]
    best_q   <- q_vec[best_idx]

    # eBIC on the cumulative model log-BF (correct formulation):
    #   eBIC(M_k) = -2 log p(y|M_k) + (k+1) (log n + 2 log m)
    # with log p(y|M_k) approximated by sum_{j=1..k} log_BF_j.
    new_cum_log_BF <- cum_log_BF + best_log_BF
    ebic_k <- -2 * new_cum_log_BF + (k_step + 1L) * (log(n) + 2 * log(m))
    if (criterion == "eBIC" && k_step >= 2L && ebic_k > ebic_path[k_step]) break
    if (criterion == "JointPosterior" && best_q < theta) break

    selected <- c(selected, best_j)
    q_at_step[k_step] <- best_q
    log_BF_at_step[k_step] <- best_log_BF
    ebic_path[k_step + 1L] <- ebic_k
    cum_log_BF <- new_cum_log_BF
    remaining <- setdiff(remaining, best_j)
  }

  K_hat <- length(selected)
  list(indices = selected,
       q_at_step = q_at_step[seq_len(K_hat)],
       log_BF_at_step = log_BF_at_step[seq_len(K_hat)],
       ebic_path = ebic_path[seq_len(K_hat + 1L)],
       K_hat = K_hat,
       variant = "joint_fast")
}


# ------------------------------------------------------------------------------
# 5. Original marginal stepwise wrapper (kept for backward compatibility)
# ------------------------------------------------------------------------------

MS_L_LMM_stepwise_fast <- function(y, X, a_n = NULL,
                                     p_prior = NULL, tau2 = 0.04,
                                     K_max = 15L,
                                     criterion = c("eBIC", "JointPosterior"),
                                     delta_eval = c("marginal", "reml"),
                                     theta = 0.95,
                                     n_nodes = 20L,
                                     prior = "half_cauchy",
                                     prior_params = list(scale = 1.0)) {
  criterion  <- match.arg(criterion)
  delta_eval <- match.arg(delta_eval)
  # REML-BF operating point: plug in the per-step REML estimate of delta
  # (balanced default for routine refinement). Dispatches to the validated
  # plug-in REML evaluation; eBIC stopping only.
  if (delta_eval == "reml") {
    fit <- plug_in_reml_eBIC_W(y, X, W = a_n, K_max = K_max, tau2 = tau2)
    return(list(indices = fit$indices,
                q_at_step = rep(NA_real_, fit$K_hat),
                log_BF_at_step = fit$log_BFs,
                ebic_path = fit$ebic_path,
                K_hat = fit$K_hat,
                variant = "reml_bf"))
  }
  n <- nrow(X); m <- ncol(X)
  # a_n may be a vector or a matrix (intercept + external covariates).
  if (is.null(a_n))           a_n <- rep(1, n)
  if (!is.matrix(a_n))        a_n <- matrix(a_n, ncol = 1L)
  if (is.null(p_prior)) p_prior <- rep(0.5, m)

  selected       <- integer(0)
  q_at_step      <- numeric(K_max)
  log_BF_at_step <- numeric(K_max)
  ebic_path      <- numeric(K_max + 1L)
  ebic_path[1L]  <- Inf
  cum_log_BF     <- 0                              # cumulative model log-BF

  remaining <- seq_len(m)

  for (k_step in seq_len(K_max)) {
    if (length(remaining) == 0L) break

    F_k <- if (length(selected) == 0L) matrix(a_n, ncol = 1L)
           else cbind(a_n, X[, selected, drop = FALSE])
    rem_mat <- X[, remaining, drop = FALSE]
    R_step  <- tcrossprod(rem_mat) / max(ncol(rem_mat) - 1L, 1L)
    R_eig   <- eigendecompose_R(R_step)

    # Batched scoring of all remaining candidates
    bf_all <- lmm_bf_marginalised_k_batched(y, rem_mat, F_k, R_eig,
                                              tau2 = tau2,
                                              prior = prior,
                                              prior_params = prior_params,
                                              n_nodes = n_nodes)
    log_BF_vec <- bf_all$log_BF
    p_vec <- p_prior[remaining]
    q_vec <- (1 - p_vec) * exp(log_BF_vec) /
             ((1 - p_vec) * exp(log_BF_vec) + p_vec + 1e-300)

    best_idx     <- which.max(log_BF_vec)
    best_j       <- remaining[best_idx]
    best_log_BF  <- log_BF_vec[best_idx]
    best_q       <- q_vec[best_idx]

    # eBIC on cumulative model log-BF :
    #   eBIC(M_k) = -2 log p(y|M_k) + (k+1)(log n + 2 log m)
    new_cum_log_BF <- cum_log_BF + best_log_BF
    ebic_k <- -2 * new_cum_log_BF + (k_step + 1L) * (log(n) + 2 * log(m))
    if (criterion == "eBIC" && k_step >= 2L && ebic_k > ebic_path[k_step]) break
    if (criterion == "JointPosterior" && best_q < theta) break

    selected <- c(selected, best_j)
    q_at_step[k_step]      <- best_q
    log_BF_at_step[k_step] <- best_log_BF
    ebic_path[k_step + 1L] <- ebic_k
    cum_log_BF <- new_cum_log_BF
    remaining <- setdiff(remaining, best_j)
  }

  K_hat <- length(selected)
  list(indices = selected,
       q_at_step = q_at_step[seq_len(K_hat)],
       log_BF_at_step = log_BF_at_step[seq_len(K_hat)],
       ebic_path = ebic_path[seq_len(K_hat + 1L)],
       K_hat = K_hat,
       variant = "marginal_fast")
}
