# ==============================================================================
# LMM_stepwise.R
# Stepwise Bayesian variable selection in the Gaussian linear mixed model.
# Companion to lmm_bayes.tex (theory document).
#
# Implements the algorithm of Section "Stepwise selection algorithm":
#   - At each step k, profile out the k-dim nuisance theta_jk
#     = (beta_0, beta_2_jk) for the selected SNPs S_{k-1}.
#   - Remaining non-selected non-candidate SNPs T_{k-1} \ {j} enter as
#     a polygenic random effect with kernel R_{jk}.
#   - Marginal Bayes factor: per-candidate slab on beta_1 only.
#   - Stepwise wrapper with eBIC and JointPosterior stopping.
#
# Depends on R/LMM_core.R (sourced via the path-resolution block below).
# ==============================================================================

suppressPackageStartupMessages({
  library(statmod)
})

# Source LMM core utilities. Resolution order:
#   1. If LMM_core.R is in the same directory, source it directly.
#   2. Else search common relative paths.
.lmm_source_core <- function() {
  candidates <- c(
    "LMM_core.R",
    "R/LMM_core.R",
    file.path(dirname(sys.frame(1)$ofile %||% "."), "LMM_core.R")
  )
  for (p in candidates) {
    if (!is.null(p) && file.exists(p)) { source(p); return(invisible(p)) }
  }
  stop("LMM_core.R not found. Please source it before sourcing LMM_stepwise.R.")
}
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x
tryCatch(.lmm_source_core(), error = function(e) {
  message("Note: LMM_core.R not auto-sourced. Source it manually before use.")
})


# ------------------------------------------------------------------------------
# 1. RSS quantities at step k, fixed delta — extends rss_at_delta to k>=2
# ------------------------------------------------------------------------------

#' Compute projected RSS quantities at step k, fixed delta.
#' Generalizes rss_at_delta() (k=1 case) by:
#'   - Projecting out the k-dim nuisance F_k = [a_n, X_{S_{k-1}}]
#'   - Working with Z_delta^{-1} computed from R_{jk} (leave-out kernel)
#'
#' @param y         n-vector
#' @param X_j       n-vector candidate
#' @param F_k       n x k nuisance (intercept + selected)
#' @param delta     variance ratio
#' @param R_eig     eigendecomp of R_{jk}
#' @param tau2      slab variance
#' @return list with rss0, rss1, Q_star, sigma2_hat, log_det_Z, D_jk
rss_at_delta_k <- function(y, X_j, F_k, delta, R_eig, tau2 = 1.0) {
  n <- length(y); k <- ncol(F_k)
  # Whiten via Z_delta^{-1/2}
  y_s   <- apply_Zdelta_invsqrt(y,   delta, R_eig)
  Xj_s  <- apply_Zdelta_invsqrt(X_j, delta, R_eig)
  Fk_s  <- apply_Zdelta_invsqrt(F_k, delta, R_eig)

  # Project out F_k_s by GLS (whitened OLS now): theta_hat = (F^T F)^{-1} F^T y
  FtF   <- crossprod(Fk_s)
  Fty_s <- crossprod(Fk_s, y_s)
  Ftxj  <- crossprod(Fk_s, Xj_s)
  L_FtF <- chol(FtF + 1e-10 * diag(k))   # ridge for numerical stability
  theta_hat <- backsolve(L_FtF, forwardsolve(t(L_FtF), Fty_s))

  # Projected residual r_perp = (I - F (F^T F)^{-1} F^T) y_s
  res_s <- as.numeric(y_s - Fk_s %*% theta_hat)
  # Projected candidate Xj_perp = (I - F(F^T F)^{-1}F^T) Xj_s
  Xj_proj_coef <- backsolve(L_FtF, forwardsolve(t(L_FtF), Ftxj))
  Xj_perp <- as.numeric(Xj_s - Fk_s %*% Xj_proj_coef)

  D_jk   <- as.numeric(crossprod(Xj_perp))                    # ≈ n(1 - rho_V)
  rss0   <- as.numeric(crossprod(res_s))
  Xjy_p  <- as.numeric(crossprod(Xj_perp, res_s))
  v      <- tau2 * D_jk
  rss1   <- rss0 - (v * Xjy_p^2) / (D_jk * (1 + v))
  sigma2_hat <- rss0 / (n - k)
  Q_star <- Xjy_p / sqrt(D_jk * sigma2_hat)

  # log |F^T F| from cholesky
  log_det_FtF <- 2 * sum(log(diag(L_FtF)))

  list(rss0 = rss0, rss1 = rss1, Q_star = Q_star,
       sigma2_hat = sigma2_hat,
       log_det_Z = log_det_Zdelta(delta, R_eig),
       log_det_FtF = log_det_FtF,
       D_jk = D_jk, v = v,
       Xjy_perp = Xjy_p)
}


# ------------------------------------------------------------------------------
# 2. Marginal BF at step k (per-candidate slab on beta_1 only)
# ------------------------------------------------------------------------------

#' Conditional Bayes factor at fixed delta for step k (eq. BF_delta_k).
#' BF(delta) = (1+v)^{-1/2} * (RSS_0 / RSS_1)^{(n-k)/2}
lmm_bf_at_delta_k <- function(y, X_j, F_k, delta, R_eig, tau2 = 1.0) {
  rss <- rss_at_delta_k(y, X_j, F_k, delta, R_eig, tau2)
  n <- length(y); k <- ncol(F_k)
  log_bf <- -0.5 * log1p(rss$v) +
             0.5 * (n - k) * (log(rss$rss0) - log(rss$rss1))
  list(log_BF = log_bf, rss = rss)
}

#' Marginalized BF at step k, integrating over delta.
lmm_bf_marginalised_k <- function(y, X_j, F_k, R_eig, tau2 = 1.0,
                                    prior        = "half_cauchy",
                                    prior_params = list(scale = 1.0),
                                    n_nodes      = 30L,
                                    delta_lo     = 1e-4,
                                    delta_hi     = 1e3) {
  n <- length(y); k <- ncol(F_k)
  grid    <- make_delta_grid(n_nodes, delta_lo, delta_hi)
  dv      <- grid$delta
  prior_d <- lmm_prior_delta(dv, prior = prior, params = prior_params)
  log_q_w <- log(grid$w_GL) + log(dv) + log(pmax(prior_d, 1e-300))
  log_q_w <- log_q_w - max(log_q_w)

  log_p_h0 <- numeric(n_nodes)
  log_p_h1 <- numeric(n_nodes)
  log_BF_at <- numeric(n_nodes)
  for (r in seq_len(n_nodes)) {
    rss <- rss_at_delta_k(y, X_j, F_k, dv[r], R_eig, tau2)
    base <- -0.5 * rss$log_det_Z - 0.5 * rss$log_det_FtF
    log_p_h0[r] <- base - 0.5 * (n - k) * log(rss$rss0)
    log_p_h1[r] <- base - 0.5 * (n - k) * log(rss$rss1) - 0.5 * log1p(rss$v)
    log_BF_at[r] <- log_p_h1[r] - log_p_h0[r]
  }

  lse <- function(x) { M <- max(x); M + log(sum(exp(x - M))) }
  log_BF <- lse(log_q_w + log_p_h1) - lse(log_q_w + log_p_h0)

  list(log_BF = log_BF,
       log_BF_at_delta = log_BF_at,
       grid = dv)
}


# ------------------------------------------------------------------------------
# 3. Joint Bayes factor at step k (Proposition prop:joint_BF in lmm_bayes.tex)
# ------------------------------------------------------------------------------
#
# Joint variant: full lambda = (beta_2_jk, beta_1_jk) ~ N(0, tau^2 sigma^2 I_k)
# rather than profiling beta_2 with a flat prior.
# Closed-form via Schur complement of the augmented information matrix.
#
# At fixed delta:
#   I_{k-1}^V = n^{-1} X_S^T Z_delta^{-1} X_S    (shared across candidates)
#   c_j^V     = n^{-1} X_S^T Z_delta^{-1} X_j
#   d_jj^V    = n^{-1} X_j^T Z_delta^{-1} X_j
#   s_jj^V    = d_jj^V - c_j^{V T} (I_{k-1}^V)^{-1} c_j^V    (Schur)
#   v_jj^V    = tau^2 (||c_j^V||^2 + (d_jj^V)^2) / s_jj^V
#   Z_j^V     = X_j^T Z^{-1}(y - F_k theta_hat) / sqrt(n s_jj^V sigma2_hat)
#   logBF     = -1/2 log(1 + v_jj^V) + v_jj^V * Z_j^V^2 / (2 (1 + v_jj^V))
#
# At k=1, X_S is empty, c_j^V = 0, s_jj^V = d_jj^V, v_jj^V = tau^2 d_jj^V.

#' Joint Bayes factor at fixed delta for step k (eq. joint_BF).
#'
#' @param y       n-vector response
#' @param X_j     n-vector candidate
#' @param F_k     n x k full nuisance design (intercept + selected SNPs)
#' @param X_S     n x (k-1) matrix of SELECTED SNPs only (NOT intercept).
#'                Pass an n x 0 matrix when k=1.
#' @param delta   variance ratio
#' @param R_eig   eigendecomp of R_{jk}
#' @param tau2    g-prior variance
#' @return list(log_BF, Z_jV, v_jj, s_jj, d_jj)
lmm_joint_bf_at_delta_k <- function(y, X_j, F_k, X_S, delta, R_eig, tau2 = 1.0) {
  # Exact joint Bayes factor at fixed delta via Sherman--Morrison--Woodbury.
  # Reference: eq.(joint_BF) of theory/lmm_bayes.tex Section 5.
  #   BF_jk^Joint(delta) = (1 + tau^2 * s_jj^{V,tau})^{-1/2}
  #                       * (RSS_0^tau / RSS_1^tau)^{(n-1)/2}
  # where s_jj^{V,tau} = d_jj^V - c_j^{V,T} (I_{k-1}^V + tau^{-2} I)^{-1} c_j^V
  # is the prior-regularized Schur complement, and RSS_0^tau, RSS_1^tau are
  # prior-augmented residual sums of squares.
  #
  # The intercept block (a_n column of F_k) keeps a flat prior; only the
  # X_S columns of F_k carry the joint Gaussian g-prior.
  n <- length(y); k_sel <- ncol(X_S); k_full <- ncol(F_k)

  # Project out intercept first (flat prior): work in F_perp metric.
  # F_k = [a_n, X_S]. Profile out a_n (column 1) under flat prior.
  Zinv_F <- apply_Zdelta_inv(F_k, delta, R_eig)
  Zinv_y <- apply_Zdelta_inv(y,   delta, R_eig)
  Zinv_Xj <- apply_Zdelta_inv(X_j, delta, R_eig)

  # Profile-likelihood residual under H_0 = {a_n} (intercept only)
  a_n <- F_k[, 1L, drop = FALSE]
  Zinv_a <- apply_Zdelta_inv(a_n, delta, R_eig)
  aZa <- as.numeric(crossprod(a_n, Zinv_a))
  alpha_hat <- as.numeric(crossprod(a_n, Zinv_y)) / aZa
  r_a <- as.numeric(y - alpha_hat * a_n)
  Zinv_ra <- apply_Zdelta_inv(r_a, delta, R_eig)
  RSS_0 <- as.numeric(crossprod(r_a, Zinv_ra))   # RSS after intercept profile

  # Information quantities in V = Z_delta metric (centred basis):
  d_jj <- as.numeric(crossprod(X_j, Zinv_Xj))    # ~ n

  if (k_sel > 0L) {
    Zinv_XS <- apply_Zdelta_inv(X_S, delta, R_eig)
    I_km1   <- crossprod(X_S, Zinv_XS)            # I_{k-1}^V
    c_j     <- as.numeric(crossprod(X_S, Zinv_Xj))  # c_j^V
    p_vec   <- as.numeric(crossprod(X_S, Zinv_ra))  # X_S^T Z^{-1} r_a

    # M_aug = I + tau^2 * I_{k-1}^V; v = M_aug^{-1} p
    M_aug   <- diag(k_sel) + tau2 * I_km1
    L_aug   <- chol(M_aug + 1e-10 * diag(k_sel))
    v_vec   <- backsolve(L_aug, forwardsolve(t(L_aug), p_vec))
    pT_v    <- as.numeric(crossprod(p_vec, v_vec))

    # u = M_aug^{-1} c_j;  bT_u = c_j^T u;  bT_v = c_j^T v
    u_vec   <- backsolve(L_aug, forwardsolve(t(L_aug), c_j))
    bT_u    <- as.numeric(crossprod(c_j, u_vec))
    bT_v    <- as.numeric(crossprod(c_j, v_vec))
  } else {
    pT_v <- 0; bT_u <- 0; bT_v <- 0
  }

  # q_j = X_j^T Z^{-1} r_a
  q_j <- as.numeric(crossprod(X_j, Zinv_ra))

  # S = 1 + tau^2 (d_jj - tau^2 c^T M_aug^{-1} c) = 1 + tau^2 s_jj^{V,tau}
  S_val   <- 1 + tau2 * (d_jj - tau2 * bT_u)
  num     <- q_j - tau2 * bT_v          # numerator of the SMW update
  RSS0_t  <- RSS_0 - tau2 * pT_v        # prior-augmented RSS_0
  RSS1_t  <- pmax(RSS0_t - tau2 * num^2 / S_val, 1e-300)

  log_bf  <- -0.5 * log(S_val) + 0.5 * (n - 1) * (log(RSS0_t) - log(RSS1_t))

  # Diagnostic Z_jV (kept for backward compatibility with V5/V6)
  s_jj    <- d_jj - if (k_sel > 0L) tau2 * bT_u else 0
  sigma2_hat <- RSS0_t / (n - 1)
  Z_jV    <- num / sqrt(pmax(s_jj * sigma2_hat, 1e-300))

  # Pieces needed for the H_0 marginal-likelihood weight in the
  # marginalization wrapper:
  #   log p(y | delta, H_0 joint) = -0.5*log|Z|
  #                                 -0.5*log|I+tau^2 I_{k-1}^V|
  #                                 -0.5*log(a^T Z^{-1} a)
  #                                 -0.5*(n-1)*log(RSS_0^tau)
  log_det_IA <- if (k_sel > 0L) {
    as.numeric(determinant(M_aug, logarithm = TRUE)$modulus)
  } else 0

  list(log_BF = log_bf, Z_jV = Z_jV, S = S_val,
       s_jj = s_jj, d_jj = d_jj,
       sigma2_hat = sigma2_hat,
       RSS_0 = RSS0_t, RSS_1 = RSS1_t,
       log_det_Z = log_det_Zdelta(delta, R_eig),
       log_det_IA = log_det_IA, log_aV = log(aZa))
}


#' Joint Bayes factor marginalized over delta.
#' Same quadrature scheme as the marginal variant, but using
#' lmm_joint_bf_at_delta_k() for each grid point and weighting by the
#' appropriate likelihood normalisation.
lmm_joint_bf_marginalised_k <- function(y, X_j, F_k, X_S, R_eig, tau2 = 1.0,
                                          prior = "half_cauchy",
                                          prior_params = list(scale = 1.0),
                                          n_nodes = 30L,
                                          delta_lo = 1e-4,
                                          delta_hi = 1e3) {
  n <- length(y); k_full <- ncol(F_k)
  grid    <- make_delta_grid(n_nodes, delta_lo, delta_hi)
  dv      <- grid$delta
  prior_d <- lmm_prior_delta(dv, prior = prior, params = prior_params)
  log_q_w <- log(grid$w_GL) + log(dv) + log(pmax(prior_d, 1e-300))
  log_q_w <- log_q_w - max(log_q_w)

  # Marginalization: BF_marg = E_{posterior(delta | y, H_0 joint)} [ BF(delta) ]
  # The H_0 joint marginal likelihood incorporates:
  #   p(y|delta, H_0_joint) ~ |Z|^{-1/2} |I+tau^2 I_{k-1}^V|^{-1/2}
  #                          (a^T Z^{-1} a)^{-1/2} (RSS_0^tau)^{-(n-1)/2}
  # (intercept under flat prior, X_S under joint Gaussian g-prior, sigma^2
  # under Jeffreys; cf. eq.(joint_BF) of theory/lmm_bayes.tex Section 5.)
  log_p_h0    <- numeric(n_nodes)   # H_0 joint marginal log-likelihood
  log_BF_at   <- numeric(n_nodes)
  for (r in seq_len(n_nodes)) {
    fit <- lmm_joint_bf_at_delta_k(y, X_j, F_k, X_S, dv[r], R_eig, tau2)
    log_BF_at[r] <- fit$log_BF
    log_p_h0[r]  <- -0.5 * fit$log_det_Z -
                    0.5 * fit$log_det_IA -
                    0.5 * fit$log_aV -
                    0.5 * (n - 1) * log(fit$RSS_0)
  }

  # Proper marginal Bayes factor: ratio of marginal likelihoods, not E[BF].
  #   BF^marg = int p(y|H_1, delta) pi(delta) d delta /
  #             int p(y|H_0, delta) pi(delta) d delta
  # log p(y|H_1, delta) = log p(y|H_0, delta) + log BF(delta)
  lse <- function(x) { M <- max(x); M + log(sum(exp(x - M))) }
  log_BF_marg <- lse(log_q_w + log_p_h0 + log_BF_at) -
                 lse(log_q_w + log_p_h0)

  log_w_post  <- exp(log_q_w + log_p_h0 - lse(log_q_w + log_p_h0))
  list(log_BF       = log_BF_marg,
       log_BF_at_delta = log_BF_at,
       posterior_w  = log_w_post,
       grid         = dv)
}


# ------------------------------------------------------------------------------
# 4. Stepwise wrapper — MarginalScore_L_LMM with eBIC stopping
# ------------------------------------------------------------------------------

#' Stepwise Bayesian LMM-MarginalScore selection.
#' @param y          n-vector response
#' @param X          n x m genotype matrix (centred-scaled)
#' @param a_n        n-vector intercept basis (default rep(1, n))
#' @param p_prior    m-vector prior probability of NULLITY
#' @param tau2       g-prior variance
#' @param K_max      maximum number of steps
#' @param criterion  "eBIC" or "JointPosterior"
#' @param theta      JointPosterior threshold (default 0.95)
#' @param n_nodes    quadrature nodes per step
#' @param prior      prior on delta
#' @param prior_params list of prior params
#' @return list(indices, q_at_step, log_BF_at_step, ebic_path, K_hat)
MS_L_LMM_stepwise <- function(y, X, a_n = NULL,
                                p_prior = NULL, tau2 = 0.04,
                                K_max = 15L,
                                criterion = c("eBIC", "JointPosterior"),
                                theta = 0.95,
                                n_nodes = 20L,
                                prior = "half_cauchy",
                                prior_params = list(scale = 1.0)) {
  criterion <- match.arg(criterion)
  n <- nrow(X); m <- ncol(X)
  if (is.null(a_n)) a_n <- rep(1, n)
  if (is.null(p_prior)) p_prior <- rep(0.5, m)

  selected <- integer(0)
  q_at_step <- numeric(K_max)
  log_BF_at_step <- numeric(K_max)
  ebic_path <- numeric(K_max + 1L)
  ebic_path[1L] <- Inf   # placeholder before any selection
  cum_log_BF <- 0        # cumulative model log-BF (correct eBIC)

  remaining <- seq_len(m)

  for (k_step in seq_len(K_max)) {
    if (length(remaining) == 0L) break

    F_k <- if (length(selected) == 0L) matrix(a_n, ncol = 1L)
           else cbind(a_n, X[, selected, drop = FALSE])

    # Build leave-out kernel for this step (Option B: candidate excluded
    # AND selected excluded). Use approximate full-remaining kernel to
    # share eigendecomposition across candidates within the step.
    rem_mat <- X[, remaining, drop = FALSE]
    R_step <- tcrossprod(rem_mat) / max(ncol(rem_mat) - 1L, 1L)
    R_eig <- eigendecompose_R(R_step)

    # Score each remaining candidate
    best_log_BF <- -Inf; best_j <- NA_integer_; best_q <- 0
    for (idx in seq_along(remaining)) {
      j_orig <- remaining[idx]
      X_j <- X[, j_orig]
      bf <- lmm_bf_marginalised_k(y, X_j, F_k, R_eig, tau2 = tau2,
                                    prior = prior,
                                    prior_params = prior_params,
                                    n_nodes = n_nodes)
      logBF_j <- bf$log_BF
      p_j <- p_prior[j_orig]
      q_j <- (1 - p_j) * exp(logBF_j) /
             ((1 - p_j) * exp(logBF_j) + p_j + 1e-300)
      if (logBF_j > best_log_BF) {
        best_log_BF <- logBF_j; best_j <- j_orig; best_q <- q_j
      }
    }

    # eBIC on CUMULATIVE model log-BF (correct):
    #   eBIC(M_k) = -2 log p(y|M_k) + (k+1)(log n + 2 log m)
    # Approximate log p(y|M_k) by sum_j log_BF_j (cumulative sum).
    new_cum_log_BF <- cum_log_BF + best_log_BF
    ebic_k <- -2 * new_cum_log_BF + (k_step + 1L) * (log(n) + 2 * log(m))
    if (criterion == "eBIC" && k_step >= 2L && ebic_k > ebic_path[k_step]) {
      break   # eBIC increased -> stop
    }
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
       variant = "marginal")
}


# ------------------------------------------------------------------------------
# 5. Stepwise wrapper — JS_L_LMM (joint variant) using Schur-form joint BF
# ------------------------------------------------------------------------------

#' Stepwise Bayesian LMM-JointScore selection.
#' Identical to MS_L_LMM_stepwise() but uses lmm_joint_bf_marginalised_k()
#' (Proposition prop:joint_BF in lmm_bayes.tex), giving a Schur-complement
#' joint Bayes factor under the augmented Gaussian prior on
#' (beta_2_jk, beta_1_jk).
#'
#' @inheritParams MS_L_LMM_stepwise
#' @return same structure as MS_L_LMM_stepwise() with variant = "joint"
JS_L_LMM_stepwise <- function(y, X, a_n = NULL,
                                p_prior = NULL, tau2 = 0.04,
                                K_max = 15L,
                                criterion = c("eBIC", "JointPosterior"),
                                theta = 0.95,
                                n_nodes = 20L,
                                prior = "half_cauchy",
                                prior_params = list(scale = 1.0)) {
  criterion <- match.arg(criterion)
  n <- nrow(X); m <- ncol(X)
  if (is.null(a_n))     a_n <- rep(1, n)
  if (is.null(p_prior)) p_prior <- rep(0.5, m)

  selected       <- integer(0)
  q_at_step      <- numeric(K_max)
  log_BF_at_step <- numeric(K_max)
  ebic_path      <- numeric(K_max + 1L)
  ebic_path[1L]  <- Inf
  cum_log_BF     <- 0      # cumulative model log-BF (correct eBIC)

  remaining <- seq_len(m)

  for (k_step in seq_len(K_max)) {
    if (length(remaining) == 0L) break

    F_k <- if (length(selected) == 0L) matrix(a_n, ncol = 1L)
           else cbind(a_n, X[, selected, drop = FALSE])
    X_S <- if (length(selected) == 0L) matrix(0, nrow = n, ncol = 0L)
           else X[, selected, drop = FALSE]

    # Random-effect kernel: same Option B convention (selected and candidate
    # both removed, but kernel shared across candidates in the step).
    rem_mat <- X[, remaining, drop = FALSE]
    R_step  <- tcrossprod(rem_mat) / max(ncol(rem_mat) - 1L, 1L)
    R_eig   <- eigendecompose_R(R_step)

    best_log_BF <- -Inf; best_j <- NA_integer_; best_q <- 0
    for (idx in seq_along(remaining)) {
      j_orig <- remaining[idx]
      X_j <- X[, j_orig]
      bf <- lmm_joint_bf_marginalised_k(y, X_j, F_k, X_S, R_eig,
                                          tau2 = tau2,
                                          prior = prior,
                                          prior_params = prior_params,
                                          n_nodes = n_nodes)
      logBF_j <- bf$log_BF
      p_j <- p_prior[j_orig]
      q_j <- (1 - p_j) * exp(logBF_j) /
             ((1 - p_j) * exp(logBF_j) + p_j + 1e-300)
      if (logBF_j > best_log_BF) {
        best_log_BF <- logBF_j; best_j <- j_orig; best_q <- q_j
      }
    }

    # eBIC on CUMULATIVE model log-BF (correct).
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
       variant = "joint")
}
