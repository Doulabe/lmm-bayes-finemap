# ==============================================================================
# LMM_core.R
# Core utilities for Bayesian variable selection in the Gaussian linear
# mixed model. 
#
# Functions provided:
#   - eigendecompose_R()     : eigendecomposition of the random-effect kernel
#                              R = X X^T / m (or leave-out variant).
#   - apply_Zdelta_invsqrt() : Z_delta^{-1/2} A via cached eigendecomp.
#   - log_det_Zdelta()       : log |Z_delta|.
#   - rss_at_delta()         : RSS_0(delta), RSS_1(delta) at step k=1.
#   - lmm_bf_at_delta()      : conditional Bayes factor at fixed delta.
#   - lmm_reml_delta()       : REML estimate of delta via Brent's method
#                              (1-D root-finding on profile log-likelihood).
#   - lmm_prior_delta()      : prior densities (half-Cauchy / log-uniform /
#                              inverse-gamma).
#   - make_delta_grid()      : Gauss-Legendre quadrature grid on log(delta).
#   - lmm_bf_marginalised()  : Bayes factor marginalised over delta.
#   - lmm_score_stat()       : LMM-conditioned Score statistic Q^*(delta).
# ==============================================================================

suppressPackageStartupMessages({
  library(statmod)   # gauss.quad
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x


# ------------------------------------------------------------------------------
# 1. Eigendecomposition (shared across all delta values)
# ------------------------------------------------------------------------------

#' Eigendecomposition of a symmetric PSD matrix R = U S U^T.
#' @param R  n x n GRM (or any symmetric PSD).
#' @param floor_eigval  clamp eigenvalues below this (numerical stability).
#' @return list(U, s) with R = U %*% diag(s) %*% t(U).
eigendecompose_R <- function(R, floor_eigval = 1e-10) {
  e <- eigen(R, symmetric = TRUE)
  list(U = e$vectors, s = pmax(e$values, floor_eigval))
}


# ------------------------------------------------------------------------------
# 2. Whitening transformation y -> y* via eigendecomposition
# ------------------------------------------------------------------------------

#' Apply Z_delta^{-1/2} to a vector or matrix using cached eigendecomp.
#' Z_delta = I + delta * R = U (I + delta * S) U^T
#' Z_delta^{-1/2} = U diag(1/sqrt(1 + delta*s)) U^T
apply_Zdelta_invsqrt <- function(A, delta, R_eig) {
  d <- 1 / sqrt(1 + delta * R_eig$s)
  U <- R_eig$U
  if (is.null(dim(A))) {
    UtA <- as.numeric(crossprod(U, A))
    as.numeric(U %*% (d * UtA))
  } else {
    UtA <- crossprod(U, A)
    U %*% (d * UtA)
  }
}

#' log |Z_delta| = sum_i log(1 + delta * s_i)
log_det_Zdelta <- function(delta, R_eig) {
  sum(log1p(delta * R_eig$s))
}

#' Apply Z_delta^{-1} to a vector or matrix using cached eigendecomp.
#' Z_delta^{-1} = U diag(1/(1 + delta*s)) U^T
#' Used by the joint Bayes-factor variant (Section "Joint Bayesian variant"
#' of lmm_bayes.tex) where we work directly with V^{-1} = Z_delta^{-1}/sigma^2.
apply_Zdelta_inv <- function(A, delta, R_eig) {
  d <- 1 / (1 + delta * R_eig$s)
  U <- R_eig$U
  if (is.null(dim(A))) {
    UtA <- as.numeric(crossprod(U, A))
    as.numeric(U %*% (d * UtA))
  } else {
    UtA <- crossprod(U, A)
    U %*% (d * UtA)
  }
}


# ------------------------------------------------------------------------------
# 3. RSS quantities at fixed delta
# ------------------------------------------------------------------------------

#' Compute RSS_0(delta), RSS_1(delta), and the LMM score Q*(delta) at fixed delta.
#'
#' Inputs are *original* (un-whitened) y, X_j, and the eigendecomposition of R.
#' The whitened quantities are computed inside.
#'
#' @param y        n-vector response
#' @param X_j      n-vector candidate genotype (centred, standardised)
#' @param a_n      n-vector all-ones (or any orthogonal-to-X intercept basis)
#' @param delta    scalar delta value
#' @param R_eig    eigendecomp of R
#' @param tau2     prior variance scaling (g-prior); v = tau2 * n
#' @return list(rss0, rss1, Q_star, sigma2_hat, log_det_Z)
rss_at_delta <- function(y, X_j, a_n, delta, R_eig, tau2 = 1.0) {
  n <- length(y)
  # Whitening
  y_s   <- apply_Zdelta_invsqrt(y,   delta, R_eig)
  Xj_s  <- apply_Zdelta_invsqrt(X_j, delta, R_eig)
  a_s   <- apply_Zdelta_invsqrt(a_n, delta, R_eig)

  # Profile out beta_0 via OLS in whitened space:
  #   beta_0_hat = (a_s^T y_s) / (a_s^T a_s)
  ata     <- as.numeric(crossprod(a_s))
  beta0   <- as.numeric(crossprod(a_s, y_s)) / ata
  res     <- y_s - a_s * beta0

  rss0    <- as.numeric(crossprod(res))                 # RSS under H0
  Xj_y    <- as.numeric(crossprod(Xj_s, res))           # X^*T y_centred
  XjtXj   <- as.numeric(crossprod(Xj_s))                # approx n; see Lemma orthogonality
  v       <- tau2 * XjtXj
  rss1    <- rss0 - (v * Xj_y^2) / (XjtXj * (1 + v))    # RSS under H1 slab
  sigma2_hat <- rss0 / n
  Q_star  <- Xj_y / sqrt(n * sigma2_hat)
  list(rss0 = rss0, rss1 = rss1, Q_star = Q_star,
       sigma2_hat = sigma2_hat,
       log_det_Z = log_det_Zdelta(delta, R_eig),
       v = v, XjtXj = XjtXj, ata = ata)
}


# ------------------------------------------------------------------------------
# 4. BF at fixed delta — eq. BF_delta
# ------------------------------------------------------------------------------

#' Closed-form Bayes factor at fixed delta (eq. \ref{eq:BF_delta} in the LaTeX).
#' BF(delta) = (1+v)^{-1/2} * (RSS_0/RSS_1)^{(n-1)/2}
#' Returns the LOG Bayes factor.
lmm_bf_at_delta <- function(y, X_j, a_n, delta, R_eig, tau2 = 1.0) {
  rss <- rss_at_delta(y, X_j, a_n, delta, R_eig, tau2)
  n <- length(y)
  log_bf <- -0.5 * log1p(rss$v) +
             0.5 * (n - 1) * (log(rss$rss0) - log(rss$rss1))
  list(log_BF = log_bf, rss = rss)
}


# ------------------------------------------------------------------------------
# 5. REML for delta (1-D Brent on profile log-likelihood)
# ------------------------------------------------------------------------------

#' REML estimate of delta under H_0 by Brent search on log(delta).
#' Profiled REML log-likelihood:
#'   ell(delta) = -0.5 log|Z_delta| - 0.5 (n-1) log RSS_0(delta) - 0.5 log(a_s^T a_s)
lmm_reml_delta <- function(y, X_j, a_n, R_eig,
                            log_delta_lo = log(1e-4),
                            log_delta_hi = log(1e3)) {
  neg_ll <- function(log_delta) {
    rss <- rss_at_delta(y, X_j, a_n, exp(log_delta), R_eig, tau2 = 1)
    n <- length(y)
    0.5 * rss$log_det_Z +
      0.5 * (n - 1) * log(rss$rss0) +
      0.5 * log(rss$ata)
  }
  opt <- optimise(neg_ll, c(log_delta_lo, log_delta_hi), tol = 1e-5)
  exp(opt$minimum)
}


# ------------------------------------------------------------------------------
# 6. Prior densities for delta
# ------------------------------------------------------------------------------

#' Prior density on delta. Half-Cauchy on sqrt(delta) by default.
lmm_prior_delta <- function(delta,
                              prior  = c("half_cauchy", "log_uniform", "inv_gamma"),
                              params = list()) {
  prior <- match.arg(prior)
  if (prior == "half_cauchy") {
    s <- params$scale %||% 1.0
    sd_ <- sqrt(pmax(delta, 0))
    dens_sd <- 2 / (pi * s * (1 + (sd_ / s)^2))
    dens_sd / (2 * pmax(sd_, 1e-12))
  } else if (prior == "log_uniform") {
    lo <- params$lo %||% 1e-4
    hi <- params$hi %||% 1e3
    ifelse(delta >= lo & delta <= hi, 1 / (delta * log(hi / lo)), 0)
  } else {
    a <- params$shape %||% 1
    b <- params$scale %||% 1
    b^a / gamma(a) * delta^(-a - 1) * exp(-b / delta)
  }
}


# ------------------------------------------------------------------------------
# 7. Quadrature grid on log(delta)
# ------------------------------------------------------------------------------

#' Build Gauss-Legendre quadrature grid on log(delta).
#' Returns nodes (delta values) and Gauss-Legendre weights on the
#' u = log(delta) scale.
make_delta_grid <- function(n_nodes = 30L,
                              delta_lo = 1e-4,
                              delta_hi = 1e3) {
  gl <- statmod::gauss.quad(n_nodes, kind = "legendre")
  u_lo <- log(delta_lo); u_hi <- log(delta_hi)
  u    <- 0.5 * (u_hi - u_lo) * gl$nodes + 0.5 * (u_hi + u_lo)
  w_GL <- 0.5 * (u_hi - u_lo) * gl$weights
  list(delta = exp(u), w_GL = w_GL)
}


# ------------------------------------------------------------------------------
# 8. Marginalized Bayes factor (k=1 main entry point)
# ------------------------------------------------------------------------------

#' Compute the marginalized Bayes factor
#'   BF^Bayes = [ int p(y|delta, H_1) pi(delta) ddelta ]
#'           / [ int p(y|delta, H_0) pi(delta) ddelta ]
#' via 1-D Gauss-Legendre quadrature on log(delta).
#'
#' For numerical stability, we use log-space and the log-sum-exp trick after
#' to combine grid points with their (highly variable) weights.
#'
#' @param y         n-vector response
#' @param X_j       n-vector candidate
#' @param a_n       n-vector intercept basis (typically all-ones)
#' @param R_eig     eigendecomp of R = X^{-j} X^{-j T}/m
#' @param tau2      prior variance scaling (g-prior)
#' @param prior     prior on delta
#' @param prior_params list of prior parameters
#' @param n_nodes   GL nodes (default 30)
#' @param delta_lo,delta_hi quadrature bounds
#' @return list(log_BF, log_BF_at_delta, weights, grid)
lmm_bf_marginalised <- function(y, X_j, a_n, R_eig, tau2 = 1.0,
                                  prior        = "half_cauchy",
                                  prior_params = list(scale = 1.0),
                                  n_nodes      = 30L,
                                  delta_lo     = 1e-4,
                                  delta_hi     = 1e3) {
  n     <- length(y)
  grid  <- make_delta_grid(n_nodes, delta_lo, delta_hi)
  dv    <- grid$delta
  prior_d <- lmm_prior_delta(dv, prior = prior, params = prior_params)
  # Quadrature weight density on log(delta), incl. Jacobian delta_r:
  log_q_w <- log(grid$w_GL) + log(dv) + log(pmax(prior_d, 1e-300))
  log_q_w <- log_q_w - max(log_q_w)   # stability shift; cancels in BF

  # log p(y | delta_r, H_h)  (un-normalized; constants cancel in BF)
  
  log_p_h0 <- numeric(n_nodes)
  log_p_h1 <- numeric(n_nodes)
  log_BF_at <- numeric(n_nodes)
  for (r in seq_len(n_nodes)) {
    rss <- rss_at_delta(y, X_j, a_n, dv[r], R_eig, tau2)
    log_p_h0[r] <- -0.5 * rss$log_det_Z - 0.5 * log(rss$ata) -
                   0.5 * (n - 1) * log(rss$rss0)
    log_p_h1[r] <- -0.5 * rss$log_det_Z - 0.5 * log(rss$ata) -
                   0.5 * (n - 1) * log(rss$rss1) -
                   0.5 * log1p(rss$v)
    log_BF_at[r] <- log_p_h1[r] - log_p_h0[r]
  }

  # log-sum-exp for numerator/denominator
  lse <- function(x) {
    M <- max(x); M + log(sum(exp(x - M)))
  }
  log_num <- lse(log_q_w + log_p_h1)
  log_den <- lse(log_q_w + log_p_h0)
  log_BF  <- log_num - log_den

  list(log_BF = log_BF,
       log_BF_at_delta = log_BF_at,
       grid = dv,
       log_quad_weights = log_q_w + log_p_h0,
       posterior_delta_weights =
         exp(log_q_w + log_p_h0 - lse(log_q_w + log_p_h0)))
}


# ------------------------------------------------------------------------------
# 9. Score statistic Q*(delta) — diagnostic helper
# ------------------------------------------------------------------------------

lmm_score_stat <- function(y, X_j, a_n, delta, R_eig) {
  rss <- rss_at_delta(y, X_j, a_n, delta, R_eig, tau2 = 1)
  rss$Q_star
}
