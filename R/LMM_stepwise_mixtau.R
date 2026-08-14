# ==============================================================================
# LMM_stepwise_mixtau.R
#
# Mixture-slab variants of the proposed framework: at each step, the
# conditional Bayes factor BF_j(delta, tau^2) is averaged over a discrete
# mixture of slab variances tau_r^2:
#
#     BF_j^mix(delta) = sum_r omega_r * BF_j(delta, tau_r^2)
#
#It addresses the polygenicity-dependent optimal tau^2 observed in the sensitivity analysis
# (Appendix F) by averaging across a range of plausible effect sizes
# instead of committing to a single tau^2 value.
#
# Implementation: thin wrappers around the existing single-tau^2 functions.
# We call the underlying lmm_*_bf_marginalised_k_batched once per tau_r^2,
# then combine the R log-BF vectors via log-sum-exp with mixture weights.
# This makes the wrapper R times slower than the single-tau^2 version
# (R = length(tau2_mix), typically 4-5) but does not require modifying
# the lower-level batched evaluators.
# ==============================================================================

#' Joint variant with mixture slab
#'
#' @param tau2_mix vector of slab variances, e.g. c(0.01, 0.04, 0.10, 0.25)
#' @param weights mixture weights summing to 1; if NULL, uniform
JS_L_LMM_stepwise_fast_mixtau <- function(y, X, a_n = NULL,
                                            p_prior = NULL,
                                            tau2_mix = c(0.01, 0.04, 0.10, 0.25),
                                            weights = NULL,
                                            K_max = 15L,
                                            criterion = c("eBIC", "JointPosterior"),
                                            theta = 0.95,
                                            n_nodes = 20L,
                                            prior = "half_cauchy",
                                            prior_params = list(scale = 1.0)) {
  criterion <- match.arg(criterion)
  n <- nrow(X); m <- ncol(X)
  R <- length(tau2_mix)
  if (is.null(weights)) weights <- rep(1/R, R)
  weights <- weights / sum(weights)
  log_w <- log(weights)
  if (is.null(a_n))    a_n <- rep(1, n)
  if (!is.matrix(a_n)) a_n <- matrix(a_n, ncol = 1L)
  if (is.null(p_prior)) p_prior <- rep(0.5, m)

  selected <- integer(0)
  q_at_step <- numeric(K_max)
  log_BF_at_step <- numeric(K_max)
  ebic_path <- numeric(K_max + 1L); ebic_path[1L] <- Inf
  cum_log_BF <- 0
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

    # Mixture: evaluate BF at each tau_r^2 and combine via log-sum-exp
    log_BFs <- matrix(0, R, length(remaining))
    for (r in seq_len(R)) {
      bf_r <- lmm_joint_bf_marginalised_k_batched(
                  y, rem_mat, F_k, X_S, R_eig,
                  tau2 = tau2_mix[r],
                  prior = prior, prior_params = prior_params,
                  n_nodes = n_nodes)
      log_BFs[r, ] <- bf_r$log_BF
    }
    # log( sum_r omega_r * exp(log_BF_r) ) per candidate
    log_BF_vec <- apply(log_BFs, 2, function(lbfs) {
      x <- lbfs + log_w
      M <- max(x)
      M + log(sum(exp(x - M)))
    })
    p_vec <- p_prior[remaining]
    q_vec <- (1 - p_vec) * exp(log_BF_vec) /
             ((1 - p_vec) * exp(log_BF_vec) + p_vec + 1e-300)

    best_idx <- which.max(log_BF_vec)
    best_j   <- remaining[best_idx]
    best_log_BF <- log_BF_vec[best_idx]
    best_q   <- q_vec[best_idx]

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
       tau2_mix = tau2_mix,
       weights = weights,
       variant = "joint_fast_mixtau")
}


#' Marginal variant with mixture slab (same idea as above)
MS_L_LMM_stepwise_fast_mixtau <- function(y, X, a_n = NULL,
                                            p_prior = NULL,
                                            tau2_mix = c(0.01, 0.04, 0.10, 0.25),
                                            weights = NULL,
                                            K_max = 15L,
                                            criterion = c("eBIC", "JointPosterior"),
                                            theta = 0.95,
                                            n_nodes = 20L,
                                            prior = "half_cauchy",
                                            prior_params = list(scale = 1.0)) {
  criterion <- match.arg(criterion)
  n <- nrow(X); m <- ncol(X)
  R <- length(tau2_mix)
  if (is.null(weights)) weights <- rep(1/R, R)
  weights <- weights / sum(weights)
  log_w <- log(weights)
  if (is.null(a_n))    a_n <- rep(1, n)
  if (!is.matrix(a_n)) a_n <- matrix(a_n, ncol = 1L)
  if (is.null(p_prior)) p_prior <- rep(0.5, m)

  selected <- integer(0)
  q_at_step <- numeric(K_max)
  log_BF_at_step <- numeric(K_max)
  ebic_path <- numeric(K_max + 1L); ebic_path[1L] <- Inf
  cum_log_BF <- 0
  remaining <- seq_len(m)

  for (k_step in seq_len(K_max)) {
    if (length(remaining) == 0L) break
    F_k <- if (length(selected) == 0L) a_n
           else cbind(a_n, X[, selected, drop = FALSE])
    rem_mat <- X[, remaining, drop = FALSE]
    R_step  <- tcrossprod(rem_mat) / max(ncol(rem_mat) - 1L, 1L)
    R_eig   <- eigendecompose_R(R_step)

    log_BFs <- matrix(0, R, length(remaining))
    for (r in seq_len(R)) {
      bf_r <- lmm_bf_marginalised_k_batched(
                  y, rem_mat, F_k, R_eig,
                  tau2 = tau2_mix[r],
                  prior = prior, prior_params = prior_params,
                  n_nodes = n_nodes)
      log_BFs[r, ] <- bf_r$log_BF
    }
    log_BF_vec <- apply(log_BFs, 2, function(lbfs) {
      x <- lbfs + log_w
      M <- max(x)
      M + log(sum(exp(x - M)))
    })
    p_vec <- p_prior[remaining]
    q_vec <- (1 - p_vec) * exp(log_BF_vec) /
             ((1 - p_vec) * exp(log_BF_vec) + p_vec + 1e-300)

    best_idx <- which.max(log_BF_vec)
    best_j   <- remaining[best_idx]
    best_log_BF <- log_BF_vec[best_idx]
    best_q   <- q_vec[best_idx]

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
       tau2_mix = tau2_mix,
       weights = weights,
       variant = "marginal_fast_mixtau")
}
