# ==============================================================================
# LMM_stepwise_pathavg.R
#
# Path-averaged variants of the framework 
#Instead of following a single greedy stepwise path, we maintain a
# beam of B parallel paths, each with a cumulative log-evidence weight.
# After K_max steps, the marginal posterior inclusion probability for each
# SNP j is computed by summing path weights over all beam paths that
# contain j:
#
#     q^marginal_j = sum_{p in beam : j in p.S} exp(p.log_w)
#                  / sum_{p in beam} exp(p.log_w)
#
# This gives an approximation of P(gamma_j = 1 | y) restricted to the
# greedy-reachable model space, addressing the conditional-vs-marginal
# critique on Brier/ECE calibration.
#
# Implementation notes:
# - Beam width B is a key hyperparameter; B=10-20 typical
# - Deduplication: two paths can reach the same S; we merge their weights
#   via log-sum-exp before keeping the top-B
# - Path weight = sum of log-BFs along the path (cumulative model evidence
#   under uniform prior over models)
# - Cost: ~B times the single-path stepwise cost
# ==============================================================================

#' Path-averaged marginal-variant stepwise selection (PA-MS-L)
#'
#' @param B beam width (number of parallel paths to maintain)
#' @param K_max max selection steps
#' @param prior_inclusion prior P(gamma_j = 1) for posterior conversion
PA_MS_L_LMM_stepwise_fast <- function(y, X, a_n = NULL,
                                          tau2 = 0.04,
                                          B = 10L,
                                          K_max = 10L,
                                          n_nodes = 15L,
                                          prior = "half_cauchy",
                                          prior_params = list(scale = 1.0),
                                          prior_inclusion = 0.5,
                                          verbose = FALSE) {
  n <- nrow(X); m <- ncol(X)
  if (is.null(a_n))    a_n <- rep(1, n)
  if (!is.matrix(a_n)) a_n <- matrix(a_n, ncol = 1L)

  # Initialize beam with empty selection at log_w = 0
  beam <- list(list(S = integer(0), log_w = 0,
                       log_BFs_path = numeric(0)))

  for (k_step in seq_len(K_max)) {
    if (verbose) cat(sprintf("  step %d: beam size = %d\n",
                                k_step, length(beam)))
    extensions <- list()

    for (path in beam) {
      F_k <- if (length(path$S) == 0L) a_n
             else cbind(a_n, X[, path$S, drop = FALSE])
      remaining <- setdiff(seq_len(m), path$S)
      if (length(remaining) == 0L) next
      rem_mat <- X[, remaining, drop = FALSE]
      R_step  <- tcrossprod(rem_mat) / max(ncol(rem_mat) - 1L, 1L)
      R_eig   <- eigendecompose_R(R_step)

      bf_all <- lmm_bf_marginalised_k_batched(y, rem_mat, F_k, R_eig,
                                                  tau2 = tau2,
                                                  prior = prior,
                                                  prior_params = prior_params,
                                                  n_nodes = n_nodes)
      log_BFs <- bf_all$log_BF

      # Generate extensions (top-B locally per path for efficiency)
      # Only keep candidates with sufficient evidence to make beam later
      local_top <- order(-log_BFs)[seq_len(min(B, length(log_BFs)))]
      for (i in local_top) {
        j <- remaining[i]
        new_S <- sort(c(path$S, j))
        new_log_w <- path$log_w + log_BFs[i]
        extensions[[length(extensions) + 1L]] <- list(
          S = new_S, log_w = new_log_w,
          log_BFs_path = c(path$log_BFs_path, log_BFs[i])
        )
      }
    }

    if (length(extensions) == 0L) break

    # Deduplicate: same S -> sum weights via log-sum-exp
    S_keys <- sapply(extensions, function(p) paste(p$S, collapse = ","))
    unique_keys <- unique(S_keys)
    deduped <- lapply(unique_keys, function(k) {
      ix <- which(S_keys == k)
      log_ws <- sapply(extensions[ix], function(p) p$log_w)
      M <- max(log_ws)
      list(
        S = extensions[[ix[1]]]$S,
        log_w = M + log(sum(exp(log_ws - M))),
        log_BFs_path = extensions[[ix[1]]]$log_BFs_path
      )
    })

    # Keep top-B by log_w
    log_ws <- sapply(deduped, function(p) p$log_w)
    keep_idx <- order(-log_ws)[seq_len(min(B, length(deduped)))]
    beam <- deduped[keep_idx]
  }

  # Marginal PIPs: q^marginal_j = sum_{p : j in p.S} exp(p.log_w) / sum exp(p.log_w)
  log_ws <- sapply(beam, function(p) p$log_w)
  M <- max(log_ws)
  w_norm <- exp(log_ws - M) / sum(exp(log_ws - M))
  q_marginal <- numeric(m)
  for (i in seq_along(beam)) {
    q_marginal[beam[[i]]$S] <- q_marginal[beam[[i]]$S] + w_norm[i]
  }
  q_marginal <- pmin(pmax(q_marginal, 0), 1)

  # Selected set via posterior threshold (default 0.5, applied to marginal PIPs)
  selected <- which(q_marginal > 0.5)
  K_hat <- length(selected)

  # Best-path indices (for backward compat with stepwise output)
  best_path <- beam[[which.max(log_ws)]]

  list(
    indices = selected,
    q_marginal = q_marginal,           # NEW: marginal PIP per SNP
    K_hat = K_hat,
    best_path_indices = best_path$S,
    best_path_log_BFs = best_path$log_BFs_path,
    beam_size = length(beam),
    variant = "path_averaged_marginal"
  )
}

# Joint variant: identical structure, uses lmm_joint_bf_marginalised_k_batched
PA_JS_L_LMM_stepwise_fast <- function(y, X, a_n = NULL,
                                          tau2 = 0.04,
                                          B = 10L,
                                          K_max = 10L,
                                          n_nodes = 15L,
                                          prior = "half_cauchy",
                                          prior_params = list(scale = 1.0),
                                          prior_inclusion = 0.5,
                                          verbose = FALSE) {
  n <- nrow(X); m <- ncol(X)
  if (is.null(a_n))    a_n <- rep(1, n)
  if (!is.matrix(a_n)) a_n <- matrix(a_n, ncol = 1L)

  beam <- list(list(S = integer(0), log_w = 0, log_BFs_path = numeric(0)))

  for (k_step in seq_len(K_max)) {
    if (verbose) cat(sprintf("  step %d: beam size = %d\n",
                                k_step, length(beam)))
    extensions <- list()
    for (path in beam) {
      F_k <- if (length(path$S) == 0L) a_n
             else cbind(a_n, X[, path$S, drop = FALSE])
      X_S <- if (length(path$S) == 0L) matrix(0, nrow = n, ncol = 0L)
             else X[, path$S, drop = FALSE]
      remaining <- setdiff(seq_len(m), path$S)
      if (length(remaining) == 0L) next
      rem_mat <- X[, remaining, drop = FALSE]
      R_step  <- tcrossprod(rem_mat) / max(ncol(rem_mat) - 1L, 1L)
      R_eig   <- eigendecompose_R(R_step)
      bf_all <- lmm_joint_bf_marginalised_k_batched(y, rem_mat, F_k, X_S, R_eig,
                                                        tau2 = tau2,
                                                        prior = prior,
                                                        prior_params = prior_params,
                                                        n_nodes = n_nodes)
      log_BFs <- bf_all$log_BF
      local_top <- order(-log_BFs)[seq_len(min(B, length(log_BFs)))]
      for (i in local_top) {
        j <- remaining[i]
        new_S <- sort(c(path$S, j))
        extensions[[length(extensions) + 1L]] <- list(
          S = new_S, log_w = path$log_w + log_BFs[i],
          log_BFs_path = c(path$log_BFs_path, log_BFs[i]))
      }
    }
    if (length(extensions) == 0L) break

    S_keys <- sapply(extensions, function(p) paste(p$S, collapse = ","))
    unique_keys <- unique(S_keys)
    deduped <- lapply(unique_keys, function(k) {
      ix <- which(S_keys == k)
      log_ws <- sapply(extensions[ix], function(p) p$log_w)
      M <- max(log_ws)
      list(S = extensions[[ix[1]]]$S,
            log_w = M + log(sum(exp(log_ws - M))),
            log_BFs_path = extensions[[ix[1]]]$log_BFs_path)
    })
    log_ws <- sapply(deduped, function(p) p$log_w)
    keep_idx <- order(-log_ws)[seq_len(min(B, length(deduped)))]
    beam <- deduped[keep_idx]
  }

  log_ws <- sapply(beam, function(p) p$log_w)
  M <- max(log_ws)
  w_norm <- exp(log_ws - M) / sum(exp(log_ws - M))
  q_marginal <- numeric(m)
  for (i in seq_along(beam)) {
    q_marginal[beam[[i]]$S] <- q_marginal[beam[[i]]$S] + w_norm[i]
  }
  q_marginal <- pmin(pmax(q_marginal, 0), 1)
  selected <- which(q_marginal > 0.5)
  K_hat <- length(selected)
  best_path <- beam[[which.max(log_ws)]]

  list(
    indices = selected,
    q_marginal = q_marginal,
    K_hat = K_hat,
    best_path_indices = best_path$S,
    best_path_log_BFs = best_path$log_BFs_path,
    beam_size = length(beam),
    variant = "path_averaged_joint"
  )
}
