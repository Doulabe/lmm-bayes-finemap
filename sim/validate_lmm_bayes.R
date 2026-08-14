# ==============================================================================
# validate_lmm_bayes.R
# Numerical validation of the Bayesian Gaussian LMM framework.
# Runs the checks V1-V6 from Section "Validation programme" of
# lmm_bayes.tex (companion theory document).
#
# V1: BF at delta=0 matches Wakefield ABF (closed-form check)
# V2: Plug-in REML vs Bayesian agreement at large n (Prop. asymptotic)
# V3: Type-I error control under the global null
# V4: Quadrature convergence in N_delta
# V5: Stepwise consistency (added when LMM_stepwise.R is sourced)
# V6: Joint vs marginal advantage under correlation (added when stepwise sourced)
#
# Usage:  Rscript sim/validate_lmm_bayes.R   (from project root)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
})

# Resolve LMM_core.R location: try same dir, then R/, then ../R/
.try_source <- function(paths) {
  for (p in paths) if (file.exists(p)) { source(p); return(invisible(p)) }
  stop("LMM_core.R not found. Run from project root with R/LMM_core.R available.")
}
.try_source(c("R/LMM_core.R", "../R/LMM_core.R", "LMM_core.R"))

# Helper: compute Wakefield ABF directly on raw y, X_j (delta=0 case)
wakefield_abf_raw <- function(y, X_j, a_n, tau2) {
  n   <- length(y)
  ata <- as.numeric(crossprod(a_n))
  beta0 <- as.numeric(crossprod(a_n, y)) / ata
  res <- y - a_n * beta0
  rss0 <- as.numeric(crossprod(res))
  Xj_y <- as.numeric(crossprod(X_j, res))
  XjtXj <- as.numeric(crossprod(X_j))
  v <- tau2 * XjtXj
  rss1 <- rss0 - v * Xj_y^2 / (XjtXj * (1 + v))
  -0.5 * log1p(v) + 0.5 * (n - 1) * (log(rss0) - log(rss1))
}

# Generate a single dataset under the LMM at k=1
gen_data <- function(n, m, sigma_g2, beta_1 = 0, intercept = 0,
                      seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- matrix(rnorm(n * m), n, m)
  X <- scale(X)
  # Polygenic random effect: u = (1/sqrt(m)) X[,-j] beta_others, beta ~ N(0, sigma_g2)
  # Approximate by drawing u ~ N(0, sigma_g2 * R) directly
  R <- tcrossprod(X[, -1L]) / (m - 1L)
  ee <- eigen(R, symmetric = TRUE)
  half <- ee$vectors %*% diag(sqrt(pmax(sigma_g2 * ee$values, 0))) %*% t(ee$vectors)
  u <- as.numeric(half %*% rnorm(n))
  eps <- rnorm(n)
  y <- intercept + X[, 1L] * beta_1 + u + eps
  list(y = y, X = X, R = R, R_eig = list(U = ee$vectors,
                                          s = pmax(ee$values, 1e-10)))
}


# ==============================================================================
# V1 — Closed-form check at delta = 0
# ==============================================================================
v1_test <- function(seed = 42L) {
  cat("\n=== V1: BF at delta=0 should match plain Wakefield ABF ===\n")
  d <- gen_data(n = 200L, m = 50L, sigma_g2 = 0, beta_1 = 0, seed = seed)
  X_j <- d$X[, 1L]; X_others <- d$X[, -1L]
  R_eig <- list(U = d$R_eig$U, s = d$R_eig$s)
  a_n <- rep(1, length(d$y))

  # Our Bayesian BF at delta = 0
  bf_lmm <- lmm_bf_at_delta(d$y, X_j, a_n, delta = 0, R_eig, tau2 = 0.04)$log_BF
  # Plain Wakefield ABF on raw y
  bf_wake <- wakefield_abf_raw(d$y, X_j, a_n, tau2 = 0.04)
  cat(sprintf("  log-BF (LMM at delta=0) = %.6f\n",  bf_lmm))
  cat(sprintf("  log-BF (Wakefield raw)  = %.6f\n",  bf_wake))
  cat(sprintf("  Difference              = %.2e  (target: < 1e-10)\n",
              abs(bf_lmm - bf_wake)))
  invisible(abs(bf_lmm - bf_wake))
}


# ==============================================================================
# V2 — Plug-in REML vs Bayesian at large n
# ==============================================================================
v2_test <- function(B = 100L, n = 1000L, m = 200L, sigma_g2 = 1.0,
                     tau2 = 0.04, seed = 123L) {
  cat(sprintf("\n=== V2: Plug-in vs Bayesian at n=%d (B=%d reps) ===\n", n, B))
  diffs <- numeric(B); plug <- numeric(B); bayes <- numeric(B); deltahat <- numeric(B)
  for (b in seq_len(B)) {
    d <- gen_data(n, m, sigma_g2, beta_1 = 0, seed = seed + b)
    X_j <- d$X[, 1L]; a_n <- rep(1, n)
    R_eig <- list(U = d$R_eig$U, s = d$R_eig$s)

    delta_hat <- lmm_reml_delta(d$y, X_j, a_n, R_eig)
    bf_plug <- lmm_bf_at_delta(d$y, X_j, a_n, delta_hat, R_eig, tau2)$log_BF
    bf_bay  <- lmm_bf_marginalised(d$y, X_j, a_n, R_eig, tau2,
                                     n_nodes = 30L)$log_BF
    plug[b]  <- bf_plug
    bayes[b] <- bf_bay
    diffs[b] <- bf_bay - bf_plug
    deltahat[b] <- delta_hat
  }
  cat(sprintf("  Mean |log-BF_Bayes - log-BF_Plug| = %.4f\n", mean(abs(diffs))))
  cat(sprintf("  Median diff                       = %.4f\n", median(diffs)))
  cat(sprintf("  REML delta_hat:  mean=%.3f  sd=%.3f  (true=%.2f)\n",
              mean(deltahat), sd(deltahat), sigma_g2))
  cat(sprintf("  Plug-in  log-BF  mean=%.3f  sd=%.3f\n", mean(plug), sd(plug)))
  cat(sprintf("  Bayesian log-BF  mean=%.3f  sd=%.3f\n", mean(bayes), sd(bayes)))
  invisible(list(diffs = diffs, plug = plug, bayes = bayes,
                 delta_hat = deltahat))
}


# ==============================================================================
# V3 — Bayesian advantage at small n (heavy structure)
# ==============================================================================
v3_test <- function(B = 100L, n = 100L, m = 200L, sigma_g2 = 1.5,
                     tau2 = 0.04, seed = 456L) {
  cat(sprintf("\n=== V3: Type-I error at small n=%d (B=%d reps under H0) ===\n",
              n, B))
  fp_plug <- 0L; fp_bay <- 0L
  threshold <- 3                              # log-BF cutoff (Kass-Raftery)
  for (b in seq_len(B)) {
    d <- gen_data(n, m, sigma_g2, beta_1 = 0, seed = seed + b)
    X_j <- d$X[, 1L]; a_n <- rep(1, n)
    R_eig <- list(U = d$R_eig$U, s = d$R_eig$s)
    delta_hat <- lmm_reml_delta(d$y, X_j, a_n, R_eig)
    bf_plug <- lmm_bf_at_delta(d$y, X_j, a_n, delta_hat, R_eig, tau2)$log_BF
    bf_bay  <- lmm_bf_marginalised(d$y, X_j, a_n, R_eig, tau2,
                                     n_nodes = 30L)$log_BF
    if (bf_plug > threshold) fp_plug <- fp_plug + 1L
    if (bf_bay  > threshold) fp_bay  <- fp_bay  + 1L
  }
  cat(sprintf("  Plug-in REML  type-I (log-BF > %.0f): %d / %d = %.2f\n",
              threshold, fp_plug, B, fp_plug / B))
  cat(sprintf("  Bayesian      type-I (log-BF > %.0f): %d / %d = %.2f\n",
              threshold, fp_bay, B, fp_bay / B))
  invisible(list(fp_plug = fp_plug, fp_bay = fp_bay))
}


# ==============================================================================
# V4 — Quadrature convergence
# ==============================================================================
v4_test <- function(seed = 789L) {
  cat("\n=== V4: Quadrature convergence in N_delta ===\n")
  d <- gen_data(n = 500L, m = 200L, sigma_g2 = 0.5, beta_1 = 0.3,
                seed = seed)
  X_j <- d$X[, 1L]; a_n <- rep(1, length(d$y))
  R_eig <- list(U = d$R_eig$U, s = d$R_eig$s)
  ref <- lmm_bf_marginalised(d$y, X_j, a_n, R_eig, tau2 = 0.04,
                              n_nodes = 100L)$log_BF
  for (N in c(5L, 10L, 20L, 30L, 50L)) {
    bf <- lmm_bf_marginalised(d$y, X_j, a_n, R_eig, tau2 = 0.04,
                                n_nodes = N)$log_BF
    cat(sprintf("  N=%3d   log-BF = %.6f   |diff| from N=100 = %.2e\n",
                N, bf, abs(bf - ref)))
  }
  invisible(NULL)
}


# ==============================================================================
# V5 — Stepwise consistency
# Inject K_true causals at known effect sizes; run JS_L_LMM_stepwise with
# JointPosterior stopping; verify selection of all true causals.
# ==============================================================================
v5_test <- function(B = 30L, n = 300L, m = 100L, sigma_g2 = 0.3,
                     beta_true = c(0.6, 0.5, 0.4), tau2 = 0.04, theta = 0.9,
                     seed = 1234L) {
  K_true <- length(beta_true)
  cat(sprintf("\n=== V5: Stepwise consistency (B=%d, K_true=%d, theta=%.2f) ===\n",
              B, K_true, theta))
  if (!exists("JS_L_LMM_stepwise_fast") || !exists("MS_L_LMM_stepwise_fast")) {
    .try_source(c("R/LMM_stepwise_fast.R", "../R/LMM_stepwise_fast.R"))
  }
  JS_L_LMM_stepwise <- JS_L_LMM_stepwise_fast
  MS_L_LMM_stepwise <- MS_L_LMM_stepwise_fast

  perfect_js <- 0L; perfect_ms <- 0L
  recall_js <- numeric(B); recall_ms <- numeric(B)
  K_hat_js  <- integer(B); K_hat_ms  <- integer(B)

  for (b in seq_len(B)) {
    set.seed(seed + b)
    X <- matrix(rnorm(n*m), n, m); X <- scale(X)
    R <- tcrossprod(X[, -seq_len(K_true)]) / (m - K_true)
    ee <- eigen(R, symmetric = TRUE)
    half <- ee$vectors %*% diag(sqrt(pmax(sigma_g2 * ee$values, 0))) %*% t(ee$vectors)
    u <- as.numeric(half %*% rnorm(n))
    y <- as.numeric(X[, seq_len(K_true)] %*% beta_true) + u + rnorm(n)

    res_js <- JS_L_LMM_stepwise(y, X, tau2 = tau2, K_max = K_true + 2L,
                                  criterion = "JointPosterior", theta = theta,
                                  n_nodes = 15L)
    res_ms <- MS_L_LMM_stepwise(y, X, tau2 = tau2, K_max = K_true + 2L,
                                  criterion = "JointPosterior", theta = theta,
                                  n_nodes = 15L)

    truth <- seq_len(K_true)
    tp_js <- length(intersect(res_js$indices, truth))
    tp_ms <- length(intersect(res_ms$indices, truth))
    recall_js[b] <- tp_js / K_true
    recall_ms[b] <- tp_ms / K_true
    K_hat_js[b]  <- res_js$K_hat
    K_hat_ms[b]  <- res_ms$K_hat
    if (tp_js == K_true && res_js$K_hat == K_true) perfect_js <- perfect_js + 1L
    if (tp_ms == K_true && res_ms$K_hat == K_true) perfect_ms <- perfect_ms + 1L
  }
  cat(sprintf("  JS_L: perfect recovery %d/%d (%.0f%%), mean recall=%.3f, mean K_hat=%.2f\n",
              perfect_js, B, 100*perfect_js/B, mean(recall_js), mean(K_hat_js)))
  cat(sprintf("  MS_L: perfect recovery %d/%d (%.0f%%), mean recall=%.3f, mean K_hat=%.2f\n",
              perfect_ms, B, 100*perfect_ms/B, mean(recall_ms), mean(K_hat_ms)))
  invisible(list(perfect_js = perfect_js, perfect_ms = perfect_ms,
                 recall_js = recall_js, recall_ms = recall_ms))
}


# ==============================================================================
# V6 — Joint vs Marginal under correlated predictors (block-AR(1) LD)
# Generate correlated X via block-AR(1); inject correlated causals; compare
# JS_L vs MS_L on F1 and recall.
# ==============================================================================
v6_test <- function(B = 30L, n = 300L, m = 100L, block_size = 10L, r = 0.6,
                     sigma_g2 = 0.3, beta_true = c(0.6, 0.5, 0.4),
                     tau2 = 0.04, theta = 0.9, seed = 5678L) {
  K_true <- length(beta_true)
  cat(sprintf("\n=== V6: Joint vs marginal under block-AR(1) LD (B=%d, r=%.1f) ===\n",
              B, r))
  if (!exists("JS_L_LMM_stepwise_fast") || !exists("MS_L_LMM_stepwise_fast")) {
    .try_source(c("R/LMM_stepwise_fast.R", "../R/LMM_stepwise_fast.R"))
  }
  JS_L_LMM_stepwise <- JS_L_LMM_stepwise_fast
  MS_L_LMM_stepwise <- MS_L_LMM_stepwise_fast
  # Block-AR(1) covariance helper
  ar1_block <- function(p, r) {
    M <- matrix(0, p, p)
    for (i in seq_len(p)) for (j in seq_len(p)) M[i,j] <- r^abs(i-j)
    M
  }
  S_block <- ar1_block(block_size, r)
  L_block <- chol(S_block + 1e-8 * diag(block_size))

  f1_js <- numeric(B); f1_ms <- numeric(B)
  for (b in seq_len(B)) {
    set.seed(seed + b)
    Z <- matrix(rnorm(n*m), n, m)
    # Block-correlate
    n_blocks <- m %/% block_size
    for (bk in seq_len(n_blocks)) {
      cols <- ((bk-1)*block_size + 1):(bk*block_size)
      Z[, cols] <- Z[, cols] %*% L_block
    }
    X <- scale(Z)
    R <- tcrossprod(X[, -seq_len(K_true)]) / (m - K_true)
    ee <- eigen(R, symmetric = TRUE)
    half <- ee$vectors %*% diag(sqrt(pmax(sigma_g2 * ee$values, 0))) %*% t(ee$vectors)
    u <- as.numeric(half %*% rnorm(n))
    y <- as.numeric(X[, seq_len(K_true)] %*% beta_true) + u + rnorm(n)

    res_js <- JS_L_LMM_stepwise(y, X, tau2 = tau2, K_max = K_true + 2L,
                                  criterion = "JointPosterior", theta = theta,
                                  n_nodes = 15L)
    res_ms <- MS_L_LMM_stepwise(y, X, tau2 = tau2, K_max = K_true + 2L,
                                  criterion = "JointPosterior", theta = theta,
                                  n_nodes = 15L)

    truth <- seq_len(K_true)
    f1 <- function(sel, K_hat) {
      tp <- length(intersect(sel, truth))
      if (K_hat == 0L || K_true == 0L) return(0)
      prec <- tp / K_hat; rec <- tp / K_true
      if (prec + rec == 0) 0 else 2*prec*rec/(prec + rec)
    }
    f1_js[b] <- f1(res_js$indices, res_js$K_hat)
    f1_ms[b] <- f1(res_ms$indices, res_ms$K_hat)
  }
  win_js <- sum(f1_js > f1_ms); tie <- sum(f1_js == f1_ms)
  cat(sprintf("  JS_L mean F1 = %.3f  (sd %.3f)\n", mean(f1_js), sd(f1_js)))
  cat(sprintf("  MS_L mean F1 = %.3f  (sd %.3f)\n", mean(f1_ms), sd(f1_ms)))
  cat(sprintf("  JS_L wins %d/%d, ties %d, MS_L wins %d/%d (under block-AR(1) r=%.1f)\n",
              win_js, B, tie, B - win_js - tie, B, r))
  invisible(list(f1_js = f1_js, f1_ms = f1_ms))
}


# ==============================================================================
# Main entry point
# ==============================================================================
if (sys.nframe() == 0L) {
  message("=== Validation: Bayesian Gaussian LMM with variance-component marginalisation ===")
  v1 <- v1_test()
  v2 <- v2_test(B = 50L)
  v3 <- v3_test(B = 50L)
  v4 <- v4_test()

  # Need fast stepwise wrappers for V5/V6 (uses exact joint BF via SMW)
  .try_source(c("R/LMM_stepwise_fast.R", "../R/LMM_stepwise_fast.R"))

  v5 <- v5_test(B = 20L)
  v6 <- v6_test(B = 20L)
  message("\nValidation complete.")
}
