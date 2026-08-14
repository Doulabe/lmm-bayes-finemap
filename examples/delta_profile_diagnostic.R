# ==============================================================================
# delta_profile_diagnostic.R
#
# Quality-control diagnostic: the restricted profile log-likelihood of the
# variance-component ratio delta = sigma_g^2 / sigma_e^2 at a region's leading
# candidates.  This is the check behind the REML-BF / MBF decision rule
# (Supplementary Notes S8.4-S8.5, Figures S12-S13):
#
#   * sharply peaked profile  -> delta is well identified
#                                -> plug-in REML-BF is the efficient default,
#                                   MBF a conservative cross-check;
#   * flat or multimodal       -> delta is weakly identified
#                                -> prefer the marginalised MBF and read the
#                                   selection as operating-point dependent.
#
# Usage:
#   Rscript examples/delta_profile_diagnostic.R
# ==============================================================================

# --- A small synthetic region (replace X, y by your own candidate region) ----
set.seed(1)
n <- 600L; m <- 150L
X <- scale(matrix(rnorm(n * m), n, m))
truth <- c(20L, 70L, 120L)
y <- as.numeric(X[, truth] %*% c(0.8, 0.4, 0.4)) + rnorm(n)
# Optional fixed covariates (intercept is added internally):
W <- matrix(1, n, 1L)            # e.g. cbind(1, sex) for sex adjustment

# --- Restricted profile log-likelihood of delta at a candidate locus ----------
# Conditions on W and the candidate; the random-effect kernel is built from the
# remaining markers (leave-candidate-out background).
profile_loglik_delta <- function(y, X, W, cand, delta_grid) {
  n <- nrow(X); m <- ncol(X)
  F_cond <- cbind(W, X[, cand])
  rest   <- setdiff(seq_len(m), cand)
  K  <- tcrossprod(X[, rest]) / max(length(rest) - 1L, 1L)
  ee <- eigen(K, symmetric = TRUE); s <- pmax(ee$values, 0); U <- ee$vectors
  Uty <- as.numeric(crossprod(U, y)); UtF <- crossprod(U, F_cond); k <- ncol(F_cond)
  vapply(delta_grid, function(d) {
    dinv <- 1 / (1 + d * s)
    FtF  <- crossprod(UtF * dinv, UtF); Fty <- crossprod(UtF * dinv, Uty)
    beta <- solve(FtF, Fty)
    rss  <- sum(dinv * (Uty - as.numeric(UtF %*% beta))^2)
    -0.5 * (sum(log(1 + d * s)) +
              as.numeric(determinant(FtF, logarithm = TRUE)$modulus) +
              (n - k) * log(rss))
  }, numeric(1))
}

# --- Pick the strongest candidate (single-SNP |z| after a mixed-model fit) -----
# Here we simply use the marginal least-squares |z| as a cheap proxy.
z <- abs(colSums(X * as.numeric(y)) / sqrt(colSums(X^2)))
lead <- which.max(z)

delta_grid <- exp(seq(log(1e-3), log(1e3), length.out = 60))
ll  <- profile_loglik_delta(y, X, W, lead, delta_grid)
llc <- ll - max(ll)                       # centre at the maximum
dhat   <- delta_grid[which.max(ll)]
in_ci  <- which(llc > -2)                 # 2-log-likelihood interval
width  <- log10(delta_grid[max(in_ci)]) - log10(delta_grid[min(in_ci)])

cat(sprintf("Lead candidate: column %d\n", lead))
cat(sprintf("delta_hat = %.3f   2-log-likelihood width = %.2f decades   -> %s\n",
            dhat, width, if (width > 3) "FLAT: prefer MBF" else "PEAKED: REML-BF is fine"))

# --- Plot (Figure S12 style) --------------------------------------------------
pdf("delta_profile_diagnostic.pdf", width = 6, height = 4)
plot(log10(delta_grid), llc, type = "l", lwd = 2, col = "#1b6ca8",
     xlab = expression(log[10](delta)),
     ylab = "profile restricted log-likelihood (centred)")
abline(h = -2, lty = 2, col = "grey40")
legend("bottomright", bty = "n", lwd = 2, col = "#1b6ca8",
       legend = sprintf("lead candidate (delta_hat = %.2f)", dhat))
dev.off()
cat("Wrote delta_profile_diagnostic.pdf\n")
