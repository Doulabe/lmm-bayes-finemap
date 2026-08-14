# ==============================================================================
# 23_delta_profile_mouse.R
#
# Identifiability of the variance-component ratio delta at the sex-adjusted
# competing loci of the mouse BMI panel (Supplementary Note S8.4, Figure on the
# delta profile).
#
# The sex-adjusted whole-genome lead shifts from chromosome X (under MBF) to
# chromosome 15 (under plug-in REML-BF). This script plots the restricted
# profile log-likelihood of delta at the leading chr-X and chr-15 markers,
# conditioning on [intercept, sex] and the candidate, with the random-effect
# kernel built from the remaining markers. A sharply peaked profile means delta
# is well identified and the divergence is not an identifiability artefact.
#
# Output: results/bench_full/delta_profile_mouse.rds
#         results/bench_full/fig_delta_profile.pdf
#
# Usage:
#   Rscript sim/bench_full/23_delta_profile_mouse.R
# ==============================================================================

Sys.setenv(OPENBLAS_NUM_THREADS = "4", OMP_NUM_THREADS = "4")
suppressPackageStartupMessages({ library(BGLR); library(rrBLUP) })

OUT_DIR <- "results/bench_full"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

data(mice)
X <- mice.X; y <- mice.pheno$Obesity.BMI
sex <- as.integer(mice.pheno$GENDER == "M")
X <- scale(X, center = TRUE, scale = TRUE); X[!is.finite(X)] <- 0
chr <- as.character(mice.map$chr); snp_id <- mice.map$snp_id
n <- nrow(X); m <- ncol(X)

# Residualise y on (intercept, sex).
Wc <- cbind(1, sex - mean(sex))
y_resid <- as.numeric(y - Wc %*% solve(crossprod(Wc), crossprod(Wc, y)))

# Profile restricted log-likelihood of delta = Vu/Ve at a candidate locus,
# conditioning on W=[1,sex] + the candidate, kernel from the rest.
profile_loglik_delta <- function(cand_idx, delta_grid) {
  F_cond <- cbind(Wc, X[, cand_idx])
  rest <- setdiff(seq_len(m), cand_idx)
  K <- tcrossprod(X[, rest]) / max(length(rest) - 1L, 1L)
  ee <- eigen(K, symmetric = TRUE); s <- pmax(ee$values, 0); U <- ee$vectors
  Uty <- as.numeric(crossprod(U, y)); UtF <- crossprod(U, F_cond)
  k <- ncol(F_cond)
  ll <- numeric(length(delta_grid))
  for (i in seq_along(delta_grid)) {
    d <- delta_grid[i]; dinv <- 1 / (1 + d * s)
    FtF <- crossprod(UtF * dinv, UtF); Fty <- crossprod(UtF * dinv, Uty)
    beta <- solve(FtF, Fty)
    resid <- Uty - as.numeric(UtF %*% beta)
    rss <- sum(dinv * resid^2)
    ll[i] <- -0.5 * (sum(log(1 + d * s)) +
                       as.numeric(determinant(FtF, logarithm = TRUE)$modulus) +
                       (n - k) * log(rss))
  }
  ll
}

# Single-SNP sex-adjusted |Z| to pick the leading SNP on each chromosome.
single_Z <- function(yv) {
  K <- tcrossprod(X) / max(m - 1L, 1L); fit <- rrBLUP::mixed.solve(y = yv, K = K)
  ee <- eigen(K, symmetric = TRUE); dh <- fit$Vu / fit$Ve
  di <- 1 / sqrt(1 + dh * pmax(ee$values, 0))
  UtX <- crossprod(ee$vectors, X); Uty <- crossprod(ee$vectors, yv)
  as.numeric(crossprod(UtX * di, Uty * di) / sqrt(colSums((UtX * di)^2) + 1e-300))
}
Z <- abs(single_Z(y_resid))
chrX_idx  <- which(chr == "X")[which.max(Z[chr == "X"])]
chr15_idx <- which(chr == "15")[which.max(Z[chr == "15"])]
cat(sprintf("ChrX  lead: %s (|Z|=%.2f)\n", snp_id[chrX_idx], Z[chrX_idx]))
cat(sprintf("Chr15 lead: %s (|Z|=%.2f)\n", snp_id[chr15_idx], Z[chr15_idx]))

delta_grid <- exp(seq(log(1e-3), log(1e3), length.out = 60))
llX  <- profile_loglik_delta(chrX_idx, delta_grid)
ll15 <- profile_loglik_delta(chr15_idx, delta_grid)

summ <- function(ll, lab) {
  llc <- ll - max(ll); imax <- which.max(ll); dhat <- delta_grid[imax]
  in_ci <- which(llc > -2)
  width <- log10(delta_grid[max(in_ci)]) - log10(delta_grid[min(in_ci)])
  cat(sprintf("%s: delta_hat=%.3f, log10-CI-width(2LL)=%.2f decades, flat=%s\n",
              lab, dhat, width, ifelse(width > 3, "YES (near-flat)", "no")))
}
summ(llX, "ChrX ")
summ(ll15, "Chr15")

saveRDS(list(delta = delta_grid, llX = llX - max(llX), ll15 = ll15 - max(ll15),
             snpX = snp_id[chrX_idx], snp15 = snp_id[chr15_idx]),
        file.path(OUT_DIR, "delta_profile_mouse.rds"))

# ---- Figure: centred profile log-likelihoods with the 2-LL reference line ----
fig <- file.path(OUT_DIR, "fig_delta_profile.pdf")
pdf(fig, width = 6.5, height = 4.2)
op <- par(mar = c(4.2, 4.2, 1.0, 1.0))
llXc <- llX - max(llX); ll15c <- ll15 - max(ll15)
yl <- range(c(llXc, ll15c, -6))
plot(log10(delta_grid), llXc, type = "l", lwd = 2, col = "#1b6ca8",
     ylim = yl, xlab = expression(log[10](delta)),
     ylab = "profile restricted log-likelihood (centred)")
lines(log10(delta_grid), ll15c, lwd = 2, col = "#c8401b", lty = 1)
abline(h = -2, lty = 2, col = "grey40")
legend("bottomright", bty = "n", lwd = 2, col = c("#1b6ca8", "#c8401b"),
       legend = c(sprintf("chr X (%s)", snp_id[chrX_idx]),
                  sprintf("chr 15 (%s)", snp_id[chr15_idx])))
par(op); dev.off()

cat(sprintf("\nSaved: %s\n       %s\n",
            file.path(OUT_DIR, "delta_profile_mouse.rds"), fig))
