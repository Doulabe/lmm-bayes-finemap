# ==============================================================================
#  17_rho_V_diagnostic.R
#  Empirical distribution of the multiple-correlation coefficient rho^V_jk(delta)
#  and of D_jk/n across replicates and candidate-step pairs.
#  Methodological choice: for the empirical distribution of rho^V we condition
#  the analysis on the TRUE causal set (rather than on a stepwise-selected
#  set), because (a) the distributional question only depends on Fk, not on
#  the path that produced Fk, and (b) this removes ~50x of compute (no
#  per-step REML refit; one eigen and one REML fit per replicate are enough).
#  Each step k corresponds to conditioning on F_k = [1, X_{truth[1:k-1]}].
#
#  Output:
#    results/bench_full/rho_V/rho_V_diagnostic.rds  -- raw per-(j,k) values
#
#  Run from project root:
#    Rscript sim/bench_full/17_rho_V_diagnostic.R --cores 6 --B 100
# ==============================================================================

suppressPackageStartupMessages({
  library(parallel)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# --- BLAS thread pinning (CRITICAL for mclapply efficiency) ------------------
Sys.setenv(OMP_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1",
           BLIS_NUM_THREADS = "1")
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
  cat("BLAS thread pinning applied via RhpcBLASctl.\n")
}

# --- Args --------------------------------------------------------------------
args  <- commandArgs(trailingOnly = TRUE)
.find <- function(flag, default) {
  i <- which(args == flag); if (length(i)) args[i + 1] else default
}
N_CORES <- as.integer(.find("--cores", 6L))
B       <- as.integer(.find("--B",     100L))

# --- Config ------------------------------------------------------------------
source("sim/bench_full/00_config.R")
if (!requireNamespace("rrBLUP", quietly = TRUE))
  stop("rrBLUP needed for REML; install via install.packages('rrBLUP')")

OUT_DIR <- "results/bench_full/rho_V"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
FIG_DIR <- "theory/overleaf_compact/figures"
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

GRID <- expand.grid(
  rho      = c(0.80, 0.95, 0.98),
  sigma_g2 = c(0.0, 0.5),
  stringsAsFactors = FALSE
)

# --- Worker: rho^V_jk and D_jk/n on the truth-conditioned path ---------------
#
# Algorithmic structure (per replicate):
#   1. Generate data ; truth indices known.
#   2. Eigendecompose full-X kernel K = X X^T / (m-1)  [ONCE per replicate]
#   3. Fit REML on the null model y ~ 1 + GRM           [ONCE per replicate]
#   4. For each step k = 1, ..., length(truth)+1:
#        F_k = [1, X[, truth[1:(k-1)]]]
#        Compute rho^V_jk and D_jk/n for ALL remaining candidates by 3
#        matrix operations (vectorised over candidates).
#
diag_one_replicate <- function(rep_id, rho, sigma_g2) {
  set.seed(20260601L + rep_id + 1000L * which(c(0.80, 0.95, 0.98) == rho)
                                 + 100L  * which(c(0.0, 0.5) == sigma_g2))
  d <- gen_dataset(n = ANCHOR$n, m = ANCHOR$m, rho = rho,
                   K_true = ANCHOR$K_true, beta_true = ANCHOR$beta_true,
                   sigma_g2 = sigma_g2, block_size = ANCHOR$block_size,
                   seed = 20260601L + rep_id)
  y <- d$y; X <- d$X; truth <- d$truth
  n <- nrow(X); m <- ncol(X)

  # Step 1: one eigendecomposition for the full-X kernel
  Kmat <- tcrossprod(X) / max(m - 1L, 1L)
  ee   <- eigen(Kmat, symmetric = TRUE)
  U <- ee$vectors; s <- pmax(ee$values, 0)

  # Step 2: one REML fit on the null model y ~ 1 + GRM
  fit <- tryCatch(rrBLUP::mixed.solve(y = y, K = Kmat),
                  error = function(e) NULL)
  delta_hat <- if (is.null(fit)) 1 else max(fit$Vu / fit$Ve, 1e-6)
  w <- 1 / (1 + delta_hat * s)

  # Step 3 (ONCE): precompute U^T X, Z^{-1} X, d_jj^V for ALL columns of X.
  # The per-step rho^V then reduces to small (p x m) ops + a subset mask.
  UtX     <- crossprod(U, X)                      # n x m (~8s once)
  Zinv_X  <- U %*% (w * UtX)                      # n x m (~8s once)
  d_jj_V_full <- colSums(X * Zinv_X)              # m  (~0.1s)

  # Step 4: per step, just subset & small linear algebra (no big matmul)
  res <- vector("list", length(truth) + 1L)
  for (k in seq_len(length(truth) + 1L)) {
    Sk_minus_1 <- if (k == 1L) integer(0) else truth[seq_len(k - 1L)]
    Tk_minus_1 <- setdiff(seq_len(m), Sk_minus_1)
    # Conditioning matrix F_k = [1, X[, S_{k-1}]] is p x n with p <= K+1
    if (length(Sk_minus_1) == 0L) {
      Fk      <- matrix(1, n, 1L)
      UtFk    <- matrix(colSums(U), ncol = 1L)   # = crossprod(U, 1)
      Zinv_Fk <- U %*% (w * UtFk)
    } else {
      Fk      <- cbind(1, X[, Sk_minus_1, drop = FALSE])
      UtFk    <- cbind(colSums(U), UtX[, Sk_minus_1, drop = FALSE])
      Zinv_Fk <- U %*% (w * UtFk)
    }
    FtZF     <- crossprod(Fk, Zinv_Fk)            # p x p
    FtZF_inv <- solve(FtZF + 1e-10 * diag(ncol(FtZF)))
    # C[, j] = F_k^T Z^{-1} X_j, computed over candidates via precomputed Zinv_X
    C_full <- crossprod(Fk, Zinv_X)               # p x m  (fast, small p)
    C_mat  <- C_full[, Tk_minus_1, drop = FALSE]
    AC     <- FtZF_inv %*% C_mat                  # p x m_k
    d_jj_V <- d_jj_V_full[Tk_minus_1]
    rho_V_raw <- colSums(C_mat * AC) / d_jj_V
    rho_V_vec <- pmin(pmax(rho_V_raw, 0), 1 - 1e-12)
    D_vec     <- d_jj_V * (1 - rho_V_vec) / n

    res[[k]] <- data.frame(rep_id = rep_id, k = k,
                           rho_V = rho_V_vec, D_over_n = D_vec)
  }
  do.call(rbind, res)
}

# --- Main loop: parallelise over (rho, sigma_g2, replicate) ------------------
cat(sprintf("Diagnostic: B=%d replicates x %d (rho, sigma_g2) cells = %d runs\n",
            B, nrow(GRID), B * nrow(GRID)))
all_rows <- list()
for (i in seq_len(nrow(GRID))) {
  rho_i   <- GRID$rho[i]
  sigg2_i <- GRID$sigma_g2[i]
  cat(sprintf(" Cell %d/%d : rho=%.2f sigma_g2=%.1f ...", i, nrow(GRID),
              rho_i, sigg2_i))
  flush.console()
  t0 <- Sys.time()
  rows <- mclapply(seq_len(B), function(b) {
    tryCatch(diag_one_replicate(b, rho = rho_i, sigma_g2 = sigg2_i),
             error = function(e) NULL)
  }, mc.cores = N_CORES)
  rows <- rows[!sapply(rows, is.null)]
  rows <- do.call(rbind, rows)
  rows$rho      <- rho_i
  rows$sigma_g2 <- sigg2_i
  all_rows[[i]] <- rows
  cat(sprintf(" done in %.1f s (%d rows)\n",
              as.numeric(Sys.time() - t0, units = "secs"), nrow(rows)))
  flush.console()
}
dat <- do.call(rbind, all_rows)
saveRDS(dat, file.path(OUT_DIR, "rho_V_diagnostic.rds"))
cat(sprintf("Saved raw data: %s  (n = %d rows)\n",
            file.path(OUT_DIR, "rho_V_diagnostic.rds"), nrow(dat)))

# --- Figure ------------------------------------------------------------------
dat$rho_label      <- factor(sprintf("rho == %.2f", dat$rho))
dat$sigma_g2_label <- factor(ifelse(dat$sigma_g2 == 0,
                                    "sigma[g]^2 == 0",
                                    "sigma[g]^2 == 0.5"))
long <- dat |>
  tidyr::pivot_longer(c(rho_V, D_over_n),
                      names_to = "quantity", values_to = "value") |>
  dplyr::mutate(quantity = factor(quantity, levels = c("rho_V", "D_over_n"),
                                  labels = c("rho[jk]^V(hat(delta))",
                                             "D[jk]/n")))

p <- ggplot(long, aes(x = value)) +
  geom_histogram(bins = 60, fill = "#1F77B4", colour = NA) +
  facet_grid(rho_label ~ quantity + sigma_g2_label,
             scales = "free", labeller = label_parsed) +
  theme_bw(base_size = 9) +
  labs(x = NULL, y = "count") +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
        panel.grid.minor = element_blank())

ggsave(file.path(FIG_DIR, "fig_rho_V_distribution.pdf"),
       p, width = 8, height = 6.5)
cat(sprintf("Figure written: %s\n",
            file.path(FIG_DIR, "fig_rho_V_distribution.pdf")))
cat("Done.\n")
