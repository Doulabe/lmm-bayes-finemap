# ==============================================================================
#  15b_mouse_reml_bf_sensitivity.R
#  Plug-in REML (REML-BF) sensitivity scan on the mouse BMI panel.
#  Mirrors the 4 analyses of 12_mouse_autosomes_only.R and
#  15_mouse_bmi_with_sex.R but evaluates delta by plug-in REML
#  (no quadrature) under the same eBIC stopping rule.
#
#  Output:
#    results/bench_full/15b_mouse_reml_bf_sensitivity/
#       wg_no_sex.rds, wg_with_sex.rds, autosomes_no_sex.rds, autosomes_with_sex.rds
#    theory/overleaf_compact/tables/tab_bmi_reml_bf_sensitivity.tex
#
#  Run from project root:
#    Rscript sim/bench_full/15b_mouse_reml_bf_sensitivity.R
# ==============================================================================

# BLAS pinning
Sys.setenv(OMP_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           BLIS_NUM_THREADS = "1")
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
}

source("sim/bench_full/00_config.R")

OUT_DIR <- "results/bench_full/15b_mouse_reml_bf_sensitivity"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
TEX_OUT <- "theory/overleaf_compact/tables/tab_bmi_reml_bf_sensitivity.tex"

if (!requireNamespace("rrBLUP", quietly = TRUE))
  stop("rrBLUP needed for REML; install via install.packages('rrBLUP')")
if (!requireNamespace("BGLR", quietly = TRUE))
  stop("BGLR needed for mouse data; install via install.packages('BGLR')")
suppressPackageStartupMessages(library(BGLR))
data(mice)

#  Plug-in REML stepwise selection with eBIC stopping 
# Supports a generic conditioning matrix W (default intercept only).
# Mirrors the algebra of plug_in_reml_eBIC() in 16_ablation.R but takes W
# so that sex (or any fixed covariate) can be entered.
plug_in_reml_eBIC_W <- function(y, X, W = NULL, K_max = 30L, tau2 = 0.04) {
  n <- nrow(X); m <- ncol(X)
  if (is.null(W)) W <- matrix(1, n, 1L)
  if (!is.matrix(W)) W <- matrix(W, ncol = 1L)
  kW <- ncol(W)

  selected <- integer(0); scores <- numeric(0); log_BFs <- numeric(0)
  remaining <- seq_len(m)
  ebic_path <- numeric(K_max + 1L); ebic_path[1L] <- Inf

  F_k <- W
  for (k_step in seq_len(K_max)) {
    if (length(remaining) == 0L) break
    rem <- X[, remaining, drop = FALSE]
    K   <- tcrossprod(rem) / max(ncol(rem) - 1L, 1L)

    # Residualise y on F_k for the null-model REML fit
    proj <- solve(crossprod(F_k), crossprod(F_k, y))
    y_r  <- y - F_k %*% proj
    fit  <- tryCatch(rrBLUP::mixed.solve(y = y_r, K = K),
                     error = function(e) NULL)
    if (is.null(fit)) break
    delta_hat <- max(fit$Vu / fit$Ve, 1e-6)

    # Whitening at delta_hat
    ee <- eigen(K, symmetric = TRUE)
    s <- pmax(ee$values, 0); U <- ee$vectors
    w_white <- 1 / sqrt(1 + delta_hat * s)
    Uty   <- crossprod(U, y) * w_white
    UtFk  <- crossprod(U, F_k) * w_white
    UtX   <- crossprod(U, rem) * w_white

    # Project out F_k in the whitened metric
    A    <- crossprod(UtFk)
    A_inv<- solve(A)
    beta <- A_inv %*% crossprod(UtFk, Uty)
    ty   <- Uty - UtFk %*% beta
    proj_X <- UtFk %*% (A_inv %*% crossprod(UtFk, UtX))
    tXj    <- UtX - proj_X

    D_jk   <- colSums(tXj^2)
    u_j    <- as.numeric(crossprod(tXj, ty))
    rss0   <- as.numeric(crossprod(ty))
    v_j    <- tau2 * D_jk
    df     <- n - ncol(F_k)
    Q2_j   <- u_j^2 / (D_jk * rss0 / df)
    log_BF <- -0.5 * log1p(v_j) - 0.5 * df *
              log1p(-v_j * Q2_j / df / (1 + v_j))
    j_best <- which.max(log_BF)
    if (!is.finite(log_BF[j_best])) break

    # eBIC stopping (selected covariates penalized; W not penalized)
    sel_new <- c(selected, remaining[j_best])
    F_new   <- cbind(F_k, X[, remaining[j_best], drop = FALSE])
    proj2   <- solve(crossprod(F_new), crossprod(F_new, y))
    rss1    <- as.numeric(crossprod(y - F_new %*% proj2))
    ll_new  <- -0.5 * (n - ncol(F_new)) * log(rss1)
    ebic_new <- -2 * ll_new + length(sel_new) * (log(n) + 2 * log(m))
    if (k_step >= 2L && ebic_new > ebic_path[k_step]) break

    selected <- sel_new
    scores   <- c(scores, log_BF[j_best])
    log_BFs  <- c(log_BFs, log_BF[j_best])
    F_k      <- F_new
    remaining <- remaining[-j_best]
    ebic_path[k_step + 1L] <- ebic_new
  }
  list(method = "REML_BF_eBIC",
       indices = selected, K_hat = length(selected),
       scores = scores, log_BFs = log_BFs)
}

#  Mouse data preparation 
y_full   <- as.numeric(mice.pheno$Obesity.BMI)
sex_full <- as.integer(mice.pheno$GENDER == "M")
keep <- !is.na(y_full) & !is.na(sex_full)
y    <- y_full[keep]
sex  <- sex_full[keep]
X    <- scale(mice.X[keep, ])
map  <- mice.map
n <- length(y); m_wg <- ncol(X)
cat(sprintf("Mouse data after filtering: n=%d, m_WG=%d\n", n, m_wg))

# Autosomes-only mask (drop chromosome X)
chr_vec  <- as.character(map$chr)
auto_idx <- which(chr_vec != "X")
X_auto   <- X[, auto_idx, drop = FALSE]
map_auto <- map[auto_idx, , drop = FALSE]
cat(sprintf("Autosomes-only: m_auto=%d (dropped %d X-linked SNPs)\n",
            ncol(X_auto), m_wg - ncol(X_auto)))

W_intercept <- matrix(1, n, 1L)
W_with_sex  <- cbind(intercept = rep(1, n), sex = sex - mean(sex))

#  Run the 4 analyses 
analyses <- list(
  list(tag = "wg_no_sex",       label = "WG, no sex",         X = X,      map = map,      W = W_intercept),
  list(tag = "wg_with_sex",     label = "WG, with sex",       X = X,      map = map,      W = W_with_sex),
  list(tag = "autosomes_no_sex",label = "Autosomes, no sex",  X = X_auto, map = map_auto, W = W_intercept),
  list(tag = "autosomes_with_sex",label = "Autosomes, with sex",X = X_auto, map = map_auto, W = W_with_sex)
)

summarise_chr <- function(map, indices) {
  if (length(indices) == 0L) return("(none)")
  chrs <- as.character(map$chr[indices])
  counts <- sort(table(chrs), decreasing = TRUE)
  paste(sprintf("%s:%d", names(counts), counts), collapse = ", ")
}

results <- list()
for (a in analyses) {
  cat(sprintf("\n=== [%s] starting REML-BF + eBIC scan ===\n", a$label))
  t0 <- Sys.time()
  fit <- plug_in_reml_eBIC_W(y, a$X, W = a$W, K_max = 30L, tau2 = 0.04)
  elapsed <- as.numeric(Sys.time() - t0, units = "secs")
  cat(sprintf(" -> K_hat = %d, elapsed = %.1f s\n", fit$K_hat, elapsed))
  chr_summary <- summarise_chr(a$map, fit$indices)
  cat(sprintf("    chromosomes: %s\n", chr_summary))

  results[[a$tag]] <- list(label = a$label, K_hat = fit$K_hat,
                            indices = fit$indices,
                            chr_summary = chr_summary,
                            elapsed = elapsed)
  saveRDS(results[[a$tag]],
          file.path(OUT_DIR, sprintf("%s.rds", a$tag)))
}

#  LaTeX table 
sink(TEX_OUT)
cat("% Auto-generated by sim/bench_full/15b_mouse_reml_bf_sensitivity.R\n")
cat("% Plug-in REML (REML-BF) sensitivity scan on the mouse BMI panel.\n")
cat("\\begin{tabular}{lrll}\n\\toprule\n")
cat("Analysis & $\\hat K$ & \\# chr.\\ & Target chromosomes (chr:count) \\\\\n")
cat("\\midrule\n")
for (tag in c("wg_no_sex", "wg_with_sex",
              "autosomes_no_sex", "autosomes_with_sex")) {
  r <- results[[tag]]
  n_chr <- if (r$K_hat == 0L) 0L
           else length(unique(as.character(
                  analyses[[which(sapply(analyses, function(a) a$tag) == tag)]]$map$chr[r$indices])))
  cat(sprintf("%s & $%d$ & $%d$ & %s \\\\\n",
              r$label, r$K_hat, n_chr, r$chr_summary))
}
cat("\\bottomrule\n\\end{tabular}\n")
sink()
cat(sprintf("\nWrote %s\n", TEX_OUT))
cat("Done.\n")
