# ==============================================================================
#  16_ablation.R
#  Ablation table decomposing the contributions of the proposed framework:
#  (A) plug-in REML + eBIC         -- delta plugged in
#  (B) marginalised delta + eBIC   -- MS_L_eBIC
#  (C) joint Schur + eBIC          -- JS_L_eBIC
#  (D) Bayesian score + JP99       -- JS_L_JP99
#  Anchor: n=1000, m=5000, rho=0.95, K_true=5, sigma_g2 in {0, 0.5}
#
#  Output:
#    results/bench_full/16_ablation/<cell>_b<rep>.rds
#    theory/overleaf_compact/tables/tab_ablation.tex
#
#  Run from project root:
#    Rscript sim/bench_full/16_ablation.R --cores 8 --B 100
# ==============================================================================
suppressPackageStartupMessages({
  library(parallel)
})

# --- BLAS thread pinning (CRITICAL for mclapply efficiency) ------------------
# Without this, each mclapply worker spawns its own BLAS thread pool, causing
# severe contention on multi-core runs (observed: ~10x slowdown on 6 cores).
# Set both via env vars (for any future child R) and via R-level setters.
Sys.setenv(OMP_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1",
           BLIS_NUM_THREADS = "1")
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
  cat("BLAS thread pinning applied via RhpcBLASctl.\n")
} else {
  message("RhpcBLASctl not installed -- relying on env vars only.\n",
          "  install.packages('RhpcBLASctl') for guaranteed pinning.")
}

source("sim/bench_full/00_config.R")

args <- commandArgs(trailingOnly = TRUE)
.arg <- function(flag, def) {
  i <- which(args == flag); if (length(i)) args[i + 1] else def
}
N_CORES <- as.integer(.arg("--cores", 4L))
B       <- as.integer(.arg("--B",     100L))

OUT_DIR <- "results/bench_full/16_ablation"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Variant A: plug-in REML + eBIC 
# Implementation: fit REML at each step on the leave-out kernel of the null
# model, get delta_hat, evaluate the closed-form BF at that delta only,
# stop by extended BIC.
plug_in_reml_eBIC <- function(y, X, K_max = 10L, tau2 = 0.04) {
  if (!requireNamespace("rrBLUP", quietly = TRUE)) {
    stop("Variant A requires rrBLUP::mixed.solve")
  }
  n <- nrow(X); m <- ncol(X)
  selected <- integer(0); scores <- numeric(0); log_BFs <- numeric(0)
  remaining <- seq_len(m)
  ll_path <- numeric(K_max + 1L); ebic_path <- numeric(K_max + 1L)
  ebic_path[1L] <- Inf

  F_k <- matrix(1, n, 1L)
  for (k_step in seq_len(K_max)) {
    if (length(remaining) == 0L) break
    rem <- X[, remaining, drop = FALSE]
    K   <- tcrossprod(rem) / max(ncol(rem) - 1L, 1L)

    # Residualize y on F_k for the null-model REML fit
    proj <- solve(crossprod(F_k), crossprod(F_k, y))
    y_r  <- y - F_k %*% proj
    fit  <- tryCatch(rrBLUP::mixed.solve(y = y_r, K = K),
                     error = function(e) NULL)
    if (is.null(fit)) break
    delta_hat <- fit$Vu / fit$Ve

    # Whitening at delta_hat
    ee <- eigen(K, symmetric = TRUE)
    s <- pmax(ee$values, 0); U <- ee$vectors
    w <- 1 / sqrt(1 + delta_hat * s)
    Uty   <- crossprod(U, y) * w
    UtFk  <- crossprod(U, F_k) * w
    UtX   <- crossprod(U, rem) * w

    # Project out F_k in the whitened metric: residuals tilde_y and tilde_Xj
    A    <- crossprod(UtFk)
    A_inv<- solve(A)
    beta <- A_inv %*% crossprod(UtFk, Uty)
    ty   <- Uty - UtFk %*% beta
    proj_X <- UtFk %*% (A_inv %*% crossprod(UtFk, UtX))
    tXj    <- UtX - proj_X

    D_jk   <- colSums(tXj^2)
    u_j    <- as.numeric(crossprod(tXj, ty))
    rss0   <- as.numeric(crossprod(ty))
    # Closed-form conditional BF (eq. for MS_L)
    v_j    <- tau2 * D_jk
    Q2_j   <- u_j^2 / (D_jk * rss0 / (n - ncol(F_k)))
    df     <- n - ncol(F_k)
    log_BF <- -0.5 * log1p(v_j) - 0.5 * df *
              log1p(-v_j * Q2_j / df / (1 + v_j))
    j_best <- which.max(log_BF)
    if (!is.finite(log_BF[j_best])) break

    # eBIC stopping: -2 * profile_loglik + k * (log n + 2 log m)
    # profile loglik at the model with k SNPs included: -0.5 * df * log(RSS_min)
    sel_new <- c(selected, remaining[j_best])
    F_new   <- cbind(F_k, X[, remaining[j_best], drop = FALSE])
    # refit on new F to get profile loglik
    proj2  <- solve(crossprod(F_new), crossprod(F_new, y))
    rss1   <- as.numeric(crossprod(y - F_new %*% proj2))
    ll_new <- -0.5 * (n - ncol(F_new)) * log(rss1)
    ebic_new <- -2 * ll_new + length(sel_new) * (log(n) + 2 * log(m))
    if (k_step >= 2L && ebic_new > ebic_path[k_step]) break

    selected <- sel_new
    scores   <- c(scores, log_BF[j_best])
    log_BFs  <- c(log_BFs, log_BF[j_best])
    F_k      <- F_new
    remaining <- remaining[-j_best]
    ebic_path[k_step + 1L] <- ebic_new
  }
  list(method = "plug_in_REML_eBIC",
       indices = selected, K_hat = length(selected),
       scores = scores, log_BFs = log_BFs)
}

# --- Run one (cell, rep) ------------------------------------------------------
run_ablation_rep <- function(rep_id, sigma_g2) {
  set.seed(20260602L + rep_id + (sigma_g2 > 0) * 100000L)
  d <- gen_dataset(n = ANCHOR$n, m = ANCHOR$m, rho = ANCHOR$rho,
                   K_true = ANCHOR$K_true, beta_true = ANCHOR$beta_true,
                   sigma_g2 = sigma_g2, block_size = ANCHOR$block_size,
                   seed = 20260602L + rep_id)
  y <- d$y; X <- d$X; truth <- d$truth

  cell_tag <- sprintf("anchor_sg%.1f", sigma_g2)
  out_fn <- file.path(OUT_DIR, sprintf("%s_b%03d.rds", cell_tag, rep_id))
  if (file.exists(out_fn)) return(invisible(NULL))

  variants <- list()

  # Variant A: plug-in REML + eBIC
  t0 <- Sys.time()
  vA <- tryCatch(plug_in_reml_eBIC(y, X, K_max = ANCHOR$K_max,
                                   tau2 = ANCHOR$tau2),
                 error = function(e) list(method = "plug_in_REML_eBIC",
                                          indices = integer(0), K_hat = 0L))
  vA$elapsed <- as.numeric(Sys.time() - t0, units = "secs")
  variants$A_plug_REML_eBIC <- vA

  # Variant B: marginalized delta + eBIC (= MS_L_eBIC)
  t0 <- Sys.time()
  vB <- MS_L_LMM_stepwise_fast(y, X, tau2 = ANCHOR$tau2,
                                K_max = ANCHOR$K_max,
                                criterion = "eBIC",
                                n_nodes = ANCHOR$N_delta)
  variants$B_marg_eBIC <- list(method = "marg_delta_eBIC",
                               indices = vB$indices,
                               K_hat = length(vB$indices),
                               elapsed = as.numeric(Sys.time() - t0, units="secs"))

  # Variant C: joint Schur + eBIC (= JS_L_eBIC)
  t0 <- Sys.time()
  vC <- JS_L_LMM_stepwise_fast(y, X, tau2 = ANCHOR$tau2,
                                K_max = ANCHOR$K_max,
                                criterion = "eBIC",
                                n_nodes = ANCHOR$N_delta)
  variants$C_joint_eBIC <- list(method = "joint_Schur_eBIC",
                                indices = vC$indices,
                                K_hat = length(vC$indices),
                                elapsed = as.numeric(Sys.time() - t0, units="secs"))

  # Variant D: joint Schur + JP99 (= JS_L_JP99)
  t0 <- Sys.time()
  vD <- JS_L_LMM_stepwise_fast(y, X, tau2 = ANCHOR$tau2,
                                K_max = ANCHOR$K_max,
                                criterion = "JointPosterior",
                                theta = ANCHOR$theta,
                                n_nodes = ANCHOR$N_delta)
  variants$D_bayes_JP99 <- list(method = "bayes_JP99",
                                indices = vD$indices,
                                K_hat = length(vD$indices),
                                elapsed = as.numeric(Sys.time() - t0, units="secs"))

  # Compute metrics
  metrics <- do.call(rbind, lapply(variants, function(v) {
    tp <- length(intersect(v$indices, truth))
    fp <- length(setdiff(v$indices, truth))
    fn <- length(setdiff(truth, v$indices))
    recall    <- if (length(truth) > 0) tp / length(truth) else NA
    precision <- if (length(v$indices) > 0) tp / length(v$indices) else NA
    f1 <- if (is.na(recall) || is.na(precision) || (recall + precision) == 0) 0
          else 2 * recall * precision / (recall + precision)
    data.frame(method = v$method, K_hat = v$K_hat, tp = tp, fp = fp,
               recall = recall, precision = precision, f1 = f1,
               elapsed = v$elapsed, rep = rep_id, sigma_g2 = sigma_g2,
               stringsAsFactors = FALSE)
  }))
  payload <- list(metrics = metrics, truth = truth, results = variants)
  saveRDS(payload, out_fn)
  invisible(payload)
}

#  Main loop 
cat(sprintf("Ablation: %d replicates x 2 sigma_g2 cells (cores=%d)\n",
            B, N_CORES))
for (sg2 in c(0.0, 0.5)) {
  cat(sprintf(" sigma_g2 = %.1f ... ", sg2))
  t0 <- Sys.time()
  mclapply(seq_len(B), function(b) {
    tryCatch(run_ablation_rep(b, sg2),
             error = function(e) message(sprintf("  rep %d (sg=%.1f): %s",
                                                  b, sg2, conditionMessage(e))))
    NULL
  }, mc.cores = N_CORES)
  cat(sprintf("done in %.1f min\n", as.numeric(Sys.time() - t0, units = "mins")))
}

#  Aggregate and write LaTeX table 
files <- list.files(OUT_DIR, "[.]rds$", full.names = TRUE)
cat(sprintf("Aggregating %d RDS files...\n", length(files)))
M <- do.call(rbind, lapply(files, function(f) readRDS(f)$metrics))

agg <- aggregate(cbind(f1, precision, recall, K_hat, elapsed) ~ method + sigma_g2,
                 data = M, FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                               se = sd(x, na.rm = TRUE) / sqrt(length(x))))
# flatten the matrix columns
flat <- data.frame(method = agg$method, sigma_g2 = agg$sigma_g2,
                   f1_mean = agg$f1[,"mean"], f1_se = agg$f1[,"se"],
                   prec_mean = agg$precision[,"mean"], prec_se = agg$precision[,"se"],
                   rec_mean = agg$recall[,"mean"], rec_se = agg$recall[,"se"],
                   K_mean = agg$K_hat[,"mean"], K_se = agg$K_hat[,"se"],
                   t_mean = agg$elapsed[,"mean"])
write.csv(flat, file.path(OUT_DIR, "ablation_summary.csv"), row.names = FALSE)

display <- c(plug_in_REML_eBIC = "Plug-in REML + eBIC",
             marg_delta_eBIC   = "Marginalised $\\delta$ + eBIC",
             joint_Schur_eBIC  = "Joint Schur + eBIC",
             bayes_JP99        = "Joint Schur + JP99")
order_methods <- c("plug_in_REML_eBIC", "marg_delta_eBIC",
                   "joint_Schur_eBIC", "bayes_JP99")

fmt <- function(m, se, dig = 3) sprintf("$%.*f \\pm %.*f$", dig, m, dig, se)

tex_path <- "theory/overleaf_compact/tables/tab_ablation.tex"
sink(tex_path)
cat("% Auto-generated by sim/bench_full/16_ablation.R\n")
cat("\\begin{tabular}{lcccc}\n\\toprule\n")
cat("Variant & $\\overline{F_1}$ & $\\overline{\\mathrm{Prec}}$ & ",
    "$\\overline{\\mathrm{Rec}}$ & $\\overline{\\hat K}$ \\\\\n")
for (sg2 in c(0.0, 0.5)) {
  cat(sprintf("\\midrule\n\\multicolumn{5}{l}{\\emph{$\\sigma_g^2 = %.1f$}}\\\\\n",
              sg2))
  for (mth in order_methods) {
    r <- flat[flat$method == mth & flat$sigma_g2 == sg2, , drop = FALSE]
    if (nrow(r) == 0L) next
    cat(sprintf("%s & %s & %s & %s & %s \\\\\n",
                display[mth],
                fmt(r$f1_mean, r$f1_se),
                fmt(r$prec_mean, r$prec_se),
                fmt(r$rec_mean, r$rec_se),
                fmt(r$K_mean, r$K_se, dig = 1)))
  }
}
cat("\\bottomrule\n\\end{tabular}\n")
sink()
cat(sprintf("Wrote %s\n", tex_path))
cat("Done.\n")
