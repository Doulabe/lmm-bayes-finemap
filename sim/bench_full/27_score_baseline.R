# ==============================================================================
# 27_score_baseline.R
# The most direct internal baseline: the SAME algorithm without the Bayes
# factor. Conditional Score-LMM ranks candidates by the conditional score
# statistic Q^2_{jk} computed with the identical adaptive kernel K_k,
# identical conditioning projection P_k^* (intercept + selected variants),
# identical plug-in REML treatment of delta, and eBIC stopping in the same
# form. Compared against CBF-LMM under the plug-in REML evaluation (so the
# ONLY difference is Bayes factor vs score ranking), on the same seeds.
#
# Also records, for the score path, the exact combinatorial eBIC stop
#   pen_exact(k) = k log n + 2 log C(m, k)
# alongside the sparse approximation k(log n + 2 log m), quantifying the
# effect of the log-binomial approximation on selection.
#
# Design: anchor cell (n=1000, m=5000, rho=0.95, medium beta) x
#         sigma_g2 in {0, 0.5}, B=100 replicates, checkpointed.
# Output: results/bench_full/27_score_baseline/{cell}_b*.rds + summary CSV
# Usage:  Rscript sim/bench_full/27_score_baseline.R [--B 100] [--cores 6]
# ==============================================================================

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressPackageStartupMessages(library(parallel))
source("sim/bench_full/00_config.R")
.try_source(c("R/LMM_reml_bf.R", "../R/LMM_reml_bf.R", "../../R/LMM_reml_bf.R"))

OUT_DIR <- "results/bench_full/27_score_baseline"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
args <- parse_args()
B  <- if (!is.null(args$B)) as.integer(args$B) else 100L
NC <- if (!is.null(args$cores)) as.integer(args$cores) else
        max(1L, detectCores() - 2L)
K_MAX <- 10L

## Conditional Score-LMM: adaptive kernel + P_k^* conditioning + Q^2 ranking.
## Mirrors the CBF-LMM path with the Bayes factor replaced by the score.
score_stepwise <- function(y, X, K_max = 10L) {
  n <- nrow(X); m <- ncol(X)
  selected <- integer(0); remaining <- seq_len(m); F_k <- matrix(1, n, 1L)
  # eBIC paths under both penalties (OLS profile likelihood, gamma = 1)
  ll_ols <- function(Fm) {
    b <- solve(crossprod(Fm), crossprod(Fm, y))
    r <- as.numeric(crossprod(y - Fm %*% b))
    -0.5 * n * log(r / n)
  }
  pen_apx <- function(k) k * (log(n) + 2 * log(m))
  pen_exa <- function(k) k * log(n) + 2 * lchoose(m, k)
  ebic_apx <- c(-2 * ll_ols(F_k) + pen_apx(0), rep(NA_real_, K_max))
  ebic_exa <- ebic_apx
  sel_apx <- NA_integer_; sel_exa <- NA_integer_
  for (k_step in seq_len(K_max)) {
    rem <- X[, remaining, drop = FALSE]
    Kk  <- tcrossprod(rem) / max(ncol(rem) - 1L, 1L)
    ee  <- eigen(Kk, symmetric = TRUE); s <- pmax(ee$values, 0); U <- ee$vectors
    Uty0 <- crossprod(U, y); UtF0 <- crossprod(U, F_k)
    # plug-in REML delta by profile likelihood (same treatment as CBF-REML)
    negll <- function(ld) {
      d <- exp(ld); w <- 1 / (1 + d * s)
      A <- crossprod(UtF0 * w, UtF0)
      ch <- tryCatch(chol(A), error = function(e) NULL)
      if (is.null(ch)) return(1e10)
      b <- backsolve(ch, forwardsolve(t(ch), crossprod(UtF0 * w, Uty0)))
      r <- as.numeric(crossprod((Uty0 - UtF0 %*% b) * sqrt(w)))
      0.5 * (sum(log1p(d * s)) + (n - ncol(F_k)) * log(r) +
             2 * sum(log(diag(ch))))
    }
    o <- tryCatch(optimize(negll, c(log(1e-4), log(1e3))),
                  error = function(e) NULL)
    dh <- if (is.null(o)) 1 else exp(o$minimum)
    w  <- 1 / sqrt(1 + dh * s)
    Uty <- Uty0 * w; UtF <- UtF0 * w; UtX <- crossprod(U, rem) * w
    A <- crossprod(UtF); Ai <- solve(A)
    ty <- Uty - UtF %*% (Ai %*% crossprod(UtF, Uty))
    tX <- UtX - UtF %*% (Ai %*% crossprod(UtF, UtX))
    D  <- colSums(tX^2); u <- as.numeric(crossprod(tX, ty))
    df <- n - ncol(F_k)
    Q2 <- u^2 / (D * as.numeric(crossprod(ty)) / df)
    jb <- which.max(Q2)
    F_new <- cbind(F_k, X[, remaining[jb], drop = FALSE])
    m2ll <- -2 * ll_ols(F_new)
    ebic_apx[k_step + 1L] <- m2ll + pen_apx(k_step)
    ebic_exa[k_step + 1L] <- m2ll + pen_exa(k_step)
    if (is.na(sel_apx) && ebic_apx[k_step + 1L] > ebic_apx[k_step])
      sel_apx <- k_step - 1L
    if (is.na(sel_exa) && ebic_exa[k_step + 1L] > ebic_exa[k_step])
      sel_exa <- k_step - 1L
    if (!is.na(sel_apx) && !is.na(sel_exa)) break
    selected <- c(selected, remaining[jb]); F_k <- F_new
    remaining <- remaining[-jb]
  }
  if (is.na(sel_apx)) sel_apx <- length(selected)
  if (is.na(sel_exa)) sel_exa <- length(selected)
  list(idx_apx = selected[seq_len(min(sel_apx, length(selected)))],
       idx_exa = selected[seq_len(min(sel_exa, length(selected)))])
}

met <- function(sel, truth, K_star = 5L) {
  tp <- length(intersect(sel, truth)); K <- length(sel)
  rec <- tp / K_star; prec <- if (K > 0) tp / K else 0
  c(K_hat = K, recall = rec, precision = prec,
    f1 = if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0)
}

one_rep <- function(sg, b) {
  ck <- file.path(OUT_DIR, sprintf("sg%.1f_b%03d.rds", sg, b))
  if (file.exists(ck) && file.size(ck) > 0) return(invisible("cached"))
  d <- gen_dataset(n = 1000L, m = 5000L, rho = 0.95, K_true = 5L,
                   beta_true = c(0.8, 0.4, 0.4, 0.2, 0.2), sigma_g2 = sg,
                   block_size = 10L, seed = 20270100L + 1000L * b +
                                            as.integer(round(100 * sg)))
  sc <- score_stepwise(d$y, d$X, K_max = K_MAX)
  cb <- MS_L_LMM_stepwise_fast(d$y, d$X, tau2 = 0.04, K_max = K_MAX,
                               criterion = "eBIC", delta_eval = "reml")
  out <- rbind(
    data.frame(method = "Score_eBIC",       t(met(sc$idx_apx, d$truth))),
    data.frame(method = "Score_eBIC_exact", t(met(sc$idx_exa, d$truth))),
    data.frame(method = "CBF_REML_eBIC",    t(met(cb$indices, d$truth))))
  out$sigma_g2 <- sg; out$rep <- b
  saveRDS(out, ck)
  invisible("done")
}

cells <- expand.grid(sg = c(0, 0.5), b = seq_len(B))
cat(sprintf("Score baseline: %d reps x 2 regimes, mc.cores=%d\n", B, NC))
t0 <- Sys.time()
res <- mclapply(seq_len(nrow(cells)), function(i)
  tryCatch(one_rep(cells$sg[i], cells$b[i]),
           error = function(e) { message("ERR sg", cells$sg[i], " b",
                                          cells$b[i], ": ",
                                          conditionMessage(e)); "error" }),
  mc.cores = NC, mc.preschedule = FALSE)
cat(sprintf("done in %.1f min: %s\n",
            as.numeric(Sys.time() - t0, units = "mins"),
            paste(names(table(unlist(res))), table(unlist(res)),
                  sep = "=", collapse = ", ")))
files <- list.files(OUT_DIR, pattern = "_b[0-9]+\\.rds$", full.names = TRUE)
all <- do.call(rbind, lapply(files, readRDS))
agg <- aggregate(cbind(f1, recall, precision, K_hat) ~ method + sigma_g2,
                 all, mean)
se  <- aggregate(f1 ~ method + sigma_g2, all,
                 function(v) sd(v) / sqrt(length(v)))
agg$f1_se <- se$f1[match(paste(agg$method, agg$sigma_g2),
                          paste(se$method, se$sigma_g2))]
write.csv(all, file.path(OUT_DIR, "score_baseline_raw.csv"), row.names = FALSE)
write.csv(agg, file.path(OUT_DIR, "score_baseline_summary.csv"),
          row.names = FALSE)
cat("\n=== Score vs CBF (plug-in REML, same seeds, anchor) ===\n")
print(format(agg[order(agg$sigma_g2, -agg$f1), ], digits = 3),
      row.names = FALSE)
