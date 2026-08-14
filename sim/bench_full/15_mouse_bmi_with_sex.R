# ==============================================================================
# 15_mouse_bmi_with_sex.R
#
# Mouse BMI analysis with explicit sex covariate adjustment.
# Uses the now-matrix-aware a_n argument of the stepwise functions:
#   F_k = [W, X_{S_{k-1}}] where W = [intercept, sex]
#
# Runs two scans:
#   * whole-genome (m=10346) with sex
#   * autosomes-only (m=10074) with sex
#
# Output: results/bench_full/15_mouse_bmi_with_sex/
#         results/bench_full/15_mouse_bmi_with_sex_summary.csv
# ==============================================================================

source("sim/bench_full/00_config.R")

OUT_DIR <- "results/bench_full/15_mouse_bmi_with_sex"
OUT_CSV <- "results/bench_full/15_mouse_bmi_with_sex_summary.csv"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

if (!requireNamespace("BGLR", quietly = TRUE))
  stop("BGLR package required")
suppressPackageStartupMessages(library(BGLR))
data(mice)

#  Phenotype + sex covariate 
y_full   <- as.numeric(mice.pheno$Obesity.BMI)
sex_full <- as.integer(mice.pheno$GENDER == "M")   # 1=M, 0=F
keep <- !is.na(y_full) & !is.na(sex_full)
y    <- y_full[keep]
sex  <- sex_full[keep]
X    <- scale(mice.X[keep, ])
map  <- mice.map
cat(sprintf("After filtering: n=%d, m=%d, sex M/F=%d/%d\n",
            length(y), ncol(X), sum(sex == 1), sum(sex == 0)))

# Conditioning matrix W = [intercept, sex]
W_with_sex <- cbind(intercept = rep(1, length(y)),
                       sex       = sex - mean(sex))   # centred sex

# Method runner: takes (y, X, W, map, scope_label) returns named list of results
run_methods_with_W <- function(y, X, W, map, scope_label,
                                  K_max = 30L, N_delta = 10L, tau2 = 0.04,
                                  theta = 0.99) {
  n <- nrow(X); m <- ncol(X)
  results <- list()

  cat(sprintf("\n[%s] running framework methods with sex covariate...\n",
              scope_label))

  # JS_L_eBIC
  t0 <- Sys.time()
  res <- tryCatch(JS_L_LMM_stepwise_fast(y, X, a_n = W,
                                            tau2 = tau2, K_max = K_max,
                                            criterion = "eBIC",
                                            theta = theta, n_nodes = N_delta),
                    error = function(e) { cat("ERR JS_L_eBIC:", conditionMessage(e), "\n"); NULL })
  results$JS_L_eBIC <- list(elapsed = as.numeric(Sys.time() - t0, units = "secs"),
                              indices = if (!is.null(res)) res$indices else integer(0),
                              K_hat = if (!is.null(res)) res$K_hat else 0L)

  # MS_L_eBIC
  t0 <- Sys.time()
  res <- tryCatch(MS_L_LMM_stepwise_fast(y, X, a_n = W,
                                            tau2 = tau2, K_max = K_max,
                                            criterion = "eBIC",
                                            theta = theta, n_nodes = N_delta),
                    error = function(e) NULL)
  results$MS_L_eBIC <- list(elapsed = as.numeric(Sys.time() - t0, units = "secs"),
                              indices = if (!is.null(res)) res$indices else integer(0),
                              K_hat = if (!is.null(res)) res$K_hat else 0L)

  # JS_L_JP99
  t0 <- Sys.time()
  res <- tryCatch(JS_L_LMM_stepwise_fast(y, X, a_n = W,
                                            tau2 = tau2, K_max = K_max,
                                            criterion = "JointPosterior",
                                            theta = theta, n_nodes = N_delta),
                    error = function(e) NULL)
  results$JS_L_JP99 <- list(elapsed = as.numeric(Sys.time() - t0, units = "secs"),
                              indices = if (!is.null(res)) res$indices else integer(0),
                              K_hat = if (!is.null(res)) res$K_hat else 0L)

  # MS_L_JP99
  t0 <- Sys.time()
  res <- tryCatch(MS_L_LMM_stepwise_fast(y, X, a_n = W,
                                            tau2 = tau2, K_max = K_max,
                                            criterion = "JointPosterior",
                                            theta = theta, n_nodes = N_delta),
                    error = function(e) NULL)
  results$MS_L_JP99 <- list(elapsed = as.numeric(Sys.time() - t0, units = "secs"),
                              indices = if (!is.null(res)) res$indices else integer(0),
                              K_hat = if (!is.null(res)) res$K_hat else 0L)

  # For competitors that don't natively accept W, we residualise y against W
  # before fitting (Frisch-Waugh-Lovell argument)
  cat(sprintf("[%s] residualising y against W for competitor methods...\n",
              scope_label))
  WtW <- crossprod(W)
  Wty <- crossprod(W, y)
  coef_W <- solve(WtW, Wty)
  y_resid <- as.numeric(y - W %*% coef_W)

  # SuSiE
  t0 <- Sys.time()
  res <- tryCatch(susieR::susie(X, y_resid, L = K_max, verbose = FALSE),
                    error = function(e) NULL)
  if (!is.null(res)) sel <- which(res$pip > 0.99)
  else sel <- integer(0)
  results$SuSiE <- list(elapsed = as.numeric(Sys.time() - t0, units = "secs"),
                          indices = sel, K_hat = length(sel))

  # BSLMM
  t0 <- Sys.time()
  res <- tryCatch({
    eta <- list(list(X = X, model = "BayesB"))
    suppressMessages(suppressWarnings(
      BGLR::BGLR(y = y_resid, ETA = eta, nIter = 3000L, burnIn = 500L,
                  verbose = FALSE)))
  }, error = function(e) NULL)
  if (!is.null(res) && !is.null(res$ETA[[1]]$d)) {
    pip <- res$ETA[[1]]$d
    sel <- which(pip > 0.5)
  } else sel <- integer(0)
  results$BSLMM <- list(elapsed = as.numeric(Sys.time() - t0, units = "secs"),
                          indices = sel, K_hat = length(sel))

  # Per-method: build chromosome distribution string
  rows <- list()
  for (mname in names(results)) {
    r <- results[[mname]]
    if (length(r$indices) > 0L) {
      chrs <- as.character(map$chr[r$indices])
      chr_tab <- table(chrs)
      chr_str <- paste(sprintf("%s:%d", names(chr_tab),
                                  as.integer(chr_tab)), collapse = ", ")
      n_chr <- length(unique(chrs))
    } else { chr_str <- "(none)"; n_chr <- 0L }
    rows[[length(rows)+1L]] <- data.frame(
      scope = scope_label, method = mname,
      K_hat = r$K_hat, n_chr = n_chr,
      chr_distribution = chr_str,
      elapsed_s = round(r$elapsed, 1),
      stringsAsFactors = FALSE)
  }
  list(rows = do.call(rbind, rows), full_results = results)
}

#  Whole-genome with sex 
cat("\n=== Whole-genome with sex ===\n")
res_wg <- run_methods_with_W(y, X, W_with_sex, map, "whole_with_sex")
print(res_wg$rows)
saveRDS(res_wg, file.path(OUT_DIR, "whole_with_sex.rds"))

#  Autosomes-only with sex 
cat("\n=== Autosomes-only with sex ===\n")
autosome_idx <- which(map$chr != "X")
X_aut   <- X[, autosome_idx]
map_aut <- map[autosome_idx, ]
res_aut <- run_methods_with_W(y, X_aut, W_with_sex, map_aut, "autosomes_with_sex")
print(res_aut$rows)
saveRDS(res_aut, file.path(OUT_DIR, "autosomes_with_sex.rds"))

# Combined CSV
df <- rbind(res_wg$rows, res_aut$rows)
write.csv(df, OUT_CSV, row.names = FALSE)
cat(sprintf("\nWrote %s\n", OUT_CSV))
print(df)
