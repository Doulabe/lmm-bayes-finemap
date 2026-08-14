# ==============================================================================
# 12_mouse_autosomes_only.R
#
# Autosomes-only re-analysis of mouse BMI as a sensitivity check on the
# chromosome-X selection reported in Section 7 of the manuscript.
#
# Motivation: chromosome X in mice is subject to dosage compensation
# (random X-inactivation in females, ~15-25% of X-linked genes escaping
# inactivation, monoallelic expression patterns).  Per-SNP effect sizes on
# chromosome X are therefore confounded with these phenomena unless sex is
# jointly modelled.  Since the sex covariate is not used in the main
# whole-genome scan (which models the polygenic background through the
# GRM only), we run a clean sensitivity check by dropping the 272 X-linked
# markers and re-running all methods on the 10,074 autosomal markers.
#
# Output: results/bench_full/12_mouse_autosomes/whole_genome_autosomes.rds
#         results/bench_full/12_mouse_autosomes_summary.csv
# ==============================================================================

source("sim/bench_full/00_config.R")

OUT_DIR <- "results/bench_full/12_mouse_autosomes"
OUT_CSV <- "results/bench_full/12_mouse_autosomes_summary.csv"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

if (!requireNamespace("BGLR", quietly = TRUE))
  stop("BGLR package required for mouse data")
suppressPackageStartupMessages(library(BGLR))
data(mice)

# Extract autosomes-only subset 
autosome_idx <- which(mice.map$chr != "X")
cat(sprintf("Autosomes-only: dropping %d X-linked markers, keeping %d\n",
            sum(mice.map$chr == "X"), length(autosome_idx)))

X_aut <- mice.X[, autosome_idx]
map_aut <- mice.map[autosome_idx, ]
y <- as.numeric(mice.pheno$Obesity.BMI)
keep <- !is.na(y)
y <- y[keep]
X_aut <- X_aut[keep, ]
cat(sprintf("After dropping missing BMI: n=%d, m_aut=%d\n", length(y),
            ncol(X_aut)))

# Standardise X column-wise
X_aut <- scale(X_aut)

#  Run all methods 
K_max <- 30L
N_delta <- 10L
theta <- 0.99
tau2 <- 0.04

methods_list <- list()

run_method <- function(name, fn) {
  cat(sprintf("  Running %s...\n", name))
  t0 <- Sys.time()
  res <- tryCatch(fn(), error = function(e) NULL)
  elapsed <- as.numeric(Sys.time() - t0, units = "secs")
  if (is.null(res)) return(list(name = name, indices = integer(0),
                                   K_hat = 0L, elapsed = elapsed))
  res$elapsed <- elapsed
  res$name    <- name
  res
}

methods_list$JS_L_eBIC <- run_method("JS_L_eBIC", function()
  JS_L_LMM_stepwise_fast(y, X_aut, tau2 = tau2, K_max = K_max,
                          criterion = "eBIC",
                          theta = theta, n_nodes = N_delta))
methods_list$MS_L_eBIC <- run_method("MS_L_eBIC", function()
  MS_L_LMM_stepwise_fast(y, X_aut, tau2 = tau2, K_max = K_max,
                          criterion = "eBIC",
                          theta = theta, n_nodes = N_delta))
methods_list$JS_L_JP99 <- run_method("JS_L_JP99", function()
  JS_L_LMM_stepwise_fast(y, X_aut, tau2 = tau2, K_max = K_max,
                          criterion = "JointPosterior",
                          theta = theta, n_nodes = N_delta))
methods_list$MS_L_JP99 <- run_method("MS_L_JP99", function()
  MS_L_LMM_stepwise_fast(y, X_aut, tau2 = tau2, K_max = K_max,
                          criterion = "JointPosterior",
                          theta = theta, n_nodes = N_delta))

methods_list$SuSiE <- run_method("SuSiE", function() {
  fit <- susieR::susie(X_aut, y, L = K_max, verbose = FALSE)
  sel <- which(fit$pip > 0.99)
  list(indices = sel, K_hat = length(sel),
        scores = fit$pip[sel])
})

methods_list$BSLMM <- run_method("BSLMM", function() {
  eta <- list(list(X = X_aut, model = "BayesB"))
  fit <- suppressMessages(suppressWarnings(
    BGLR::BGLR(y = y, ETA = eta, nIter = 3000L, burnIn = 500L,
                verbose = FALSE)))
  pip <- fit$ETA[[1]]$d
  sel <- which(pip > 0.5)
  list(indices = sel, K_hat = length(sel), scores = pip[sel])
})

methods_list$fastlmm <- run_method("fastlmm", function() {
  m <- ncol(X_aut); n <- nrow(X_aut)
  K <- tcrossprod(X_aut) / max(m - 1L, 1L)
  fit <- rrBLUP::mixed.solve(y = y, K = K)
  ee <- eigen(K, symmetric = TRUE)
  delta_hat <- fit$Vu / fit$Ve
  d_inv_sqrt <- 1 / sqrt(1 + delta_hat * ee$values)
  UtX <- crossprod(ee$vectors, X_aut)
  Uty <- crossprod(ee$vectors, y)
  Z   <- as.numeric(crossprod(UtX * d_inv_sqrt, Uty * d_inv_sqrt) /
                       sqrt(colSums((UtX * d_inv_sqrt)^2) + 1e-300))
  z_crit <- qnorm(1 - 0.05 / (2 * m))
  sel <- which(abs(Z) > z_crit)
  list(indices = sel, K_hat = length(sel), scores = abs(Z)[sel])
})

#  Map back to chromosomes 
to_chr <- function(indices) {
  if (length(indices) == 0L) return(character(0))
  as.character(map_aut$chr[indices])
}

rows <- list()
for (mname in names(methods_list)) {
  r <- methods_list[[mname]]
  chrs <- to_chr(r$indices)
  chr_tab <- if (length(chrs) > 0L) table(chrs) else integer(0)
  chr_str <- if (length(chr_tab) > 0L)
    paste(sprintf("%s:%d", names(chr_tab), as.integer(chr_tab)),
          collapse = ", ") else "(none)"
  n_chr <- length(unique(chrs))
  rows[[length(rows)+1L]] <- data.frame(
    method = mname, K_hat = r$K_hat, n_chr = n_chr,
    chr_distribution = chr_str, elapsed_s = round(r$elapsed, 1),
    stringsAsFactors = FALSE)
}
df <- do.call(rbind, rows)
print(df)

write.csv(df, OUT_CSV, row.names = FALSE)
saveRDS(list(methods = methods_list, map = map_aut,
              autosome_idx = autosome_idx),
         file.path(OUT_DIR, "whole_genome_autosomes.rds"))
cat(sprintf("\nWrote %s and whole_genome_autosomes.rds\n", OUT_CSV))
