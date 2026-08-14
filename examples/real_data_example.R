# ==============================================================================
# examples/real_data_example.R
#
# Minimal real-data analysis: Valdar et al. (2006) heterogeneous-stock mouse
# panel, phenotype = Body Mass Index, with optional sex-covariate adjustment.
#
# Demonstrates the F_k = [W, X_{S_{k-1}}] generalization: passing a covariate
# matrix W via the a_n argument adjusts for sex (or any other fixed effect)
# without modifying the Bayes-factor algebra.
#
# Run from the repository root:
#   Rscript examples/real_data_example.R
# ==============================================================================

if (!requireNamespace("BGLR", quietly = TRUE))
  stop("This example requires BGLR for the mouse dataset.\n",
        "  install.packages(\"BGLR\")")

source("R/LMM_core.R")
source("R/LMM_stepwise_fast.R")

suppressPackageStartupMessages(library(BGLR))
data(mice)

# Phenotype and SNP matrix
y_full <- as.numeric(mice.pheno$Obesity.BMI)
keep   <- !is.na(y_full)
y      <- y_full[keep]
X      <- scale(mice.X[keep, ])
map    <- mice.map
n      <- length(y); m <- ncol(X)
cat(sprintf("After dropping NA BMI: n = %d, m = %d markers\n", n, m))

# -- Analysis 1: no covariates (intercept only)
cat("\n=== Analysis 1: intercept only (no sex covariate) ===\n")
t0 <- Sys.time()
res1 <- JS_L_LMM_stepwise_fast(y, X,
                                   tau2 = 0.04, K_max = 30L,
                                   criterion = "eBIC",
                                   theta = 0.99, n_nodes = 10L)
elapsed1 <- as.numeric(Sys.time() - t0, units = "secs")
chrs1 <- if (length(res1$indices) > 0L)
           paste(sort(unique(as.character(map$chr[res1$indices]))),
                   collapse = ", ") else "(none)"
cat(sprintf("  K_hat = %d, target chromosomes: %s, elapsed %.1f s\n",
              res1$K_hat, chrs1, elapsed1))

# -- Analysis 2: adjust for sex
cat("\n=== Analysis 2: sex covariate in W ===\n")
sex <- as.integer(mice.pheno$GENDER[keep] == "M")
W <- cbind(intercept = rep(1, n),
            sex       = sex - mean(sex))    # centred

t0 <- Sys.time()
res2 <- JS_L_LMM_stepwise_fast(y, X, a_n = W,
                                   tau2 = 0.04, K_max = 30L,
                                   criterion = "eBIC",
                                   theta = 0.99, n_nodes = 10L)
elapsed2 <- as.numeric(Sys.time() - t0, units = "secs")
chrs2 <- if (length(res2$indices) > 0L)
           paste(sort(unique(as.character(map$chr[res2$indices]))),
                   collapse = ", ") else "(none)"
cat(sprintf("  K_hat = %d, target chromosomes: %s, elapsed %.1f s\n",
              res2$K_hat, chrs2, elapsed2))

# -- Compare: are the same chromosomes selected ?
cat("\n=== Comparison ===\n")
cat(sprintf("  no-sex   selection on chromosomes: %s\n", chrs1))
cat(sprintf("  with-sex selection on chromosomes: %s\n", chrs2))
if (identical(chrs1, chrs2)) {
  cat("  The chromosome-level selection is robust to sex adjustment.\n")
} else {
  cat("  The chromosome-level selection changes with sex adjustment;\n")
  cat("  this indicates either a sex-by-chromosome interaction or a\n")
  cat("  sex-confounded marginal signal.\n")
}
