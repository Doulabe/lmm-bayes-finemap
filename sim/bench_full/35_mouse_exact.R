# ==============================================================================
# 35_mouse_exact.R
# Primary autosomal mouse BMI analysis with the exact candidate-specific
# engine (K_jk + profile-ML eBIC), replacing the shared-kernel CBF-LMM row of
# 12_mouse_autosomes_only.R. Comparator rows are unchanged (their thresholds
# were not crossed in the original analysis).
# Output: results/bench_full_exact/12_mouse_autosomes/mouse_exact.rds
# Usage:  Rscript sim/bench_full/35_mouse_exact.R
# ==============================================================================

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
source("sim/bench_full/00_config.R")
source("R/CBF_LMM_exact.R")
suppressPackageStartupMessages(library(BGLR))

OUT_DIR <- "results/bench_full_exact/12_mouse_autosomes"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
ck <- file.path(OUT_DIR, "mouse_exact.rds")
if (file.exists(ck)) { cat("[cached]\n"); print(readRDS(ck)$summary); quit() }

data(mice)
autosome_idx <- which(mice.map$chr != "X")
X_aut <- scale(mice.X[, autosome_idx]); X_aut[!is.finite(X_aut)] <- 0
y <- as.numeric(scale(mice.pheno$Obesity.BMI))
map_aut <- mice.map[autosome_idx, ]
cat(sprintf("mouse autosomes: n=%d, m=%d\n", nrow(X_aut), ncol(X_aut)))

t0  <- Sys.time()
fit <- CBF_LMM_stepwise_exact(y, X_aut, tau2 = ANCHOR$tau2,
                                K_max = ANCHOR$K_max,
                                n_nodes = ANCHOR$N_delta)
el  <- as.numeric(Sys.time() - t0, units = "mins")

sel_info <- if (fit$K_hat > 0) {
  data.frame(marker = colnames(mice.X)[autosome_idx][fit$indices],
             chr = map_aut$chr[fit$indices],
             logBF = fit$log_BF_at_step)
} else data.frame()

out <- list(summary = list(K_hat = fit$K_hat, indices = fit$indices,
                            selected = sel_info,
                            ebic_path = fit$ebic_path,
                            elapsed_min = el),
            variant = fit$variant)
saveRDS(out, ck)
cat(sprintf("[done] exact mouse: K_hat=%d (%.1f min)\n", fit$K_hat, el))
if (fit$K_hat > 0) print(sel_info)
