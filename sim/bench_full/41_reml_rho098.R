# BF-REML plug-in arm on the rho=0.98 cells (for the same-delta ranking
# ablation), same seeds as axis 02. Checkpointed.
Sys.setenv(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1",
           VECLIB_MAXIMUM_THREADS="1", MKL_NUM_THREADS="1")
suppressPackageStartupMessages(library(parallel))
source("sim/bench_full/00_config.R"); source("R/CBF_LMM_exact.R")
OUT <- "results/bench_full_exact/41_reml_rho098"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
one <- function(sg, b) {
  ck <- file.path(OUT, sprintf("rho098_sg%.1f_b%03d.rds", sg, b))
  if (file.exists(ck) && file.size(ck) > 0) return("cached")
  seed <- 20260425L + 1000L*(b-1L) + as.integer(round(1000*0.98)) +
          as.integer(round(100*sg))
  d <- gen_dataset(n=1000L, m=ANCHOR$m, rho=0.98, K_true=ANCHOR$K_true,
                   beta_true=ANCHOR$beta_true, sigma_g2=sg,
                   block_size=ANCHOR$block_size, seed=seed)
  fit <- CBF_LMM_stepwise_exact(d$y, d$X, tau2=ANCHOR$tau2, K_max=ANCHOR$K_max,
                                delta_eval="reml", n_nodes=ANCHOR$N_delta)
  tp <- length(intersect(fit$indices, d$truth))
  saveRDS(list(sg=sg, b=b, K_hat=fit$K_hat, tp=tp,
               recall=tp/length(d$truth),
               precision=if (fit$K_hat>0) tp/fit$K_hat else 0), ck)
  "done"
}
g <- expand.grid(sg=c(0,0.5), b=1:100)
invisible(mclapply(seq_len(nrow(g)), function(i)
  tryCatch(one(g$sg[i], g$b[i]), error=function(e) message("FAIL ", i)),
  mc.cores=5L, mc.preschedule=FALSE))
message("done.")
