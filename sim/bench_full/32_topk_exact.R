# ==============================================================================
# 32_topk_exact.R
# Recompute the CBF-LMM top-K* ranking rows of the threshold-independent
# benchmark with the exact engine: the stepwise path is forced to K* = 5
# steps (no stopping), and the K* selected indices are the top-K* ranking.
# Datasets are rebuilt with the original seed formulas of
# 10_threshold_independent.R (offsets 0/1/2) and 19_threshold_full_grid.R
# (offsets 10/20/30/40/50). Payloads are cloned with the MS_L_eBIC entry
# replaced (old kept as MS_L_eBIC_sharedGk).
# Output: results/bench_full_exact/topk/<same filenames>
# Usage:  Rscript sim/bench_full/32_topk_exact.R [--cores 5]
# ==============================================================================

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressPackageStartupMessages(library(parallel))
source("sim/bench_full/00_config.R")
source("R/CBF_LMM_exact.R")

IN_DIRS <- c("results/bench_full/10_threshold_indep",
             "results/bench_full/19_threshold_full_grid")
OUT_DIR <- "results/bench_full_exact/topk"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(flag, default) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) as.integer(args[i + 1L]) else default
}
N_CORES <- arg_val("--cores", 5L)

CELLS <- list(
  anchor    = list(n = 1000L, m = 5000L,  rho = 0.95,
                   beta = c(0.8, 0.4, 0.4, 0.2, 0.2), off = 0L),
  n3000     = list(n = 3000L, m = 5000L,  rho = 0.95,
                   beta = c(0.8, 0.4, 0.4, 0.2, 0.2), off = 1L),
  rho098    = list(n = 1000L, m = 5000L,  rho = 0.98,
                   beta = c(0.8, 0.4, 0.4, 0.2, 0.2), off = 2L),
  n500      = list(n = 500L,  m = 5000L,  rho = 0.95,
                   beta = c(0.8, 0.4, 0.4, 0.2, 0.2), off = 10L),
  rho080    = list(n = 1000L, m = 5000L,  rho = 0.80,
                   beta = c(0.8, 0.4, 0.4, 0.2, 0.2), off = 20L),
  sigweak   = list(n = 1000L, m = 5000L,  rho = 0.95,
                   beta = c(0.4, 0.2, 0.2, 0.1, 0.1), off = 30L),
  sigstrong = list(n = 1000L, m = 5000L,  rho = 0.95,
                   beta = c(1.6, 0.8, 0.8, 0.4, 0.4), off = 40L),
  m10000    = list(n = 1000L, m = 10000L, rho = 0.95,
                   beta = c(0.8, 0.4, 0.4, 0.2, 0.2), off = 50L))

work <- do.call(rbind, lapply(IN_DIRS, function(dd) {
  fs <- list.files(dd, pattern = "^[a-z0-9]+_sg[0-9.]+_b[0-9]+\\.rds$")
  if (!length(fs)) return(NULL)
  data.frame(dir = dd, file = fs, stringsAsFactors = FALSE)
}))
# put the n3000 cells last (heaviest)
work <- work[order(grepl("^n3000", work$file)), ]

one_file <- function(dd, f) {
  out_fn <- file.path(OUT_DIR, f)
  if (file.exists(out_fn) && file.size(out_fn) > 0) return("cached")
  tag <- sub("_sg.*$", "", f)
  cc  <- CELLS[[tag]]
  if (is.null(cc)) return(paste("skip unknown cell", tag))
  sg  <- as.numeric(sub("^[a-z0-9]+_sg([0-9.]+)_b.*$", "\\1", f))
  b   <- as.integer(sub(".*_b([0-9]+)\\.rds$", "\\1", f))
  seed <- 20260425L + 1000L * (b - 1L) + cc$off

  d <- gen_dataset(n = cc$n, m = cc$m, rho = cc$rho, K_true = 5L,
                   beta_true = cc$beta, sigma_g2 = sg, block_size = 10L,
                   seed = seed)

  orig <- NULL
  for (att in 1:6) {
    orig <- tryCatch(readRDS(file.path(dd, f)), error = function(e) NULL)
    if (!is.null(orig)) break
    Sys.sleep(2 * att)
  }
  if (is.null(orig)) stop("unreadable original: ", f)
  stopifnot(identical(sort(orig$truth), sort(d$truth)))

  fit <- CBF_LMM_stepwise_exact(d$y, d$X, tau2 = ANCHOR$tau2,
                                  n_nodes = ANCHOR$N_delta, force_K = 5L)
  rec <- length(intersect(fit$indices, d$truth)) / 5

  orig$metrics$MS_L_eBIC_sharedGk <- orig$metrics$MS_L_eBIC
  orig$metrics$MS_L_eBIC$top_K_recall <- rec
  saveRDS(orig, out_fn)
  "done"
}

message(sprintf("topK exact: %d replicates on %d cores", nrow(work), N_CORES))
invisible(mclapply(seq_len(nrow(work)), function(i)
  tryCatch(one_file(work$dir[i], work$file[i]), error = function(e)
    message("FAIL ", work$file[i], ": ", conditionMessage(e))),
  mc.cores = N_CORES, mc.preschedule = FALSE))

fs <- list.files(OUT_DIR, pattern = "rds$", full.names = TRUE)
res <- do.call(rbind, lapply(fs, function(f) {
  p <- readRDS(f)
  data.frame(cell = p$cell_tag, sg = p$sigma_g2,
             exact = p$metrics$MS_L_eBIC$top_K_recall,
             shared = p$metrics$MS_L_eBIC_sharedGk$top_K_recall,
             susie = p$metrics$SuSiE$top_K_recall,
             bslmm = p$metrics$BSLMM$top_K_recall)
}))
cat("\n════ Top-K* recall poolé ════\n")
print(aggregate(cbind(exact, shared, susie, bslmm) ~ sg, res,
                function(x) mean(x, na.rm = TRUE)), digits = 3)
write.csv(res, file.path(OUT_DIR, "topk_exact_raw.csv"), row.names = FALSE)
message("done.")
