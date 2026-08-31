# ==============================================================================
# 37_exact_summaries_figures.R
# Build per-cell summary CSVs for the exact rerun (same format as
# 99_aggregate.R) and regenerate the two manuscript figures via
# make_cbf_figures.R logic pointed at the exact tree.
# Usage: Rscript sim/bench_full/37_exact_summaries_figures.R
# ==============================================================================

EX <- "results/bench_full_exact"
AXES <- c("01_scaling_n", "02_scaling_rho", "03_arch", "04_scaling_m")

for (ax in AXES) {
  fs <- list.files(file.path(EX, ax), pattern = "_b[0-9]+[.]rds$",
                   full.names = TRUE)
  rows <- do.call(rbind, lapply(fs, function(f) {
    m <- readRDS(f)$metrics
    if (!"signal" %in% names(m) && grepl("^sig", basename(f)))
      m$signal <- sub("^sig([a-z]+)_.*", "\\1", basename(f))
    m
  }))
  keys <- intersect(c("n", "sigma_g2", "rho", "m", "signal", "method"),
                    names(rows))
  agg <- aggregate(cbind(f1, recall, precision, K_hat, elapsed) ~ .,
                   data = rows[, c(keys, "f1", "recall", "precision",
                                   "K_hat", "elapsed")],
                   FUN = mean)
  names(agg)[names(agg) == "f1"]        <- "mean_f1"
  names(agg)[names(agg) == "recall"]    <- "mean_recall"
  names(agg)[names(agg) == "precision"] <- "mean_precision"
  names(agg)[names(agg) == "K_hat"]     <- "mean_K_hat"
  names(agg)[names(agg) == "elapsed"]   <- "mean_elapsed"
  sdv <- aggregate(f1 ~ ., data = rows[, c(keys, "f1")], FUN = sd)
  names(sdv)[names(sdv) == "f1"] <- "sd_f1"
  agg <- merge(agg, sdv, by = keys)
  agg$n_reps <- 100L
  write.csv(agg, file.path(EX, paste0(ax, "_summary.csv")), row.names = FALSE)
  cat(ax, "summary:", nrow(agg), "lignes\n")
}

# --- figures: run a path-swapped copy of make_cbf_figures.R ------------------
src <- readLines("sim/bench_full/make_cbf_figures.R")
src <- gsub("results/bench_full", EX, src, fixed = TRUE)
tmp <- file.path(EX, "make_cbf_figures_exact.R")
writeLines(src, tmp)
source(tmp, chdir = FALSE)
cat("figures régénérées\n")
