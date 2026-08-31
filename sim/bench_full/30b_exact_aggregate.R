# ==============================================================================
# 30b_exact_aggregate.R
# Build raw CSVs from the exact-rerun payloads (results/bench_full_exact/) and
# print the pooled aggregate (n >= 1000, REML_stepwise excluded — same pooling
# rule as make_aggregate_table.R), comparing the exact CBF-LMM row with the
# retired shared-kernel row and the external comparators.
#
# Usage: Rscript sim/bench_full/30b_exact_aggregate.R
# ==============================================================================

ROOT <- "results/bench_full_exact"
AXES <- c("01_scaling_n", "02_scaling_rho", "03_arch", "04_scaling_m")

all_rows <- list()
for (ax in AXES) {
  fs <- list.files(file.path(ROOT, ax), pattern = "_b[0-9]+\\.rds$",
                   full.names = TRUE)
  rows <- lapply(fs, function(f) readRDS(f)$metrics)
  df <- do.call(rbind, lapply(rows, function(r) {
    r[, intersect(names(r), c("method","K_hat","tp","fp","recall","precision",
                              "f1","elapsed","n","sigma_g2","rep","rho","m",
                              "signal")), drop = FALSE]
  }))
  write.csv(df, file.path(ROOT, paste0(ax, "_raw.csv")), row.names = FALSE)
  df$axis <- ax
  all_rows[[ax]] <- df
  cat(sprintf("%s: %d fichiers, %d lignes\n", ax, length(fs), nrow(df)))
}

common <- Reduce(intersect, lapply(all_rows, names))
pool <- do.call(rbind, lapply(all_rows, function(d) d[, common, drop = FALSE]))
if ("n" %in% names(pool)) pool <- pool[pool$n >= 1000, , drop = FALSE]
pool <- pool[pool$method != "REML_stepwise", , drop = FALSE]

agg <- aggregate(cbind(K_hat, recall, precision, f1) ~ method + sigma_g2,
                 data = pool, FUN = mean)
agg <- agg[order(agg$sigma_g2, -agg$f1), ]
cat("\n════ AGRÉGAT POOLÉ (n ≥ 1000) ════\n")
print(agg, digits = 3, row.names = FALSE)
write.csv(agg, file.path(ROOT, "aggregate_pooled.csv"), row.names = FALSE)
