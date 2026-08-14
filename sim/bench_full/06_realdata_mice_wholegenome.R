# ==============================================================================
# 06_realdata_mice_wholegenome.R
# Whole-genome joint analysis on the BGLR mouse panel (Obesity.BMI, n=1814,
# m_total = 10,346 SNPs after MAF >0.05 filter). All 8 methods run on the
# *entire* SNP matrix — no chromosome-wise tiling.
#
# This complements the chromosome-wise analysis (`05_realdata_mice_full.R`):
# rather than constraining the algorithm to "1 SNP per chromosome", we let
# each method decide K_total adaptively across the whole genome.
#
# Output: results/bench_full/06_realdata_mice_wholegenome/
#   - whole_genome.rds
#   - whole_genome_summary.csv
#
# CLI (called by run_all.R or standalone):
#   Rscript 06_realdata_mice_wholegenome.R [--K_max 30] [--include-bslmm]
# ==============================================================================

# Working-directory-robust source
.try_source_first <- function(candidates) {
  for (p in candidates) if (file.exists(p)) { source(p); return(p) }
  stop("Could not locate any of: ", paste(candidates, collapse = ", "),
       "\n  CWD = ", getwd())
}
BENCH_DIR <- .try_source_first(c(
  "sim/bench_full/00_config.R",
  "bench_full/00_config.R",
  "00_config.R"
))
BENCH_DIR    <- normalizePath(dirname(BENCH_DIR))
PROJECT_ROOT <- normalizePath(file.path(BENCH_DIR, "..", ".."))


run_wholegenome_mice <- function(K_max     = 30L,
                                    theta     = 0.99,
                                    N_delta   = 10L,
                                    tau2      = 0.04,
                                    include_bslmm = TRUE,
                                    out_dir   = file.path(PROJECT_ROOT,
                                       "results/bench_full/06_realdata_mice_wholegenome")) {

  if (!requireNamespace("BGLR", quietly = TRUE))
    stop("BGLR required.")

  message("=== Whole-genome mouse BMI analysis ===")
  message("Loading + QC ...")

  # -- Load and minimal QC -------------------------------------------------
  .env <- new.env()
  data(mice, package = "BGLR", envir = .env)
  X_full <- .env$mice.X; pheno <- .env$mice.pheno; map <- .env$mice.map

  y_full <- pheno$Obesity.BMI
  ok     <- which(is.finite(y_full))
  X_ok   <- X_full[ok, ]; y_ok <- y_full[ok]
  af     <- colMeans(X_ok, na.rm = TRUE) / 2
  maf    <- pmin(af, 1 - af)
  keep   <- which(is.finite(maf) & maf > 0.05)
  X      <- scale(X_ok[, keep])
  map_kept <- map[keep, , drop = FALSE]

  n <- nrow(X); m <- ncol(X)
  message(sprintf("  n=%d, m=%d (whole-genome, %d chromosomes pooled)",
                  n, m, length(unique(map_kept$chr))))
  message(sprintf("  K_max=%d, theta=%.2f, N_delta=%d", K_max, theta, N_delta))

  # -- Run all methods on full (X, y) -------------------------------------
  d <- list(X = X, y = y_ok, truth = NULL)
  methods <- default_methods(include_bslmm = include_bslmm)
  message("Methods: ", paste(methods, collapse = ", "))
  message("Running on full panel (this will take ~30-60 min) ...")

  res <- run_methods(d, methods, K_max = K_max, theta = theta,
                       N_delta = N_delta, tau2 = tau2)

  # Decorate selections with SNP names, chromosomes, and per-chromosome counts
  CHR_LEVELS <- c(as.character(1:19), "X")
  for (m_name in names(res)) {
    r <- res[[m_name]]
    sel <- r$indices
    if (length(sel) > 0L) {
      res[[m_name]]$snp_names  <- map_kept$snp_id[sel]
      res[[m_name]]$snp_chrs   <- as.character(map_kept$chr[sel])
      # per-chromosome counts vector aligned to CHR_LEVELS (zeros included)
      cnt <- table(factor(res[[m_name]]$snp_chrs, levels = CHR_LEVELS))
      res[[m_name]]$chr_counts <- as.integer(cnt)
      names(res[[m_name]]$chr_counts) <- CHR_LEVELS
      # explicit target chromosomes (those with K>=1)
      res[[m_name]]$target_chrs <- CHR_LEVELS[cnt > 0]
    } else {
      res[[m_name]]$snp_names  <- character(0)
      res[[m_name]]$snp_chrs   <- character(0)
      res[[m_name]]$chr_counts <- setNames(rep(0L, length(CHR_LEVELS)),
                                              CHR_LEVELS)
      res[[m_name]]$target_chrs <- character(0)
    }
  }

  # Cross-method consensus: chromosomes hit by >=k methods
  chr_matrix <- do.call(cbind, lapply(res, function(r) r$chr_counts))
  rownames(chr_matrix) <- CHR_LEVELS
  consensus_count <- rowSums(chr_matrix > 0)  # # methods hitting each chrom
  consensus_chr_2 <- CHR_LEVELS[consensus_count >= 2]
  consensus_chr_3 <- CHR_LEVELS[consensus_count >= 3]
  consensus_chr_4 <- CHR_LEVELS[consensus_count >= 4]

  # -- Save ----------------------------------------------------------------
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  payload <- list(
    n = n, m = m,
    chrom_distribution = table(map_kept$chr),
    map = map_kept,
    chr_levels = CHR_LEVELS,
    chr_matrix = chr_matrix,                 # (chr × method) selection counts
    consensus = list(
      hit_count_per_chr = setNames(consensus_count, CHR_LEVELS),
      chr_hit_by_at_least_2 = consensus_chr_2,
      chr_hit_by_at_least_3 = consensus_chr_3,
      chr_hit_by_at_least_4 = consensus_chr_4),
    results = res)
  saveRDS(payload, file.path(out_dir, "whole_genome.rds"))

  # Summary CSV (one row per method)
  summary_rows <- list()
  for (m_name in names(res)) {
    r <- res[[m_name]]
    cc <- r$chr_counts
    chr_counts_str <- paste(sprintf("chr%s:%d", names(cc[cc > 0]), cc[cc > 0]),
                              collapse = ",")
    target_chr_str <- paste(r$target_chrs, collapse = ",")
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      method   = m_name,
      K_hat    = r$K_hat,
      n_chrom_with_signal = length(r$target_chrs),
      target_chromosomes  = target_chr_str,
      top_snp  = if (length(r$snp_names) > 0L) r$snp_names[1] else NA_character_,
      top_chr  = if (length(r$snp_chrs)  > 0L) r$snp_chrs[1]  else NA_character_,
      elapsed_s = r$elapsed,
      chr_distribution = chr_counts_str,
      stringsAsFactors = FALSE)
  }
  summary_df <- do.call(rbind, summary_rows)
  write.csv(summary_df, file.path(out_dir, "whole_genome_summary.csv"),
            row.names = FALSE)

  # Long-format CSV (one row per (method, chr) with selection count) for plots
  long_rows <- list()
  for (m_name in names(res)) {
    cc <- res[[m_name]]$chr_counts
    long_rows[[length(long_rows) + 1L]] <- data.frame(
      method = m_name, chr = names(cc), K_hat = as.integer(cc),
      stringsAsFactors = FALSE)
  }
  long_df <- do.call(rbind, long_rows)
  write.csv(long_df,
            file.path(out_dir, "whole_genome_per_chromosome.csv"),
            row.names = FALSE)

  # Print compact summary highlighting target chromosomes
  message("\n=== Whole-genome results: target chromosomes ===")
  for (m_name in names(res)) {
    r <- res[[m_name]]
    if (r$K_hat > 0L) {
      tgt_str <- paste(r$target_chrs, collapse = ",")
      message(sprintf("  %-15s K=%4d on %2d chr [%s]  (%.1fs)",
                      m_name, r$K_hat, length(r$target_chrs),
                      tgt_str, r$elapsed))
    } else {
      message(sprintf("  %-15s K=   0  (none)  (%.1fs)",
                      m_name, r$elapsed))
    }
  }
  message("\n=== Cross-method consensus ===")
  message(sprintf("  chr hit by >=2 methods (%d): %s",
                  length(consensus_chr_2), paste(consensus_chr_2, collapse=",")))
  message(sprintf("  chr hit by >=3 methods (%d): %s",
                  length(consensus_chr_3), paste(consensus_chr_3, collapse=",")))
  message(sprintf("  chr hit by >=4 methods (%d): %s",
                  length(consensus_chr_4), paste(consensus_chr_4, collapse=",")))

  invisible(payload)
}


# --- CLI entry point -----------------------------------------------------
if (sys.nframe() == 0L) {
  args <- parse_args()
  K_max         <- if (!is.null(args$K_max)) as.integer(args$K_max) else 30L
  include_bslmm <- !isTRUE(args[["no-bslmm"]])  # bslmm by default

  saved_wd <- getwd()
  setwd(PROJECT_ROOT)
  on.exit(setwd(saved_wd), add = TRUE)

  run_wholegenome_mice(K_max = K_max, include_bslmm = include_bslmm)
}
