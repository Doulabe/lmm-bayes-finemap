# ==============================================================================
# 05_realdata_mice_full.R
# Chromosome-wise fine-mapping on the BGLR mouse dataset (Obesity.BMI),
# using the full marker panel (n = 1814 mice, m_total ~ 10346 SNPs across
# 20 chromosomes).
#
# For each chromosome we run all enabled methods on the within-chromosome
# X submatrix (m_chr ~ 250-880 SNPs). Results are saved per-chromosome to
# enable parallel execution and resumption.
#
# CLI:
#   Rscript 05_realdata_mice_full.R                  # run all 20 chromosomes
#   Rscript 05_realdata_mice_full.R --chr 1          # one chromosome only
#   Rscript 05_realdata_mice_full.R --chr 1 --include-bslmm
# ==============================================================================

source("sim/bench_full/00_config.R")

OUT_DIR <- "results/bench_full/05_realdata_mice"


# ------------------------------------------------------------------------------
# Load + minimal QC
# ------------------------------------------------------------------------------
load_mice <- function(maf_thresh = 0.05) {
  if (!requireNamespace("BGLR", quietly = TRUE))
    stop("BGLR package required; install with install.packages('BGLR').")
  data(mice, package = "BGLR", envir = environment())
  X_full <- get("mice.X")
  pheno  <- get("mice.pheno")
  map    <- get("mice.map")

  y_full <- pheno$Obesity.BMI
  ok <- which(is.finite(y_full))
  X_ok <- X_full[ok, , drop = FALSE]
  y_ok <- y_full[ok]

  af  <- colMeans(X_ok, na.rm = TRUE) / 2
  maf <- pmin(af, 1 - af)
  keep <- which(is.finite(maf) & maf > maf_thresh)
  X    <- scale(X_ok[, keep])
  map_kept <- map[keep, , drop = FALSE]

  cat(sprintf("  After QC + MAF>%.2f: n=%d, m=%d, %d chromosomes\n",
              maf_thresh, nrow(X), ncol(X),
              length(unique(map_kept$chr))))
  list(X = X, y = y_ok, map = map_kept)
}


# ------------------------------------------------------------------------------
# Run one chromosome
# ------------------------------------------------------------------------------
run_chromosome <- function(chr_id, dat, methods, K_max, theta, N_delta, tau2) {
  out_file <- file.path(OUT_DIR, sprintf("chr%s.rds", chr_id))
  if (file.exists(out_file) && file.size(out_file) > 0) {
    cat(sprintf("  [skip] chr%s already done\n", chr_id))
    return(invisible(NULL))
  }

  cols <- which(dat$map$chr == chr_id)
  if (length(cols) < 5L) {
    cat(sprintf("  [skip] chr%s has only %d SNPs\n", chr_id, length(cols)))
    return(invisible(NULL))
  }
  X_chr  <- dat$X[, cols, drop = FALSE]
  snp_id <- dat$map$snp_id[cols]
  cat(sprintf("\n--- Chromosome %s: m=%d SNPs ---\n", chr_id, ncol(X_chr)))

  d <- list(X = X_chr, y = dat$y, truth = NULL)
  res <- run_methods(d, methods, K_max = K_max, theta = theta,
                     N_delta = N_delta, tau2 = tau2)

  # Decorate with SNP names
  for (m_name in names(res)) {
    sel <- res[[m_name]]$indices
    res[[m_name]]$snp_names <- if (length(sel) > 0L) snp_id[sel] else character(0)
  }

  payload <- list(chr = chr_id, m = ncol(X_chr), n = nrow(X_chr),
                    snp_id = snp_id, results = res)
  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
  saveRDS(payload, out_file)

  cat(sprintf("  [done] chr%s:\n", chr_id))
  for (m_name in names(res)) {
    r <- res[[m_name]]
    sel_str <- if (length(r$snp_names) > 0L)
      paste(head(r$snp_names, 3L), collapse = ",") else "(none)"
    if (length(r$snp_names) > 3L) sel_str <- paste0(sel_str, ",...")
    cat(sprintf("    %s: K_hat=%d  selected=[%s]  time=%.1fs\n",
                m_name, r$K_hat, sel_str, r$elapsed))
  }
  invisible(out_file)
}


# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------
if (sys.nframe() == 0L) {
  args <- parse_args()
  include_bslmm <- isTRUE(args[["include-bslmm"]])
  methods <- default_methods(include_bslmm)

  K_max   <- if (!is.null(args$K_max)) as.integer(args$K_max) else 5L
  theta   <- if (!is.null(args$theta)) as.numeric(args$theta) else 0.95
  N_delta <- if (!is.null(args$N_delta)) as.integer(args$N_delta) else 10L

  cat("=== 05_realdata_mice_full: BGLR mice BMI, chromosome-wise ===\n")
  cat(sprintf("Methods: %s\n", paste(methods, collapse = ",")))
  cat("Loading + QC...\n")
  dat <- load_mice(maf_thresh = 0.05)

  chr_grid <- if (!is.null(args$chr)) as.character(args$chr)
              else as.character(c(1:19, "X"))

  for (chr_id in chr_grid) {
    run_chromosome(chr_id, dat, methods,
                     K_max = K_max, theta = theta,
                     N_delta = N_delta, tau2 = ANCHOR$tau2)
  }
  cat("\nDone.\n")
}
