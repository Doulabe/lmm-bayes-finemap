# ==============================================================================
# 99_aggregate.R
# Collect all per-(cell, rep) RDS files and produce summary CSVs and figures.
#
# Outputs (all under results/bench_full/):
#   - 01_scaling_n_summary.csv
#   - 02_scaling_rho_summary.csv
#   - 03_arch_summary.csv
#   - 04_scaling_m_summary.csv
#   - 05_realdata_mice_summary.csv
#   - fig_scaling_n.pdf, fig_scaling_rho.pdf, fig_arch.pdf, fig_scaling_m.pdf
#   - fig_realdata_mice_chr.pdf
#
# CLI:
#   Rscript 99_aggregate.R                    # all sims + figures
#   Rscript 99_aggregate.R --skip-figures     # just CSVs
# ==============================================================================

source("sim/bench_full/00_config.R")

suppressPackageStartupMessages({
  library(ggplot2)
})

method_colors <- c(
  "JS_L_eBIC"     = "#D55E00",   # vermillion
  "JS_L_JP99"     = "#A04000",   # darker vermillion
  "MS_L_eBIC"     = "#E69F00",   # orange
  "MS_L_JP99"     = "#A66B00",   # darker orange
  # legacy names (in case old RDS files exist)
  "JS_L"          = "#D55E00",
  "MS_L"          = "#E69F00",
  "SuSiE"         = "#0072B2",
  "BayesR"        = "#009E73",
  "BSLMM"         = "#CC79A7",
  "REML_stepwise" = "#56B4E9",
  "fastlmm"       = "#7570B3"    # purple — LMM-frequentist baseline
)

# Sample-size filter for n-scaling cell (skip n=500 — undersized regime,
# per binary GWAS companion paper). Apply only when --filter-n is set.
filter_n_floor <- function(df, n_min = 1000L) {
  if (is.null(df) || nrow(df) == 0L) return(df)
  if ("n" %in% names(df)) df <- df[df$n >= n_min, , drop = FALSE]
  df
}

theme_pub <- function() {
  theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.3, colour = "grey92"),
          strip.background = element_rect(fill = "grey95", colour = "grey80"),
          strip.text       = element_text(face = "bold"),
          legend.position  = "bottom",
          axis.title       = element_text(face = "bold"))
}


# 
# Helpers
# 
collect_sim <- function(dir) {
  files <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(files) == 0L) return(NULL)
  rows <- lapply(files, function(f) {
    p <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(p) || is.null(p$metrics)) return(NULL)
    p$metrics
  })
  do.call(rbind, rows)
}


summarise_sim <- function(df, group_vars) {
  if (is.null(df) || nrow(df) == 0L) return(NULL)
  df %>%
    group_by(!!!syms(c(group_vars, "method"))) %>%
    summarise(n_reps = n(),
              mean_f1 = mean(f1, na.rm = TRUE),
              sd_f1 = sd(f1, na.rm = TRUE),
              mean_recall = mean(recall, na.rm = TRUE),
              mean_precision = mean(precision, na.rm = TRUE),
              mean_K_hat = mean(K_hat, na.rm = TRUE),
              mean_elapsed = mean(elapsed, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(!!!syms(group_vars), desc(mean_f1))
}


plot_axis <- function(df_summary, x_var, label_x, out_pdf, n_facet_rows = 1L) {
  if (is.null(df_summary) || nrow(df_summary) == 0L) return(invisible(NULL))
  df <- df_summary %>%
    mutate(method = factor(method, levels = names(method_colors)),
           sigma_g2 = factor(sigma_g2,
                              levels = sort(unique(sigma_g2)),
                              labels = sprintf("sigma_g^2 = %.1f",
                                                 sort(unique(sigma_g2)))))
  p <- ggplot(df, aes(x = factor(!!sym(x_var)), y = mean_f1,
                       colour = method, group = method)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 2.0) +
    geom_errorbar(aes(ymin = pmax(0, mean_f1 - sd_f1 / sqrt(n_reps) * 1.96),
                       ymax = pmin(1, mean_f1 + sd_f1 / sqrt(n_reps) * 1.96)),
                   width = 0.15, linewidth = 0.4) +
    facet_wrap(~ sigma_g2, nrow = n_facet_rows) +
    scale_colour_manual(values = method_colors, name = "Method") +
    coord_cartesian(ylim = c(0, 1.05)) +
    labs(x = label_x, y = expression("Mean " * F[1])) +
    theme_pub()
  ggsave(out_pdf, p, width = 9, height = 4.5, units = "in")
  message("Wrote ", out_pdf)
  invisible(p)
}


# 
# Real-data mouse aggregator
# 
collect_mice <- function(dir = "results/bench_full/05_realdata_mice") {
  files <- list.files(dir, pattern = "^chr.*\\.rds$", full.names = TRUE)
  if (length(files) == 0L) return(NULL)
  rows <- list()
  for (f in files) {
    p <- readRDS(f)
    chr <- p$chr
    for (m_name in names(p$results)) {
      r <- p$results[[m_name]]
      rows[[length(rows) + 1L]] <- data.frame(
        chr = chr, method = m_name,
        K_hat = r$K_hat, m_chr = p$m,
        elapsed = r$elapsed,
        top_snps = paste(head(r$snp_names, 5L), collapse = ","),
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}


plot_mice_chr <- function(df, out_pdf) {
  if (is.null(df) || nrow(df) == 0L) return(invisible(NULL))
  chr_levels <- c(as.character(1:19), "X")
  df <- df %>% mutate(chr = factor(chr, levels = chr_levels),
                       method = factor(method, levels = names(method_colors)))
  p <- ggplot(df, aes(x = chr, y = K_hat, fill = method)) +
    geom_col(position = position_dodge(0.85), width = 0.8,
              colour = "black", linewidth = 0.15) +
    scale_fill_manual(values = method_colors, name = "Method") +
    labs(x = "Chromosome", y = expression(hat(K))) +
    theme_pub()
  ggsave(out_pdf, p, width = 11, height = 4, units = "in")
  message("Wrote ", out_pdf)
  invisible(p)
}


# 
# Whole-genome mice aggregator
# 
collect_mice_whole <- function(
    dir = "results/bench_full/06_realdata_mice_wholegenome") {
  fn <- file.path(dir, "whole_genome.rds")
  if (!file.exists(fn)) return(NULL)
  p <- readRDS(fn)
  rows <- list()
  for (m_name in names(p$results)) {
    r <- p$results[[m_name]]
    cc <- r$chr_counts
    chr_str <- if (length(cc) > 0L && any(cc > 0))
      paste(sprintf("chr%s:%d", names(cc[cc > 0]), cc[cc > 0]),
              collapse = ",") else ""
    target_str <- paste(r$target_chrs, collapse = ",")
    rows[[length(rows) + 1L]] <- data.frame(
      method   = m_name,
      K_hat    = r$K_hat,
      n_chrom_with_signal = length(r$target_chrs),
      target_chromosomes  = target_str,
      top_snp  = if (length(r$snp_names) > 0L) r$snp_names[1] else NA_character_,
      top_chr  = if (length(r$snp_chrs)  > 0L) r$snp_chrs[1]  else NA_character_,
      elapsed_s = r$elapsed,
      chr_distribution = chr_str,
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

# Whole-genome chromosome × method long-format (for figures)
collect_mice_whole_long <- function(
    dir = "results/bench_full/06_realdata_mice_wholegenome") {
  fn <- file.path(dir, "whole_genome.rds")
  if (!file.exists(fn)) return(NULL)
  p <- readRDS(fn)
  if (is.null(p$chr_matrix)) return(NULL)
  long <- as.data.frame(as.table(p$chr_matrix))
  names(long) <- c("chr", "method", "K_hat")
  long$chr <- as.character(long$chr)
  long$method <- as.character(long$method)
  long
}


# 
# Main
# 
if (sys.nframe() == 0L) {
  args <- parse_args()
  skip_figs <- isTRUE(args[["skip-figures"]])
  filter_n  <- isTRUE(args[["filter-n"]])
  n_floor   <- if (!is.null(args[["n-floor"]])) as.integer(args[["n-floor"]]) else 1000L

  cat("=== 99_aggregate ===\n")
  if (filter_n) cat(sprintf("  [filter-n] dropping cells with n < %d\n", n_floor))

  # --- 01 scaling n ---
  df <- collect_sim("results/bench_full/01_scaling_n")
  if (!is.null(df)) {
    if (filter_n) df <- filter_n_floor(df, n_floor)
    s <- summarise_sim(df, c("n", "sigma_g2"))
    write.csv(df, "results/bench_full/01_scaling_n_raw.csv", row.names = FALSE)
    write.csv(s,  "results/bench_full/01_scaling_n_summary.csv", row.names = FALSE)
    cat(sprintf("01_scaling_n:    %d rows -> summary %d cells\n", nrow(df), nrow(s)))
    if (!skip_figs)
      plot_axis(s, "n", "n", "results/bench_full/fig_scaling_n.pdf")
  } else cat("01_scaling_n:    (no data)\n")

  # --- 02 scaling rho ---
  df <- collect_sim("results/bench_full/02_scaling_rho")
  if (!is.null(df)) {
    s <- summarise_sim(df, c("rho", "sigma_g2"))
    write.csv(df, "results/bench_full/02_scaling_rho_raw.csv", row.names = FALSE)
    write.csv(s,  "results/bench_full/02_scaling_rho_summary.csv", row.names = FALSE)
    cat(sprintf("02_scaling_rho:  %d rows -> summary %d cells\n", nrow(df), nrow(s)))
    if (!skip_figs)
      plot_axis(s, "rho", expression(rho), "results/bench_full/fig_scaling_rho.pdf")
  } else cat("02_scaling_rho:  (no data)\n")

  # --- 03 arch ---
  df <- collect_sim("results/bench_full/03_arch")
  if (!is.null(df)) {
    s <- summarise_sim(df, c("signal", "sigma_g2"))
    write.csv(df, "results/bench_full/03_arch_raw.csv", row.names = FALSE)
    write.csv(s,  "results/bench_full/03_arch_summary.csv", row.names = FALSE)
    cat(sprintf("03_arch:         %d rows -> summary %d cells\n", nrow(df), nrow(s)))
    if (!skip_figs)
      plot_axis(s, "signal", "Signal level", "results/bench_full/fig_arch.pdf")
  } else cat("03_arch:         (no data)\n")

  # --- 04 scaling m ---
  df <- collect_sim("results/bench_full/04_scaling_m")
  if (!is.null(df)) {
    s <- summarise_sim(df, c("m", "sigma_g2"))
    write.csv(df, "results/bench_full/04_scaling_m_raw.csv", row.names = FALSE)
    write.csv(s,  "results/bench_full/04_scaling_m_summary.csv", row.names = FALSE)
    cat(sprintf("04_scaling_m:    %d rows -> summary %d cells\n", nrow(df), nrow(s)))
    if (!skip_figs)
      plot_axis(s, "m", "m", "results/bench_full/fig_scaling_m.pdf")
  } else cat("04_scaling_m:    (no data)\n")

  # --- 05 realdata mice (chromosome-wise) ---
  df_mice <- collect_mice()
  if (!is.null(df_mice)) {
    write.csv(df_mice, "results/bench_full/05_realdata_mice_summary.csv",
              row.names = FALSE)
    cat(sprintf("05_realdata_mice: %d rows (chr x method)\n", nrow(df_mice)))
    if (!skip_figs)
      plot_mice_chr(df_mice, "results/bench_full/fig_realdata_mice_chr.pdf")
  } else cat("05_realdata_mice: (no data)\n")

  # --- 06 realdata mice (whole-genome joint) ---
  df_mice_w <- collect_mice_whole()
  if (!is.null(df_mice_w)) {
    write.csv(df_mice_w,
              "results/bench_full/06_realdata_mice_wholegenome_summary.csv",
              row.names = FALSE)
    cat(sprintf("06_realdata_mice_whole: %d methods, K_hat range [%d, %d]\n",
                nrow(df_mice_w),
                min(df_mice_w$K_hat), max(df_mice_w$K_hat)))
    cat("  Target chromosomes per method:\n")
    for (i in seq_len(nrow(df_mice_w))) {
      cat(sprintf("    %-15s K=%4d  chr=[%s]\n",
                  df_mice_w$method[i], df_mice_w$K_hat[i],
                  df_mice_w$target_chromosomes[i]))
    }
    df_long <- collect_mice_whole_long()
    if (!is.null(df_long)) {
      write.csv(df_long,
                "results/bench_full/06_realdata_mice_wholegenome_per_chr.csv",
                row.names = FALSE)
    }
  } else cat("06_realdata_mice_whole: (no data)\n")

  cat("\nDone.\n")
}
