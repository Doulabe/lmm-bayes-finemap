# ============================================================
# make_pub_figures.R
# Publication-grade figures from results-2/bench_full v2 data.
#
# Outputs to results/bench_full/ (and theory/ for manuscript inclusion):
#   - fig_pub_axis_n.pdf       (F1 vs n × σg² × method)
#   - fig_pub_axis_rho.pdf     (F1 vs ρ × σg² × method)
#   - fig_pub_axis_arch.pdf    (F1 vs signal × σg² × method)
#   - fig_pub_axis_m.pdf       (F1 vs m × σg² × method)
#   - fig_pub_recall_fdr.pdf   (recall–FDR scatter, all 8 methods)
#   - fig_pub_mice_K.pdf       (K_hat per chromosome × method)
#
# Style: Wong palette, 7-method panel, faceted by σg², prior-style colour.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(scales)
})

OUT <- "results/bench_full"

# -- Method labels: stable ordering, custom colours ------------------
METHOD_ORDER <- c("JS_L_eBIC", "MS_L_eBIC",
                   "JS_L_JP99", "MS_L_JP99",
                   "BayesR", "BSLMM",
                   "SuSiE", "fastlmm")
METHOD_LABS <- c(
  JS_L_eBIC     = "JS-LMM (eBIC)",
  MS_L_eBIC     = "MS-LMM (eBIC)",
  JS_L_JP99     = "JS-LMM (JP99)",
  MS_L_JP99     = "MS-LMM (JP99)",
  BayesR        = "BayesR",
  BSLMM         = "BSLMM",
  SuSiE         = "SuSiE",
  fastlmm       = "FaST-LMM"
)
METHOD_COLS <- c(
  JS_L_eBIC     = "#D55E00",   # vermillion (primary)
  MS_L_eBIC     = "#E69F00",   # orange
  JS_L_JP99     = "#A04000",   # dark vermillion
  MS_L_JP99     = "#A66B00",   # dark orange
  BayesR        = "#009E73",   # green
  BSLMM         = "#CC79A7",   # pink
  SuSiE         = "#0072B2",   # blue
  fastlmm       = "#7570B3"    # purple — LMM-frequentist baseline
)
# Use plotmath strings parseable via label_parsed for clean math rendering
# in PDF (avoids Unicode-rendering issues with subscript g).
SG_LABS <- c(`0`   = "sigma[g]^2 == 0",
              `0.5` = "sigma[g]^2 == 0.5")

# Drop n=500 cells in axis-n figures (undersized regime).
N_FLOOR <- 1000L

theme_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.3, colour = "grey92"),
          strip.background = element_rect(fill = "grey95", colour = "grey80"),
          strip.text       = element_text(face = "bold"),
          legend.position  = "bottom",
          legend.title     = element_blank(),
          axis.title       = element_text(face = "bold"))
}

# -- Common helper: load a summary CSV, format columns ---------------
load_axis <- function(name) {
  df <- read.csv(file.path(OUT, paste0(name, "_summary.csv")))
  df <- df[df$method %in% METHOD_ORDER, , drop = FALSE]   # drop varbvs/BL/REML rows
  df$method <- factor(df$method, levels = METHOD_ORDER)
  df$sg <- factor(df$sigma_g2, levels = c(0, 0.5),
                   labels = SG_LABS[c("0","0.5")])
  df
}

# Single-metric plot per axis: faceted by sigma_g2, no error bars.
# `metric` ∈ {"F1", "Precision"}; emits ONE PDF per (axis, metric).
# 8 shapes for the 8 reported methods
SHAPE_VEC <- c(15, 17, 0, 2, 16, 18, 6, 4)

plot_axis_metric <- function(df, x_var, label_x, metric, out_name,
                                log_x = FALSE) {
  df <- df %>%
    mutate(method = factor(method, levels = METHOD_ORDER)) %>%
    filter(method != "REML_stepwise")
  y_col   <- if (metric == "F1") "mean_f1" else "mean_precision"
  y_label <- if (metric == "F1") expression("Mean " * F[1])
             else expression("Mean Precision")
  p <- ggplot(df, aes(x = !!sym(x_var), y = !!sym(y_col),
                       colour = method, shape = method, group = method)) +
    geom_line(linewidth = 0.6, alpha = 0.9) +
    geom_point(size = 2.8) +
    facet_wrap(~ sg, nrow = 1, labeller = label_parsed) +
    scale_colour_manual(values = METHOD_COLS, labels = METHOD_LABS,
                         drop = FALSE) +
    scale_shape_manual(values = SHAPE_VEC,
                        labels = METHOD_LABS, drop = FALSE) +
    coord_cartesian(ylim = c(0, 1.05)) +
    labs(x = label_x, y = y_label) +
    theme_pub() +
    guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))
  if (log_x) p <- p + scale_x_continuous(trans = "log10",
                                            breaks = unique(df[[x_var]]))
  ggsave(file.path(OUT, out_name), p, width = 9.5, height = 4.2, units = "in")
  message("Wrote ", out_name)
  p
}

# Convenience: emit BOTH F1 and Precision figures for one axis
plot_axis <- function(df, x_var, label_x, axis_tag, log_x = FALSE) {
  plot_axis_metric(df, x_var, label_x, "F1",
    paste0("fig_pub_axis_", axis_tag, "_F1.pdf"), log_x = log_x)
  plot_axis_metric(df, x_var, label_x, "Precision",
    paste0("fig_pub_axis_", axis_tag, "_Precision.pdf"), log_x = log_x)
}


# --- Figures 1a, 1b: scaling in n (F1, Precision) -----------------
# Skip n < N_FLOOR cells (undersized regime, framework not competitive there).
df_n <- load_axis("01_scaling_n")
df_n <- df_n[df_n$n >= N_FLOOR, , drop = FALSE]
plot_axis(df_n, "n", expression("Sample size " * n),
          axis_tag = "n", log_x = TRUE)


# --- Figures 2a, 2b: scaling in rho (F1, Precision) ---------------
df_rho <- load_axis("02_scaling_rho")
plot_axis(df_rho, "rho", expression("Block-AR(1) correlation " * rho),
          axis_tag = "rho")


# --- Figures 3a, 3b: signal architecture (F1 / Precision separate) ---
df_arch <- load_axis("03_arch")
df_arch$signal <- factor(df_arch$signal, levels = c("weak","medium","strong"))
df_arch$x <- as.integer(df_arch$signal)

plot_arch_metric <- function(metric, out_name) {
  y_col   <- if (metric == "F1") "mean_f1" else "mean_precision"
  y_label <- if (metric == "F1") expression("Mean " * F[1])
             else expression("Mean Precision")
  p <- ggplot(df_arch %>% filter(method != "REML_stepwise"),
               aes(x = x, y = !!sym(y_col), colour = method,
                    shape = method, group = method)) +
    geom_line(linewidth = 0.6, alpha = 0.9) +
    geom_point(size = 2.8) +
    facet_wrap(~ sg, nrow = 1, labeller = label_parsed) +
    scale_x_continuous(breaks = 1:3, labels = c("weak","medium","strong")) +
    scale_colour_manual(values = METHOD_COLS, labels = METHOD_LABS, drop = FALSE) +
    scale_shape_manual(values = SHAPE_VEC,
                        labels = METHOD_LABS, drop = FALSE) +
    coord_cartesian(ylim = c(0, 1.05)) +
    labs(x = "Signal strength", y = y_label) +
    theme_pub() +
    guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))
  ggsave(file.path(OUT, out_name), p, width = 9.5, height = 4.2,
          units = "in")
  message("Wrote ", out_name)
}
plot_arch_metric("F1",        "fig_pub_axis_arch_F1.pdf")
plot_arch_metric("Precision", "fig_pub_axis_arch_Precision.pdf")


# --- Figures 4a, 4b: scaling in m (F1, Precision) -----------------
df_m <- load_axis("04_scaling_m")
plot_axis(df_m, "m", expression("Number of SNPs " * m),
          axis_tag = "m", log_x = TRUE)


# --- Figure 5: recall–FDR scatter (combined view) ------------------
all_df <- bind_rows(
  df_n   %>% mutate(axis = "n",      cell_id = paste0("n=",      n)),
  df_rho %>% mutate(axis = "rho",    cell_id = paste0("rho=",    rho)),
  df_arch %>% mutate(axis = "signal",cell_id = signal),
  df_m   %>% mutate(axis = "m",      cell_id = paste0("m=",      m))
) %>% filter(method != "REML_stepwise") %>%
  mutate(mean_fdr = 1 - mean_precision)

p_rf <- ggplot(all_df, aes(x = mean_fdr, y = mean_recall,
                              colour = method, shape = method)) +
  geom_point(size = 2.4, alpha = 0.75) +
  facet_wrap(~ sg, nrow = 1, labeller = label_parsed) +
  scale_colour_manual(values = METHOD_COLS, labels = METHOD_LABS,
                       drop = FALSE) +
  scale_shape_manual(values = SHAPE_VEC,
                      labels = METHOD_LABS, drop = FALSE) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3, color = "grey60",
               linetype = "dashed") +
  labs(x = "FDR", y = "Recall") +
  theme_pub() +
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))
ggsave(file.path(OUT, "fig_pub_recall_fdr.pdf"), p_rf,
       width = 9.5, height = 5, units = "in")
message("Wrote fig_pub_recall_fdr.pdf")

# --- Figure 5b: recall–FDR scatter per axis (4 panels × 2 sg cols) -
AXIS_LABS <- c(n = "n (sample size)",
                rho = "rho (LD strength)",
                signal = "signal architecture",
                m = "m (SNPs)")
all_df$axis_lab <- factor(AXIS_LABS[all_df$axis],
                              levels = AXIS_LABS[c("n","rho","signal","m")])
p_rf_axis <- ggplot(all_df, aes(x = mean_fdr, y = mean_recall,
                                    colour = method, shape = method)) +
  geom_point(size = 2.0, alpha = 0.78) +
  facet_grid(axis_lab ~ sg,
              labeller = labeller(sg = label_parsed,
                                    axis_lab = label_value)) +
  scale_colour_manual(values = METHOD_COLS, labels = METHOD_LABS,
                       drop = FALSE) +
  scale_shape_manual(values = SHAPE_VEC,
                      labels = METHOD_LABS, drop = FALSE) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3, color = "grey60",
               linetype = "dashed") +
  labs(x = "FDR", y = "Recall") +
  theme_pub() +
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))
ggsave(file.path(OUT, "fig_pub_recall_fdr_per_axis.pdf"), p_rf_axis,
       width = 9.5, height = 12, units = "in")
message("Wrote fig_pub_recall_fdr_per_axis.pdf")


# --- Figure 6: mice K_hat per chromosome × method ------------------
mice_csv <- file.path(OUT, "05_realdata_mice_summary.csv")
if (!file.exists(mice_csv)) {
  message("(skipped fig_pub_mice_K: 05_realdata_mice_summary.csv not found)")
} else {
df_mice <- read.csv(mice_csv)
df_mice <- df_mice[df_mice$method %in% METHOD_ORDER, , drop = FALSE]
df_mice$method <- factor(df_mice$method, levels = METHOD_ORDER)
df_mice$chr <- factor(df_mice$chr, levels = c(as.character(1:19), "X"))
df_mice <- df_mice %>% filter(method != "REML_stepwise")

p_mice <- ggplot(df_mice, aes(x = chr, y = pmin(K_hat, 30), fill = method)) +
  geom_col(position = position_dodge2(preserve = "single", padding = 0.05),
            width = 0.85, colour = "black", linewidth = 0.15) +
  scale_fill_manual(values = METHOD_COLS, labels = METHOD_LABS, drop = FALSE) +
  labs(x = "Chromosome", y = expression(hat(K) * " (capped at 30)"),
        caption = "BSLMM K capped at 30 for visibility (true totals up to 548 on chr 11)") +
  theme_pub() +
  theme(plot.caption = element_text(size = 8, color = "grey40")) +
  guides(fill = guide_legend(nrow = 2))
ggsave(file.path(OUT, "fig_pub_mice_K.pdf"), p_mice,
       width = 11, height = 4.5, units = "in")
message("Wrote fig_pub_mice_K.pdf")
}  # end of mice_csv if-exists block


# --- Figure 7: Whole-genome target chromosomes per method ----------
WHOLE_FN <- file.path(OUT, "06_realdata_mice_wholegenome/whole_genome.rds")
if (file.exists(WHOLE_FN)) {
  p_whole <- readRDS(WHOLE_FN)
  long <- as.data.frame(as.table(p_whole$chr_matrix))
  names(long) <- c("chr", "method", "K_hat")
  long$chr <- factor(long$chr, levels = c(as.character(1:19), "X"))
  long <- long[long$method %in% METHOD_ORDER, , drop = FALSE]
  long$method <- factor(long$method, levels = METHOD_ORDER)
  long <- long %>% filter(method != "REML_stepwise") %>%
    mutate(K_hat = as.integer(K_hat))

  # Heatmap: chromosome × method, fill = K_hat (capped for visibility)
  long$K_capped <- pmin(long$K_hat, 30)
  long$K_label  <- ifelse(long$K_hat > 0, as.character(long$K_hat), "")

  p_whole_heat <- ggplot(long, aes(x = chr, y = method, fill = K_capped)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = K_label), size = 2.6, colour = "black") +
    scale_fill_gradient(low = "white", high = "#D55E00",
                         limits = c(0, 30), oob = scales::squish,
                         name = expression(hat(K) * " (capped 30)")) +
    scale_y_discrete(labels = METHOD_LABS) +
    labs(x = "Chromosome", y = NULL,
          caption = "Whole-genome joint analysis on mouse BMI; cell label is the
                     true selection count (some BSLMM cells exceed 30 — colour saturates).") +
    theme_pub() +
    theme(plot.caption = element_text(size = 8, color = "grey40"),
           legend.position = "right")
  ggsave(file.path(OUT, "fig_pub_mice_whole_heatmap.pdf"), p_whole_heat,
         width = 11, height = 4.5, units = "in")
  message("Wrote fig_pub_mice_whole_heatmap.pdf")

  # Consensus bar: chromosomes hit by ≥ k methods
  n_total <- dplyr::n_distinct(long$method)
  consensus <- long %>%
    group_by(chr) %>%
    summarise(n_methods = sum(K_hat > 0), .groups = "drop")
  p_consensus <- ggplot(consensus, aes(x = chr, y = n_methods)) +
    geom_col(fill = "#D55E00", colour = "black", linewidth = 0.2) +
    geom_text(aes(label = n_methods), vjust = -0.3, size = 3) +
    scale_y_continuous(breaks = 0:n_total, limits = c(0, n_total + 0.5),
                        expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Chromosome",
          y = sprintf("# methods (out of %d) selecting at least one SNP", n_total),
          caption = sprintf("Consensus across the %d benchmark methods.", n_total)) +
    theme_pub() +
    theme(plot.caption = element_text(size = 8, color = "grey40"))
  ggsave(file.path(OUT, "fig_pub_mice_whole_consensus.pdf"), p_consensus,
         width = 10, height = 3.8, units = "in")
  message("Wrote fig_pub_mice_whole_consensus.pdf")
} else {
  message("(skipped fig_pub_mice_whole_*: whole_genome.rds not found yet)")
}


# --- Copy all to the manuscript figures/ dir for direct LaTeX inclusion --
FIG_DEST <- "results/bench_full/figures"
dir.create(FIG_DEST, showWarnings = FALSE, recursive = TRUE)
fig_files <- c(
  # 4 axes × 2 metrics = 8 figures
  "fig_pub_axis_n_F1.pdf",     "fig_pub_axis_n_Precision.pdf",
  "fig_pub_axis_rho_F1.pdf",   "fig_pub_axis_rho_Precision.pdf",
  "fig_pub_axis_arch_F1.pdf",  "fig_pub_axis_arch_Precision.pdf",
  "fig_pub_axis_m_F1.pdf",     "fig_pub_axis_m_Precision.pdf",
  # other figures
  "fig_pub_recall_fdr.pdf",
  "fig_pub_recall_fdr_per_axis.pdf",
  "fig_pub_mice_K.pdf",
  "fig_pub_mice_whole_heatmap.pdf",
  "fig_pub_mice_whole_consensus.pdf"
)
for (f in fig_files) {
  src <- file.path(OUT, f)
  if (file.exists(src))
    file.copy(src, file.path(FIG_DEST, f), overwrite = TRUE)
}
message("\nAll available figures copied to ", FIG_DEST)
