# Two publication figures for the CBF-LMM restructured manuscript.
# 5-method roster: CBF-LMM (= validated MS_L_eBIC, MBF evaluation) + comparators.
# Style follows sim/bench_full/make_pub_figures.R.
# Run from the project root.
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(scales) })

OUT <- "results/bench_full/figures_cbf"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

ORDER <- c("MS_L_eBIC", "BSLMM", "BayesR", "SuSiE", "fastlmm")
LABS  <- c(MS_L_eBIC = "CBF-LMM", BSLMM = "BSLMM", BayesR = "BayesR",
           SuSiE = "SuSiE", fastlmm = "FaST-LMM")
COLS  <- c(MS_L_eBIC = "#D55E00", BSLMM = "#CC79A7", BayesR = "#009E73",
           SuSiE = "#0072B2", fastlmm = "#7570B3")
SHAPES <- c(15, 18, 16, 6, 4)
SG_LABS <- c(`0` = "sigma[g]^2 == 0", `0.5` = "sigma[g]^2 == 0.5")

theme_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.3, colour = "grey92"),
          strip.background = element_rect(fill = "grey95", colour = "grey80"),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom", legend.title = element_blank(),
          axis.title = element_text(face = "bold"))
}
load_axis <- function(name) {
  df <- read.csv(file.path("results/bench_full", paste0(name, "_summary.csv")))
  df <- df[df$method %in% ORDER, , drop = FALSE]
  df$method <- factor(df$method, levels = ORDER)
  df$sg <- factor(df$sigma_g2, levels = c(0, 0.5), labels = SG_LABS[c("0","0.5")])
  df
}

## Figure 1: mean F1 vs n (n >= 1000), faceted by regime -----------------------
df_n <- load_axis("01_scaling_n"); df_n <- df_n[df_n$n >= 1000, ]
p1 <- ggplot(df_n, aes(n, mean_f1, colour = method, shape = method, group = method)) +
  geom_line(linewidth = 0.6, alpha = 0.9) + geom_point(size = 2.8) +
  facet_wrap(~ sg, nrow = 1, labeller = label_parsed) +
  scale_colour_manual(values = COLS, labels = LABS, drop = FALSE) +
  scale_shape_manual(values = SHAPES, labels = LABS, drop = FALSE) +
  scale_x_continuous(trans = "log10", breaks = unique(df_n$n)) +
  coord_cartesian(ylim = c(0, 1.05)) +
  labs(x = expression("Sample size " * n), y = expression("Mean " * F[1])) +
  theme_pub() + guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))
ggsave(file.path(OUT, "fig_cbf_axis_n_F1.pdf"), p1, width = 9.5, height = 4.2, units = "in")
message("Wrote fig_cbf_axis_n_F1.pdf")

## Figure 2: recall-FDR scatter over all 16 cells ------------------------------
all_df <- bind_rows(
  load_axis("01_scaling_n")   %>% mutate(axis = "n"),
  load_axis("02_scaling_rho") %>% mutate(axis = "rho"),
  load_axis("03_arch")        %>% mutate(axis = "signal"),
  load_axis("04_scaling_m")   %>% mutate(axis = "m")
) %>% mutate(mean_fdr = 1 - mean_precision)
p2 <- ggplot(all_df, aes(mean_fdr, mean_recall, colour = method, shape = method)) +
  geom_point(size = 2.4, alpha = 0.75) +
  facet_wrap(~ sg, nrow = 1, labeller = label_parsed) +
  scale_colour_manual(values = COLS, labels = LABS, drop = FALSE) +
  scale_shape_manual(values = SHAPES, labels = LABS, drop = FALSE) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3, color = "grey60",
              linetype = "dashed") +
  labs(x = "FDR", y = "Recall") + theme_pub() +
  guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))
ggsave(file.path(OUT, "fig_cbf_recall_fdr.pdf"), p2, width = 9.5, height = 5, units = "in")
message("Wrote fig_cbf_recall_fdr.pdf")
