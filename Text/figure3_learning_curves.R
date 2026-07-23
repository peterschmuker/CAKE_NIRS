# ==============================================================
# Figure 3: Learning curves of average predicted-vs-reference
#           correlation as a function of training size per class,
#           faceted by evaluation class (rows) and model (cols),
#           colored by training source (Matching/Opposing/Both).
# Manuscript: Near-infrared prediction of sponge cake volume in
#             common and club soft wheats
# Repository: https://github.com/peterschmuker/CAKE_NIRS
#
# Input:  growing_sample_size_validation.csv
#         (size-structured validation output; one row per
#          replicate x eval_class x training_source x s_per_class,
#          with columns cor_pls, cor_enet, cor_local)
# Output: figures/figure3_learning_curves.tiff  (500 dpi, LZW)
# ==============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(scales); library(grid)
})

# ---- Paths (edit if your files live elsewhere) ----
data_dir     <- "."
results_file <- "growing_sample_size_validation.csv"
fig_dir      <- "figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

results_df <- read.csv(file.path(data_dir, results_file))

mlevels  <- c("LASSO", "LOCAL PLS", "PLS")
slevels  <- c("Matching", "Opposing", "Both")
src_cols <- c("Matching" = "#4E79A7", "Opposing" = "#9C9C9C", "Both" = "#E15759")
src_lty  <- c("Matching" = "solid",   "Opposing" = "dotted",  "Both" = "dashed")

results_long <- results_df %>%
  pivot_longer(c(cor_pls, cor_enet, cor_local),
               names_to = "model", values_to = "correlation") %>%
  mutate(model = dplyr::recode(model, cor_pls = "PLS", cor_enet = "LASSO", cor_local = "LOCAL PLS"),
         model = factor(model, levels = mlevels),
         training_source = factor(training_source, levels = slevels),
         eval_class = factor(eval_class, levels = c("Common", "Club")))

# mean +/- 95% CI
mean_ci <- function(x) {
  m <- mean(x, na.rm = TRUE); n <- sum(!is.na(x))
  ci <- 1.96 * sd(x, na.rm = TRUE) / sqrt(n)
  data.frame(y = m, ymin = m - ci, ymax = m + ci)
}

p_curves <- ggplot(results_long, aes(s_per_class, correlation,
                         color = training_source, fill = training_source,
                         linetype = training_source)) +
  stat_summary(fun.data = mean_ci, geom = "ribbon", alpha = 0.12, color = NA) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(fun = mean, geom = "point", size = 1.8) +
  facet_grid(eval_class ~ model, scales = "free_y") +
  scale_color_manual(values = src_cols) +
  scale_fill_manual(values = src_cols) +
  scale_linetype_manual(values = src_lty) +
  labs(x = "Training size per class (s)", y = "Average Correlation",
       color = "Training source", fill = "Training source", linetype = "Training source") +
  theme_minimal(base_size = 13) +
  theme(panel.border    = element_rect(color = "black", linewidth = 0.6, fill = NA),
        panel.spacing.y = unit(0.8, "lines"),
        panel.spacing.x = unit(0.8, "lines"),
        strip.text      = element_text(face = "bold"),
        legend.position = "bottom")

print(p_curves)

# ---- Save (500 dpi TIFF, LZW) ----
ggsave(file.path(fig_dir, "figure3_learning_curves.tiff"),
       p_curves, width = 10, height = 6.5, units = "in", dpi = 500, compression = "lzw")
