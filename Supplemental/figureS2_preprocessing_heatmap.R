# ==============================================================
# Supplemental Figure S2: Effect of spectral preprocessing on
#           calibration RMSE, shown as a signed delta-RMSE heatmap
#           (blue = lower RMSE than raw; red = higher), faceted by
#           evaluation class (rows) and training source (cols).
# Manuscript: Near-infrared prediction of sponge cake volume in
#             common and club soft wheats
# Repository: https://github.com/peterschmuker/CAKE_NIRS
#
# Input:  preprocessing_results.csv
#         (per-replicate RMSE by eval_class x training_source x
#          preproc, with columns rmse_caret [PLS], rmse_lasso,
#          rmse_local, and a "raw" preproc level as the reference)
# Output: figures/figureS2_preprocessing_heatmap.tiff (500 dpi, LZW)
#
# dRMSE = mean(preproc) - mean(raw) within
#         (class x model x training x replicate), then averaged
#         across replicates. Set label_with_abs = TRUE to print
#         absolute mean RMSE in the cells instead of dRMSE.
# ==============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(forcats); library(ggplot2)
})

# ---- Paths (edit if your files live elsewhere) ----
data_dir     <- "."
preproc_file <- "preprocessing_results.csv"
fig_dir      <- "figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Preprocessing label lookup (raw -> display)
preproc_labels <- c(
  "raw" = "Raw", "snv" = "SNV", "msc" = "MSC", "d1" = "D1",
  "snvd1" = "SNVD1", "d2" = "D2", "snvd2" = "SNVD2"
)

# Print absolute mean RMSE in cells instead of dRMSE
label_with_abs <- FALSE

# ---- 1) Load + build long RMSE ----
df <- read.csv(file.path(data_dir, preproc_file), check.names = FALSE)

df <- df %>%
  mutate(
    eval_class      = factor(eval_class, levels = c("Common", "Club")),
    training_source = factor(training_source, levels = c("Matching", "Opposing", "Both")),
    preproc         = fct_relevel(factor(preproc), "raw", "snv", "msc", "d1", "snvd1", "d2", "snvd2")
  )

long_rmse <- df %>%
  transmute(eval_class, training_source, preproc, replicate,
            rmse_pls = rmse_caret, rmse_lasso = rmse_lasso, rmse_local = rmse_local) %>%
  pivot_longer(starts_with("rmse_"), names_to = "model", values_to = "rmse") %>%
  mutate(model = recode(model,
                        rmse_pls = "PLS", rmse_lasso = "LASSO", rmse_local = "LOCAL PLS")) %>%
  filter(is.finite(rmse))

# Per-replicate dRMSE vs raw
delta_df <- long_rmse %>%
  group_by(eval_class, model, training_source, replicate) %>%
  mutate(rmse_raw = rmse[preproc == "raw"][1]) %>%
  ungroup() %>%
  filter(!is.na(rmse_raw), preproc != "raw") %>%
  mutate(delta_rmse = rmse - rmse_raw)

# ---- 2) Summaries ----
delta_sum <- delta_df %>%
  group_by(eval_class, training_source, model, preproc) %>%
  summarise(mean_delta = mean(delta_rmse, na.rm = TRUE),
            sd = sd(delta_rmse, na.rm = TRUE),
            n = dplyr::n(), .groups = "drop")

# Absolute mean RMSE per cell (for optional labeling)
rmse_abs <- long_rmse %>%
  filter(preproc != "raw") %>%
  group_by(eval_class, training_source, model, preproc) %>%
  summarise(mean_rmse = mean(rmse, na.rm = TRUE),
            sd_rmse   = sd(rmse, na.rm = TRUE), .groups = "drop")

# ---- 3) Assemble heatmap data ----
heat <- delta_sum %>%
  left_join(rmse_abs, by = c("eval_class", "training_source", "model", "preproc")) %>%
  mutate(
    eval_class      = factor(eval_class, levels = c("Common", "Club")),
    training_source = factor(training_source, levels = c("Matching", "Opposing", "Both")),
    preproc         = recode(as.character(preproc), !!!preproc_labels),
    preproc         = fct_relevel(preproc, "SNV", "MSC", "D1", "SNVD1", "D2", "SNVD2"),
    preproc         = fct_rev(preproc),
    model           = factor(model, levels = c("LASSO", "LOCAL PLS", "PLS")),
    mean_label      = if (label_with_abs) sprintf("%.1f", mean_rmse)
                      else sprintf("%+.1f", mean_delta),
    sd_label        = if (label_with_abs) sprintf("\u00B1%.1f", sd_rmse)
                      else sprintf("\u00B1%.1f", sd)
  )

# Symmetric color limit centered at 0
lim <- max(abs(heat$mean_delta), na.rm = TRUE)

# ---- 4) Plot ----
p_heat <- ggplot(heat, aes(x = model, y = preproc, fill = mean_delta)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = mean_label), position = position_nudge(y = 0.14),
            size = 2.7) +
  geom_text(aes(label = sd_label), position = position_nudge(y = -0.16),
            size = 2.0, color = "grey30") +
  facet_grid(eval_class ~ training_source) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-lim, lim),
    name = "Mean \u0394RMSE\nvs raw (mL)"
  ) +
  labs(x = NULL, y = "Preprocessing") +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid  = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    strip.text  = element_text(face = "bold"),
    legend.position = "right"
  )

print(p_heat)

# ---- 5) Save (500 dpi TIFF, LZW) ----
ggsave(file.path(fig_dir, "figureS2_preprocessing_heatmap.tiff"),
       p_heat, width = 10, height = 6, units = "in", dpi = 500, compression = "lzw")
