# ==============================================================
# Figure 1: Distribution of sponge cake volume for common, club,
#           and pooled (All) samples (violin plots with boxplots).
# Manuscript: Near-infrared prediction of sponge cake volume in
#             common and club soft wheats
# Repository: https://github.com/peterschmuker/CAKE_NIRS
#
# Input:  pheno_data_shared.csv  (CLASS coded Common / Club)
# Output: figures/figure1_cake_volume_violin.tiff  (500 dpi, LZW)
# No randomness; purely descriptive.
# ==============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# ---- Paths (edit if your files live elsewhere) ----
data_dir   <- "."
pheno_file <- "pheno_data_shared.csv"
fig_dir    <- "figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1) Load phenotype data ----
pheno <- read.csv(file.path(data_dir, pheno_file), check.names = FALSE)

# ---- 2) Column for cake volume; normalize class labels ----
cake_col <- "CAVOL"
if (!"CLASS" %in% names(pheno)) stop("Column 'CLASS' not found in pheno.")
if (!cake_col %in% names(pheno)) stop(sprintf("Column '%s' not found in pheno.", cake_col))

# CLASS is coded Common / Club.
pheno2 <- pheno %>%
  mutate(CLASS = as.character(CLASS)) %>%
  filter(CLASS %in% c("Common", "Club")) %>%
  filter(is.finite(.data[[cake_col]])) %>%
  mutate(CLASS = factor(CLASS, levels = c("Common", "Club")))

# ---- 3) Plotting frame with an "All" category ----
df_plot <- bind_rows(
  pheno2 %>% transmute(Group = CLASS, CakeVol = .data[[cake_col]]),
  pheno2 %>% transmute(Group = factor("All", levels = c("Common", "Club", "All")),
                       CakeVol = .data[[cake_col]])
) %>%
  mutate(Group = factor(Group, levels = c("Common", "Club", "All")))

# ---- 4) Stats (Mean, SD, CV) by group ----
stats_by_group <- df_plot %>%
  group_by(Group) %>%
  summarise(
    n       = sum(is.finite(CakeVol)),
    mean    = mean(CakeVol, na.rm = TRUE),
    sd      = sd(CakeVol,   na.rm = TRUE),
    cv_perc = 100 * sd / mean,
    .groups = "drop"
  )

print(
  stats_by_group %>%
    transmute(Group, n, Mean = round(mean, 1), SD = round(sd, 1), `CV (%)` = round(cv_perc, 1))
)

# ---- 5) Fixed label heights ----
fixed_y   <- 1500
gap_above <- 35

labels_mean <- stats_by_group %>%
  transmute(Group, label_y = fixed_y, lab = sprintf("Mean = %.1f cm\u00B3", mean))

labels_cv <- stats_by_group %>%
  transmute(Group, label_y = fixed_y + gap_above, lab = sprintf("CV = %.1f%%", cv_perc))

y_min <- min(df_plot$CakeVol, na.rm = TRUE)
y_max <- max(df_plot$CakeVol, na.rm = TRUE)
upper_limit <- max(y_max, fixed_y + gap_above + 20)

# ---- 6) Violin plot with labels above ----
p_violin <- ggplot(df_plot, aes(x = Group, y = CakeVol, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.85, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", alpha = 0.95) +
  stat_summary(fun = mean, geom = "point",
               shape = 21, size = 2.8, color = "black", fill = "white", stroke = 0.6) +
  geom_label(
    data = labels_mean,
    aes(x = Group, y = label_y, label = lab),
    inherit.aes = FALSE,
    label.size = 0, fill = "white", color = "black", size = 3.6
  ) +
  geom_label(
    data = labels_cv,
    aes(x = Group, y = label_y, label = lab),
    inherit.aes = FALSE,
    label.size = 0, fill = "white", color = "black", size = 3.6
  ) +
  scale_fill_manual(values = c("Common" = "#E57373", "Club" = "#4FC3F7", "All" = "gray70")) +
  labs(
    x = "Market class",
    y = "Sponge Cake Volume (cm\u00B3)"
  ) +
  coord_cartesian(ylim = c(y_min, upper_limit), clip = "off") +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.margin = margin(12, 40, 12, 40)
  )

print(p_violin)

# ---- 7) Save (500 dpi TIFF, LZW) ----
ggsave(file.path(fig_dir, "figure1_cake_volume_violin.tiff"),
       p_violin, width = 6, height = 6, units = "in", dpi = 500, compression = "lzw")
