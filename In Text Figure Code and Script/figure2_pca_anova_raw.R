# ==============================================================
# Figure 2: PCA of unprocessed (raw) spectra for common vs club,
#           with PERMANOVA and PERMDISP tests of class structure.
# Manuscript: Near-infrared prediction of sponge cake volume in
#             common and club soft wheats
# Repository: https://github.com/peterschmuker/CAKE_NIRS
#
# Inputs: spectral_data_shared.csv  (Bcode + 155 wavelength columns)
#         pheno_data_shared.csv     (CLASS coded Common / Club)
# Output: figures/figure2_pca_anova_raw.tiff  (500 dpi, LZW)
# PERMANOVA/PERMDISP use 500 permutations; set.seed(42) for reproducibility.
# ==============================================================

suppressPackageStartupMessages({
  library(vegan)     # adonis2 (PERMANOVA), betadisper (PERMDISP)
  library(ggplot2)
  library(dplyr)
})

# ---- Paths (edit if your files live elsewhere) ----
data_dir      <- "."
spectral_file <- "spectral_data_shared.csv"
pheno_file    <- "pheno_data_shared.csv"
fig_dir       <- "figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1) Load and align by Bcode ----
scans <- read.csv(file.path(data_dir, spectral_file), check.names = FALSE)
pheno <- read.csv(file.path(data_dir, pheno_file),    check.names = FALSE)
stopifnot("Bcode" %in% names(scans), "Bcode" %in% names(pheno))
scans <- scans[match(pheno$Bcode, scans$Bcode), , drop = FALSE]

# ---- 2) Keep target classes; wavelength matrix = every X-column ----
eval_classes <- c("Common", "Club")
keep_idx <- which(pheno$CLASS %in% eval_classes)
pheno2 <- pheno[keep_idx, , drop = FALSE]
X <- scans[keep_idx, grepl("^X", names(scans)), drop = FALSE]
X <- as.data.frame(X)
stopifnot(all(vapply(X, is.numeric, logical(1))))
class <- factor(pheno2$CLASS, levels = c("Common", "Club"))

# ---- 3) PCA (PC1 vs PC2) with class centroids; raw spectra as the basis ----
Xs  <- scale(as.matrix(X))          # center + scale by wavelength
pca <- prcomp(Xs, center = FALSE, scale. = FALSE)

var_exp <- (pca$sdev^2) / sum(pca$sdev^2) * 100
ve1 <- round(var_exp[1], 1)
ve2 <- round(var_exp[2], 1)

pdat <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], class = class)
centers <- pdat %>%
  group_by(class) %>%
  summarize(PC1 = mean(PC1, na.rm = TRUE), PC2 = mean(PC2, na.rm = TRUE), .groups = "drop")

p_pca <- ggplot(pdat, aes(PC1, PC2, color = class)) +
  geom_point(alpha = 0.70, size = 1.9) +
  stat_ellipse(type = "norm", linetype = 2, linewidth = 0.6, show.legend = FALSE) +
  geom_text(data = centers,
            aes(x = PC1, y = PC2, label = as.character(class)),
            inherit.aes = FALSE, color = "black", fontface = "bold", size = 4.2) +
  coord_equal() +
  theme_minimal(base_size = 13) +
  theme(aspect.ratio = 1, legend.position = "bottom",
        panel.grid.minor = element_blank(), plot.margin = margin(8, 8, 8, 8)) +
  labs(x = sprintf("PC1 (%.1f%% variance)", ve1),
       y = sprintf("PC2 (%.1f%% variance)", ve2), color = "Class")

print(p_pca)
ggsave(file.path(fig_dir, "figure2_pca_anova_raw.tiff"),
       p_pca, width = 6, height = 6, units = "in", dpi = 500, compression = "lzw")

# ---- 4) PERMANOVA (location differences) on raw ----
D <- dist(Xs, method = "euclidean")
set.seed(42)
adon <- adonis2(D ~ class, permutations = 500, by = "terms")
cat("\n=== PERMANOVA (adonis2) results (RAW) ===\n")
print(adon)

# ---- 5) PERMDISP (dispersion homogeneity) on raw ----
bd <- betadisper(D, group = class, type = "centroid")
bd_anova <- anova(bd)
bd_tuk   <- tryCatch(TukeyHSD(bd), error = function(e) NULL)

cat("\n=== PERMDISP (betadisper) ANOVA on distances to centroid (RAW) ===\n")
print(bd_anova)
if (!is.null(bd_tuk)) {
  cat("\n--- Pairwise differences in dispersion (Tukey HSD) ---\n")
  print(bd_tuk)
}

ss_groups <- bd_anova$`Sum Sq`[1]
ss_resid  <- bd_anova$`Sum Sq`[2]
eta2_permdisp <- ss_groups / (ss_groups + ss_resid)
cat(sprintf("\nPERMDISP effect size (eta^2 ~ R^2): %.4f\n", eta2_permdisp))

# ---- 6) Quick summary ----
alpha  <- 0.05
perm_p <- adon$`Pr(>F)`[1]
disp_p <- bd_anova$`Pr(>F)`[1]

cat("\n=== Quick summary ===\n")
cat(sprintf("PERMANOVA p-value (class effect): %.4g\n", perm_p))
cat(sprintf("PERMDISP  p-value (dispersion):   %.4g\n", disp_p))

if (perm_p < alpha && disp_p >= alpha) {
  cat("Interpretation: classes differ in multivariate location (centroids) with similar dispersion.\n")
} else if (perm_p < alpha && disp_p < alpha) {
  cat("Interpretation: classes differ, but dispersion heterogeneity is present; results may reflect both location and spread.\n")
} else if (perm_p >= alpha) {
  cat("Interpretation: no strong evidence of centroid separation at the chosen alpha.\n")
}
