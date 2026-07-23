# ==============================================================
# Supplemental Figure S1: PCA of SNV-transformed spectra for
#           common vs club, with PERMANOVA and PERMDISP tests.
# Manuscript: Near-infrared prediction of sponge cake volume in
#             common and club soft wheats
# Repository: https://github.com/peterschmuker/CAKE_NIRS
#
# Inputs: spectral_data_shared.csv  (Bcode + 155 wavelength columns)
#         pheno_data_shared.csv     (CLASS coded Common / Club)
# Output: figures/figureS1_pca_anova_snv.tiff  (500 dpi, LZW)
#
# SNV scatter correction is applied with the waves package
# (waves::pretreat_spectra, pretreatment = 2), which calls
# prospectr::standardNormalVariate (rowwise SNV). PERMANOVA/PERMDISP
# use 500 permutations; set.seed(42) for reproducibility.
# ==============================================================

suppressPackageStartupMessages({
  library(vegan)     # adonis2 (PERMANOVA), betadisper (PERMDISP)
  library(ggplot2)
  library(dplyr)
  library(waves)     # SNV via pretreat_spectra()
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

# ---- 2) Keep target classes ----
eval_classes <- c("Common", "Club")
keep_idx <- which(pheno$CLASS %in% eval_classes)
pheno2 <- pheno[keep_idx, , drop = FALSE]
scans2 <- scans[keep_idx, , drop = FALSE]   # Bcode (metadata) + X wavelength columns
class  <- factor(pheno2$CLASS, levels = c("Common", "Club"))

# ---- 3) SNV via waves ----
# waves treats non-"X" columns (here Bcode) as metadata and applies
# rowwise SNV to the X-named spectral columns; pretreatment = 2 = SNV.
snv_df <- waves::pretreat_spectra(df = scans2, pretreatment = 2)
stopifnot(nrow(snv_df) == nrow(pheno2))     # no rows dropped -> still aligned
X_snv <- as.matrix(snv_df[, grepl("^X", names(snv_df)), drop = FALSE])
Xs    <- scale(X_snv)                        # column-standardize for PCA and distance

# ---- 4) PCA (PC1 vs PC2) with class centroids ----
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
ggsave(file.path(fig_dir, "figureS1_pca_anova_snv.tiff"),
       p_pca, width = 6, height = 6, units = "in", dpi = 500, compression = "lzw")

# ---- 5) PERMANOVA (location differences) on SNV ----
D <- dist(X_snv, method = "euclidean")
set.seed(42)
adon <- adonis2(D ~ class, permutations = 500, by = "terms")
cat("\n=== PERMANOVA (adonis2) on SNV spectra ===\n")
print(adon)

# ---- 6) PERMDISP (dispersion homogeneity) on SNV ----
bd <- betadisper(D, group = class, type = "centroid")
bd_anova <- anova(bd)
bd_tuk   <- tryCatch(TukeyHSD(bd), error = function(e) NULL)

cat("\n=== PERMDISP (betadisper) on SNV spectra ===\n")
print(bd_anova)
ss_groups <- bd_anova$`Sum Sq`[1]
ss_resid  <- bd_anova$`Sum Sq`[2]
eta2_permdisp <- ss_groups / (ss_groups + ss_resid)
cat(sprintf("\nPERMDISP effect size (eta^2 ~ R^2): %.4f\n", eta2_permdisp))

if (!is.null(bd_tuk)) {
  cat("\n--- Pairwise differences in dispersion (Tukey HSD) ---\n")
  print(bd_tuk)
}

# ---- 7) Quick summary ----
alpha  <- 0.05
perm_p <- adon$`Pr(>F)`[1]
disp_p <- bd_anova$`Pr(>F)`[1]

cat("\n=== Quick summary (SNV) ===\n")
cat(sprintf("PERMANOVA p-value (class effect): %.4g\n", perm_p))
cat(sprintf("PERMDISP  p-value (dispersion):   %.4g\n", disp_p))

if (perm_p < alpha && disp_p >= alpha) {
  cat("Interpretation: classes differ in multivariate location (centroids) with similar dispersion (SNV data).\n")
} else if (perm_p < alpha && disp_p < alpha) {
  cat("Interpretation: classes differ, but dispersion heterogeneity is present; results may reflect both location and spread (SNV data).\n")
} else if (perm_p >= alpha) {
  cat("Interpretation: no strong evidence of centroid separation at the chosen alpha (SNV data).\n")
}
