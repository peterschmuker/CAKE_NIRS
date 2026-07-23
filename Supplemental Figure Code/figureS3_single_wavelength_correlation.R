# ==============================================================
# Supplemental Figure S3: Pearson correlation between each single
#           wavelength and sponge cake volume, computed within the
#           common class, within the club class, and pooled.
# Manuscript: Near-infrared prediction of sponge cake volume in
#             common and club soft wheats
# Repository: https://github.com/peterschmuker/CAKE_NIRS
#
# Inputs: spectral_data_shared.csv  (Bcode + 155 wavelength columns)
#         pheno_data_shared.csv     (CLASS coded Common / Club; CAVOL)
# Output: figures/figureS3_single_wavelength_correlation.tiff (500 dpi, LZW)
#         single_wavelength_correlation.csv (per-wavelength r values)
#
# Spectra are SNV scatter-corrected (rowwise), then each wavelength
# is centered and scaled, matching the calibration preprocessing.
# Descriptive per-wavelength Pearson correlation; no model, no
# resampling. Raw spectra are dominated by a shared baseline/scatter
# component, so SNV is used to expose per-wavelength structure.
# ==============================================================

suppressPackageStartupMessages({ library(ggplot2) })

# ---- Paths (edit if your files live elsewhere) ----
data_dir      <- "."
spectral_file <- "spectral_data_shared.csv"
pheno_file    <- "pheno_data_shared.csv"
fig_dir       <- "figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

response_col <- "CAVOL"
eval_classes <- c("Common", "Club")

# ---- Load and align the two files by Bcode ----
scans <- read.csv(file.path(data_dir, spectral_file), check.names = FALSE)
pheno <- read.csv(file.path(data_dir, pheno_file),    check.names = FALSE)
stopifnot("Bcode" %in% names(scans), "Bcode" %in% names(pheno))

# order spectra to match pheno rows via the shared key
scans <- scans[match(pheno$Bcode, scans$Bcode), , drop = FALSE]

# wavelength matrix = every spectral column except the Bcode key
wl_cols <- setdiff(names(scans), "Bcode")
scans   <- scans[, wl_cols, drop = FALSE]

# ---- keep classes of interest, complete response ----
keep_idx <- which(pheno$CLASS %in% eval_classes)
pheno2 <- pheno[keep_idx, , drop = FALSE]
scans2 <- as.data.frame(scans[keep_idx, , drop = FALSE])
if (!response_col %in% names(pheno2)) stop("response_col not found in pheno")
ok <- is.finite(pheno2[[response_col]])
pheno2 <- pheno2[ok, , drop = FALSE]; scans2 <- scans2[ok, , drop = FALSE]
stopifnot(all(vapply(scans2, is.numeric, logical(1))))
wl <- suppressWarnings(as.numeric(gsub("[^0-9.]+", "", colnames(scans2))))
if (anyNA(wl)) wl <- seq_len(ncol(scans2))

# ---- SNV (rowwise), matching the calibration preprocessing ----
pp_snv <- function(X) {
  X <- as.matrix(X); mu <- rowMeans(X); s <- apply(X, 1, sd)
  s[s == 0 | !is.finite(s)] <- 1
  sweep(sweep(X, 1, mu, "-"), 1, s, "/")
}

# ---- groups: two within-class, plus the pooled set ----
groups <- list(
  Common = which(pheno2$CLASS == "Common"),
  Club   = which(pheno2$CLASS == "Club"),
  Pooled = seq_len(nrow(pheno2))
)

# ---- per-wavelength correlation per group (SNV, then center + scale) ----
per_grp <- list(); summ <- list()
for (g in names(groups)) {
  idx <- groups[[g]]
  Xc  <- pp_snv(scans2[idx, , drop = FALSE])   # SNV scatter correction
  Xc  <- scale(Xc)                             # center + scale wavelengths
  yc  <- pheno2[[response_col]][idx]
  r   <- suppressWarnings(apply(Xc, 2, function(col) cor(col, yc, use = "complete.obs")))
  per_grp[[g]] <- data.frame(group = g, wavelength = wl, r = as.numeric(r))
  i <- which.max(abs(r))
  summ[[g]] <- data.frame(group = g, n = length(idx),
    max_abs_r       = round(max(abs(r), na.rm = TRUE), 3),
    peak_wavelength = wl[i],
    median_abs_r    = round(median(abs(r), na.rm = TRUE), 3),
    pct_above_0.30  = round(mean(abs(r) > 0.30, na.rm = TRUE) * 100, 1))
}
cor_df     <- do.call(rbind, per_grp)
summary_df <- do.call(rbind, summ); rownames(summary_df) <- NULL

cat("\n=== Single-wavelength correlation with CAVOL (SNV; within-class and pooled) ===\n\n")
print(summary_df, row.names = FALSE)
cat("\nFor reference, the PLS Matching models reached r = 0.61 (Common) and r = 0.69 (Club).\n")
write.csv(cor_df, "single_wavelength_correlation.csv", row.names = FALSE)

# ---- plot (no title; the figure caption carries the description) ----
cor_df$group <- factor(cor_df$group, levels = c("Common", "Club", "Pooled"))
p <- ggplot(cor_df, aes(wavelength, r, color = group)) +
  geom_hline(yintercept = 0, color = "grey75") +
  geom_hline(yintercept = c(-0.3, 0.3), linetype = "dashed", color = "grey60") +
  geom_line(linewidth = 0.8) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(x = "Wavelength (nm)", y = "Pearson r with cake volume", color = NULL) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")
print(p)

# ---- Save (500 dpi TIFF, LZW) ----
ggsave(file.path(fig_dir, "figureS3_single_wavelength_correlation.tiff"),
       p, width = 8, height = 5, units = "in", dpi = 500, compression = "lzw")
cat("\nSaved figure to", file.path(fig_dir, "figureS3_single_wavelength_correlation.tiff"), "\n")
