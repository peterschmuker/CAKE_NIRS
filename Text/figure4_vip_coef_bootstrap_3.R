# ============================================================
# Figure 4 -- Bootstrapped correlations of VIP scores and PLS
# regression coefficients, within and between market classes.
# ------------------------------------------------------------

# ============================================================

library(caret)
library(pls)
library(dplyr)
library(ggplot2)

set.seed(123)

# ---- Load data (same as flat_test_100_samples_v3.R) ----
scans <- read.csv("Z:/Grain Spectroscopy/spectral data.csv", check.names = FALSE)
pheno <- read.csv("Z:/Grain Spectroscopy/pheno data.csv",   check.names = FALSE)

drop_cols <- intersect(c(1:7, 163, 164), seq_len(ncol(scans)))
if (length(drop_cols) > 0) scans <- scans[, -drop_cols, drop = FALSE]

# Combine SWW and SWS into "Common"
pheno$CLASS[pheno$CLASS %in% c("SWW", "SWS")] <- "Common"

response_col <- "CAVOL"
eval_classes <- c("Common", "CLUB")

keep_idx <- which(pheno$CLASS %in% eval_classes)
pheno2 <- pheno[keep_idx, , drop = FALSE]
scans2 <- as.data.frame(scans[keep_idx, , drop = FALSE])

ok <- is.finite(pheno2[[response_col]])
pheno2 <- pheno2[ok, , drop = FALSE]
scans2 <- scans2[ok, , drop = FALSE]

# ============================================================
# CONFIG  -- reconstruction assumptions; adjust to match original.
# ============================================================
PREPROC      <- "snv"   # spectral preprocessing for the PLS model: "snv" or "raw".
                        # The caption's "raw / 3-point" refers to the PROFILE, not the
                        # spectra, so this is a separate, unstated choice. SNV is the
                        # project default; switch to "raw" if 0.27 is not reproduced.
NCOMP_METHOD <- "cv"    # ncomp selection: "cv" (caret 10-fold) or "fixed"
NCOMP_FIXED  <- 10      # used only if NCOMP_METHOD == "fixed"
MAX_NCOMP    <- 30      # upper bound for CV tuning
N_BOOT         <- 50    # iterations (Methods: "repeated this procedure 50 times")
SUBSAMPLE_FRAC <- 0.5   # Methods: "randomly selected half of the available samples"
BIN_K        <- 3       # "3-point adjacent" = groups of 3 adjacent wavelengths (~13.5 nm)
out_dir      <- "Z:/CAVOL/TIFFs"

cols <- c("Within club"     = "#E69F00",   # orange
          "Within common"   = "#009E73",   # green
          "Between classes" = "#0072B2")   # blue

# ---- SNV (rowwise), identical to validation script ----
pp_snv <- function(X) {
  X <- as.matrix(X)
  mu  <- rowMeans(X)
  sdv <- apply(X, 1, sd); sdv[sdv == 0 | !is.finite(sdv)] <- 1
  sweep(sweep(X, 1, mu, "-"), 1, sdv, "/")
}
preprocess <- function(X) if (PREPROC == "snv") pp_snv(X) else as.matrix(X)

# ---- 3-point adjacent averaging (non-overlapping groups of BIN_K) ----
bin_profile <- function(v, k = BIN_K) {
  g <- ceiling(seq_along(v) / k)
  as.numeric(tapply(v, g, mean))
}

# ---- caret fit at a fixed ncomp; return VIP + coef (raw and binned) ----
# VIP via caret::varImp (weighted |coef| by per-component SS reduction).
# varImp and coef are returned in predictor (wavelength) order, so club
# and common profiles stay aligned by wavelength for the correlations.
fit_profiles <- function(X, y, ncomp) {
  Xm <- as.matrix(X)
  m  <- train(x = Xm, y = y, method = "pls",
              tuneGrid = data.frame(ncomp = ncomp),
              trControl = trainControl(method = "none"))
  vip <- as.numeric(caret::varImp(m, scale = TRUE)$importance[[1]])   # 0-100, top wl = 100 (per Methods); scaling does not affect Pearson r
  cf  <- as.numeric(coef(m$finalModel, ncomp = ncomp))
  list(vip_raw = vip,             coef_raw = cf,
       vip_bin = bin_profile(vip), coef_bin = bin_profile(cf))
}

# ---- Choose ncomp per class on full data (caret 10-fold CV) ----
pick_ncomp <- function(X, y) {
  if (NCOMP_METHOD == "fixed") return(min(NCOMP_FIXED, ncol(X)))
  m <- train(x = as.matrix(X), y = y, method = "pls",
             tuneGrid = data.frame(ncomp = 1:min(MAX_NCOMP, ncol(X) - 1)),
             trControl = trainControl(method = "cv", number = 10))
  m$bestTune$ncomp
}

# ============================================================
# Full-data models per class
# ============================================================
split_class <- function(cl) {
  idx <- which(pheno2$CLASS == cl)
  list(X = preprocess(scans2[idx, , drop = FALSE]),
       y = pheno2[[response_col]][idx])
}
club   <- split_class("CLUB")
common <- split_class("Common")

nc_club   <- pick_ncomp(club$X,   club$y)
nc_common <- pick_ncomp(common$X, common$y)
cat(sprintf("Selected ncomp -- CLUB: %d, Common: %d (PREPROC = %s)\n",
            nc_club, nc_common, PREPROC))

full_club   <- fit_profiles(club$X,   club$y,   nc_club)
full_common <- fit_profiles(common$X, common$y, nc_common)

# Full-data between-class correlations (dashed reference lines)
ref_lines <- data.frame(
  metric    = c("VIP", "VIP", "Coefficient", "Coefficient"),
  treatment = c("Raw", "3-point", "Raw", "3-point"),
  ref = c(cor(full_club$vip_raw,  full_common$vip_raw),
          cor(full_club$vip_bin,  full_common$vip_bin),
          cor(full_club$coef_raw, full_common$coef_raw),
          cor(full_club$coef_bin, full_common$coef_bin))
)
print(ref_lines)   # check the VIP / Raw value against ~0.27

# ============================================================
# Bootstrap
# ============================================================
boot_rows <- vector("list", N_BOOT)
for (b in seq_len(N_BOOT)) {
  ib_club   <- sample.int(length(club$y),   size = floor(SUBSAMPLE_FRAC * length(club$y)),   replace = FALSE)
  ib_common <- sample.int(length(common$y), size = floor(SUBSAMPLE_FRAC * length(common$y)), replace = FALSE)

  bc <- fit_profiles(club$X[ib_club, , drop = FALSE],     club$y[ib_club],     nc_club)
  bo <- fit_profiles(common$X[ib_common, , drop = FALSE], common$y[ib_common], nc_common)

  boot_rows[[b]] <- data.frame(
    iter = b,
    metric     = rep(c("VIP", "Coefficient"), each = 2, times = 3),
    treatment  = rep(c("Raw", "3-point"),      times = 6),
    comparison = rep(c("Within club", "Within common", "Between classes"), each = 4),
    cor = c(
      # Within club: boot vs full within club
      cor(bc$vip_raw,  full_club$vip_raw),  cor(bc$vip_bin,  full_club$vip_bin),
      cor(bc$coef_raw, full_club$coef_raw), cor(bc$coef_bin, full_club$coef_bin),
      # Within common: boot vs full within common
      cor(bo$vip_raw,  full_common$vip_raw),  cor(bo$vip_bin,  full_common$vip_bin),
      cor(bo$coef_raw, full_common$coef_raw), cor(bo$coef_bin, full_common$coef_bin),
      # Between classes: club-boot vs common-boot at this iteration
      cor(bc$vip_raw,  bo$vip_raw),  cor(bc$vip_bin,  bo$vip_bin),
      cor(bc$coef_raw, bo$coef_raw), cor(bc$coef_bin, bo$coef_bin)
    )
  )
}
boot <- bind_rows(boot_rows)

boot$metric     <- factor(boot$metric,     levels = c("VIP", "Coefficient"))
boot$treatment  <- factor(boot$treatment,  levels = c("Raw", "3-point"))
boot$comparison <- factor(boot$comparison,
                          levels = c("Within club", "Within common", "Between classes"))
ref_lines$metric    <- factor(ref_lines$metric,    levels = levels(boot$metric))
ref_lines$treatment <- factor(ref_lines$treatment, levels = levels(boot$treatment))

# ---- Mean +/- SD per facet and comparison ----
summ <- boot %>%
  group_by(metric, treatment, comparison) %>%
  summarise(mean = mean(cor), sd = sd(cor), .groups = "drop")
print(summ)

# ============================================================
# Plot
# ============================================================
p_fig4 <- ggplot(summ, aes(x = mean, y = comparison, color = comparison)) +
  geom_vline(data = ref_lines, aes(xintercept = ref),
             linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd),
                 height = 0.25, linewidth = 0.6) +
  geom_point(size = 2.4) +
  facet_grid(metric ~ treatment) +
  scale_color_manual(values = cols, name = NULL) +
  scale_y_discrete(limits = rev(levels(summ$comparison))) +
  labs(x = "Bootstrap Pearson correlation", y = NULL) +
  coord_cartesian(xlim = c(min(0, min(summ$mean - summ$sd)), 1)) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text  = element_text(color = "black"),
    legend.position = "bottom",
    strip.background = element_rect(fill = "white", color = "grey60"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )

print(p_fig4)

# ---- 300 dpi TIFF export (LZW) for submission ----
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out_dir, "figure4_vip_bootstrap.tiff"), p_fig4,
       width = 8, height = 6, units = "in", dpi = 300, compression = "lzw")
