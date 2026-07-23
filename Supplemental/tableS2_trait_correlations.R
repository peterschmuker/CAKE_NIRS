# ==============================================================
# Supplemental Table S2: Pearson correlations between sponge cake
#           volume (CAVOL) and secondary quality traits, by market
#           class, with sample sizes and p-values.
# Manuscript: Near-infrared prediction of sponge cake volume in
#             common and club soft wheats
# Repository: https://github.com/peterschmuker/CAKE_NIRS
#
# Input:  pheno_data_shared.csv  (CLASS coded Common / Club)
# Output: tables/tableS2_trait_correlations.csv
# No randomness; purely descriptive.
# ==============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# ---- Paths (edit if your files live elsewhere) ----
data_dir   <- "."
pheno_file <- "pheno_data_shared.csv"
tab_dir    <- "tables"
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

pheno <- read.csv(file.path(data_dir, pheno_file), check.names = FALSE)
pheno <- pheno %>% filter(CLASS %in% c("Common", "Club"))

# Secondary traits to correlate with cake volume (manuscript row order)
traits <- c("FPROT", "FSDS", "FSR_S", "FSR_W", "FSR_C", "BKFYELD")
trait_labels <- c(
  FPROT   = "Flour protein",
  FSDS    = "Sedimentation volume",
  FSR_S   = "Sucrose SRC",
  FSR_W   = "Water SRC",
  FSR_C   = "Calcium SRC",
  BKFYELD = "Break flour yield"
)
ordered_labels <- unname(trait_labels[traits])

# Guard against any -9 missing-value sentinel (none expected in shared data)
pheno_clean <- pheno %>%
  mutate(across(all_of(c("CAVOL", traits)), ~ na_if(., -9)))

# Correlation, p-value, and n for each trait x class (complete pairs only)
cor_results <- lapply(traits, function(trait) {
  pheno_clean %>%
    group_by(CLASS) %>%
    summarise(
      Trait = trait_labels[trait],
      n     = sum(!is.na(CAVOL) & !is.na(.data[[trait]])),
      r     = if (n > 2) cor.test(CAVOL, .data[[trait]], use = "complete.obs")$estimate else NA_real_,
      p     = if (n > 2) cor.test(CAVOL, .data[[trait]], use = "complete.obs")$p.value  else NA_real_,
      .groups = "drop"
    )
}) %>%
  bind_rows() %>%
  mutate(r = round(r, 2), p = round(p, 3),
         Trait = factor(Trait, levels = ordered_labels)) %>%
  select(CLASS, Trait, n, r, p) %>%
  arrange(Trait, CLASS)

# Wide format for the supplemental table
cor_wide <- cor_results %>%
  pivot_wider(
    names_from  = CLASS,
    values_from = c(n, r, p),
    names_glue  = "{CLASS}_{.value}"
  ) %>%
  select(Trait, starts_with("Common"), starts_with("Club")) %>%
  arrange(Trait)

print(cor_wide)

out_csv <- file.path(tab_dir, "tableS2_trait_correlations.csv")
write.csv(cor_wide, out_csv, row.names = FALSE)
cat("\nSaved:", out_csv, "\n")
