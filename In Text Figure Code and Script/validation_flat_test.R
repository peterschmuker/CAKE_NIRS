# ==============================================================
# Validation: preprocessing sweep for Table 1 and Supplemental Figure S2.
# For each market class, replicate (25), training source (Matching/
# Opposing/Both) and preprocessing method (raw, snv, d1, d2, snvd1,
# snvd2, msc), fits PLS, LASSO, and LOCAL PLS on all remaining rows and
# scores a 100-sample held-out test set (cor, RMSE, bias, slope).
# Manuscript: Near-infrared prediction of sponge cake volume in
#             common and club soft wheats
# Repository: https://github.com/peterschmuker/CAKE_NIRS
#
# Output: preprocessing_results.csv
#   backs Supplemental Figure S2 (via figureS2_preprocessing_heatmap.R)
#   and Table 1 (the best-preprocessing row per model x scenario).
# set.seed(123) plus per-replicate seeds make the run reproducible.
# ------------------------------------------------------------

# ---- Packages ----
# install.packages(c("caret","pls","glmnet","resemble","dplyr","tidyr","ggplot2","purrr","prospectr"))
library(caret)
library(pls)
library(glmnet)
library(resemble)
library(dplyr)
library(tidyr)
library(prospectr)
library(waves)      # SNV via pretreat_spectra()

set.seed(123)

# ---- Paths (edit if your files live elsewhere) ----
data_dir      <- "."
spectral_file <- "spectral_data_shared.csv"
pheno_file    <- "pheno_data_shared.csv"

# ---- Load data and align by Bcode ----
scans <- read.csv(file.path(data_dir, spectral_file), check.names = FALSE)
pheno <- read.csv(file.path(data_dir, pheno_file),    check.names = FALSE)
stopifnot("Bcode" %in% names(scans), "Bcode" %in% names(pheno))
scans <- scans[match(pheno$Bcode, scans$Bcode), , drop = FALSE]
scans <- scans[, grepl("^X", names(scans)), drop = FALSE]   # wavelength matrix (X columns)

# CLASS is coded Common / Club in the shared data.

# ---- CONFIG ----
response_col <- "CAVOL"             # <-- change to "FYELD" or other target as needed
eval_classes <- c("Common", "Club") # the two market classes
test_size    <- 100                 # test size per eval class
n_reps       <- 25                  # number of replicates for robustness
k_grid       <- seq(50, 500, by = 50)   # neighborhood sizes (k) tuned by NNv (capped at n_train-1)

# Savitzky–Golay settings for derivatives (adjust as needed; w must be odd)
sg_window <- 11   # typical: 9, 11, 15, 21 (odd)
sg_poly   <- 2    # typical: 2 or 3

# ---- resemble controls ----
my_control <- mbl_control(validation_type = c("NNv"))   # nearest-neighbor validation
my_waplsr  <- local_fit_wapls(min_pls_c = 2, max_pls_c = 30)

# ---- Keep only rows in classes of interest and ensure complete cases ----
keep_idx <- which(pheno$CLASS %in% eval_classes)
pheno2   <- pheno[keep_idx, , drop = FALSE]
scans2   <- as.data.frame(scans[keep_idx, , drop = FALSE])

# Response column must exist and be numeric
if (!response_col %in% names(pheno2)) stop("response_col not found in pheno.")
if (!is.numeric(pheno2[[response_col]])) {
  stop(paste0("'", response_col, "' must be numeric."))
}

# Drop rows with missing response
ok <- is.finite(pheno2[[response_col]])
pheno2 <- pheno2[ok, , drop = FALSE]
scans2 <- scans2[ok, , drop = FALSE]

# Ensure spectra are numeric (caret wants DF; resemble wants matrix)
if (!all(vapply(scans2, is.numeric, logical(1)))) {
  stop("Some spectral columns are not numeric. Please convert them before modeling.")
}

# ============================================================
# Spectral preprocessing: waves for SNV and the SG derivatives; prospectr for
# SG-of-SNV and MSC (neither is a single waves pretreatment).
# waves codes: 2 = SNV, 11 = SG 1st deriv (w=11, p=2), 13 = SG 2nd deriv.
# ============================================================
apply_preproc <- function(name, Xtr, Xte) {
  wv <- function(X, code) as.matrix(waves::pretreat_spectra(as.data.frame(X, check.names = FALSE), pretreatment = code))
  sg <- function(X, m)    prospectr::savitzkyGolay(as.matrix(X), m = m, p = sg_poly, w = sg_window)
  switch(tolower(name),
    raw   = list(Xtr = as.matrix(Xtr),    Xte = as.matrix(Xte)),
    snv   = list(Xtr = wv(Xtr, 2),        Xte = wv(Xte, 2)),
    d1    = list(Xtr = wv(Xtr, 11),       Xte = wv(Xte, 11)),
    d2    = list(Xtr = wv(Xtr, 13),       Xte = wv(Xte, 13)),
    snvd1 = list(Xtr = sg(wv(Xtr, 2), 1), Xte = sg(wv(Xte, 2), 1)),
    snvd2 = list(Xtr = sg(wv(Xtr, 2), 2), Xte = sg(wv(Xte, 2), 2)),
    msc   = { r <- colMeans(as.matrix(Xtr)); list(Xtr = prospectr::msc(as.matrix(Xtr), r), Xte = prospectr::msc(as.matrix(Xte), r)) },
    stop("unknown preprocessing: ", name)
  )
}

# Drop columns that are NA/Inf or near-zero variance in TRAIN; keep alignment for TEST
finalize_columns <- function(Xtr, Xte) {
  Xtr <- as.matrix(Xtr); Xte <- as.matrix(Xte)
  bad_na <- which(!is.finite(colSums(Xtr)) | !is.finite(colSums(Xte)))
  if (length(bad_na) > 0) {
    keep <- setdiff(seq_len(ncol(Xtr)), bad_na)
    Xtr <- Xtr[, keep, drop = FALSE]
    Xte <- Xte[, keep, drop = FALSE]
  }
  nzv_idx <- caret::nearZeroVar(Xtr, saveMetrics = FALSE)
  if (length(nzv_idx) > 0) {
    keep <- setdiff(seq_len(ncol(Xtr)), nzv_idx)
    Xtr <- Xtr[, keep, drop = FALSE]
    Xte <- Xte[, keep, drop = FALSE]
  }
  list(Xtr = Xtr, Xte = Xte)
}

# ============================================================
# Robust helpers (unchanged)
# ============================================================

clean_pred <- function(x) {
  x_num <- as.numeric(unlist(x, use.names = FALSE))
  x_num[!is.finite(x_num)] <- NA_real_
  x_num
}

safe_cor <- function(y, yhat) {
  y <- as.numeric(y); yhat <- as.numeric(yhat)
  keep <- is.finite(y) & is.finite(yhat)
  if (sum(keep) < 2L) return(NA_real_)
  y <- y[keep]; yhat <- yhat[keep]
  if (sd(yhat) < .Machine$double.eps) return(NA_real_)
  suppressWarnings(cor(y, yhat))
}

rmse <- function(y, yhat) {
  y <- as.numeric(y); yhat <- as.numeric(yhat)
  keep <- is.finite(y) & is.finite(yhat)
  if (sum(keep) == 0L) return(NA_real_)
  y <- y[keep]; yhat <- yhat[keep]
  sqrt(mean((y - yhat)^2))
}

bias <- function(y, yhat) {
  y <- as.numeric(y); yhat <- as.numeric(yhat)
  keep <- is.finite(y) & is.finite(yhat)
  if (sum(keep) == 0L) return(NA_real_)
  y <- y[keep]; yhat <- yhat[keep]
  mean(yhat - y)
}

slope <- function(y, yhat) {
  y <- as.numeric(y); yhat <- as.numeric(yhat)
  keep <- is.finite(y) & is.finite(yhat)
  if (sum(keep) < 2L) return(NA_real_)
  y <- y[keep]; yhat <- yhat[keep]
  v <- var(yhat)
  if (is.na(v) || v < .Machine$double.eps) return(NA_real_)
  cov(y, yhat) / v
}

metricize <- function(y, yhat) {
  c(
    cor   = safe_cor(y, yhat),
    rmse  = rmse(y, yhat),
    bias  = bias(y, yhat),
    slope = slope(y, yhat)
  )
}

DEBUG_DIAG <- FALSE
diag_metrics <- function(name, y, yhat) {
  if (!DEBUG_DIAG) return(invisible(NULL))
  y <- as.numeric(y); yhat <- as.numeric(yhat)
  keep <- is.finite(y) & is.finite(yhat)
  cat(sprintf(
    "%s: n=%d | finite_pairs=%d | sd(pred)=%.6f | anyNA(pred)=%s | anyNA(obs)=%s\n",
    name, length(y), sum(keep), sd(yhat[keep]), anyNA(yhat), anyNA(y)
  ))
}

# ============================================================
# Results container (add 'preproc')
# ============================================================
results_df <- data.frame(
  eval_class            = character(),
  replicate             = integer(),
  training_source       = character(), # Matching / Opposing / Both
  preproc               = character(),
  n_train               = integer(),
  # PLS
  best_ncomp_caret      = integer(),
  cor_caret             = numeric(),
  rmse_caret            = numeric(),
  bias_caret            = numeric(),
  slope_caret           = numeric(),
  # LASSO
  best_lambda_lasso     = numeric(),
  n_nonzero_lasso       = integer(),
  cor_lasso             = numeric(),
  rmse_lasso            = numeric(),
  bias_lasso            = numeric(),
  slope_lasso           = numeric(),
  # Local
  cor_local             = numeric(),
  rmse_local            = numeric(),
  bias_local            = numeric(),
  slope_local           = numeric(),
  stringsAsFactors      = FALSE
)

# Define the preprocessing methods to evaluate
preproc_methods <- c("raw", "snv", "d1", "d2", "snvd1", "snvd2", "msc")

# ============================================================
# Main loop
# ============================================================
for (eval_class in eval_classes) {
  opp_class <- setdiff(eval_classes, eval_class)
  
  # Indices by class
  idx_eval <- which(pheno2$CLASS == eval_class)
  idx_opp  <- which(pheno2$CLASS == opp_class)
  
  if (length(idx_eval) < test_size) {
    stop("Not enough rows in '", eval_class, "' to draw ", test_size, " test samples.")
  }
  
  message(sprintf("\n==== Evaluating class '%s' (%d candidates) ====\n",
                  eval_class, length(idx_eval)))
  
  for (rep_i in seq_len(n_reps)) {
    set.seed(1000 + rep_i)
    test_idx <- sample(idx_eval, size = test_size, replace = FALSE)
    
    # Test sets (kept raw; preprocessing applied inside the inner loop)
    test_x_raw_df <- scans2[test_idx, , drop = FALSE]
    test_x_raw_m  <- as.matrix(test_x_raw_df)
    test_y        <- pheno2[[response_col]][test_idx]
    
    # Training pools
    pool_eval <- setdiff(idx_eval, test_idx)  # same-class remaining
    pool_opp  <- idx_opp                      # all opposing class
    
    # Three scenarios: Matching, Opposing, Both (use ALL remaining in each)
    scenarios <- list(
      Matching = pool_eval,
      Opposing = pool_opp,
      Both     = union(pool_eval, pool_opp)   # all rows except the 100 test
    )
    
    for (scenario_name in names(scenarios)) {
      tr_idx <- scenarios[[scenario_name]]
      if (length(tr_idx) == 0) next
      
      train_x_raw_df <- scans2[tr_idx, , drop = FALSE]
      train_x_raw_m  <- as.matrix(train_x_raw_df)
      train_y        <- pheno2[[response_col]][tr_idx]
      
      for (pp in preproc_methods) {
        
        # -------- Apply preprocessing (train-derived where needed) --------
        # (Derivatives & SNV operate per row; MSC uses training mean as reference.)
        ppd <- apply_preproc(pp, train_x_raw_m, test_x_raw_m)
        aligned <- finalize_columns(ppd$Xtr, ppd$Xte)
        Xtr_m <- aligned$Xtr
        Xte_m <- aligned$Xte
        
        # If all predictors were removed, skip this combo
        if (ncol(Xtr_m) < 2L) {
          warning(sprintf("Skipping: %s / %s / %s : <2 predictors after preprocessing.",
                          eval_class, scenario_name, pp))
          next
        }
        
        # For caret models we need data.frames with matching column names
        colnames(Xtr_m) <- make.names(colnames(Xtr_m), unique = TRUE)
        colnames(Xte_m) <- colnames(Xtr_m)
        
        train_x_df <- as.data.frame(Xtr_m)
        test_x_df  <- as.data.frame(Xte_m)
        
        # ==== Global PLS (caret) ====
        caret_cor <- caret_rmse <- caret_bias <- caret_slope <- NA_real_
        best_ncomp <- NA_integer_
        caret_ok <- TRUE
        
        pls_fit <- tryCatch({
          set.seed(11 + rep_i)
          train(
            x = train_x_df, y = train_y,
            method     = "pls",
            preProcess = c("center", "scale"),
            tuneLength = 30,
            trControl  = trainControl(method = "cv", number = 10),
            metric     = "RMSE"
          )
        }, error = function(e) { caret_ok <<- FALSE; NULL })
        
        if (caret_ok && !is.null(pls_fit)) {
          preds_caret <- clean_pred(predict(pls_fit, newdata = test_x_df))
          m <- metricize(test_y, preds_caret)
          caret_cor   <- m["cor"]
          caret_rmse  <- m["rmse"]
          caret_bias  <- m["bias"]
          caret_slope <- m["slope"]
          best_ncomp  <- pls_fit$bestTune$ncomp
          diag_metrics("PLS", test_y, preds_caret)
        }
        
        # ==== Global LASSO (caret + glmnet) ====
        lasso_cor <- lasso_rmse <- lasso_bias <- lasso_slope <- NA_real_
        best_lambda <- NA_real_
        n_nonzero   <- NA_integer_
        lasso_ok <- TRUE
        
        lambda_grid <- 10^seq(-6, 2, length.out = 120)
        
        lasso_fit <- tryCatch({
          set.seed(22 + rep_i)
          train(
            x = train_x_df, y = train_y,
            method     = "glmnet",
            preProcess = c("center", "scale"),
            tuneGrid   = expand.grid(alpha = 1, lambda = lambda_grid),  # alpha=1 => LASSO
            trControl  = trainControl(method = "cv", number = 10),
            metric     = "RMSE"
          )
        }, error = function(e) { lasso_ok <<- FALSE; NULL })
        
        if (lasso_ok && !is.null(lasso_fit)) {
          preds_lasso <- clean_pred(predict(lasso_fit, newdata = test_x_df))
          m2 <- metricize(test_y, preds_lasso)
          lasso_cor   <- m2["cor"]
          lasso_rmse  <- m2["rmse"]
          lasso_bias  <- m2["bias"]
          lasso_slope <- m2["slope"]
          best_lambda <- lasso_fit$bestTune$lambda
          
          coefs <- tryCatch({
            as.matrix(coef(lasso_fit$finalModel, s = best_lambda))
          }, error = function(e) NULL)
          if (!is.null(coefs)) n_nonzero <- sum(abs(coefs) > 0) - 1L
          
          diag_metrics("LASSO", test_y, preds_lasso)
        }
        
        # ==== Local WAPLS (resemble::mbl) ====
        local_cor  <- local_rmse <- local_bias <- local_slope <- NA_real_
        
        n_tr <- length(train_y)
        k_use <- k_grid[k_grid <= (n_tr - 1L)]          # tune k over 50-500, capped at n_tr-1
        if (length(k_use) == 0L) k_use <- max(2L, n_tr - 1L)
        
        if (n_tr >= 3L) {
          local_fit <- tryCatch({
            mbl(
              Xr = Xtr_m,
              Yr = train_y,
              Xu = Xte_m,
              k  = k_use,
              method      = my_waplsr,
              diss_method = "pls",
              diss_usage  = "predictors",
              control     = my_control,   # NNv selects the best neighborhood size
              scale       = TRUE
            )
          }, error = function(e) NULL)
          
          if (!is.null(local_fit)) {
            preds_local <- tryCatch({
              # choose the neighborhood size (k) with the lowest NNv RMSE
              vr <- local_fit$validation_results$nearest_neighbor_validation
              bk <- vr$k[which.min(vr$rmse)]
              pm <- as.matrix(get_predictions(local_fit))
              cn <- colnames(pm)
              ci <- if (!is.null(cn) && as.character(bk) %in% cn) which(cn == as.character(bk)) else
                    if (!is.null(cn) && paste0("k=", bk) %in% cn) which(cn == paste0("k=", bk)) else
                    which(vr$k == bk)
              clean_pred(pm[, ci[1]])
            }, error = function(e) rep(NA_real_, length(test_y)))
            m3 <- metricize(test_y, preds_local)
            local_cor   <- m3["cor"]
            local_rmse  <- m3["rmse"]
            local_bias  <- m3["bias"]
            local_slope <- m3["slope"]
            diag_metrics("LOCAL", test_y, preds_local)
          }
        }
        
        # ---- Save results ----
        results_df <- rbind(
          results_df,
          data.frame(
            eval_class            = eval_class,
            replicate             = rep_i,
            training_source       = scenario_name,
            preproc               = pp,
            n_train               = length(tr_idx),
            # PLS
            best_ncomp_caret      = best_ncomp,
            cor_caret             = as.numeric(caret_cor),
            rmse_caret            = as.numeric(caret_rmse),
            bias_caret            = as.numeric(caret_bias),
            slope_caret           = as.numeric(caret_slope),
            # LASSO
            best_lambda_lasso     = as.numeric(best_lambda),
            n_nonzero_lasso       = as.integer(n_nonzero),
            cor_lasso             = as.numeric(lasso_cor),
            rmse_lasso            = as.numeric(lasso_rmse),
            bias_lasso            = as.numeric(lasso_bias),
            slope_lasso           = as.numeric(lasso_slope),
            # Local
            cor_local             = as.numeric(local_cor),
            rmse_local            = as.numeric(local_rmse),
            bias_local            = as.numeric(local_bias),
            slope_local           = as.numeric(local_slope),
            stringsAsFactors      = FALSE
          )
        )
        
        cat(sprintf(
          "Class=%-6s Rep=%02d | %-8s | preproc=%-6s | n_train=%4d | PLS (r=%.3f rmse=%.3f bias=%.2f slope=%.3f) | LASSO (r=%.3f rmse=%.3f bias=%.2f slope=%.3f; λ=%.4g, nz=%s) | Local (r=%.3f rmse=%.3f bias=%.2f slope=%.3f)\n",
          eval_class, rep_i, scenario_name, pp, length(tr_idx),
          ifelse(is.na(caret_cor), NA, round(caret_cor,3)),
          ifelse(is.na(caret_rmse), NA, round(caret_rmse,3)),
          ifelse(is.na(caret_bias), NA, round(caret_bias,2)),
          ifelse(is.na(caret_slope), NA, round(caret_slope,3)),
          ifelse(is.na(lasso_cor), NA, round(lasso_cor,3)),
          ifelse(is.na(lasso_rmse), NA, round(lasso_rmse,3)),
          ifelse(is.na(lasso_bias), NA, round(lasso_bias,2)),
          ifelse(is.na(lasso_slope), NA, round(lasso_slope,3)),
          ifelse(is.na(best_lambda), NA, signif(best_lambda, 4)),
          ifelse(is.na(n_nonzero), NA, n_nonzero),
          ifelse(is.na(local_cor), NA, round(local_cor,3)),
          ifelse(is.na(local_rmse), NA, round(local_rmse,3)),
          ifelse(is.na(local_bias), NA, round(local_bias,2)),
          ifelse(is.na(local_slope), NA, round(local_slope,3))
        ))
      } # end preproc
    }   # end scenario
  }     # end replicate
}       # end eval_class

write.csv(results_df, file.path(data_dir, "preprocessing_results.csv"), row.names = FALSE)
cat("\nSaved: preprocessing_results.csv\n")


# ---- Summaries ----
print(head(results_df, 12))

summary_df <- results_df %>%
  group_by(eval_class, training_source, preproc) %>%
  summarise(
    mean_n_train      = mean(n_train),
    # PLS
    mean_cor_caret    = mean(cor_caret, na.rm = TRUE),
    sd_cor_caret      = sd(cor_caret, na.rm = TRUE),
    mean_rmse_caret   = mean(rmse_caret, na.rm = TRUE),
    sd_rmse_caret     = sd(rmse_caret, na.rm = TRUE),
    mean_bias_caret   = mean(bias_caret, na.rm = TRUE),
    mean_slope_caret  = mean(slope_caret, na.rm = TRUE),
    # LASSO
    mean_cor_lasso    = mean(cor_lasso, na.rm = TRUE),
    sd_cor_lasso      = sd(cor_lasso, na.rm = TRUE),
    mean_rmse_lasso   = mean(rmse_lasso, na.rm = TRUE),
    sd_rmse_lasso     = sd(rmse_lasso, na.rm = TRUE),
    mean_bias_lasso   = mean(bias_lasso, na.rm = TRUE),
    mean_slope_lasso  = mean(slope_lasso, na.rm = TRUE),
    mean_lambda_lasso = mean(best_lambda_lasso, na.rm = TRUE),
    mean_nz_lasso     = mean(n_nonzero_lasso, na.rm = TRUE),
    # Local
    mean_cor_local    = mean(cor_local, na.rm = TRUE),
    sd_cor_local      = sd(cor_local, na.rm = TRUE),
    mean_rmse_local   = mean(rmse_local, na.rm = TRUE),
    sd_rmse_local     = sd(rmse_local, na.rm = TRUE),
    mean_bias_local   = mean(bias_local, na.rm = TRUE),
    mean_slope_local  = mean(slope_local, na.rm = TRUE),
    n_reps            = dplyr::n(),
    .groups = "drop"
  )

print(summary_df)
