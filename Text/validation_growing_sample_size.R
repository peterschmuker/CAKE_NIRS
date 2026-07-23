## ==============================================================
## Validation: learning curves for Figure 3 (growing training size).
## For each market class, replicate (25), per-class training size s
## (100..700), and training source (Matching/Opposing/Both), fits PLS,
## LASSO/ENet, and LOCAL PLS using the best preprocessing per scenario
## and records the held-out test-set correlation.
## Models: PLS (caret), LASSO/ENet (caret+glmnet), LOCAL PLS (resemble::mbl)
## Manuscript: Near-infrared prediction of sponge cake volume in
##             common and club soft wheats
## Repository: https://github.com/peterschmuker/CAKE_NIRS
##
## Output: growing_sample_size_validation.csv
##   backs Figure 3 (via figure3_learning_curves.R).
## set.seed(123) plus per-replicate seeds make the run reproducible.
## ==============================================================

# ---- Paths (edit if your files live elsewhere) ----
data_dir      <- "."
spectral_file <- "spectral_data_shared.csv"
pheno_file    <- "pheno_data_shared.csv"

# ---- Data & packages ----
scans <- read.csv(file.path(data_dir, spectral_file), check.names = FALSE)
pheno <- read.csv(file.path(data_dir, pheno_file),    check.names = FALSE)
stopifnot("Bcode" %in% names(scans), "Bcode" %in% names(pheno))
scans <- scans[match(pheno$Bcode, scans$Bcode), , drop = FALSE]
scans <- scans[, grepl("^X", names(scans)), drop = FALSE]   # wavelength matrix (X columns)

# CLASS is coded Common / Club in the shared data.

# Core packages
suppressPackageStartupMessages({
  library(caret)
  library(pls)
  library(glmnet)
  library(resemble)
  library(dplyr)
  library(tidyr)
  library(doParallel)
  library(waves)      # SNV via pretreat_spectra()
})

## ---- Parallel backend ----
n_cores <- max(1L, parallel::detectCores() - 1L)   # leave one core free
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)
cat(sprintf("\nRegistered parallel backend with %d workers.\n", n_cores))

set.seed(123)  # overall seed

# Local PLS controls
my_control <- mbl_control(validation_type = c("NNv"))
my_waplsr  <- local_fit_wapls(min_pls_c = 2, max_pls_c = 30)

# ----------------------------
# Parameters
# ----------------------------
response_col <- "CAVOL"                  # <-- set to "FYELD" or "CAVOL"
eval_classes <- c("Common", "Club")      # two market classes to evaluate
s_values     <- seq(100, 700, by = 100)  # per-class training sizes
test_size    <- 100                      # test size from the evaluated class
n_reps       <- 25                       # number of replicates
quiet        <- FALSE                    # set TRUE to reduce console chatter

# ----------------------------
# Keep only the two classes and align rows
# ----------------------------
keep_idx <- which(pheno$CLASS %in% eval_classes)
pheno2   <- pheno[keep_idx, , drop = FALSE]
scans2   <- as.data.frame(scans[keep_idx, , drop = FALSE])  # caret-friendly

# ----------------------------
# (A) PREPROCESSING UTILITIES
# ----------------------------
# Spectral preprocessing via waves. Codes: raw=1, snv=2, d1=5, d2=6,
# snvd1=3, snvd2=4 (derivatives by adjacent differencing). MSC not used here.
apply_preproc <- function(X_train, X_test, method = "raw") {
  code <- c(raw = 1, snv = 2, d1 = 5, d2 = 6, snvd1 = 3, snvd2 = 4)[tolower(method)]
  if (is.na(code)) stop("unsupported preprocessing: ", method)
  wv <- function(X) as.matrix(waves::pretreat_spectra(as.data.frame(X, check.names = FALSE), pretreatment = code))
  list(train = wv(X_train), test = wv(X_test))
}

# ----------------------------
# (B) Best preprocessing per scenario (selected from the preprocessing sweep)
# ----------------------------
best_preproc_map <- tibble::tribble(
  ~eval_class, ~model,      ~training_source, ~preproc,
  "Common",    "LASSO",     "Matching",       "snv",
  "Common",    "LASSO",     "Opposing",       "snv",
  "Common",    "LASSO",     "Both",           "snv",
  "Common",    "LOCAL PLS", "Matching",       "raw",
  "Common",    "LOCAL PLS", "Opposing",       "raw",
  "Common",    "LOCAL PLS", "Both",           "raw",
  "Common",    "PLS",       "Matching",       "raw",
  "Common",    "PLS",       "Opposing",       "snvd2",
  "Common",    "PLS",       "Both",           "snvd1",
  "Club",      "LASSO",     "Matching",       "snv",
  "Club",      "LASSO",     "Opposing",       "d2",
  "Club",      "LASSO",     "Both",           "snv",
  "Club",      "LOCAL PLS", "Matching",       "raw",
  "Club",      "LOCAL PLS", "Opposing",       "snvd1",
  "Club",      "LOCAL PLS", "Both",           "raw",
  "Club",      "PLS",       "Matching",       "snv",
  "Club",      "PLS",       "Opposing",       "d2",
  "Club",      "PLS",       "Both",           "snv"
)

fallback_preproc <- "raw"  # if a combo isn't in the map, use this

choose_preproc <- function(cls, mdl, src) {
  hit <- best_preproc_map %>% dplyr::filter(eval_class == cls, model == mdl, training_source == src)
  if (nrow(hit) == 1) return(hit$preproc[1])
  fallback_preproc
}

# ----------------------------
# (C) caret trainControl (parallel-aware)
# ----------------------------
trctrl_pls  <- trainControl(method = "cv", number = 10, allowParallel = TRUE)
trctrl_enet <- trainControl(method = "cv", number = 10, allowParallel = TRUE)

# ----------------------------
# Results container
# ----------------------------
results_df <- data.frame(
  eval_class = character(),
  replicate = integer(),
  training_source = character(),   # Matching / Opposing / Both
  s_per_class = integer(),
  requested_total_train = integer(),
  # Which preprocessing used per model (traceability)
  preproc_pls   = character(),
  preproc_enet  = character(),
  preproc_local = character(),
  # PLS (caret)
  best_ncomp_pls = integer(),
  cor_pls        = numeric(),
  # LASSO / ENet (caret+glmnet)
  best_alpha_enet  = numeric(),
  best_lambda_enet = numeric(),
  n_nonzero_enet   = integer(),
  cor_enet         = numeric(),
  # Local PLS (mbl)
  cor_local      = numeric(),
  stringsAsFactors = FALSE
)

# ----------------------------
# Helper: safe correlation
# ----------------------------
safe_cor <- function(y, yhat) {
  y <- as.numeric(y); yhat <- as.numeric(yhat)
  keep <- is.finite(y) & is.finite(yhat)
  if (sum(keep) < 2L) return(NA_real_)
  suppressWarnings(cor(y[keep], yhat[keep]))
}

# ----------------------------
# (D) PROGRESS UTILITIES
# ----------------------------
total_sources <- length(eval_classes) * n_reps * length(s_values) * 3L  # approx (Matching/Opposing/Both)
progress_count <- 0L
start_time <- Sys.time()
pb <- utils::txtProgressBar(min = 0, max = total_sources, style = 3)

note <- function(msg, ...) {
  if (!quiet) cat(sprintf(msg, ...))
}
tick <- function(incr = 1L, cls, rep_i, s, src, ppls, pen, ploc) {
  # advance
  progress_count <<- progress_count + incr
  utils::setTxtProgressBar(pb, min(progress_count, total_sources))
  # ETA
  elapsed <- Sys.time() - start_time
  eta <- if (progress_count > 0)
    (as.numeric(elapsed, units = "secs") * (total_sources - progress_count) / progress_count) else NA
  eta_str <- if (is.na(eta)) "NA" else paste0(round(eta/60, 1), " min")
  pct <- 100 * progress_count / total_sources
  note("  ✓ %-6s | rep %02d | s=%3d | %-8s | PLS=%-6s ENET=%-6s LOCAL=%-6s | %5d/%5d (%.1f%%) | elapsed %s | ETA %s\n",
       cls, rep_i, s, src, ppls, pen, ploc, progress_count, total_sources, pct,
       format(elapsed, digits = 2), eta_str)
}

# ----------------------------
# Main loop (nested sampling) with scenario-specific preprocessing + PROGRESS
# ----------------------------
tryCatch({
  
  for (c_idx in seq_along(eval_classes)) {
    eval_class <- eval_classes[c_idx]
    opp_class  <- setdiff(eval_classes, eval_class)
    
    idx_eval <- which(pheno2$CLASS == eval_class)
    idx_opp  <- which(pheno2$CLASS == opp_class)
    
    if (length(idx_eval) < test_size) {
      stop("Not enough rows in '", eval_class, "' to draw ", test_size, " test samples.")
    }
    
    note("\n=== CLASS %-6s (%d/%d) | n_eval=%d | n_opp=%d ===\n",
         eval_class, c_idx, length(eval_classes), length(idx_eval), length(idx_opp))
    
    for (rep_i in seq_len(n_reps)) {
      set.seed(1000 + rep_i)  # reproducible test split
      
      test_idx  <- sample(idx_eval, size = test_size, replace = FALSE)
      X_testraw <- as.matrix(scans2[test_idx, , drop = FALSE])
      y_test    <- pheno2[[response_col]][test_idx]
      
      pool_eval <- setdiff(idx_eval, test_idx)
      pool_opp  <- idx_opp
      
      set.seed(2000 + rep_i)  # nested orderings
      nested_eval_order <- sample(pool_eval, length(pool_eval), replace = FALSE)
      nested_opp_order  <- sample(pool_opp,  length(pool_opp),  replace = FALSE)
      
      note("-- Start rep %02d | test=%d | pool_eval=%d | pool_opp=%d\n",
           rep_i, length(test_idx), length(pool_eval), length(pool_opp))
      
      for (s in s_values) {
        
        fit_and_score <- function(tr_idx, training_source_label, requested_total_train_val) {
          
          X_train_raw <- as.matrix(scans2[tr_idx, , drop = FALSE])
          y_train     <- pheno2[[response_col]][tr_idx]
          
          ## ---------------- PLS ----------------
          preproc_pls <- choose_preproc(eval_class, "PLS", training_source_label)
          pp_pls      <- apply_preproc(X_train_raw, X_testraw, preproc_pls)
          set.seed(11 + rep_i)
          pls_fit <- train(
            x = as.data.frame(pp_pls$train), y = y_train,
            method = "pls",
            preProcess = c("center", "scale"),
            tuneLength = 30,
            trControl  = trctrl_pls,
            metric     = "RMSE"
          )
          preds_pls <- as.numeric(predict(pls_fit, newdata = as.data.frame(pp_pls$test)))
          r_pls     <- safe_cor(y_test, preds_pls)
          
          ## --------------- LASSO / ENet ----------------
          preproc_enet <- choose_preproc(eval_class, "LASSO", training_source_label)
          pp_enet      <- apply_preproc(X_train_raw, X_testraw, preproc_enet)
          set.seed(22 + rep_i)
          enet_grid <- expand.grid(
            alpha  = seq(0, 1, by = 0.1),
            lambda = 10^seq(-6, 1, length.out = 60)
          )
          enet_fit <- train(
            x = as.data.frame(pp_enet$train), y = y_train,
            method = "glmnet",
            preProcess = c("center", "scale"),
            tuneGrid   = enet_grid,
            trControl  = trctrl_enet,
            metric     = "RMSE"
          )
          preds_enet <- as.numeric(predict(enet_fit, newdata = as.data.frame(pp_enet$test)))
          r_enet     <- safe_cor(y_test, preds_enet)
          best_alpha  <- enet_fit$bestTune$alpha
          best_lambda <- enet_fit$bestTune$lambda
          nnz <- NA_integer_
          coefs <- tryCatch(as.matrix(coef(enet_fit$finalModel, s = best_lambda)), error = function(e) NULL)
          if (!is.null(coefs)) nnz <- max(0L, sum(abs(coefs) > 0) - 1L)
          
          ## ---------------- LOCAL PLS ----------------
          preproc_local <- choose_preproc(eval_class, "LOCAL PLS", training_source_label)
          pp_local      <- apply_preproc(X_train_raw, X_testraw, preproc_local)
          
          k_use <- seq(50, 500, by = 50)                       # tune k over 50-500
          k_use <- k_use[k_use <= (nrow(pp_local$train) - 1L)] # capped at n_train-1
          if (length(k_use) == 0L) k_use <- max(2L, nrow(pp_local$train) - 1L)
          local_fit <- mbl(
            Xr = pp_local$train,  Yr = y_train,
            Xu = pp_local$test,
            k  = k_use,
            method      = my_waplsr,
            diss_method = "pls",
            diss_usage  = "predictors",
            control     = my_control,   # NNv for validation-driven k selection
            scale       = TRUE
          )
          
          # ---- Pick best k by smallest NNv RMSE and extract predictions for that k (vignette pattern) ----
          vr <- local_fit$validation_results$nearest_neighbor_validation
          if (is.null(vr) || !all(c("k", "rmse") %in% names(vr))) {
            stop("LOCAL PLS: nearest_neighbor_validation results are missing (no 'k'/'rmse').")
          }
          bki <- which.min(vr$rmse)
          bk  <- vr$k[bki]
          
          # All predictions (n_test x n_k); columns map to neighbor sizes
          pred_mat <- as.matrix(get_predictions(local_fit))
          
          # Column alignment: prefer exact name, else align by position using vr$k
          cn <- colnames(pred_mat)
          col_idx <- NULL
          if (!is.null(cn)) {
            if (as.character(bk) %in% cn) {
              col_idx <- which(cn == as.character(bk))
            } else if (paste0("k=", bk) %in% cn) {
              col_idx <- which(cn == paste0("k=", bk))
            }
          }
          if (is.null(col_idx)) {
            col_idx <- which(vr$k == bk)
            if (length(col_idx) != 1L) {
              stop("LOCAL PLS: could not match best k to a predictions column.")
            }
          }
          
          preds_local <- as.numeric(pred_mat[, col_idx])
          if (length(preds_local) != length(y_test)) {
            stop(sprintf("LOCAL PLS: prediction length mismatch (%d != %d).",
                         length(preds_local), length(y_test)))
          }
          
          r_local <- safe_cor(y_test, preds_local)
          
          # progress tick for this single (class,rep,s,source)
          tick(1L, eval_class, rep_i, s, training_source_label, preproc_pls, preproc_enet, preproc_local)
          
          # Row result
          data.frame(
            eval_class = eval_class,
            replicate  = rep_i,
            training_source = training_source_label,
            s_per_class = s,
            requested_total_train = requested_total_train_val,
            preproc_pls   = preproc_pls,
            preproc_enet  = preproc_enet,
            preproc_local = preproc_local,
            # PLS
            best_ncomp_pls = pls_fit$bestTune$ncomp,
            cor_pls        = r_pls,
            # LASSO / ENet
            best_alpha_enet  = best_alpha,
            best_lambda_enet = best_lambda,
            n_nonzero_enet   = nnz,
            cor_enet         = r_enet,
            # Local
            cor_local        = r_local,
            stringsAsFactors = FALSE
          )
        } # fit_and_score
        
        rows_list <- list()
        
        if (length(nested_eval_order) >= s) {
          rows_list[["Matching"]] <- fit_and_score(nested_eval_order[1:s], "Matching", s)
        }
        if (length(nested_opp_order) >= s) {
          rows_list[["Opposing"]] <- fit_and_score(nested_opp_order[1:s], "Opposing", s)
        }
        if (length(nested_eval_order) >= s && length(nested_opp_order) >= s) {
          rows_list[["Both"]] <- fit_and_score(c(nested_eval_order[1:s], nested_opp_order[1:s]),
                                               "Both", 2L * s)
        }
        
        if (length(rows_list)) {
          results_df <- bind_rows(results_df, dplyr::bind_rows(rows_list))
        }
      } # s
    }   # rep_i
  }     # eval_class
  
}, finally = {
  try(parallel::stopCluster(cl), silent = TRUE)
  doParallel::registerDoSEQ()
  utils::setTxtProgressBar(pb, total_sources); close(pb)
  cat("\nParallel backend stopped. Re-registered to sequential.\n")
})

# ----------------------------
# Quick look + summary
# ----------------------------
print(head(results_df, 12))

summary_df <- results_df %>%
  group_by(eval_class, training_source, s_per_class) %>%
  summarise(
    mean_cor_pls   = mean(cor_pls,  na.rm = TRUE),
    sd_cor_pls     = sd(cor_pls,    na.rm = TRUE),
    mean_cor_enet  = mean(cor_enet, na.rm = TRUE),
    sd_cor_enet    = sd(cor_enet,   na.rm = TRUE),
    mean_cor_local = mean(cor_local,na.rm = TRUE),
    sd_cor_local   = sd(cor_local,  na.rm = TRUE),
    n              = dplyr::n(),
    .groups = "drop"
  )

print(summary_df)

# ----------------------------
# Write per-replicate results for the learning-curve figure (Figure 3)
# ----------------------------
write.csv(results_df, file.path(data_dir, "growing_sample_size_validation.csv"), row.names = FALSE)
cat("\nSaved: growing_sample_size_validation.csv\n")
