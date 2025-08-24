# ============================ Setup ============================================
library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(glmnet)     # LASSO (inner)
library(logistf)    # Firth logistic (outer fits)
library(pROC)       # ROC/AUC
library(caret)      # confusionMatrix

set.seed(123)

# ============================ Load & prep data =================================
PCL_CAPS <- read_excel("Desktop/caps/PCL_CAPS_ALL.xlsx")

# keep rows with PCL items present; drop any remaining NA rows
PCL_CAPS <- PCL_CAPS[!is.na(PCL_CAPS$pcl_q1), ]
PCL_CAPS$caps_ptsd1 <- as.numeric(as.character(as.factor(PCL_CAPS$caps_ptsd1)))
PCL_CAPS <- na.omit(PCL_CAPS)

# ensure PCL items are numeric
for (i in 1:20) {
  var <- paste0("pcl_q", i)
  PCL_CAPS[[var]] <- as.numeric(PCL_CAPS[[var]])
}

# vector of the 20 PCL items
item_names <- paste0("pcl_q", 1:20)

# Precompute baseline totals (predictors only; not leakage)
PCL_CAPS$pcl_all_total <- rowSums(PCL_CAPS[, item_names, drop = FALSE], na.rm = TRUE)

# ============================ Nested bootstrap plan ============================
B_outer   <- 1000
threshold <- 0.20   # operating point for thresholded metrics

n <- nrow(PCL_CAPS)

# AUC storage
auc_full_vec <- c()
auc_sum_vec  <- c()
auc_mul_vec  <- c()

# thresholded metric storage (per model)
metrics_full_list <- list()
metrics_sum_list  <- list()
metrics_mul_list  <- list()

# also store vectors for CIs
sens_full_boot <- c(); spec_full_boot <- c(); prec_full_boot <- c(); recall_full_boot <- c(); f1_full_boot <- c()
sens_sum_boot  <- c(); spec_sum_boot  <- c(); prec_sum_boot  <- c(); recall_sum_boot  <- c(); f1_sum_boot  <- c()
sens_mul_boot  <- c(); spec_mul_boot  <- c(); prec_mul_boot  <- c(); recall_mul_boot  <- c(); f1_mul_boot  <- c()

# pooled OOB for ROC plotting
all_oob_labels      <- c()
all_oob_preds_full  <- c()
all_oob_preds_sum   <- c()
all_oob_preds_multi <- c()

# --- Importance tracking across outer loops (for Reduced-Multi) ---
sel_counts       <- setNames(rep(0L, length(item_names)), item_names)   # selection frequency
coef_logOR_store <- lapply(item_names, function(.) numeric(0))          # Firth log-ORs (training fit)
perm_drop_store  <- lapply(item_names, function(.) numeric(0))          # OOB permutation AUC drops

outer_success <- 0L   # number of outer iterations successfully used for Reduced-Multi

# ============================ Outer bootstrap loop =============================
for (b in 1:B_outer) {
  print(b)
  # cat("Outer iter:", b, "\n")
  boot_idx <- sample(1:n, replace = TRUE)
  oob_idx  <- setdiff(1:n, unique(boot_idx))
  if (length(oob_idx) < 5) next  # too few OOB
  
  train_df <- PCL_CAPS[boot_idx, , drop = FALSE]
  oob_df   <- PCL_CAPS[oob_idx,  , drop = FALSE]
  
  # class sanity in OOB
  labels <- oob_df$caps_ptsd1
  if (length(unique(labels)) < 2) next
  lab_fac <- factor(labels, levels = c(0, 1))
  
  # ---------------- Inner selection (LASSO on TRAIN only) ----------------------
  set.seed(123 + b*17)  # reproducible inner CV
  X_tr <- as.matrix(train_df[, item_names, drop = FALSE])
  y_tr <- train_df$caps_ptsd1
  
  # if degenerate in TRAIN (single class), skip
  if (length(unique(y_tr)) < 2) next
  
  cvfit <- tryCatch(cv.glmnet(X_tr, y_tr, alpha = 1, family = "binomial",
                              nfolds = 5, standardize = TRUE),
                    error = function(e) NULL)
  if (is.null(cvfit)) next
  
  coef_mat <- as.matrix(coef(cvfit, s = "lambda.min"))
  nonzero  <- rownames(coef_mat)[coef_mat[,1] != 0]
  selected_vars <- setdiff(nonzero, "(Intercept)")
  
  # fallback if nothing selected (rare)
  if (length(selected_vars) == 0) {
    # select the single strongest (by absolute coef at 1-SE) or just the highest variance item
    coef_mat_1se <- as.matrix(coef(cvfit, s = "lambda.1se"))
    nz_1se <- setdiff(rownames(coef_mat_1se)[coef_mat_1se[,1] != 0], "(Intercept)")
    selected_vars <- if (length(nz_1se)) nz_1se else item_names[1]
  }
  
  # build a summed feature from SELECTED vars (train + OOB, predictors only; no outcome)
  train_df$pcl_selected_total <- rowSums(train_df[, selected_vars, drop = FALSE], na.rm = TRUE)
  oob_df$pcl_selected_total   <- rowSums(oob_df[,   selected_vars, drop = FALSE], na.rm = TRUE)
  
  # ---------------- Fit 3 models on TRAIN (Firth logistic) ---------------------
  # Full-SUM model: single predictor = sum of all 20 items
  fit_full <- tryCatch(
    logistf(caps_ptsd1 ~ pcl_all_total, data = train_df,
            control   = logistf.control(maxit = 1000, maxstep = 0.5),
            plcontrol = logistpl.control(maxit = 2000)),
    error = function(e) NULL
  )
  if (is.null(fit_full)) next
  
  # Reduced-SUM model: single predictor = sum of selected items
  fit_sum <- tryCatch(
    logistf(caps_ptsd1 ~ pcl_selected_total, data = train_df,
            control   = logistf.control(maxit = 1000, maxstep = 0.5),
            plcontrol = logistpl.control(maxit = 2000)),
    error = function(e) NULL
  )
  if (is.null(fit_sum)) next
  
  # Reduced-MULTI model: each selected item individually
  fm_mul <- as.formula(paste("caps_ptsd1 ~", paste(selected_vars, collapse = " + ")))
  fit_mul <- tryCatch(
    logistf(fm_mul, data = train_df,
            control   = logistf.control(maxit = 1000, maxstep = 0.5),
            plcontrol = logistpl.control(maxit = 3000)),
    error = function(e) NULL
  )
  if (is.null(fit_mul)) next
  
  # ---------------- OOB predictions & ROC/AUC ----------------------------------
  preds_full <- tryCatch(predict(fit_full, newdata = oob_df, type = "response"), error = function(e) NULL)
  preds_sum  <- tryCatch(predict(fit_sum,  newdata = oob_df, type = "response"), error = function(e) NULL)
  preds_mul  <- tryCatch(predict(fit_mul,  newdata = oob_df, type = "response"), error = function(e) NULL)
  if (any(vapply(list(preds_full, preds_sum, preds_mul), is.null, logical(1)))) next
  
  r_full <- tryCatch(pROC::roc(lab_fac, preds_full, levels = c(0,1), direction = "auto", quiet = TRUE), error = function(e) NULL)
  r_sum  <- tryCatch(pROC::roc(lab_fac, preds_sum,  levels = c(0,1), direction = "auto", quiet = TRUE), error = function(e) NULL)
  r_mul  <- tryCatch(pROC::roc(lab_fac, preds_mul,  levels = c(0,1), direction = "auto", quiet = TRUE), error = function(e) NULL)
  if (any(vapply(list(r_full, r_sum, r_mul), is.null, logical(1)))) next
  
  auc_full_vec <- c(auc_full_vec, as.numeric(pROC::auc(r_full)))
  auc_sum_vec  <- c(auc_sum_vec,  as.numeric(pROC::auc(r_sum)))
  auc_mul_vec  <- c(auc_mul_vec,  as.numeric(pROC::auc(r_mul)))
  
  # pooled OOB storage for later ROC curves
  all_oob_labels      <- c(all_oob_labels,      as.numeric(labels))
  all_oob_preds_full  <- c(all_oob_preds_full,  as.numeric(preds_full))
  all_oob_preds_sum   <- c(all_oob_preds_sum,   as.numeric(preds_sum))
  all_oob_preds_multi <- c(all_oob_preds_multi, as.numeric(preds_mul))
  
  # ---------------- Thresholded metrics at fixed operating point ---------------
  # Factor-encode predictions at threshold
  bin_full_fac <- factor(ifelse(preds_full >= threshold, 1, 0), levels = c(0, 1))
  bin_sum_fac  <- factor(ifelse(preds_sum  >= threshold, 1, 0), levels = c(0, 1))
  bin_mul_fac  <- factor(ifelse(preds_mul  >= threshold, 1, 0), levels = c(0, 1))
  
  cm_full <- caret::confusionMatrix(bin_full_fac, lab_fac, positive = "1")
  cm_sum  <- caret::confusionMatrix(bin_sum_fac,  lab_fac, positive = "1")
  cm_mul  <- caret::confusionMatrix(bin_mul_fac,  lab_fac, positive = "1")
  
  mnames <- c("Sensitivity","Specificity","Precision","Recall","F1")
  metrics_full_list[[length(metrics_full_list) + 1]] <- cm_full$byClass[mnames]
  metrics_sum_list [[length(metrics_sum_list)  + 1]] <- cm_sum$byClass[mnames]
  metrics_mul_list [[length(metrics_mul_list)  + 1]] <- cm_mul$byClass[mnames]
  
  sens_full_boot <- c(sens_full_boot, unname(cm_full$byClass["Sensitivity"]))
  spec_full_boot <- c(spec_full_boot, unname(cm_full$byClass["Specificity"]))
  prec_full_boot <- c(prec_full_boot, unname(cm_full$byClass["Precision"]))
  recall_full_boot <- c(recall_full_boot, unname(cm_full$byClass["Recall"]))
  f1_full_boot <- c(f1_full_boot, unname(cm_full$byClass["F1"]))
  
  sens_sum_boot <- c(sens_sum_boot, unname(cm_sum$byClass["Sensitivity"]))
  spec_sum_boot <- c(spec_sum_boot, unname(cm_sum$byClass["Specificity"]))
  prec_sum_boot <- c(prec_sum_boot, unname(cm_sum$byClass["Precision"]))
  recall_sum_boot <- c(recall_sum_boot, unname(cm_sum$byClass["Recall"]))
  f1_sum_boot <- c(f1_sum_boot, unname(cm_sum$byClass["F1"]))
  
  sens_mul_boot <- c(sens_mul_boot, unname(cm_mul$byClass["Sensitivity"]))
  spec_mul_boot <- c(spec_mul_boot, unname(cm_mul$byClass["Specificity"]))
  prec_mul_boot <- c(prec_mul_boot, unname(cm_mul$byClass["Precision"]))
  recall_mul_boot <- c(recall_mul_boot, unname(cm_mul$byClass["Recall"]))
  f1_mul_boot <- c(f1_mul_boot, unname(cm_mul$byClass["F1"]))
  
  # ---------------- Importance for Reduced-Multi -------------------------------
  # 1) Stability: which items were selected by inner LASSO?
  sel_counts[selected_vars] <- sel_counts[selected_vars] + 1L
  
  # 2) Effect sizes: training log-ORs from Firth (skip intercept)
  coefs_mul <- coef(fit_mul)
  vars_in_model <- intersect(names(coefs_mul), item_names)
  for (vn in vars_in_model) {
    coef_logOR_store[[vn]] <- c(coef_logOR_store[[vn]], unname(coefs_mul[vn]))
  }
  
  # 3) OOB permutation AUC drop (contribution) for each selected item
  base_auc_mul <- as.numeric(pROC::auc(r_mul))
  for (vn in selected_vars) {
    oob_perm <- oob_df
    oob_perm[[vn]] <- sample(oob_perm[[vn]])  # permute one variable only
    pred_perm <- tryCatch(predict(fit_mul, newdata = oob_perm, type = "response"), error = function(e) NULL)
    if (!is.null(pred_perm)) {
      r_perm <- tryCatch(pROC::roc(lab_fac, pred_perm, levels = c(0,1), direction = "auto", quiet = TRUE), error = function(e) NULL)
      if (!is.null(r_perm)) {
        perm_drop_store[[vn]] <- c(perm_drop_store[[vn]], base_auc_mul - as.numeric(pROC::auc(r_perm)))
      }
    }
  }
  
  outer_success <- outer_success + 1L
}

cat("\nOuter bootstrap iterations kept:", outer_success, "out of", B_outer, "\n")

# ============================ Aggregate metrics & CIs ==========================
avg_tbl <- function(lst) {
  if (!length(lst)) return(tibble(Sensitivity=NA, Specificity=NA, Precision=NA, Recall=NA, F1=NA))
  bind_rows(lst) %>% summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))
}

avg_full <- avg_tbl(metrics_full_list)
avg_sum  <- avg_tbl(metrics_sum_list)
avg_mul  <- avg_tbl(metrics_mul_list)

report_model <- function(name, avg_tbl, auc_vec) {
  cat(paste0("\n====== ", name, " (averaged across bootstraps) ======\n"))
  print(round(avg_tbl, 3))
  cat("AUC (mean):", round(mean(auc_vec, na.rm = TRUE), 3), "\n")
  cat("AUC 95% CI:", round(quantile(auc_vec, c(0.025, 0.975), na.rm = TRUE), 3), "\n")
}

report_model("Full-Sum",       avg_full, auc_full_vec)
report_model("Reduced-Sum",    avg_sum,  auc_sum_vec)
report_model("Reduced-Multi",  avg_mul,  auc_mul_vec)




safe_mean <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_quant <- function(x, probs = c(0.025, 0.975)) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(rep(NA_real_, length(probs)))
  as.numeric(quantile(x, probs = probs, na.rm = TRUE, names = FALSE))
}

mk_ci_df <- function(sens, spec, prec, rec, f1, auc, label) {
  mns <- c(
    safe_mean(sens),
    safe_mean(spec),
    safe_mean(prec),
    safe_mean(rec),
    safe_mean(f1),
    safe_mean(auc)
  )
  
  lo <- c(
    safe_quant(sens, 0.025),
    safe_quant(spec, 0.025),
    safe_quant(prec, 0.025),
    safe_quant(rec,  0.025),
    safe_quant(f1,   0.025),
    safe_quant(auc,  0.025)
  )
  
  hi <- c(
    safe_quant(sens, 0.975),
    safe_quant(spec, 0.975),
    safe_quant(prec, 0.975),
    safe_quant(rec,  0.975),
    safe_quant(f1,   0.975),
    safe_quant(auc,  0.975)
  )
  
  data.frame(
    Model  = label,
    Metric = c("Sensitivity","Specificity","Precision","Recall","F1","AUC"),
    Mean   = round(mns, 3),
    CI_Low = round(lo, 3),
    CI_High= round(hi, 3),
    row.names = NULL,
    check.names = FALSE
  )
}

ci_full <- mk_ci_df(sens_full_boot, spec_full_boot, prec_full_boot, recall_full_boot, f1_full_boot, auc_full_vec, "Full-Sum")
ci_sum  <- mk_ci_df(sens_sum_boot,  spec_sum_boot,  prec_sum_boot,  recall_sum_boot,  f1_sum_boot,  auc_sum_vec,  "Reduced-Sum")
ci_mul  <- mk_ci_df(sens_mul_boot,  spec_mul_boot,  prec_mul_boot,  recall_mul_boot,  f1_mul_boot,  auc_mul_vec,  "Reduced-Multi")
ci_all  <- bind_rows(ci_full, ci_sum, ci_mul)

cat("\n====== Metric Means and 95% CIs (all models) ======\n")
print(ci_all)

# ============================ ROC Curves (pooled OOB) ==========================
if (length(all_oob_labels) > 1 &&
    length(all_oob_labels) == length(all_oob_preds_full) &&
    length(all_oob_labels) == length(all_oob_preds_sum)  &&
    length(all_oob_labels) == length(all_oob_preds_multi)) {
  
  lab_all <- factor(all_oob_labels, levels = c(0,1))
  roc_full  <- pROC::roc(lab_all, all_oob_preds_full,  levels = c(0,1), direction = "auto", quiet = TRUE)
  roc_sum   <- pROC::roc(lab_all, all_oob_preds_sum,   levels = c(0,1), direction = "auto", quiet = TRUE)
  roc_multi <- pROC::roc(lab_all, all_oob_preds_multi, levels = c(0,1), direction = "auto", quiet = TRUE)
  
  auc_full <- as.numeric(pROC::auc(roc_full))
  auc_sum  <- as.numeric(pROC::auc(roc_sum))
  auc_mul  <- as.numeric(pROC::auc(roc_multi))
  
  roc_df <- dplyr::bind_rows(
    data.frame(FPR = 1 - roc_full$specificities,  TPR = roc_full$sensitivities,  Model = "Full-Sum"),
    data.frame(FPR = 1 - roc_sum$specificities,   TPR = roc_sum$sensitivities,   Model = "Reduced-Sum"),
    data.frame(FPR = 1 - roc_multi$specificities, TPR = roc_multi$sensitivities, Model = "Reduced-Multi")
  )
  
  ggplot(roc_df, aes(x = FPR, y = TPR, color = Model)) +
    geom_line(size = 1.3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
      title = "Pooled OOB ROC Curves (3 Models)",
      x = "1 - Specificity (False Positive Rate)",
      y = "Sensitivity (True Positive Rate)",
      subtitle = paste0(
        "AUC Full-Sum = ", sprintf("%.3f", auc_full),
        " | AUC Reduced-Sum = ", sprintf("%.3f", auc_sum),
        " | AUC Reduced-Multi = ", sprintf("%.3f", auc_mul)
      )
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position      = c(0.98, 0.02),
      legend.justification = c(1, 0),
      legend.background    = element_rect(fill = scales::alpha("white", 0.7), color = NA),
      legend.title         = element_blank()
    )
  
  cat("\nPairwise DeLong tests (pooled OOB):\n")
  print(pROC::roc.test(roc_sum,   roc_full,  method = "delong", paired = TRUE))
  print(pROC::roc.test(roc_multi, roc_full,  method = "delong", paired = TRUE))
  print(pROC::roc.test(roc_multi, roc_sum,   method = "delong", paired = TRUE))
}

# ============================ Item Importance Summary ==========================
# Pre-reqs expected from earlier code:
# - item_names: character vector of item names to summarize (e.g., selected_vars or pcl_items)
# - outer_success: number of outer resamples that finished successfully
# - sel_counts: named integer vector of selection counts per item (names = item_names)
# - coef_logOR_store: named list; each [[item]] is a numeric vector of log-ORs across outer runs
# - perm_drop_store: named list; each [[item]] is a numeric vector of AUC drops across outer runs
# If any of the above may be missing, this code guards gracefully.

suppressWarnings({
  if (!exists("item_names")) {
    item_names <- if (exists("selected_vars")) selected_vars else colnames(PCL_CAPS)[grepl("^pcl_q", colnames(PCL_CAPS))]
  }
  item_names <- intersect(item_names, colnames(PCL_CAPS))  # keep only valid columns
})

# Denominator for selection frequency = number of outer successes
kept_B <- if (exists("outer_success")) as.integer(outer_success) else NA_integer_

# Ensure sel_counts exists and is complete over item_names
if (!exists("sel_counts") || is.null(sel_counts)) sel_counts <- setNames(integer(0), character(0))
sel_counts_full <- setNames(integer(length(item_names)), item_names)
if (length(sel_counts)) {
  idx <- intersect(names(sel_counts), item_names)
  sel_counts_full[idx] <- as.integer(sel_counts[idx])
}

# Safe extractors for per-item logs/permutation drops
get_logs  <- function(vn) {
  if (exists("coef_logOR_store") && !is.null(coef_logOR_store[[vn]])) {
    as.numeric(coef_logOR_store[[vn]])
  } else numeric(0)
}
get_perms <- function(vn) {
  if (exists("perm_drop_store") && !is.null(perm_drop_store[[vn]])) {
    as.numeric(perm_drop_store[[vn]])
  } else numeric(0)
}

summarize_item <- function(vn) {
  logs  <- get_logs(vn)
  perms <- get_perms(vn)
  c(
    sel_freq      = if (!is.na(kept_B) && kept_B > 0) sel_counts_full[vn] / kept_B else NA_real_,
    logOR_median  = if (length(logs)) stats::median(logs, na.rm = TRUE) else NA_real_,
    logOR_q1      = if (length(logs)) stats::quantile(logs, 0.25, na.rm = TRUE) else NA_real_,
    logOR_q3      = if (length(logs)) stats::quantile(logs, 0.75, na.rm = TRUE) else NA_real_,
    OR_median     = if (length(logs)) exp(stats::median(logs, na.rm = TRUE)) else NA_real_,
    perm_drop_med = if (length(perms)) stats::median(perms, na.rm = TRUE) else NA_real_,
    perm_drop_q1  = if (length(perms)) stats::quantile(perms, 0.25, na.rm = TRUE) else NA_real_,
    perm_drop_q3  = if (length(perms)) stats::quantile(perms, 0.75, na.rm = TRUE) else NA_real_
  )
}

# Build importance table
imp_mat <- t(vapply(item_names, summarize_item,
                    FUN.VALUE = c(sel_freq = 0, logOR_median = 0, logOR_q1 = 0, logOR_q3 = 0,
                                  OR_median = 0, perm_drop_med = 0, perm_drop_q1 = 0, perm_drop_q3 = 0)))
imp_df  <- as.data.frame(imp_mat, stringsAsFactors = FALSE)
imp_df$Item <- rownames(imp_mat)

# Ensure numeric columns are numeric; keep Item as character
num_cols <- setdiff(colnames(imp_df), "Item")
imp_df[num_cols] <- lapply(imp_df[num_cols], function(x) as.numeric(x))

# Order by stability (sel_freq) then contribution (perm_drop_med)
imp_df <- imp_df |>
  dplyr::arrange(dplyr::desc(sel_freq), dplyr::desc(perm_drop_med)) |>
  dplyr::relocate(Item)

# ------- Print Top N table (round numeric columns only) -------
topN <- min(10L, nrow(imp_df))
cat("\n====== Top Items by Stability & Contribution ======\n")
if (topN > 0) {
  top_tbl <- imp_df |>
    dplyr::slice(1:topN) |>
    dplyr::select(Item, sel_freq, OR_median, perm_drop_med) |>
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3)))
  print(top_tbl, row.names = FALSE)
} else {
  cat("(No items to display)\n")
}

# ------- Scatter: selection stability vs. contribution -------
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  warning("Package 'ggrepel' not installed; plotting without repelled labels.")
  ggplot(imp_df, aes(x = sel_freq, y = perm_drop_med, label = Item)) +
    geom_point() +
    geom_text(size = 3, vjust = -0.6) +
    labs(x = "Selection frequency (stability)",
         y = "Median AUC drop on OOB (permutation importance)",
         title = "Item Importance: stability vs. contribution (Reduced-Multi)") +
    theme_minimal(base_size = 13)
} else {
  ggplot(imp_df, aes(x = sel_freq, y = perm_drop_med, label = Item)) +
    geom_point() +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 20) +
    labs(x = "Selection frequency (stability)",
         y = "Median AUC drop on OOB (permutation importance)",
         title = "Item Importance: stability vs. contribution (Reduced-Multi)") +
    theme_minimal(base_size = 13)
}


# ============================ Distribution of Responses per PCL-5 Question ==========================
pcl_long <- PCL_CAPS %>%
  select(all_of(item_names)) %>%
  pivot_longer(cols = everything(), names_to = "Item", values_to = "Response") %>%
  filter(!is.na(Response))

pcl_long$Item <- factor(pcl_long$Item, levels = item_names)

value_counts <- pcl_long %>%
  group_by(Item, Response) %>%
  summarise(Count = n(), .groups = "drop")

ggplot(value_counts, aes(x = Item, y = Count, fill = factor(Response))) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count),
            position = position_stack(vjust = 0.5),
            size = 3, color = "black") +
  scale_fill_brewer(palette = "YlGnBu", name = "Response Value") +
  labs(title = "Distribution of Responses per PCL-5 Question",
       x = "PCL-5 Item", y = "Number of Responses") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#=========== model with 6 items ================================================
# --- Libraries ----------------------------------------------------------------
library(readxl)
library(dplyr)
library(logistf)
library(pROC)
library(caret)

set.seed(1234)

# --- Load & prep ---------------------------------------------------------------
PCL_CAPS <- read_excel("Desktop/caps/PCL_CAPS_ALL.xlsx") %>%
  filter(!is.na(pcl_q1))

# outcome as 0/1 numeric
PCL_CAPS$caps_ptsd1 <- as.numeric(as.character(as.factor(PCL_CAPS$caps_ptsd1)))

# ensure all PCL items numeric
for (i in 1:20) {
  vi <- paste0("pcl_q", i)
  PCL_CAPS[[vi]] <- as.numeric(PCL_CAPS[[vi]])
}
PCL_CAPS <- na.omit(PCL_CAPS)

# --- Your chosen items ---------------------------------------------------------
selected_vars <- c("pcl_q10","pcl_q2","pcl_q7","pcl_q15","pcl_q4","pcl_q9")

# --- Train/test split (stratified; enforce both classes in each) ---------------
make_split <- function(dat, p = 0.7, tries = 100) {
  for (t in 1:tries) {
    idx <- caret::createDataPartition(dat$caps_ptsd1, p = p, list = FALSE)
    tr <- dat[idx, , drop = FALSE]
    te <- dat[-idx, , drop = FALSE]
    if (length(unique(tr$caps_ptsd1)) == 2 && length(unique(te$caps_ptsd1)) == 2) {
      return(list(train = tr, test = te))
    }
  }
  stop("Could not make a split with both classes in train and test.")
}
spl <- make_split(PCL_CAPS, p = 0.7)
train_df <- spl$train
test_df  <- spl$test

# --- Fit Firth logistic on train ----------------------------------------------
fm <- as.formula(paste("caps_ptsd1 ~", paste(selected_vars, collapse = " + ")))
fit <- logistf(fm, data = train_df,
               control   = logistf.control(maxit = 1000, maxstep = 0.5),
               plcontrol = logistpl.control(maxit = 1000))
cat("\nCoefficients (Firth on train):\n")
print(round(coef(fit), 3))

# --- Scores & train thresholds -------------------------------------------------
train_scores <- predict(fit, newdata = train_df, type = "response")
test_scores  <- predict(fit, newdata  = test_df,  type = "response")

roc_train <- pROC::roc(factor(train_df$caps_ptsd1, levels = c(0,1)), train_scores,
                       levels = c(0,1), direction = "auto", quiet = TRUE)
roc_test  <- pROC::roc(factor(test_df$caps_ptsd1,  levels = c(0,1)), test_scores,
                       levels = c(0,1), direction = "auto", quiet = TRUE)

# AUCs with 95% CI (DeLong)
auc_train <- as.numeric(pROC::auc(roc_train))
ci_train  <- as.numeric(pROC::ci.auc(roc_train, method = "delong"))
auc_test  <- as.numeric(pROC::auc(roc_test))
ci_test   <- as.numeric(pROC::ci.auc(roc_test, method = "delong"))

cat(sprintf("\nTrain AUC = %.3f (95%% CI %.3f–%.3f)\n", auc_train, ci_train[1], ci_train[3]))
cat(sprintf(" Test  AUC = %.3f (95%% CI %.3f–%.3f)\n", auc_test,  ci_test[1],  ci_test[3]))

# Thresholds: fixed 0.20, Youden on train, and high-sensitivity (~0.90) on train
thr_fixed <- 0.20
thr_youden <- as.numeric(pROC::coords(roc_train, x = "best", best.method = "youden",
                                      ret = "threshold", transpose = FALSE))
thr_sens90 <- tryCatch({
  as.numeric(pROC::coords(roc_train, x = 0.90, input = "sensitivity",
                          ret = "threshold", transpose = FALSE))
}, error = function(e) NA_real_)

cat(sprintf("\nChosen thresholds (from TRAIN): fixed=%.2f | Youden=%.3f | Sens90=%s\n",
            thr_fixed, thr_youden, ifelse(is.na(thr_sens90), "NA", sprintf("%.3f", thr_sens90))))

# --- Helper to compute classification metrics ---------------------------------
get_metrics <- function(scores, labels, thr) {
  pred <- factor(ifelse(scores >= thr, 1, 0), levels = c(0,1))
  lab  <- factor(labels, levels = c(0,1))
  cm   <- caret::confusionMatrix(pred, lab, positive = "1")
  c(Accuracy    = unname(cm$overall["Accuracy"]),
    Sensitivity = unname(cm$byClass["Sensitivity"]),
    Specificity = unname(cm$byClass["Specificity"]),
    Precision   = unname(cm$byClass["Precision"]),
    Recall      = unname(cm$byClass["Recall"]),
    F1          = unname(cm$byClass["F1"]))
}

# --- Evaluate on TEST at each threshold ---------------------------------------
m_fixed  <- get_metrics(test_scores, test_df$caps_ptsd1, thr_fixed)
m_youden <- get_metrics(test_scores, test_df$caps_ptsd1, thr_youden)
m_sens90 <- if (!is.na(thr_sens90)) get_metrics(test_scores, test_df$caps_ptsd1, thr_sens90) else rep(NA_real_, 6)

# --- Report --------------------------------------------------------------------
print_block <- function(name, m) {
  cat(paste0("\n== Test metrics @ ", name, " ==\n"))
  print(round(m, 3))
}
print_block(sprintf("fixed %.2f", thr_fixed), m_fixed)
print_block(sprintf("Youden %.3f (train)", thr_youden), m_youden)
if (!is.na(thr_sens90[1])) print_block(sprintf("Sens≈0.90 (train) %.3f", thr_sens90), m_sens90)
