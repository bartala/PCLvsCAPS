library(readxl)
library(Metrics)
library(logistf)
library(pROC)
library(glmnet)
library(ggplot2)
library(dplyr)

# load the data
PCL_CAPS <- read_excel("path_to_data/CAPS_deidentified_t3retestadded_9.21_SM.xlsx")

keep <- c("pcl_q1",
          "pcl_q2",	
          "pcl_q3",
          "pcl_q4",
          "pcl_q5",	
          "pcl_q6",
          "pcl_q7",
          "pcl_q8",
          "pcl_q9",
          "pcl_q10",
          "pcl_q11",
          "pcl_q12",
          "pcl_q13",
          "pcl_q14",
          "pcl_q15",
          "pcl_q16",
          "pcl_q17",
          "pcl_q18",
          "pcl_q19",
          "pcl_q20",
          "caps_ptsd1"
          )

PCL_CAPS <- PCL_CAPS[,keep]
PCL_CAPS <- PCL_CAPS[!is.na(PCL_CAPS$pcl_q1),]

PCL_CAPS$caps_ptsd1 <- as.factor(PCL_CAPS$caps_ptsd1)
PCL_CAPS$caps_ptsd1 <- as.numeric(as.character(PCL_CAPS$caps_ptsd1))

# Make sure all predictors are numeric
for (i in 1:20) {
  var <- paste0("pcl_q", i)
  PCL_CAPS[[var]] <- as.numeric(PCL_CAPS[[var]])
}



##------------------------------------------------------------------------------
## STEP 1: Use LASSO logistic regression to select the most important items (PCL5 questions)
##------------------------------------------------------------------------------

# Prepare data: assume PCL_CAPS is your data.frame
# Outcome is be binary (0/1)
y <- PCL_CAPS$caps_ptsd1
X <- as.matrix(PCL_CAPS[, grep("^pcl_q", names(PCL_CAPS))])

set.seed(123)

# Find the optimal lambda (penalty) using cross-validation.
cv_lasso <- cv.glmnet(X, y, alpha = 1, family = "binomial", nfolds = 5)

nonzero_coefs <- coef(cv_lasso, s = "lambda.min")

# Extract selected variables (non-zero coefficients)
selected_vars <- rownames(nonzero_coefs)[which(nonzero_coefs != 0)]
selected_vars <- setdiff(selected_vars, "(Intercept)")

cat("Selected PCL items by LASSO:\n")
print(selected_vars)

#  "pcl_q2"  "pcl_q4"  "pcl_q7"  "pcl_q9"  "pcl_q10" "pcl_q13" "pcl_q14" "pcl_q15"

#-------------------------------------------------------------------------------
## STEP 2: Bootstrap for Reduced and Full models
#-------------------------------------------------------------------------------
library(logistf)
library(pROC)

# Storage for bootstrap AUCs and pooled OOB predictions
auc_boot_full <- c()
auc_boot_reduced <- c()
all_oob_preds_full <- c()
all_oob_preds_reduced <- c()
all_oob_labels <- c()

set.seed(123)
n_boot <- 1000
n <- nrow(PCL_CAPS)

for (i in 1:n_boot) {
  boot_idx <- sample(1:n, replace = TRUE)
  oob_idx <- setdiff(1:n, unique(boot_idx))
  if (length(oob_idx) < 5) next
  
  try({
    # Fit full model
    model_full <- logistf(caps_ptsd1 ~ ., 
                          data = PCL_CAPS[boot_idx, ],
                          control = logistf.control(maxit = 1000, maxstep = 0.5),
                          plcontrol = logistpl.control(maxit = 1000))
    preds_full <- predict(model_full, newdata = PCL_CAPS[oob_idx, ], type = "response")
    
    # Fit reduced model
    model_reduced <- logistf(
      formula = as.formula(paste("caps_ptsd1 ~", paste(selected_vars, collapse = " + "))),
      data = PCL_CAPS[boot_idx, ],
      control = logistf.control(maxit = 1000, maxstep = 0.5),
      plcontrol = logistpl.control(maxit = 1000)
    )
    preds_reduced <- predict(model_reduced, newdata = PCL_CAPS[oob_idx, ], type = "response")
    
    labels <- PCL_CAPS$caps_ptsd1[oob_idx]
    
    # Store AUCs
    auc_boot_full <- c(auc_boot_full, auc(labels, preds_full))
    auc_boot_reduced <- c(auc_boot_reduced, auc(labels, preds_reduced))
    
    # Accumulate predictions for pooled ROC
    all_oob_preds_full <- c(all_oob_preds_full, preds_full)
    all_oob_preds_reduced <- c(all_oob_preds_reduced, preds_reduced)
    all_oob_labels <- c(all_oob_labels, labels)
    
  }, silent = TRUE)
}

# ----- Boxplot -----
df_auc <- data.frame(
  AUC = c(auc_boot_full, auc_boot_reduced),
  Model = factor(c(rep("Full Model", length(auc_boot_full)),
                   rep("Reduced Model", length(auc_boot_reduced))))
)

boxplot(AUC ~ Model, data = df_auc,
        col = c("lightblue", "lightgreen"),
        main = "Bootstrap AUC Comparison",
        ylab = "AUC",
        notch = TRUE,
        border = "gray30")
grid()


# Wilcoxon signed-rank test
wilcox.test(auc_boot_reduced, auc_boot_full, paired = TRUE)


# --- Compute ROC curves ---
roc_full <- roc(all_oob_labels, all_oob_preds_full, direction = "<", quiet = TRUE)
roc_reduced <- roc(all_oob_labels, all_oob_preds_reduced, direction = "<", quiet = TRUE)

# Extract ROC data for ggplot
roc_full_df <- data.frame(
  FPR = 1 - roc_full$specificities,
  TPR = roc_full$sensitivities,
  Model = "Full Model"
)

roc_reduced_df <- data.frame(
  FPR = 1 - roc_reduced$specificities,
  TPR = roc_reduced$sensitivities,
  Model = "Reduced Model"
)

# Combine data
roc_df <- bind_rows(roc_full_df, roc_reduced_df)

# Create plot
ggplot(roc_df, aes(x = FPR, y = TPR, color = Model)) +
  geom_line(size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("Full Model" = "blue", "Reduced Model" = "darkgreen")) +
  labs(
    title = "ROC Curves from Pooled OOB Predictions",
    x = "1 - Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)",
    color = "Model"
  ) +
  annotate("text", x = 0.65, y = 0.2,
           label = paste0("AUC Full: ", round(auc(roc_full), 3), 
                          "\nAUC Reduced: ", round(auc(roc_reduced), 3)),
           hjust = 0, size = 4) +
  theme_minimal(base_size = 14) +
  coord_equal()

