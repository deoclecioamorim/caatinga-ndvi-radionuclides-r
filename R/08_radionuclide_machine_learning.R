#Manuscript: Section 2.3.3 (Random Forest and GBM); Results 3.3; Fig. 8 and Fig. A.5
#Classification of the areas (A x B) from natural radionuclides and variable importance
#Input: data/processed/soil_radionuclides.xlsx
#
#Design (avoids information leakage between subsamples of the same plot):
#- training/test split at the PLOT level (70/30, stratified by area): the 3 subsamples of a plot
#  are all in the training set or all in the test set
#- cross-validation grouped by plot (5 folds x 20 repeats) on the training set, for hyperparameter
#  tuning and model comparison
#- comparison of RF and GBM by paired Wilcoxon test (both models are evaluated on the same folds)
#- test-set performance with the 95% confidence interval of accuracy
#- permutation importance (Random Forest, ranger) and SHAP values

#Control options --------------------------------------------------------------------------------
options(OutDec = ".")
rm(list = ls(all = TRUE))


# Packages -------------------------------------------------------------------------------------
for (p in c("readxl", "writexl", "dplyr", "tibble", "caret", "ranger", "gbm", "fastshap", "shapviz", "ggplot2", "patchwork")) {
  if (!require(p, character.only = TRUE)) { install.packages(p, dependencies = TRUE); library(p, character.only = TRUE) }
}
source("R/00_publication_theme.R") #figure standard (font and print size)

predictors <- c("ra228", "th228", "ra226", "k40")
radionuclide_labels <- c(ra228 = expression({}^{228}*Ra), th228 = expression({}^{228}*Th),
                         ra226 = expression({}^{226}*Ra), k40 = expression({}^{40}*K))


# Data -----------------------------------------------------------------------------------------
dataset <- read_excel("data/processed/soil_radionuclides.xlsx", sheet = "samples") %>%
  select(area, plot, subsample, all_of(predictors)) %>%
  mutate(area = as.factor(area), plot = as.factor(plot), subsample = as.factor(subsample),
         plot_id = interaction(area, plot, drop = TRUE))
cat("Plots:", length(unique(dataset$plot_id)), "| subsamples per plot:",
    paste(unique(table(dataset$plot_id)), collapse = ","), "\n")
print(nearZeroVar(dataset[, predictors], saveMetrics = TRUE))


# Training/test split at the plot level --------------------------------------------------------
set.seed(1653)
plots <- dataset %>% distinct(plot_id, area)
train_plots <- plots %>% group_by(area) %>% slice_sample(prop = 0.7) %>% pull(plot_id)
training <- dataset %>% filter(plot_id %in% train_plots)
testing <- dataset %>% filter(!plot_id %in% train_plots)
stopifnot(length(intersect(training$plot_id, testing$plot_id)) == 0)
split_summary <- data.frame(set = c("training", "test"),
                            plots_A = c(sum(grepl("^A", unique(training$plot_id))), sum(grepl("^A", unique(testing$plot_id)))),
                            plots_B = c(sum(grepl("^B", unique(training$plot_id))), sum(grepl("^B", unique(testing$plot_id)))),
                            subsamples = c(nrow(training), nrow(testing)))
print(split_summary)


# Cross-validation grouped by plot (5 folds x 20 repeats) --------------------------------------
create_grouped_folds <- function(data, group_var, class_var, k = 5, repeats = 20, seed = 123) {
  set.seed(seed)
  groups <- data %>% distinct(group = .data[[group_var]], class = .data[[class_var]])
  index <- list(); indexOut <- list(); counter <- 1
  for (r in seq_len(repeats)) {
    folds <- createFolds(groups$class, k = k, returnTrain = FALSE)
    for (f in seq_along(folds)) {
      test_groups <- groups$group[folds[[f]]]
      train_groups <- setdiff(groups$group, test_groups)
      index[[counter]] <- which(data[[group_var]] %in% train_groups)
      indexOut[[counter]] <- which(data[[group_var]] %in% test_groups)
      names(index)[counter] <- names(indexOut)[counter] <- paste0("Repeat", r, "_Fold", f)
      counter <- counter + 1
    }
  }
  list(index = index, indexOut = indexOut)
}
grouped_folds <- create_grouped_folds(training, "plot_id", "area", k = 5, repeats = 20, seed = 123)
fitControl_grouped <- trainControl(method = "cv", number = 5, index = grouped_folds$index,
                                   indexOut = grouped_folds$indexOut, returnResamp = "all", savePredictions = "final")


# Random Forest (ranger) -----------------------------------------------------------------------
set.seed(1653)
mrf <- train(area ~ ra228 + th228 + ra226 + k40, data = training, method = "ranger",
             trControl = fitControl_grouped, metric = "Accuracy", importance = "permutation", num.trees = 500,
             tuneGrid = expand.grid(mtry = 1:4, splitrule = "gini", min.node.size = c(1, 3, 5)))
print(mrf$bestTune)

# GBM ------------------------------------------------------------------------------------------
set.seed(1653)
mgbm <- train(area ~ ra228 + th228 + ra226 + k40, data = training, method = "gbm",
              trControl = fitControl_grouped, metric = "Accuracy", verbose = FALSE,
              tuneGrid = expand.grid(n.trees = c(50, 100, 250, 500), interaction.depth = 1:3,
                                     shrinkage = c(0.01, 0.05, 0.1), n.minobsinnode = c(5, 10)))
print(mgbm$bestTune)


# Model comparison (cross-validation, best tune of each model) ---------------------------------
res <- resamples(list(RF = mrf, GBM = mgbm))
cv <- res$values
comp <- data.frame(Resample = cv$Resample,
                   RF_Accuracy = cv$`RF~Accuracy`, GBM_Accuracy = cv$`GBM~Accuracy`,
                   RF_Kappa = cv$`RF~Kappa`, GBM_Kappa = cv$`GBM~Kappa`)
cv_summary <- data.frame(model = c("RF", "GBM"),
                         accuracy_mean = c(mean(comp$RF_Accuracy), mean(comp$GBM_Accuracy)),
                         accuracy_sd = c(sd(comp$RF_Accuracy), sd(comp$GBM_Accuracy)),
                         kappa_mean = c(mean(comp$RF_Kappa), mean(comp$GBM_Kappa)),
                         kappa_sd = c(sd(comp$RF_Kappa), sd(comp$GBM_Kappa)))
#Paired Wilcoxon test (same folds for both models)
w_acc <- wilcox.test(comp$RF_Accuracy, comp$GBM_Accuracy, paired = TRUE, exact = FALSE)
w_kap <- wilcox.test(comp$RF_Kappa, comp$GBM_Kappa, paired = TRUE, exact = FALSE)
model_tests <- data.frame(metric = c("Accuracy", "Kappa"),
                          p_paired_wilcoxon = c(w_acc$p.value, w_kap$p.value))
cat("\n===== Grouped cross-validation (100 resamples) =====\n"); print(cv_summary, digits = 3)
print(model_tests, digits = 3)


# Test-set performance (plots not used for fitting) --------------------------------------------
test_results <- do.call(rbind, lapply(list(RF = mrf, GBM = mgbm), function(m) {
  cm <- confusionMatrix(predict(m, newdata = testing), testing$area)
  data.frame(accuracy = cm$overall[["Accuracy"]], ci95_lower = cm$overall[["AccuracyLower"]],
             ci95_upper = cm$overall[["AccuracyUpper"]], kappa = cm$overall[["Kappa"]],
             sensitivity_A = cm$byClass[["Sensitivity"]], specificity_A = cm$byClass[["Specificity"]])
}))
test_results <- data.frame(model = rownames(test_results), test_results, row.names = NULL)
cat("\n===== Test set (", nrow(testing), "subsamples from", length(unique(testing$plot_id)), "plots) =====\n")
print(test_results, digits = 3)
print(confusionMatrix(predict(mrf, newdata = testing), testing$area)$table)


# Variable importance (permutation, Random Forest) ---------------------------------------------
imp_raw <- mrf$finalModel$variable.importance #mean decrease in accuracy by permutation
imp <- data.frame(variable = names(imp_raw), importance = unname(imp_raw),
                  importance_rel = 100 * unname(imp_raw) / max(imp_raw)) %>% arrange(desc(importance))
cat("\n===== Permutation importance (RF) =====\n"); print(imp, digits = 3)


# SHAP values (final RF model with the best hyperparameters) -----------------------------------
set.seed(1177)
rf_final <- ranger(area ~ ra228 + th228 + ra226 + k40, data = training, num.trees = 500,
                   mtry = mrf$bestTune$mtry, min.node.size = mrf$bestTune$min.node.size,
                   probability = TRUE, importance = "permutation")
pred_fun <- function(object, newdata) predict(object, data = newdata)$predictions[, "B"] #probability of Area B
set.seed(123)
shap_values <- fastshap::explain(rf_final, X = as.data.frame(training[, predictors]), pred_wrapper = pred_fun,
                                 nsim = 1000, adjust = TRUE)
shap_imp <- data.frame(variable = colnames(shap_values), mean_abs_shap = colMeans(abs(as.matrix(shap_values)))) %>%
  arrange(desc(mean_abs_shap))
cat("\n===== SHAP importance (mean |SHAP|, probability of Area B) =====\n"); print(shap_imp, digits = 3)


# Figures --------------------------------------------------------------------------------------
#Fig. 8: (A) cross-validation accuracy; (B) permutation importance
acc_long <- rbind(data.frame(model = "RF", accuracy = comp$RF_Accuracy),
                  data.frame(model = "GBM", accuracy = comp$GBM_Accuracy))
acc_long$model <- factor(acc_long$model, levels = c("RF", "GBM"))
gA <- ggplot(acc_long, aes(model, accuracy)) +
  geom_boxplot(fill = "grey85", outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.12, height = 0, size = 0.6, alpha = 0.4) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 2.6, colour = "darkred") +
  labs(x = NULL, y = "Accuracy (grouped cross-validation)", title = "A") +
  theme_pub() + theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"))
imp$variable <- factor(imp$variable, levels = rev(imp$variable))
gB <- ggplot(imp, aes(importance_rel, variable)) +
  geom_segment(aes(x = 0, xend = importance_rel, yend = variable), colour = "grey40") +
  geom_point(size = 2.2) +
  scale_y_discrete(labels = radionuclide_labels) +
  scale_x_continuous(limits = c(0, 105), expand = c(0, 0)) +
  labs(x = "Permutation importance (%)", y = NULL, title = "B") +
  theme_pub() + theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"))
save_pub("Fig8", gA + gB, width_double, 75)

#Fig. A.5: SHAP values (beeswarm)
sv <- shapviz(shap_values, X = as.data.frame(training[, predictors]))
g_shap <- sv_importance(sv, kind = "beeswarm") + theme_pub() + labs(x = "SHAP value (probability of Area B)") +
  scale_y_discrete(labels = radionuclide_labels)
save_pub("FigA5", g_shap, width_mid, 80)


# Export results -------------------------------------------------------------------------------
write_xlsx(list(split = split_summary,
                best_hyperparameters = data.frame(model = c("RF", "GBM"),
                  parameters = c(paste(names(mrf$bestTune), mrf$bestTune, sep = "=", collapse = "; "),
                                 paste(names(mgbm$bestTune), mgbm$bestTune, sep = "=", collapse = "; "))),
                cross_validation = cv_summary, model_tests = model_tests, cv_resamples = comp,
                test_set = test_results, permutation_importance = imp, shap_importance = shap_imp),
           "data/processed/ml_results.xlsx")
