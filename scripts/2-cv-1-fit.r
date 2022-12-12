#' Compare model specifications and engines at AAT genotype prediction
#'
#' The purpose of this script is to evaluate several optimized predictive models
#' of abnormal AAT genotype. The script compares several model specifications
#' (lung and liver history alone or together with either gender or smoking
#' history) and a subset of the engines considered in a previous script
#' (logistic regression, random forest, and nearest neighbor). The same data set
#' is used for all evaluations and comparisons. Models are tuned and evaluated
#' using folded cross-validation on these data. The partitions are stratified by
#' genotype class (abnormal versus normal). Each comparison also includes the
#' current guideline of COPD as the sole indication for screening.

library(tidyverse)
library(tidymodels)

#' Setup

source(here::here("code/settings.r"))

# genotypes to include in analysis
genotype_incl <- read_rds(here::here("data/genotype-incl.rds"))

#' Evaluation measures

# measures to calculate
aatd_met <- metric_set(accuracy, roc_auc, pr_auc)

#' Model specifications

# tuning restrictions
read_rds(here::here("data/aatd-pred.rds")) %>%
  # all predictors from any specification
  select(unique(unlist(sapply(vars_predictors, eval)))) %>%
  ncol() ->
  n_pred
trees_values <- round(10 ^ seq(0, 2, by = .5))
mtry_values <- sort(unique(c(
  round(n_pred ^ (1/3)),
  round(sqrt(n_pred)),
  floor(n_pred / 2)
)))
mtry_values <- unique(mtry_values)
# 10^4 requires too much memory
neighbors_values <- as.integer(10 ^ seq(1L, 3L, by = 1L))

# logistic regression
logistic_reg(penalty = tune(), mixture = 0) %>%
  set_engine("glmnet") %>%
  set_mode("classification") ->
  aatd_lr_spec
grid_regular(penalty(), levels = 6L) ->
  aatd_lr_grid

# random forest
rand_forest(mtry = tune(), trees = tune()) %>%
  set_engine("randomForest") %>%
  set_mode("classification") ->
  aatd_rf_spec
crossing(trees = trees_values, mtry = mtry_values) ->
  aatd_rf_grid

# nearest neighbor
nearest_neighbor(neighbors = tune(), weight_func = tune()) %>%
  set_engine("kknn") %>%
  set_mode("classification") ->
  aatd_nn_spec
grid_regular(weight_func(), levels = 2L) %>%
  mutate(weight_func = fct_inorder(weight_func)) %>%
  crossing(neighbors = neighbors_values) ->
  aatd_nn_grid

# # quadratic SVM model
# svm_poly(degree = 2L, cost = tune(), scale_factor = tune()) %>%
#   set_engine("kernlab") %>%
#   set_mode("classification") ->
#   aatd_svm2_mod
# grid_regular(cost(), scale_factor()) ->
#   aatd_svm2_grid

ii <- if (file.exists(here::here("data/aatd-2-cv-ii.rds"))) {
  read_rds(here::here("data/aatd-2-cv-ii.rds"))
} else {
  c(0L, 0L)
}
aatd_metrics <- if (file.exists(here::here("data/aatd-2-eval-ml.rds"))) {
  read_rds(here::here("data/aatd-2-eval-ml.rds"))
} else {
  tibble()
}

for (i_pred in seq_along(vars_predictors)) {#LOOP
for (i_resp in seq_along(vars_response)) {#LOOP

if (i_pred < ii[[1L]] || (i_pred == ii[[1L]] && i_resp < ii[[2L]])) next

pred <- names(vars_predictors)[[i_pred]]
resp <- names(vars_response)[[i_resp]]

print("--------------------------------")
print(str_c("Predictors: ", pred))
print(str_c("Response: ", resp))
print("--------------------------------")

#' Pre-process data

# load and subset predictors data
read_rds(here::here("data/aatd-pred.rds")) %>%
  # predictors from current specification
  select(record_id, eval(vars_predictors[[i_pred]])) %>%
  # restrict to training set
  semi_join(read_rds(here::here("data/aatd-2-train.rds")), by = "record_id") %>%
  # remove empty factor levels
  mutate(across(where(is.factor), fct_drop)) ->
  aatd_data

# join in responses
read_rds(here::here("data/aatd-resp.rds")) %>%
  select(record_id, genotype) %>%
  semi_join(genotype_incl, by = "genotype") %>%
  # choice of response
  mutate(
    geno_class = ifelse(
      eval(vars_response[[i_resp]]),
      "Abnormal", "Normal"
    ),
    # make abnormal genotype the first factor level (outcome)
    geno_class = factor(geno_class, c("Abnormal", "Normal"))
  ) %>%
  # remove missing responses
  drop_na() %>%
  inner_join(aatd_data, by = "record_id") ->
  aatd_data

#' Specify pre-processing recipes

# prepare regression recipe
recipe(aatd_data, geno_class ~ .) %>%
  # stop treating the ID as a predictor
  #update_role(record_id, new_role = "id variable") %>%
  step_rm(record_id, genotype) %>%
  # drop cases with missing values
  # step_naomit() %>%
  # remove redundant lung & liver categories
  step_rm(ends_with("_none")) %>%
  # remove any variables that are constant within classes (zero variance)
  step_zv(all_predictors()) %>%
  # linearly independent encoding of factors
  step_dummy(all_nominal_predictors(), one_hot = FALSE) %>%
  # binary encoding of logicals
  # step_mutate_at(has_type(match = "logical"), fn = ~ . * 2L - 1L) %>%
  step_mutate_at(has_type(match = "logical"), fn = as.integer) %>%
  prep() ->
  aatd_reg_rec

# prepare numeric recipe
recipe(aatd_data, geno_class ~ .) %>%
  # stop treating the ID as a predictor
  #update_role(record_id, new_role = "id variable") %>%
  step_rm(record_id, genotype) %>%
  # drop cases with missing values
  # step_naomit() %>%
  # one-hot encoding of factors
  step_dummy(all_nominal_predictors(), one_hot = TRUE) %>%
  # binary encoding of logicals
  # step_mutate_at(has_type(match = "logical"), fn = ~ . * 2L - 1L) %>%
  step_mutate_at(has_type(match = "logical"), fn = as.integer) %>%
  prep() ->
  aatd_num_rec

#' Folds

# folds for cross-validation evaluation
aatd_cv <- vfold_cv(aatd_data, v = n_folds_2, strata = geno_class)

#' Logistic regression

# # tune hyperparameters
# workflow() %>%
#   #add_formula(geno_class ~ .) %>%
#   add_recipe(aatd_reg_rec) %>%
#   add_model(aatd_lr_spec) %>%
#   tune_grid(resamples = aatd_cv, grid = aatd_lr_grid, metrics = aatd_met)
# 
# # fit model
# workflow() %>%
#   #add_formula(geno_class ~ .) %>%
#   add_recipe(aatd_reg_rec) %>%
#   add_model(aatd_lr_spec) %>%
#   fit_resamples(resamples = aatd_cv) ->
#   aatd_lr_fit
# 
# # average evaluation
# aatd_lr_fit %>%
#   collect_metrics() %>%
#   mutate(predictors = pred, response = resp, model = "logistic regression") %>%
#   relocate(predictors, response, model) ->
#   aatd_lr_metrics
# 
# # augment results
# aatd_metrics %>%
#   bind_rows(aatd_lr_metrics) ->
#   aatd_metrics

# tune hyperparameters & record metrics (manually, since *tidymodels* won't)
crossing(aatd_lr_grid, transmute(aatd_cv, id, fold = row_number())) %>%
  mutate(metrics = vector(mode = "list", length = nrow(.))) ->
  aatd_lr_met
pb <- progress::progress_bar$new(
  format = "Logistic regression [:bar] :percent",
  total = nrow(aatd_lr_met), clear = FALSE
)
for (i in rev(seq(nrow(aatd_lr_met)))) {
  fold_i <- aatd_lr_met$fold[[i]]
  penalty_i <- aatd_lr_met$penalty[[i]]
  # mixture_i <- aatd_lr_met$mixture[[i]]
  train_i <- bake(aatd_reg_rec, new_data = training(aatd_cv$splits[[fold_i]]))
  test_i <- bake(aatd_reg_rec, new_data = testing(aatd_cv$splits[[fold_i]]))
  logistic_reg(penalty = penalty_i, mixture = 0) %>%
    set_engine("glmnet") %>%
    set_mode("classification") %>%
    fit(geno_class ~ ., train_i) ->
    fit_i
  bind_cols(
    select(test_i, class = geno_class),
    predict(fit_i, new_data = test_i),
    predict(fit_i, new_data = test_i, type = "prob")
  ) %>%
    metrics(truth = class, estimate = .pred_class, .pred_Abnormal) ->
    met_i
  aatd_lr_met$metrics[[i]] <- met_i
  pb$tick()
}
aatd_lr_met %>%
  mutate(
    model = "logistic regression",
    predictors = pred, response = resp
  ) %>%
  select(-fold) %>%
  unnest(metrics) %>%
  nest(hyperparameters = c(penalty)) ->
  aatd_lr_met
aatd_metrics <- bind_rows(aatd_metrics, aatd_lr_met)

#' Random forests

# # fit model
# workflow() %>%
#   #add_formula(geno_class ~ .) %>%
#   add_recipe(aatd_num_rec) %>%
#   add_model(aatd_rf_spec) %>%
#   fit_resamples(resamples = aatd_cv) ->
#   aatd_rf_fit
# 
# # average evaluation
# aatd_rf_fit %>%
#   collect_metrics() %>%
#   mutate(predictors = pred, response = resp, model = "random forest") %>%
#   relocate(predictors, response, model) ->
#   aatd_rf_metrics
# 
# # augment results
# aatd_metrics %>%
#   bind_rows(aatd_rf_metrics) ->
#   aatd_metrics

# tune hyperparameters & record metrics
crossing(aatd_rf_grid, transmute(aatd_cv, id, fold = row_number())) %>%
  mutate(metrics = vector(mode = "list", length = nrow(.))) ->
  aatd_rf_met
pb <- progress::progress_bar$new(
  format = "Random forest [:bar] :percent",
  total = nrow(aatd_rf_met), clear = FALSE
)
for (i in rev(seq(nrow(aatd_rf_met)))) {
  fold_i <- aatd_rf_met$fold[[i]]
  mtry_i <- aatd_rf_met$mtry[[i]]
  trees_i <- aatd_rf_met$trees[[i]]
  train_i <- bake(aatd_num_rec, new_data = training(aatd_cv$splits[[fold_i]]))
  test_i <- bake(aatd_num_rec, new_data = testing(aatd_cv$splits[[fold_i]]))
  rand_forest(mtry = mtry_i, trees = trees_i) %>%
    set_engine("randomForest") %>%
    set_mode("classification") %>%
    fit(geno_class ~ ., train_i) ->
    fit_i
  # TODO: debug `predict()` calls for `i >= 13L`
  bind_cols(
    select(test_i, class = geno_class),
    predict(fit_i, new_data = test_i),
    predict(fit_i, new_data = test_i, type = "prob")
  ) %>%
    metrics(truth = class, estimate = .pred_class, .pred_Abnormal) ->
    met_i
  aatd_rf_met$metrics[[i]] <- met_i
  pb$tick()
}
aatd_rf_met %>%
  mutate(
    model = "random forest",
    predictors = pred, response = resp
  ) %>%
  select(-fold) %>%
  unnest(metrics) %>%
  nest(hyperparameters = c(mtry, trees)) ->
  aatd_rf_met
aatd_metrics <- bind_rows(aatd_metrics, aatd_rf_met)

#' Nearest neighbors

# # fit model
# workflow() %>%
#   #add_formula(geno_class ~ .) %>%
#   add_recipe(aatd_num_rec) %>%
#   add_model(aatd_nn_spec) %>%
#   fit_resamples(resamples = aatd_cv) ->
#   aatd_nn_fit
# 
# # average evaluation
# aatd_nn_fit %>%
#   collect_metrics() %>%
#   mutate(predictors = pred, response = resp, model = "nearest neighbor") %>%
#   relocate(predictors, response, model) ->
#   aatd_nn_metrics
# 
# # augment results
# aatd_metrics %>%
#   bind_rows(aatd_nn_metrics) ->
#   aatd_metrics

# tune hyperparameters & record metrics
crossing(aatd_nn_grid, transmute(aatd_cv, id, fold = row_number())) %>%
  mutate(metrics = vector(mode = "list", length = nrow(.))) ->
  aatd_nn_met
pb <- progress::progress_bar$new(
  format = "Nearest neighbors [:bar] :percent",
  total = nrow(aatd_nn_met), clear = FALSE
)
for (i in rev(seq(nrow(aatd_nn_met)))) {
  fold_i <- aatd_nn_met$fold[[i]]
  neighbors_i <- aatd_nn_met$neighbors[[i]]
  weight_func_i <- as.character(aatd_nn_met$weight_func[[i]])
  train_i <- bake(aatd_num_rec, new_data = training(aatd_cv$splits[[fold_i]]))
  test_i <- bake(aatd_num_rec, new_data = testing(aatd_cv$splits[[fold_i]]))
  nearest_neighbor(neighbors = neighbors_i, weight_func = weight_func_i) %>%
    set_engine("kknn") %>%
    set_mode("classification") %>%
    fit(geno_class ~ ., train_i) ->
    fit_i
  bind_cols(
    select(test_i, class = geno_class),
    predict(fit_i, new_data = test_i),
    predict(fit_i, new_data = test_i, type = "prob")
  ) %>%
    metrics(truth = class, estimate = .pred_class, .pred_Abnormal) ->
    met_i
  aatd_nn_met$metrics[[i]] <- met_i
  pb$tick()
}
aatd_nn_met %>%
  mutate(
    model = "nearest neighbors",
    predictors = pred, response = resp
  ) %>%
  select(-fold) %>%
  unnest(metrics) %>%
  nest(hyperparameters = c(neighbors, weight_func)) ->
  aatd_nn_met
aatd_metrics <- bind_rows(aatd_metrics, aatd_nn_met)

write_rds(aatd_metrics, here::here("data/aatd-2-eval-ml.rds"))
write_rds(c(i_pred, i_resp), here::here("data/aatd-2-cv-ii.rds"))

}#LOOP
}#LOOP

file.remove(here::here("data/aatd-2-cv-ii.rds"))
