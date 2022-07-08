library(tidyverse)
library(tidymodels)

#' 0. Setup

# hyperparameters
p_data <- 1/6
n_train <- 2/3
n_folds <- 6L
#n_folds <- 12L

# load and subset predictors data
read_rds(here::here("data/aatd-pred.rds")) %>%
  # exclude variables never used as predictors
  select(-race) %>%
  sample_frac(size = p_data) %>%
  print() ->
  aatd_data

# genotypes to include in analysis
read_rds(here::here("data/aatd-resp.rds")) %>%
  count(genotype) %>%
  filter(n > 120L) %T>% print() %>%
  select(genotype) %>%
  drop_na() ->
  genotype_incl

# join in responses
read_rds(here::here("data/aatd-resp.rds")) %>%
  select(record_id, genotype) %>%
  semi_join(genotype_incl, by = "genotype") %>%
  mutate(
    geno_class = ifelse(
      #genotype == "ZZ",
      genotype == "ZZ" | genotype == "SZ",
      #str_detect(genotype, "Z$"),
      "Abnormal", "Normal"
    ),
    # make abnormal genotype the first factor level (outcome)
    geno_class = factor(geno_class, c("Abnormal", "Normal"))
  ) %>%
  # remove missing responses
  drop_na() %>%
  inner_join(aatd_data, by = "record_id") ->
  aatd_data

#' 0.1 Model specifications

# logistic regression, fix parameters
logistic_reg(penalty = 1) %>%
  set_engine("glm") %>%
  set_mode("classification") ->
  aatd_lr_spec

# logistic regression, tune parameters
logistic_reg(penalty = tune()) %>%
  set_engine("glm") %>%
  set_mode("classification") ->
  aatd_lr_spec_tune

# random forest, fix parameters
rand_forest(mtry = NULL, trees = 120L) %>%
  set_engine("randomForest") %>%
  set_mode("classification") ->
  aatd_rf_spec

# linear SVM, fix parameters
svm_linear() %>%
  #set_engine("LiblineaR") %>%
  set_engine("kernlab") %>%
  set_mode("classification") ->
  aatd_svm1_spec

# quadratic SVM model specification
svm_poly(degree = 2L) %>%
  set_engine("kernlab") %>%
  set_mode("classification") ->
  aatd_svm2_spec

# cubic SVM model specification
svm_poly(degree = 3L) %>%
  set_engine("kernlab") %>%
  set_mode("classification") ->
  aatd_svm3_spec

#' 1. Single-partition, fixed-parameter

# initial partition
aatd_split <- initial_split(aatd_data, prop = n_train, strata = geno_class)

# prepare regression recipe
recipe(training(aatd_split), geno_class ~ .) %>%
  # stop treating the ID as a predictor
  #update_role(record_id, new_role = "id variable") %>%
  step_rm(record_id, genotype) %>%
  prep() %>%
  print() ->
  aatd_reg_rec

# prepare numeric recipe
recipe(training(aatd_split), geno_class ~ .) %>%
  # stop treating the ID as a predictor
  #update_role(record_id, new_role = "id variable") %>%
  step_rm(record_id, genotype) %>%
  # one-hot encoding of factors
  step_dummy(all_nominal_predictors(), one_hot = TRUE) %>%
  # binary encoding of logicals
  step_mutate_at(has_type(match = "logical"), fn = ~ . * 2L - 1L) %>%
  prep() %>%
  print() ->
  aatd_num_rec

# summarize results from model fit
fit_results <- function(fit) {
  bind_cols(
    predict(fit, bake(aatd_num_rec, testing(aatd_split))),
    # -+- manually add any missing classes -+-
    predict(fit, bake(aatd_num_rec, testing(aatd_split)), type = "prob"),
    select(testing(aatd_split), geno_class)
  )
}

#' 1.1. Logistic regression

# fit model
aatd_lr_spec %>%
  fit(geno_class ~ ., bake(aatd_reg_rec, NULL)) %>%
  print() ->
  aatd_lr_fit

# evaluate model
aatd_lr_res <- fit_results(aatd_lr_fit)
aatd_lr_res %>%
  count(.pred_class, geno_class)
aatd_lr_res %>%
  metrics(truth = geno_class, estimate = .pred_class, .pred_Abnormal)

#' 1.2. Random forest classification

# fit model
aatd_rf_spec %>%
  fit(geno_class ~ ., bake(aatd_num_rec, NULL)) %>%
  print() ->
  aatd_rf_fit

# evaluate model
aatd_rf_res <- fit_results(aatd_rf_fit)
aatd_rf_res %>%
  count(.pred_class, geno_class)
aatd_rf_res %>%
  metrics(truth = geno_class, estimate = .pred_class, .pred_Abnormal)

#' 1.3. Support vector machine classification

# fit model
aatd_svm1_spec %>%
  fit(geno_class ~ ., bake(aatd_num_rec, NULL)) %>%
  print() ->
  aatd_svm1_fit

# evaluate model
aatd_svm1_res <- fit_results(aatd_svm1_fit)
aatd_svm1_res %>%
  count(.pred_class, geno_class)
aatd_svm1_res %>%
  metrics(truth = geno_class, estimate = .pred_class, .pred_Abnormal)

# fit model
aatd_svm2_spec %>%
  fit(geno_class ~ ., bake(aatd_num_rec, NULL)) %>%
  print() ->
  aatd_svm2_fit

# evaluate model
aatd_svm2_res <- fit_results(aatd_svm2_fit)
aatd_svm2_res %>%
  count(.pred_class, geno_class)
aatd_svm2_res %>%
  metrics(truth = geno_class, estimate = .pred_class, .pred_Abnormal)

# fit model
aatd_svm3_spec %>%
  fit(geno_class ~ ., bake(aatd_num_rec, NULL)) %>%
  print() ->
  aatd_svm3_fit

# evaluate model
aatd_svm3_res <- fit_results(aatd_svm3_fit)
aatd_svm3_res %>%
  count(.pred_class, geno_class)
aatd_svm3_res %>%
  metrics(truth = geno_class, estimate = .pred_class, .pred_Abnormal)

# compare ROC curves
list(
  logistic_regression = aatd_lr_res,
  random_forest = aatd_rf_res,
  linear_svm = aatd_svm1_res,
  quadratic_svm = aatd_svm2_res,
  cubic_svm = aatd_svm3_res
) %>%
  enframe(name = "model", value = "results") %>%
  mutate(model = fct_inorder(model)) %>%
  unnest(results) %>%
  group_by(model) %>%
  roc_curve(truth = geno_class, estimate = .pred_Abnormal) %>%
  autoplot() ->
  aatd_roc
print(aatd_roc)
ggsave(here::here("fig/aatd-roc.png"), aatd_roc)

# compare PR curves
list(
  logistic_regression = aatd_lr_res,
  random_forest = aatd_rf_res,
  linear_svm = aatd_svm1_res,
  quadratic_svm = aatd_svm2_res,
  cubic_svm = aatd_svm3_res
) %>%
  enframe(name = "model", value = "results") %>%
  mutate(model = fct_inorder(model)) %>%
  unnest(results) %>%
  group_by(model) %>%
  pr_curve(truth = geno_class, estimate = .pred_Abnormal) %>%
  autoplot() ->
  aatd_pr
print(aatd_pr)
ggsave(here::here("fig/aatd-pr.png"), aatd_pr)

# compare models at optimal trade-off between sensitivity and specificity?

stop()

#' 2. Single-partition, CV-optimized



#' 3. CV-evaluated, CV-optimized

# partition for machine learning
aatd_cv_outer <- vfold_cv(aatd_data, v = n_folds, strata = geno_class)

# across folds
i <- 1L

# prepare recipe
recipe(training(aatd_cv_outer$splits[[i]]), geno_class ~ .) %>%
  # stop treating the ID as a predictor
  #update_role(record_id, new_role = "id variable") %>%
  step_rm(record_id) %>%
  prep() %>%
  print() ->
  aatd_rec

# resamples
training(aatd_cv_outer$splits[[i]]) %>%
  vfold_cv(v = n_folds, strata = geno_class) ->
  aatd_cv_inner
# check numbers of cases in training sets
aatd_cv_inner %>%
  mutate(cases = map_int(
    splits,
    ~ nrow(filter(analysis(.), geno_class == "Abnormal"))
  ))

# logistic regression parameter ranges
grid_regular(
  penalty(),
  levels = 6L
) ->
  aatd_lr_grid

# tune logistic regression model
workflow() %>%
  add_model(aatd_lr_spec_tune) %>%
  #add_formula(geno_class ~ .) %>%
  add_recipe(aatd_rec) %>%
  print() ->
  aatd_lr_wf
aatd_lr_wf %>%
  tune_grid(resamples = aatd_cv_inner, grid = aatd_lr_grid) ->
  aatd_lr_tune
aatd_lr_tune %>%
  collect_metrics() ->
  aatd_lr_metrics



# random forest parameter ranges
n_pred <- ncol(bake(aatd_rec, new_data = NULL)) - 1L
mtry_values <- c(
  round(n_pred ^ (1/3)),
  round(sqrt(n_pred)),
  round(n_pred / 2),
  n_pred
)
grid_regular(
  trees(),
  levels = 6L
) %>%
  crossing(mtry = mtry_values) ->
  aatd_rf_grid

# support vector machine parameter ranges
grid_regular(
  cost(),
  levels = 6L
) %>%
  crossing(degree = c(1L, 2L, 3L, 4L)) ->
  aatd_svm_grid
