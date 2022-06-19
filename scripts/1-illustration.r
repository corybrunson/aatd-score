library(tidyverse)
library(tidymodels)

#' 0. Setup

# hyperparameters
p_data <- 1/10
n_folds <- 6L
#n_folds <- 12L

# load and subset predictors data
read_rds(here::here("data/aatd-pred.rds")) %>%
  select(-race) %>%
  sample_frac(size = p_data) %>%
  print() ->
  aatd_data
# genotypes to include in analysis
read_rds(here::here("data/aatd-resp.rds")) %>%
  count(genotype) %>%
  filter(n > 120L) %>%
  select(genotype) %>%
  drop_na() ->
  genotype_incl
# join in responses
read_rds(here::here("data/aatd-resp.rds")) %>%
  select(record_id, genotype) %>%
  semi_join(genotype_incl) %>%
  mutate(
    genotype = ifelse(
      #genotype == "ZZ",
      genotype == "ZZ" | genotype == "SZ",
      "Abnormal", "Normal"
    ),
    # make abnormal genotype the first factor level (outcome)
    genotype = fct_rev(fct_infreq(genotype))
  ) %>%
  # remove missing responses
  drop_na() %>%
  inner_join(aatd_data, by = "record_id") ->
  aatd_data

#' 1. Single-partition, fixed-parameter

# initial partition
aatd_split <- initial_split(aatd_data, prop = 2/3, strata = genotype)

# prepare regression recipe
recipe(training(aatd_split), genotype ~ .) %>%
  # stop treating the ID as a predictor
  #update_role(record_id, new_role = "id variable") %>%
  step_rm(record_id) %>%
  prep() %>%
  print() ->
  aatd_reg_rec

# prepare numeric recipe
recipe(training(aatd_split), genotype ~ .) %>%
  # stop treating the ID as a predictor
  #update_role(record_id, new_role = "id variable") %>%
  step_rm(record_id) %>%
  # one-hot encoding of factors
  step_dummy(all_nominal_predictors(), one_hot = TRUE) %>%
  # binary encoding of logicals
  step_mutate_at(has_type(match = "logical"), fn = ~ . * 2L - 1L) %>%
  prep() %>%
  print() ->
  aatd_num_rec

#' 1.1. Logistic regression

# logistic regression model specification
logistic_reg(penalty = 1) %>%
  set_engine("glm") %>%
  set_mode("classification") ->
  aatd_lr_spec

# fit model
aatd_lr_spec %>%
  fit(genotype ~ ., bake(aatd_reg_rec, NULL)) %>%
  print() ->
  aatd_lr_fit

# evaluate model
bind_cols(
  predict(aatd_lr_fit, bake(aatd_reg_rec, testing(aatd_split))),
  predict(aatd_lr_fit, bake(aatd_reg_rec, testing(aatd_split)), type = "prob")
) %>%
  bind_cols(select(testing(aatd_split), genotype)) %>%
  print() ->
  aatd_lr_res
aatd_lr_res %>%
  count(.pred_class, genotype)
aatd_lr_res %>%
  metrics(truth = genotype, estimate = .pred_class, .pred_Abnormal)

#' 1.2. Random forest classification

# random forest model specification
rand_forest(mtry = NULL, trees = 120L) %>%
  set_engine("randomForest") %>%
  set_mode("classification") ->
  aatd_rf_spec

# fit model
aatd_rf_spec %>%
  fit(genotype ~ ., bake(aatd_num_rec, NULL)) %>%
  print() ->
  aatd_rf_fit

# evaluate model
bind_cols(
  predict(aatd_rf_fit, bake(aatd_num_rec, testing(aatd_split))),
  predict(aatd_rf_fit, bake(aatd_num_rec, testing(aatd_split)), type = "prob")
) %>%
  bind_cols(select(testing(aatd_split), genotype)) %>%
  print() ->
  aatd_rf_res
aatd_rf_res %>%
  count(.pred_class, genotype)
aatd_rf_res %>%
  metrics(truth = genotype, estimate = .pred_class, .pred_Abnormal)

#' 1.3. Support vector machine classification

# linear SVM model specification
svm_linear() %>%
  #set_engine("LiblineaR") %>%
  set_engine("kernlab") %>%
  set_mode("classification") ->
  aatd_svm1_spec

# fit model
aatd_svm1_spec %>%
  fit(genotype ~ ., bake(aatd_num_rec, NULL)) %>%
  print() ->
  aatd_svm1_fit

# evaluate model
bind_cols(
  predict(aatd_svm1_fit, bake(aatd_num_rec, testing(aatd_split))),
  predict(aatd_svm1_fit, bake(aatd_num_rec, testing(aatd_split)), type = "prob")
) %>%
  bind_cols(select(testing(aatd_split), genotype)) %>%
  print() ->
  aatd_svm1_res
aatd_svm1_res %>%
  count(.pred_class, genotype)
aatd_svm1_res %>%
  metrics(truth = genotype, estimate = .pred_class, .pred_Abnormal)

# quadratic SVM model specification
svm_poly(degree = 2L) %>%
  set_engine("kernlab") %>%
  set_mode("classification") ->
  aatd_svm2_spec

# fit model
aatd_svm2_spec %>%
  fit(genotype ~ ., bake(aatd_num_rec, NULL)) %>%
  print() ->
  aatd_svm2_fit

# evaluate model
bind_cols(
  predict(aatd_svm2_fit, bake(aatd_num_rec, testing(aatd_split))),
  predict(aatd_svm2_fit, bake(aatd_num_rec, testing(aatd_split)), type = "prob")
) %>%
  bind_cols(select(testing(aatd_split), genotype)) %>%
  print() ->
  aatd_svm2_res
aatd_svm2_res %>%
  count(.pred_class, genotype)
aatd_svm2_res %>%
  metrics(truth = genotype, estimate = .pred_class, .pred_Abnormal)

# cubic SVM model specification
svm_poly(degree = 3L) %>%
  set_engine("kernlab") %>%
  set_mode("classification") ->
  aatd_svm3_spec

# fit model
aatd_svm3_spec %>%
  fit(genotype ~ ., bake(aatd_num_rec, NULL)) %>%
  print() ->
  aatd_svm3_fit

# evaluate model
bind_cols(
  predict(aatd_svm3_fit, bake(aatd_num_rec, testing(aatd_split))),
  predict(aatd_svm3_fit, bake(aatd_num_rec, testing(aatd_split)), type = "prob")
) %>%
  bind_cols(select(testing(aatd_split), genotype)) %>%
  print() ->
  aatd_svm3_res
aatd_svm3_res %>%
  count(.pred_class, genotype)
aatd_svm3_res %>%
  metrics(truth = genotype, estimate = .pred_class, .pred_Abnormal)

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
  roc_curve(truth = genotype, estimate = .pred_Abnormal) %>%
  autoplot() ->
  aatd_roc
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
  pr_curve(truth = genotype, estimate = .pred_Abnormal) %>%
  autoplot()

stop()

#' 2. Single-partition, CV-optimized



#' 3. CV-evaluated, CV-optimized

# partition for machine learning
aatd_cv_outer <- vfold_cv(aatd_data, v = n_folds, strata = genotype)

# across folds
i <- 1L

# prepare recipe
recipe(training(aatd_cv_outer$splits[[i]]), genotype ~ .) %>%
  # stop treating the ID as a predictor
  #update_role(record_id, new_role = "id variable") %>%
  step_rm(record_id) %>%
  prep() %>%
  print() ->
  aatd_rec

# resamples
training(aatd_cv_outer$splits[[i]]) %>%
  vfold_cv(v = n_folds, strata = genotype) ->
  aatd_cv_inner
# check numbers of cases in training sets
aatd_cv_inner %>%
  mutate(cases = map_int(
    splits,
    ~ nrow(filter(analysis(.), genotype == "Abnormal"))
  ))

# logistic regression model specification
logistic_reg(penalty = tune()) %>%
  set_engine("glm") %>%
  set_mode("classification") ->
  aatd_lr_spec

# logistic regression parameter ranges
grid_regular(
  penalty(),
  levels = 6L
) ->
  aatd_lr_grid

# tune logistic regression model
workflow() %>%
  add_model(aatd_lr_spec) %>%
  #add_formula(genotype ~ .) %>%
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
