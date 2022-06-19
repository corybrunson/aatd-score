library(tidyverse)
library(tidymodels)

# load and subset predictors data
read_rds(here::here("data/aatd-pred.rds")) %>%
  select(-race) %>%
  sample_frac(size = p_data) %>%
  print() ->
  aatd_data
# join in responses
read_rds(here::here("data/aatd-resp.rds")) %>%
  select(record_id, genotype) %>%
  mutate(
    genotype = ifelse(
      #genotype == "ZZ",
      genotype == "ZZ" | genotype == "SZ",
      "Abnormal", "Normal"
    ),
    genotype = fct_infreq(genotype)
  ) %>%
  # remove missing responses
  drop_na() %>%
  inner_join(aatd_data, by = "record_id") ->
  aatd_data

# initial partition
aatd_split <- initial_split(aatd_data, prop = 2/3, strata = genotype)

# prepare recipe
recipe(training(aatd_split), genotype ~ .) %>%
  # stop treating the ID as a predictor
  #update_role(record_id, new_role = "id variable") %>%
  step_rm(record_id) %>%
  prep() %>%
  print() ->
  aatd_rec

# logistic regression model specification
logistic_reg(penalty = 1) %>%
  set_engine("glm") %>%
  set_mode("classification") ->
  aatd_lr_spec

# fit model
aatd_lr_spec %>%
  fit(genotype ~ ., bake(aatd_rec, NULL)) %>%
  print() ->
  aatd_lr_fit

# evaluate model
bind_cols(
  predict(aatd_lr_fit, bake(aatd_rec, testing(aatd_split))),
  predict(aatd_lr_fit, bake(aatd_rec, testing(aatd_split)), type = "prob")
) %>%
  bind_cols(select(testing(aatd_split), genotype)) %>%
  print() ->
  aatd_lr_res
aatd_lr_res %>%
  count(.pred_class, genotype)
aatd_lr_res %>%
  metrics(truth = genotype, estimate = .pred_class, .pred_Normal)

stop()

# hyperparameters
p_data <- 1/100
n_folds <- 6L
#n_folds <- 12L

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
