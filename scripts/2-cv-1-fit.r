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

# hyperparameters
# p_data <- 1/10
# p_data <- 1/6
p_data <- 1
n_nests <- 6L
# n_folds <- 3L
n_folds <- 6L

# genotypes to include in analysis
genotype_incl <- read_rds(here::here("data/genotype-incl.rds"))

# model specifications as tidy selections
vars_predictors <- list(
  Dx = expr(c(starts_with("lung_"), starts_with("liver_"))),
  # `Dx+age` =
  #   expr(c(contains("age_guess"), contains("receipt_date"),
  #          starts_with("lung_"), starts_with("liver_"))),
  `Dx+gender` =
    expr(c(contains("gender"), starts_with("lung_"), starts_with("liver_"))),
  `Dx+tobacco` =
    expr(c(contains("smoking_history_cigarette"),
           # contains("any_tobacco_exposure"),
           starts_with("lung_"), starts_with("liver_")))
)

# response variables as logical tests
vars_response <- list(
  ZZ = expr(genotype == "ZZ"),
  SZ = expr(genotype == "SZ" | genotype == "ZZ"),
  # MZ = expr(genotype == "MZ" | genotype == "SZ" | genotype == "ZZ"),
  # Z = expr(grepl(".*Z$", genotype)),
  Ab = expr(genotype != "MM")
)

#' Subset data

read_rds(here::here("data/aatd-pred.rds")) %>%
  sample_frac(size = p_data) %>%
  # all predictors from any specification
  select(record_id, unique(unlist(sapply(vars_predictors, eval)))) %>%
  # filter missing gender
  filter(gender != "(Missing)") %>%
  # drop any cases with missing values
  drop_na() %>%
  # store `record_id`
  select(record_id) ->
  elig_ids

#' Model specifications

# logistic regression
logistic_reg(penalty = tune()) %>%
  set_engine("glm") %>%
  set_mode("classification") ->
  aatd_lr_spec

# random forest
rand_forest(mtry = NULL, trees = 120L) %>%
  set_engine("randomForest") %>%
  set_mode("classification") ->
  aatd_rf_spec

# nearest neighbor
nearest_neighbor(neighbors = 360L, weight_func = "triangular") %>%
  set_engine("kknn") %>%
  set_mode("classification") ->
  aatd_nn_spec

# progress bar
n_loop <- length(vars_predictors) * length(vars_response)
pb <- progress::progress_bar$new(total = n_loop)

ii <- if (file.exists(here::here("data/aatd-cv-ii.rds"))) {
  read_rds(here::here("data/aatd-cv-ii.rds"))
} else {
  c(0L, 0L)
}
aatd_metrics <- if (file.exists(here::here("data/aatd-eval.rds"))) {
  read_rds(here::here("data/aatd-eval.rds"))
} else {
  tibble()
}

for (i_pred in seq_along(vars_predictors)) {#LOOP
for (i_resp in seq_along(vars_response)) {#LOOP

pb$tick()
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
  sample_frac(size = p_data) %>%
  # predictors from current specification
  select(record_id, eval(vars_predictors[[i_pred]])) %>%
  # eligible records
  semi_join(elig_ids, by = "record_id") %>%
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

#' Specify models

# prepare regression recipe
recipe(aatd_data, geno_class ~ .) %>%
  # stop treating the ID as a predictor
  #update_role(record_id, new_role = "id variable") %>%
  step_rm(record_id, genotype, ends_with("_none")) %>%
  prep() ->
  aatd_reg_rec

# prepare numeric recipe
recipe(aatd_data, geno_class ~ .) %>%
  # stop treating the ID as a predictor
  #update_role(record_id, new_role = "id variable") %>%
  step_rm(record_id, genotype) %>%
  # one-hot encoding of factors
  step_dummy(all_nominal_predictors(), one_hot = TRUE) %>%
  # binary encoding of logicals
  step_mutate_at(has_type(match = "logical"), fn = ~ . * 2L - 1L) %>%
  prep() ->
  aatd_num_rec

#' Folds

# folds for cross-validation evaluation
aatd_cv <- vfold_cv(aatd_data, v = n_folds, strata = geno_class)

#' Logistic regression

# tune hyperparameters
tune_grid(aatd_lr_spec, aatd_reg_rec, resamples = aatd_cv)
workflow() %>%
  #add_formula(geno_class ~ .) %>%
  add_recipe(aatd_reg_rec) %>%
  add_model(aatd_lr_spec) %>%
  tune_grid(resamples = aatd_cv)

# fit model
workflow() %>%
  #add_formula(geno_class ~ .) %>%
  add_recipe(aatd_reg_rec) %>%
  add_model(aatd_lr_spec) %>%
  fit_resamples(resamples = aatd_cv) ->
  aatd_lr_fit

# average evaluation
aatd_lr_fit %>%
  collect_metrics() %>%
  mutate(predictors = pred, response = resp, model = "logistic regression") %>%
  relocate(predictors, response, model) ->
  aatd_lr_metrics

# augment results
aatd_metrics %>%
  bind_rows(aatd_lr_metrics) ->
  aatd_metrics

#' Random forests

# fit model
workflow() %>%
  #add_formula(geno_class ~ .) %>%
  add_recipe(aatd_num_rec) %>%
  add_model(aatd_rf_spec) %>%
  fit_resamples(resamples = aatd_cv) ->
  aatd_rf_fit

# average evaluation
aatd_rf_fit %>%
  collect_metrics() %>%
  mutate(predictors = pred, response = resp, model = "random forest") %>%
  relocate(predictors, response, model) ->
  aatd_rf_metrics

# augment results
aatd_metrics %>%
  bind_rows(aatd_rf_metrics) ->
  aatd_metrics

#' Nearest neighbors

# fit model
workflow() %>%
  #add_formula(geno_class ~ .) %>%
  add_recipe(aatd_num_rec) %>%
  add_model(aatd_nn_spec) %>%
  fit_resamples(resamples = aatd_cv) ->
  aatd_nn_fit

# average evaluation
aatd_nn_fit %>%
  collect_metrics() %>%
  mutate(predictors = pred, response = resp, model = "nearest neighbor") %>%
  relocate(predictors, response, model) ->
  aatd_nn_metrics

# augment results
aatd_metrics %>%
  bind_rows(aatd_nn_metrics) ->
  aatd_metrics

write_rds(aatd_metrics, here::here("data/aatd-eval.rds"))
write_rds(c(i_pred, i_resp), here::here("data/aatd-cv-ii.rds"))

}#LOOP
}#LOOP

