library(tidyverse)
library(tidymodels)
library(reticulate)
np <- import("numpy")
fr <- import("fasterrisk.fasterrisk")

#' Setup

source(here::here("code/settings.r"))

# genotypes to include in analysis
genotype_incl <- read_rds(here::here("data/genotype-incl.rds"))

#' Evaluation measures

# measures to calculate
aatd_met <- metric_set(accuracy, roc_auc, pr_auc)

#' Model specifications

# `FasterRisk` settings
# n_terms <- 7L
n_models <- 2L
# abs_bound <- 5
n_retain <- 8L
n_mults <- 16L

# tuning restrictions
read_rds(here::here("data/aatd-pred.rds")) %>%
  # all predictors from any specification
  select(unique(unlist(sapply(vars_predictors, eval)))) %>%
  ncol() ->
  n_pred
ns_terms <- seq(5L, ceiling(3/4 * n_pred), length.out = 3L)
abs_bounds <- c(2, 6, 10)

# read in existing data
aatd_metrics <- if (file.exists(here::here("data/aatd-2-eval-fr.rds"))) {
  here::here("data/aatd-2-eval-fr.rds") %>%
    read_rds() %>%
    group_by(predictors, response, .metric) %>%
    add_count(name = "count") %>%
    ungroup() %>%
    filter(count == max(count))
} else {
  tibble()
}

# check assumption
if (nrow(aatd_metrics) > 0L)
  stopifnot(all(aatd_metrics$count ==
                  n_folds_2 * n_models * length(ns_terms) * length(abs_bounds)))

for (i_pred in seq_along(vars_predictors)) {#LOOP
for (i_resp in seq_along(vars_response)) {#LOOP

pred <- names(vars_predictors)[[i_pred]]
resp <- names(vars_response)[[i_resp]]

# skip this loop if already done
done <- if (nrow(aatd_metrics) == 0L) FALSE else {
  aatd_metrics %>%
    filter(predictors == pred & response == resp) %>%
    nrow() %>%
    as.logical()
}
if (done) next

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

# prepare binarization recipe
recipe(aatd_data, geno_class ~ .) %>%
  # stop treating the ID as a predictor
  #update_role(record_id, new_role = "id variable") %>%
  step_rm(record_id, genotype) %>%
  # remove redundant lung & liver categories
  step_rm(ends_with("_none")) %>%
  # linearly independent encoding of factors
  step_dummy(all_nominal_predictors(), one_hot = FALSE) %>%
  # binary encoding of logicals
  step_mutate_at(has_type(match = "logical"), fn = as.integer) %>%
  # -1/1 encoding of response
  step_mutate_at(
    has_role(match = "outcome"),
    fn = ~ ifelse(. == "Abnormal", -1L, 1L)
  ) %>%
  prep() ->
  aatd_int_rec

#' Folds

# folds for cross-validation evaluation
aatd_cv <- vfold_cv(aatd_data, v = n_folds_2, strata = geno_class)

#' FasterRisk

transmute(aatd_cv, id) %>%
  mutate(fold = seq(n_folds_2)) %>%
  crossing(terms = ns_terms, bound = abs_bounds) %>%
  mutate(metrics = vector(mode = "list", length = nrow(.))) ->
  aatd_fr_met

for (r in seq(nrow(aatd_fr_met))) {#LOOP
# r <- 1L

# set parameters from data frame
i_fold <- aatd_fr_met$fold[[r]]
n_terms <- aatd_fr_met$terms[[r]]
abs_bound <- aatd_fr_met$bound[[r]]

# obtain training and testing sets
aatd_train <- bake(aatd_int_rec, training(aatd_cv$splits[[i_fold]]))
aatd_test <- bake(aatd_int_rec, testing(aatd_cv$splits[[i_fold]]))

# convert to `numpy` arrays
aatd_train %>%
  select(-geno_class) %>%
  np$asarray() %>%
  np_array(dtype = "int") ->
  X_train
aatd_train %>%
  pull(geno_class) %>%
  np$asarray() %>%
  np_array(dtype = "int") ->
  y_train
aatd_test %>%
  select(-geno_class) %>%
  np$asarray() %>%
  np_array(dtype = "int") ->
  X_test
aatd_test %>%
  pull(geno_class) %>%
  np$asarray() %>%
  np_array(dtype = "int") ->
  y_test

# specify model
aatd_rso <- fr$RiskScoreOptimizer(
  X = X_train, y = y_train,
  k = n_terms, select_top_m = n_models,
  lb = -abs_bound, ub = abs_bound,
  parent_size = n_retain, num_ray_search = n_mults
)

# optimize model
aatd_rso$optimize()

# save results and destroy optimizer
aatd_rso_res <- aatd_rso$get_models()
rm(aatd_rso)

#' Evaluation and comparison of models

aatd_fr_fold_met <- tibble()
for (i_model in seq(n_models)) {
  
  # build classifier
  aatd_rsc <- fr$RiskScoreClassifier(
    multiplier = aatd_rso_res[[1L]][[i_model]],
    intercept = aatd_rso_res[[2L]][i_model],
    coefficients = np$asarray(aatd_rso_res[[3L]][i_model, , drop = TRUE])
  )
  
  # summarize results from model fit
  bind_cols(
    .pred_class = as.integer(aatd_rsc$predict(X_test)),
    .pred_Normal = aatd_rsc$predict_prob(X_test),
    .pred_Abnormal = 1 - aatd_rsc$predict_prob(X_test),
    select(aatd_test, geno_class)
  ) %>%
    # restore factor levels
    mutate(across(
      c(.pred_class, geno_class),
      ~ factor(ifelse(. == 1L, "Normal", "Abnormal"), c("Abnormal", "Normal"))
    )) ->
    aatd_fr_res
  
  # metric tables
  aatd_fr_res %>%
    metrics(truth = geno_class, estimate = .pred_class, .pred_Abnormal) %>%
    mutate(
      model = "FasterRisk",
      predictors = pred, response = resp,
      number = i_model
    ) ->
    aatd_fr_fold_met_i
  aatd_fr_fold_met <- bind_rows(aatd_fr_fold_met, aatd_fr_fold_met_i)
  
}

aatd_fr_met$metrics[[r]] <- aatd_fr_fold_met

}#LOOP

aatd_fr_met %>%
  unnest(metrics) %>%
  select(model, predictors, response, terms, id, number, everything()) ->
  aatd_fr_met

aatd_metrics <- bind_rows(aatd_metrics, aatd_fr_met)

write_rds(aatd_metrics, here::here("data/aatd-2-eval-fr.rds"))

}#LOOP
}#LOOP
