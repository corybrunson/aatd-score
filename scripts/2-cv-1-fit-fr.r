
library(tidyverse)
library(tidymodels)
library(reticulate)
np <- import("numpy")
fr <- import("fasterrisk.fasterrisk")

#' Setup

source(here::here("code/settings.r"))

# genotypes to include in analysis
genotype_incl <- read_rds(here::here("data/genotype-incl.rds"))

# model specifications as tidy selections
vars_predictors <- list(
  Dx = expr(c(starts_with("lung_"), starts_with("liver_"))),
  # `Dx+age` =
  #   expr(c(contains("age_guess"), contains("receipt_date"),
  #          starts_with("lung_"), starts_with("liver_"))),
  # `Dx+tobacco` =
  #   expr(c(contains("smoking_history_cigarette"),
  #          # contains("any_tobacco_exposure"),
  #          starts_with("lung_"), starts_with("liver_"))),
  `Dx+smoke hx` = expr(c(contains("smoking_hx"),
                         starts_with("lung_"), starts_with("liver_"))),
  `Dx+smoke use` = expr(c(contains("smoking_use"),
                          starts_with("lung_"), starts_with("liver_"))),
  `Dx+gender` =
    expr(c(contains("gender"), starts_with("lung_"), starts_with("liver_")))
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
  sample_frac(size = p_data_2) %>%
  # all predictors from any specification
  select(record_id, unique(unlist(sapply(vars_predictors, eval)))) %>%
  # drop any cases with missing values
  drop_na() %>%
  # filter missing gender
  filter(gender != "(Missing)") %>%
  # store `record_id`
  select(record_id) ->
  elig_ids

#' Evaluation measures

# measures to calculate
aatd_met <- metric_set(accuracy, roc_auc, pr_auc)

#' Model specifications

# `FasterRisk` settings
# TODO: ask how the following parameters are used
n_terms <- 7L
n_models <- 6L
abs_bound <- 5
n_retain <- 12L
n_mults <- 24L

aatd_metrics <- if (file.exists(here::here("data/aatd-2-fr-eval.rds"))) {
  read_rds(here::here("data/aatd-2-fr-eval.rds"))
} else {
  tibble()
}

for (i_pred in seq_along(vars_predictors)) {#LOOP
for (i_resp in seq_along(vars_response)) {#LOOP

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

#' Specify pre-processing recipes

# prepare binarization recipe
recipe(aatd_data, geno_class ~ .) %>%
  # stop treating the ID as a predictor
  #update_role(record_id, new_role = "id variable") %>%
  step_rm(record_id, genotype) %>%
  # one-hot encoding of factors
  step_dummy(all_nominal_predictors(), one_hot = TRUE) %>%
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
  mutate(metrics = vector(mode = "list", length = nrow(.))) ->
  aatd_fr_met

for (i_fold in seq(n_folds_2)) {#LOOP
# i_fold <- 1L

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
rso <- fr$RiskScoreOptimizer(
  X = X_train, y = y_train,
  k = n_terms, select_top_m = n_models,
  lb = -abs_bound, ub = abs_bound,
  parent_size = n_retain, num_ray_search = n_mults
)

# optimize model
rso$optimize()

# save results
aatd_rso_res <- rso$get_models()

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
      ~ factor(ifelse(. == 1L, "Normal", "Abnormal"))
    )) %>%
    mutate(across(
      c(.pred_class, geno_class),
      ~ factor(., c("Normal", "Abnormal"))
    )) ->
    aatd_fr_res
  
  # metric tables
  aatd_fr_res %>%
    metrics(truth = geno_class, estimate = .pred_class, .pred_Abnormal) %>%
    mutate(
      model = "FasterRisk", number = i_model,
      predictors = pred, response = resp
    ) ->
    aatd_fr_fold_met_i
  aatd_fr_fold_met <- bind_rows(aatd_fr_fold_met, aatd_fr_fold_met_i)
  
}

aatd_fr_met$metrics[[i_fold]] <- aatd_fr_fold_met

}#LOOP

aatd_fr_met %>%
  unnest(metrics) %>%
  write_rds(here::here("data/aatd-2-fr-eval.rds"))

}#LOOP
}#LOOP