#' Evaluated selected models (risk scores) on testing data

library(tidyverse)
library(tidymodels)
library(reticulate)
np <- import("numpy")
fr <- import("fasterrisk.fasterrisk")

#' Setup

source(here::here("code/settings.r"))

# genotypes to include in analysis
genotype_incl <- read_rds(here::here("data/genotype-incl.rds"))

# fixed hyperparameters
mixture_fix <- 0
n_models_fix <- 2L
n_retain_fix <- 12L
n_mults_fix <- 24L

# optimized hyperparameters
stop("Choose these hyperparameter values based on the CV results.")
penalty_opt <- 10^(-10)
# n_terms_opt <- 6L
# abs_bound_opt <- 5

#' Evaluation measures

# measures to calculate
aatd_met <- metric_set(accuracy, roc_auc, pr_auc)

for (i_pred in seq_along(vars_predictors)) {#LOOP
for (i_resp in seq_along(vars_response)) {#LOOP

pred <- names(vars_predictors)[[i_pred]]
resp <- names(vars_response)[[i_resp]]

#' Pre-process partitioned data

# load and subset predictors data
read_rds(here::here("data/aatd-pred.rds")) %>%
  # predictors from current specification
  select(record_id, eval(vars_predictors[[i_pred]])) %>%
  # restrict to training set
  semi_join(read_rds(here::here("data/aatd-2-train.rds")), by = "record_id") %>%
  # remove empty factor levels
  mutate(across(where(is.factor), fct_drop)) ->
  aatd_train

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
  inner_join(aatd_train, by = "record_id") ->
  aatd_train

# load and subset predictors data
read_rds(here::here("data/aatd-pred.rds")) %>%
  # predictors from current specification
  select(record_id, eval(vars_predictors[[i_pred]])) %>%
  # restrict to testing set
  semi_join(read_rds(here::here("data/aatd-2-test.rds")), by = "record_id") %>%
  # remove empty factor levels
  mutate(across(where(is.factor), fct_drop)) ->
  aatd_test

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
  inner_join(aatd_test, by = "record_id") ->
  aatd_test

#' Specify pre-processing recipes

# prepare regression recipe
recipe(aatd_train, geno_class ~ .) %>%
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

# prepare binarization recipe
recipe(aatd_train, geno_class ~ .) %>%
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

#' Evaluate ML-based risk score

# fit logistic regression model
logistic_reg(penalty = penalty_opt, mixture = mixture_fix) %>%
  set_engine("glmnet") %>%
  set_mode("classification") %>%
  fit(geno_class ~ ., bake(aatd_reg_rec, new_data = aatd_train)) ->
  aatd_lr_fit

# get referent risk factor profile (non-Dx, non-smoke, female)
bake(aatd_reg_rec, new_data = aatd_train) %>%
  # remove response
  select(-geno_class) %>%
  # set reference values
  summarize(across(everything(), min)) ->
  aatd_ref

# get coefficients (referent profile and distances from it)
aatd_lr_fit %>%
  tidy() %>%
  select(term, estimate) %>%
  deframe() ->
  aatd_coef

# separate intercept as baseline value
aatd_int <- aatd_coef[names(aatd_coef) == "(Intercept)"]
aatd_coef <- aatd_coef[names(aatd_coef) != "(Intercept)"]

# determine weights (divide maximum effect by maximum weight)
wt_max <- 25
coef_min <- max(abs(aatd_coef)) / wt_max
# aatd_wt <- round(aatd_coef[abs(aatd_coef) > coef_min / 2] / coef_min)
aatd_wt <- round(aatd_coef / coef_min)

# calculate risk estimates
aatd_risk_fun <- function(x) {
  stopifnot(names(x) == names(aatd_wt))
  aatd_int + sum(aatd_wt * unlist(x))
}
aatd_test %>%
  rowwise() %>%
  mutate(score = aatd_int + sum(aatd_wt * c_across(all_of(names(aatd_wt))))) %>%
  ungroup() %>%
  mutate(risk = 1 / (1 + exp(score))) %>%
  select(score, risk) ->
  aatd_test_risk

# evaluate logistic regression-based risk score
bind_cols(
  select(aatd_test, class = geno_class),
  transmute(
    aatd_test_risk,
    .pred_class = factor(ifelse(risk > .5, "Abnormal", "Normal")),
    .pred_Abnormal = risk, .pred_Normal = 1 - risk
  )
) %>%
  metrics(truth = class, estimate = .pred_class, .pred_Abnormal) ->
  aatd_lr_risk_eval

#' Evaluate FR risk score

# convert to `numpy` arrays
bake(aatd_int_rec, new_data = aatd_train) %>%
  select(-geno_class) %>%
  np$asarray() %>%
  np_array(dtype = "int") ->
  X_train
bake(aatd_int_rec, new_data = aatd_train) %>%
  pull(geno_class) %>%
  np$asarray() %>%
  np_array(dtype = "int") ->
  y_train
bake(aatd_int_rec, new_data = aatd_test) %>%
  select(-geno_class) %>%
  np$asarray() %>%
  np_array(dtype = "int") ->
  X_test
bake(aatd_int_rec, new_data = aatd_test) %>%
  pull(geno_class) %>%
  np$asarray() %>%
  np_array(dtype = "int") ->
  y_test

# specify model
aatd_rso <- fr$RiskScoreOptimizer(
  X = X_train, y = y_train,
  k = n_terms_opt, select_top_m = n_models_fix,
  lb = -abs_bound_opt, ub = abs_bound_opt,
  parent_size = n_retain_fix, num_ray_search = n_mults_fix
)

# optimize model
aatd_rso$optimize()

# save results and destroy optimizer
aatd_rso_res <- aatd_rso$get_models()
rm(aatd_rso)



stop("Extract risk score!")



}#LOOP
}#LOOP
