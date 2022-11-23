#' Evaluate integer programming optimizer(s)
#'
#' This script exports CSV files formatted for FasterRisk, currently the only
#' integer programming optimizer for which code is available and working. The
#' steps should match those of '2-cv-1-fit.r'.

library(tidyverse)
library(tidymodels)

#' Setup

# hyperparameters
# p_data <- 1/10
# p_data <- 1/6
p_data <- 1
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
  # `Dx+tobacco` =
  #   expr(c(contains("smoking_history_cigarette"),
  #          # contains("any_tobacco_exposure"),
  #          starts_with("lung_"), starts_with("liver_"))),
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

set.seed(385635L)
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

for (i_pred in seq_along(vars_predictors)) {###
for (i_resp in seq_along(vars_response)) {###
# TODO: for loop
# i_pred <- 2L
# i_resp <- 1L

pred <- names(vars_predictors)[[i_pred]]
resp <- names(vars_response)[[i_resp]]

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

#' Specify pre-processing recipe

# prepare binarization recipe
recipe(aatd_data, geno_class ~ .) %>%
  # stop treating the ID as a predictor
  #update_role(record_id, new_role = "id variable") %>%
  step_rm(record_id, genotype) %>%
  # one-hot encoding of factors
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
set.seed(385635L)
aatd_cv <- vfold_cv(aatd_data, v = n_folds, strata = geno_class)

#' Export data

for (i_fold in seq_along(aatd_cv$id)) {###
# TODO: for loop
# i_fold <- 1L

# pre-process training and testing data
aatd_int_rec %>%
  bake(new_data = training(aatd_cv$splits[[i_fold]])) %>%
  select(geno_class, everything()) %>%
  rename_with(~ str_replace_all(., "\\.", "\\_")) %>%
  write_csv(here::here(str_c(
    "data/aatd-pred-", i_pred, "resp-", i_resp, "-train-", i_fold, ".csv"
  )))
aatd_int_rec %>%
  bake(new_data = testing(aatd_cv$splits[[i_fold]])) %>%
  select(geno_class, everything()) %>%
  rename_with(~ str_replace_all(., "\\.", "\\_")) %>%
  write_csv(here::here(str_c(
    "data/aatd-pred-", i_pred, "resp-", i_resp, "-test-", i_fold, ".csv"
  )))

}###
}###

}###
