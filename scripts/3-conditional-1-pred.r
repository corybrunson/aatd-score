#' This document will calculate P( abnormal ) and P( ZZ | abnormal ) using a
#' Dx-only whole-population model and the following conditional models:
#' 1. Dx-only
#' 2. age-informed
#' 3,4. age- and smoking-informed
#' These models will be fitted to the subset of data with all of these elements.
#' They will be evaluated for their predictive performance, tested for goodness
#' of fit, and compared for their interpretability.

library(tidyverse)
library(tidymodels)

#' Setup

set.seed(78280L)

# hyperparameters
p_data <- 1/6
p_data <- 1
n_train <- 2/3

# genotypes to include in analysis
genotype_incl <- read_rds(here::here("data/genotype-incl.rds"))

# illustrative prediction
n_demo <- 4L

#' Pre-process data

# load and subset predictors data
read_rds(here::here("data/aatd-pred.rds")) %>%
  sample_frac(size = p_data) %>%
  # choice of predictors
  select(
    record_id,
    c(starts_with("lung_"), starts_with("liver_")),
    age_guess, smoking_history_cigarette, any_tobacco_exposure
  ) %>%
  filter(smoking_history_cigarette != "(Missing)" &
           any_tobacco_exposure != "(Missing)") ->
  aatd_data
# inspect rates of tobacco variables
aatd_data %>%
  count(smoking_history_cigarette, any_tobacco_exposure)

# join in responses
read_rds(here::here("data/aatd-resp.rds")) %>%
  select(record_id, genotype) %>%
  semi_join(genotype_incl, by = "genotype") %>%
  # remove missing responses
  drop_na() %>%
  # responses
  mutate(
    # abnormal
    geno_class = ifelse(
      genotype == "SZ" | genotype == "ZZ",
      "Abnormal", "Normal"
    ),
    geno_class = factor(geno_class, c("Abnormal", "Normal")),
    # ZZ specifically
    geno_abtype = ifelse(geno_class == "Abnormal", genotype, NA_character_),
    geno_abtype = factor(geno_abtype, c("ZZ", "SZ"))
  ) %>%
  inner_join(aatd_data, by = "record_id") ->
  aatd_data
# inspect responses
aatd_data %>%
  select(starts_with("geno")) %>%
  distinct() %>%
  arrange(geno_abtype, geno_class)

#' Model specifications

# logistic regression, fix parameters
logistic_reg(penalty = 1) %>%
  set_engine("glm") %>%
  set_mode("classification") ->
  aatd_lr_spec

# initial partition
aatd_split <- initial_split(aatd_data, prop = n_train, strata = geno_class)

#' Whole-population model (Dx only)

# prepare regression recipe
recipe(training(aatd_split), geno_class ~ .) %>%
  # remove extraneous variables
  step_rm(record_id, genotype, ends_with("_none")) %>%
  # classify secondary response and predictors as non-predictors
  remove_role(
    geno_abtype, age_guess, smoking_history_cigarette, any_tobacco_exposure,
    old_role = "predictor"
  ) %>%
  prep() ->
  aatd_dx_rec

# prepare workflow
workflow() %>%
  add_model(aatd_lr_spec) %>%
  add_recipe(aatd_dx_rec) %>%
  fit(training(aatd_split)) ->
  aatd_dx_fit

# sample for illustration, balanced by normality, abnormal allelle, and tobacco
testing(aatd_split) %>%
  nest(cohort = -c(geno_class, geno_abtype, any_tobacco_exposure)) %>%
  mutate(cohort = map(cohort, ~ sample_n(., size = n_demo))) %>%
  arrange(geno_class, geno_abtype, any_tobacco_exposure) %>%
  unnest(cohort) ->
  test_aatd_sample

# predict abnormal genotype
aatd_dx_fit %>%
  predict(test_aatd_sample, type = "prob") %>%
  print() ->
  aatd_dx_abnormal_pred

# subset to abnormal genotypes
training(aatd_split) %>%
  filter(geno_class == "Abnormal") ->
  train_aatd_abnormal
testing(aatd_split) %>%
  filter(geno_class == "Abnormal") ->
  test_aatd_abnormal

# predict ZZ conditional on abnormal
train_aatd_abnormal %>%
  count(geno_abtype) %>%
  summarize(
    .pred_ZZ = sum(n * (geno_abtype == "ZZ")) / sum(n),
    .pred_SZ = sum(n * (geno_abtype == "SZ")) / sum(n)
  ) %>%
  print() ->
  aatd_dx_zz_pred
# write to file
write_rds(aatd_dx_zz_pred, here::here("data/aatd-dx-zz-pred.rds"))

# summarize results from model fit
fit_results <- function(wf) {
  bind_cols(
    predict(wf, test_aatd_abnormal),
    # -+- manually add any missing classes -+-
    predict(wf, test_aatd_abnormal, type = "prob"),
    select(test_aatd_abnormal, genotype)
  )
}

#' Conditional Dx-only model

# prepare regression recipe
recipe(train_aatd_abnormal, geno_abtype ~ .) %>%
  # remove extraneous variables
  step_rm(record_id, genotype, ends_with("_none")) %>%
  # classify secondary response and predictors as non-predictors
  remove_role(
    geno_class, age_guess, smoking_history_cigarette, any_tobacco_exposure,
    old_role = "predictor"
  ) %>%
  prep() ->
  aatd_cond_dx_rec

# prepare workflow
workflow() %>%
  add_model(aatd_lr_spec) %>%
  add_recipe(aatd_cond_dx_rec) %>%
  fit(train_aatd_abnormal) ->
  aatd_cond_dx_fit

# predict abnormal genotype
aatd_cond_dx_fit %>%
  predict(test_aatd_sample, type = "prob") %>%
  print() ->
  aatd_cond_dx_pred

#' Conditional Dx/age model

# prepare regression recipe
recipe(train_aatd_abnormal, geno_abtype ~ .) %>%
  # remove extraneous variables
  step_rm(record_id, genotype, ends_with("_none")) %>%
  # classify secondary response and predictors as non-predictors
  remove_role(
    geno_class, smoking_history_cigarette, any_tobacco_exposure,
    old_role = "predictor"
  ) %>%
  prep() ->
  aatd_cond_age_rec

# prepare workflow
workflow() %>%
  add_model(aatd_lr_spec) %>%
  add_recipe(aatd_cond_age_rec) %>%
  fit(train_aatd_abnormal) ->
  aatd_cond_age_fit

# predict abnormal genotype
aatd_cond_age_fit %>%
  predict(test_aatd_sample, type = "prob") %>%
  print() ->
  aatd_cond_age_pred

#' Conditional Dx/age/smoking model

# prepare regression recipe
recipe(train_aatd_abnormal, geno_abtype ~ .) %>%
  # remove extraneous variables
  step_rm(record_id, genotype, ends_with("_none")) %>%
  # classify secondary response and predictors as non-predictors
  remove_role(
    geno_class, smoking_history_cigarette,
    old_role = "predictor"
  ) %>%
  prep() ->
  aatd_cond_cig_rec

# prepare workflow
workflow() %>%
  add_model(aatd_lr_spec) %>%
  add_recipe(aatd_cond_cig_rec) %>%
  fit(train_aatd_abnormal) ->
  aatd_cond_cig_fit

# predict abnormal genotype
aatd_cond_cig_fit %>%
  predict(test_aatd_sample, type = "prob") %>%
  print() ->
  aatd_cond_cig_pred

#' Conditional Dx/age/exposure model

# prepare regression recipe
recipe(train_aatd_abnormal, geno_abtype ~ .) %>%
  # remove extraneous variables
  step_rm(record_id, genotype, ends_with("_none")) %>%
  # classify secondary response and predictors as non-predictors
  remove_role(
    geno_class, any_tobacco_exposure,
    old_role = "predictor"
  ) %>%
  prep() ->
  aatd_cond_exp_rec

# prepare workflow
workflow() %>%
  add_model(aatd_lr_spec) %>%
  add_recipe(aatd_cond_exp_rec) %>%
  fit(train_aatd_abnormal) ->
  aatd_cond_exp_fit

# predict abnormal genotype
aatd_cond_exp_fit %>%
  predict(test_aatd_sample, type = "prob") %>%
  print() ->
  aatd_cond_exp_pred

#' Analysis

# combine predictions
bind_rows(
  mutate(fit_results(aatd_cond_dx_fit), predictors = "Dx"),
  mutate(fit_results(aatd_cond_age_fit), predictors = "Dx,age"),
  mutate(fit_results(aatd_cond_cig_fit), predictors = "Dx,age,cigarette"),
  mutate(fit_results(aatd_cond_exp_fit), predictors = "Dx,age,exposure")
) %>%
  mutate(predictors = fct_inorder(predictors)) ->
  aatd_cond_pred
# write to file
write_rds(aatd_cond_pred, here::here("data/aatd-cond-pred.rds"))

#' Illustrations

# bind predictions for analysis
test_aatd_sample %>%
  bind_cols(aatd_dx_abnormal_pred) %>%
  bind_cols(rename_with(aatd_dx_zz_pred, ~ str_replace(., "_", "_1_dx_"))) %>%
  bind_cols(rename_with(aatd_cond_dx_pred, ~ str_replace(., "_", "_2_dx_"))) %>%
  bind_cols(rename_with(aatd_cond_age_pred, ~ str_replace(., "_", "_2_age_"))) %>%
  bind_cols(rename_with(aatd_cond_cig_pred, ~ str_replace(., "_", "_2_cig_"))) %>%
  bind_cols(rename_with(aatd_cond_exp_pred, ~ str_replace(., "_", "_2_exp_"))) %>%
  print() ->
  aatd_pred

# write to CSV
aatd_pred %>%
  mutate(across(starts_with(".pred_"), round, digits = 8L)) %>%
  select(
    record_id,
    genotype, starts_with("geno_"),
    starts_with(".pred"),
    starts_with("lung_"), starts_with("liver_"),
    age_guess, any_tobacco_exposure, smoking_history_cigarette
  ) %>%
  write_csv(here::here("data/aatd-conditional-predictions.csv"))
