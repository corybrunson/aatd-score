#' This document will calculate P( abnormal ) and P( ZZ | abnormal ) using the
#' following two modeling approaches:
#' 1. a single, Dx-only whole-sample model
#' 2. a Dx-only whole-sample model and an age- and smoking-informed conditional
#' model
#' These models will be evaluated qualitatively in terms of their
#' interpretability.

library(tidyverse)
library(tidymodels)

#' Setup

# hyperparameters
p_data <- 1/10
p_data <- 1/6
#p_data <- 1
n_train <- 2/3

# genotypes to include in analysis
genotype_incl <- read_rds(here::here("data/genotype-incl.rds"))

# illustrative prediction
i <- 1L

#' Pre-process data

# load and subset predictors data
read_rds(here::here("data/aatd-pred.rds")) %>%
  sample_frac(size = p_data) %>%
  # choice of predictors
  select(
    record_id,
    c(starts_with("lung_"), starts_with("liver_")),
    age_guess, any_tobacco_exposure
  ) ->
  aatd_data

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

#' Model fits

#' 1. Single whole-population model, Dx only

# prepare regression recipe
recipe(training(aatd_split), geno_class ~ .) %>%
  # remove extraneous variables
  step_rm(record_id, genotype, ends_with("_none")) %>%
  # classify secondary response and predictors as non-predictors
  remove_role(
    geno_abtype, age_guess, any_tobacco_exposure,
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

# predict abnormal genotype
predict(aatd_dx_fit, slice(testing(aatd_split), i), type = "prob")

# subset to abnormal genotypes
testing(aatd_split) %>%
  filter(geno_class == "Abnormal") ->
  test_aatd_abnormal

# predict ZZ conditional on abnormal
test_aatd_abnormal %>%
  count(geno_abtype) %>%
  summarize(
    .pred_ZZ = sum(n * (geno_abtype == "ZZ")) / sum(n),
    .pred_SZ = sum(n * (geno_abtype == "SZ")) / sum(n)
  )

#' 2. + conditional Dx-only model

# prepare regression recipe
recipe(test_aatd_abnormal, geno_abtype ~ .) %>%
  # remove extraneous variables
  step_rm(record_id, genotype, ends_with("_none")) %>%
  # classify secondary response and predictors as non-predictors
  remove_role(
    geno_class, age_guess, any_tobacco_exposure,
    old_role = "predictor"
  ) %>%
  prep() ->
  aatd_cond_rec

# prepare workflow
workflow() %>%
  add_model(aatd_lr_spec) %>%
  add_recipe(aatd_cond_rec) %>%
  fit(training(aatd_split)) ->
  aatd_cond_fit

# predict abnormal genotype
predict(aatd_cond_fit, slice(testing(aatd_split), i), type = "prob")

#' 3. + conditional Dx/age/smoking model

# prepare regression recipe
recipe(test_aatd_abnormal, geno_abtype ~ .) %>%
  # remove extraneous variables
  step_rm(record_id, genotype, ends_with("_none")) %>%
  # classify secondary response and predictors as non-predictors
  remove_role(
    geno_class,
    old_role = "predictor"
  ) %>%
  prep() ->
  aatd_aug_rec

# prepare workflow
workflow() %>%
  add_model(aatd_lr_spec) %>%
  add_recipe(aatd_aug_rec) %>%
  fit(training(aatd_split)) ->
  aatd_aug_fit

# predict abnormal genotype
predict(aatd_aug_fit, slice(testing(aatd_split), i), type = "prob")
