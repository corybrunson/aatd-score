
#' Setup

library(tidyverse)
library(tidymodels)
library(reticulate)
fr <- import("fasterrisk.fasterrisk")
os <- import("os.path")
np <- import("numpy")
pd <- import("pandas")
tm <- import("time")

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

#' Result tables

aatd_rf_res_metric <- if (file.exists(here::here("data/aatd-1-rf-metric.rds"))) {
  read_rds(here::here("data/aatd-1-rf-metric.rds"))
} else {
  tibble()
}
aatd_rf_res_pred <- if (file.exists(here::here("data/aatd-1-rf-pred.rds"))) {
  read_rds(here::here("data/aatd-1-rf-pred.rds"))
} else {
  tibble()
}

# load data and subset to a fixed stratified sample for all experiments
set.seed(seed)
read_rds(here::here("data/aatd-pred.rds")) %>%
  inner_join(read_rds(here::here("data/aatd-resp.rds")), by = "record_id") %>%
  mutate(stratum = case_when(
    genotype == "SZ" | genotype == "ZZ" | genotype == "MM" ~ genotype,
    TRUE ~ "Ab"
  )) %>%
  group_by(stratum) %>%
  slice_sample(prop = p_data_1) %>%
  ungroup() %>%
  select(record_id) ->
  aatd_subset

# for (i_pred in seq_along(vars_predictors)) {#LOOP
# for (i_resp in seq_along(vars_response)) {#LOOP

if (i_pred < ii[[1L]] || (i_pred == ii[[1L]] && i_resp <= ii[[2L]])) next

pred <- names(vars_predictors)[[i_pred]]
resp <- names(vars_response)[[i_resp]]

print("--------------------------------")
print(str_c("Predictors: ", pred))
print(str_c("Response: ", resp))
print("--------------------------------")

#' Pre-process data

# load predictors data
read_rds(here::here("data/aatd-pred.rds")) %>%
  # choice of predictors
  select(record_id, eval(vars_predictors[[i_pred]])) %>%
  # subset
  semi_join(aatd_subset, by = "record_id") ->
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

#' Prepare recipes

# initial partition
set.seed(seed)
aatd_split <- initial_split(aatd_data, prop = n_train_1, strata = geno_class)

# TODO: copy recipe from other script

#' FasterRisk

# save as `aatd_train` and `aatd_test`
# TODO: replace below with in-script `bake()`
here::here("pkg/FasterRisk/uf_alpha1/aatd_fr_train.csv") %>%
  read_csv(col_types = str_c(rep("i", 23), collapse = "")) ->
  aatd_train
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
here::here("pkg/FasterRisk/uf_alpha1/aatd_fr_test.csv") %>%
  read_csv(col_types = str_c(rep("i", 23), collapse = "")) ->
  aatd_test
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

# settings
n_terms <- 7L
abs_bound <- 5
# TODO: ask how the following parameters are used
n_retain <- 12L
n_models <- 24L

# specify model
rso <- fr$RiskScoreOptimizer(
  X = X_train, y = y_train,
  k = n_terms, select_top_m = n_models,
  lb = -abs_bound, ub = abs_bound, parent_size = n_retain
)

# optimize model
rso$optimize()

# save results
aatd_rso_res <- rso$get_models()

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
    )) ->
    aatd_fr_res
  
  # metric tables
  aatd_fr_res %>%
    metrics(truth = geno_class, estimate = .pred_class, .pred_Abnormal) %>%
    mutate(model = "FasterRisk", predictors = pred, response = resp) %>%
    select(response, predictors, model, everything()) ->
    aatd_rf_res_metric_i
  aatd_rf_res_metric <- bind_rows(aatd_rf_res_metric, aatd_rf_res_metric_i)
  # write to data file
  write_rds(aatd_rf_res_metric, here::here("data/aatd-1-rf-metric.rds"))
  
  # prediction tables
  aatd_fr_res %>%
    mutate(model = "FasterRisk", predictors = pred, response = resp) %>%
    select(response, predictors, model, everything()) ->
    aatd_rf_res_pred_i
  aatd_rf_res_pred <- bind_rows(aatd_rf_res_pred, aatd_rf_res_pred_i)
  # write to data file
  write_rds(aatd_rf_res_pred, here::here("data/aatd-1-rf-pred.rds"))
  
}

#' Comparison of models

}#LOOP
}#LOOP
