#' This script uses several families of models to predict abnormal genotype from
#' liver and lung disease histories and other personal factors. Each combination
#' of predictors, response, and model is conducted, and for each combination of
#' predictors and response the models are compared using performance curves.
#' Fixed parameter values are used for all models. On each plot of curves is
#' superimposed the performance of COPD as a screening indication. Evaluations
#' are performed using a designated subset of the data, partitioned into
#' training and testing sets in a designated proportion.
#' 
#' See here for model families available through {parsnip}:
#' https://www.tidymodels.org/find/parsnip/

#' Setup

library(tidyverse)
library(discrim)
library(tidymodels)

# hyperparameters
# p_data <- 1/10
# p_data <- 1/6
p_data <- 1
n_train <- 2/3

# genotypes to include in analysis
genotype_incl <- read_rds(here::here("data/genotype-incl.rds"))

# model specifications as tidy selections
vars_predictors <- list(
  Dx = expr(c(starts_with("lung_"), starts_with("liver_"))),
  `Dx+age` =
    expr(c(contains("age_guess"), contains("receipt_date"),
           starts_with("lung_"), starts_with("liver_"))),
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

#' Model specifications

# logistic regression, fix parameters
logistic_reg(penalty = 1) %>%
  set_engine("glm") %>%
  set_mode("classification") ->
  aatd_lr_spec

# linear discriminant analysis, fix parameters
discrim_linear(penalty = 1) %>%
  set_engine("MASS") %>%
  set_mode("classification") ->
  aatd_lda_spec

# random forest, fix parameters
rand_forest(mtry = NULL, trees = 120L) %>%
  set_engine("randomForest") %>%
  set_mode("classification") ->
  aatd_rf_spec

# nearest neighbor, fix parameters
nearest_neighbor(neighbors = 360L, weight_func = "triangular") %>%
  set_engine("kknn") %>%
  set_mode("classification") ->
  aatd_nn_spec

# linear SVM, fix parameters
svm_linear() %>%
  #set_engine("LiblineaR") %>%
  set_engine("kernlab") %>%
  set_mode("classification") ->
  aatd_svm1_spec

# quadratic SVM model specification
svm_poly(degree = 2L) %>%
  set_engine("kernlab") %>%
  set_mode("classification") ->
  aatd_svm2_spec

# cubic SVM model specification
svm_poly(degree = 3L) %>%
  set_engine("kernlab") %>%
  set_mode("classification") ->
  aatd_svm3_spec

#' Result tables

aatd_copd_res <- if (file.exists(here::here("data/aatd-1-copd.rds"))) {
  read_rds(here::here("data/aatd-1-copd.rds"))
} else {
  tibble()
}
aatd_mod_res_count <- if (file.exists(here::here("data/aatd-1-count.rds"))) {
  read_rds(here::here("data/aatd-1-count.rds"))
} else {
  tibble()
}
aatd_mod_res_metric <- if (file.exists(here::here("data/aatd-1-metric.rds"))) {
  read_rds(here::here("data/aatd-1-metric.rds"))
} else {
  tibble()
}
aatd_mod_res_pred <- if (file.exists(here::here("data/aatd-1-pred.rds"))) {
  read_rds(here::here("data/aatd-1-pred.rds"))
} else {
  tibble()
}

ii <- if (file.exists(here::here("data/aatd-1-split-ii.rds"))) {
  read_rds(here::here("data/aatd-1-split-ii.rds"))
} else {
  c(0L, 0L)
}

for (i_pred in seq_along(vars_predictors)) {#LOOP
for (i_resp in seq_along(vars_response)) {#LOOP

if (i_pred < ii[[1L]] || (i_pred == ii[[1L]] && i_resp <= ii[[2L]])) next

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
  # choice of predictors
  select(record_id, eval(vars_predictors[[i_pred]])) ->
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

#' COPD indication

aatd_data %>%
  transmute(screen = lung_hx_copd, test = geno_class == "Abnormal") %>%
  count(screen, test, name = "count") %>%
  mutate(quadrant = case_when(
    screen & test ~ "TP",
    ! screen & ! test ~ "TN",
    screen & ! test ~ "FP",
    ! screen & test ~ "FN"
  )) %>%
  pivot_wider(id_cols = c(), names_from = quadrant, values_from = count) %>%
  transmute(
    sensitivity = TP / (TP + FN),
    specificity = TN / (TN + FP),
    precision = TP / (TP + FP),
    recall = sensitivity
  ) %>%
  mutate(predictors = pred, response = resp) ->
  aatd_copd_res_i
aatd_copd_res <- bind_rows(aatd_copd_res, aatd_copd_res_i)
# write to data file
write_rds(aatd_copd_res, here::here("data/aatd-1-copd.rds"))

#' Prepare recipes

# initial partition
aatd_split <- initial_split(aatd_data, prop = n_train, strata = geno_class)

# prepare regression recipe
recipe(training(aatd_split), geno_class ~ .) %>%
  # stop treating the ID as a predictor
  #update_role(record_id, new_role = "id variable") %>%
  step_rm(record_id, genotype, ends_with("_none")) %>%
  prep() ->
  aatd_reg_rec

# prepare numeric recipe
recipe(training(aatd_split), geno_class ~ .) %>%
  # stop treating the ID as a predictor
  #update_role(record_id, new_role = "id variable") %>%
  step_rm(record_id, genotype) %>%
  # one-hot encoding of factors
  step_dummy(all_nominal_predictors(), one_hot = TRUE) %>%
  # binary encoding of logicals
  step_mutate_at(has_type(match = "logical"), fn = ~ . * 2L - 1L) %>%
  prep() ->
  aatd_num_rec

# summarize results from model fit
fit_results <- function(fit, rec) {
  bind_cols(
    predict(fit, bake(rec, testing(aatd_split))),
    # -+- manually add any missing classes -+-
    predict(fit, bake(rec, testing(aatd_split)), type = "prob"),
    select(testing(aatd_split), geno_class)
  )
}

#' Logistic regression

# fit model
aatd_lr_spec %>%
  fit(geno_class ~ ., bake(aatd_reg_rec, NULL)) ->
  aatd_lr_fit

# evaluate model
aatd_lr_res <- fit_results(aatd_lr_fit, aatd_reg_rec)

#' Linear discriminant analysis

# fit model
aatd_lda_spec %>%
  fit(geno_class ~ ., bake(aatd_reg_rec, NULL)) ->
  aatd_lda_fit

# evaluate model
aatd_lda_res <- fit_results(aatd_lda_fit, aatd_reg_rec)

#' Random forest classification

# fit model
aatd_rf_spec %>%
  fit(geno_class ~ ., bake(aatd_num_rec, NULL)) ->
  aatd_rf_fit

# evaluate model
aatd_rf_res <- fit_results(aatd_rf_fit, aatd_num_rec)

#' Nearest neighbors classification

# fit model
aatd_nn_spec %>%
  fit(geno_class ~ ., bake(aatd_num_rec, NULL)) ->
  aatd_nn_fit

# evaluate model
aatd_nn_res <- fit_results(aatd_nn_fit, aatd_num_rec)

#' Support vector machine classification

# # fit model
# aatd_svm1_spec %>%
#   fit(geno_class ~ ., bake(aatd_num_rec, NULL)) ->
#   aatd_svm1_fit
# 
# # evaluate model
# aatd_svm1_res <- fit_results(aatd_svm1_fit, aatd_num_rec)

# # fit model
# aatd_svm2_spec %>%
#   fit(geno_class ~ ., bake(aatd_num_rec, NULL)) ->
#   aatd_svm2_fit
# 
# # evaluate model
# aatd_svm2_res <- fit_results(aatd_svm2_fit, aatd_num_rec)

# # fit model
# aatd_svm3_spec %>%
#   fit(geno_class ~ ., bake(aatd_num_rec, NULL)) ->
#   aatd_svm3_fit
# 
# # evaluate model
# aatd_svm3_res <- fit_results(aatd_svm3_fit, aatd_num_rec)

#' Comparison of models

# count tables
list(
  `logistic regression` = aatd_lr_res
  , `linear discriminant` = aatd_lda_res
  , `random forest` = aatd_rf_res
  , `nearest neighbor` = aatd_nn_res
  # , `linear svm` = aatd_svm1_res
  # , `quadratic svm` = aatd_svm2_res
  # , `cubic svm` = aatd_svm3_res
) %>%
  map(count, .pred_class, geno_class, name = "count") %>%
  enframe(name = "model", value = "count") %>%
  mutate(model = fct_inorder(model), predictors = pred, response = resp) %>%
  unnest(count) %>%
  select(response, predictors, model, everything(), count) ->
  aatd_mod_res_count_i
aatd_mod_res_count <- bind_rows(aatd_mod_res_count, aatd_mod_res_count_i)
# write to data file
write_rds(aatd_mod_res_count, here::here("data/aatd-1-count.rds"))

# metric tables
list(
  `logistic regression` = aatd_lr_res
  , `linear discriminant` = aatd_lda_res
  , `random forest` = aatd_rf_res
  , `nearest neighbor` = aatd_nn_res
  # , `linear svm` = aatd_svm1_res
  # , `quadratic svm` = aatd_svm2_res
  # , `cubic svm` = aatd_svm3_res
) %>%
  map(metrics, truth = geno_class, estimate = .pred_class, .pred_Abnormal) %>%
  enframe(name = "model", value = "metrics") %>%
  mutate(model = fct_inorder(model), predictors = pred, response = resp) %>%
  unnest(metrics) %>%
  select(response, predictors, model, everything()) ->
  aatd_mod_res_metric_i
aatd_mod_res_metric <- bind_rows(aatd_mod_res_metric, aatd_mod_res_metric_i)
# write to data file
write_rds(aatd_mod_res_metric, here::here("data/aatd-1-metric.rds"))

# prediction tables
list(
  `logistic regression` = aatd_lr_res
  , `linear discriminant` = aatd_lda_res
  , `random forest` = aatd_rf_res
  , `nearest neighbor` = aatd_nn_res
  # , `linear svm` = aatd_svm1_res
  # , `quadratic svm` = aatd_svm2_res
  # , `cubic svm` = aatd_svm3_res
) %>%
  enframe(name = "model", value = "results") %>%
  mutate(model = fct_inorder(model), predictors = pred, response = resp) %>%
  unnest(results) ->
  aatd_mod_res_pred_i
aatd_mod_res_pred <- bind_rows(aatd_mod_res_pred, aatd_mod_res_pred_i)
# write to data file
write_rds(aatd_mod_res_pred, here::here("data/aatd-1-pred.rds"))

# keep track of step
write_rds(c(i_pred, i_resp), here::here("data/aatd-1-split-ii.rds"))

}#LOOP
}#LOOP
