#' Consider predictive model engines on AAT genotype data
#'
#' The purpose of this script is to quickly compare the performance of several
#' predictive modeling engines (algorithms) and decide their suitability for a
#' more rigorous test. The script uses several predictive algorithms to predict
#' abnormal genotype from liver and lung disease histories and, optionally,
#' other personal factors. Each combination of predictors, response, and model
#' is used, and for each combination of predictors and response the models are
#' compared using performance curves. Fixed parameter values are used for all
#' models. On each plot of curves is superimposed the performance of COPD as a
#' screening indication. Evaluations are performed using a designated subset of
#' the data, partitioned into training and testing sets in a designated
#' proportion.
#'
#' See here for model families available through {parsnip}:
#' https://www.tidymodels.org/find/parsnip/

#' Setup

library(tidyverse)
library(tidymodels)
library(discrim)
library(rules)

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

#' Model specifications

# logistic regression, fix parameters
logistic_reg(penalty = 1) %>%
  set_engine("glm") %>%
  set_mode("classification") ->
  aatd_lr_mod

# linear discriminant analysis, fix parameters
discrim_linear(penalty = 1) %>%
  set_engine("MASS") %>%
  set_mode("classification") ->
  aatd_lda_mod

# decision rules, fix parameters
C5_rules(trees = 1L) %>%
  set_engine("C5.0") %>%
  set_mode("classification") ->
  aatd_dr_mod

# decision tree, fix parameters
decision_tree(tree_depth = 6L) %>%
  set_engine("rpart") %>%
  set_mode("classification") ->
  aatd_dt_mod

# random forest, fix parameters
rand_forest(mtry = NULL, trees = 120L) %>%
  set_engine("randomForest") %>%
  set_mode("classification") ->
  aatd_rf_mod

# nearest neighbor, fix parameters
nearest_neighbor(neighbors = 360L, weight_func = "triangular") %>%
  set_engine("kknn") %>%
  set_mode("classification") ->
  aatd_nn_mod

# linear SVM, fix parameters
svm_linear() %>%
  #set_engine("LiblineaR") %>%
  set_engine("kernlab") %>%
  set_mode("classification") ->
  aatd_svm1_mod

# quadratic SVM model specification
svm_poly(degree = 2L) %>%
  set_engine("kernlab") %>%
  set_mode("classification") ->
  aatd_svm2_mod

# cubic SVM model specification
svm_poly(degree = 3L) %>%
  set_engine("kernlab") %>%
  set_mode("classification") ->
  aatd_svm3_mod

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

#' Load data

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

# load predictors data
read_rds(here::here("data/aatd-pred.rds")) %>%
  # choice of predictors
  select(record_id, eval(vars_predictors[[i_pred]])) %>%
  # remove missing values
  drop_na() %>%
  # subset
  semi_join(aatd_subset, by = "record_id") ->
  aatd_data

# join in responses
read_rds(here::here("data/aatd-resp.rds")) %>%
  select(record_id, genotype) %>%
  # remove missing responses
  drop_na() %>%
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
set.seed(seed)
aatd_split <- initial_split(aatd_data, prop = n_train_1, strata = geno_class)

# prepare regression recipe
recipe(training(aatd_split), geno_class ~ .) %>%
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
  step_mutate_at(has_type(match = "logical"), fn = as.integer) %>%
  prep() ->
  aatd_reg_rec

# prepare numeric recipe
recipe(training(aatd_split), geno_class ~ .) %>%
  # stop treating the ID as a predictor
  #update_role(record_id, new_role = "id variable") %>%
  step_rm(record_id, genotype) %>%
  # drop cases with missing values
  # step_naomit() %>%
  # one-hot encoding of factors
  step_dummy(all_nominal_predictors(), one_hot = TRUE) %>%
  # binary encoding of logicals
  step_mutate_at(has_type(match = "logical"), fn = as.integer) %>%
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
aatd_lr_mod %>%
  fit(geno_class ~ ., bake(aatd_reg_rec, NULL)) ->
  aatd_lr_fit

# evaluate model
aatd_lr_res <- fit_results(aatd_lr_fit, aatd_reg_rec)

#' Linear discriminant analysis

# fit model
aatd_lda_mod %>%
  fit(geno_class ~ ., bake(aatd_reg_rec, NULL)) ->
  aatd_lda_fit

# evaluate model
aatd_lda_res <- fit_results(aatd_lda_fit, aatd_reg_rec)

#' Decision rule classification

# fit model
aatd_dr_mod %>%
  fit(geno_class ~ ., bake(aatd_num_rec, NULL)) ->
  aatd_dr_fit

# evaluate model
aatd_dr_res <- fit_results(aatd_dr_fit, aatd_num_rec)

#' Decision tree classification

# fit model
aatd_dt_mod %>%
  fit(geno_class ~ ., bake(aatd_num_rec, NULL)) ->
  aatd_dt_fit

# evaluate model
aatd_dt_res <- fit_results(aatd_dt_fit, aatd_num_rec)

#' Random forest classification

# fit model
aatd_rf_mod %>%
  fit(geno_class ~ ., bake(aatd_num_rec, NULL)) ->
  aatd_rf_fit

# evaluate model
aatd_rf_res <- fit_results(aatd_rf_fit, aatd_num_rec)

#' Nearest neighbors classification

# fit model
aatd_nn_mod %>%
  fit(geno_class ~ ., bake(aatd_num_rec, NULL)) ->
  aatd_nn_fit

# evaluate model
aatd_nn_res <- fit_results(aatd_nn_fit, aatd_num_rec)

#' Support vector machine classification

# fit model
aatd_svm1_mod %>%
  fit(geno_class ~ ., bake(aatd_num_rec, NULL)) ->
  aatd_svm1_fit

# evaluate model
aatd_svm1_res <- fit_results(aatd_svm1_fit, aatd_num_rec)

# fit model
aatd_svm2_mod %>%
  fit(geno_class ~ ., bake(aatd_num_rec, NULL)) ->
  aatd_svm2_fit

# evaluate model
aatd_svm2_res <- fit_results(aatd_svm2_fit, aatd_num_rec)

# fit model
aatd_svm3_mod %>%
  fit(geno_class ~ ., bake(aatd_num_rec, NULL)) ->
  aatd_svm3_fit

# evaluate model
aatd_svm3_res <- fit_results(aatd_svm3_fit, aatd_num_rec)

#' Comparison of models

# count tables
list(
  `logistic regression` = aatd_lr_res
  , `linear discriminant` = aatd_lda_res
  , `decision rules` = aatd_dr_res
  , `decision tree` = aatd_dt_res
  , `random forest` = aatd_rf_res
  , `nearest neighbor` = aatd_nn_res
  , `linear svm` = aatd_svm1_res
  , `quadratic svm` = aatd_svm2_res
  , `cubic svm` = aatd_svm3_res
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
  , `decision rules` = aatd_dr_res
  , `decision tree` = aatd_dt_res
  , `random forest` = aatd_rf_res
  , `nearest neighbor` = aatd_nn_res
  , `linear svm` = aatd_svm1_res
  , `quadratic svm` = aatd_svm2_res
  , `cubic svm` = aatd_svm3_res
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
  , `decision rules` = aatd_dr_res
  , `decision tree` = aatd_dt_res
  , `random forest` = aatd_rf_res
  , `nearest neighbor` = aatd_nn_res
  , `linear svm` = aatd_svm1_res
  , `quadratic svm` = aatd_svm2_res
  , `cubic svm` = aatd_svm3_res
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

file.remove(here::here("data/aatd-1-split-ii.rds"))
