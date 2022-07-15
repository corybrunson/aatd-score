library(tidyverse)
library(tidymodels)

#' Setup

# hyperparameters
p_data <- 1/10
p_data <- 1/6
#p_data <- 1
n_train <- 2/3

# genotypes to include in analysis
read_rds(here::here("data/aatd-resp.rds")) %>%
  count(genotype) %>%
  filter(n > 120L) %>%
  select(genotype) %>%
  drop_na() ->
  genotype_incl

# model specifications as tidy selections
vars_predictors <- list(
  Dx = expr(c(starts_with("lung_"), starts_with("liver_"))),
  `Dx+age` =
    expr(c(contains("age_guess"), starts_with("lung_"), starts_with("liver_"))),
  `Dx+gender` =
    expr(c(contains("gender"), starts_with("lung_"), starts_with("liver_"))),
  `Dx+tobacco` =
    expr(c(contains("smoking_history_cigarette"),
           contains("any_tobacco_exposure"),
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

#' 1. Single-partition, fixed-parameter

#' 1.0 Model specifications

# logistic regression, fix parameters
logistic_reg(penalty = 1) %>%
  set_engine("glm") %>%
  set_mode("classification") ->
  aatd_lr_spec

# random forest, fix parameters
rand_forest(mtry = NULL, trees = 120L) %>%
  set_engine("randomForest") %>%
  set_mode("classification") ->
  aatd_rf_spec

# linear SVM, fix parameters
svm_linear() %>%
  #set_engine("LiblineaR") %>%
  set_engine("kernlab") %>%
  set_mode("classification") ->
  aatd_svm1_spec

# # quadratic SVM model specification
# svm_poly(degree = 2L) %>%
#   set_engine("kernlab") %>%
#   set_mode("classification") ->
#   aatd_svm2_spec

# cubic SVM model specification
svm_poly(degree = 3L) %>%
  set_engine("kernlab") %>%
  set_mode("classification") ->
  aatd_svm3_spec

ii <- if (file.exists(here::here("data/aatd-curve-ii.rds"))) {
  read_rds(here::here("data/aatd-curve-ii.rds"))
} else {
  c(0L, 0L)
}

for (i_pred in seq_along(vars_predictors))
for (i_resp in seq_along(vars_response)) {#LOOP

if (i_pred < ii[[1L]] || (i_pred == ii[[1L]] && i_resp < ii[[2L]])) next

pred <- names(vars_predictors)[[i_pred]]
resp <- names(vars_response)[[i_resp]]

#' 1.1. Pre-process data

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
  mutate(geno_class = ifelse(
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

#' 1.2 COPD indication

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
  ) ->
  copd_res

#' 1.3 Prepare recipes

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

#' 1.4. Logistic regression

# fit model
aatd_lr_spec %>%
  fit(geno_class ~ ., bake(aatd_reg_rec, NULL)) ->
  aatd_lr_fit

# evaluate model
aatd_lr_res <- fit_results(aatd_lr_fit, aatd_reg_rec)
aatd_lr_res %>%
  count(.pred_class, geno_class)
aatd_lr_res %>%
  metrics(truth = geno_class, estimate = .pred_class, .pred_Abnormal)

#' 1.5. Random forest classification

# fit model
aatd_rf_spec %>%
  fit(geno_class ~ ., bake(aatd_num_rec, NULL)) ->
  aatd_rf_fit

# evaluate model
aatd_rf_res <- fit_results(aatd_rf_fit, aatd_num_rec)
aatd_rf_res %>%
  count(.pred_class, geno_class)
aatd_rf_res %>%
  metrics(truth = geno_class, estimate = .pred_class, .pred_Abnormal)

#' 1.6. Support vector machine classification

# fit model
aatd_svm1_spec %>%
  fit(geno_class ~ ., bake(aatd_num_rec, NULL)) ->
  aatd_svm1_fit

# evaluate model
aatd_svm1_res <- fit_results(aatd_svm1_fit, aatd_num_rec)
aatd_svm1_res %>%
  count(.pred_class, geno_class)
aatd_svm1_res %>%
  metrics(truth = geno_class, estimate = .pred_class, .pred_Abnormal)

# # fit model
# aatd_svm2_spec %>%
#   fit(geno_class ~ ., bake(aatd_num_rec, NULL)) ->
#   aatd_svm2_fit

# # evaluate model
# aatd_svm2_res <- fit_results(aatd_svm2_fit, aatd_num_rec)
# aatd_svm2_res %>%
#   count(.pred_class, geno_class)
# aatd_svm2_res %>%
#   metrics(truth = geno_class, estimate = .pred_class, .pred_Abnormal)

# fit model
aatd_svm3_spec %>%
  fit(geno_class ~ ., bake(aatd_num_rec, NULL)) ->
  aatd_svm3_fit

# evaluate model
aatd_svm3_res <- fit_results(aatd_svm3_fit, aatd_num_rec)
aatd_svm3_res %>%
  count(.pred_class, geno_class)
aatd_svm3_res %>%
  metrics(truth = geno_class, estimate = .pred_class, .pred_Abnormal)

#' 1.7. Comparison of models

# compare ROC curves
list(
  `logistic regression` = aatd_lr_res
  , `random forest` = aatd_rf_res
  , `linear svm` = aatd_svm1_res
  # , `quadratic svm` = aatd_svm2_res
  , `cubic svm` = aatd_svm3_res
) %>%
  enframe(name = "model", value = "results") %>%
  mutate(model = fct_inorder(model)) %>%
  unnest(results) %>%
  group_by(model) %>%
  roc_curve(truth = geno_class, estimate = .pred_Abnormal) %>%
  autoplot() +
  # account for silliness of autoplot
  geom_point(data = mutate(copd_res, specificity = 1 - specificity),
             aes(x = specificity, y = sensitivity)) +
  ggtitle(str_c(pred, "-based screen for ", resp)) ->
  aatd_roc
ggsave(
  here::here(str_c("fig/aatd-", tolower(pred), "-", tolower(resp), "-roc.png")),
  aatd_roc
)

# compare PR curves
list(
  `logistic regression` = aatd_lr_res
  , `random forest` = aatd_rf_res
  , `linear svm` = aatd_svm1_res
  # , `quadratic svm` = aatd_svm2_res
  , `cubic svm` = aatd_svm3_res
) %>%
  enframe(name = "model", value = "results") %>%
  mutate(model = fct_inorder(model)) %>%
  unnest(results) %>%
  group_by(model) %>%
  pr_curve(truth = geno_class, estimate = .pred_Abnormal) %>%
  autoplot() +
  geom_point(data = copd_res, aes(x = recall, y = precision)) +
  scale_x_continuous(trans = "log") + scale_y_continuous(trans = "log") +
  ggtitle(str_c(pred, "-based screen for ", resp)) ->
  aatd_pr
ggsave(
  here::here(str_c("fig/aatd-", tolower(pred), "-", tolower(resp), "-pr.png")),
  aatd_pr
)

write_rds(c(i_pred, i_resp), here::here("data/aatd-curve-ii.rds"))

}#LOOP
