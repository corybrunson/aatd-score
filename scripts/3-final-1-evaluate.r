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
warning("Choose these hyperparameter values based on the CV results.")
penalty_opt <- 10 ^ c(-4, 0)
n_terms_opt <- c(7L, 10L)
abs_bound_opt <- c(6L, 25L, 100L)

#' Model specifications

# restrict to models of interest
vars_predictors <- vars_predictors[c("Dx", "Dx+smoke use")]
vars_response <- vars_response[c("ZZ")]

# initialize results tibble
eval_data <- tibble()

#' Evaluation measures

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
    # make abnormal genotype the second factor level
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
    # make abnormal genotype the second factor level
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
  # reorder factors for fit
  step_mutate_at(all_outcomes(), fn = ~ factor(., c("Normal", "Abnormal"))) %>%
  prep() ->
  aatd_reg_rec

# prepare binarization recipe
recipe(aatd_train, geno_class ~ .) %>%
  # stop treating the ID as a predictor
  #update_role(record_id, new_role = "id variable") %>%
  step_rm(record_id, genotype) %>%
  # remove redundant lung & liver categories
  step_rm(ends_with("_none")) %>%
  # one-hot encoding of factors
  step_dummy(all_nominal_predictors(), one_hot = FALSE) %>%
  # binary encoding of logicals
  step_mutate_at(has_type(match = "logical"), fn = as.integer) %>%
  # -1/1 encoding of response
  step_mutate_at(
    has_role(match = "outcome"),
    # (CHANGED)
    fn = ~ ifelse(. == "Abnormal", 1L, -1L)
  ) %>%
  prep() ->
  aatd_int_rec

#' Evaluate COPD model

# calculate proportions of abnormal genotypes by COPD
aatd_train %>%
  group_by(geno_class, lung_hx_copd) %>%
  count(name = "count") %>%
  ungroup() %>%
  pivot_wider(lung_hx_copd, names_from = geno_class, values_from = count) %>%
  rowwise() %>%
  transmute(
    score = as.integer(lung_hx_copd),
    .pred_Normal = Normal / (Abnormal + Normal),
    .pred_Abnormal = Abnormal / (Abnormal + Normal)
  ) %>%
  mutate(risk = .pred_Abnormal) %>%
  ungroup() ->
  aatd_copd_pred_tab

# evaluate COPD as a predictor
aatd_train %>%
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
  aatd_copd_table

# predictions for testing data
aatd_test %>%
  transmute(class = geno_class, score = as.integer(lung_hx_copd)) %>%
  left_join(aatd_copd_pred_tab, by = "score") ->
  aatd_copd_pred

aatd_copd_pred %>%
  # pull(class) %>% levels()
  # NB: specify which level is the 'event'
  roc_curve(truth = class, .pred_Abnormal, event_level = "first") ->
  aatd_copd_roc
aatd_copd_roc %>%
  ggplot(aes(x = specificity, y = sensitivity)) +
  coord_equal() +
  geom_path() +
  geom_abline(intercept = 1, slope = -1, lty = 3)
aatd_copd_roc %>%
  slice_min(abs(specificity - aatd_copd_table$specificity), n = 1L) %>%
  as_tibble() ->
  aatd_copd_cut

# metric tables
aatd_copd_pred %>%
  mutate(.pred_class = factor(
    ifelse(score == 0, "Normal", "Abnormal"),
    levels = levels(class)
  )) %>%
  metrics(truth = class, estimate = .pred_class, .pred_Abnormal) ->
  aatd_copd_met

# append to evaluation data
eval_data <- bind_rows(eval_data, tibble(
  model = "COPD", predictors = "COPD", response = resp,
  hyperparameters = list(tibble()),
  referent_risk = 0,
  point_vals = list(c(lung_hx_copd = 1)),
  score_risk = list(select(aatd_copd_pred_tab, score, risk = .pred_Abnormal)),
  predictions = list(aatd_copd_pred),
  cutoff = list(aatd_copd_cut),
  metrics = list(aatd_copd_met)
))

#' Evaluate logistic regression-based risk score

aatd_train_reg <- bake(aatd_reg_rec, new_data = aatd_train)
aatd_test_reg <- bake(aatd_reg_rec, new_data = aatd_test)

for (penalty in penalty_opt) {

# fit logistic regression model
logistic_reg(penalty = penalty, mixture = mixture_fix) %>%
  set_engine("glmnet") %>%
  set_mode("classification") %>%
  fit(geno_class ~ ., aatd_train_reg) ->
  aatd_lr_fit

# get referent risk factor profile (non-Dx, non-smoke, female)
aatd_train_reg %>%
  # remove response
  select(-geno_class) %>%
  # set reference values
  summarize(across(everything(), min)) ->
  aatd_lr_ref

# get coefficients (referent profile and distances from it)
aatd_lr_fit %>%
  tidy() %>%
  select(term, estimate) %>%
  deframe() ->
  aatd_lr_coef

# separate intercept as baseline value
aatd_lr_int <- aatd_lr_coef[names(aatd_lr_coef) == "(Intercept)"]
aatd_lr_coef <- aatd_lr_coef[names(aatd_lr_coef) != "(Intercept)"]

for (abs_bound in abs_bound_opt) {

# determine weights (divide maximum effect by maximum weight)
coef_denom <- max(abs(aatd_lr_coef)) / abs_bound
aatd_lr_pt <- round(aatd_lr_coef / coef_denom)

# calculate probability conversion
# tibble(score = seq(0, sum(aatd_lr_pt))) %>%
#   mutate(risk = 1 / (1 + exp(- aatd_lr_int - score))) ->
#   aatd_lr_tab
tibble(score = seq(sum(pmin(aatd_lr_pt, 0)), sum(pmax(aatd_lr_pt, 0)))) %>%
  mutate(score_adj = aatd_lr_int + score) %>%
  mutate(risk = 1 / (1 + exp(- score_adj))) ->
  aatd_lr_tab

# calculate scores and risk estimates
aatd_test_reg %>%
  rowwise() %>%
  mutate(score = sum(aatd_lr_pt * c_across(all_of(names(aatd_lr_pt))))) %>%
  ungroup() %>%
  mutate(risk = 1 / (1 + exp(- aatd_lr_int - score))) %>%
  select(score, risk) ->
  aatd_lr_risk

# evaluate logistic regression-based risk score
bind_cols(
  select(aatd_test_reg, class = geno_class),
  mutate(
    aatd_lr_risk,
    .pred_Normal = 1 - risk,
    .pred_Abnormal = risk
  )
) ->
  aatd_lr_pred

# cutoffs for specificity near that of COPD recommendation
aatd_lr_pred %>%
  # pull(class) %>% levels()
  # NB: specify which level is the 'event'
  roc_curve(truth = class, .pred_Abnormal, event_level = "second") ->
  aatd_lr_roc
aatd_lr_roc %>%
  ggplot(aes(x = specificity, y = sensitivity)) +
  coord_equal() +
  geom_path() +
  geom_abline(intercept = 1, slope = -1, lty = 3)
aatd_lr_roc %>%
  slice_min(abs(specificity - aatd_copd_table$specificity), n = 1L) %>%
  as_tibble() ->
  aatd_lr_cut

# metric tables
aatd_lr_pred %>%
  mutate(
    .pred_class = factor(
      ifelse(.pred_Abnormal > aatd_lr_cut$.threshold, "Abnormal", "Normal"),
      levels = levels(aatd_test_reg$geno_class)
    )
  ) %>%
  # NB: correct order of factors for `metrics()`
  mutate(across(c(class, .pred_class), fct_rev)) %>%
  metrics(truth = class, estimate = .pred_class, .pred_Abnormal) ->
  aatd_lr_met

# append to evaluation data
eval_data <- bind_rows(eval_data, tibble(
  model = "logistic regression", predictors = pred, response = resp,
  hyperparameters = list(tibble(penalty = penalty, abs_bound = abs_bound)),
  referent_risk = aatd_lr_int,
  point_vals = list(aatd_lr_pt),
  score_risk = list(aatd_lr_tab),
  predictions = list(aatd_lr_pred),
  cutoff = list(aatd_lr_cut),
  metrics = list(aatd_lr_met)
))

}

}

#' Evaluate FasterRisk risk score

aatd_train_int <- bake(aatd_int_rec, new_data = aatd_train)
aatd_test_int <- bake(aatd_int_rec, new_data = aatd_test)

# convert to `numpy` arrays
aatd_train_int %>%
  select(-geno_class) %>%
  np$asarray() %>%
  np_array(dtype = "int") ->
  X_train
aatd_train_int %>%
  pull(geno_class) %>%
  np$asarray() %>%
  np_array(dtype = "int") ->
  y_train
aatd_test_int %>%
  select(-geno_class) %>%
  np$asarray() %>%
  np_array(dtype = "int") ->
  X_test
aatd_test_int %>%
  pull(geno_class) %>%
  np$asarray() %>%
  np_array(dtype = "int") ->
  y_test

for (n_terms in n_terms_opt) {
for (abs_bound in abs_bound_opt) {

# specify model
aatd_rso <- fr$RiskScoreOptimizer(
  X = X_train, y = y_train,
  k = n_terms, select_top_m = n_models_fix,
  lb = -abs_bound, ub = abs_bound,
  parent_size = n_retain_fix, num_ray_search = n_mults_fix
)

# optimize model
aatd_rso$optimize()

# save results and destroy optimizer
# https://github.com/jiachangliu/FasterRisk/blob/main/src/fasterrisk/fasterrisk.py#L102
aatd_rso_res <- aatd_rso$get_models()
rm(aatd_rso)

for (i_model in seq(n_models_fix)) {

# build classifier
aatd_rsc <- fr$RiskScoreClassifier(
  multiplier = aatd_rso_res[[1L]][[i_model]],
  intercept = aatd_rso_res[[2L]][i_model],
  coefficients = np$asarray(aatd_rso_res[[3L]][i_model, , drop = TRUE])
)

# recover coefficients
aatd_fr_pt <- aatd_rsc$coefficients
names(aatd_fr_pt) <- names(aatd_lr_pt)

# calculate probability conversion
# https://github.com/jiachangliu/FasterRisk/blob/main/src/fasterrisk/fasterrisk.py#L147
# https://github.com/jiachangliu/FasterRisk/blob/main/src/fasterrisk/fasterrisk.py#L164
# NB: this is the integer score, without the intercept or multiplier
# tibble(score = seq(0, sum(aatd_fr_pt))) %>%
#   mutate(risk = 1 / (1 + exp(- aatd_rsc$intercept - score))) ->
#   aatd_fr_tab
tibble(score = seq(
  aatd_rsc$intercept + sum(pmin(aatd_fr_pt, 0)),
  aatd_rsc$intercept + sum(pmax(aatd_fr_pt, 0))
)) %>%
  mutate(score_adj = score / aatd_rsc$multiplier) %>%
  mutate(risk = 1 / (1 + exp(- score_adj))) ->
  aatd_fr_tab

# calculate scores and risk estimates
aatd_test_int %>%
  rowwise() %>%
  mutate(score = sum(aatd_fr_pt * c_across(all_of(names(aatd_fr_pt))))) %>%
  ungroup() %>%
  mutate(risk = 1 /
           (1 + exp(- (aatd_rsc$intercept + score) / aatd_rsc$multiplier))) %>%
  select(score, risk) ->
  aatd_fr_risk

# summarize results from model fit
aatd_fr_risk %>%
  # validate manual calculations above with internal calculations here
  bind_cols(
  # .pred_class = as.integer(aatd_rsc$predict(X_test)),
  .pred_Normal = as.vector(1 - aatd_rsc$predict_prob(X_test)),
  .pred_Abnormal = as.vector(aatd_rsc$predict_prob(X_test)),
  select(aatd_test_int, geno_class)
) %>%
  # restore factor levels
  mutate(across(
    any_of(c(".pred_class", "geno_class")),
    ~ factor(ifelse(. == 1L, "Abnormal", "Normal"), c("Normal", "Abnormal"))
  )) ->
  aatd_fr_pred

# cutoffs for specificity near that of COPD recommendation
aatd_fr_pred %>%
  # pull(geno_class) %>% levels()
  # NB: specify which level is the 'event'
  roc_curve(truth = geno_class, .pred_Abnormal, event_level = "second") ->
  aatd_fr_roc
aatd_fr_roc %>%
  ggplot(aes(x = specificity, y = sensitivity)) +
  coord_equal() +
  geom_path() +
  geom_abline(intercept = 1, slope = -1, lty = 3)
aatd_fr_roc %>%
  slice_min(abs(specificity - aatd_copd_table$specificity), n = 1L) %>%
  as_tibble() ->
  aatd_fr_cut

# metric tables
aatd_fr_pred %>%
  mutate(
    .pred_class = factor(
      ifelse(.pred_Abnormal > aatd_fr_cut$.threshold, "Abnormal", "Normal"),
      levels = levels(aatd_fr_pred$geno_class)
    )
  ) %>%
  metrics(truth = geno_class, estimate = .pred_class, .pred_Normal) ->
  aatd_fr_met

# append to evaluation data
eval_data <- bind_rows(eval_data, tibble(
  model = "FasterRisk", predictors = pred, response = resp,
  hyperparameters = list(tibble(n_terms = n_terms, abs_bound = abs_bound)),
  referent_risk = aatd_rsc$intercept,
  point_vals = list(aatd_fr_pt),
  score_risk = list(aatd_fr_tab),
  predictions = list(aatd_fr_pred),
  cutoff = list(aatd_fr_cut),
  metrics = list(aatd_fr_met)
))

}

}
}

write_rds(eval_data, here::here("data/aatd-2-eval-final.rds"))

}#LOOP
}#LOOP
