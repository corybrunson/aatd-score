
#' Setup

library(tidyverse)
library(tidymodels)
sample_denom <- "6"

# read in model summary data
read_rds(here::here(str_c("data/aatd-1-copd-", sample_denom, ".rds"))) %>%
  filter(predictors != "Dx+age") %>%
  mutate(across(c(predictors, response), fct_inorder)) ->
  aatd_copd_res
read_rds(here::here(str_c("data/aatd-1-count-", sample_denom, ".rds"))) %>%
  filter(predictors != "Dx+age") %>%
  mutate(across(c(predictors, response), fct_inorder)) ->
  aatd_mod_res_count
read_rds(here::here(str_c("data/aatd-1-metric-", sample_denom, ".rds"))) %>%
  filter(predictors != "Dx+age") %>%
  mutate(across(c(predictors, response), fct_inorder)) ->
  aatd_mod_res_metric
read_rds(here::here(str_c("data/aatd-1-pred-", sample_denom, ".rds"))) %>%
  filter(predictors != "Dx+age") %>%
  mutate(across(c(predictors, response), fct_inorder)) ->
  aatd_mod_res_pred
read_rds(here::here(str_c("data/aatd-1-fr-pred-", sample_denom, ".rds"))) %>%
  filter(predictors != "Dx+age") %>%
  mutate(across(c(predictors, response), fct_inorder)) ->
  aatd_fr_res_pred

pred <- "Dx"
resp <- "ZZ"

# compare ROC curves of ML models
aatd_mod_res_pred %>%
  filter(response == resp, predictors == pred) %>%
  select(-response, -predictors) %>%
  group_by(model) %>%
  roc_curve(truth = geno_class, estimate = .pred_Abnormal) %>%
  ggplot(aes(x = specificity, y = sensitivity)) +
  coord_equal() +
  geom_path(aes(color = model)) +
  geom_abline(intercept = 1, slope = -1, lty = 3) +
  geom_point(
    data = filter(aatd_copd_res, response == resp, predictors == pred)
  ) +
  ggtitle(str_c(pred, "-based screen for ", resp)) ->
  aatd_roc
print(aatd_roc)
ggsave(
  here::here(str_c("fig/aatd-", tolower(pred), "-", tolower(resp), "-roc.pdf")),
  aatd_roc
)

# compare ROC curves of FR models
aatd_fr_res_pred %>%
  filter(response == resp, predictors == pred) %>%
  select(-response, -predictors) %>%
  group_by(model) %>%
  roc_curve(truth = geno_class, estimate = .pred_Abnormal) %>%
  ggplot(aes(x = specificity, y = sensitivity)) +
  coord_equal() +
  geom_path(aes(color = model)) +
  geom_abline(intercept = 1, slope = -1, lty = 3) +
  geom_point(
    data = filter(aatd_copd_res, response == resp, predictors == pred)
  ) +
  ggtitle(str_c(pred, "-based screen for ", resp)) ->
  aatd_roc
print(aatd_roc)
ggsave(
  here::here(str_c(
    "fig/aatd-", tolower(pred), "-", tolower(resp), "-fr-roc.pdf"
  )),
  aatd_roc
)

# compare PR curves of ML models
aatd_mod_res_pred %>%
  filter(response == resp, predictors == pred) %>%
  select(-response, -predictors) %>%
  group_by(model) %>%
  pr_curve(truth = geno_class, estimate = .pred_Abnormal) %>%
  ggplot(aes(x = recall, y = precision)) +
  coord_equal() +
  geom_path(aes(color = model)) +
  geom_point(
    data = filter(aatd_copd_res, response == resp, predictors == pred)
  ) +
  # lims(x = c(0, 1), y = c(0, 1)) +
  scale_x_continuous(trans = "log", breaks = breaks_log(base = 10)) +
  scale_y_continuous(trans = "log", breaks = breaks_log(base = 10)) +
  ggtitle(str_c(pred, "-based screen for ", resp)) ->
  aatd_pr
print(aatd_pr)
ggsave(
  here::here(str_c("fig/aatd-", tolower(pred), "-", tolower(resp), "-pr.pdf")),
  aatd_pr
)

# compare PR curves of FR models
aatd_fr_res_pred %>%
  filter(response == resp, predictors == pred) %>%
  select(-response, -predictors) %>%
  group_by(model) %>%
  pr_curve(truth = geno_class, estimate = .pred_Abnormal) %>%
  ggplot(aes(x = recall, y = precision)) +
  coord_equal() +
  geom_path(aes(color = model)) +
  geom_point(
    data = filter(aatd_copd_res, response == resp, predictors == pred)
  ) +
  # lims(x = c(0, 1), y = c(0, 1)) +
  scale_x_continuous(trans = "log", breaks = breaks_log(base = 10)) +
  scale_y_continuous(trans = "log", breaks = breaks_log(base = 10)) +
  ggtitle(str_c(pred, "-based screen for ", resp)) ->
  aatd_pr
print(aatd_pr)
ggsave(
  here::here(str_c(
    "fig/aatd-", tolower(pred), "-", tolower(resp), "-fr-pr.pdf"
  )),
  aatd_pr
)

# select best ML and best FR models
stop("This step must be manually curated based on the above ROC & PR curves.")
best_models <- c(
  "logistic regression", "nearest neighbor",
  "FasterRisk 5", "FasterRisk 6"
)

# compare ROC curves of best ML and best FR models
aatd_mod_res_pred %>%
  bind_rows(aatd_fr_res_pred) %>%
  filter(model %in% best_models) %>%
  mutate(model = fct_inorder(model)) %>%
  filter(response == resp, predictors == pred) %>%
  select(-response, -predictors) %>%
  group_by(model) %>%
  roc_curve(truth = geno_class, estimate = .pred_Abnormal) %>%
  ggplot(aes(x = specificity, y = sensitivity)) +
  coord_equal() +
  geom_path(aes(color = model)) +
  geom_abline(intercept = 1, slope = -1, lty = 3) +
  geom_point(
    data = filter(aatd_copd_res, response == resp, predictors == pred)
  ) +
  ggtitle(str_c(pred, "-based screen for ", resp)) ->
  aatd_roc
print(aatd_roc)
ggsave(
  here::here(str_c(
    "fig/aatd-", tolower(pred), "-", tolower(resp), "-best-roc.pdf"
  )),
  aatd_roc
)

# compare PR curves of best ML and best FR models
aatd_mod_res_pred %>%
  bind_rows(aatd_fr_res_pred) %>%
  filter(model %in% best_models) %>%
  mutate(model = fct_inorder(model)) %>%
  filter(response == resp, predictors == pred) %>%
  select(-response, -predictors) %>%
  group_by(model) %>%
  pr_curve(truth = geno_class, estimate = .pred_Abnormal) %>%
  ggplot(aes(x = recall, y = precision)) +
  coord_equal() +
  geom_path(aes(color = model)) +
  geom_point(
    data = filter(aatd_copd_res, response == resp, predictors == pred)
  ) +
  # lims(x = c(0, 1), y = c(0, 1)) +
  scale_x_continuous(trans = "log", breaks = breaks_log(base = 10)) +
  scale_y_continuous(trans = "log", breaks = breaks_log(base = 10)) +
  ggtitle(str_c(pred, "-based screen for ", resp)) ->
  aatd_pr
print(aatd_pr)
ggsave(
  here::here(str_c(
    "fig/aatd-", tolower(pred), "-", tolower(resp), "-best-pr.pdf"
  )),
  aatd_pr
)
