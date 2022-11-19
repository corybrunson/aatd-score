library(tidyverse)

# read in evaluation data
read_rds(here::here("data/aatd-2-eval-ml.rds")) %>%
  mutate(across(c(model, predictors, response), fct_inorder)) %>%
  mutate(number = 1L) %>%
  unnest(hyperparameters) ->
  aatd_ml_metrics
read_rds(here::here("data/aatd-2-eval-fr.rds")) %>%
  select(-fold) %>%
  mutate(across(c(model, predictors, response), fct_inorder)) ->
  aatd_fr_metrics

# check for duplicates
aatd_ml_metrics %>%
  count(predictors, response, model, .metric, .estimator, name = "count") %>%
  count(predictors, response, model, count, name = "n")
# remind self of numbers of models obtained
aatd_fr_metrics %>%
  distinct(predictors, response, terms, number) %>%
  count(predictors, response, terms)

# compare accuracy across hyperparameter settings
aatd_ml_metrics %>%
  bind_rows(aatd_fr_metrics) %>%
  group_by_at(vars(-id, -number, -.estimate)) %>%
  summarize(mean = mean(.estimate), .groups = "drop") %>%
  filter(.metric == "accuracy") %>%
  mutate(formula = interaction(predictors, response, sep = " -> ")) %>%
  mutate(formula = fct_rev(formula)) %>%
  ggplot(aes(x = mean, y = formula)) +
  facet_grid(model ~ .) +
  geom_boxplot() +
  labs(x = "Accuracy", y = "Predictors -> Response")
# choice of response is determinative; need to compare models of same response

# compare AUROC across hyperparameter settings
aatd_ml_metrics %>%
  bind_rows(aatd_fr_metrics) %>%
  mutate(model = fct_inorder(model)) %>%
  group_by_at(vars(-id, -number, -.estimate)) %>%
  summarize(mean = mean(.estimate), .groups = "drop") %>%
  filter(.metric == "roc_auc") %>%
  mutate(formula = interaction(predictors, response, sep = " -> ")) %>%
  mutate(formula = fct_rev(formula)) %>%
  ggplot(aes(x = mean, y = formula)) +
  facet_grid(model ~ .) +
  geom_boxplot() +
  labs(x = "AUROC", y = "Predictors -> Response")

# compare accuracy across hyperparameter settings: ZZ (greatest accuracy)
aatd_ml_metrics %>%
  bind_rows(aatd_fr_metrics) %>%
  filter(response == "ZZ") %>%
  mutate(model = fct_inorder(model)) %>%
  group_by_at(vars(-id, -number, -.estimate)) %>%
  summarize(mean = mean(.estimate)) %>%
  filter(.metric == "accuracy") %>%
  mutate(predictors = fct_rev(predictors)) %>%
  ggplot(aes(x = mean, y = predictors)) +
  facet_grid(model ~ .) +
  geom_jitter(width = 0, height = 1/6, alpha = .5) +
  labs(x = "Accuracy", y = "Predictors") +
  ggtitle("Accuracy of predictive models of ZZ genotype")
# logistic regression is more reliable, but random forest achieves greater
# performance

# best parameter settings for each model and formula by two measures
aatd_ml_metrics %>%
  bind_rows(aatd_fr_metrics) %>%
  filter(.metric == "accuracy" | .metric == "roc_auc") %>%
  mutate(model = fct_inorder(model)) %>%
  group_by_at(vars(-id, -number, -.estimate)) %>%
  summarize(mean = mean(.estimate), sd = sd(.estimate), .groups = "drop") %>%
  group_by(model, predictors, response, .metric) %>%
  filter(mean == max(mean)) %>%
  mutate(formula = interaction(predictors, response, sep = " -> ")) %>%
  mutate(formula = fct_rev(formula)) %>%
  ggplot(aes(x = mean, xmin = mean - 2*sd, xmax = mean + 2*sd, y = formula)) +
  facet_grid(model ~ .metric) +
  geom_pointrange()
# logistic regression achieves more robust (AUROC) performance

# track logistic regression performance across penalties
aatd_ml_metrics %>%
  filter(model == "logistic regression" & .metric == "roc_auc") %>%
  group_by_at(vars(-id, -number, -.estimate)) %>%
  summarize(mean = mean(.estimate), sd = sd(.estimate), .groups = "drop") %>%
  group_by(model, predictors, response, .metric) %>%
  mutate(formula = interaction(predictors, response, sep = " -> ")) %>%
  mutate(formula = fct_rev(formula)) %>%
  ggplot(aes(x = penalty, y = mean, ymin = mean - 2*sd, ymax = mean + 2*sd,
             color = predictors, linetype = response, group = formula)) +
  geom_line() +
  geom_point(
    data = ~ filter(filter(., mean == max(mean)), penalty == max(penalty))
  ) +
  scale_x_continuous(trans = "log", breaks = 10 ^ seq(-9, 0))
# tobacco use obtains greatest predictive value for each response

# compare FasterRisk models
aatd_fr_metrics %>%
  filter(.metric == "roc_auc") %>%
  group_by_at(vars(-id, -number, -.estimate)) %>%
  summarize(mean = mean(.estimate), sd = sd(.estimate), .groups = "drop") %>%
  mutate(hyperparameters = str_c("terms = ", terms, ", bound = ", bound)) %>%
  mutate(hyperparameters = fct_rev(fct_inorder(hyperparameters))) %>%
  ggplot(aes(x = mean, xmin = mean - 2*sd, xmax = mean + 2*sd,
             y = hyperparameters)) +
  facet_grid(rows = vars(predictors), cols = vars(response)) +
  geom_pointrange()
