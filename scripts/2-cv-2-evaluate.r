library(tidyverse)

# read in evaluation data
read_rds(here::here("data/aatd-eval.rds")) %>%
  mutate(across(c(model, predictors, response), fct_inorder)) %>%
  print() ->
  aatd_metrics

# check for duplicates
aatd_metrics %>%
  unnest(hyperparameters) %>%
  count(predictors, response, model, .metric, .estimator, name = "count") %>%
  count(predictors, response, model, count, name = "n")

# compare accuracy across hyperparameter settings
aatd_metrics %>%
  unnest(hyperparameters) %>%
  group_by_at(vars(-id, -.estimate)) %>%
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
aatd_metrics %>%
  unnest(hyperparameters) %>%
  group_by_at(vars(-id, -.estimate)) %>%
  summarize(mean = mean(.estimate), .groups = "drop") %>%
  filter(.metric == "roc_auc") %>%
  mutate(formula = interaction(predictors, response, sep = " -> ")) %>%
  mutate(formula = fct_rev(formula)) %>%
  ggplot(aes(x = mean, y = formula)) +
  facet_grid(model ~ .) +
  geom_boxplot() +
  labs(x = "AUROC", y = "Predictors -> Response")

# compare accuracy across hyperparameter settings: ZZ (greatest accuracy)
aatd_metrics %>%
  filter(response == "ZZ") %>%
  unnest(hyperparameters) %>%
  group_by_at(vars(-id, -.estimate)) %>%
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
aatd_metrics %>%
  filter(.metric == "accuracy" | .metric == "roc_auc") %>%
  group_by_at(vars(-id, -.estimate)) %>%
  summarize(mean = mean(.estimate), sd = sd(.estimate), .groups = "drop") %>%
  group_by(model, predictors, response, .metric) %>%
  filter(mean == max(mean)) %>%
  unnest(hyperparameters) %>%
  mutate(formula = interaction(predictors, response, sep = " -> ")) %>%
  mutate(formula = fct_rev(formula)) %>%
  ggplot(aes(x = mean, xmin = mean - 2*sd, xmax = mean + 2*sd, y = formula)) +
  facet_grid(model ~ .metric) +
  geom_pointrange()
# logistic regression achieves more robust (AUROC) performance

# track logistic regression performance across penalties
aatd_metrics %>%
  filter(model == "logistic regression" & .metric == "roc_auc") %>%
  group_by_at(vars(-id, -.estimate)) %>%
  summarize(mean = mean(.estimate), sd = sd(.estimate), .groups = "drop") %>%
  group_by(model, predictors, response, .metric) %>%
  # ungroup() %>%
  unnest(hyperparameters) %>%
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

stop("Below lies code for an older version of the results data.")

# compare accuracy
aatd_metrics %>%
  filter(.metric == "accuracy") %>%
  mutate(model = fct_rev(model)) %>%
  ggplot(aes(x = mean, xmin = mean - 2 * std_err, xmax = mean + 2 * std_err)) +
  facet_grid(. ~ response, scales = "free_x", space = "free_x") +
  geom_pointrange(
    aes(y = predictors, color = model),
    position = position_dodge(width = 2/3)
  ) +
  scale_y_discrete(limits = rev(levels(aatd_metrics$predictors))) +
  scale_color_brewer(
    type = "qual", palette = 2L, direction = -1,
    limits = rev(levels(aatd_metrics$model))
  ) +
  guides(color = guide_legend(reverse = TRUE)) +
  theme(
    legend.position = "bottom", legend.direction = "horizontal",
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) ->
  aatd_acc_eval
print(aatd_acc_eval)
ggsave(
  here::here(str_c("fig/aatd-cv-acc-eval.png")), aatd_acc_eval,
  width = 6, height = 4
)

# compare AUC
aatd_metrics %>%
  filter(.metric == "roc_auc") %>%
  mutate(predictors = fct_rev(predictors)) %>%
  ggplot(aes(x = mean, xmin = mean - 2 * std_err, xmax = mean + 2 * std_err)) +
  facet_grid(model ~ .) +
  geom_vline(xintercept = .5, linetype = "dotted") +
  geom_pointrange(
    aes(y = response, color = predictors),
    position = position_dodge(width = 2/3)
  ) +
  scale_y_discrete(limits = rev(levels(aatd_metrics$response))) +
  scale_color_brewer(
    type = "qual", palette = 2L, direction = -1,
    limits = rev(levels(aatd_metrics$predictors))
  ) +
  guides(color = guide_legend(reverse = TRUE)) +
  theme(legend.position = "bottom", legend.direction = "horizontal") ->
  aatd_roc_eval
print(aatd_roc_eval)
ggsave(
  here::here(str_c("fig/aatd-cv-roc-eval.png")), aatd_roc_eval,
  width = 6, height = 8
)
