library(tidyverse)

# read in evaluation data
read_rds(here::here("data/aatd-eval.rds")) %>%
  mutate(across(c(model, predictors, response), fct_inorder)) %>%
  print() ->
  aatd_metrics

# compare accuracy
aatd_metrics %>%
  filter(.metric == "accuracy") %>%
  ggplot(aes(x = mean, xmin = mean - 2 * std_err, xmax = mean + 2 * std_err)) +
  facet_grid(. ~ response, scales = "free_x", space = "free_x") +
  geom_pointrange(
    aes(y = predictors, color = model),
    position = position_dodge(width = 2/3)
  ) +
  scale_y_discrete(limits = rev(levels(aatd_metrics$predictors))) +
  scale_color_discrete(limits = rev(levels(aatd_metrics$model))) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.text.x = element_text(angle = 30, hjust = 1)) ->
  aatd_acc_eval
print(aatd_acc_eval)
ggsave(
  here::here(str_c("fig/aatd-cv-acc-eval.png")), aatd_acc_eval,
  width = 6, height = 4
)

# compare AUC
aatd_metrics %>%
  filter(.metric == "roc_auc") %>%
  ggplot(aes(x = mean, xmin = mean - 2 * std_err, xmax = mean + 2 * std_err)) +
  facet_grid(model ~ .) +
  geom_vline(xintercept = .5, linetype = "dotted") +
  geom_pointrange(
    aes(y = response, color = predictors),
    position = position_dodge(width = 2/3)
  ) +
  scale_y_discrete(limits = rev(levels(aatd_metrics$response))) +
  scale_color_discrete(limits = rev(levels(aatd_metrics$predictors))) +
  theme(legend.position = "bottom", legend.direction = "horizontal") ->
  aatd_roc_eval
print(aatd_roc_eval)
ggsave(
  here::here(str_c("fig/aatd-cv-roc-eval.png")), aatd_roc_eval,
  width = 6, height = 8
)
