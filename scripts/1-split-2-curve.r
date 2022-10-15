
#' Setup

library(tidyverse)
library(tidymodels)

# read in model summary data
aatd_copd_res <- read_rds(here::here("data/aatd-1-copd.rds"))
aatd_mod_res_count <- read_rds(here::here("data/aatd-1-count.rds"))
aatd_mod_res_metric <- read_rds(here::here("data/aatd-1-metric.rds"))
aatd_mod_res_pred <- read_rds(here::here("data/aatd-1-pred.rds"))

pred <- "Dx"
resp <- "ZZ"

# compare ROC curves
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
ggsave(
  here::here(str_c("fig/aatd-", tolower(pred), "-", tolower(resp), "-roc.pdf")),
  aatd_roc
)

# compare PR curves
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
ggsave(
  here::here(str_c("fig/aatd-", tolower(pred), "-", tolower(resp), "-pr.pdf")),
  aatd_pr
)
