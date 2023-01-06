#' Report results using testing data

library(tidyverse)
library(tidymodels)

# read in evaluations on testing data
eval_data <- read_rds(here::here("data/aatd-3-eval-final.rds"))

# receiver operating characteristic curves
eval_data %>%
  filter(model != "COPD") %>%
  rename(family = model) %>%
  mutate(specification = map_chr(
    hyperparameters,
    ~ str_c(names(.), unlist(.), sep = " = ", collapse = ", ")
  )) %>%
  unite(col = "model", family, specification, sep = ": ", remove = FALSE) %>%
  mutate(model = ifelse(
    str_detect(model, "FasterRisk"),
    str_c(model, (row_number() - 1L) %% 2L + 1L, sep = ", #"),
    model
  )) %>%
  mutate(model = fct_inorder(model)) %>%
  # OOPS
  mutate(predictions = map(
    predictions,
    ~ rename_with(., function(x) str_remove(x, "^geno_"))
  )) %>%
  # restrict to best-performing models
  unnest(metrics) %>% filter(.metric == "roc_auc") %>%
  # pull(.estimate) %>% hist()
  filter(.estimate > 0.55) %>%
  unnest(hyperparameters) %>%
  unnest(predictions) %>%
  group_by(predictors, response, family, abs_bound, model) %>%
  roc_curve(truth = class, estimate = .pred_Abnormal,
            event_level = "second") %>%
  ggplot(aes(x = specificity, y = sensitivity, group = model)) +
  facet_grid(cols = vars(predictors)) +
  coord_equal() +
  geom_path(aes(color = family, linetype = factor(abs_bound))) +
  geom_abline(intercept = 1, slope = -1, lty = 3) +
  labs(color = "Model family", linetype = "Abs. bound") +
  theme(legend.position = "bottom")

# compare AUROC
eval_data %>%
  mutate(model = fct_inorder(model)) %>%
  unnest(metrics) %>%
  filter(.metric == "roc_auc") ->
  eval_auroc
eval_auroc %>%
  filter(model != "COPD") %>%
  mutate(specification = map_chr(
    hyperparameters,
    ~ str_c(names(.), unlist(.), sep = " = ", collapse = ", ")
  )) %>%
  ggplot(aes(x = .estimate, y = specification, color = model)) +
  facet_grid(rows = vars(predictors)) +
  geom_vline(
    data = select(distinct(filter(eval_auroc, model == "COPD")), -predictors),
    aes(xintercept = .estimate)
  ) +
  geom_point() +
  labs(x = "AUROC", y = NULL, color = NULL)

# compare histograms of scores
eval_data %>%
  filter(model != "COPD" & predictors == "Dx") %>%
  mutate(model = fct_inorder(model)) %>%
  mutate(number = case_when(
    model == "FasterRisk" ~ as.character((row_number() - 1L) %% 2L + 1L),
    TRUE ~ NA_character_
  )) %>%
  unnest(hyperparameters) %>%
  unnest(predictions) %>%
  ggplot(aes(x = score)) +
  # facet_grid(rows = vars(model), cols = vars(abs_bound), scales = "free_x",
  #            labeller = label_both) +
  facet_wrap(facets = vars(model, abs_bound), nrow = 2L, scales = "free_x",
             labeller = labeller(.cols = label_both, .multi_line = FALSE)) +
  geom_histogram(aes(alpha = number), position = "dodge",
                 binwidth = function(x) max(1L, diff(range(x)) / 24L)) +
  scale_y_continuous(trans = "log", breaks = 10 ^ seq(0L, 5L),
                     labels = scales::comma) +
  scale_alpha_manual(values = c(.5, .9), guide = "none") +
  labs(x = "Score", y = NULL)

# compare coefficients
eval_data %>%
  filter(model != "COPD" & predictors == "Dx") %>%
  mutate(model = fct_inorder(model)) %>%
  # distinguish FR models (necessary for pivoting below)
  mutate(number = case_when(
    model == "FasterRisk" ~ as.character((row_number() - 1L) %% 2L + 1L),
    TRUE ~ NA_character_
  )) %>%
  select(-predictors, -response) %>%
  mutate(point_vals = map(
    point_vals,
    enframe, name = "item", value = "points"
  )) %>%
  unnest(point_vals) ->
  eval_points
# in a table
eval_points %>%
  mutate(specification = map_chr(
    hyperparameters,
    ~ str_c(names(.), unlist(.), sep = " = ", collapse = ", ")
  )) %>%
  unite(col = "model", model, specification, sep = ": ") %>%
  mutate(points = as.integer(points)) %>%
  pivot_wider(c(model, number), names_from = item, values_from = points) %>%
  select(-number) %>%
  rename_with(~ str_remove(., "(lung|liver)\\_hx\\_")) %>%
  knitr::kable()
# in a plot
eval_points %>%
  mutate(model = fct_inorder(ifelse(
    model == "FasterRisk",
    str_c(model, " #", number), as.character(model)
  ))) %>%
  mutate(model = fct_inorder(model)) %>%
  unnest(hyperparameters) %>%
  mutate(sign = factor(sign(points), levels = c("-1", "1"))) %>%
  ggplot(aes(x = points, y = item, fill = sign)) +
  facet_grid(rows = vars(model), cols = vars(abs_bound), scales = "free_x") +
  geom_vline(xintercept = 0) +
  geom_col() +
  labs(x = "Point value", y = NULL, fill = NULL)
