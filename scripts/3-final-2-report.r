#' Report results using testing data

library(tidyverse)

# read in evaluations on testing data
eval_data <- read_rds(here::here("data/aatd-2-eval-final.rds"))

# bar plot of AUROC
eval_data %>%
  filter(model != "COPD") %>%
  mutate(model = fct_inorder(model)) %>%
  unnest(metrics) %>%
  filter(.metric == "roc_auc") %>%
  mutate(specification = map_chr(
    hyperparameters,
    ~ str_c(names(.), unlist(.), sep = " = ", collapse = ", ")
  )) %>%
  ggplot(aes(x = .estimate, y = specification, color = model)) +
  facet_grid(rows = vars(predictors)) +
  scale_color_manual(values = c("#e41a1c", "#377eb8", "darkgrey")) +
  geom_point() +
  labs(x = "AUROC", y = NULL, color = NULL)
