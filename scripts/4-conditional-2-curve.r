
library(tidyverse)
library(tidymodels)

aatd_dx_zz_pred <- read_rds(here::here("data/aatd-dx-zz-pred.rds"))
aatd_cond_pred <- read_rds(here::here("data/aatd-cond-pred.rds"))

aatd_cond_pred %>%
  filter(genotype == "SZ" | genotype == "ZZ") %>%
  mutate(genotype = factor(genotype)) %>%
  group_by(predictors) %>%
  roc_curve(truth = genotype, estimate = .pred_SZ) %>%
  autoplot() ->
  aatd_cond_roc
ggsave(
  here::here(str_c("fig/aatd-cond-roc.pdf")),
  aatd_cond_roc
)
