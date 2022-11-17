#' Compare model specifications and engines at AAT genotype prediction
#'
#' The purpose of this script is to evaluate several optimized predictive models
#' of abnormal AAT genotype. The script compares several model specifications
#' (lung and liver history alone or together with either gender or smoking
#' history) and a subset of the engines considered in a previous script
#' (logistic regression, random forest, and nearest neighbor). The same data set
#' is used for all evaluations and comparisons. Models are tuned and evaluated
#' using folded cross-validation on these data. The partitions are stratified by
#' genotype class (abnormal versus normal). Each comparison also includes the
#' current guideline of COPD as the sole indication for screening.

library(tidyverse)
library(tidymodels)

#' Setup

source(here::here("code/settings.r"))

# genotypes to include in analysis
genotype_incl <- read_rds(here::here("data/genotype-incl.rds"))

#' Split data

set.seed(804001L)
read_rds(here::here("data/aatd-pred.rds")) %>%
  # all predictors from any specification
  select(record_id, unique(unlist(sapply(vars_predictors, eval)))) %>%
  # filter missing gender
  filter(gender != "(Missing)") %>%
  # drop any cases with missing values
  drop_na() %>%
  # join in responses
  inner_join(read_rds(here::here("data/aatd-resp.rds")), by = "record_id") %>%
  # partition by genotype group
  mutate(stratum = case_when(
    genotype == "SZ" | genotype == "ZZ" | genotype == "MM" ~ genotype,
    TRUE ~ "Ab"
  )) %>%
  initial_split(prop = p_train_2, strata = stratum, breaks = 6L) ->
  aatd_split

# save training and testing IDs
training(aatd_split) %>%
  select(record_id) %>%
  write_rds(here::here("data/aatd-2-train.rds"))
testing(aatd_split) %>%
  select(record_id) %>%
  write_rds(here::here("data/aatd-2-test.rds"))
