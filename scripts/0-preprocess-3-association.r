#' This document will assess associations among predictors based on the training
#' data.

#' Setup

library(tidyverse)
library(latentcor)
library(ordr)
library(jSDM)
library(sjSDM)

source(here::here("code/settings.r"))

# genotypes to include in analysis
genotype_incl <- read_rds(here::here("data/genotype-incl.rds"))

# load training data
read_rds(here::here("data/aatd-pred.rds")) %>%
  semi_join(read_rds(here::here("data/aatd-2-train.rds")), by = "record_id") %>%
  inner_join(read_rds(here::here("data/aatd-resp.rds")), by = "record_id") %>%
  mutate(stratum = case_when(
    genotype == "SZ" | genotype == "ZZ" | genotype == "MM" ~ genotype,
    TRUE ~ "Ab"
  )) %>%
  select(record_id) ->
  aatd_subset

# load predictors data
read_rds(here::here("data/aatd-pred.rds")) %>%
  # subset
  semi_join(aatd_subset, by = "record_id") ->
  aatd_data

#' Correlations

# latent correlations by variable type
aatd_data %>%
  select(
    gender,
    starts_with("lung_"), starts_with("liver_"),
    smoking_history_cigarette, any_tobacco_exposure
  ) ->
  aatd_corr
aatd_corr %>%
  map_chr(~ if (is.factor(.)) "ter" else if (is.logical(.)) "bin" else "con") ->
  aatd_type
aatd_corr %>%
  mutate(across(where(is.factor), as.integer)) %>%
  mutate(across(where(is.logical), ~ as.integer(.) + 1L)) %>%
  drop_na() %>%
  as.matrix() %>%
  latentcor(types = aatd_type, showplot = TRUE) ->
  aatd_latcorr
aatd_latcorr$Rpointwise %>%
  as.data.frame() %>%
  set_names(pred_names[match(names(.), names(pred_names))]) %>%
  `rownames<-`(pred_names[match(rownames(.), names(pred_names))]) %>%
  mutate(across(everything(), ~ format(round(., digits = 3L), nsmall = 3L))) %>%
  rename_with(abbreviate, minlength = 4L) %>%
  knitr::kable()
aatd_latcorr$plotR

# correlation monoplot
aatd_latcorr$Rpointwise %>%
  as.data.frame() %>%
  set_names(pred_names[match(names(.), names(pred_names))]) %>%
  `rownames<-`(pred_names[match(rownames(.), names(pred_names))]) %>%
  as.matrix() %>%
  eigen_ord() %>%
  as_tbl_ord() %>%
  augment_ord() %>%
  confer_inertia(1) ->
  aatd_lateigen
aatd_lateigen %>%
  ggbiplot() +
  theme_void() +
  geom_unit_circle(segments = 120L) +
  geom_rows_vector() +
  geom_rows_text_radiate(aes(label = name)) +
  lims(x = c(-1, 1), y = c(-1, 1))

#' Joint interaction-distribution model

aatd_corr %>%
  select(contains("_hx")) %>%
  mutate(across(everything(), as.integer)) ->
  aatd_dx
aatd_corr %>%
  select(gender, smoking_history_cigarette) %>%
  fastDummies::dummy_cols(
    c("gender", "smoking_history_cigarette"),
    remove_first_dummy = TRUE,
    remove_selected_columns = TRUE
  ) ->
  aatd_cov

# https://cran.r-project.org/web/packages/jSDM/vignettes/jSDM.html
aatd_jidm <- jSDM_binomial_probit(
  burnin = 1000L, mcmc = 1000L, thin = 1,
  presence_data = aatd_dx, site_formula = ~ ., site_data = aatd_cov,
  n_latent = 2L, site_effect = "fixed"
)
plot_residual_cor(aatd_jidm)
aatd_jidmcor <- get_residual_cor(aatd_jidm)
aatd_jidmcor$cov.median %>%
  as.data.frame() %>%
  set_names(pred_names[match(names(.), names(pred_names))]) %>%
  `rownames<-`(pred_names[match(rownames(.), names(pred_names))]) %>%
  as.matrix() %>%
  eigen_ord() %>%
  as_tbl_ord() %>%
  augment_ord() %>%
  confer_inertia(1) ->
  aatd_jidmcor_eigen
aatd_jidmcor_eigen %>%
  ggbiplot() +
  theme_void() +
  geom_unit_circle(segments = 120L) +
  geom_rows_vector() +
  geom_rows_text_radiate(aes(label = name)) +
  lims(x = c(-1, 1), y = c(-1, 1))

# https://theoreticalecology.github.io/s-jSDM/
# aatd_sjidm <- sjSDM(
#   Y = aatd_dx, env = linear(data = aatd_cov, formula = ~ .),
#   family = binomial("probit"),
#   sampling = 100L
# )
