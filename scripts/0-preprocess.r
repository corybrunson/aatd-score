library(tidyverse)
library(lubridate)
library(magrittr)
library(visdat)

here::here("data-raw/AAT Detection Database_deidenti.xlsx") %>%
  readxl::read_xlsx() %>%
  rename_with(snakecase::to_snake_case, everything()) %>%
  # format dates as dates
  mutate(across(contains("date"), as_date, format = "%m/%d/%Y")) %>%
  # coerce checkbox and yes/no responses to logical
  mutate(across(
    where(~ all(. %in% c("Unchecked", "Checked"))),
    ~ . == "Checked"
  )) %>%
  mutate(across(
    where(~ all(. %in% c("No", "Yes", NA_character_))),
    ~ . == "Yes"
  )) %>%
  # coerce categorical variables to factors
  mutate(
    gender = factor(gender),
    # select the reference category
    race = fct_relevel(factor(race), "White"),
    smoking_history_cigarette = factor(
      smoking_history_cigarette,
      c("None", "Passive Smoke as a Child", "Passive Smoke as am Adult",
        "Past", "Current")
    )
  ) %>%
  rename_with(~ str_replace(., "liver_disease_history_choice", "liver_hx")) %>%
  rename_with(~ str_replace(., "lung_disease_history_choice", "lung_hx")) %>%
  mutate(
    physician_specialty = case_when(
      physician_specialty == "0" ~ NA_character_,
      str_detect(physician_specialty, "Allergist") ~ "Allergist",
      TRUE ~ physician_specialty
    ),
    physician_specialty = fct_explicit_na(physician_specialty)
  ) ->
  aatd_data

# check variables for completeness
aatd_data %>%
  mutate(across(everything(), is.na)) %>%
  summarize(across(everything(), ~ sum(.) / n())) %>%
  t() %>% as.data.frame() %>%
  set_names("Missingness") %>%
  rownames_to_column("Variable") %>%
  mutate(Missingness = str_c(round(Missingness * 100, digits = 0L), "%")) %>%
  knitr::kable()

# age distribution
ggplot(aatd_data, aes(x = age_calculated)) + geom_histogram()
# distribution of dates received
aatd_data %>%
  filter(! is.na(date_received)) %>%
  summarize(across(
    date_received,
    list(min = min, max = max, mean = mean, median = median)
  )) %>%
  glimpse()
aatd_data %>%
  ggplot(aes(x = date_received)) +
  geom_histogram()
# date received earlier than date of birth
aatd_data %>%
  mutate(row = row_number()) %>%
  filter(is.na(age_calculated)) %>%
  filter(! is.na(date_of_birth) & ! is.na(date_received)) %>%
  mutate(age_calc = interval(date_of_birth, date_received) / years(1)) %>%
  select(row, date_of_birth, date_received, age_calc) %>%
  print() %T>% { hist(.$age_calc) } %>%
  mutate(contiguous_group = cumsum(is.na(lag(row)) | row != lag(row) + 1L)) %>%
  select(row, contiguous_group) %>%
  group_by(contiguous_group) %>%
  nest(rows = row) %>%
  pull(rows) %>%
  map(~ unname(unlist(.)))
# guess missing ages
aatd_data %>%
  mutate(age_avail = case_when(
    ! is.na(age_calculated) ~ "already calculated",
    ! is.na(date_of_birth) & ! is.na(date_received) &
      interval(date_of_birth, date_received) / years(1) < -3/4 ~ "implausible",
    ! is.na(date_of_birth) & ! is.na(date_received) &
      interval(date_of_birth, date_received) / years(1) <= 0 ~ "in utero?",
    ! is.na(date_of_birth) & ! is.na(date_received) ~ "birth & receipt dates",
    ! is.na(date_of_birth) ~ "birth date only",
    TRUE ~ "nothing"
  )) %>%
  # count(age_avail)
  mutate(age_guess = case_when(
    ! is.na(age_calculated) ~ age_calculated,
    ! is.na(date_of_birth) & ! is.na(date_received) ~
      interval(date_of_birth, date_received) / years(1),
    ! is.na(date_of_birth) ~ interval(date_of_birth, "2009-01-01") / years(1),
    TRUE ~ NA_real_
  )) %>%
  select(age_avail, age_guess) %>%
  ggplot(aes(x = age_guess)) +
  facet_wrap(~ age_avail, scales = "free") +
  geom_histogram()
# plot date of receipt against row
aatd_data %>%
  transmute(row = row_number(), date_received) %>%
  ggplot(aes(x = row, y = date_received)) +
  geom_point()
aatd_data %>%
  transmute(row = row_number(), date_received) %>%
  mutate(batch = cut(row, c(0, 5e4, 2e5, Inf))) %>%
  ggplot(aes(x = batch, y = date_received)) +
  geom_boxplot()
# missingness plot
aatd_data %>%
  select(-record_id, -genotype_simple, -starts_with("aat_")) %>%
  vis_dat(warn_large_data = FALSE)

# transform variables in preparation for predictive modeling
aatd_data %>%
  # impute reception date guesses
  mutate(receipt_date = ! is.na(date_received)) %>%
  mutate(date_received_impute = if_else(
    is.na(date_received),
    date_received[sample(
      which(! is.na(date_received)),
      nrow(aatd_data), replace = TRUE
    )],
    date_received
  )) %>%
  # guess ages
  mutate(age_guess = case_when(
    ! is.na(age_calculated) ~ age_calculated,
    ! is.na(date_of_birth) & ! is.na(date_received) ~ NA_real_,
    ! is.na(date_of_birth) ~
      interval(date_of_birth, date_received_impute) / years(1),
    TRUE ~ NA_real_
  )) %>%
  # note: not enough `date_received` values for temporal validation
  select(
    -date_of_birth, -date_received, -date_received_impute
  ) %>%
  # drop response variables
  select(-starts_with("genotype"), -starts_with("aat_")) %>%
  # multiple packs per day by years
  mutate(
    pack_years_current = packs_per_day_current * for_how_many_years_current,
    pack_years_past = packs_per_day_past * for_how_many_years_past
  ) %>%
  select(-starts_with("packs_per_day"), -starts_with("for_how_many_years")) %>%
  # drop remaining non-predictors
  select(
    -patient_received_what_is_alpha_1_brochure,
    -physician_address_zipcode
  ) %>%
  # drop or explicate missing values
  select(-pack_years_current, -pack_years_past) %>%
  mutate(across(
    c(gender, race, smoking_history_cigarette, any_tobacco_exposure),
    ~ fct_explicit_na(factor(.))
  )) %>%
  # allow missing values of calculated age
  drop_na(-age_calculated) %>%
  print() ->
  aatd_pred

# age guess distribution
ggplot(aatd_pred, aes(x = age_guess)) + geom_histogram()
aatd_pred %>%
  ggplot(aes(x = age_guess)) +
  facet_grid(receipt_date ~ ., scales = "free_y") +
  geom_histogram()

# join response variables back in (regardless of missingness)
aatd_data %>%
  select(record_id, starts_with("genotype"), starts_with("aat_")) %>%
  semi_join(aatd_pred, by = "record_id") %>%
  print() ->
  aatd_resp

# save patient-level data
write_rds(aatd_pred, here::here("data/aatd-pred.rds"))
write_rds(aatd_resp, here::here("data/aatd-resp.rds"))

# genotypes to include in analysis
read_rds(here::here("data/aatd-resp.rds")) %>%
  count(genotype) %>%
  filter(n > 120L) %>%
  select(genotype) %>%
  drop_na() ->
  genotype_incl

# save genotypes to include
write_rds(genotype_incl, here::here("data/genotype-incl.rds"))
