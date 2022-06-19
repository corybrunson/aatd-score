library(tidyverse)
library(lubridate)

here::here("data-raw/AAT Detection Database_deidenti.xlsx") %>%
  readxl::read_xlsx() %>%
  rename_with(snakecase::to_snake_case, everything()) %>%
  # remove extraneous variables
  # select(
  #   -date_of_birth, -genotype_simple,
  #   -patient_received_what_is_alpha_1_brochure, -physician_address_zipcode
  # ) %>%
  # format dates as dates
  mutate(across(contains("date"), as_date, format = "%m/%d/%y")) %>%
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
  glimpse()

# transform variables in preparation for predictive modeling
aatd_data %>%
  # guess ages
  mutate(age_guess = case_when(
    ! is.na(age_calculated) ~ age_calculated,
    ! is.na(date_of_birth) & ! is.na(date_received) ~
      interval(date_of_birth, date_received) / years(1),
    ! is.na(date_of_birth) ~ interval(date_of_birth, "2020-06-01") / years(1),
    TRUE ~ NA_real_
  )) %>%
  select(-date_of_birth, -date_received, -age_calculated) %>%
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
  drop_na() %>%
  print() ->
  aatd_pred

# join response variables back in (regardless of missingness)
aatd_data %>%
  select(record_id, starts_with("genotype"), starts_with("aat_")) %>%
  semi_join(aatd_pred, by = "record_id") ->
  aatd_resp

# save 
write_rds(aatd_pred, here::here("data/aatd-pred.rds"))
write_rds(aatd_resp, here::here("data/aatd-resp.rds"))
