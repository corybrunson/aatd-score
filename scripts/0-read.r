library(tidyverse)

here::here("data-raw/AAT Detection Database_deidenti.xlsx") %>%
  readxl::read_xlsx() %>%
  rename_with(snakecase::to_snake_case, everything()) %>%
  # format dates as dates
  mutate(across(contains("date"), lubridate::as_date, format = "%m/%d/%y")) %>%
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
  # -+- need to pre-process physician specialty -+-
  print() ->
  aat_data

write_rds(aat_data, here::here("data/aat-data.rds"))
