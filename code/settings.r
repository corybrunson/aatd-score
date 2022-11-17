
# standard seed
seed <- 38345L

# hyperparameters for first phase

# p_data_1 <- 1/10
p_data_1 <- 1/6
# p_data_1 <- 1
n_train_1 <- 2/3

# hyperparameters for second phase

# p_train_2 <- 1/2
p_train_2 <- 4/5
# n_folds_2 <- 3L
n_folds_2 <- 6L

# model specifications as tidy selections
vars_predictors <- list(
  Dx = expr(c(starts_with("lung_"), starts_with("liver_"))),
  # `Dx+age` =
  #   expr(c(contains("age_guess"), contains("receipt_date"),
  #          starts_with("lung_"), starts_with("liver_"))),
  # `Dx+tobacco` =
  #   expr(c(contains("smoking_history_cigarette"),
  #          # contains("any_tobacco_exposure"),
  #          starts_with("lung_"), starts_with("liver_"))),
  `Dx+smoke hx` = expr(c(contains("smoking_hx"),
                         starts_with("lung_"), starts_with("liver_"))),
  `Dx+smoke use` = expr(c(contains("smoking_use"),
                          starts_with("lung_"), starts_with("liver_"))),
  `Dx+gender` =
    expr(c(contains("gender"), starts_with("lung_"), starts_with("liver_")))
)

# response variables as logical tests
vars_response <- list(
  ZZ = expr(genotype == "ZZ"),
  SZ = expr(genotype == "SZ" | genotype == "ZZ"),
  # MZ = expr(genotype == "MZ" | genotype == "SZ" | genotype == "ZZ"),
  # Z = expr(grepl(".*Z$", genotype)),
  Ab = expr(genotype != "MM")
)
