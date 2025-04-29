# PACKAGES ---- 
# Required packages
library(tidymodels)
library(readxl)
library(tidyverse)
library(here)
library(performance)
library(janitor)

# Read the dataset into R
ch4 <- read_xlsx(here("data", "DNA methylation data.xlsm"), sheet = 1)

# Explore the first few rows
head(ch4)

# Cleaning the data ---- 
colnames(ch4)

ch4 <- janitor::clean_names(ch4)

colnames(ch4)

# Checking for duplicates and deleting duplicates  
ch4 |> 
  duplicated() |>
  sum()

# [1] 0 (meaning there are no duplication's)

# Deleting missing data 
ch4 <- ch4 |>
  drop_na(cp_g_gria2_1,cp_g_gria2_2,aspa_1) %>% 
  select(!age_category & !sample)

# The names have been cleaned and any missing data has been deleted as it can interfere with the model. 

# Define the linear regression model ---- 
linear_model <- linear_reg() |> 
  set_engine("lm")  # We are using the "lm" engine for linear regression

# Preprocess the Data with a Recipe ---- 
age_recipe <- recipe(age ~ ., data = ch4) |> 
  step_center(all_predictors())  # Centering the predictors

# Create a workflow ----
workflow_model <- workflow() |> 
  add_model(linear_model) |> 
  add_recipe(age_recipe)

# Fit the model ----
fit_model <- fit(workflow_model, data = ch4)

# Evaluate the model ---- 
fit_model_summary <- tidy(fit_model)
fit_model_summary 

# We can visualize how well our predictions match the actual ages: 
# Get predictions on the training data

predictions_lm <- augment(fit_model, new_data = ch4)

# Plot observed vs. predicted values
ggplot(data = ch4, aes(x = age, y = predictions_lm$.pred)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Observed vs Predicted Bat Age") 

# Performance metrics ---- 
mse <- mean((ch4$age - predictions_lm$.pred)^2)

# A cleaner function for calculating MSE
mse_impl <- function(model, data, predictor) {
  augment(model, new_data = data) |> 
    mutate(squared_error = (.pred - {{predictor}})^2) |> 
    summarise(mse = mean(squared_error)) |> 
    pull(mse)
}

mse <-  mse_impl(fit_model, ch4, age)

rmse <- sqrt(mse)
rmse
# rmse = 2.087959 

# Get more comprehensive model statistics
glance(fit_model)

# Check model assumptions ---- 
fit_model |> 
  extract_fit_engine() |> 
  check_model()



