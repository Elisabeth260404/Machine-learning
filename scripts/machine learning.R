# PACKAGES ---- 
# Packages must be installed and ran before they will work! 
library(tidymodels)
library(readxl)
library(tidyverse)
library(here)
library(performance)
library(janitor)
library(yardstick) 
library(patchwork)

#__________________________________________________________________________ ----
# I have added a line in between each section so that I can keep track easier of each section as there is a lot of coding. 

# Read the data set into R
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

#__________________________________________________________________________ ----

# Define the linear regression model ---- 
linear_model <- linear_reg() |> 
  set_engine("lm")  # We are using the "lm" engine for linear regression

#__________________________________________________________________________ ----

# Preprocess the Data with a Recipe ---- 
age_recipe <- recipe(age ~ ., data = ch4) |> 
  step_center(all_predictors())  # Centering the predictors

#__________________________________________________________________________ ----

# Create a workflow ----
workflow_model <- workflow() |> 
  add_model(linear_model) |> 
  add_recipe(age_recipe)

#__________________________________________________________________________ ----

# Fit the model ----
fit_model <- fit(workflow_model, data = ch4)

#__________________________________________________________________________ ----

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

#__________________________________________________________________________ ----

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

#__________________________________________________________________________ ---- 

# Check model assumptions ---- 
fit_model |> 
  extract_fit_engine() |> 
  check_model()

#__________________________________________________________________________ ----

# Training and testing sets ---- 
# A crucial step in any machine learning project is to evaluate the model on data it hasn’t seen during training. 
## This helps assess how well it will generalise to new samples.
# A typical split is 80% for training and 20% for testing

# Split the data into 80% training and 20% testing
set.seed(123)  # For reproducibility
split <- initial_split(ch4, prop = 0.8)

# Extract training and testing sets
train_data <- training(split)
test_data <- testing(split)
# Check the dimensions of the splits
glimpse(train_data)
glimpse(test_data)
# Setting a seed ensures reproducibility of your random split. Always set and document your seed value! 

# Fit the model on the training data
lm_fit <- fit(workflow_model, data = train_data)
# View the model summary
tidy(lm_fit)

# Let's evaluate our model's performance on the test set: 
# Calculate performance metrics for linear regression
cat("RMSE is:", sqrt(mse_impl(lm_fit, test_data, age)))

####

# Make predictions on the test set using the linear model
lm_predictions <- predict(lm_fit, new_data = test_data)

# Combine predictions with actual values (truth) into a tibble
results <- tibble(
  truth = test_data$age,  # Actual values (truth)
  estimate = lm_predictions$.pred  # Predicted values (estimate)
)

# Now use yardstick's rsq function to calculate R-squared
rsq_result <- rsq(results, truth = truth, 
                  estimate = estimate)

# Print the R-squared result
rsq_result
# rsq results = rsq standard 0.301 

#__________________________________________________________________________ ----

# Cross-Validation for Robust Evaluation (K-folds) ----
# Perform 10-fold cross-validation
#v/k same thing

folds <- vfold_cv(train_data, v = 10)

# Fit models using cross-validation
cv_results <- fit_resamples(workflow_model, 
                            resamples = folds)

# Collect and visualize metrics
cv_metrics <- collect_metrics(cv_results)

ggplot(cv_metrics, aes(x = .metric, y = mean, color = .metric)) +
  geom_boxplot() +
  labs(title = "Cross-Validation Performance")

#__________________________________________________________________________ ----

# Exploring different fold numbers ---- 
fold_numbers <- seq(2, 10, 1)

cv_results_by_fold <- map_dfr(fold_numbers, function(k_value) {
  # Create cross-validation with k folds
  k_folds <- vfold_cv(train_data, v = k_value)
  
  # Fit models using cross-validation
  cv_results <- fit_resamples(workflow_model, resamples = k_folds)
  
  # Collect metrics and add k value
  cv_metrics <- collect_metrics(cv_results)
  cv_metrics$.k <- k_value
  
  return(cv_metrics)
})

# Plot performance by number of folds
ggplot(cv_results_by_fold, aes(x = .k, y = mean, color = .metric)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ .metric, scales = "free_y") +
  labs(title = "Performance by Number of Cross-Validation Folds",
       x = "Number of Folds (k)",
       y = "Mean Performance")

#__________________________________________________________________________ ----

# Final model evaluation ---- 
# Use 5-fold cross-validation for final assessment
final_k <- vfold_cv(train_data, v = 5)

# Fit with prediction saving
final_cv_results <- fit_resamples(
  workflow_model, 
  resamples = final_k,
  control = control_resamples(save_pred = TRUE)
)

# Examine performance metrics
collect_metrics(final_cv_results)

# Get all predictions across folds for visualization
cv_predictions <- collect_predictions(final_cv_results)

# Visualize predictions vs actual values across all folds
ggplot(cv_predictions, aes(x = age, y = .pred)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "Cross-Validation Predictions vs Actual Age",
       x = "Actual Age (years)",
       y = "Predicted Age (years)")

#__________________________________________________________________________ ----

# Regularised Regression: Ridge and Lasso ----
## Lasso Regression (L1 Regularisation) ----
# Define a Lasso regression model with a specific penalty
lasso_model <- linear_reg(penalty = 0.1, mixture = 1) |> 
  set_engine("glmnet")

# Create a workflow with the Lasso model
workflow_model_lasso <- workflow() |> 
  add_model(lasso_model) |> 
  add_recipe(age_recipe)

# Fit the Lasso model
fit_lasso <- fit(workflow_model_lasso, data = ch4)

# Check model performance
mse_impl(fit_lasso, ch4, age)
# IMPORTANT NOTE - The penalty parameter controls the strength of regularization. 
# Higher values lead to more coefficients being zeroed out. 
# In practice, this parameter should be tuned using cross-validation. 

## Ridge Regression (L2 Regularization) ----
# Define a Ridge regression model
ridge_model <- linear_reg(penalty = 0.1, mixture = 0) |> 
  set_engine("glmnet")

# Create a workflow with the Ridge model
workflow_model_ridge <- workflow() |> 
  add_model(ridge_model) |> 
  add_recipe(age_recipe)

# Fit the Ridge model
fit_ridge <- fit(workflow_model_ridge, data = ch4)

# Check model performance
mse_impl(fit_ridge, ch4, age)

#__________________________________________________________________________ ----

# K-Nearest Neighbors (KNN) Regression ----
# Define a KNN model with K=1
knn_model <- nearest_neighbor(neighbors = 1) |> 
  set_engine("kknn") |> 
  set_mode("regression")

# Create a workflow with the KNN model
workflow_model_knn <- workflow() |> 
  add_model(knn_model) |> 
  add_recipe(age_recipe)

# Fit the KNN model
fit_knn <- fit(workflow_model_knn, data = ch4)

# Check model performance
mse_impl(fit_knn, ch4, age)
# WARNING! When K=1, we’re only using the single most similar sample to predict age. 
# This is highly flexible but prone to over fitting, especially with noisy methylation data. 

# Let’s explore how different values of K affect the model performance:
# Try different K values
k_values <- c(1, 5, 10, 25, 50, 100)

k_results <- map_dfr(k_values, function(k) {
  # Define KNN model with current K value
  knn_model <- nearest_neighbor(neighbors = k) |>  
    set_engine("kknn") |>  
    set_mode("regression")
  
  # Create workflow
  workflow_model_knn <- workflow() |>  
    add_model(knn_model) |>  
    add_recipe(age_recipe)
  
  # Fit model
  fit_knn <- fit(workflow_model_knn, data = ch4)
  
  # Calculate MSE
  mse_value <- mse_impl(fit_knn, ch4, age)
  
  # Return results as a data frame row
  tibble(k = k, mse = mse_value, rmse = sqrt(mse_value))
})

# Plot results
ggplot(k_results, aes(x = k, y = rmse)) +
  geom_line() +
  geom_point() +
  labs(title = "KNN Performance by Number of Neighbors (K)",
       x = "Number of Neighbors (K)",
       y = "Root Mean Squared Error (years)") +
  theme_minimal()

#__________________________________________________________________________ ----

# Train and Test Comparison of All Models ---- 
# Make sure we've created train/test split first
set.seed(123)
split <- initial_split(ch4, prop = 0.8)
train_data <- training(split)
test_data <- testing(split)

# Train models on training data
fit_lm <- fit(workflow_model, data = train_data)
fit_ridge <- fit(workflow_model_ridge, data = train_data)
fit_lasso <- fit(workflow_model_lasso, data = train_data)
fit_knn <- fit(workflow_model_knn, data = train_data)

# Generate predictions on test data
predictions <- list(
  "Linear" = augment(fit_lm, new_data = test_data)$.pred,
  "Ridge" = augment(fit_ridge, new_data = test_data)$.pred,
  "Lasso" = augment(fit_lasso, new_data = test_data)$.pred,
  "KNN" = augment(fit_knn, new_data = test_data)$.pred
)

# Create visualization of predictions vs actual values
plots_list <- map2(predictions, names(predictions), function(preds, model_name) {
  ggplot(data = test_data, aes(x = age, y = preds)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(title = model_name,
         x = "Actual Age (years)",
         y = "Predicted Age (years)") +
    theme_minimal()+
    scale_y_continuous(limits = c(0,9))
})

# Display plots in a grid
wrap_plots(plots_list)

## Now, let’s quantitatively compare the models using R²: ----
# Calculate R-squared for each model
model_metrics <- map2_dfr(predictions, names(predictions), 
                          function(preds, model_name) {
                            results <- tibble(
                              truth = test_data$age,  # Actual values
                              estimate = preds        # Predicted values
                            )
                            
                            # Calculate R-squared
                            rsq_val <- yardstick::rsq(results, truth = truth, estimate = estimate)
                            
                            # Add model name and return
                            rsq_val |>  mutate(model = model_name)
                          })

# Create a bar plot comparing R-squared values
ggplot(model_metrics, aes(x = model, y = .estimate, fill = model)) +
  geom_col() +
  labs(title = "Model Comparison: R-squared on Test Data",
       x = "Model",
       y = "R-squared") +
  theme_minimal() +
  theme(legend.position = "none")

## Let’s also compare the models using RMSE: ----
# Calculate RMSE for each model
rmse_metrics <- map2_dfr(predictions, names(predictions), 
                         function(preds, model_name) {
                           results <- tibble(
                             truth = test_data$age,
                             estimate = preds
                           )
                           
                           # Calculate RMSE
                           rmse_val <- yardstick::rmse(results, truth = truth, estimate = estimate)
                           
                           # Add model name and return
                           rmse_val |>  mutate(model = model_name)
                         })

# Create a bar plot comparing RMSE values
ggplot(rmse_metrics, aes(x = model, y = .estimate, fill = model)) +
  geom_col() +
  labs(title = "Model Comparison: RMSE on Test Data",
       x = "Model",
       y = "Root Mean Squared Error (years)") +
  theme_minimal() +
  theme(legend.position = "none")

# _________________________________________________________________________ ----

# Proper Model Tuning with Cross-Validation ----
# Define parameter grid for Lasso regression
lasso_grid <- tibble(penalty = 10^seq(-3, 0, length.out = 10))

# Create tunable Lasso model
lasso_tune <- linear_reg(penalty = tune(), mixture = 1) |> 
  set_engine("glmnet")

# Create workflow
lasso_workflow <- workflow() |> 
  add_model(lasso_tune) |> 
  add_recipe(age_recipe)

# Set up cross-validation folds
set.seed(234)
folds <- vfold_cv(train_data, v = 5)

# Tune model
lasso_results <- tune_grid(
  lasso_workflow,
  resamples = folds,
  grid = lasso_grid
)

# Visualize tuning results
autoplot(lasso_results)


# Select best penalty value
best_penalty <- select_best(lasso_results, metric = "rmse")

# Finalize workflow with best parameters
final_lasso <- finalize_workflow(lasso_workflow, best_penalty)

# Fit final model on training data
final_fit <- fit(final_lasso, data = train_data)

# Evaluate on test data
final_results <- augment(final_fit, new_data = test_data)
rmse(final_results, truth = age, estimate = .pred)

# _________________________________________________________________________ ----

# Which Features Matter Most? ----
# Extract coefficients from the Lasso model
lasso_coefs <- tidy(final_fit) |> 
  filter(term != "(Intercept)") |> 
  mutate(abs_estimate = abs(estimate)) |> 
  arrange(desc(abs_estimate))

# Plot the most important CpG sites
top_n_sites <- 10
ggplot(head(lasso_coefs, top_n_sites), 
       aes(x = reorder(term, abs_estimate), y = estimate, fill = estimate > 0)) +
  geom_col() +
  coord_flip() +
  labs(title = paste("Top", top_n_sites, "CpG Sites for Age Prediction"),
       x = "CpG Site",
       y = "Coefficient Value",
       fill = "Increases with Age") +
  theme_minimal()





