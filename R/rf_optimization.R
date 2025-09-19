# Random Forest Optimization Functions
# ============================================================================

#' Simple cLHS sampling
#'
#' @param covariates SpatRaster object containing environmental covariates
#' @param n_samples Number of samples to generate
#' @param seed Random seed for reproducibility (default: NULL)
#'
#' @return Data frame with sample locations and covariate values
#' @export
simple_clhs_sampling <- function(covariates, n_samples, seed = NULL) {
  cat("Starting cLHS sampling...\n")
  
  # Set seed for cLHS reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
    cat("cLHS seed set to", seed, "\n")
  }
  
  # Convert to data.frame
  cov_df <- as.data.frame(covariates, xy = TRUE, na.rm = TRUE)
  cov_df <- na.omit(cov_df)
  
  # Prepare data for clhs (exclude x, y columns)
  clhs_data <- cov_df[, -c(1, 2), drop = FALSE]
  
  # Clean and ensure numeric data
  names(clhs_data) <- make.names(names(clhs_data))
  clhs_data <- as.data.frame(lapply(clhs_data, function(x) as.numeric(as.character(x))))
  
  # Remove constant variables
  valid_cols <- sapply(clhs_data, function(x) {
    if (all(is.na(x))) return(FALSE)
    var_val <- var(x, na.rm = TRUE)
    return(!is.na(var_val) && var_val > 1e-10)
  })
  
  clhs_data <- clhs_data[, valid_cols, drop = FALSE]
  
  cat("Using", ncol(clhs_data), "covariates for cLHS\n")
  
  # Perform cLHS
  sample_indices <- clhs::clhs(clhs_data, size = n_samples, progress = FALSE)
  
  return(cov_df[sample_indices, ])
}

#' Calculate Random Forest MSE for sampling design evaluation
#'
#' @param sample_points Data frame containing sample points with coordinates and covariates
#' @param full_data Data frame containing full population data
#' @param n_folds Number of cross-validation folds (default: 5)
#'
#' @return Numeric value representing mean squared error
#' @keywords internal
calculate_rf_mse <- function(sample_points, full_data, n_folds = 5) {
  # Prepare data (exclude coordinates for RF)
  sample_data <- sample_points[, -c(1,2)]  # Remove x, y coordinates
  
  # Use the third environmental variable as target (or first if not enough)
  if (ncol(sample_data) >= 3) {
    target_var <- names(sample_data)[3]
    predictors <- names(sample_data)[-3]
  } else {
    target_var <- names(sample_data)[1]
    predictors <- names(sample_data)[-1]
  }
  
  # Cross-validation to estimate MSE
  n_samples <- nrow(sample_data)
  fold_size <- ceiling(n_samples / n_folds)
  mse_values <- numeric(n_folds)
  
  for (i in 1:n_folds) {
    test_indices <- ((i-1) * fold_size + 1):min(i * fold_size, n_samples)
    train_data <- sample_data[-test_indices, ]
    test_data <- sample_data[test_indices, ]
    
    if (nrow(test_data) == 0) next
    
    # Train Random Forest
    rf_model <- randomForest::randomForest(
      as.formula(paste(target_var, "~", paste(predictors, collapse = "+"))),
      data = train_data, ntree = 100
    )
    
    # Predict and calculate MSE
    predictions <- predict(rf_model, test_data)
    mse_values[i] <- mean((test_data[[target_var]] - predictions)^2, na.rm = TRUE)
  }
  
  return(mean(mse_values, na.rm = TRUE))
}

#' Random Forest optimization with simulated annealing
#'
#' @param covariates SpatRaster object containing environmental covariates
#' @param n_samples Number of samples to optimize
#' @param n_iterations Number of optimization iterations (default: 500)
#' @param seed Random seed for reproducibility (default: NULL)
#'
#' @return List containing optimization results and sample locations
#' @export
rf_optimize_sampling <- function(covariates, n_samples, n_iterations = 500, seed = NULL) {
  cat("Starting RF optimization...\n")
  
  # Set seed for RF optimization reproducibility
  if (!is.null(seed)) {
    set.seed(seed + 1)  # Slightly different seed for RF
    cat("RF optimization seed set to", seed + 1, "\n")
  }
  
  # Convert to data.frame
  cov_df <- as.data.frame(covariates, xy = TRUE, na.rm = TRUE)
  cov_df <- na.omit(cov_df)
  
  # Start with cLHS design
  initial_samples <- simple_clhs_sampling(covariates, n_samples, seed = seed)
  current_sample <- initial_samples
  best_sample <- current_sample
  
  # Calculate initial MSE
  current_mse <- calculate_rf_mse(current_sample, cov_df)
  best_mse <- current_mse
  initial_mse <- current_mse
  
  # Simulated annealing parameters
  temperature <- 1000
  cooling_rate <- 0.95
  
  cat("Initial MSE:", round(current_mse, 4), "\n")
  
  for (i in 1:n_iterations) {
    # Generate candidate by swapping a random point
    candidate_sample <- current_sample
    replace_idx <- sample(nrow(candidate_sample), 1)
    
    # Select random replacement
    remaining_points <- setdiff(seq_len(nrow(cov_df)), 
                               match(paste(candidate_sample$x, candidate_sample$y), 
                                     paste(cov_df$x, cov_df$y)))
    new_point_idx <- sample(remaining_points, 1)
    candidate_sample[replace_idx, ] <- cov_df[new_point_idx, ]
    
    # Calculate MSE for candidate
    candidate_mse <- calculate_rf_mse(candidate_sample, cov_df)
    
    # Accept or reject candidate
    delta <- candidate_mse - current_mse
    
    if (delta < 0 || runif(1) < exp(-delta / temperature)) {
      current_sample <- candidate_sample
      current_mse <- candidate_mse
      
      if (current_mse < best_mse) {
        best_sample <- current_sample
        best_mse <- current_mse
        cat("Iteration", i, "- New best MSE:", round(best_mse, 4), "\n")
      }
    }
    
    # Cool down
    temperature <- temperature * cooling_rate
    
    # Progress update
    if (i %% 100 == 0) {
      cat("Iteration", i, "- Current MSE:", round(current_mse, 4), 
          "- Best MSE:", round(best_mse, 4), "\n")
    }
  }
  
  improvement <- ((initial_mse - best_mse) / initial_mse) * 100
  cat("Optimization completed. Improvement:", round(improvement, 2), "%\n")
  
  return(list(
    optimized_samples = best_sample,
    initial_samples = initial_samples,
    initial_mse = initial_mse,
    final_mse = best_mse,
    improvement = improvement
  ))
}
