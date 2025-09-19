# Core Optimization Functions for cLHS Sampling
# ============================================================================

# Silence R CMD check notes for non-standard evaluation
utils::globalVariables(c(
  "sample_size", "kl_divergence", "mean_kl", "sd_kl", "fitted_kl", "cdf", 
  ".data", "x", "y", "value"
))

#' Calculate KL divergence between sample and population distributions
calculate_kl_divergence <- function(population_data, sample_data, n_bins = 25) {
  kl_values <- numeric()
  
  # Calculate KL divergence for each continuous variable
  for (col in names(population_data)) {
    if (is.numeric(population_data[[col]])) {
      
      # Create histogram bins based on population data range
      pop_range <- range(population_data[[col]], na.rm = TRUE)
      breaks <- seq(pop_range[1], pop_range[2], length.out = n_bins + 1)
      
      # Calculate population distribution (Ei)
      pop_hist <- hist(population_data[[col]], breaks = breaks, plot = FALSE)
      pop_density <- pop_hist$counts / sum(pop_hist$counts)
      
      # Calculate sample distribution (Oi)
      sample_hist <- hist(sample_data[[col]], breaks = breaks, plot = FALSE)
      sample_density <- sample_hist$counts / sum(sample_hist$counts)
      
      # Avoid log(0) by adding small constant
      pop_density <- pmax(pop_density, 1e-10)
      sample_density <- pmax(sample_density, 1e-10)
      
      # Calculate KL divergence: KL = Σ Oi * log(Oi) - log(Ei)
      kl <- sum(sample_density * (log(sample_density) - log(pop_density)))
      kl_values <- c(kl_values, kl)
    }
  }
  
  # Return mean KL divergence across all variables
  return(mean(kl_values, na.rm = TRUE))
}

#' Optimize sample size using KL divergence approach (Malone et al. 2019)
optimize_sample_size <- function(population_data, 
                               min_samples = 10, 
                               max_samples = 500, 
                               step_size = 10,
                               n_replicates = 10,
                               n_bins = 25,
                               probability_threshold = 0.95) {
  
  # Generate sequence of sample sizes
  sample_sizes <- seq(min_samples, max_samples, by = step_size)
  
  # Initialize results storage
  results <- data.frame(
    sample_size = integer(),
    replicate = integer(),
    kl_divergence = numeric(),
    stringsAsFactors = FALSE
  )
  
  cat("Running cLHS optimization...\n")
  
  # Main optimization loop
  for (n_samples in sample_sizes) {
    cat(paste("Testing sample size:", n_samples, "\n"))
    
    for (rep in 1:n_replicates) {
      # Run cLHS sampling with error handling
      tryCatch({
        clhs_result <- clhs::clhs(population_data, size = n_samples, iter = 10000)
        sample_data <- population_data[clhs_result, ]
        
        # Calculate KL divergence
        kl_div <- calculate_kl_divergence(population_data, sample_data, n_bins)
        
        # Store results
        results <- rbind(results, data.frame(
          sample_size = n_samples,
          replicate = rep,
          kl_divergence = kl_div,
          stringsAsFactors = FALSE
        ))
        
      }, error = function(e) {
        cat(paste("Error at sample size", n_samples, "replicate", rep, ":", e$message, "\n"))
      })
    }
  }
  
  # Check if we have any results
  if (nrow(results) == 0) {
    cat("No successful results obtained\n")
    return(list(
      raw_results = results,
      summary_results = data.frame(),
      fitted_model = NULL,
      fitted_curve = NULL,
      optimal_sample_size = NA,
      plot_kl = NULL,
      plot_cdf = NULL
    ))
  }
  
  # Calculate mean KL divergence for each sample size
  summary_results <- results %>%
    dplyr::group_by(.data$sample_size) %>%
    dplyr::summarise(
      mean_kl = mean(.data$kl_divergence, na.rm = TRUE),
      sd_kl = sd(.data$kl_divergence, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Initialize variables for model fitting
  exp_model <- NULL
  fitted_curve <- NULL
  cdf_values <- NULL
  optimal_sample_size <- NA
  
  # Fit exponential decay function: y = b1 * exp(-k*x) + b0
  cat("Fitting exponential decay function...\n")
  
  if (nrow(summary_results) > 3) {
    tryCatch({
      # Initial parameter estimates
      start_params <- list(
        b0 = min(summary_results$mean_kl), 
        b1 = max(summary_results$mean_kl) - min(summary_results$mean_kl),
        k = 0.01
      )
      
      # Fit exponential model
      exp_model <- minpack.lm::nlsLM(
        mean_kl ~ b1 * exp(-k * sample_size) + b0,
        data = summary_results,
        start = start_params
      )
      
      # Generate fitted curve for plotting
      fitted_curve <- data.frame(
        sample_size = sample_sizes,
        fitted_kl = predict(exp_model, newdata = data.frame(sample_size = sample_sizes))
      )
      
      # Calculate cumulative density function of (1 - KL divergence)
      max_improvement <- max(fitted_curve$fitted_kl) - min(fitted_curve$fitted_kl)
      if (max_improvement > 1e-10) {
        cdf_values <- (max(fitted_curve$fitted_kl) - fitted_curve$fitted_kl) / max_improvement
        
        # Find optimal sample size where CDF reaches threshold
        optimal_idx <- which(cdf_values >= probability_threshold)[1]
        optimal_sample_size <- ifelse(is.na(optimal_idx), max_samples, sample_sizes[optimal_idx])
      } else {
        optimal_sample_size <- max_samples
      }
      
      cat(paste("Optimal sample size:", optimal_sample_size, "\n"))
      
    }, error = function(e) {
      cat("Error fitting exponential model:", e$message, "\n")
      optimal_sample_size <- max_samples
    })
  } else {
    cat("Not enough data points to fit exponential model\n")
    optimal_sample_size <- max_samples
  }
  
  # Create visualizations
  cat("Creating visualizations...\n")
  
  # Plot 1: KL divergence vs sample size
  p1 <- NULL
  tryCatch({
    p1 <- ggplot2::ggplot(summary_results, ggplot2::aes(x = .data$sample_size, y = .data$mean_kl)) +
      ggplot2::geom_point(size = 2, color = "blue") +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = pmax(.data$mean_kl - .data$sd_kl, 0), 
                        ymax = .data$mean_kl + .data$sd_kl), 
                    width = step_size/2, alpha = 0.7) +
      ggplot2::labs(title = "KL Divergence vs Sample Size",
           x = "Number of Samples",
           y = "KL Divergence") +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    
    # Add fitted curve if available
    if (!is.null(fitted_curve)) {
      p1 <- p1 + ggplot2::geom_line(data = fitted_curve, 
                           ggplot2::aes(x = .data$sample_size, y = .data$fitted_kl), 
                           color = "red", size = 1, linetype = "solid")
    }
    
  }, error = function(e) {
    cat("Error creating KL plot:", e$message, "\n")
    p1 <- NULL
  })
  
  # Plot 2: Cumulative density function
  p2 <- NULL
  if (!is.null(cdf_values) && !is.na(optimal_sample_size)) {
    tryCatch({
      cdf_data <- data.frame(sample_size = sample_sizes, cdf = cdf_values)
      
      p2 <- ggplot2::ggplot(cdf_data, ggplot2::aes(x = .data$sample_size, y = .data$cdf)) +
        ggplot2::geom_line(size = 1, color = "darkgreen") +
        ggplot2::geom_hline(yintercept = probability_threshold, color = "red", linetype = "dashed") +
        ggplot2::geom_vline(xintercept = optimal_sample_size, color = "red", linetype = "dashed") +
        ggplot2::labs(title = "Cumulative Density Function of (1 - KL Divergence)",
             x = "Number of Samples", 
             y = "CDF of (1 - KL divergence)") +
        ggplot2::theme_minimal() +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::annotate("text", x = optimal_sample_size + step_size * 2, y = probability_threshold - 0.05,
                 label = paste("Optimal size:", optimal_sample_size), color = "red")
      
    }, error = function(e) {
      cat("Error creating CDF plot:", e$message, "\n")
      p2 <- NULL
    })
  }
  
  # Return results
  return(list(
    raw_results = results,
    summary_results = summary_results,
    fitted_model = exp_model,
    fitted_curve = fitted_curve,
    optimal_sample_size = optimal_sample_size,
    plot_kl = p1,
    plot_cdf = p2
  ))
}
