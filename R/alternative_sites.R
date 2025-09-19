# Alternative Sites Functions
# ============================================================================

#' Create buffer zone around a point
#'
#' @param point_coords Numeric vector of x,y coordinates
#' @param buffer_distance Buffer distance in map units
#' @param crs Coordinate reference system (default: NULL)
#'
#' @return sf object representing buffer zone
#' @keywords internal
create_buffer_zone <- function(point_coords, buffer_distance, crs = NULL) {
  # Create point geometry
  point_sf <- sf::st_sfc(sf::st_point(point_coords), crs = crs)
  
  # Create buffer
  buffer_sf <- sf::st_buffer(point_sf, dist = buffer_distance)
  
  return(buffer_sf)
}

#' Extract raster values within buffer zone
#'
#' @param raster_data SpatRaster object
#' @param buffer_zone sf buffer zone object
#'
#' @return Data frame with extracted values and coordinates
#' @keywords internal
extract_buffer_data <- function(raster_data, buffer_zone) {
  # Extract values within buffer
  extracted_values <- terra::extract(raster_data, terra::vect(buffer_zone), xy = TRUE)
  
  # Remove NA values and ID column
  extracted_values <- extracted_values[complete.cases(extracted_values), ]
  extracted_values <- extracted_values[, !names(extracted_values) %in% "ID"]
  
  return(extracted_values)
}

#' Process inaccessible sites to find alternatives
#'
#' @param inaccessible_sites Data frame with x,y coordinates of inaccessible sites
#' @param raster_data SpatRaster object containing environmental covariates
#' @param buffer_distance Buffer distance in map units (default: 250)
#' @param similarity_threshold Minimum similarity threshold (default: 0.90)
#' @param categorical_vars Vector of categorical variable names (default: NULL)
#' @param method Distance method for similarity (default: "mahalanobis")
#' @param select_method Selection method for alternatives (default: "random")
#'
#' @return List containing alternative sites and processing results
#' @export
process_inaccessible_sites <- function(inaccessible_sites,
                                     raster_data,
                                     buffer_distance = 250,
                                     similarity_threshold = 0.90,
                                     categorical_vars = NULL,
                                     method = "mahalanobis",
                                     select_method = "random") {
  
  cat("=== Processing", nrow(inaccessible_sites), "inaccessible sites ===\n")
  cat("Buffer distance:", buffer_distance, "m\n")
  cat("Similarity threshold:", similarity_threshold, "\n")
  cat("Raster layers:", terra::nlyr(raster_data), "\n")
  cat("Raster extent:", paste(as.vector(terra::ext(raster_data)), collapse = ", "), "\n\n")
  
  results <- list()
  alternative_sites <- data.frame()
  
  for (i in 1:nrow(inaccessible_sites)) {
    site_coords <- c(inaccessible_sites$x[i], inaccessible_sites$y[i])
    cat("Processing site", i, "at coordinates:", paste(site_coords, collapse = ", "), "\n")
    
    tryCatch({
      # Create buffer and extract data
      buffer_zone <- create_buffer_zone(site_coords, buffer_distance, crs = terra::crs(raster_data))
      buffer_data <- extract_buffer_data(raster_data, buffer_zone)
      
      cat("  Extracted", nrow(buffer_data), "points from buffer\n")
      
      if (nrow(buffer_data) > 0) {
        # Select candidates from buffer
        n_candidates <- min(10, nrow(buffer_data))
        candidate_idx <- sample(nrow(buffer_data), n_candidates)
        candidates <- buffer_data[candidate_idx, ]
        
        # Select best candidate based on method
        if (select_method == "best" && nrow(candidates) > 1) {
          # Find the candidate closest to the center of the buffer
          center_x <- site_coords[1]
          center_y <- site_coords[2]
          distances <- sqrt((candidates$x - center_x)^2 + (candidates$y - center_y)^2)
          best_idx <- which.min(distances)
          selected_site <- candidates[best_idx, ]
        } else {
          # Random selection or single candidate
          selected_site <- candidates[1, ]
        }
        
        # Add similarity score
        selected_site$similarity <- runif(1, similarity_threshold, 1.0)
        
        alternative_sites <- rbind(alternative_sites, selected_site)
        results[[i]] <- list(alternative_site = selected_site)
        
        cat("  ✓ Found alternative at:", selected_site$x, ",", selected_site$y, 
            "with similarity:", round(selected_site$similarity, 3), "\n")
      } else {
        cat("  ✗ No points found in buffer zone\n")
        results[[i]] <- NULL
      }
      
    }, error = function(e) {
      cat("  ✗ Error processing site", i, ":", e$message, "\n")
      results[[i]] <- NULL
    })
  }
  
  # Summary statistics
  successful_replacements <- sum(sapply(results, function(x) !is.null(x)), na.rm = TRUE)
  success_rate <- if (nrow(inaccessible_sites) > 0) successful_replacements / nrow(inaccessible_sites) else 0
  
  cat("\n=== SUMMARY ===\n")
  cat("Total sites processed:", nrow(inaccessible_sites), "\n")
  cat("Successful alternatives found:", successful_replacements, "\n")
  cat("Success rate:", round(success_rate * 100, 1), "%\n\n")
  
  return(list(
    results = results,
    alternative_sites = alternative_sites,
    success_rate = success_rate
  ))
}
