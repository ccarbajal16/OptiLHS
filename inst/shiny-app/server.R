# Server logic for cLHS Optimization Shiny App

server <- function(input, output, session) {
  # Reactive values to store data
  values <- reactiveValues(
    rasters = NULL,
    kl_results = NULL,
    rf_results = NULL,
    alternative_results = NULL,
    final_design = NULL,
    uploaded_inaccessible = NULL,
    current_sample_data = NULL
  )
  
  # File upload reactive
  observeEvent(input$raster_files, {
    req(input$raster_files)
    
    # Load and stack rasters
    raster_paths <- input$raster_files$datapath
    names(raster_paths) <- tools::file_path_sans_ext(input$raster_files$name)
    
    tryCatch({
      values$rasters <- terra::rast(raster_paths)
      names(values$rasters) <- names(raster_paths)
      showNotification("Raster files loaded successfully!", type = "message")
    }, error = function(e) {
      showNotification(paste("Error loading rasters:", e$message), type = "error")
    })
  })
  
  # Check if files are uploaded
  output$files_uploaded <- reactive({
    !is.null(values$rasters)
  })
  outputOptions(output, "files_uploaded", suspendWhenHidden = FALSE)
  
  # Raster info table
  output$raster_info <- DT::renderDataTable({
    req(values$rasters)
    
    data.frame(
      Layer = names(values$rasters),
      Rows = nrow(values$rasters),
      Columns = ncol(values$rasters),
      Resolution = paste(terra::res(values$rasters), collapse = " x "),
      CRS = as.character(terra::crs(values$rasters))
    )
  }, options = list(pageLength = 10, scrollX = TRUE))
  
  # Raster preview plot
  output$raster_preview <- renderPlot({
    req(values$rasters)
    
    # Create preview of first raster layer
    raster_df <- as.data.frame(values$rasters[[1]], xy = TRUE, na.rm = TRUE)
    names(raster_df)[3] <- "value"
    
    ggplot2::ggplot(raster_df, ggplot2::aes(x = x, y = y, fill = value)) +
      ggplot2::geom_raster() +
      ggplot2::scale_fill_gradientn(colors = terrain.colors(255)) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = paste("Preview:", names(values$rasters)[1]),
           x = "X Coordinate", y = "Y Coordinate") +
      ggplot2::theme(aspect.ratio = 1)
  })
  
  # Module 1: KL Optimization
  observeEvent(input$run_module1, {
    req(values$rasters)
    
    shinyWidgets::updateProgressBar(session, "pb_module1", value = 10)
    
    tryCatch({
      shinyWidgets::updateProgressBar(session, "pb_module1", value = 30)
      
      # Extract population data
      population_values <- terra::values(values$rasters, na.rm = TRUE)
      population_data <- as.data.frame(population_values)
      population_data <- population_data[complete.cases(population_data), ]
      
      shinyWidgets::updateProgressBar(session, "pb_module1", value = 50)
      
      # Subsample if too large
      if (nrow(population_data) > 100000) {
        idx <- sample(nrow(population_data), 100000)
        population_data <- population_data[idx, ]
      }
      
      shinyWidgets::updateProgressBar(session, "pb_module1", value = 70)
      
      # Run optimization
      values$kl_results <- optimize_sample_size(
        population_data = population_data,
        min_samples = input$min_samples,
        max_samples = input$max_samples,
        step_size = input$step_size,
        n_replicates = input$n_replicates,
        probability_threshold = input$prob_threshold
      )
      
      shinyWidgets::updateProgressBar(session, "pb_module1", value = 100)
      showNotification("Module 1 completed successfully!", type = "message")
      
    }, error = function(e) {
      showNotification(paste("Error in Module 1:", e$message), type = "error")
    })
  })
  
  # Check if Module 1 is complete
  output$module1_complete <- reactive({
    !is.null(values$kl_results)
  })
  outputOptions(output, "module1_complete", suspendWhenHidden = FALSE)
  
  # Optimal samples value box
  output$optimal_samples <- shinydashboard::renderValueBox({
    req(values$kl_results)
    shinydashboard::valueBox(
      value = values$kl_results$optimal_sample_size,
      subtitle = "Optimal Sample Size",
      icon = shiny::icon("bullseye"),
      color = "green"
    )
  })
  
  # KL summary table
  output$kl_summary_table <- DT::renderDataTable({
    req(values$kl_results)
    values$kl_results$summary_results
  }, options = list(pageLength = 10, scrollX = TRUE))
  
  # KL plot
  output$kl_plot <- plotly::renderPlotly({
    req(values$kl_results, values$kl_results$plot_kl)
    plotly::ggplotly(values$kl_results$plot_kl)
  })
  
  # CDF plot  
  output$cdf_plot <- plotly::renderPlotly({
    req(values$kl_results, values$kl_results$plot_cdf)
    plotly::ggplotly(values$kl_results$plot_cdf)
  })
  
  # Display optimal size in Module 2
  output$optimal_size_display <- renderText({
    req(values$kl_results)
    paste("Optimal Size:", values$kl_results$optimal_sample_size, "samples")
  })
  
  # Module 2: RF Optimization
  observeEvent(input$run_module2, {
    req(values$rasters, values$kl_results)
    
    shinyWidgets::updateProgressBar(session, "pb_module2", value = 10)
    
    tryCatch({
      shinyWidgets::updateProgressBar(session, "pb_module2", value = 30)
      
      # Run RF optimization
      values$rf_results <- rf_optimize_sampling(
        covariates = values$rasters,
        n_samples = values$kl_results$optimal_sample_size,
        n_iterations = input$rf_iterations,
        seed = input$rf_seed
      )
      
      shinyWidgets::updateProgressBar(session, "pb_module2", value = 100)
      showNotification("Module 2 completed successfully!", type = "message")
      
    }, error = function(e) {
      showNotification(paste("Error in Module 2:", e$message), type = "error")
    })
  })
  
  # Check if Module 2 is complete
  output$module2_complete <- reactive({
    !is.null(values$rf_results)
  })
  outputOptions(output, "module2_complete", suspendWhenHidden = FALSE)
  
  # RF comparison table
  output$rf_comparison_table <- DT::renderDataTable({
    req(values$rf_results)
    
    data.frame(
      Method = c("Initial cLHS", "RF Optimized"),
      MSE = c(round(values$rf_results$initial_mse, 4), 
              round(values$rf_results$final_mse, 4)),
      Improvement = c("0%", paste0(round(values$rf_results$improvement, 2), "%"))
    )
  }, options = list(pageLength = 5, scrollX = TRUE))
  
  # RF improvement plot
  output$rf_improvement_plot <- plotly::renderPlotly({
    req(values$rf_results)
    
    improvement_data <- data.frame(
      Method = c("Initial cLHS", "RF Optimized"),
      MSE = c(values$rf_results$initial_mse, values$rf_results$final_mse)
    )
    
    p <- ggplot2::ggplot(improvement_data, ggplot2::aes(x = Method, y = MSE, fill = Method)) +
      ggplot2::geom_col() +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "MSE Comparison", y = "Mean Squared Error") +
      ggplot2::theme(legend.position = "none")
    
    plotly::ggplotly(p)
  })
  
  # Initialize reactive values for selected difficult sites
  selected_difficult_sites <- reactiveVal(numeric(0))

  # Interactive sampling locations map
  output$sampling_locations_map <- leaflet::renderLeaflet({
    req(values$rf_results)

    # Prepare sample data with indices
    sample_data <- values$rf_results$optimized_samples
    sample_data$Point_ID <- 1:nrow(sample_data)
    sample_data$Site_Index <- paste("Site", 1:nrow(sample_data))

    # Convert coordinates if they're in projected system
    if (mean(sample_data$x, na.rm = TRUE) > 1000) {
      sample_sf <- sf::st_as_sf(sample_data, coords = c("x", "y"),
                           crs = sf::st_crs(values$rasters))
      sample_wgs84 <- sf::st_transform(sample_sf, 4326)
      coords_wgs84 <- sf::st_coordinates(sample_wgs84)
      sample_data$lon <- coords_wgs84[, 1]
      sample_data$lat <- coords_wgs84[, 2]
    } else {
      sample_data$lon <- sample_data$x
      sample_data$lat <- sample_data$y
    }

    values$current_sample_data <- sample_data

    # Create leaflet map
    map <- leaflet::leaflet(sample_data) %>%
      leaflet::addProviderTiles("Esri.WorldImagery", group = "Satellite") %>%
      leaflet::addProviderTiles("OpenTopoMap", group = "Terrain") %>%
      leaflet::addProviderTiles("OpenStreetMap", group = "Streets") %>%
      leaflet::addLayersControl(
        baseGroups = c("Satellite", "Terrain", "Streets"),
        options = leaflet::layersControlOptions(collapsed = FALSE)
      ) %>%
      leaflet::addCircleMarkers(
        lng = ~lon, lat = ~lat,
        radius = 10,
        color = "orange",
        fillColor = "yellow",
        fillOpacity = 0.8,
        stroke = TRUE,
        weight = 3,
        layerId = ~Point_ID,
        popup = ~paste(
          "<b>", Site_Index, "</b><br>",
          "Point ID: ", Point_ID, "<br>",
          "Coordinates: ", round(lon, 6), ", ", round(lat, 6), "<br>",
          "<hr>",
          "<small>🎯 Click to select as difficult access site</small>"
        ),
        label = ~Site_Index
      ) %>%
      leaflet::fitBounds(
        lng1 = min(sample_data$lon), lat1 = min(sample_data$lat),
        lng2 = max(sample_data$lon), lat2 = max(sample_data$lat)
      )

    return(map)
  })

  # Handle map clicks for selecting difficult access sites
  observeEvent(input$sampling_locations_map_marker_click, {
    click_info <- input$sampling_locations_map_marker_click
    if (!is.null(click_info$id)) {
      clicked_point_id <- as.numeric(click_info$id)
      current_selection <- selected_difficult_sites()

      # Toggle selection
      if (clicked_point_id %in% current_selection) {
        new_selection <- current_selection[current_selection != clicked_point_id]
      } else {
        new_selection <- c(current_selection, clicked_point_id)
      }

      selected_difficult_sites(new_selection)
      updateMapMarkers()
    }
  })

  # Function to update map marker colors
  updateMapMarkers <- function() {
    req(values$current_sample_data)
    sample_data <- values$current_sample_data
    selected_ids <- selected_difficult_sites()

    for (i in 1:nrow(sample_data)) {
      point_id <- sample_data$Point_ID[i]
      if (point_id %in% selected_ids) {
        leaflet::leafletProxy("sampling_locations_map") %>%
          leaflet::removeMarker(layerId = point_id) %>%
          leaflet::addCircleMarkers(
            lng = sample_data$lon[i], lat = sample_data$lat[i],
            radius = 10, color = "darkred", fillColor = "red",
            fillOpacity = 0.9, stroke = TRUE, weight = 3,
            layerId = point_id,
            popup = paste("<b>🔴 SELECTED: ", sample_data$Site_Index[i], "</b>"),
            label = paste("🔴", sample_data$Site_Index[i])
          )
      } else {
        leaflet::leafletProxy("sampling_locations_map") %>%
          leaflet::removeMarker(layerId = point_id) %>%
          leaflet::addCircleMarkers(
            lng = sample_data$lon[i], lat = sample_data$lat[i],
            radius = 10, color = "orange", fillColor = "yellow",
            fillOpacity = 0.8, stroke = TRUE, weight = 3,
            layerId = point_id,
            popup = paste("<b>", sample_data$Site_Index[i], "</b>"),
            label = sample_data$Site_Index[i]
          )
      }
    }
  }

  # Clear all selections
  observeEvent(input$clear_difficult_selection, {
    selected_difficult_sites(numeric(0))
    updateMapMarkers()
  })

  # Selected sites count
  output$selected_sites_count <- reactive({
    length(selected_difficult_sites())
  })
  outputOptions(output, "selected_sites_count", suspendWhenHidden = FALSE)

  output$selected_sites_count_display <- renderText({
    selected_count <- length(selected_difficult_sites())
    if (selected_count > 0) {
      paste0(selected_count, " sites selected")
    } else {
      "0 sites selected"
    }
  })

  # Selected sites summary
  output$selected_sites_summary <- renderText({
    selected_count <- length(selected_difficult_sites())
    if (selected_count > 0) {
      paste0("You have selected ", selected_count, " difficult access site(s) for alternative location analysis.")
    } else {
      "No sites selected yet. Click on map points to select difficult access sites."
    }
  })

  # Difficult access sites table
  output$difficult_access_table <- DT::renderDataTable({
    req(values$current_sample_data)
    selected_ids <- selected_difficult_sites()

    if (length(selected_ids) > 0) {
      selected_sites <- values$current_sample_data[selected_ids, ]
      selected_sites$Site_ID <- selected_ids
      selected_sites[, c("Site_ID", "x", "y", setdiff(names(selected_sites), c("Site_ID", "x", "y")))]
    } else {
      # Return empty data frame with same structure
      empty_df <- values$current_sample_data[0, ]
      empty_df$Site_ID <- integer(0)
      empty_df[, c("Site_ID", "x", "y", setdiff(names(empty_df), c("Site_ID", "x", "y")))]
    }
  }, options = list(pageLength = 10, scrollX = TRUE, dom = 't'))

  # Sample locations table
  output$sample_locations_table <- DT::renderDataTable({
    req(values$rf_results)
    
    sample_data <- values$rf_results$optimized_samples
    sample_data$Point_ID <- 1:nrow(sample_data)
    sample_data[, c("Point_ID", "x", "y", setdiff(names(sample_data), c("Point_ID", "x", "y")))]
  }, options = list(pageLength = 10, scrollX = TRUE))
  
  # CSV upload handling for Module 3
  observeEvent(input$inaccessible_csv, {
    req(input$inaccessible_csv)
    
    tryCatch({
      values$uploaded_inaccessible <- read.csv(input$inaccessible_csv$datapath, stringsAsFactors = FALSE)
      
      required_columns <- c("x", "y")
      if (!all(required_columns %in% names(values$uploaded_inaccessible))) {
        showNotification("CSV must contain columns: x, y", type = "error")
        values$uploaded_inaccessible <- NULL
        return()
      }
      
      if (!is.numeric(values$uploaded_inaccessible$x) || !is.numeric(values$uploaded_inaccessible$y)) {
        showNotification("x and y columns must contain numeric coordinates", type = "error")
        values$uploaded_inaccessible <- NULL
        return()
      }
      
      showNotification(paste("Uploaded", nrow(values$uploaded_inaccessible), "inaccessible sites successfully!"), type = "message")
      
    }, error = function(e) {
      showNotification(paste("Error reading CSV file:", e$message), type = "error")
      values$uploaded_inaccessible <- NULL
    })
  })
  
  # Check if CSV is uploaded
  output$csv_uploaded <- reactive({
    !is.null(values$uploaded_inaccessible)
  })
  outputOptions(output, "csv_uploaded", suspendWhenHidden = FALSE)
  
  # Display uploaded inaccessible sites preview
  output$inaccessible_preview <- DT::renderDataTable({
    req(values$uploaded_inaccessible)
    
    preview_data <- values$uploaded_inaccessible
    preview_data$Site_ID <- 1:nrow(preview_data)
    preview_data[, c("Site_ID", "x", "y")]
  }, options = list(pageLength = 5, scrollX = TRUE, dom = 't'))
  
  # Module 3: Alternative Sites
  observeEvent(input$run_module3, {
    if (is.null(values$rasters)) {
      showNotification("Please upload raster files first", type = "error")
      return()
    }
    
    if (is.null(values$uploaded_inaccessible)) {
      showNotification("Please upload CSV file with inaccessible sites first", type = "error")
      return()
    }
    
    showNotification("Starting alternative site processing...", type = "message", duration = 2)
    
    tryCatch({
      inaccessible_sites <- values$uploaded_inaccessible[, c("x", "y")]
      
      # Find alternatives
      values$alternative_results <- process_inaccessible_sites(
        inaccessible_sites = inaccessible_sites,
        raster_data = values$rasters,
        buffer_distance = input$buffer_distance,
        similarity_threshold = input$similarity_threshold,
        method = input$distance_method,
        select_method = input$selection_method
      )
      
      # Create final design with both original and alternative sites
      final_design_data <- data.frame()

      # Add original inaccessible sites
      if (!is.null(values$uploaded_inaccessible)) {
        original_sites <- values$uploaded_inaccessible[, c("x", "y")]
        original_sites$Type <- "Original (Inaccessible)"
        original_sites$Site_Source <- "Original"
        final_design_data <- rbind(final_design_data, original_sites)
      }

      # Add alternative sites
      if (nrow(values$alternative_results$alternative_sites) > 0) {
        alternative_sites <- values$alternative_results$alternative_sites[, c("x", "y")]
        alternative_sites$Type <- "Alternative"
        alternative_sites$Site_Source <- "Alternative"
        final_design_data <- rbind(final_design_data, alternative_sites)
      }

      # Add Point_ID for all sites
      if (nrow(final_design_data) > 0) {
        final_design_data$Point_ID <- 1:nrow(final_design_data)
        values$final_design <- final_design_data
      }
      
      showNotification("Module 3 completed successfully!", type = "message")
      
    }, error = function(e) {
      showNotification(paste("Error in Module 3:", e$message), type = "error")
    })
  })
  
  # Check if Module 3 is complete
  output$module3_complete <- reactive({
    !is.null(values$alternative_results)
  })
  outputOptions(output, "module3_complete", suspendWhenHidden = FALSE)
  
  # Success rate value box
  output$success_rate <- shinydashboard::renderValueBox({
    req(values$alternative_results)
    
    success_rate <- round(values$alternative_results$success_rate * 100, 1)
    
    shinydashboard::valueBox(
      value = paste0(success_rate, "%"),
      subtitle = "Success Rate",
      icon = shiny::icon("check-circle"),
      color = if(success_rate >= 80) "green" else if(success_rate >= 60) "yellow" else "red"
    )
  })
  
  # Alternative summary table
  output$alternative_summary_table <- DT::renderDataTable({
    req(values$alternative_results)
    
    if (nrow(values$alternative_results$alternative_sites) > 0) {
      alt_data <- values$alternative_results$alternative_sites
      alt_data$Point_ID <- 1:nrow(alt_data)
      alt_data[, c("Point_ID", "x", "y", "similarity")]
    } else {
      data.frame(Message = "No alternative sites found")
    }
  }, options = list(pageLength = 10, scrollX = TRUE))

  # Alternative Sites Satellite Map
  output$alternatives_satellite_map <- leaflet::renderLeaflet({
    req(values$alternative_results)

    # Combine original inaccessible sites with alternative sites
    all_sites <- data.frame()

    # Add original inaccessible sites (red)
    if (!is.null(values$uploaded_inaccessible)) {
      original_sites <- values$uploaded_inaccessible[, c("x", "y")]
      original_sites$site_type <- "inaccessible"
      original_sites$Point_ID <- paste0("Original_", 1:nrow(original_sites))
      all_sites <- rbind(all_sites, original_sites)
    }

    # Add alternative sites (blue)
    if (!is.null(values$alternative_results$alternative_sites) &&
        nrow(values$alternative_results$alternative_sites) > 0) {
      alt_sites <- values$alternative_results$alternative_sites[, c("x", "y")]
      alt_sites$site_type <- "alternative"
      alt_sites$Point_ID <- paste0("Alternative_", 1:nrow(alt_sites))
      all_sites <- rbind(all_sites, alt_sites)
    }

    if (nrow(all_sites) == 0) {
      # Return empty map if no sites
      leaflet::leaflet() %>%
        leaflet::addProviderTiles("OpenStreetMap") %>%
        leaflet::setView(lng = 0, lat = 0, zoom = 2)
    } else {
      # Convert coordinates if they're in projected system
      if (mean(all_sites$x, na.rm = TRUE) > 1000) {
        sites_sf <- sf::st_as_sf(all_sites, coords = c("x", "y"),
                               crs = sf::st_crs(values$rasters))
        sites_wgs84 <- sf::st_transform(sites_sf, 4326)
        coords_wgs84 <- sf::st_coordinates(sites_wgs84)
        all_sites$lon <- coords_wgs84[, 1]
        all_sites$lat <- coords_wgs84[, 2]
      } else {
        all_sites$lon <- all_sites$x
        all_sites$lat <- all_sites$y
      }

      # Create color palette
      color_pal <- c("inaccessible" = "red", "alternative" = "blue")

      # Create leaflet map
      map <- leaflet::leaflet(all_sites) %>%
        leaflet::addProviderTiles("Esri.WorldImagery", group = "Satellite") %>%
        leaflet::addProviderTiles("OpenTopoMap", group = "Terrain") %>%
        leaflet::addProviderTiles("OpenStreetMap", group = "Streets") %>%
        leaflet::addLayersControl(
          baseGroups = c("Satellite", "Terrain", "Streets"),
          options = leaflet::layersControlOptions(collapsed = FALSE)
        )

      # Add markers for each site type with explicit colors
      if ("inaccessible" %in% all_sites$site_type) {
        inaccessible_data <- all_sites[all_sites$site_type == "inaccessible", ]
        map <- map %>%
          leaflet::addCircleMarkers(
            data = inaccessible_data,
            lng = ~lon, lat = ~lat,
            radius = 8,
            color = "darkred",
            fillColor = "red",
            fillOpacity = 0.8,
            weight = 2,
            popup = ~paste0("<b>Original (Inaccessible) Site</b><br>",
                          "ID: ", Point_ID, "<br>",
                          "Coordinates: ", round(x, 2), ", ", round(y, 2)),
            group = "Original"
          )
      }

      if ("alternative" %in% all_sites$site_type) {
        alternative_data <- all_sites[all_sites$site_type == "alternative", ]
        map <- map %>%
          leaflet::addCircleMarkers(
            data = alternative_data,
            lng = ~lon, lat = ~lat,
            radius = 8,
            color = "darkblue",
            fillColor = "blue",
            fillOpacity = 0.8,
            weight = 2,
            popup = ~paste0("<b>Alternative Site</b><br>",
                          "ID: ", Point_ID, "<br>",
                          "Coordinates: ", round(x, 2), ", ", round(y, 2)),
            group = "Alternative"
          )
      }

      # Add legend
      map <- map %>%
        leaflet::addLegend(
          position = "bottomright",
          colors = c("red", "blue"),
          labels = c("Original (Inaccessible)", "Alternative"),
          title = "Site Types",
          opacity = 1
        )

      map
    }
  })


  # Final design table
  output$final_design_table <- DT::renderDataTable({
    req(values$final_design)
    values$final_design
  }, options = list(pageLength = 10, scrollX = TRUE))
  
  # Summary value boxes
  output$summary_optimal_size <- shinydashboard::renderValueBox({
    req(values$kl_results)
    shinydashboard::valueBox(
      value = values$kl_results$optimal_sample_size,
      subtitle = "Optimal Sample Size",
      icon = shiny::icon("bullseye"),
      color = "blue"
    )
  })
  
  output$summary_rf_improvement <- shinydashboard::renderValueBox({
    req(values$rf_results)
    shinydashboard::valueBox(
      value = paste0(round(values$rf_results$improvement, 1), "%"),
      subtitle = "RF Improvement",
      icon = shiny::icon("arrow-up"),
      color = "green"
    )
  })
  
  output$summary_alternatives_found <- shinydashboard::renderValueBox({
    req(values$alternative_results)
    shinydashboard::valueBox(
      value = nrow(values$alternative_results$alternative_sites),
      subtitle = "Alternatives Found",
      icon = shiny::icon("map-marker-alt"),
      color = "orange"
    )
  })
  
  # Final summary plot
  output$final_summary_plot <- renderPlot({
    req(values$final_design, values$rasters)
    
    # Background raster
    raster_df <- as.data.frame(values$rasters[[1]], xy = TRUE, na.rm = TRUE)
    names(raster_df)[3] <- "value"
    
    ggplot2::ggplot() +
      ggplot2::geom_raster(data = raster_df, ggplot2::aes(x = x, y = y, fill = value), alpha = 0.7) +
      ggplot2::geom_point(data = values$final_design, ggplot2::aes(x = x, y = y, color = Type), size = 3) +
      ggplot2::scale_fill_gradientn(colors = terrain.colors(255), name = names(values$rasters)[1]) +
      ggplot2::scale_color_manual(values = c("Original (Inaccessible)" = "red", "Alternative" = "blue")) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Final Sampling Design",
           subtitle = paste("Total sites:", nrow(values$final_design)),
           x = "X Coordinate", y = "Y Coordinate") +
      ggplot2::theme(aspect.ratio = 1, legend.position = "bottom")
  })

  # Download Handlers
  # =================

  # Download raster preview plot
  output$download_raster_plot <- downloadHandler(
    filename = function() {
      paste0("raster_preview_", Sys.Date(), ".png")
    },
    content = function(file) {
      req(values$rasters)
      png(file, width = 3000, height = 2400, res = 300)
      plot(values$rasters[[1]], main = paste("Preview:", names(values$rasters)[1]))
      dev.off()
    }
  )

  # Download KL divergence plot
  output$download_kl_plot <- downloadHandler(
    filename = function() {
      paste0("kl_divergence_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      req(values$kl_results$plot_kl)
      ggplot2::ggsave(file, values$kl_results$plot_kl, width = 10, height = 6, dpi = 300)
    }
  )

  # Download CDF plot
  output$download_cdf_plot <- downloadHandler(
    filename = function() {
      paste0("cdf_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      req(values$kl_results$plot_cdf)
      ggplot2::ggsave(file, values$kl_results$plot_cdf, width = 10, height = 6, dpi = 300)
    }
  )

  # Download selected difficult sites
  output$download_difficult_sites <- downloadHandler(
    filename = function() {
      paste0("difficult_sites_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(values$current_sample_data)
      selected_ids <- selected_difficult_sites()

      if (length(selected_ids) > 0) {
        selected_sites <- values$current_sample_data[selected_ids, ]
        selected_sites$Site_ID <- selected_ids
        write.csv(selected_sites, file, row.names = FALSE)
      } else {
        # Create empty file with headers
        empty_df <- values$current_sample_data[0, ]
        empty_df$Site_ID <- integer(0)
        write.csv(empty_df, file, row.names = FALSE)
      }
    }
  )

  # Download locations plot
  output$download_locations_plot <- downloadHandler(
    filename = function() {
      paste0("sample_locations_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      req(values$current_sample_data, values$rasters)

      raster_df <- as.data.frame(values$rasters[[1]], xy = TRUE)
      names(raster_df)[3] <- "value"

      p <- ggplot2::ggplot() +
        ggplot2::geom_raster(data = raster_df, ggplot2::aes(x = x, y = y, fill = value), alpha = 0.7) +
        ggplot2::geom_point(data = values$current_sample_data, ggplot2::aes(x = x, y = y),
                           color = "red", size = 2) +
        ggplot2::scale_fill_gradientn(colors = terrain.colors(255), name = names(values$rasters)[1]) +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = "Sample Locations", x = "X Coordinate", y = "Y Coordinate") +
        ggplot2::theme(aspect.ratio = 1)

      ggplot2::ggsave(file, p, width = 10, height = 8, dpi = 300)
    }
  )

  # Download sample locations CSV
  output$download_samples_csv <- downloadHandler(
    filename = function() {
      paste0("sample_locations_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(values$current_sample_data)
      write.csv(values$current_sample_data, file, row.names = FALSE)
    }
  )

  # Download alternatives plot
  output$download_alternatives_plot <- downloadHandler(
    filename = function() {
      paste0("alternatives_map_", Sys.Date(), ".png")
    },
    content = function(file) {
      req(values$alternative_results, values$rasters)

      raster_df <- as.data.frame(values$rasters[[1]], xy = TRUE)
      names(raster_df)[3] <- "value"

      all_sites <- values$alternative_results$alternative_sites
      if (!is.null(all_sites) && nrow(all_sites) > 0) {
        p <- ggplot2::ggplot() +
          ggplot2::geom_raster(data = raster_df, ggplot2::aes(x = x, y = y, fill = value), alpha = 0.7) +
          ggplot2::geom_point(data = all_sites, ggplot2::aes(x = x, y = y, color = site_type), size = 3) +
          ggplot2::scale_fill_gradientn(colors = terrain.colors(255), name = names(values$rasters)[1]) +
          ggplot2::scale_color_manual(values = c("inaccessible" = "red", "alternative" = "blue")) +
          ggplot2::theme_minimal() +
          ggplot2::labs(title = "Alternative Sites Analysis", x = "X Coordinate", y = "Y Coordinate") +
          ggplot2::theme(aspect.ratio = 1, legend.position = "bottom")

        ggplot2::ggsave(file, p, width = 10, height = 8, dpi = 300)
      }
    }
  )

  # Download final design CSV
  output$download_final_design <- downloadHandler(
    filename = function() {
      paste0("final_sampling_design_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(values$final_design)
      write.csv(values$final_design, file, row.names = FALSE)
    }
  )

  # Download summary plot
  output$download_summary_plot <- downloadHandler(
    filename = function() {
      paste0("summary_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      req(values$final_design, values$rasters)

      raster_df <- as.data.frame(values$rasters[[1]], xy = TRUE)
      names(raster_df)[3] <- "value"

      p <- ggplot2::ggplot() +
        ggplot2::geom_raster(data = raster_df, ggplot2::aes(x = x, y = y, fill = value), alpha = 0.7) +
        ggplot2::geom_point(data = values$final_design, ggplot2::aes(x = x, y = y, color = Type), size = 3) +
        ggplot2::scale_fill_gradientn(colors = terrain.colors(255), name = names(values$rasters)[1]) +
        ggplot2::scale_color_manual(values = c("Original (Inaccessible)" = "red", "Alternative" = "blue")) +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = "Final Sampling Design",
             subtitle = paste("Total sites:", nrow(values$final_design)),
             x = "X Coordinate", y = "Y Coordinate") +
        ggplot2::theme(aspect.ratio = 1, legend.position = "bottom")

      ggplot2::ggsave(file, p, width = 12, height = 10, dpi = 300)
    }
  )

  # Download all results as ZIP
  output$download_all_results <- downloadHandler(
    filename = function() {
      paste0("optilhs_results_", Sys.Date(), ".zip")
    },
    content = function(file) {
      # Create temporary directory
      temp_dir <- tempdir()

      # Save all available data
      files_to_zip <- c()

      # Save KL results if available
      if (!is.null(values$kl_results)) {
        kl_file <- file.path(temp_dir, "kl_divergence_results.csv")
        write.csv(values$kl_results$summary_results, kl_file, row.names = FALSE)
        files_to_zip <- c(files_to_zip, kl_file)

        # Save KL plot
        if (!is.null(values$kl_results$plot_kl)) {
          kl_plot_file <- file.path(temp_dir, "kl_divergence_plot.png")
          ggplot2::ggsave(kl_plot_file, values$kl_results$plot_kl, width = 10, height = 6, dpi = 300)
          files_to_zip <- c(files_to_zip, kl_plot_file)
        }
      }

      # Save RF results if available
      if (!is.null(values$rf_results)) {
        rf_file <- file.path(temp_dir, "rf_optimization_results.csv")
        write.csv(values$rf_results$optimized_samples, rf_file, row.names = FALSE)
        files_to_zip <- c(files_to_zip, rf_file)
      }

      # Save current sample data if available
      if (!is.null(values$current_sample_data)) {
        samples_file <- file.path(temp_dir, "sample_locations.csv")
        write.csv(values$current_sample_data, samples_file, row.names = FALSE)
        files_to_zip <- c(files_to_zip, samples_file)
      }

      # Save final design if available
      if (!is.null(values$final_design)) {
        final_file <- file.path(temp_dir, "final_sampling_design.csv")
        write.csv(values$final_design, final_file, row.names = FALSE)
        files_to_zip <- c(files_to_zip, final_file)
      }

      # Save alternative results if available
      if (!is.null(values$alternative_results$alternative_sites)) {
        alt_file <- file.path(temp_dir, "alternative_sites.csv")
        write.csv(values$alternative_results$alternative_sites, alt_file, row.names = FALSE)
        files_to_zip <- c(files_to_zip, alt_file)
      }

      # Create ZIP file
      if (length(files_to_zip) > 0) {
        zip(file, files_to_zip, flags = "-j")
      }
    }
  )
}

# Return the server function
server
