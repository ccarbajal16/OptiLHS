# cLHS Sampling Optimization Shiny App
# Three-module application for soil sampling optimization

# Load required libraries
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(DT)
library(plotly)
library(leaflet)
library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(clhs)
library(randomForest)
library(cluster)
library(minpack.lm)
library(gridExtra)

# Get the directory where this app.R file is located
if (exists("app_dir") && !is.null(app_dir)) {
  # Use provided app_dir if available
  source_dir <- app_dir
} else {
  # Try to determine the directory of this script
  source_dir <- tryCatch({
    # Try different methods to get script location
    if (exists("rstudioapi") && requireNamespace("rstudioapi", quietly = TRUE)) {
      dirname(rstudioapi::getSourceEditorContext()$path)
    } else {
      getwd()
    }
  }, error = function(e) getwd())
}

# Source optimization functions with correct paths
source(file.path(source_dir, "optimization_functions.R"))
source(file.path(source_dir, "rf_optimization.R")) 
source(file.path(source_dir, "alternative_sites.R"))

# Source UI and server
source(file.path(source_dir, "ui.R"))
source(file.path(source_dir, "server.R"))

# Run the application
shinyApp(ui = ui, server = server)