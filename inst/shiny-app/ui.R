# UI for cLHS Sampling Optimization Shiny App

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

# Get the directory where this ui.R file is located
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

# Define UI
ui <- dashboardPage(
  dashboardHeader(
    title = tags$div(
      style = "display: flex; align-items: center; gap: 10px;",
      tags$img(
        src = "logo.svg",
        height = "40px",
        width = "auto",
        style = "vertical-align: middle;"
      ),
      "OptiLHS"
    )
  ),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data Upload", tabName = "upload", icon = icon("upload")),
      menuItem("Module 1: Sample Size", tabName = "module1", icon = icon("calculator")),
      menuItem("Module 2: RF Optimization", tabName = "module2", icon = icon("tree")),
      menuItem("Module 3: Alternative Sites", tabName = "module3", icon = icon("map-marker-alt")),
      menuItem("Results Summary", tabName = "summary", icon = icon("chart-bar")),
      menuItem("Documentation", tabName = "documentation", icon = icon("book"))
    )
  ),
  
  dashboardBody(
    tags$head(
      # Viewport meta tag for mobile responsiveness
      tags$meta(name = "viewport", content = "width=device-width, initial-scale=1.0"),

      tags$style(HTML("
        .content-wrapper, .right-side {
          background-color: #f4f4f4;
        }
        .box {
          margin-bottom: 20px;
        }
        .progress-bar {
          background-color: #3c8dbc;
        }

        /* Responsive header and logo */
        @media (max-width: 768px) {
          .main-header .navbar-custom-menu > .navbar-nav > li > .dropdown-menu {
            left: auto;
            right: 0;
          }

          /* Adjust logo size on mobile */
          .main-header .navbar-brand img {
            height: 30px !important;
          }

          /* Ensure header text doesn't wrap */
          .main-header .navbar-brand {
            white-space: nowrap;
            overflow: hidden;
          }
        }

        @media (max-width: 480px) {
          /* Even smaller logo for very small screens */
          .main-header .navbar-brand img {
            height: 25px !important;
          }

          /* Adjust header title font size */
          .main-header .navbar-brand {
            font-size: 16px;
          }
        }

        /* Responsive boxes and content */
        @media (max-width: 768px) {
          .box {
            margin-bottom: 15px;
          }

          /* Stack columns on mobile */
          .col-md-6, .col-md-4, .col-md-8 {
            margin-bottom: 15px;
          }
        }

        /* Responsive tables */
        @media (max-width: 768px) {
          .dataTables_wrapper {
            overflow-x: auto;
          }
        }

        /* Responsive plots */
        @media (max-width: 768px) {
          .plotly, .leaflet {
            height: 300px !important;
          }
        }
      "))
    ),
    
    tabItems(
      # Data Upload Tab
      tabItem(tabName = "upload",
        fluidRow(
          box(
            title = "Upload Environmental Raster Files", 
            status = "primary", 
            solidHeader = TRUE,
            width = 12,
            
            fileInput("raster_files", 
                     "Choose Raster Files (.tif)",
                     multiple = TRUE,
                     accept = c(".tif", ".tiff")),
            
            conditionalPanel(
              condition = "output.files_uploaded",
              h4("Uploaded Rasters:"),
              DT::dataTableOutput("raster_info"),
              br(),
              plotOutput("raster_preview", height = "400px"),
              br(),
              downloadButton("download_raster_plot", "Download Preview (300 DPI)", 
                           class = "btn-primary")
            )
          )
        )
      ),
      
      # Module 1: Sample Size Optimization
      tabItem(tabName = "module1",
        fluidRow(
          box(
            title = "KL Divergence Sample Size Optimization",
            status = "primary",
            solidHeader = TRUE,
            width = 4,
            
            conditionalPanel(
              condition = "!output.files_uploaded",
              div(style = "text-align: center; padding: 20px;",
                  h4("Please upload raster files first", style = "color: #dd4b39;"),
                  icon("exclamation-triangle", style = "font-size: 48px; color: #dd4b39;")
              )
            ),
            
            conditionalPanel(
              condition = "output.files_uploaded",
              numericInput("min_samples", "Minimum Samples:", value = 10, min = 5, max = 100),
              numericInput("max_samples", "Maximum Samples:", value = 500, min = 50, max = 1000),
              numericInput("step_size", "Step Size:", value = 10, min = 5, max = 50),
              numericInput("n_replicates", "Number of Replicates:", value = 10, min = 5, max = 20),
              sliderInput("prob_threshold", "Probability Threshold:", 
                         min = 0.85, max = 0.99, value = 0.95, step = 0.01),
              br(),
              actionButton("run_module1", "Run Optimization", 
                          class = "btn-success btn-block", icon = icon("play")),
              br(),
              conditionalPanel(
                condition = "input.run_module1 > 0",
                div(id = "progress_module1", 
                    h4("Processing..."),
                    progressBar("pb_module1", value = 0, display_pct = TRUE)
                )
              )
            )
          ),
          
          box(
            title = "Results",
            status = "success",
            solidHeader = TRUE,
            width = 8,
            
            conditionalPanel(
              condition = "output.module1_complete",
              tabsetPanel(
                tabPanel("Optimal Sample Size",
                  br(),
                  valueBoxOutput("optimal_samples", width = 12),
                  br(),
                  DT::dataTableOutput("kl_summary_table")
                ),
                tabPanel("KL Divergence Plot",
                  br(),
                  plotlyOutput("kl_plot", height = "500px"),
                  br(),
                  downloadButton("download_kl_plot", "Download Plot (300 DPI)", 
                               class = "btn-primary")
                ),
                tabPanel("CDF Plot",
                  br(),
                  plotlyOutput("cdf_plot", height = "500px"),
                  br(),
                  downloadButton("download_cdf_plot", "Download Plot (300 DPI)", 
                               class = "btn-primary")
                )
              )
            )
          )
        )
      ),
      
      # Module 2: RF Optimization
      tabItem(tabName = "module2",
        fluidRow(
          box(
            title = "Random Forest Sample Optimization",
            status = "primary",
            solidHeader = TRUE,
            width = 4,
            
            conditionalPanel(
              condition = "!output.module1_complete",
              div(style = "text-align: center; padding: 20px;",
                  h4("Complete Module 1 first", style = "color: #dd4b39;"),
                  icon("exclamation-triangle", style = "font-size: 48px; color: #dd4b39;")
              )
            ),
            
            conditionalPanel(
              condition = "output.module1_complete",
              h4("Use Optimal Sample Size:"),
              textOutput("optimal_size_display"),
              br(),
              numericInput("rf_iterations", "RF Iterations:", 
                          value = 500, min = 100, max = 2000),
              numericInput("rf_seed", "Random Seed:", 
                          value = 123, min = 1, max = 9999),
              br(),
              actionButton("run_module2", "Run RF Optimization", 
                          class = "btn-success btn-block", icon = icon("tree")),
              br(),
              conditionalPanel(
                condition = "input.run_module2 > 0",
                div(id = "progress_module2",
                    h4("Processing..."),
                    progressBar("pb_module2", value = 0, display_pct = TRUE)
                )
              )
            )
          ),
          
          box(
            title = "Results",
            status = "success", 
            solidHeader = TRUE,
            width = 8,
            
            conditionalPanel(
              condition = "output.module2_complete",
              tabsetPanel(
                tabPanel("Comparison",
                  br(),
                  DT::dataTableOutput("rf_comparison_table"),
                  br(),
                  plotlyOutput("rf_improvement_plot", height = "300px")
                ),
                tabPanel("Sample Locations",
                  br(),
                  leafletOutput("sampling_locations_map", height = "500px"),
                  br(),
                  h5("🎯 Interactive Site Selection"),
                  wellPanel(
                    style = "background-color: #e8f4fd; border: 1px solid #bee5eb;",
                    h6("📍 Click-to-Select Method:", style = "margin-top: 0; color: #0c5460;"),
                    p("• Click directly on map points to select/deselect difficult access sites",
                      style = "margin: 5px 0; color: #0c5460;"),
                    p("• Selected sites appear in the table below with coordinates",
                      style = "margin: 5px 0; color: #0c5460;"),
                    p("• 🔴 Red points = Selected as difficult access",
                      style = "margin: 5px 0; color: #0c5460;"),
                    p("• 🟡 Yellow points = Available for selection",
                      style = "margin: 5px 0; color: #0c5460;")
                  ),

                  fluidRow(
                    column(12,
                      div(style = "text-align: center;",
                        h5("🎮 Interactive Controls"),
                        actionButton("clear_difficult_selection", "Clear All Selections",
                                   class = "btn-warning btn-sm", style = "margin-right: 10px;"),
                        span(style = "margin: 0 15px; color: #666;", "|"),
                        strong("Selected Sites: "),
                        textOutput("selected_sites_count_display", inline = TRUE)
                      )
                    )
                  ),

                  conditionalPanel(
                    condition = "output.selected_sites_count > 0",
                    hr(),
                    h5("🗺️ Selected Difficult Access Sites:"),
                    div(id = "selected_sites_info",
                      textOutput("selected_sites_summary"),
                      style = "background-color: #f8f9fa; padding: 10px; border-radius: 4px; margin-bottom: 10px;"
                    ),
                    DT::dataTableOutput("difficult_access_table"),
                    br(),
                    downloadButton("download_difficult_sites", "Download Selected Sites (CSV)",
                                   class = "btn-primary")
                  ),
                  br(),
                  downloadButton("download_locations_plot", "Download Plot (300 DPI)", 
                               class = "btn-primary")
                ),
                tabPanel("Sample Data",
                  br(),
                  DT::dataTableOutput("sample_locations_table"),
                  br(),
                  downloadButton("download_samples_csv", "Download Sample Locations (CSV)", 
                               class = "btn-primary")
                )
              )
            )
          )
        )
      ),
      
      # Module 3: Alternative Sites
      tabItem(tabName = "module3",
        fluidRow(
          box(
            title = "Alternative Site Selection",
            status = "primary",
            solidHeader = TRUE,
            width = 4,
            
            conditionalPanel(
              condition = "!output.module2_complete && !output.files_uploaded",
              div(style = "text-align: center; padding: 20px;",
                  h4("Upload raster files first", style = "color: #dd4b39;"),
                  icon("exclamation-triangle", style = "font-size: 48px; color: #dd4b39;")
              )
            ),
            
            conditionalPanel(
              condition = "output.files_uploaded",
              h4("Upload Inaccessible Sites:"),
              
              fileInput("inaccessible_csv", 
                       "Choose CSV File:",
                       accept = c(".csv")),
              helpText("CSV file should contain columns: x, y (coordinates of inaccessible sites)"),
              
              conditionalPanel(
                condition = "output.csv_uploaded",
                h5("Uploaded Inaccessible Sites:"),
                DT::dataTableOutput("inaccessible_preview", height = "200px")
              ),
              br(),
              h4("Parameters:"),
              sliderInput("buffer_distance", "Buffer Distance (m):", 
                         min = 100, max = 500, value = 250, step = 50),
              sliderInput("similarity_threshold", "Similarity Threshold:", 
                         min = 0.70, max = 0.95, value = 0.90, step = 0.05),
              selectInput("distance_method", "Distance Method:",
                         choices = list("Mahalanobis" = "mahalanobis",
                                       "Gower" = "gower"),
                         selected = "mahalanobis"),
              selectInput("selection_method", "Selection Method:",
                         choices = list("Best Similarity" = "best",
                                       "Nearest to Original" = "nearest", 
                                       "Random" = "random"),
                         selected = "best"),
              br(),
              actionButton("run_module3", "Find Alternatives", 
                          class = "btn-success btn-block", icon = icon("search")),
              br(),
              conditionalPanel(
                condition = "input.run_module3 > 0",
                div(id = "progress_module3",
                    h4("Processing alternative sites..."),
                    p("Check notifications for progress updates.")
                )
              )
            )
          ),
          
          box(
            title = "Results",
            status = "success",
            solidHeader = TRUE, 
            width = 8,
            
            conditionalPanel(
              condition = "output.module3_complete",
              tabsetPanel(
                tabPanel("Summary",
                  br(),
                  valueBoxOutput("success_rate", width = 12),
                  br(),
                  DT::dataTableOutput("alternative_summary_table")
                ),
                tabPanel("Alternative Sites Map",
                  br(),
                  # Full-width map display
                  fluidRow(
                    column(12,
                      div(
                        style = "border: 1px solid #ddd; border-radius: 4px;",
                        leafletOutput("alternatives_satellite_map", height = "600px")
                      )
                    )
                  ),
                  br(),
                  fluidRow(
                    column(6,
                      downloadButton("download_alternatives_plot", "Download Map (300 DPI)",
                                   class = "btn-primary")
                    ),
                    column(6,
                      div(class = "pull-right",
                        p(style = "margin-top: 7px; color: #666;",
                          "🛰️ Use controls to manage point visibility • Red: Original • Blue: Alternative"
                        )
                      )
                    )
                  )
                ),
                tabPanel("Final Design",
                  br(),
                  DT::dataTableOutput("final_design_table"),
                  br(),
                  downloadButton("download_final_design", "Download Final Design (CSV)", 
                               class = "btn-primary")
                )
              )
            )
          )
        )
      ),
      
      # Results Summary
      tabItem(tabName = "summary",
        fluidRow(
          box(
            title = "Project Summary",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            
            conditionalPanel(
              condition = "output.module3_complete",
              h3("Complete Workflow Results"),
              br(),
              fluidRow(
                column(4,
                  valueBoxOutput("summary_optimal_size", width = NULL)
                ),
                column(4,
                  valueBoxOutput("summary_rf_improvement", width = NULL)  
                ),
                column(4,
                  valueBoxOutput("summary_alternatives_found", width = NULL)
                )
              ),
              br(),
              tabsetPanel(
                tabPanel("Final Sampling Design",
                  br(),
                  plotOutput("final_summary_plot", height = "600px"),
                  br(),
                  downloadButton("download_summary_plot", "Download Summary Plot (300 DPI)", 
                               class = "btn-primary")
                ),
                tabPanel("Export All Results",
                  br(),
                  h4("Download Complete Results Package:"),
                  br(),
                  downloadButton("download_all_results", "Download All Results (ZIP)", 
                               class = "btn-success btn-lg",
                               icon = icon("download"))
                )
              )
            ),
            
            conditionalPanel(
              condition = "!output.module3_complete",
              div(style = "text-align: center; padding: 40px;",
                  h4("Complete all three modules to see summary", style = "color: #999;"),
                  icon("info-circle", style = "font-size: 48px; color: #999;")
              )
            )
          )
        )
      ),
      
      # Documentation
      tabItem(tabName = "documentation",
        fluidRow(
          box(
            title = "Documentation and References",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            
            h3("Main References"),
            br(),
            
            h4("1. Wadoux, Alexandre M.J-C.; Brus, Dick J.; Heuvelink, Gerard B.M. (2019)"),
            p("Sampling design optimization for soil mapping with random forest."),
            p(strong("Journal:"), "Geoderma, 355, 113913."),
            p(strong("DOI:"), tags$a(href = "https://doi.org/10.1016/j.geoderma.2019.113913", 
                                    "https://doi.org/10.1016/j.geoderma.2019.113913", target = "_blank")),
            p("This paper provides foundational methods for digital soil mapping using R, including sampling optimization techniques that form the basis of the cLHS approach implemented in this application."),
            br(),

            h4("2. Malone BP, Minansy B, Brungard C. (2019)"),
            p("Some methods to improve the utility of conditioned Latin hypercube sampling."),
            p(strong("Journal:"), "PeerJ, 3, e6451."),
            p(strong("DOI:"), tags$a(href = "https://doi.org/10.7717/peerj.6451", 
                                    "https://doi.org/10.7717/peerj.6451", target = "_blank")),
            p("This study addresses solutions to problems that field scientists face when using cLHS. These problems include optimizing the sample size, re-locating sites when an original site is deemed inaccessible."),
            
            h4("Application Overview"),
            p("This Shiny application implements a three-module workflow for optimizing soil sampling designs using conditioned Latin Hypercube Sampling (cLHS) and Random Forest optimization:"),
            tags$ul(
              tags$li(strong("Module 1:"), " Determines optimal sample size using KL divergence analysis"),
              tags$li(strong("Module 2:"), " Refines sampling locations using Random Forest optimization"),
              tags$li(strong("Module 3:"), " Finds alternative sites for inaccessible locations")
            ),
            br(),
            
            h4("Technical Implementation"),
            p("The application uses R packages including: terra, sf, clhs, randomForest, ggplot2, plotly, leaflet, and shinydashboard for spatial data processing, optimization, visualization, and interactive web interface."),
            br(),
            
            h4("Contact Information"),
            p("For questions or feedback about this application, please refer to the original research papers or contact the development team.")
          )
        )
      )
    )
  )
)