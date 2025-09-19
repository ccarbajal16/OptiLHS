#' Launch OptiLHS Optimization Shiny App
#'
#' @description
#' Launches the cLHS (conditioned Latin Hypercube Sampling) Optimization Shiny application.
#' This app implements a three-module workflow for optimizing soil sampling designs:
#' \itemize{
#'   \item Module 1: Sample size optimization using KL divergence analysis
#'   \item Module 2: Random Forest-based sampling optimization  
#'   \item Module 3: Alternative site selection for inaccessible locations
#' }
#'
#' @param ... Additional arguments passed to \code{\link[shiny]{runApp}}
#' @param port The port number for the Shiny app (default: NULL for automatic)
#' @param host The host IP address (default: "127.0.0.1")
#' @param launch.browser Whether to launch the app in browser (default: TRUE)
#'
#' @details
#' The application requires environmental raster files (.tif format) as input and provides
#' comprehensive tools for:
#' \itemize{
#'   \item Determining optimal sample sizes using KL divergence methodology (Malone et al. 2019)
#'   \item Refining sampling locations using Random Forest optimization techniques
#'   \item Finding alternative sites when original locations are inaccessible
#'   \item Interactive visualization and export of results
#' }
#'
#' @references
#' Malone BP, Minansy B, Brungard C. (2019). Some methods to improve the utility of 
#' conditioned Latin hypercube sampling. PeerJ, 3, e6451. 
#' \doi{10.7717/peerj.6451}
#'
#' Wadoux, Alexandre M.J-C.; Brus, Dick J.; Heuvelink, Gerard B.M. (2019).
#' Sampling design optimization for soil mapping with random forest.
#' Geoderma, 355, 113913. \doi{10.1016/j.geoderma.2019.113913}
#'
#' @return No return value, launches the Shiny application
#' @export
#'
#' @examples
#' \dontrun{
#' # Launch the app with default settings
#' launch_optilhs_app()
#' 
#' # Launch on specific port
#' launch_optilhs_app(port = 3838)
#' 
#' # Launch without opening browser
#' launch_optilhs_app(launch.browser = FALSE)
#' }
launch_optilhs_app <- function(..., port = NULL, host = "127.0.0.1", launch.browser = TRUE) {
  
  # Load required libraries
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required but not installed.")
  }
  if (!requireNamespace("shinydashboard", quietly = TRUE)) {
    stop("Package 'shinydashboard' is required but not installed.")
  }
  
  # Get the app directory
  app_dir <- system.file("shiny-app", package = "OptiLHS")
  
  # If package is not installed, try to find the app in development mode
  if (app_dir == "") {
    # Try to find the app directory in the current package structure
    possible_paths <- c(
      file.path(getwd(), "inst", "shiny-app"),
      file.path(dirname(getwd()), "inst", "shiny-app"),
      file.path("..", "inst", "shiny-app"),
      file.path("inst", "shiny-app")
    )
    
    for (path in possible_paths) {
      if (dir.exists(path) && file.exists(file.path(path, "app.R"))) {
        app_dir <- path
        break
      }
    }
    
    if (app_dir == "") {
      stop("Could not find Shiny app directory. Please ensure you are in the package directory or install the package.")
    }
  }
  
  # Run the app
  shiny::runApp(
    appDir = app_dir,
    port = port,
    host = host,
    launch.browser = launch.browser,
    ...
  )
}
