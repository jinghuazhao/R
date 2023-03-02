#' Start shinygap
#'
#' @param ... Additional arguments passed to the 'runApp' function from the 'shiny' package.
#' @return These are design specific.
#' @export
#'
#' @details
#' This function starts the interactive 'shinygap' shiny web application that allows for flexible model specification.
#'
#' The 'shiny' based web application allows for flexible model specification for the implemented study designs.

runshinygap <- function (...)
{
    if (requireNamespace("shiny", quietly = TRUE)) {
        message("Starting the EnsDb shiny web app. Use Ctrl-C to stop.")
        shiny::runApp(appDir = system.file("shinygap",
            package = "gap"), ...)
    }
    else {
        stop("Package shiny not installed!")
    }
}

