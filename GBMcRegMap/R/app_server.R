#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
#'


# Main function
app_server <- function( input, output, session ) {
    session$onSessionEnded(stopApp)
    print("SERVER INFO :")
    print(paste("SHINY_PORT :",Sys.getenv("SHINY_PORT")))
    print(paste("SHINY_SERVER_VERSION :",Sys.getenv("SHINY_SERVER_VERSION")))
    timeout <- 60
    options(timeout = timeout)
    options(shiny.maxRequestSize = -1)
    r <- reactiveValues()
    r$mchRMO <- GBMRegMap::mchRMO
    r$ccleRMO <- GBMRegMap::ccleRMO
    r$cregmapVersionUpper <- "GBM"
    r$cregmapCancerName <- "glioblastoma"
    mod_Data_sets_server("Data_sets_1", r=r)
    mod_heterogeneityUserSelectionOriginExclusiveChoice_server("heterogeneityUserSelectionOriginExclusiveChoice_1", r=r)
    mod_heterogeneityUserSelectionClusteringChoice_server("heterogeneityUserSelectionClusteringChoice_ui_1", r=r)
    mod_heterogeneityUserSelectionInfluenceProjection_server("heterogeneityUserSelectionInfluenceProjection_ui_1", r=r)
    mod_nertwork_viz_server("nertwork_viz_ui_1", r=r)
    mod_FAQ_server("FAQ_1", r=r)
    waiter::waiter_hide()
}
