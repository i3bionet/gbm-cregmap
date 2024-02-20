library(shiny)
library(waiter)

#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request) {
    tagList(
        golem_add_external_resources(),
        waiter::useWaiter(),
        waiter::waiterShowOnLoad(html = waiter::spin_fading_circles()),
        div(
            id = "app_content",
            navbarPage(
                position = "static-top",
                theme = shinythemes::shinytheme("journal"),
                title="",
                tabPanel(title = "",
                    icon = div(img(src='img/logo/logo_gbm_cRegMap_resized.jpg',
                        style="float:left; margin-left: 5px; margin-right: 5px; margin-top: -15px")),
                    mod_welcomePage_ui("welcomePage_1")
                ),
                tabPanel("Data sets",
                    mod_Data_sets_ui("Data_sets_1")
                ),
                tabPanel("CoRegNet",
                    fluidRow(
                        fluidRow(style = "padding: 20px",
                            column(6,
                                mod_heterogeneityUserSelectionOriginExclusiveChoice_ui("heterogeneityUserSelectionOriginExclusiveChoice_1"),
                                align = "center"),
                            column(6,
                                mod_heterogeneityUserSelectionClusteringChoice_ui("heterogeneityUserSelectionClusteringChoice_ui_1"),
                                align = "center"
                            )),
                        fluidRow(style = "padding-left: 40px",
                            mod_heterogeneityUserSelectionInfluenceProjection_ui("heterogeneityUserSelectionInfluenceProjection_ui_1"),
                            align = "center"
                        ),
                        fluidRow(style = "border: 4px solid #3E3E3A; padding: 24px; margin:10px",
                            column(11,
                                mod_nertwork_viz_ui("nertwork_viz_ui_1")
                            ),
                            column(1,
                                br(),
                                br(),
                                br(),
                                br(),
                                br(),
                                br(),
                                br(),
                                br(),
                                div(
                                    img(src = "img/netVizLegend.jpg",
                                        style = "width: 100px; vertical-align:middle")
                                ),
                                align = "center",
                                tags$style(HTML('.verticalcenter {
                                                          vertical-align: middle;}')
                                )
                            )
                        ),
                    )
                ),
                tabPanel("Tutorial",
                    mod_FAQ_ui("FAQ_1")),
            )
        )
    )
}

#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function(){

    golem::add_resource_path(
        'www', '/inst/app/www'
    )

    golem::add_resource_path(
        'img', '/inst/app/img'
        # 'img', '/Volumes/Cruxial_X8/Thesis/projects/RegMap/gbm.cregmap/inst/app/img/'
    )

    tags$head(
        favicon(ext="png"),
        bundle_resources(
            path = app_sys('app/www'),
            app_title = 'RegMap'
        )
        # Add here other external resources
        # for example, you can add shinyalert::useShinyalert()
    )
}

