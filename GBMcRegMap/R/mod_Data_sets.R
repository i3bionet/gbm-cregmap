#' Data_sets UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_Data_sets_ui <- function(id){
    ns <- NS(id)
    tagList(
        uiOutput(ns("datasetUI"))

    )
}

#' Data_sets Server Functions
#'
#' @noRd
mod_Data_sets_server <- function(id, r){
    moduleServer( id, function(input, output, session){
        ns <- session$ns
        ## -------------------------------------------------------------------------- ##
        # render UI
        output$datasetUI <- renderUI({
            fluidRow(
                p(
                    span(paste0("Transcriptomic datasets used to generate ",tolower(isolate(r$cregmapVersionUpper)), ".cregmap reference compounds ",
                                isolate(r$cregmapVersionUpper),
                                "-CoRegNet and ",isolate(r$cregmapVersionUpper),"-CoRegMap")
                         , style = "padding:30px;font-weight:bold;font-size:large"),
                    style = "padding:30px"),
                p(h6(paste0("Meta-cohort of ",isolate(r$cregmapVersionUpper), ", tumors  (",
                            length(unique(isolate(r$mchRMO@descriptiveData$Dataset))),
                            " studies, n=",
                            ncol(isolate(r$mchRMO@influence)),")"),
                     style = "padding-left:10%"),
                  div(
                      # DT::DTOutput(outputId = ns("mchSummaryTable")),
                      # dataTableOutput(outputId = ns("mchSummaryTable")),
                      fluidRow(align = "center", tableOutput(ns('mchSummaryTable'))),
                      style = "font-size:80%;padding-left:10%; padding-right:10%"),
                  br(),
                ),
                p(h6(paste0("CCLE ",isolate(r$cregmapVersionUpper)," cell lines (n=",
                            nrow(unique(isolate(r$ccleRMO@descriptiveData))),
                            ")"),
                     style = "padding-left:10%"),
                  div(
                      # dataTableOutput(outputId = ns("ccleSummaryTable")),
                      fluidRow(align = "center", tableOutput(ns('ccleSummaryTable'))),
                      style = "font-size:80%;padding-left:10%; padding-right:10%"),
                  br(),
                ),



            )
        })

        # store datasets
        # mchInf <- BRCARegMap::mchRMO@influence
        mchRMO <- isolate(r$mchRMO)
        # ccleInf <- BRCARegMap::ccleRMO@influence
        ccleRMO <- isolate(r$ccleRMO)
        ## -------------------------------------------------------------------------- ##
        # main Panel
        # output$mchSummaryTable <- renderDataTable(BRCARegMap::mchRMO@descriptiveData, options = list(pageLength = 5))
        # output$ccleSummaryTable <- renderDataTable(BRCARegMap::ccleRMO@descriptiveData, options = list(pageLength = 5))
        output$mchSummaryTable <- renderTable(mchRMO@descriptiveData, options = list(pageLength = 5))
        output$ccleSummaryTable <- renderTable(ccleRMO@descriptiveData, options = list(pageLength = 5))



    })
}

## To be copied in the UI
# mod_Data_sets_ui("Data_sets_1")

## To be copied in the server
# mod_Data_sets_server("Data_sets_1", r=r)
