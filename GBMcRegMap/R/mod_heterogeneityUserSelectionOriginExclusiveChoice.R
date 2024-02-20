#' heterogeneityUserSelectionOriginExclusiveChoice UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_heterogeneityUserSelectionOriginExclusiveChoice_ui <- function(id){
    ns <- NS(id)
    tagList(
        shinyWidgets::awesomeRadio(
            inputId = ns("exclusiveDataChoice"),
            label = h4("1. Select data "),
            choices = c("Meta-cohort of GBM tumors  (16 studies, n=1612)",
                        "CCLE GBM cell lines (n=42)"),
            status = "success"
        )

    )
}

#' heterogeneityUserSelectionOriginExclusiveChoice Server Functions
#'
#' @noRd
mod_heterogeneityUserSelectionOriginExclusiveChoice_server <- function(id, r){
    moduleServer( id, function(input, output, session){
        ns <- session$ns

        originChoices <- c("Meta-cohort of GBM tumors  (16 studies, n=1612)",
                           "CCLE GBM cell lines (n=42)")

        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # Update choices
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ## Update origin choices if user upload data or not
        # observeEvent(r$leftPanel_UserInput, {
        #     ### If user uploaded some RNAseq data
        #     if (!is.null(r$leftPanel_UserInput )){
        #         updateChoices <- c(originChoices,"Uploaded data")
        #         r$oldOrigins <- updateChoices
        #         # print("updated oldOrigins when uploaded data")
        #         # print(r$oldOrigins)
        #         shinyWidgets::updateAwesomeRadio(
        #             session = session,
        #             inputId = "exclusiveDataChoice",
        #             choices = updateChoices)
        #     }
        #
        # })
        # ## Update origin if user has selected a toyset
        # observeEvent(r$toysets,{
        #     ### If a toyset has been selected
        #     if (!is.null(r$toysets)){
        #         ### If dataOrigin has already been updated
        #         if (!is.null(r$oldOrigins)){
        #             ### If in data oriign there is "Uploaded data"
        #             if (!any(grepl(pattern = "Uploaded", r$oldOrigins))){
        #                 updateOrigins <- c(originChoices, r$toysets)
        #             } else {
        #                 updateOrigins <- c(originChoices,"Uploaded data", r$toysets)
        #             }
        #         } else {
        #             updateOrigins <- c(originChoices, r$toysets)
        #         }
        #         ### Store updates to keep in memory previous state
        #         r$oldOrigins <- updateOrigins
        #         shinyWidgets::updateAwesomeRadio(
        #             session = session,
        #             inputId = "exclusiveDataChoice",
        #             choices = updateOrigins)
        #     } else { ### If no toyset has been selected
        #         # If r$oldOrigins doesn't exists that means that no data has been
        #         # uploaded
        #         if (!is.null(r$oldOrigins)){
        #             ### If no data has been uploaded
        #             if (!any(grepl(pattern = "Uploaded", r$oldOrigins))){
        #                 updateOrigins <- originChoices
        #             } else {
        #                 updateOrigins <- c(originChoices,"Uploaded data")
        #             }
        #         } else {
        #             updateOrigins <- originChoices
        #         }
        #         ### Store updates to keep in memory previous state
        #         r$oldOrigins <- updateOrigins
        #         shinyWidgets::updateAwesomeRadio(
        #             session = session,
        #             inputId = "exclusiveDataChoice",
        #             choices = updateOrigins)
        #     }
        # }, ignoreNULL = FALSE)

        # Update origin according to value changing inside a list
        # data possibility :
        #       - NHU
        #       - Uploaded data
        #       - GEO dataset

        observe({
            supplumentaryDataUMAPNames <- c("Uploaded", "Toyset", "GSE")
            supplumentaryDataUMAPLabel <- list("Uploaded data", r$toysets$toysetName, r$gseName)
            names(supplumentaryDataUMAPLabel) <- supplumentaryDataUMAPNames

            filterSupplementaryData <- c(!is.null(r$leftPanel_UserInput), # Uploaded data
                                         !(is.null(r$toysets)), # Toyset
                                         !is.null(r$gseName)) # GSE dataset
            names(filterSupplementaryData) <- c("Uploaded", "Toyset", "GSE")
            newChoices <- c(originChoices, unname(unlist(supplumentaryDataUMAPLabel[filterSupplementaryData])))

            # Update data choices
            shinyWidgets::updateAwesomeRadio(
                session = session,
                inputId = "exclusiveDataChoice",
                choices = newChoices)


        })


        ##  Store selected data
        observeEvent(input$exclusiveDataChoice,{
            # # Check input value to assign correct value to reactive value
            # ## Put old value to not break the rest of code
            # if (input$exclusiveDataChoice == "Meta-cohort of MBIC tumors  (11 studies, n = 1530)"){
            #     r$heterogeneityUserSelection_exclusiveData <- "MIBC Tumor metacohort (1530)"
            # } else if (input$exclusiveDataChoice == "CCLE BLCA cell lines (n=36)"){
            #     r$heterogeneityUserSelection_exclusiveData <- "CCLE (36)"
            # } else {
            #     r$heterogeneityUserSelection_exclusiveData <- input$exclusiveDataChoice
            # }
            r$heterogeneityUserSelection_exclusiveData <- input$exclusiveDataChoice

        }, ignoreNULL = FALSE)
    })
}

## To be copied in the UI
# mod_heterogeneityUserSelectionOriginExclusiveChoice_ui("heterogeneityUserSelectionOriginExclusiveChoice_1")

## To be copied in the server
# mod_heterogeneityUserSelectionOriginExclusiveChoice_server("heterogeneityUserSelectionOriginExclusiveChoice_1", r=r)
