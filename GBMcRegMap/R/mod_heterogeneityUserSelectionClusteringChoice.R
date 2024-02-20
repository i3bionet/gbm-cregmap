#' heterogeneityUserSelectionClusteringChoice UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
#' @importFrom shinyWidgets switchInput
mod_heterogeneityUserSelectionClusteringChoice_ui <- function(id){
    ns <- NS(id)
    tagList(
        # uiOutput(ns("clustersChoices"))
        # Reder desired UI
        uiOutput(ns("CoRegNetSubtypesSelection")),


        # conditionalPanel(condition = "output.chosenData == 'noCCLE'",
        #                  ns = ns,
        #                  shinyWidgets::awesomeRadio(ns("clustersChoices"),
        #                                             h4("2. Choose annotation"),
        #                                             # /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\
        #                                             # For now it is a static choice but in the x
        #                                             # future if the app use other this input will
        #                                             # be dynamicaly updated.
        #                                             # /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\
        #                                             choices = c("Consensus subtypes (Kamoun et al. 2020)",
        #                                                         "BLCA RegMap subtypes"),
        #                                             status = "default")
        # ),
        # conditionalPanel(condition = "output.chosenData == 'CCLEselected'",
        #                  ns = ns,
        #                  fluidRow()
        # )

    )
}

#' heterogeneityUserSelectionClusteringChoice Server Functions
#'
#' @noRd
mod_heterogeneityUserSelectionClusteringChoice_server <- function(id, r){
    moduleServer( id, function(input, output, session){
        ns <- session$ns
        observeEvent(r$heterogeneityUserSelection_exclusiveData,{
            ### check toyset and gse names
            toysetsName <- ifelse(is.null(r$toysets), yes = "noToyset", no = r$toysets$toysetName)
            gseName <- ifelse(is.null(r$gseName), yes = "noGSE", no = r$gseName)
            dataOrigin <-r$heterogeneityUserSelection_exclusiveData
            # chosenData <- gsub(r$heterogeneityUserSelection_exclusiveData, pattern = " ", replacement = "_")
            # print(chosenData)
            # chosenData <- gsub(chosenData,pattern = "[()]", replacement = "_")


            if (dataOrigin== "Meta-cohort of GBM tumors  (16 studies, n=1612)"){
                output$CoRegNetSubtypesSelection <- renderUI({
                    ### For classification choice we have Molecular subtypes, PAM clustering or user's classification
                    shinyWidgets::awesomeRadio(ns("clustersChoices"),
                        h4("2. Choose annotation"),
                        # /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\
                        # For now it is a static choice but in the x
                        # future if the app use other this input will
                        # be dynamicaly updated.
                        # /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\
                        choices = c("cRegMap.k7",
                            "Verhaak"
                        ),
                        status = "default")
                })
            } else if (dataOrigin== "CCLE GBM cell lines (n=42)"){
                output$CoRegNetSubtypesSelection <- renderUI({
                    ### For classification choice we have Molecular subtypes, PAM clustering or user's classification
                    shinyWidgets::awesomeRadio(ns("clustersChoices"),
                        h4("2. Choose annotation"),
                        # /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\
                        # For now it is a static choice but in the x
                        # future if the app use other this input will
                        # be dynamicaly updated.
                        # /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\
                        choices = c("cRegMap.k7",
                            "No classification, per sample"),
                        status = "default")
                })
            } else if (dataOrigin == "Uploaded data"){
                if (is.null(r$leftPanel_UserClassification)){
                    output$CoRegNetSubtypesSelection <- renderUI({
                        ### For classification choice we have Molecular subtypes, PAM clustering or user's classification
                        shinyWidgets::awesomeRadio(ns("clustersChoices"),
                            h4("2. Choose annotation"),
                            # /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\
                            # For now it is a static choice but in the x
                            # future if the app use other this input will
                            # be dynamicaly updated.
                            # /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\
                            choices = c("cRegMap.k7",
                                # "cRegMap.k12",
                                "No classification, per sample"),
                            status = "default")
                    })
                } else {
                    output$CoRegNetSubtypesSelection <- renderUI({
                        ### For classification choice we have Molecular subtypes, PAM clustering or user's classification
                        shinyWidgets::awesomeRadio(ns("clustersChoices"),
                            h4("2. Choose annotation"),
                            # /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\
                            # For now it is a static choice but in the x
                            # future if the app use other this input will
                            # be dynamicaly updated.
                            # /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\
                            choices = c("cRegMap.k7",
                                # "cRegMap.k12",
                                "No classification, per sample",
                                "Your classification"),
                            status = "default")
                    })
                }
            } else if (dataOrigin == toysetsName){
                output$CoRegNetSubtypesSelection <- renderUI({
                    ### For classification choice we have Molecular subtypes, PAM clustering or user's classification
                    shinyWidgets::awesomeRadio(ns("clustersChoices"),
                        h4("2. Choose annotation"),
                        # /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\
                        # For now it is a static choice but in the x
                        # future if the app use other this input will
                        # be dynamicaly updated.
                        # /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\
                        choices = c("cRegMap.k7",
                            # "cRegMap.k12",
                            "No classification, per sample"),
                        status = "default")
                })
            } else if (dataOrigin == gseName){
                output$CoRegNetSubtypesSelection <- renderUI({
                    ### For classification choice we have Molecular subtypes, PAM clustering or user's classification
                    shinyWidgets::awesomeRadio(ns("clustersChoices"),
                        h4("2. Choose annotation"),
                        # /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\
                        # For now it is a static choice but in the x
                        # future if the app use other this input will
                        # be dynamicaly updated.
                        # /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\
                        choices = c("cRegMap.k7",
                            # "cRegMap.k12",
                            "No classification, per sample"),
                        status = "default")
                })
            }

            #output$chosenData <-  renderText({chosenData})





            # If user's data is selected without user's classification
            # else if ((!is.null(r$leftPanel_UserInput)) & is.null(r$leftPanel_UserClassification)) {
            #     shinyWidgets::updateAwesomeRadio(session = session,
            #                                      inputId = "clustersChoices",
            #                                      choices = c("BLCA RegMap subtypes"),
            #                                      status = "default",
            #                                      selected = "Pam clustexring")
            #     # If user updated his data and classification
            # } else if (!is.null(r$leftPanel_UserInput) & !is.null(r$leftPanel_UserClassification))  {
            #     shinyWidgets::updateAwesomeRadio(session = session,
            #                                      inputId = "clustersChoices",
            #                                      choices = c("BLCA RegMap subtypes","Your classification"),
            #                                      status = "default",
            #                                      selected = "BLCA RegMap subtypes")
            # }
        })

        ## Cluster choice
        observeEvent(input$clustersChoices,{
            # Change clustersChoice reactive Value
            # if (input$clustersChoices == "Consensus") {
            #     r$heterogeneityUserSelection_clustersChoices <- "Consensus subtypes (Kamoun et al. 2020)"
            # } else if (input$clustersChoices == "BLCA-cRegMap") {
            #     r$heterogeneityUserSelection_clustersChoices <- "BLCA RegMap subtypes"
            # } else
            if (input$clustersChoices == "No classification, per sample") {
                r$heterogeneityUserSelection_clustersChoices <- "Samples"
            }
            else {
                r$heterogeneityUserSelection_clustersChoices <- input$clustersChoices
            }
            # Reset sub-group projection reactive value
            r$heterogeneityUserSelection_subGroup <- NULL
        })
    })
}

## To be copied in the UI
# mod_heterogeneityUserSelectionClusteringChoice_ui("heterogeneityUserSelectionClusteringChoice_ui_1")

## To be copied in the server
# mod_heterogeneityUserSelectionClusteringChoice_server("heterogeneityUserSelectionClusteringChoice_ui_1", r=r)
