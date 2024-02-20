#' Function to output some classification
#'
#'
#'
#' @noRd
renderClassification <- function(RMOobject, classification, inputID){
    y <-  renderUI({
        #### It is possible to select more than one regulator to project influence.
        #### In this case the mean of influences will be displayed
        shinyWidgets::pickerInput(inputId = inputID,
            h4("3. Choose a subtype :"),
            choices = c("None",
                sort(unique(attr(RMOobject,classification)))
            ),
            selected = "None")

    })

    return(y)
}



#' heterogeneityUserSelectionInfluenceProjection UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @importFrom shiny NS tagList
#'
#' @noRd
mod_heterogeneityUserSelectionInfluenceProjection_ui <- function(id){
    ns <- NS(id)
    tagList(
        # Add pickerInput to select clustering subgroup
        uiOutput(ns("info"))
    )
}

#' heterogeneityUserSelectionInfluenceProjection Server Functions
#'
#' @noRd
mod_heterogeneityUserSelectionInfluenceProjection_server <- function(id, r){
    moduleServer( id, function(input, output, session){
        ns <- session$ns
        # Store RegMap variable
        mchRMO <- GBMRegMap::mchRMO
        ccleRMO <- GBMRegMap::ccleRMO
        stemDiffRMO <- GBMRegMap::stemDiffRMO
        u87mg <- GBMRegMap::u87RMO
        # Observe the user's selection on choice of clustering
        observe({
            whereami::cat_where( whereami::whereami())
            print("r$toysets In heterogeneityUserSelectionInfluenceProjection")
            print(r$toysets)
            print(ifelse(is.null(r$toysets), yes = "noToyset", no = r$toysets$toysetName))
            ### check toyset and gse names
            toysetsName <- ifelse(is.null(r$toysets), yes = "noToyset", no = r$toysets$toysetName)
            gseName <- ifelse(is.null(r$gseName), yes = "noGSE", no = r$gseName)

            # Store reactive variables
            dataOrigin <- r$heterogeneityUserSelection_exclusiveData
            projection <- r$heterogeneityUserSelection_clustersChoices
            print("heterogeneityUserSelectionInfluenceProjection - dataOrigin")
            print(dataOrigin)
            print("heterogeneityUserSelectionInfluenceProjection - projection")
            print(projection)

            # If variable are NULL render raw metacohort network with Consensus subtypes (Kamoun et al. 2020) selected
            if (is.null(projection)){
                r$heterogeneityUserSelection_clustersChoices <- "Verhaak"
                projection <- r$heterogeneityUserSelection_clustersChoices
            }
            #r$heterogeneityUserSelection_exclusiveData

            # First if level : dataOrigin
            ## - MIBC Tumor metacohort (1530)
            ## - CCLE (36)
            ## - NHU
            ## - Uploaded data
            # Second if level : projection
            ## - Consensus subtypes (Kamoun et al. 2020)
            ## - BLCA RegMap subtypes
            ## - Your classification
            if (dataOrigin == "Meta-cohort of GBM tumors  (16 studies, n=1612)"){
                if (projection == "Verhaak"){
                    output$info <- renderUI({
                        ### For classification choice we have Consensus subtypes (Kamoun et al. 2020), PAM clustering or user's classification
                        shinyWidgets::pickerInput(inputId = ns("subProjection"),
                            h4("3. Choose a subtype"),
                            choices = c("None",
                                sort(as.character(unique(mchRMO@classificationSota)))
                            ),
                            selected = "None")

                    })
                } else if (projection == "cRegMap.k7"){
                    # output$info <- renderClassification(RMOobject = mchRMO,
                    #     classification = "classificationInfluence1",
                    #     inputID = ns("subProjection"))
                    output$info <- renderUI({
                        shinyWidgets::pickerInput(inputId = ns("subProjection"),
                            h4("3. Choose a subtype"),
                            choices = c("None",
                                as.character(levels(mchRMO@classificationInfluence1))
                            ),
                            selected = "None")
                    })
                }
            } else if (dataOrigin == "CCLE GBM cell lines (n=42)"){
                # Change according to projection
                if (projection == "cRegMap.k7"){
                    output$info <- renderClassification(RMOobject = ccleRMO,
                        classification = "classificationInfluence1",
                        inputID = ns("subProjection"))
                } else if (projection == "Samples"){
                    output$info <- renderUI({
                        shinyWidgets::pickerInput(
                            inputId = ns("subProjection"),
                            label = h4("3. Choose a sample"),
                            choices = c("None",colnames(ccleRMO@influence)),
                            selected = "None",
                            options = list(
                                `live-search` = TRUE)
                        )
                    })
                }
            } else if (dataOrigin == toysetsName){
                if (projection == "cRegMap.k7"){
                    output$info <- renderUI({
                        shinyWidgets::pickerInput(
                            inputId = ns("subProjection"),
                            label = h4("3. Choose a subtype"),
                            choices = c("None",sort(unique(r$toysets$rmoObject@classificationInfluence1))),
                            selected = "None",
                            options = list(
                                `live-search` = TRUE)
                        )
                    })
                } else if (projection == "Samples"){
                    output$info <- renderUI({
                        shinyWidgets::pickerInput(
                            inputId = ns("subProjection"),
                            label = h4("3. Choose a sample"),
                            choices = c("None",colnames(r$toysets$rmoObject@influence)),
                            selected = "None",
                            options = list(
                                `live-search` = TRUE)
                        )
                    })
                }
            } else if (dataOrigin == gseName){
                if (projection == "cRegMap.k7"){
                    output$info <- renderUI({
                        shinyWidgets::pickerInput(
                            inputId = ns("subProjection"),
                            label = h4("3. Choose a subtype"),
                            choices = c("None",sort(unique(r$gseRMO@classificationInfluence1))),
                            selected = "None",
                            options = list(
                                `live-search` = TRUE)
                        )
                    })
                } else if (projection == "cRegMap.k12"){
                    output$info <- renderUI({
                        shinyWidgets::pickerInput(
                            inputId = ns("subProjection"),
                            label = h4("3. Choose a subtype"),
                            choices = c("None",sort(unique(r$gseRMO@classificationInfluence2))),
                            selected = "None",
                            options = list(
                                `live-search` = TRUE)
                        )
                    })
                } else if (projection == "Samples"){
                    output$info <- renderUI({
                        shinyWidgets::pickerInput(
                            inputId = ns("subProjection"),
                            label = h4("3. Choose a sample"),
                            choices = c("None",colnames(r$gseRMO@influence)),
                            selected = "None",
                            options = list(
                                `live-search` = TRUE)
                        )
                    })
                }
            } else if (dataOrigin == "Uploaded data"){
                # Change according to projection
                if (projection == "cRegMap.k7"){
                    output$info <- output$info <- renderUI({
                        shinyWidgets::pickerInput(
                            inputId = ns("subProjection"),
                            label = h4("3. Choose a sample"),
                            choices = c("None", sort(unique(r$leftPanel_UserInput_pamClassification1))),
                            selected = "None",
                            options = list(
                                `live-search` = TRUE)
                        )
                    })
                } else if (projection == "cRegMap.k12"){
                    output$info <- output$info <- renderUI({
                        shinyWidgets::pickerInput(
                            inputId = ns("subProjection"),
                            label = h4("3. Choose a sample"),
                            choices = c("None", sort(unique(r$leftPanel_UserInput_pamClassification2))),
                            selected = "None",
                            options = list(
                                `live-search` = TRUE)
                        )
                    })
                } else if (projection == "Samples"){
                    output$info <- renderUI({
                        shinyWidgets::pickerInput(
                            inputId = ns("subProjection"),
                            label = h4("3. Choose a sample"),
                            choices = c("None",colnames(r$leftPanel_UserInput)),
                            selected = "None",
                            options = list(
                                `live-search` = TRUE)
                        )
                    })
                }else if (projection == "Your classification"){
                    output$info <- renderUI({
                        shinyWidgets::pickerInput(
                            inputId = ns("subProjection"),
                            label = h4("3. Choose a sample"),
                            choices = c("None",sort(unique(r$leftPanel_UserClassification))),
                            selected = "None",
                            options = list(
                                `live-search` = TRUE)
                        )
                    })
                } else {
                    print("Setting cRegMap.k7")
                    projection <- "cRegMap.k7"
                    output$info <- renderUI({
                        shinyWidgets::pickerInput(
                            inputId = ns("subProjection"),
                            label = h4("3. Choose a subtype"),
                            choices = c("None",sort(unique(r$leftPanel_UserInput_pamClassification1))),
                            selected = "None",
                            options = list(
                                `live-search` = TRUE)
                        )
                    })
                }
            }
        })

        ## Store Sub-group
        observeEvent(input$subProjection,{
            r$heterogeneityUserSelection_subGroup <- input$subProjection
        }, ignoreNULL = FALSE)
    })
}

## To be copied in the UI
# mod_heterogeneityUserSelectionInfluenceProjection_ui("heterogeneityUserSelectionInfluenceProjection_ui_1")

## To be copied in the server
# mod_heterogeneityUserSelectionInfluenceProjection_server("heterogeneityUserSelectionInfluenceProjection_ui_1")
