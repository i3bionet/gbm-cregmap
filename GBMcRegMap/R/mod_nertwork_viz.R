# Script inspired by script developped by Aurelien Dispot, Julia Puig

#' Function to return min or max of metacohort
#'
#'
#' @noRd
infSubgroup <- function(mchRMO_Object,
    ccleRMO_Object,
    # nhuRMO_Object,
    dataOrigin,
    projection,
    fun = min,
    r,
    nodes){
    ## Store toyset name
    toysetsName <- ifelse(is.null(r$toysets), yes = "noToyset", no = r$toysets$toysetName)
    ## GSE name
    gseName <- ifelse(is.null(r$gseName), yes = "noGSE", no = r$gseName)
    # Set influence variable
    if (dataOrigin == "Meta-cohort of GBM tumors  (16 studies, n=1612)"){
        rmoObject<- mchRMO_Object
    } else if (dataOrigin == "CCLE GBM cell lines (n=42)"){
        rmoObject <- ccleRMO_Object
    } else if (dataOrigin == toysetsName){
        rmoObject <- r$toysets$rmoObject
    } else if (dataOrigin == gseName){
        rmoObject <- r$gseRMO
    } else if(dataOrigin == "Uploaded data"){
        tmpClasses1 <- as.numeric(r$leftPanel_UserInput_pamClassification1)
        names(tmpClasses1) <- names(r$leftPanel_UserInput_pamClassification1)
        tmpClasses2 <- as.numeric(r$leftPanel_UserInput_pamClassification2)
        names(tmpClasses2) <- names(r$leftPanel_UserInput_pamClassification2)
        rmoObject <- RegMapObject(influence = r$leftPanel_UserInput,
            UMAP2D = r$leftPanel_UserInput_UMAP2D,
            classificationInfluence1 = tmpClasses1,
            classificationInfluence2 = tmpClasses2)
    }
    # Restrict influence data to nodes in network
    rmoObject@influence <- rmoObject@influence[intersect(nodes$id, rownames(rmoObject@influence)),]
    # Compute min or max
    if (projection == "Verhaak"){
        y <- sapply(sort(unique(rmoObject@classificationSota)), function(x) {
            tmpInf <- rmoObject@influence[,names(rmoObject@classificationSota)[rmoObject@classificationSota == x]]
            y <- apply(tmpInf, 1, mean)
            y <- fun(y)
        })
    } else if (projection == "cRegMap.k7"){
        # Need to restore and rename classification (we had a unnaming problem with gseRMO)
        tmpClassif <- rmoObject@classificationInfluence1
        names(tmpClassif) <- colnames(rmoObject@influence)
        y <- sapply(sort(unique(tmpClassif)), function(x) {
            tmpInf <- rmoObject@influence[,names(tmpClassif)[tmpClassif == x]]
            if (is.null(dim(tmpInf))){
                y <- fun(tmpInf)
            } else {
                y <- apply(tmpInf, 1, mean)
                y <- fun(y)
            }
        })
    }  else if (projection == "Your classification"){
        # Need to restore and rename classification (we had a unnaming problem with gseRMO)
        tmpClassif <- r$leftPanel_UserClassification
        names(tmpClassif) <- colnames(rmoObject@influence)
        y <- sapply(sort(unique(tmpClassif)), function(x) {
            tmpInf <- rmoObject@influence[,names(tmpClassif)[tmpClassif == x]]
            if (is.null(dim(tmpInf))){
                y <- fun(tmpInf)
            } else {
                y <- apply(tmpInf, 1, mean)
                y <- fun(y)
            }
        })
    } else (
        stop(paste("Projection :", projection, "is not supported inside infSubgroup function. Please refer this message to contact."))
    )
    return(y)
}


#' nertwork_viz UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
#' @importFrom visNetwork visGetPositions
mod_nertwork_viz_ui <- function(id){

    ns <- NS(id)
    tagList(

        # Main panel
        #     - Title
        #     - Network Vizualisation + legend
        fluidRow(#h4("Network"),
            #actionButton("layout", "Layout"),,
            visNetwork::visNetworkOutput(ns("network"), height = "700"),
            br(),
            br(),
            br(),
            br(),
            br(),
            fluidRow(column(3,offset = 3,
                actionButton(ns("store_position"), "Store current network")
            ),
                column(3, offset = 2,
                    downloadButton(ns("downloadNetwork"), "Download stored network as HTML")
                )))


    )
}

#' nertwork_viz Server Functions
#'
#'
#' OUTPUTS :
#' - output$maxGRN : maximum of shared GRN inside specified GRN.
#' - output$clusterChoices : list of cluster choices for classification results.
#'
#'
#' @importFrom visNetwork visNetworkOutput renderVisNetwork visUpdateEdges visUpdateNodes
#' @importFrom RColorBrewer  brewer.pal
#' @noRd
mod_nertwork_viz_server <- function(id, r){
    moduleServer( id, function(input, output, session){
        ns <- session$ns
        # Store RegMap variable
        mchRMO <- GBMRegMap::mchRMO
        ccleRMO <- GBMRegMap::ccleRMO
        nodesPositions <- GBMRegMap::nodesPositions

        # stemDiffRMO <- GBMRegMap::stemDiffRMO
        # u87mg <- GBMRegMap::u87RMO


        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # Load node and edge data
        networkData = computeNetworkEdgesNodes(influence = GBMRegMap::mchRMO@influence, grn = GBMRegMap::refinedGRN)
        # -----------------------------
        # /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\
        # For now it is a static choice but in the
        # future if the app use other this input will
        # be dynamicaly updated.
        # /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\
        data.edges <- networkData$networkEdges
        # data about edges :
        # id : <string> id of the edge
        #   'p1', 'p2', ... : interaction predicted by CoRegNet
        #   'e1', 'e2', ... : evidence innteration
        #       some have arrow "to" (cooperative evidence)
        # from : <string> node id
        # to : <string> node id
        # nGRN : <number>
        #   number of GRN share by 2 nodes
        #   of infinite Inf (for edges that were not predicted by hLICORN)
        # title : <string> popup when the mouse over the edge
        # value : <number> scale edge width
        #   range from 0 to 1 : ]0; 1]
        #   depend of nGRN value
        # arrow : <string>
        #   'none' if no arrow
        #   'to'
        # hidden : <boolean>
        #   TRUE : edge is not drawn but still par of the physic simulation
        # color : <string>
        #   'grey' for predicted edges
        #   other color for evidence edges

        # /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\
        # For now it is a static choice but in the
        # future if the app use other this input will
        # be dynamicaly updated.
        # /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\
        data.nodes <- networkData$networkNodes
        # data about nodes :
        # id : <string> id of the node == label
        # label : <string> label of the node that is drawn on the network == id
        # shape : <string>
        #   'circle'
        # title : <string> popup when the mouse over the node
        # value : <number> scale node width
        #   range from 0 to 1 : ]0; 1]
        #   depends on how many genes are regulated by this node
        # hidden : <boolean>
        #   TRUE : node not drawn but still part of the physic simulation
        # color : <string>
        #   from blue - white - red : depends on the influence value
        #   grey : not influence calculated

        # data.TFA <- read.table("limitsTFA.data", sep = "\t", stringsAsFactors = FALSE, header = FALSE)
        # minTFA <- data.TFA[1,1]
        # maxTFA <- data.TFA[2,1]

        # For the network legend
        # -----------------------------
        # nodes data.frame for legend
        lnodes <- data.frame(
            label = c("Low influence", "High influence"),
            shape = c("circle", "circle"),
            color = c("blue", "red"),
            title = "Informations",
            id = 1:2
        )
        # edges data.frame for legend
        ledges <- data.frame(
            color = c("red", "blue"),
            label = c("Regulatory_evidence", "Coregulatory_evidence"),
            arrows = c("to","none")
        )
        # For the minimum and maximum influecne value by classes (usefull below)
        # -----------------------------

        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # observe Events
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        # selectInput - clsutersChoices : selecting clustering results
        # -----------------------------
        # clustersChoices.r <- reactive({input$clustersChoices})
        #
        # observeEvent(clustersChoices.r(),{
        #     # Store cluster selected
        #     selected.cluster <- input$clustersChoices
        #     ## Get clustering results according to cluster choices
        #     clusteringResults <- metacohortData$clustering[[selected.cluster]]
        #     ## Get clusters names plus "None" option
        #     clustersList <- clusteringPossibilities(clusteringResults)
        #     ## Add 'None' option
        #     clustersList <- c(clustersList, "None")
        #     ## Update selectInput
        #     updateSelectInput(inputId = "info", choices = clustersList)
        # }, ignoreNULL = FALSE)


        # To store reactive values
        # -----------------------------
        reactives <- reactiveValues()

        # The following part will re-run whenever info or mingrn change
        # --------------------------------------------------------------------------
        observe({
            whereami::cat_where( whereami::whereami() )
            # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            # Store reactiveValues
            ## Data
            dataOrigin <- r$heterogeneityUserSelection_exclusiveData
            ### If dataOrigin is null ()
            ## Projection / Clustering choice
            projection <- r$heterogeneityUserSelection_clustersChoices
            ## Sub-projection
            sub_projection <- r$heterogeneityUserSelection_subGroup
            ## Store toyset name
            toysetsName <- ifelse(is.null(r$toysets), yes = "noToyset", no = r$toysets$toysetName)
            ## GSE name
            gseName <- r$gseName
            gseName <- ifelse(is.null(gseName), yes = "noGSE", no = gseName)
            ## Sub-group selection
            ## Check if sub_projection is null.
            ## At the starting of the app fix subgroup to 'None'
            if (is.null(sub_projection)){
                sub_projection <- "None"
            }
            ## Check inputs

            ## If gseName has been selected and user change GSE dataset we
            ## need to check the dataOrigin and gseName
            if (grepl(pattern = "GSE", x = dataOrigin)){
                if (dataOrigin != gseName){
                    dataOrigin <- "Meta-cohort of GBM tumors  (16 studies, n=1612)"
                }
            }
            ## Create a list for correct projection for each dataOrigin
            ### check toyset and gse names
            toysetsName <- ifelse(is.null(r$toysets), yes = "noToyset", no = r$toysets$toysetName)
            safeguardProjectionList <- list(c("Verhaak",
                "cRegMap.k7"),
                c("cRegMap.k7",
                    "Samples"),
                c("cRegMap.k7",
                    "Samples",
                    "Your classification"),
                c("cRegMap.k7",
                    "cRegMap.k12",
                    "Samples"),
                c("cRegMap.k7",
                    "cRegMap.k12",
                    "Samples"))
            names(safeguardProjectionList) <- c("Meta-cohort of GBM tumors  (16 studies, n=1612)",
                "CCLE GBM cell lines (n=42)",
                "Uploaded data",
                toysetsName,
                gseName)


            ### Check and correct projection
            if (!(projection %in% safeguardProjectionList[[dataOrigin]])){
                projection <- safeguardProjectionList[[dataOrigin]][1]
                sub_projection <- "None"
            }
            # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            # Select correct influence data according to user's choices
            ## Save regulators influence as TFA
            if(dataOrigin == "Meta-cohort of GBM tumors  (16 studies, n=1612)"){
                req(dataOrigin == "Meta-cohort of GBM tumors  (16 studies, n=1612)")
                TFA <- mchRMO@influence
            } else if (dataOrigin == "CCLE GBM cell lines (n=42)"){
                TFA <- ccleRMO@influence
            } else if (dataOrigin == toysetsName){
                TFA <- r$toysets$rmoObject@influence
            } else if (dataOrigin == gseName){
                TFA <- r$gseRMO@influence
            } else if (dataOrigin == "Uploaded data"){
                TFA <- r$leftPanel_UserInput
            }

            # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            #  Select correct classification according to user's choices

            if (dataOrigin == "Meta-cohort of GBM tumors  (16 studies, n=1612)"){
                if (projection == "Verhaak"){
                    clinicalData <- mchRMO@classificationSota
                } else if (projection == "cRegMap.k7"){
                    clinicalData <- mchRMO@classificationInfluence1
                } else {
                    clinicalData <- NA
                }
            } else if (dataOrigin == "CCLE GBM cell lines (n=42)"){
                if (projection == "cRegMap.k7"){
                    clinicalData <- ccleRMO@classificationInfluence1
                } else if (projection == "Samples"){
                    clinicalData <- colnames(ccleRMO@influence)
                    names(clinicalData) <- colnames(ccleRMO@influence)
                } else {
                    clinicalData <- NA
                }
            } else if (dataOrigin == toysetsName){
                if (projection == "cRegMap.k7"){
                    clinicalData <-r$toysets$rmoObject@classificationInfluence1
                } else if (projection == "cRegMap.k12"){
                    clinicalData <- r$toysets$rmoObject@classificationInfluence2
                } else if( projection == "Samples"){
                    clinicalData <- colnames(r$toysets$rmoObject@influence)
                    names(clinicalData) <- colnames(r$toysets$rmoObject@influence)
                }
            } else if (dataOrigin == gseName){
                if (projection == "cRegMap.k7"){
                    clinicalData <- r$gseRMO@classificationInfluence1
                    # We need to rename classification (Don't know why is has been unnamed)
                    names(clinicalData) <- colnames(r$gseRMO@influence)
                } else if (projection == "cRegMap.k12"){
                    clinicalData <- r$gseRMO@classificationInfluence2
                    # We need to rename classification (Don't know why is has been unnamed)
                    names(clinicalData) <- colnames(r$gseRMO@influence)
                } else if( projection == "Samples"){
                    clinicalData <- colnames(r$gseRMO@influence)
                    names(clinicalData) <- colnames(r$gseRMO@influence)
                }
            } else if (dataOrigin == "Uploaded data"){
                if (projection == "cRegMap.k7"){
                    clinicalData <- r$userRMO@classificationInfluence1
                    # We need to rename classification (Don't know why is has been unnamed)
                    names(clinicalData) <- colnames(r$userRMO@influence)
                } else if (projection == "cRegMap.k12"){
                    clinicalData <- r$userRMO@classificationInfluence2
                    # We need to rename classification (Don't know why is has been unnamed)
                    names(clinicalData) <- colnames(r$userRMO@influence)
                } else if (projection == "Samples"){
                    clinicalData <- colnames(r$userRMO@influence)
                    names(clinicalData) <- colnames(r$userRMO@influence)
                } else if (projection == "Your classification"){
                    clinicalData <- r$leftPanel_UserClassification
                }
            } else {
                stop("Uncorrect clustering choice parameter")
            }
            # which edges to take into account
            data.edgesOK <- data.edges
            # which nodes to take into account (i.e. nodes that are present in data.edgesOK)
            nodes.with_sufficient_nGRN.vector = unique(unlist(data.edgesOK[c("from", "to")]))
            data.nodesOK <- data.nodes[which(data.nodes$id %in% nodes.with_sufficient_nGRN.vector), ]

            # node coloring
            # nodes can be colored: (1) all in None, (2) all with a palette, (3) according to clinical data
            if (sub_projection == "None") {
                reactives$TFA <- TFA
                reactives$clinicalData <- clinicalData
            } else {
                selected_samples <- names(clinicalData)[clinicalData == sub_projection & !is.na(clinicalData)]
                reactives$TFA <- TFA[,selected_samples]
                reactives$clinicalData <- clinicalData[selected_samples]
            }
            # take mean influence of each node, and 0 if no influence
            reacTFA <- as.matrix(reactives$TFA)
            idxGrey <- NULL
            ## Update TFA value
            data.nodesOK$TFA <- sapply(
                data.nodesOK$id,
                function(id){
                    if ((id %in% rownames(reacTFA)) & ncol(TFA) > 1) {
                        return(mean(reacTFA[id, ]))
                    } else if (id %in% rownames(reacTFA)) {
                        return(reacTFA[id, ])
                    } else {
                        return(-100)
                    }
                }
            )
            idxGrey <- which(data.nodesOK$TFA == -100)
            data.nodesOK$TFA[idxGrey] <- 0
            # assign node color with the palette function
            if (sub_projection == "None") {
                data.nodesOK$color <- rep("lightblue", length(data.nodesOK$id))
            } else {
                # version without rescaling
                #reactives$breaks <- seq(min(data.nodesOK$TFA), max(data.nodesOK$TFA), by = (max(data.nodesOK$TFA) - min(data.nodesOK$TFA)) / 75)
                #pal <- leaflet::colorNumeric(colorRampPalette(c("blue","honeydew","red"))(75), domain = reactives$breaks, reverse = FALSE)
                #data.nodesOK$color <- pal(data.nodesOK$TFA)

                # Comput max and min influence for rescaling
                # maxTFA <- max(data.nodesOK$TFA) + 1e-5
                # minTFA <- min(data.nodesOK$TFA) - 1e-5


                # If user want to projection "Samples" or "Samples" influence
                ## Possible use of case :
                ##              - CCLE (36)
                ##              - r$toysets
                ##              - r$gseName
                ##              - Uploaded data

                if (projection %in% c("Samples")) {
                    # maxTFA <- max(reacTFA)
                    # minTFA <- min(reacTFA)
                    maxTFA <- max(data.nodesOK$TFA)
                    minTFA <- min(data.nodesOK$TFA)
                } else if ((dataOrigin == "Uploaded data") & (projection == "Samples")) {
                    maxTFA <- max(r$userRMO@influence)
                    minTFA <- min(r$userRMO@influence)
                } else if ((dataOrigin == toysetsName) & (projection == "Samples")){
                    maxTFA <- max(r$toysets$rmoObject@influence)
                    minTFA <- min(r$toysets$rmoObject@influence)
                } else if ((dataOrigin == gseName) & (projection == "Samples")){
                    maxTFA <- max(r$gseRMO@influence)
                    minTFA <- min(r$gseRMO@influence)
                } else {
                    # Take min and max of average influence per groups
                    # minInfGroups <- sapply(sort(unique(mchRMO@classificationSota)), function(x) {
                    #     tmpInf <- mchRMO@influence[,names(mchRMO@classificationSota)[mchRMO@classificationSota == x]]
                    #     y <- apply(tmpInf, 1, mean)
                    #     return(min(y))
                    # })
                    minInfGroups <- try({
                        infSubgroup(mchRMO_Object = mchRMO,
                            ccleRMO_Object = ccleRMO,
                            # nhuRMO_Object = nhuRMO,
                            dataOrigin = dataOrigin,
                            projection = projection,
                            fun = min,
                            r=r, nodes = data.nodesOK)
                    })
                    # maxInfGroups <- sapply(sort(unique(mchRMO@classificationSota)), function(x) {
                    #     tmpInf <- mchRMO@influence[,names(mchRMO@classificationSota)[mchRMO@classificationSota == x]]
                    #     y <- apply(tmpInf, 1, mean)
                    #     return(max(y))
                    # })
                    maxInfGroups <- try({
                        infSubgroup(mchRMO_Object = mchRMO,
                            ccleRMO_Object = ccleRMO,
                            # nhuRMO_Object = nhuRMO,
                            dataOrigin = dataOrigin,
                            projection = projection, fun = max, r=r, nodes = data.nodesOK)
                    })

                    maxTFA <- max(maxInfGroups)
                    minTFA <- min(minInfGroups)
                }

                # version with rescaling
                reactives$breaks <- seq(minTFA, maxTFA, by = (maxTFA - minTFA) / 75)
                pal <- leaflet::colorNumeric(colorRampPalette(c("blue","honeydew","red"))(75),
                    domain = reactives$breaks, reverse = FALSE)
                data.nodesOK$color <- pal(data.nodesOK$TFA)
                data.nodesOK$color[idxGrey] <- "grey"

            }
            # update nodes and edges
            reactives$data.nodes <- data.nodesOK
            reactives$data.edges <- data.edgesOK

            ## Store data.nodes and data.edges inside r
            #r$networkViz_data.nodes <- reactives$data.nodes
            #r$networkViz_data.edges <- reactives$data.edges
            # update network according to node and edge changes
            try({visNetwork::visNetworkProxy(ns("network"))  %>%
                    visNetwork::visRemoveNodes(id = setdiff(data.nodes$id, data.nodesOK$id)) %>%
                    visNetwork::visRemoveEdges(id = setdiff(data.edges$id, data.edgesOK$id)) %>%
                    visNetwork::visUpdateNodes(nodes = data.nodesOK) %>%
                    visNetwork::visUpdateEdges(edges = data.edgesOK) %>%
                    visNetwork::visOptions(
                        # highlight nearest when clicking a node
                        highlightNearest = list(enabled = TRUE, hover = FALSE),
                        # create an html select element
                        nodesIdSelection = list(enabled = TRUE, values = sort(data.nodes$id)))})

            # prepare ComplexHeatmap annotation
            # ha_subtype = HeatmapAnnotation(
            #     df = data.frame("subtype" = reactives$clinicalData),
            #     col = list(subtype = col_subtype)
            # )
            # reactives$heatmap = ComplexHeatmap::Heatmap(
            #     reactives$TFA,
            #     top_annotation = ha_subtype,
            #     name = "influence", row_title = "genes", column_title = "samples"
            # )
        }, label = "Network module")




        # Nodes positions for download
        # -----------------------------
        # get position info
        observeEvent(input$store_position, {
            # Carefull we need to add 'ns()' to communicate with ui part.
            visNetwork::visNetworkProxy(ns("network")) %>%
                visNetwork::visGetPositions()
        })
        # format positions
        nodes_positions <- reactive({
            positions <- input$network_positions
            if(!is.null(positions)){
                nodes_positions <- do.call("rbind", lapply(positions, function(x){ data.frame(x = x$x, y = x$y)}))
                nodes_positions$id <- names(positions)
                nodes_positions
            } else {
                NULL
            }
        })



        # # For the heatmap legend
        # # -----------------------------
        # # Adapt the palette
        # palette(RColorBrewer::brewer.pal(length (levels(clinicalData)),"Set3"))
        # palette()
        # # Prepare a named vector with a color attached to each subtype
        # col_subtype <- palette()[as.factor( levels(clinicalData) ) ]
        # names(col_subtype) <- levels(clinicalData)
        # end of observe
        # --------------------------------------------------------------------------




        # Nodes positions for download
        # -----------------------------
        # # format positions
        # nodes_positions <- reactive({
        #     positions <- input$network_positions
        #     if(!is.null(positions)){
        #         nodes_positions <- do.call("rbind", lapply(positions, function(x){ data.frame(x = x$x, y = x$y)}))
        #         nodes_positions$id <- names(positions)
        #         nodes_positions
        #     } else {
        #         NULL
        #     }
        # })

        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # OUTPUTS
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # Network
        # -----------------------------

        # # DEBUG : save nodes and edges to plot article figures
        # if (exists("data.nodes")){
        #     networkElements <- list(data.nodes, data.edges)
        #     names(networkElements) <- c("nodes","edges")
        #     save(networkElements, file = "../../../articles/2023/neuroOncology/scripts/inputs/gbmNetwork/nodesEdges.RData")
        #     print("network saved !")
        # }

        # Increase label size
        # data.nodes$font.size <- 90


        # --------------------------- #
        # get stabilized nodes positions
        coords <- nodesPositions[data.nodes$id,c("x","y")] %>% as.matrix()
        data.nodes <- cbind(data.nodes, coords)

        output$network <- visNetwork::renderVisNetwork({
            visNetwork::visNetwork(data.nodes, data.edges, main = "GBM Co-Regulatory Network") %>%
                visNetwork::visPhysics(
                    # whether to wait for node stabilisation for visualisation
                    stabilization = FALSE,
                    # method
                    solver = "barnesHut",
                    # gravitationalConstant: gravitation power
                    # centralGravity: gravitation power towards the network center (0,1)
                    # avoidOverlap: avoid nodes to overlap with each other (0,1)
                    barnesHut = list(gravitationalConstant = -5000, centralGravity = 0.1, avoidOverlap = 0.1)) %>% #-3000
                visNetwork::visEdges(
                    # smooth edges
                    smooth = list(enabled = TRUE, type = "dynamic")) %>%
                # visNetwork::visLegend(width = 0.1, position = "right", main = "Network legend",
                #           addEdges = ledges, addNodes = lnodes, useGroups = FALSE) %>%
                visNetwork::visOptions(
                    # highlight nearest when clicking a node
                    highlightNearest = list(enabled = TRUE, hover = FALSE),
                    # create an html select element
                    nodesIdSelection = list(enabled = TRUE, values = sort(data.nodes$id))) %>%
                visNetwork::visLayout(randomSeed = 123) %>%
                visNetwork::visInteraction(navigationButtons = TRUE) %>%
                visNetwork::visEvents(startStabilizing = "function() {
            this.moveTo({scale:0.15})}", type = "once")
            # %>%
            #     visNetwork::visIgraphLayout(layout = "layout.norm", layoutMatrix = coords)
        })


        # Select a node in the network
        # -----------------------------
        output$network_selected = renderPrint({
            input$network_selected
        })

        # if (length(input$network_selected) == 1){
        #   output$connectedNodes <- renderPrint({
        #     visNetworkProxy("network") %>%
        #       #visNearestNodes(target = input$network_selected)$id
        #       visGetConnectedNodes(id = input$network_selected, input = "connected_nodes")
        #   })
        # }

        # TO DO: print neighbors when we select a node
        # output$connectedNodes = renderPrint({
        #   input$connectedNodes
        # })
        # observe({
        #     if (length(input$network_selected) == 1){
        #         visNetworkProxy("network") %>%
        #             visGetConnectedNodes(id = input$network_selected, input = "connected_nodes")
        #     }
        # })



        # Download current network
        # -----------------------------
        output$downloadNetwork <- downloadHandler(

            # download a still moving network
            # filename = "network.html",
            # content = function(con) {
            #   visNetwork(data.nodes, data.edges, main = "Co-Regulatory network") %>%
            #     visPhysics(
            #       stabilization = FALSE,      # (better for performance) : do not wait for node stabilization to start showing the network
            #       solver = "barnesHut", # either 'barnesHut' 'repulsion', 'hierarchicalRepulsion' (if hierachical layout), 'forceAtlas2Based' or your own solver .
            #       barnesHut = list(gravitationalConstant = -4000, centralGravity = 0.1, avoidOverlap = 0.1) # avoid nodes to overlap with each other (from 0 to 1)
            #     ) %>%
            #     visEdges(smooth = list(enabled = TRUE, type = "dynamic")) %>% visSave(con)
            # }

            # download a static but modifiable network
            filename = paste("object.html", sep=""),
            content = function(con) {
                # store nodes current positions
                nodes_positions <- nodes_positions()
                # nodesPositions <-
                # save(nodesPositions, file = "data/nodesPositions.rda")
                if(!is.null(nodes_positions)){
                    nodes_save <- merge(reactives$data.nodes, nodes_positions, by = "id", all = T)
                } else  {
                    nodes_save <- reactives$data.nodes
                }
                # network to download
                visNetwork::visNetwork(nodes_save, reactives$data.edges, width = "1300", height = "1000") %>%
                    visNetwork::visEdges(smooth = list(enabled = TRUE)) %>%
                    visNetwork::visOptions(highlightNearest = list(enabled = TRUE, hover = FALSE),
                        nodesIdSelection = list(enabled = TRUE, values = sort(data.nodes$id))) %>%
                    visNetwork::visExport() %>%
                    visNetwork::visPhysics(enabled = FALSE) %>%
                    visNetwork::visSave(con)
            }
        )

    })
}

## To be copied in the UI
# mod_nertwork_viz_ui("nertwork_viz_ui_1")

## To be copied in the server
# mod_nertwork_viz_server("nertwork_viz_ui_1")
