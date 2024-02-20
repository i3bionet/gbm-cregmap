#' A function to return maximum number of designed GRN
#' @importFrom CoRegNet coregulators
#' @param grn GRN used for projection
#' @noRd
numberGRN_params <- function(grn){
    # Extract coregulators from GRN
    coregs <- grn@coRegulators
    # Get maximum of shared GRN among coregulators
    maxGRN <- max(coregs$nGRN)
    return(maxGRN)
}


#' Function to compute network's nodes and edges
#' ## Inspired by code developped by Aurelien Dispot, Julia Puig
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# INITIATE VALUES
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#'
#'
#' @noRd
computeNetworkEdgesNodes <- function(influence, grn){
    ## code to comoute network's edges and nodes
    # TFA = metacohortData$Influence
    TFA <- influence
    nmsSquare=NULL
    nmsTriangle=NULL
    nmsDmd=NULL
    keepNodes=TRUE
    minGRNdefault=5
    alphaThresh=0.05

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # COMPUTATION
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # refined and completed GRN
    coregnetwork <- grn

    # network's edges and nodes
    # coregulators
    # ----------------
    coregs = coregnetwork@coRegulators[coregnetwork@coRegulators$adjustedPvalue < 0.05,]
    coregs.with_sufficient_nGRN = coregs[which(coregs$nGRN >= minGRNdefault), ]

    # SUBNET (keep exclusively regs with an influence)
    idx <- unname(apply(coregs.with_sufficient_nGRN, 1, function(x) isTRUE(x[1] %in% rownames(TFA) & x[2] %in% rownames(TFA)) ))
    coregs.with_sufficient_nGRN <- coregs.with_sufficient_nGRN[which(idx == TRUE),]

    colnames(coregs.with_sufficient_nGRN)[ colnames(coregs.with_sufficient_nGRN)  == "Reg1" ]  = "from"
    colnames(coregs.with_sufficient_nGRN)[ colnames(coregs.with_sufficient_nGRN)  == "Reg2" ]  = "to"
    coregs.with_sufficient_nGRN = data.frame(
        "id" = paste0("p", 1:nrow(coregs.with_sufficient_nGRN)),
        "title" = paste(coregs.with_sufficient_nGRN$from, coregs.with_sufficient_nGRN$to),
        "from" = coregs.with_sufficient_nGRN$from,
        "to" = coregs.with_sufficient_nGRN$to,
        "arrows" = rep("none", nrow(coregs.with_sufficient_nGRN)),
        "color" = rep("#bdbdbd", nrow(coregs.with_sufficient_nGRN)),
        "nGRN" = coregs.with_sufficient_nGRN$nGRN,
        "value" = coregs.with_sufficient_nGRN$nGRN / max(coregs.with_sufficient_nGRN$nGRN),
        "hidden" = rep(FALSE, nrow(coregs.with_sufficient_nGRN)),
        "physics" = rep(TRUE, nrow(coregs.with_sufficient_nGRN)),
        stringsAsFactors = FALSE
    )
    coregs.with_sufficient_nGRN.vector = unique(unlist(coregs.with_sufficient_nGRN[c("from", "to")]))

    # evidences
    # ----------------
    if (length(coregnetwork@evidences) > 0) {
        coregnetwork@evidences
        palette.evidence <- RColorBrewer::brewer.pal(length(coregnetwork@evidences), "Set1")
        names(palette.evidence) <- names(coregnetwork@evidences)
        additionalInts <- do.call(rbind,
                                  lapply(names(coregnetwork@evidences),
                                         function(evname) {
                                             ev = coregnetwork@evidences[[evname]]
                                             colnames(ev) = c("from", "to")

                                             #-----
                                             # keep only the edges that was already found by hlicorn :

                                             # Prefilter using evidences edges with 2 diffrent nodes in coregnetwork
                                             to_be_add = ev[
                                                 which(  ev[, 1] %in% coregs.with_sufficient_nGRN.vector # node in the coregnetwork
                                                         & ev[, 2] %in% coregs.with_sufficient_nGRN.vector # node in the coregnetwork
                                                         & ev[, 1] != ev[, 2] # evidence between 2 different nodes
                                                 ),
                                             ]
                                             all( c(to_be_add$from, to_be_add$to) %in% coregs.with_sufficient_nGRN.vector)
                                             # Keep only the evidence edges predicted by CoRegNet :
                                             to_be_add = to_be_add[ apply(
                                                 to_be_add,
                                                 1,
                                                 function(r){
                                                     w =  which(   (r["from"] == coregs.with_sufficient_nGRN$from | r["from"] == coregs.with_sufficient_nGRN$to)
                                                                   & (r["to"] == coregs.with_sufficient_nGRN$from | r["to"] == coregs.with_sufficient_nGRN$to)
                                                     )
                                                     #if ((r["from"] == "STAT4")&(r["to"] == "IRF4")){ browser()}
                                                     return(length(w) != 0) # TRUE -> will be kept
                                                 }
                                             ), ]

                                             #-----
                                             # Add informations for visualisation :
                                             to_be_add = data.frame(
                                                 "id" = paste0("e-",evname, 1:nrow(to_be_add)),
                                                 "from" = to_be_add$from,
                                                 "to" = to_be_add$to,
                                                 "nGRN" = ifelse(
                                                     coregnetwork@evidenceDescription[evname, "evidenceType"] == "coregulatory",
                                                     rep(Inf, nrow(to_be_add)), # Always show ppi data
                                                     rep(Inf, nrow(to_be_add)) # TODO : show TFBS conditionally to number of nGRN
                                                 ),
                                                 "title" =  paste(to_be_add$from, to_be_add$to),
                                                 "value" = rep(0.5, nrow(to_be_add)),
                                                 "arrows" = rep(
                                                     ifelse(
                                                         coregnetwork@evidenceDescription[evname, "evidenceType"] == "coregulatory",
                                                         "none",
                                                         "to"
                                                     ), nrow(to_be_add)),
                                                 "hidden" = rep(FALSE, nrow(to_be_add)),
                                                 "color" = rep(palette.evidence[evname], nrow(to_be_add)),
                                                 "physics" = rep(TRUE, nrow(to_be_add)),
                                                 stringsAsFactors = FALSE
                                             )
                                         }
                                  ))
    }
    # message("TODO : assert that added evidence correspond to predicted edges.")

    # edges data
    # ----------------
    regulatorsGRN <- sapply(names(GBMRegMap::refinedGRN@adjacencyList$bytf), function(x){
        tmp <- GBMRegMap::refinedGRN@adjacencyList$bytf[[x]]
        tmp_y <- length(unlist(tmp))
        return(tmp_y)
    }) %>% sort(decreasing = TRUE)


    to_save.edges = gtools::smartbind(coregs.with_sufficient_nGRN, additionalInts)

    # nodes data
    # ----------------
    nodes = regulatorsGRN[coregs.with_sufficient_nGRN.vector]
    # nodes name to show
    highlighted <- names(nodes)
    nameSize <- 50
    if (!is.null(keepNodes)) {
        idx <- setdiff(highlighted,keepNodes)
        highlighted[which(highlighted %in% idx)] <- ""
        nameSize <- 50
    }
    # node shape
    shapes <- rep("dot", length(nodes))
    if (!is.null(nmsSquare)) {
        idx <- which(names(nodes) %in% nmsSquare)
        shapes[idx] <- "square"
    }
    if (!is.null(nmsTriangle)) {
        idx <- which(names(nodes) %in% nmsTriangle)
        shapes[idx] <- "dot"
    }
    if (!is.null(nmsDmd)) {
        idx <- which(names(nodes) %in% nmsDmd)
        shapes[idx] <- "diamond"
    }

    # named vector (regulators) and how many gene they influence on
    to_save.nodes = data.frame(
        id = names(nodes),
        label = highlighted,
        #label = names(nodes), #when all are highlighted
        value = nodes / max(nodes), # how many target genes they regulate,
        # label inside version
        # shape = rep("circle", length(nodes)),
        # font = list(strokeWidth = 1, strokeColor = "white", color = "white"),
        # color = list(background = "purple", border = "purple"),
        # title = paste0("<p><b>", names(nodes),"</b><br>Node !</p>"), # tooltip (html or character), when the mouse is above
        # shadow = rep(FALSE, length(nodes)),
        # #scaling = list(min = 70, max = 100, label = list(min = 20, max = 50)) # node size
        # #scaling = list(min = 70, max = 150, label = list(min = 14, max = 30)) #if we do not put label, the node size depends on the name length
        # #scaling = list(min = 40, max = 70, label = list(min = 40, max = 70)) # long labels are given bigger edges than short ones...
        # scaling = list(min = 300, max = 1000, label = list(min = 40, max = 100)) # long labels are given bigger edges than short ones...
        # label outside version
        shape = shapes,
        font = list(strokeWidth = 4, strokeColor = "black", color = "black", size = nameSize), #80), when all are highlighted
        color = list(background = "lightblue", border = "lightblue"),
        title = paste0("<p><b>", names(nodes),"</b><br>Node !</p>"), # tooltip (html or character), when the mouse is above
        shadow = rep(FALSE, length(nodes)),
        scaling = list(min = 20, max = 80) #, label = list(min = 50, max = 40))
    )

    to_save.nodes = data.frame(to_save.nodes, TFA = sapply(to_save.nodes$id,
                                                           function(id) {
                                                               if (id %in% rownames(TFA)) {
                                                                   return(mean(TFA[as.character(id), ]))
                                                               } else {
                                                                   return(0)
                                                               }
                                                           }))

    minTFA <- min(to_save.nodes$TFA)
    maxTFA <- max(to_save.nodes$TFA)
    breaks <- seq(minTFA, maxTFA, by = (maxTFA - minTFA) / 75)
    pal <- leaflet::colorNumeric(colorRampPalette(c("blue","honeydew","red"))(75), domain = breaks, reverse = FALSE)

    to_save.nodes <- data.frame(to_save.nodes, color = pal(to_save.nodes$TFA))

    to_save.nodes$label <- to_save.nodes$id
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # RETURN VALUES
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    networkData <- list(to_save.nodes, to_save.edges)
    names(networkData) <- c("networkNodes","networkEdges")

    return(networkData)
}


