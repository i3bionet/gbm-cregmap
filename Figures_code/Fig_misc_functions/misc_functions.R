library(dplyr)
library(ggplot2)


recodeTCGA <- function(names){
#' recodeTCGA
#'
#' @param names vector ....
#'
#' @return vector
#' @export
#'
#' @examples
#'
    sapply(names, function(x){
        if (grepl("TCGA", x)){
            y <- gsub("\\.", "-", x = x)
            y <- substring(y,1,12)
            return(y)
        } else {
            return(x)
        }
    })
}

#' geneConversionHGNC
#'
#' Function to proceed gene conversion from gene conversion list
#'
#' @param data dataframe to convert with gene as rownames and samples as colnames
#' @param conversionList list of conversion/correspondance table
#'
#' @return
#' @export
#'
#' @examples
#'
geneConversionHGNC <- function(data, conversionList){
    ## Remove extra blanck space in rownames (that was the cas for GSE214814)
    tmpRownames <- gsub(pattern = " ", replacement = "", rownames(data))
    rownames(data) <- tmpRownames
    # Parse conversion lists and return number of matched genes
    nMatch <- sapply(names(conversionList), function(x){
        convTable <- conversionList[[x]]
        numberMatches <- sum(rownames(data) %in% convTable[,2])
    })
    maxMatch <- names(nMatch[nMatch == max(nMatch)])
    ## Check if multiple conversion list have same matched genes
    if (length(maxMatch) > 1){
        maxMatch <- maxMatch[[1]]
    }
    if (max(nMatch) < 8000) {
        stop(paste("Not enought gene can be converted :", max(nMatch)))
    }
    print(paste("Maximum of matches for :",maxMatch, "with ",nMatch[maxMatch], " genes converted."))
    ## Select table for wich there is a maximum of gene matched
    tmpConvTable <- conversionList[[maxMatch]]
    ## Convert probes IDs into characters
    tmpConvTable[,2] <- as.character(tmpConvTable[,2])
    ## Remove NAs
    tmpConvTable <- tmpConvTable[!is.na(tmpConvTable[,2]),]
    ## Keep gene that are inside conversion table
    tmpConvTable <- tmpConvTable[tmpConvTable[,2] %in% rownames(data),]
    ## If there is duplicated gene HGNC IDs keep first
    tmpConvTable <- tmpConvTable[!duplicated(tmpConvTable[,1]),]



    # ## Aligne genes
    # tmpConvTable <- tmpConvTable[na.omit(match(rownames(data), tmpConvTable[,2])),]
    commonGenes <- intersect(tmpConvTable[,2], rownames(data))
    tmpConvTable <- tmpConvTable[match(commonGenes, tmpConvTable[,2]),]
    ## reorder data
    data <- data[tmpConvTable[,2],]
    ## Reorder tmpConvTable
    # tmpConvTable <- tmpConvTable[na.omit(match(rownames(data),intersect(rownames(data), tmpConvTable[,2]))),]
    ## Change rownames
    rownames(data) <- tmpConvTable[,1]

    return(data)
}


#' oneVsAllMarkers
#' Function to have gene marker for cluster
#'
#' @param data influence or normalized expression matrix
#' @param classification
#' @param only.pos
#' @param logfc.threshold
#'
#' @return
#' @export
#'
#' @examples
oneVsAllMarkers <- function(data, classification, only.pos = TRUE, logfc.threshold=0.25){
    #Append cell type to infleunce colnames to specify groups in seurat object
    classification <- classification[colnames(data)]
    classificationSeurat <- classification
    # Substitute "_" character in data colnames (will cause problem in differential analysis)
    colnames(data) <- gsub(pattern = "_",replacement = ".",x = colnames(data))
    # Change previous "_" in classification with another character
    tmpNames <- names(classificationSeurat)
    classificationSeurat <- gsub(classificationSeurat, pattern = "_", replacement = ".")
    names(classificationSeurat) <- tmpNames
    dataSeurat <- data
    colnames(dataSeurat) <- paste(colnames(dataSeurat), classificationSeurat, sep = "_")

    dataSeuratObject <- Seurat::CreateSeuratObject(dataSeurat, names.field = 2)
    # Identify gene markers
    all_markers <- Seurat::FindAllMarkers(dataSeuratObject,
                                          only.pos = only.pos,
                                          min.pct = 0,
                                          min.cells.group = 1,
                                          verbose = FALSE,
                                          logfc.threshold = logfc.threshold)
    # ## Change p_val_adj values by adjusting p_val with Benjamini and Hochberg method
    # ## This method is too conservative for our analysis
    # all_markers$p_val_adj <- p.adjust(all_markers$p_val, method = "fdr")
    # Condition if markers dataframe is not empty
    if (!(nrow(all_markers) == 0)) {
        # Create a new variable to specifiy if regulator are UP, DOWN influent or
        # not influent
        all_markers$diffInfluent <- apply(all_markers, 1, function(x){
            tmpLog2 <- as.numeric(x["avg_log2FC"])
            tmpPVal <- as.numeric(x["p_val_adj"])
            if ((tmpLog2 < 0) & (tmpPVal< 0.05)){
                return("DOWN")
            } else if ((tmpLog2 > 0) & (tmpPVal < 0.05)){
                return("UP")
            } else {
                return("NO")
            }
        })
        # Create variable to specifiy labelisation on volcanoplots
        # If not influent -> no label
        all_markers$label <- apply(all_markers,1, function(x){
            label = x[["gene"]]
            if (x["diffInfluent"] == "NO"){
                return(NULL)
            } else {
                return(label)
            }
        })

        # Create new variable 'adj_pVal_volcano' where there is no p_value equal to
        # zero (-log10(0) = Inf -> ylim = Inf)
        ## Get the lower pValue after 0
        tmpMinPval <- min(all_markers$p_val_adj[!all_markers$p_val_adj == 0])
        ## Set a new lowest pValue equal
        newLowestPval <- tmpMinPval/10
        all_markers$adj_pVal_volcano <- ifelse(test = all_markers$p_val_adj == 0,
                                               yes = newLowestPval,
                                               no = all_markers$p_val_adj)
        # Add influence info



        # reformat names(classification)
        names(classification) <- gsub(pattern = "_", replacement = ".", x = names(classification))



        tmpInfColumns <- apply(all_markers, 1, function(x){
            gene <- x[["gene"]]
            cluster <- x[["cluster"]]
            clusterInfluence <- mean(data[gene, names(classification)[classification == cluster]])
            restInfluence <- mean(data[gene, names(classification)[classification != cluster]])
            return(c(clusterInfluence, restInfluence))
        }) %>% t() %>% data.frame()
        colnames(tmpInfColumns) <- c("clusterMeanInfluence", "restMeanInfluence")

        all_markers <- cbind(all_markers, tmpInfColumns)

    }
    # ---------------------------------------------------------------------------- %
    # add a "name" column with gene's name without supplementary '.X' from rownames
    # --------------------------- #
    ## Create "name" column
    all_markers$name <- rownames(all_markers)
    # --------------------------- #
    # remove supplementary ".X"
    all_markers$name <- sapply(all_markers$name, function(x){
        strsplit(x, split = "\\.")[[1]][1]
    }) %>% unname()


    return(all_markers)
}

#' quantile_df
#'
#' @param x tibble object
#' @param probs vector of probabilities
#'
#' @return
#' @export
#'
#' @examples
quantile_df <- function(x, probs = c(0, 0.25, 0.5, 0.75, 1)) {
    tibble(quantile = probs, value = quantile(x, probs))
}

# load_rdata------------------------------------------------
#' load_rdata
#'
#' @param file_path
#'
#' @return
#' @export
#'
#' @examples
load_rdata <- function(file_path) {
    res <- local({
        load(file_path)
        return(get(ls()))
    })
    return(res)
}

# object to format classification levels
re_format <- c(
    "PN"     = "C3",
    "PN-L"   = "C7",
    "NL"     = "C5",
    "CL-A"   = "C2",
    "CL-B"   = "C6",
    "CL-C"   = "C4",
    "MES"    = "C1"
)

# object to Declare genes associated to BP for heatmaps
Mesenchymal <- c("DKK1", "BCL3",
                 "TGFB1","TGFBI","SERPINE1","TIMP1","RELB",
                 "CHI3L1","TRADD","MET","NF1", "SNAI2")

Classical <- c("MEOX2","CDH4","NES","PDGFA","KCNF1","NOTCH3",
               "JAG1","SMO","GAS1","GLI2")

Proneural <- c("DLL3","NKX2-2","OLIG2","PDGFRA")
Proliferation <- c("MKI67","DLGAP5")

LateCellCycle <- c("CCNE2","CDC25A","BUB1","CDC20",
                   "CCNA2","CCNB1")
Stemness <- c("KLF4","SALL4","POU5F1","NANOG", "SERPINB2","SERPINB3")

extraCellMatrixRemod <- c("COL1A1","COL5A1","COL6A2",
                          "LOX","LOXL1","TAGLN","CD44")
egfrMapkPathways <- c("GAB1", "EGFR")
Neural <- c("NEFL","GABRA1", "GABRA2","GABRA4","GABRA5", "KCNJ3","SYT1")
coagulation <- c("C8A","FGA","FGB")
gliogenesis <- c("ZNF488","NOG","DLL1","BMP2", "MAPT","CNTN1")
cilium <- c("EFCAB6", "DNAH7",  "DRC7", "CCDC65", "CFAP206","SPEF1","RSPH4A" )
FAM <- c("ELOVL2", "ACSL3", "PLA2G5")
kerat <- c("KRT6A", "KDF1", "SPRR1B", "KRT16", "KRT6C")
stem_pr <- c("CDX2", "PON1", "ADH4", "RBP2", "LCAT", "GAL", "LPA", "CYP19A1",
             "GAST", "SST", "UGT1A8")
# Define markers
markers <- list(
    Proliferation,       #"PN"
    LateCellCycle,       #"PN"
    Proneural,           #"PN"
    gliogenesis,
    Neural,              #"NL"

    Classical,           #"CL-A"
    egfrMapkPathways,    #"CL-A"
    FAM,

    Stemness,            #"CL-B"
    coagulation,         #"CL-B"
    kerat,               #"CL-B"

    cilium,              #"CL-C"

    Mesenchymal,         #"MES"
    extraCellMatrixRemod #"MES"
)

names(markers) <- c(
    "Proliferation",
    "Late cell cycle",
    "Proneural",
    "Gliogenesis",
    "Neural",

    "Classical",
    "MAPK pathway",
    "Fatty Acid Synthesis",

    "Stemness",
    "Coagulation",
    "Keratinization",
    "Cillia",

    "Mesenchymal",
    "Extracellular Matrix\nRemodeling"
)

markers <- markers %>%
    unlist() %>%
    tibble::enframe(name = "Signature", value = "gene") %>%
    dplyr::mutate(Signature = stringr::str_remove(Signature, "[0-9]{1,2}$") %>%
                      forcats::as_factor())

classes_to_colors_Verhak <- c(
    "Classical" = "#A3A500",
    "Mesenchymal" = "#F8766D",
    "Neural"    = "#00B6EB",
    "Proneural" = "#53B400"
)

levels_gbm_Verhaak <-c("Proneural" = "G-CIMP")
size_p <- 24
pub_theme_purre <- list(
    theme(
        title        = element_text(size = size_p,
                                    colour = "black", face = "bold"),
        legend.title = element_text(size = size_p,
                                    colour = "black", face = "bold"),
        legend.text  = element_text(size = size_p,
                                    colour = "black", face = "bold"),
        axis.text.x  = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y  = element_text(size = size_p,
                                    colour = "black", face = "bold"),
        axis.title.y = element_text(size = size_p,
                                    colour = "black", face = "bold"),
        strip.text   = element_text(size = size_p, colour = "black",
                                    face = "bold"),
        axis.line    = element_line(color='black', linewidth = 1.3),

        axis.ticks.x     = element_blank(),
        #plot.background  = element_blank(),
        legend.position="bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border     = element_rect(color = "black", linewidth = 1.3)
    )
)

