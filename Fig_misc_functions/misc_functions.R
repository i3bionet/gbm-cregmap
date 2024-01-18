recodeTCGA <- function(names){
#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
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
#' @param conversionList a list of conversion/correspondance table
#'
#' @return
#' @export
#'
#' @examples
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

# Declare genes associated to BP for heatmaps
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
    enframe(name = "Signature", value = "gene") %>%
    mutate(Signature = str_remove(Signature, "[0-9]{1,2}$") %>% as_factor())

classes_to_colors_Verhak <- c(
    "Classical" = "#A3A500",
    "Mesenchymal" = "#F8766D",
    "Neural"    = "#00B6EB",
    "Proneural" = "#53B400"
)
