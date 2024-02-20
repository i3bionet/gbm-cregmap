# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Script to create a RegMap S4 object
# This object will store all necessary data for RegMap application :
#   - descriptiveData : used in front page of application
#   - influence : infleunce matrix
#   - UMAP2D : umap object storing umap::umap results for 2D UMAP
#   - seuratObject : Seurat object of considered 'cohort'
#   - classificationSota : named character vector
#   - classificationInfluence : named integer vector storing cluster::pam results or assignSVM results
#   - markersSota : differentially influence analysis on State of the Art (Sota)
#                   classification, performed by oneVsAllMarkers function (using Seurat::FindAllMarkers)
#   - markersPam : differentially influence analysis on pam clustering classification,
#                   performed by oneVsAllMarkers function (using Seurat::FindAllMarkers)
#                   or limmaOneVsAllDiffAnalysis function (using limma::lmFit and retreaving more
#                   differentially influent regulators)
#   - crisprCas9Data : CRISPR/Cas9 data from depmap project for CCLE data

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Set S4 class for RegMap object
setClass("RegMapObject",
         slots = c(
             descriptiveData = "ANY",
             influence = "matrix",
             expression = "ANY",
             UMAP2D = "ANY",
             seuratObject = "ANY",
             classificationSota = "ANY",
             classificationInfluence1 = "ANY",
             classificationInfluence1_table = "ANY",
             classificationInfluence2 = "ANY",
             markersSota = "data.frame",
             markersPam1 = "data.frame",
             markersPam2 = "data.frame",
             cnvData = "matrix",
             somaticMutations = "ANY",
             clinicalData = "data.frame",
             crisprCas9Data = "data.frame"
                   ),
         prototype = list(
             descriptiveData = NA,
             influence = matrix(),
             expression = NA,
             UMAP2D = NA,
             seuratObject = NA,
             classificationSota = NA,
             classificationInfluence1 = NA,
             classificationInfluence1_table = NA,
             classificationInfluence2 = NA,
             markersSota = data.frame(),
             markersPam1 = data.frame(),
             markersPam2 = data.frame(),
             cnvData = matrix(),
             somaticMutations = NA,
             clinicalData = data.frame(),
             crisprCas9Data = data.frame()
         ))

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Construct helper
RegMapObject <- function(descriptiveData = data.frame(),
                         influence,
                         expression = NA,
                         UMAP2D,
                         seuratObject = NA,
                         classificationSota = NA,
                         classificationInfluence1 = NA,
                         classificationInfluence1_table = NA,
                         classificationInfluence2 = NA,
                         markersSota = data.frame(),
                         markersPam1 = data.frame(),
                         markersPam2 = data.frame(),
                         cnvData = matrix(),
                         somaticMutations = NA,
                         clinicalData = data.frame(),
                         crisprCas9Data = data.frame()){
    new("RegMapObject",
        descriptiveData = descriptiveData,
        influence = influence,
        expression = expression,
        UMAP2D = UMAP2D,
        seuratObject = seuratObject,
        classificationSota = classificationSota,
        classificationInfluence1 = classificationInfluence1,
        classificationInfluence1_table = classificationInfluence1_table,
        classificationInfluence2 = classificationInfluence2,
        markersSota = markersSota,
        markersPam1 = markersPam1,
        markersPam2 = markersPam2,
        cnvData = cnvData,
        somaticMutations = somaticMutations,
        clinicalData = clinicalData,
        crisprCas9Data = crisprCas9Data
        )
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # Set validity controls
setValidity("RegMapObject",function(object){
    if (!all(names(object@classificationInfluence1) %in% colnames(object@influence))){
        "@classificationInfluence names must be in influence colnames"
    } else {
        TRUE
    }

    if (!all(names(object@classificationSota) %in% colnames(object@influence))){
        "@classificationSota names must be in influence colnames"
    } else {
        TRUE
    }
})









