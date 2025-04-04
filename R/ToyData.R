## ---------------------------
##
## Script name: HelperFunctions
##
## Purpose of script: General helper functions to check function input and save results
##
## Author: Christina Schmidt
##
## Date Created: 2023-06-14
##
## Copyright (c) Christina Schmidt
## Email:
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

################################################################################################
### ### ### Example Data ### ### ###
################################################################################################

#' Access built-in example data
#'
#' @param Dataset Character: name of a built-in dataset:
#'     \itemize{
#'         \item{\code{"IntraCells_Raw"}: }
#'         \item{\code{"IntraCells_DMA"}: }
#'         \item{\code{"CultureMedia_Raw"}: }
#'         \item{\code{"Cells_MetaData"}: }
#'         \item{\code{"Tissue_Norm"}: }
#'         \item{\code{"Tissue_MetaData"}: }
#'         \item{\code{"Tissue_DMA"}: }
#'         \item{\code{"Tissue_DMA_Old"}: }
#'         \item{\code{"Tissue_DMA_Young"}: }
#'         \item{\code{"Tissue_TvN_Proteomics"}: }
#'         \item{\code{"Tissue_TvN_RNAseq"}: }
#'         \item{\code{"EquivalentFeatures"}: }
#'         \item{\code{"BiocratesFeatureTable"}: }
#'     }
#'
#' @return A data frame containing the toy data.
#'
#' @description Import and process .csv file to create toy data DF.
#'
#' @examples
#' Intra <- MetaProViz::ToyData("IntraCells_Raw")
#'
#' @importFrom readr read_csv cols
#' @importFrom magrittr %>% extract2
#' @importFrom tibble column_to_rownames
#' @importFrom logger log_trace
#'
#' @export
#'
ToyData <- function(Dataset) {
  ## ------------ Create log file ----------- ##
  MetaProViz_Init()

  #Available Datasets:
  datasets <- list(
    IntraCells_Raw = "MS55_RawPeakData.csv.gz",
    IntraCells_DMA = "MS55_DMA_786M1A_vs_HK2.csv.gz",
    CultureMedia_Raw = "MS51_RawPeakData.csv.gz",
    Cells_MetaData = "MappingTable_SelectPathways.csv.gz",
    Tissue_Norm = "Hakimi_ccRCC-Tissue_Data.csv.gz",
    Tissue_MetaData = "Hakimi_ccRCC-Tissue_FeatureMetaData.csv.gz",
    Tissue_DMA = "Hakimi_ccRCC-Tissue_DMA_TvsN.csv.gz",
    Tissue_DMA_Old ="Hakimi_ccRCC-Tissue_DMA_TvsN-Old.csv.gz",
    Tissue_DMA_Young ="Hakimi_ccRCC-Tissue_DMA_TvsN-Young.csv.gz",
    Tissue_TvN_Proteomics ="ccRCC-Tissue_TvN_Proteomics.csv.gz",
    Tissue_TvN_RNAseq = "ccRCC-Tissue_TvN_RNAseq.csv.gz",
    AlaninePathways = "AlaninePathways.csv.gz",
    EquivalentFeatures = "EquivalentFeatureTable.csv.gz",
    BiocratesFeatureTable = "BiocratesFeatureTable.csv.gz"
  )

  rncols <- c("Code", "Metabolite")

  #Load dataset:
  if (!Dataset %in% names(datasets)) {
    message <- sprintf("No such dataset: `%s`. Available datasets: %s", Dataset, paste(names(datasets), collapse = ", "))
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }

  datasets %>%
  magrittr::extract2(Dataset) %>%
    system.file("extdata", ., package = "MetaProViz") %>%
    readr::read_csv(col_types = readr::cols()) %>%
    {`if`(
    (rncol <- names(.) %>% intersect(rncols)) %>% length,
    tibble::column_to_rownames(., rncol),
    .
    )}
}

