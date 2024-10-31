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
#'     }
#'
#' @return A data frame containing the toy data.
#'
#' @description Import and process .csv file to create toy data.
#'
#' @examples
#' Intra <- ToyData("IntraCells_Raw")
#'
#' @importFrom readr read_csv cols
#' @importFrom magrittr %>% extract2
#' @importFrom tibble column_to_rownames
#'
#' @export
#'
ToyData <- function(Dataset) {
  ## ------------ Create log file ----------- ##
  MetaProViz_Init()

  #Available Datasets:
  datasets <- list(
    IntraCells_Raw = "MS55_RawPeakData.csv",
    IntraCells_DMA = "MS55_DMA_786M1A_vs_HK2.csv",
    CultureMedia_Raw = "MS51_RawPeakData.csv",
    Cells_MetaData = "MappingTable_SelectPathways.csv",
    Tissue_Norm = "Hakimi_ccRCC-Tissue_Data.csv",
    Tissue_MetaData = "Hakimi_ccRCC-Tissue_FeatureMetaData.csv",
    Tissue_DMA = "Hakimi_ccRCC-Tissue_DMA_TvsN.csv",
    Tissue_DMA_Old ="Hakimi_ccRCC-Tissue_DMA_TvsN-Old.csv",
    Tissue_DMA_Young ="Hakimi_ccRCC-Tissue_DMA_TvsN-Young.csv"
  )

  rncols <- c("Code", "Metabolite")

  #Load dataset:
  if (!Dataset %in% names(datasets)) {
    stop(sprintf(
      "No such dataset: `%s`. Available datasets: %s",
      Dataset,
      paste(names(datasets), collapse = ", ")
    ))
  }

  datasets %>%
  magrittr::extract2(Dataset) %>%
    system.file("data", ., package = "MetaProViz") %>%
    readr::read_csv(col_types = cols()) %>%
    {`if`(
    (rncol <- names(.) %>% intersect(rncols)) %>% length,
    tibble::column_to_rownames(., rncol),
    .
    )}
}

