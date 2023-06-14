## ---------------------------
##
## Script name: HelperFunctions
##
## Purpose of script: Metabolomics (raw ion counts) pre-processing, normalisation and outlier detection
##
## Author: Dimitrios Prymidis and Christina Schmidt
##
## Date Created: 2023-06-14
##
## Copyright (c) Dimitrios Prymidis and Christina Schmidt
## Email:
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------


#' Imports toy data into environment
#'
#' @title Toy Data Import
#' @description Import and process .csv file to create toy data.
#' @importFrom utils read.csv
#' @return A data frame containing the toy data.
#' @export
toy_data <- function() {
  # Read the .csv files
  Intra <- system.file("data", "MS55_RawPeakData.csv", package = "MetaProViz")
  Intra<- read.csv(Intra, check.names=FALSE)
  
  Media <- system.file("data", "MS51_RawPeakData.csv", package = "MetaProViz")
  Media<- read.csv(Media, check.names=FALSE)
  
  Pathways <-system.file("data", "MappingTable_SelectPathways.csv", package = "MetaProViz")
  Pathways<- read.csv(Pathways, check.names=FALSE)
  
  # Return the toy data into environment
  assign("Intra", Intra, envir=.GlobalEnv)
  assign("Media", Intra, envir=.GlobalEnv)
  assign("Pathways", Intra, envir=.GlobalEnv)
}



