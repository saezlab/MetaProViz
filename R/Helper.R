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
#' @param data Either "Standard", "CoRe" or "Pathways" depending which data you would like to load
#' @title Toy Data Import
#' @description Import and process .csv file to create toy data.
#' @importFrom utils read.csv
#' @return A data frame containing the toy data.
#' @export
#'
toy_data <- function(data) {
  # Read the .csv files
  Intra <- system.file("data", "MS55_RawPeakData.csv", package = "MetaProViz")
  Intra<- read.csv(Intra, check.names=FALSE)

  Media <- system.file("data", "MS51_RawPeakData.csv", package = "MetaProViz")
  Media<- read.csv(Media, check.names=FALSE)

  Pathways <-system.file("data", "MappingTable_SelectPathways.csv", package = "MetaProViz")
  Pathways<- read.csv(Pathways, check.names=FALSE)

  # Return the toy data into environment
  if(data=="Standard"){
    assign("Intra", Intra, envir=.GlobalEnv)
  } else if(data=="CoRe"){
    assign("Media", Media, envir=.GlobalEnv)
  } else if(data=="Pathways"){
    assign("Pathways", Pathways, envir=.GlobalEnv)
  } else{
    warning("Please choose a toy dataset you would like to use: Standard, CoRe, Pathways")
  }
}


#' Imports MCA regulatory rules into environment
#'
#' @param data Either "2Cond" or "CoRe" depending which regulatory rules you would like to load
#' @title MCA regulatory rules Import
#' @description Import and process .csv file to create toy data.
#' @importFrom utils read.csv
#' @return A data frame containing the toy data.
#' @export
#'
MCA_rules<- function(data) {
  # Read the .csv files
  Cond <- system.file("data", "MCA_2Cond.csv", package = "MetaProViz")
  Cond<- read.csv( Cond, check.names=FALSE)

  CoRe <- system.file("data", "MCA_CoRe.csv", package = "MetaProViz")
  CoRe<- read.csv(CoRe, check.names=FALSE)

  # Return the toy data into environment
  if(data=="2Cond"){
    assign("MCA_2Cond", Cond, envir=.GlobalEnv)
  } else if(data=="CoRe"){
    assign("MCA_CoRe", CoRe, envir=.GlobalEnv)
  } else{
    warning("Please choose the MCA regulatory rules you would like to load: 2Cond, CoRe")
  }
}

