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


#' Imports KEGG pathways into the environment
#'
#' @title MCA regulatory rules Import
#' @description Import and process .csv file to create toy data.
#' @importFrom utils read.csv
#' @return A data frame containing the KEGG pathways for ORA.
#' @export
#'
#'
Load_KEGG<- function(){
  #Get the package:
  RequiredPackages <- c("KEGGREST", "tidyverse")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)

  suppressMessages(library(KEGGREST))
  suppressMessages(library(tidyverse))

  # 1. Make a list of all available human pathways in KEGG
  Pathways_H <- as.data.frame(keggList("pathway", "hsa"))  # hsa = human

  # 2. Initialize the result data frame
  KEGG_H <- data.frame(KEGGPathway = character(nrow(Pathways_H)),
                       PathID = 1:nrow(Pathways_H),
                       Compound = 1:nrow(Pathways_H),
                       KEGG_CompoundID = 1:nrow(Pathways_H),
                       stringsAsFactors = FALSE)

  # 3. Iterate over each pathway and extract the needed information
  for (k in 1:nrow(Pathways_H)) {
    path <- rownames(Pathways_H)[k]

    tryCatch({#try-catch block is used to catch any errors that occur during the query process
      # Query the pathway information
      query <- keggGet(path)

      # Extract the necessary information and store it in the result data frame
      KEGG_H[k, "KEGGPathway"] <- Pathways_H[k,]
      KEGG_H[k, "PathID"] <- path
      KEGG_H[k, "Compound"] <- paste(query[[1]]$COMPOUND, collapse = ";")
      KEGG_H[k, "KEGG_CompoundID"] <- paste(names(query[[1]]$COMPOUND), collapse = ";")
    }, error = function(e) {
      # If an error occurs, store "error" in the corresponding row and continue to the next query
      KEGG_H[k, "KEGGPathway"] <- "error"
      message(paste("`Error in .getUrl(url, .flatFileParser) : Not Found (HTTP 404).` for pathway", path, "- Skipping and continuing to the next query."))
    })
  }

  # 3. Remove the pathways that do not have any metabolic compounds associated to them
  KEGG_H_Select <- subset(KEGG_H, !KEGG_H$KEGG_CompoundID=="")

  # 4. Make the Metabolite DF
  KEGG_Metabolite <- separate_longer_delim(KEGG_H_Select[,-5], c(Compound, KEGG_CompoundID), delim = ";")

  # 5. Remove Metabolites
  ### 5.1. Ions should be removed
  Remove_Ions <- c("Calcium cation","Potassium cation","Sodium cation","H+","Cl-", "Fatty acid", "Superoxide","H2O", "CO2", "Copper", "Fe2+", "Magnesium cation", "Fe3+",  "Zinc cation", "Nickel", "NH4+")
  ### 5.2. Unspecific small molecules
  Remove_Small <- c("Nitric oxide","Hydrogen peroxide", "Superoxide","H2O", "CO2", "Hydroxyl radical", "Ammonia", "HCO3-",  "Oxygen", "Diphosphate", "Reactive oxygen species", "Nitrite", "Nitrate", "Hydrogen", "RX", "Hg")

  KEGG_Metabolite <- KEGG_Metabolite[!(KEGG_Metabolite$Compound %in% c(Remove_Ions, Remove_Small)), ]

  #Return into environment
  assign("KEGG_Pathways", KEGG_Metabolite, envir=.GlobalEnv)

}
