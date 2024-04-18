## ---------------------------
##
## Script name: Make_GeneMetabSet
##
## Purpose of script: Create gene-metabolite sets for pathway enrichment analysis.
##
## Author: Christina Schmidt
##
## Date Created: 2024-01-21
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

##########################################################################################
### ### ### Get KEGG prior knowledge ### ### ###
##########################################################################################

#' Imports KEGG pathways into the environment
#'
#' @title KEGG
#' @description Import and process KEGG.
#' @importFrom utils read.csv
#' @return A data frame containing the KEGG pathways for ORA.
#' @export
#'
#'
LoadKEGG <- function(){
  #Get the package:
  RequiredPackages <- c("rappdirs")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)

  #------------------------------------------------------------------
  #Get the directory and filepath of cache results of R
  directory <- rappdirs::user_cache_dir()#get chache directory
  File_path <-paste(directory, "/KEGG_Metabolite.rds", sep="")

  if(file.exists(File_path)==TRUE){# First we will check the users chache directory and weather there are rds files with KEGG_pathways already:
    KEGG_Metabolite <- readRDS(File_path)
    message("Cached file loaded from: ", File_path)
  }else{# load from KEGG
    RequiredPackages <- c("KEGGREST", "tidyverse")
    new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
    if(length(new.packages)) install.packages(new.packages)

    suppressMessages(library(KEGGREST))
    suppressMessages(library(tidyverse))

    #--------------------------------------------------------------------------------------------
    # 1. Make a list of all available human pathways in KEGG
    Pathways_H <- as.data.frame(keggList("pathway", "hsa"))  # hsa = human

    # 2. Initialize the result data frame
    KEGG_H <- data.frame(KEGGPathway = character(nrow(Pathways_H)),
                         #PathID = 1:nrow(Pathways_H),
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
        #KEGG_H[k, "PathID"] <- path
        KEGG_H[k, "Compound"] <- paste(query[[1]]$COMPOUND, collapse = ";")
        KEGG_H[k, "KEGG_CompoundID"] <- paste(names(query[[1]]$COMPOUND), collapse = ";")
      }, error = function(e) {
        # If an error occurs, store "error" in the corresponding row and continue to the next query
        KEGG_H[k, "KEGGPathway"] <- "error"
        message(paste("`Error in .getUrl(url, .flatFileParser) : Not Found (HTTP 404).` for pathway", path, "- Skipping and continuing to the next query."))
      })
    }

    # 3. Remove the pathways that do not have any metabolic compounds associated to them
    KEGG_H_Select <-KEGG_H%>%
      subset(!KEGG_CompoundID=="")%>%
      subset(!KEGGPathway=="")

    # 4. Make the Metabolite DF
    KEGG_Metabolite <- separate_longer_delim(KEGG_H_Select[,-5], c(Compound, KEGG_CompoundID), delim = ";")

    # 5. Remove Metabolites
    ### 5.1. Ions should be removed
    Remove_Ions <- c("Calcium cation","Potassium cation","Sodium cation","H+","Cl-", "Fatty acid", "Superoxide","H2O", "CO2", "Copper", "Fe2+", "Magnesium cation", "Fe3+",  "Zinc cation", "Nickel", "NH4+")
    ### 5.2. Unspecific small molecules
    Remove_Small <- c("Nitric oxide","Hydrogen peroxide", "Superoxide","H2O", "CO2", "Hydroxyl radical", "Ammonia", "HCO3-",  "Oxygen", "Diphosphate", "Reactive oxygen species", "Nitrite", "Nitrate", "Hydrogen", "RX", "Hg")

    KEGG_Metabolite <- KEGG_Metabolite[!(KEGG_Metabolite$Compound %in% c(Remove_Ions, Remove_Small)), ]

    #Change syntax as required for ORA
    KEGG_Metabolite <- KEGG_Metabolite%>%
      dplyr::rename("term"=1,
                    "Metabolite"=2,
                    "MetaboliteID"=3)
    KEGG_Metabolite$Description <- KEGG_Metabolite$term

    #Save the results as an RDS file in the Cache directory of R
    if(!dir.exists(directory)) {dir.create(directory)}
    saveRDS(KEGG_Metabolite, file = paste(directory, "/KEGG_Metabolite.rds", sep=""))
  }

  #Return into environment
  assign("KEGG_Pathways", KEGG_Metabolite, envir=.GlobalEnv)
}



##########################################################################################
### ### ### Get Metabolite Pathways using Cosmos prior knowledge ### ### ###
##########################################################################################
#'
#' Function to add metabolite HMDB IDs to existing genesets based on cosmosR prior knowledge
#'
#' @param Input_GeneSet Dataframe with two columns for source (=term) and target (=gene), e.g. Hallmarks.
#' @param target \emph{Optional: }  Column name of target in Input_GeneSet. \strong{Default = "gene}
#' @param OutputName {Optional:} String which is added to the output files name.\strong(Default = NULL)
#' @param Folder_Name {Optional:} String which is added to the resulting folder name \strong(Default = NULL)
#'
#' @export

Make_GeneMetabSet <- function(Input_GeneSet,
                              target=c("gene"),
                              OutputName=NULL,
                              Folder_Name = NULL){

  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse","cosmosR")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  suppressMessages(library(tidyverse))

  ## ------------ Check Input files ----------- ##
  # 1. The input data:
  if(is.data.frame(Input_GeneSet)==FALSE){
    stop("`Input_GeneSet` must be of class data.frame with columns for source (=term) and target (=gene). Please check your input")
  }
  # 2. target:

  if(is.character(target)==FALSE | (length(target) ==1)==FALSE){
    stop("`target` must be of class character withonly one entry. Please check your input")
  }

  ## ------------ Folder ----------- ##
  if(is.null(Folder_Name)){
    name <- paste("MetaProViz_Results",Sys.Date(),sep = "_" )
  }else{
    if(grepl('[^[:alnum:]]', Folder_Name)){
      stop("The 'Folder_Name' must not contain any special character.")
    }else{
      name <- paste("MetaProViz_Results",Sys.Date(),Folder_Name,sep = "_" )
    }
  }
  WorkD <- getwd()
  Results_folder <- file.path(WorkD, name)
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)} # Make Results folder
  Results_folder_GeneMetabSet_folder = file.path(Results_folder, "GeneMetabSets")  # This searches for a folder called "GeneMetabSet" within the "Results" folder in the current working directory and if its not found it creates one
  if (!dir.exists(Results_folder_GeneMetabSet_folder )) {dir.create(Results_folder_GeneMetabSet_folder)}  # check and create folder

  ######################################################
  ##-------------- Cosmos PKN
  #load the network from cosmos
  data("meta_network", package = "cosmosR")
  meta_network <- meta_network[which(meta_network$source != meta_network$target),]

  #adapt to our needs extracting the metabolites:
  meta_network_metabs <- meta_network[grepl("Metab__HMDB", meta_network$source) | grepl("Metab__HMDB", meta_network$target),-2]#extract entries with metabolites in source or target
  meta_network_metabs <- meta_network_metabs[grepl("Gene", meta_network_metabs$source) | grepl("Gene", meta_network_metabs$target),]#extract entries with genes in source or target

  #Get reactant and product
  meta_network_metabs_reactant <-  meta_network_metabs[grepl("Metab__HMDB", meta_network_metabs$source),]%>% dplyr::rename("metab"=1, "gene"=2)
  meta_network_metabs_products <-  meta_network_metabs[grepl("Metab__HMDB", meta_network_metabs$target),]%>% dplyr::rename("gene"=1, "metab"=2)

  meta_network_metabs <- as.data.frame(rbind(meta_network_metabs_reactant, meta_network_metabs_products))
  meta_network_metabs$gene <- gsub("Gene.*__","",meta_network_metabs$gene)
  meta_network_metabs$metab <- gsub("_[a-z]$","",meta_network_metabs$metab)
  meta_network_metabs$metab <- gsub("Metab__","",meta_network_metabs$metab)

  ##--------------metalinks transporters
  #Add metalinks transporters to Cosmos PKN

  ##-------------- Combine with Input_GeneSet
  #add pathway names --> File that can be used for metabolite pathway analysis
  MetabSet <- merge(meta_network_metabs,Input_GeneSet, by.x="gene", by.y=target)

  #combine with pathways --> File that can be used for combined pathway analysis (metabolites and gene t.vals)
  GeneMetabSet <- unique(as.data.frame(rbind(Input_GeneSet%>%dplyr::rename("feature"=target), MetabSet[,-1]%>%dplyr::rename("feature"=1))))

  ##-------------- Save and return
  if(is.null(OutputName)){
    write.csv(GeneMetabSet, paste(Results_folder_GeneMetabSet_folder, "/GeneMetabSet.csv", sep = ""), row.names = FALSE)
  }else{
    write.csv(GeneMetabSet, paste(Results_folder_GeneMetabSet_folder, "/GeneMetabSet_",   OutputName ,".csv", sep = ""), row.names = FALSE)
  }

  return(invisible(GeneMetabSet))
}



##########################################################################################
### ### ### Load MetaLinksDB prior knowledge ### ### ###
##########################################################################################
#'
#' Function to
#'
#' @param x test
#'
#' @export

LoadMetalinks <- function(){
  metalinks_db_url <- "https://figshare.com/ndownloader/files/44624707?private_link=4744950f8768d5c8f68c"

  # Download the Metalinks database file
  download.file(metalinks_db_url, destfile = "metalinks.db")

  # Connect to the SQLite database
  con <- DBI::dbConnect(RSQLite::SQLite(), "metalinks.db")
  tables <- DBI::dbListTables(con)

  query <- "SELECT * FROM metalinks"
  metalinks_data <- DBI::dbGetQuery(con, query)


  columns_query <- "PRAGMA table_info(metalinks)"
  columns_info <-  DBI::dbGetQuery(con, columns_query)


}
