## ---------------------------
##
## Script name: GetPriorKnowledge
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
  ## ------------ Create log file ----------- ##
  MetaProViz_Init()

  logger::log_info("Load KEGG.")

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
### ### ### Load Hallmark prior knowledge ### ### ###
##########################################################################################
#'
#' @title Toy Data Import
#' @description Import and process .csv file to create toy data.
#' @importFrom utils read.csv
#' @return A data frame containing the toy data.
#' @export
#'
LoadHallmarks <- function() {
  ## ------------ Create log file ----------- ##
  MetaProViz_Init()

  # Read the .csv files
  Hallmark <- system.file("data", "Hallmarks.csv", package = "MetaProViz")
  Hallmark <- read.csv(Hallmark, check.names=FALSE)

  # Return into environment
  assign("Hallmark_Pathways", Hallmark, envir=.GlobalEnv)
}



##########################################################################################
### ### ### Get Metabolite Pathways using Cosmos prior knowledge ### ### ###
##########################################################################################
#'
#' Function to add metabolite HMDB IDs to existing genesets based on cosmosR prior knowledge
#'
#' @param Input_GeneSet Dataframe with two columns for source (=term) and Target (=gene), e.g. Hallmarks.
#' @param SettingsInfo \emph{Optional: }  Column name of Target in Input_GeneSet. \strong{Default = c(Target="gene")}
#' @param PKName \emph{Optional: } Name of the prior knowledge resource. \strong{default: NULL}
#' @param SaveAs_Table \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt". \strong{Default = "csv"}
#' @param FolderPath {Optional:} String which is added to the resulting folder name \strong{default: NULL}
#'
#' @export

Make_GeneMetabSet <- function(Input_GeneSet,
                              SettingsInfo=c(Target="gene"),
                              PKName= NULL,
                              SaveAs_Table = "csv",
                              FolderPath = NULL){

  ## ------------ Create log file ----------- ##
  MetaProViz_Init()

  logger::log_info("Make_GeneMetabSet.")

  ## ------------ Check Input files ----------- ##
  # 1. The input data:
  if(is.data.frame(Input_GeneSet)==FALSE){
    stop("`Input_GeneSet` must be of class data.frame with columns for source (=term) and Target (=gene). Please check your input")
  }
  # 2. Target:
  if("Target" %in% names(SettingsInfo)){
    if(SettingsInfo[["Target"]] %in% colnames(Input_GeneSet)== FALSE){
      stop("The ", SettingsInfo[["Target"]], " column selected as Conditions in SettingsInfo was not found in Input_GeneSet. Please check your input.")
    }
  }else{
    stop("Please provide a column name for the Target in SettingsInfo.")
  }

  if(is.null(PKName)){
    PKName <- "GeneMetabSet"
  }

  ## ------------ Folder ----------- ##
  if(is.null(SaveAs_Table)==FALSE){
    Folder <- SavePath(FolderName= "PriorKnowledge",
                                    FolderPath=FolderPath)

    SubFolder <- file.path(Folder, "MetabSet")
    if (!dir.exists(SubFolder)) {dir.create(SubFolder)}
  }


  ######################################################
  ##-------------- Cosmos PKN
  #load the network from cosmos
  data("meta_network", package = "cosmosR")
  meta_network <- meta_network[which(meta_network$source != meta_network$target),]

  #adapt to our needs extracting the metabolites:
  meta_network_metabs <- meta_network[grepl("Metab__", meta_network$source) | grepl("Metab__HMDB", meta_network$target),-2]#extract entries with metabolites in source or Target
  meta_network_metabs <- meta_network_metabs[grepl("Gene", meta_network_metabs$source) | grepl("Gene", meta_network_metabs$target),]#extract entries with genes in source or Target

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
  MetabSet <- merge(meta_network_metabs,Input_GeneSet, by.x="gene", by.y=SettingsInfo[["Target"]])

  #combine with pathways --> File that can be used for combined pathway analysis (metabolites and gene t.vals)
  GeneMetabSet <- unique(as.data.frame(rbind(Input_GeneSet%>%dplyr::rename("feature"=SettingsInfo[["Target"]]), MetabSet[,-1]%>%dplyr::rename("feature"=1))))


  ##------------ Select metabolites only
  MetabSet <-  GeneMetabSet %>%
    filter(grepl("HMDB", feature))


  ##-------------- Save and return
  DF_List <- list("GeneMetabSet"=GeneMetabSet,
                  "MetabSet"=MetabSet)
  suppressMessages(suppressWarnings(
    SaveRes(InputList_DF= DF_List,#This needs to be a list, also for single comparisons
                         InputList_Plot= NULL,
                         SaveAs_Table=SaveAs_Table,
                         SaveAs_Plot=NULL,
                         FolderPath= SubFolder,
                         FileName= PKName,
                         CoRe=FALSE,
                         PrintPlot=FALSE)))

  return(invisible(DF_List))
}



##########################################################################################
### ### ### Load MetaLinksDB prior knowledge ### ### ###
##########################################################################################
#'
#' Function to
#'
#' @param types Desired edge types. Options are: "lr", "pd", where 'lr' stands for 'ligand-receptor' and 'pd' stands for 'production-degradation'.\strong{default: NULL}
#' @param cell_location Desired metabolite cell locations. Pass selection using c("Select1", "Select2", "Selectn"). View options setting "?". Options are: "Cytoplasm", "Endoplasmic reticulum", "Extracellular", "Lysosome" , "Mitochondria", "Peroxisome", "Membrane", "Nucleus", "Golgi apparatus" , "Inner mitochondrial membrane". \strong{default: NULL}
#' @param tissue_location Desired metabolite tissue locations. Pass selection using c("Select1", "Select2", "Selectn"). View options setting "?". Options are: "Placenta", "Adipose Tissue","Bladder", "Brain", "Epidermis","Kidney", "Liver", "Neuron", "Pancreas", "Prostate", "Skeletal Muscle", "Spleen", "Testis", "Thyroid Gland", "Adrenal Medulla", "Erythrocyte","Fibroblasts", "Intestine", "Ovary", "Platelet", "All Tissues", "Semen", "Adrenal Gland", "Adrenal Cortex", "Heart", "Lung", "Hair", "Eye Lens", "Leukocyte", Retina", "Smooth Muscle", "Gall Bladder", "Bile",  "Bone Marrow", "Blood", "Basal Ganglia", "Cartilage". \strong{default: NULL}
#' @param biospecimen_location Desired metabolite biospecimen locations.Pass selection using c("Select1", "Select2", "Selectn").View options setting "?".  "Blood", "Feces", "Saliva", "Sweat", "Urine", "Breast Milk", "Cellular Cytoplasm", "Cerebrospinal Fluid (CSF)", "Amniotic Fluid" , "Aqueous Humour", "Ascites Fluid", "Lymph", "Tears", "Breath", "Bile", "Semen", "Pericardial Effusion".\strong{default: NULL}
#' @param disease Desired metabolite diseases.Pass selection using c("Select1", "Select2", "Selectn"). View options setting "?". \strong{default: NULL}
#' @param pathway Desired metabolite pathways.Pass selection using c("Select1", "Select2", "Selectn"). View options setting "?".\strong{default: NULL}
#' @param hmdb_ids Desired HMDB IDs.Pass selection using c("Select1", "Select2", "Selectn"). View options setting "?".\strong{default: NULL}
#' @param uniprot_ids Desired UniProt IDs.Pass selection using c("Select1", "Select2", "Selectn"). View options setting "?".\strong{default: NULL}
#' @param SaveAs_Table \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt". \strong{Default = "csv"}
#' @param FolderPath \emph{Optional:} Path to the folder the results should be saved at. \strong{default: NULL}

#' @export

LoadMetalinks <- function(types = NULL,
                          cell_location = NULL,
                          tissue_location = NULL,
                          biospecimen_location = NULL,
                          disease = NULL,
                          pathway = NULL,
                          hmdb_ids = NULL,
                          uniprot_ids = NULL,
                          SaveAs_Table="csv",
                          FolderPath=NULL){
  ## ------------ Create log file ----------- ##
  MetaProViz_Init()

  logger::log_info("MetaLinksDB.")


  #------------------------------------------------------------------
  if((any(c(types, cell_location, tissue_location, biospecimen_location, disease, pathway, hmdb_ids, uniprot_ids)=="?"))==FALSE){
  #Check Input parameters



    }

  #Python version enables the user to add their own link to the database dump (probably to obtain a specific version. Lets check how the link was generated and see if it would make sense for us to do the same.)
  # --> At the moment arbitrary!
  # We could provide the user the ability to point to their own path were they already dumpled/stored QA version of metalinks they like to use!

  ## ------------ Folder ----------- ##
  if(is.null(SaveAs_Table)==FALSE){
    Folder <- SavePath(FolderName= "PriorKnowledge",
                                    FolderPath=FolderPath)

    SubFolder <- file.path(Folder, "MetabSet")
    if (!dir.exists(SubFolder)) {dir.create(SubFolder)}
  }

  #------------------------------------------------------------------
  #Get the directory and filepath of cache results of R
  directory <- rappdirs::user_cache_dir()#get chache directory
  File_path <-paste(directory, "/metalinks.db", sep="")

  if(file.exists(File_path)==TRUE){# First we will check the users chache directory and weather there are rds files with KEGG_pathways already:
    # Connect to the SQLite database
    con <- DBI::dbConnect(RSQLite::SQLite(), File_path, synchronous = NULL)
    message("Cached file loaded from: ", File_path)
  }else{# load from API
    RequiredPackages <- c("tidyverse", "RSQLite", "DBI")
    new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
    if(length(new.packages)) install.packages(new.packages)

    suppressMessages(library(tidyverse))

    #--------------------------------------------------------------------------------------------
    #Python availability via Liana: https://github.com/saezlab/liana-py/blob/main/liana/resource/get_metalinks.py
    metalinks_db_url <- "https://figshare.com/ndownloader/files/47567597"
    # Download the Metalinks database file and save where the cache is stored
    download.file(metalinks_db_url, destfile =  File_path , mode = "wb")#WB: This mode is used for writing binary files. It opens the destination file for writing in binary mode.
    message("Metalinks database downloaded and saved to: ", File_path)

    # Connect to the SQLite database
    con <- DBI::dbConnect(RSQLite::SQLite(), File_path, synchronous = NULL)
  }

  #------------------------------------------------------------------
  #Query the database for a specific tables
  tables <- DBI::dbListTables(con)

  TablesList <- list()
  for(table in tables){
    query <- paste("SELECT * FROM", table)
    data <- DBI::dbGetQuery(con, query)
    TablesList[[table]] <- data
  }

  # Close the connection
  DBI::dbDisconnect(con)

  MetalinksDB <- TablesList[["edges"]]#extract the edges table
  #------------------------------------------------------------------
  # Answer questions about the database
  #if any parameter is ? then return the data
  if(any(c(types, cell_location, tissue_location, biospecimen_location, disease, pathway, hmdb_ids, uniprot_ids)=="?")){
    Questions <- which(c(types, cell_location, tissue_location, biospecimen_location, disease, pathway, hmdb_ids, uniprot_ids)=="?")
    #Check tables where the user has questions
    if(length(Questions)>0){
      for(i in Questions){
        if(i==1){
          print("Types:")
          print(unique(MetalinksDB$type))
        }
        if(i==2){
          print("Cell Location:")
          CellLocation <- TablesList[["cell_location"]]
          print(unique(CellLocation$cell_location))
        }
        if(i==3){
          print("Tissue Location:")
          TissueLocation <- TablesList[["tissue_location"]]
          print(unique(TissueLocation$tissue_location))
        }
        if(i==4){
          print("Biospecimen Location:")
          BiospecimenLocation <- TablesList[["biospecimen_location"]]
          print(unique(BiospecimenLocation$biospecimen_location))
        }
        if(i==5){
          print("Disease:")
          Disease <- TablesList[["disease"]]
          print(unique(Disease$disease))
        }
        if(i==6){
          print("Pathway:")
          Pathway <- TablesList[["pathway"]]
          print(unique(Pathway$pathway))
        }
        if(i==7){
          print("HMDB IDs:")
          print(unique(MetalinksDB$hmdb))
        }
        if(i==8){
          print("UniProt IDs:")
          print(unique(MetalinksDB$uniprot))
        }
      }
    }
    message("No result is returned unless correct options for your selections are used. `?` is not a valid option, but only returns you the list of options.")
    return()
  }

  #------------------------------------------------------------------
  #Extract specific connections based on parameter settings. If any parameter is not NULL, filter the data:
  ## types
  if(!is.null(types)){
    MetalinksDB <- MetalinksDB[MetalinksDB$type %in% types,]
  }

  ## cell_location
  if(!is.null(cell_location)){
    CellLocation <- TablesList[["cell_location"]]
    CellLocation <- CellLocation[CellLocation$cell_location %in% cell_location,]#Filter the cell location

    #Get unique HMDB IDs
    CellLocation_HMDB <- unique(CellLocation$hmdb)

    #Only keep selected HMDB IDs
    MetalinksDB <- MetalinksDB[MetalinksDB$hmdb %in%  CellLocation_HMDB,]
  }

  ## tissue_location
  if(!is.null(tissue_location)){# "All Tissues"?
    TissueLocation <- TablesList[["tissue_location"]]
    TissueLocation <- TissueLocation[TissueLocation$tissue_location %in% tissue_location,]#Filter the tissue location

    #Get unique HMDB IDs
    TissueLocation_HMDB <- unique(TissueLocation$hmdb)

    #Only keep selected HMDB IDs
    MetalinksDB <- MetalinksDB[MetalinksDB$hmdb %in%  TissueLocation_HMDB,]
  }

  ## biospecimen_location
  if(!is.null(biospecimen_location)){
    BiospecimenLocation <- TablesList[["biospecimen_location"]]
    BiospecimenLocation <- BiospecimenLocation[BiospecimenLocation$biospecimen_location %in% biospecimen_location,]#Filter the biospecimen location

    #Get unique HMDB IDs
    BiospecimenLocation_HMDB <- unique(BiospecimenLocation$hmdb)

    #Only keep selected HMDB IDs
    MetalinksDB <- MetalinksDB[MetalinksDB$hmdb %in%  BiospecimenLocation_HMDB,]
  }

  ## disease
  if(!is.null(disease)){
    Disease <- TablesList[["disease"]]
    Disease <- Disease[Disease$disease %in% disease,]#Filter the disease

    #Get unique HMDB IDs
    Disease_HMDB <- unique(Disease$hmdb)

    #Only keep selected HMDB IDs
    MetalinksDB <- MetalinksDB[MetalinksDB$hmdb %in%  Disease_HMDB,]
  }

  ## pathway
  if(!is.null(pathway)){
    Pathway <- TablesList[["pathway"]]
    Pathway <- Pathway[Pathway$pathway %in% pathway,]#Filter the pathway

    #Get unique HMDB IDs
    Pathway_HMDB <- unique(Pathway$hmdb)

    #Only keep selected HMDB IDs
    MetalinksDB <- MetalinksDB[MetalinksDB$hmdb %in%  Pathway_HMDB,]
  }

  ## hmdb_ids
  if(!is.null(hmdb_ids)){
    #Only keep selected HMDB IDs
    MetalinksDB <- MetalinksDB[MetalinksDB$hmdb %in%  hmdb_ids,]
  }

  ## uniprot_ids
  if(!is.null(uniprot_ids)){
    #Only keep selected UniProt IDs
    MetalinksDB <- MetalinksDB[MetalinksDB$uniprot %in%  uniprot_ids,]
  }

  #------------------------------------------------------------------
  #Add other ID types:
  ## Metabolite Name
  MetalinksDB <- merge(MetalinksDB, TablesList[["metabolites"]], by="hmdb", all.x=TRUE)

  ## Gene Name
  MetalinksDB <- merge(MetalinksDB, TablesList[["proteins"]], by="uniprot", all.x=TRUE)


  ## Rearrange columns:
  MetalinksDB <- MetalinksDB[,c(2,10:12, 1, 13:14, 3:9)]%>%
    dplyr::mutate(type = dplyr::case_when(
      type == "lr" ~ "Ligand-Receptor",
      type == "pd" ~ "Production-Degradation",
      TRUE ~ type  # this keeps the original value if it doesn't match any condition
    ))%>%
    dplyr::mutate(mode_of_regulation = dplyr::case_when(
    mor == -1 ~ "Inhibiting",
    mor == 1 ~ "Activating",
    mor == 0 ~ "Binding",
    TRUE ~ as.character(mor)  # this keeps the original value if it doesn't match any condition
  ))
  #--------------------------------------------------------
  #Remove metabolites that are not detectable by mass spectrometry


  #------------------------------------------------------------------
  #Decide on useful selections term-metabolite for MetaProViz.
  #MetalinksDB_Pathways <- merge(MetalinksDB[, c(1:3)], TablesList[["pathway"]], by="hmdb", all.x=TRUE)

  MetalinksDB_Type <-MetalinksDB[, c(1:3, 7,13)]%>%
    subset(!is.na(protein_type))%>%
    dplyr::mutate(term = gsub('\"', '', protein_type))%>%
    unite(term_specific, c("term", "type"), sep = "_", remove = FALSE)%>%
    dplyr::select(-"type", -"protein_type")

  #------------------------------------------------------------------
  #Save results in folder
  ##-------------- Save and return
  DF_List <- list("MetalinksDB"=MetalinksDB,
                  "MetalinksDB_Type"=MetalinksDB_Type)
  suppressMessages(suppressWarnings(
    SaveRes(InputList_DF= DF_List,#This needs to be a list, also for single comparisons
                         InputList_Plot= NULL,
                         SaveAs_Table=SaveAs_Table,
                         SaveAs_Plot=NULL,
                         FolderPath= SubFolder,
                         FileName= "MetaLinksDB",
                         CoRe=FALSE,
                         PrintPlot=FALSE)))

  #Return into environment
  return(invisible(DF_List))

}

##########################################################################################
### ### ### Helper: Load compound lists of ions, xenobiotics, cofactors (not detectable by LC-MS), etc ### ### ###
##########################################################################################

# This is needed to remove ions, xenobiotics, cofactors, etc. from the metabolite list for any of the functions above.

