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
#' @param data Either "Standard", "Standard_DMA", "CoRe" or "MappingInfo" depending which data you would like to load
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

  Intra_DMA <- system.file("data", "MS55_DMA_786M1A_vs_HK2.csv", package = "MetaProViz")
  Intra_DMA<- read.csv(Intra_DMA, check.names=FALSE)

  Media <- system.file("data", "MS51_RawPeakData.csv", package = "MetaProViz")
  Media<- read.csv(Media, check.names=FALSE)

  Pathways <-system.file("data", "MappingTable_SelectPathways.csv", package = "MetaProViz")
  Pathways<- read.csv(Pathways, check.names=FALSE)

  # Return the toy data into environment
  if(data=="Standard"){
    assign("Intra", Intra, envir=.GlobalEnv)
  } else if(data=="Standard_DMA"){
    assign("Intra_DMA_786M1A_vs_HK2",  Intra_DMA, envir=.GlobalEnv)
  }else if(data=="CoRe"){
    assign("Media", Media, envir=.GlobalEnv)
  } else if(data=="MappingInfo"){
    assign("MappingInfo", Pathways, envir=.GlobalEnv)
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
MCA_rules<- function(data){
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
#' @title KEGG
#' @description Import and process KEGG.
#' @importFrom utils read.csv
#' @return A data frame containing the KEGG pathways for ORA.
#' @export
#'
#'
Load_KEGG<- function(){
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



#' SavePath is the helper function to create the folder structure and path
#'
#' @param FolderName Name of the folder, which can not contain any special characters.
#'
#'
#' @keywords
#' @noRd
#'

SavePath<- function(FolderName, FolderPath){

  #Check if FolderName includes special characters that are not allowed
  cleaned_FolderName <- gsub("[^a-zA-Z0-9 ]", "", FolderName)
  if (FolderName != cleaned_FolderName){
    message("Special characters were removed from `FolderName`.")
  }

  #Check if FolderPath exist
  if(is.null(FolderPath)){
    FolderPath <- getwd()
    FolderPath <- file.path(FolderPath, "MetaProViz_Results")
    if(!dir.exists(FolderPath)){dir.create(FolderPath)}
  }else{
    if(dir.exists(FolderPath)==FALSE){
      FolderPath <- getwd()
      message("Provided `FolderPath` does not exist and hence results are saved here: ", FolderPath, sep="")
    }
  }

  #Create the folder name
  Results_folder <- file.path(FolderPath, cleaned_FolderName)
  if(!dir.exists(Results_folder)){dir.create(Results_folder)}

  #Return the folder path:
  return(invisible(Results_folder))
}


#' SaveRes
#'
#' @param FolderPath
#' @param CoRe
#' @param List TRUE or FALSE
#'
#' @keywords
#' @noRd
#'

SaveRes<- function(InputList_DF, InputList_Plot, SaveAs_Table,SaveAs_Plot, FolderPath, OutputFileName, CoRe, PrintPlot){

  if(is.null(SaveAs_Table)==FALSE){
    # Excel File: One file with multiple sheets:
    if(SaveAs_Table == "xlsx"){
      #Make FileName
      if(CoRe==FALSE){
        FileName <- paste0(FolderPath,"/" , OutputFileName, "_",Sys.Date(), sep = "")
      }else{
        FileName <- paste0(FolderPath,"/CoRe_" , OutputFileName,"_",Sys.Date(), sep = "")
      }
      #Save Excel
      writexl::write_xlsx(InputList_DF, paste0(FileName,".xlsx", sep = "") , col_names = TRUE)
    }else{
      for(DF in names(InputList_DF)){
        #Make FileName
        if(CoRe==FALSE){
          FileName <- paste0(FolderPath,"/" , OutputFileName, "_", DF ,"_",Sys.Date(), sep = "")
        }else{
          FileName <- paste0(FolderPath,"/CoRe_" , OutputFileName, "_", DF ,"_",Sys.Date(), sep = "")
        }

        #Save table
        if (SaveAs_Table == "csv"){
          write.csv(InputList_DF[[DF]], paste0(FileName,".csv", sep = ""))
          }else if (SaveAs_Table == "txt"){
            write.table(InputList_DF[[DF]], paste0(FileName,".txt", sep = "") , col.names = TRUE, row.names = FALSE)
          }
        }
    }
  }

  if(is.null(SaveAs_Plot)==FALSE){
    for(Plot in names(InputList_Plot)){
      #Make FileName
      if(CoRe==FALSE){
        FileName <- paste0(FolderPath,"/" , OutputFileName,"_", Plot , "_",Sys.Date(), sep = "")
      }else{
        FileName <- paste0(FolderPath,"/CoRe_" , OutputFileName,"_", Plot ,"_",Sys.Date(), sep = "")
      }

      #Save
      ggsave(filename = paste0(FileName, ".",SaveAs_Plot, sep=""), plot = InputList_Plot[[Plot]], width = 10,  height = 8)

      if(PrintPlot==TRUE){
        suppressMessages(suppressWarnings(plot(InputList_Plot[[Plot]])))
      }
    }
  }
}



#' Check input parameters
#'
#' @param Function Name of the MetaProViz Function that is checked.
#' @param InputList
#'
#'
#' @keywords
#' @noRd
#'
#'

CheckInput <- function(InputData,
                       SettingsFile_Sample,
                       SettingsFile_Metab,
                       SettingsInfo,
                       SaveAs_Plot,
                       SaveAs_Table,
                       CoRe,
                       PrintPlot){
  ############## Parameters valid for multiple MetaProViz functions
  #-------------InputData
  if(class(InputData) != "data.frame"){
    stop("InputData should be a data.frame. It's currently a ", paste(class(InputData), ".",sep = ""))
  }
  if(any(duplicated(row.names(InputData)))==TRUE){
    stop("Duplicated row.names of InputData, whilst row.names must be unique")
  }

  Test_num <- apply(InputData, 2, function(x) is.numeric(x))
  if((any(Test_num) ==  FALSE) ==  TRUE){
      stop("InputData needs to be of class numeric")
  }

  #-------------SettingsFile
  if(is.null(SettingsFile_Sample)==FALSE){
    Test_match <- merge(SettingsFile_Sample, InputData, by = "row.names", all =  FALSE)
    if(nrow(Test_match) ==  0){
        stop("row.names InputData need to match row.names SettingsFile_Sample.")
      }
  }

  if(is.null(SettingsFile_Metab)==FALSE){
    Test_match <- merge(SettingsFile_Metab, as.data.frame(InputData), by = "row.names", all =  FALSE)
    if(nrow(Test_match) ==  0){
      stop("col.names InputData need to match col.names SettingsFile_Metab.")
    }
  }

  #-------------SettingsInfo
  if(is.vector(SettingsInfo)==FALSE & is.null(SettingsInfo)==FALSE){
    stop("SettingsInfo should be NULL or a vector. It's currently a ", paste(class(SettingsInfo), ".", sep = ""))
  }

  if((is.null(SettingsInfo)==TRUE & is.null(SettingsFile_Sample)==TRUE) & (is.null(SettingsInfo)==TRUE & is.null(SettingsFile_Metab)==TRUE)){
    message("No Input_Settings have been added.")
  }



  #-------------SaveAs
  Save_as_Plot_options <- c("svg","pdf", "png")
  if((SaveAs_Plot %in% Save_as_Plot_options == FALSE) | (is.null(SaveAs_Plot)==TRUE)){
    stop("Check input. The selected SaveAs_Plot option is not valid. Please select one of the folowwing: ",paste(Save_as_Plot_options,collapse = ", ")," or set to NULL if no plots should be saved." )
  }

  SaveAs_Table_options <- c("txt","csv", "xlsx", "RData")#RData = SummarizedExperiment (?)
  if(is.null(SaveAs_Table)==FALSE){
    if((SaveAs_Table %in% SaveAs_Table_options == FALSE)| (is.null(SaveAs_Table)==TRUE)){
      stop("Check input. The selected SaveAs_Table option is not valid. Please select one of the folowwing: ",paste(SaveAs_Table_options,collapse = ", "),"." )
    }
  }

  #-------------CoRe
  if(is.logical(CoRe) == FALSE | (is.null(CoRe)==TRUE)){
    stop("Check input. The CoRe value should be either =TRUE for preprocessing of Consuption/Release experiment or =FALSE if not.")
  }

  #------------- general
  if(is.logical(PrintPlot) == FALSE){
    stop("Check input. PrintPlot should be either =TRUE or =FALSE.")
  }
}
