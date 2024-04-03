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


#' Imports toy data into environment
#'
#' @param data Either "Standard", "Standard_DMA", "CoRe" or "MappingInfo" depending which data you would like to load
#' @title Toy Data Import
#' @description Import and process .csv file to create toy data.
#' @importFrom utils read.csv
#' @return A data frame containing the toy data.
#' @export
#'
ToyData <- function(data) {
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


################################################################################################
### ### ### General helper function: Save folder path ### ### ###
################################################################################################

#' SavePath is the helper function to create the folder structure and path
#'
#' @param FolderName Name of the folder, which can not contain any special characters. Created within  the individual MetaProViz functions and can not be changed by the user.
#' @param FolderPath Passed to main function by the user
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


################################################################################################
### ### ### General helper function: Save tables and plot ### ### ###
################################################################################################

#' SaveRes is the helper function to create the folder structure and path
#'
#' @param InputList_DF Generated within the MetaProViz function. Contains named DFs.
#' @param InputList_Plot Generated within the MetaProViz function. Contains named Plots.
#' @param SaveAs_Table Passed to main function by the user. If not avalailable can be set to NULL.
#' @param SaveAs_Plot Passed to main function by the user.If not avalailable can be set to NULL.
#' @param FolderPath Passed to main function by the user
#' @param FileName Passed to main function by the user
#' @param CoRe Passed to main function by the user. If not avalailable can be set to NULL.
#' @param PrintPlot Passed to main function by the user. If not avalailable can be set to NULL.
#' @param PlotHeight Parameter for ggsave.
#' @param PlotWidth Parameter for ggsave.
#' @param PlotUnit Parameter for ggsave.
#'
#' @keywords
#' @noRd
#'

SaveRes<- function(InputList_DF,
                   InputList_Plot,
                   SaveAs_Table,
                   SaveAs_Plot=SaveAs_Plot,
                   FolderPath,
                   FileName,
                   CoRe=FALSE,
                   PrintPlot=TRUE,
                   PlotHeight=NULL,
                   PlotWidth=NULL,
                   PlotUnit=NULL){

  RequiredPackages <- c("tidyverse", "writexl")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  suppressMessages(library(tidyverse))


  if(is.null(SaveAs_Table)==FALSE){
    # Excel File: One file with multiple sheets:
    if(SaveAs_Table == "xlsx"){
      #Make FileName
      if(CoRe==FALSE){
        FileName <- paste0(FolderPath,"/" , FileName, "_",Sys.Date(), sep = "")
      }else{
        FileName <- paste0(FolderPath,"/CoRe_" , FileName,"_",Sys.Date(), sep = "")
      }
      #Save Excel
      writexl::write_xlsx(InputList_DF, paste0(FileName,".xlsx", sep = "") , col_names = TRUE)
    }else{
      for(DF in names(InputList_DF)){
        #Make FileName
        if(CoRe==FALSE){
          FileName <- paste0(FolderPath,"/" , FileName, "_", DF ,"_",Sys.Date(), sep = "")
        }else{
          FileName <- paste0(FolderPath,"/CoRe_" , FileName, "_", DF ,"_",Sys.Date(), sep = "")
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
        FileName <- paste0(FolderPath,"/" , FileName,"_", Plot , "_",Sys.Date(), sep = "")
      }else{
        FileName <- paste0(FolderPath,"/CoRe_" , FileName,"_", Plot ,"_",Sys.Date(), sep = "")
      }

      #Save
      if(is.null(PlotHeight)){
        PlotHeight <- 8
      }
      if(is.null(PlotWidth)){
        PlotWidth <- 10
      }
      if(is.null(PlotUnit)){
        PlotUnit <- "cm"
      }

      ggsave(filename = paste0(FileName, ".",SaveAs_Plot, sep=""), plot = InputList_Plot[[Plot]], width = PlotWidth,  height = PlotHeight, unit=PlotUnit)

      if(PrintPlot==TRUE){
        suppressMessages(suppressWarnings(plot(InputList_Plot[[Plot]])))
      }
    }
  }
}



################################################################################################
### ### ### Helper function: Internal Function to check function input ### ### ###
################################################################################################

#' Check input parameters
#'
#' @param InputData Passed to main function MetaProViz::PreProcessing()
#' @param SettingsFile_Sample Passed to main function MetaProViz::PreProcessing()
#' @param SettingsFile_Metab Passed to main function MetaProViz::PreProcessing(). If not avaliable can be set to NULL.
#' @param SettingsInfo Passed to main function MetaProViz::PreProcessing()
#' @param SaveAs_Plot Passed to main function MetaProViz::PreProcessing(). If not avaliable can be set to NULL.
#' @param SaveAs_Table Passed to main function MetaProViz::PreProcessing(). If not avaliable can be set to NULL.
#' @param CoRe Passed to main function MetaProViz::PreProcessing(). If not avaliable can be set to NULL.
#' @param PrintPlot Passed to main function MetaProViz::PreProcessing(). If not avaliable can be set to NULL.
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
  RequiredPackages <- c("tidyverse")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  suppressMessages(library(tidyverse))

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

  if(sum(duplicated(colnames(InputData))) > 0){
    doublons <- as.character(colnames(InputData)[duplicated(colnames(InputData))])#number of duplications
    #data <-data[!duplicated(colnames(InputData)),]#remove duplications
    stop("InputData contained duplicates column names, whilst col.names must be unique.")
  }

  #-------------SettingsFile
  if(is.null(SettingsFile_Sample)==FALSE){
    Test_match <- merge(SettingsFile_Sample, InputData, by = "row.names", all =  FALSE)
    if(nrow(Test_match) ==  0){
        stop("row.names InputData need to match row.names SettingsFile_Sample.")
      }
  }

  if(is.null(SettingsFile_Metab)==FALSE){
    Test_match <- merge(SettingsFile_Metab, as.data.frame(t(InputData)), by = "row.names", all =  FALSE)
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

  if(is.null(SettingsInfo)==FALSE){
    #Conditions
    if("Conditions" %in% names(SettingsInfo)){
      if(SettingsInfo[["Conditions"]] %in% colnames(SettingsFile_Sample)== FALSE){
        stop("The ", SettingsInfo[["Conditions"]], " column selected as Conditions in SettingsInfo was not found in SettingsFile. Please check your input.")
      }
    }

    #Biological replicates
    if("Biological_Replicates" %in% names(SettingsInfo)){
      if(SettingsInfo[["Biological_Replicates"]] %in% colnames(SettingsFile_Sample)== FALSE){
        stop("The ",SettingsInfo[["Biological_Replicates"]], " column selected as Biological_Replicates in SettingsInfo was not found in SettingsFile_Sample. Please check your input.")
      }
    }

    #Numerator
    if("Numerator" %in% names(SettingsInfo)==TRUE){
      if(SettingsInfo[["Numerator"]] %in% SettingsFile_Sample[[SettingsInfo[["Conditions"]]]]==FALSE){
        stop("The ",SettingsInfo[["Numerator"]], " column selected as numerator in SettingsInfo was not found in SettingsFile_Sample. Please check your input.")
      }
    }

   #Denominator
    if("Denominator" %in% names(SettingsInfo)==TRUE){
      if(SettingsInfo[["Denominator"]] %in% SettingsFile_Sample[[SettingsInfo[["Conditions"]]]]==FALSE){
        stop("The ",SettingsInfo[["Denominator"]], " column selected as denominator in SettingsInfo was not found in SettingsFile_Sample. Please check your input.")
      }
    }

    #Denominator & Numerator
    if("Denominator" %in% names(SettingsInfo)==FALSE  & "Numerator" %in% names(SettingsInfo) ==TRUE){
      stop("Check input. The selected denominator option is empty while ",paste(SettingsInfo[["Numerator"]])," has been selected as a numerator. Please add a denominator for 1-vs-1 comparison or remove the numerator for all-vs-all comparison." )
    }

    #Plot colour
    if("color" %in% names(SettingsInfo)){
      if(SettingsInfo[["color"]] %in% colnames(SettingsFile_Sample)== FALSE){
        stop("The ",SettingsInfo[["color"]], " column selected as color in SettingsInfo was not found in SettingsFile_Sample. Please check your input.")
      }
    }

    #Plot colour sample
    if("color_Sample" %in% names(SettingsInfo)){
      if(SettingsInfo[["color_Sample"]] %in% colnames(SettingsFile_Sample)== FALSE){
        stop("The ",SettingsInfo[["color_Sample"]], " column selected as color_Sample in SettingsInfo was not found in SettingsFile_Sample. Please check your input.")
      }
    }

    #Plot colour Metab
    if("color_Metab" %in% names(SettingsInfo)){
      if(SettingsInfo[["color_Metab"]] %in% colnames(SettingsFile_Metab)== FALSE){
        stop("The ",SettingsInfo[["color_Metab"]], " column selected as color_Metab in SettingsInfo was not found in SettingsFile_Metab. Please check your input.")
      }
      if(sum(colnames(InputData) %in% SettingsFile_Metab$Metabolite) < length(InputData)  ){
        warning("The InputData contains metabolites not found in SettingsFile_Metab.")
      }
    }

    #Plot shape
    if("shape" %in% names(SettingsInfo)){
      if(SettingsInfo[["shape"]] %in% colnames(SettingsFile_Sample)== FALSE){
        stop("The ",SettingsInfo[["shape"]], " column selected as shape in SettingsInfo was not found in SettingsFile_Sample. Please check your input.")
      }
    }

    #Plot individual
    if("individual" %in% names(SettingsInfo)){
      if(SettingsInfo[["individual"]] %in% colnames(SettingsFile_Sample)== FALSE){
        stop("The ",SettingsInfo[["individual"]], " column selected as individual in SettingsInfo was not found in SettingsFile_Sample. Please check your input.")
      }
    }

    #Plot individual_Metab
    if("individual_Metab" %in% names(SettingsInfo)){
      if(SettingsInfo[["individual_Metab"]] %in% colnames(SettingsFile_Metab)== FALSE){
        stop("The ",SettingsInfo[["individual_Metab"]], " column selected as individual_Metab in SettingsInfo was not found in SettingsFile_Metab. Please check your input.")
      }
    }

    #Plot individual_Sample
    if("individual_Sample" %in% names(SettingsInfo)){
      if(SettingsInfo[["individual_Sample"]] %in% colnames(SettingsFile_Sample)== FALSE){
        stop("The ",SettingsInfo[["individual_Sample"]], " column selected as individual_Sample in SettingsInfo was not found in SettingsFile_Sample. Please check your input.")
      }
    }



  }

  #-------------SaveAs
  Save_as_Plot_options <- c("svg","pdf", "png")
  if(is.null(SaveAs_Plot)==FALSE){
    if(SaveAs_Plot %in% Save_as_Plot_options == FALSE){
    stop("Check input. The selected SaveAs_Plot option is not valid. Please select one of the folowwing: ",paste(Save_as_Plot_options,collapse = ", ")," or set to NULL if no plots should be saved." )
  }
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
