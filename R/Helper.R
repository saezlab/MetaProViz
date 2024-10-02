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
#' intra <- ToyData("IntraCells_Raw")
#'
#' @importFrom readr read_csv cols
#' @importFrom magrittr %>% extract2
#' @importFrom tibble column_to_rownames
#' @export
ToyData <- function(Dataset) {

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
    read_csv(col_types = cols()) %>%
    {`if`(
    (rncol <- names(.) %>% intersect(rncols)) %>% length,
    column_to_rownames(., rncol),
    .
    )}

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
#' @keywords Create folder and path
#' @noRd
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


#' Make sure the results directory exists
#'
#' @importFrom magrittr %>%
#' @noRd
ResultsDir <- function(path = 'MetaProViz_Results') {

  # TODO: options?
  path %>% {`if`(!dir.exists(.), {dir.create(.); .}, .) }

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
#' @keywords Save
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
      if(CoRe==FALSE | is.null(CoRe)==TRUE){
        FileName <- paste0(FolderPath,"/" , FileName, "_",Sys.Date(), sep = "")
      }else{
        FileName <- paste0(FolderPath,"/CoRe_" , FileName,"_",Sys.Date(), sep = "")
      }
      #Save Excel
      writexl::write_xlsx(InputList_DF, paste0(FileName,".xlsx", sep = "") , col_names = TRUE, rownames = FALSE)
    }else{
      for(DF in names(InputList_DF)){
        #Make FileName
        if(CoRe==FALSE | is.null(CoRe)==TRUE){
          FileName_Save <- paste0(FolderPath,"/" , FileName, "_", DF ,"_",Sys.Date(), sep = "")
        }else{
          FileName_Save <- paste0(FolderPath,"/CoRe_" , FileName, "_", DF ,"_",Sys.Date(), sep = "")
        }

        #Save table
        if (SaveAs_Table == "csv"){
          write.csv(InputList_DF[[DF]], paste0(FileName_Save,".csv", sep = ""), row.names = FALSE)
        }else if (SaveAs_Table == "txt"){
          write.table(InputList_DF[[DF]], paste0(FileName_Save,".txt", sep = "") , col.names = TRUE, row.names = FALSE)
        }
      }
    }
  }

  if(is.null(SaveAs_Plot)==FALSE){
    for(Plot in names(InputList_Plot)){
      #Make FileName
      if(CoRe==FALSE | is.null(CoRe)==TRUE){
        FileName_Save <- paste0(FolderPath,"/" , FileName,"_", Plot , "_",Sys.Date(), sep = "")
      }else{
        FileName_Save <- paste0(FolderPath,"/CoRe_" , FileName,"_", Plot ,"_",Sys.Date(), sep = "")
      }

      #Save
      if(is.null(PlotHeight)){
        PlotHeight <- 12
      }
      if(is.null(PlotWidth)){
        PlotWidth <- 16
      }
      if(is.null(PlotUnit)){
        PlotUnit <- "cm"
      }

      ggsave(filename = paste0(FileName_Save, ".",SaveAs_Plot, sep=""), plot = InputList_Plot[[Plot]], width = PlotWidth,  height = PlotHeight, unit=PlotUnit)

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
# REFACT: Description of each argument should start with its type; e.g.
# "@param x Character: name of the variable mapped to the x axis."
#' @param InputData Passed to main function MetaProViz::Function()
#' @param InputData_Num  \emph{Optional: } If InputData must be numeric \strong{Default = TRUE}
#' @param SettingsFile_Sample Passed to main function MetaProViz::Function()
#' @param SettingsFile_Metab Passed to main function MetaProViz::Function(). If not avaliable can be set to NULL.
#' @param SettingsInfo Passed to main function MetaProViz::Function()
#' @param SaveAs_Plot Passed to main function MetaProViz::Function(). If not avaliable can be set to NULL.
#' @param SaveAs_Table Passed to main function MetaProViz::Function(). If not avaliable can be set to NULL.
#' @param CoRe \emph{Optional: } Passed to main function MetaProViz::Function(). If not avaliable can be set to NULL. \strong{Default = FALSE}
#' @param PrintPlot Passed to main function MetaProViz::Function(). If not avaliable can be set to NULL.
#' @param Theme \emph{Optional: } Passed to main function MetaProViz::Function(). If not avaliable can be set to NULL.  \strong{Default = NULL}
#' @param PlotSettings \emph{Optional: } Needs to be set for MetaProViz::VizX functions. Options are "Sample", "Feature", Both". This refers to SettingsInfo color, shape, individual as for some plots we have both feature and sample settings. \strong{Default = NULL}
#'
#'
#' @keywords Input check
#' @noRd
#'
#'

CheckInput <- function(InputData,
                       InputData_Num=TRUE,
                       SettingsFile_Sample,
                       SettingsFile_Metab,
                       SettingsInfo,
                       SaveAs_Plot,
                       SaveAs_Table,
                       CoRe=FALSE,
                       PrintPlot,
                       Theme=NULL,
                       PlotSettings=NULL){
  ############## Parameters valid for multiple MetaProViz functions

  #-------------InputData
  # REFACT: this will fail with "condition has length > 1" error, and also
  # gives wrong result. Use is.data.frame() instead.
  if(class(InputData) != "data.frame"){
    stop("InputData should be a data.frame. It's currently a ", paste(class(InputData), ".",sep = ""))
  }
  if(any(duplicated(row.names(InputData)))==TRUE){
    stop("Duplicated row.names of InputData, whilst row.names must be unique")
  }

  if(InputData_Num==TRUE){
     Test_num <- apply(InputData, 2, function(x) is.numeric(x))
     if((any(Test_num) ==  FALSE) ==  TRUE){
       stop("InputData needs to be of class numeric")
       }
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
      stop("col.names InputData need to match row.names SettingsFile_Metab.")
    }
  }

  #-------------SettingsInfo
  if(is.vector(SettingsInfo)==FALSE & is.null(SettingsInfo)==FALSE){
    stop("SettingsInfo should be NULL or a vector. It's currently a ", paste(class(SettingsInfo), ".", sep = ""))
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

    #Superplot
    if("Superplot" %in% names(SettingsInfo)){
      if(SettingsInfo[["Superplot"]] %in% colnames(SettingsFile_Sample)== FALSE){
        stop("The ",SettingsInfo[["Superplot"]], " column selected as Superplot column in SettingsInfo was not found in SettingsFile_Sample. Please check your input.")
      }
    }

    if(is.null(PlotSettings)==FALSE){
      if(PlotSettings== "Sample"){
        #Plot colour
        if("color" %in% names(SettingsInfo)){
          if(SettingsInfo[["color"]] %in% colnames(SettingsFile_Sample)== FALSE){
            stop("The ",SettingsInfo[["color"]], " column selected as color in SettingsInfo was not found in SettingsFile_Sample. Please check your input.")
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
      }else if(PlotSettings== "Feature"){
        if("color" %in% names(SettingsInfo)){
          if(SettingsInfo[["color"]] %in% colnames(SettingsFile_Metab)== FALSE){
            stop("The ",SettingsInfo[["color"]], " column selected as color in SettingsInfo was not found in SettingsFile_Metab. Please check your input.")
          }
        }

        #Plot shape
        if("shape" %in% names(SettingsInfo)){
          if(SettingsInfo[["shape"]] %in% colnames(SettingsFile_Metab)== FALSE){
            stop("The ",SettingsInfo[["shape"]], " column selected as shape in SettingsInfo was not found in SettingsFile_Metab. Please check your input.")
          }
        }

        #Plot individual
        if("individual" %in% names(SettingsInfo)){
          if(SettingsInfo[["individual"]] %in% colnames(SettingsFile_Metab)== FALSE){
            stop("The ",SettingsInfo[["individual"]], " column selected as individual in SettingsInfo was not found in SettingsFile_Metab. Please check your input.")
          }
        }
      }else if(PlotSettings== "Both"){
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

       # Plot shape_metab
        if("shape_Metab" %in% names(SettingsInfo)){
          if(SettingsInfo[["shape_Metab"]] %in% colnames(SettingsFile_Metab)== FALSE){
            stop("The ",SettingsInfo[["shape_Metab"]], " column selected as shape_Metab in SettingsInfo was not found in SettingsFile_Metab. Please check your input.")
          }
        }

        # Plot shape_metab
        if("shape_Sample" %in% names(SettingsInfo)){
          if(SettingsInfo[["shape_Sample"]] %in% colnames(SettingsFile_Metab)== FALSE){
            stop("The ",SettingsInfo[["shape_Sample"]], " column selected as shape_Metab in SettingsInfo was not found in SettingsFile_Sample. Please check your input.")
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
  if(is.logical(CoRe) == FALSE){
    stop("Check input. The CoRe value should be either =TRUE for preprocessing of Consuption/Release experiment or =FALSE if not.")
  }

  #-------------Theme
  if(is.null(Theme)==FALSE){
    Theme_options <- c("theme_grey()", "theme_gray()", "theme_bw()", "theme_linedraw()", "theme_light()", "theme_dark()", "theme_minimal()", "theme_classic()", "theme_void()", "theme_test()")
    if (Theme %in% Theme_options == FALSE){
      stop("Theme option is incorrect. You can check for complete themes here: https://ggplot2.tidyverse.org/reference/ggtheme.html. Options are the following: ",paste(Theme_options, collapse = ", "),"." )
    }
  }
  #------------- general
  if(is.logical(PrintPlot) == FALSE){
    stop("Check input. PrintPlot should be either =TRUE or =FALSE.")
  }
}
