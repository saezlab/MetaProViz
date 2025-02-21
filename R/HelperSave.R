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
### ### ### General helper function: Save folder path ### ### ###
################################################################################################

#' SavePath is the helper function to create the folder structure and path
#'
#' @param FolderName Name of the folder, which can not contain any special characters. Created within  the individual MetaProViz functions and can not be changed by the user.
#' @param FolderPath Passed to main function by the user
#'
#' @description Create folder and path
#'
#' @keywords folder, path
#'
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

#' SaveRes is the helper function to save the plots and tables
#'
#' @param InputList_DF \emph{Optional: } Generated within the MetaProViz function. Contains named DFs. If not avalailable can be set to NULL.\strong{Default = NULL}
#' @param InputList_Plot \emph{Optional: } Generated within the MetaProViz function. Contains named Plots. If not avalailable can be set to NULL.\strong{Default = NULL}
#' @param SaveAs_Table \emph{Optional: } Passed to main function by the user. If not avalailable can be set to NULL.\strong{Default = NULL}
#' @param SaveAs_Plot \emph{Optional: } Passed to main function by the user. If not avalailable can be set to NULL. \strong{Default = NULL}
#' @param FolderPath Passed to main function by the user.
#' @param FileName Passed to main function by the user.
#' @param CoRe \emph{Optional: } Passed to main function by the user. If not avalailable can be set to NULL.\strong{Default = FALSE}
#' @param PrintPlot \emph{Optional: } Passed to main function by the user. If not avalailable can be set to NULL.\strong{Default = TRUE}
#' @param PlotHeight \emph{Optional: } Parameter for ggsave.\strong{Default = NULL}
#' @param PlotWidth \emph{Optional: } Parameter for ggsave. \strong{Default = NULL}
#' @param PlotUnit \emph{Optional: } Parameter for ggsave. \strong{Default = NULL}
#'
#' @keywords Save
#' @noRd
#'

SaveRes<- function(InputList_DF= NULL,
                   InputList_Plot= NULL,
                   SaveAs_Table = NULL,
                   SaveAs_Plot = NULL,
                   FolderPath,
                   FileName,
                   CoRe=FALSE,
                   PrintPlot=TRUE,
                   PlotHeight=NULL,
                   PlotWidth=NULL,
                   PlotUnit=NULL){

  ################ Save Tables:
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
      writexl::write_xlsx(InputList_DF, paste0(FileName,".xlsx", sep = "") , col_names = TRUE)
    }else{
      for(DF in names(InputList_DF)){
        #Make FileName
        if(CoRe==FALSE | is.null(CoRe)==TRUE){
          FileName_Save <- paste0(FolderPath,"/" , FileName, "_", DF ,"_",Sys.Date(), sep = "")
        }else{
          FileName_Save <- paste0(FolderPath,"/CoRe_" , FileName, "_", DF ,"_",Sys.Date(), sep = "")
        }

        #unlist DF columns if needed
        InputList_DF[[DF]] <- InputList_DF[[DF]]%>%
          mutate(
            across(
              where(is.list),
              ~map_chr(.x, ~ paste(sort(unique(.x)), collapse = "; "))
            )
          )
        #Save table
        if (SaveAs_Table == "csv"){
          InputList_DF[[DF]]%>%
            readr::write_csv(paste0(FileName_Save,".csv", sep = ""))
        }else if (SaveAs_Table == "txt"){
          InputList_DF[[DF]]%>%
            readr::write_delim(paste0(FileName_Save,".csv", sep = ""))
        }
      }
    }
  }

  ################ Save Plots:
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

      ggplot2::ggsave(filename = paste0(FileName_Save, ".",SaveAs_Plot, sep=""), plot = InputList_Plot[[Plot]], width = PlotWidth,  height = PlotHeight, unit=PlotUnit)

      if(PrintPlot==TRUE){
        suppressMessages(suppressWarnings(plot(InputList_Plot[[Plot]])))
      }
    }
  }
}
