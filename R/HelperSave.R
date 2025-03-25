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
SavePath<- function(FolderName, FolderPath) {
  
    ## check if FolderName includes special characters that are not allowed
    cleaned_FolderName <- gsub("[^a-zA-Z0-9 ]", "", FolderName)
    if (FolderName != cleaned_FolderName){
        message("Special characters were removed from `FolderName`.")
    }

    ## check if FolderPath exist
    if (is.null(FolderPath)) {
        FolderPath <- getwd()
        FolderPath <- file.path(FolderPath, "MetaProViz_Results")
        if (!dir.exists(FolderPath)) {
            dir.create(FolderPath)
        }
    } else {
        if (!dir.exists(FolderPath)) {
            FolderPath <- getwd()
            message("Provided `FolderPath` does not exist and hence results are saved here: ", 
                FolderPath, sep = "")
        }
    }

    ## create the folder name
    Results_folder <- file.path(FolderPath, cleaned_FolderName)
    if (!dir.exists(Results_folder)) {
        dir.create(Results_folder)
    }

    ## return the folder path
    invisible(Results_folder)
}


#' Make sure the results directory exists
#'
#' @importFrom magrittr %>%
#' @noRd
ResultsDir <- function(path = 'MetaProViz_Results') {

  # TODO: options?
  path %>% {`if`(!dir.exists(.), {dir.create(.); .}, .) } ## EDIT: is this easily understandable? Is the classical, non-"tidy" notation easier to read?

}


################################################################################################
### ### ### General helper function: Save tables and plot ### ### ###
################################################################################################

#' SaveRes is the helper function to save the plots and tables
#'
#' @param InputList_DF \emph{Optional: } Generated within the MetaProViz function. Contains named DFs. If not availalable can be set to NULL.\strong{Default = NULL}
#' @param InputList_Plot \emph{Optional: } Generated within the MetaProViz function. Contains named Plots. If not availalable can be set to NULL.\strong{Default = NULL}
#' @param SaveAs_Table \emph{Optional: } Passed to main function by the user. If not availalable can be set to NULL.\strong{Default = NULL}
#' @param SaveAs_Plot \emph{Optional: } Passed to main function by the user. If not availalable can be set to NULL. \strong{Default = NULL}
#' @param FolderPath Passed to main function by the user.
#' @param FileName Passed to main function by the user.
#' @param CoRe \emph{Optional: } Passed to main function by the user. If not availalable can be set to NULL.\strong{Default = FALSE}
#' @param PrintPlot \emph{Optional: } Passed to main function by the user. If not availalable can be set to NULL.\strong{Default = TRUE}
#' @param PlotHeight \emph{Optional: } Parameter for ggsave.\strong{Default = NULL}
#' @param PlotWidth \emph{Optional: } Parameter for ggsave. \strong{Default = NULL}
#' @param PlotUnit \emph{Optional: } Parameter for ggsave. \strong{Default = NULL}
#'
#' @keywords Save
#' @noRd
#'

SaveRes <- function(data = NULL,
    plot = NULL,
    SaveAs_Table = NULL, ## EDIT: name the options here and use match.arg ## EDIT: should have a different name, e.g. saveAsFormat ? ## EDIT: camel case and snake case notation should not be mixed
    SaveAs_Plot = NULL,
    FolderPath,
    FileName,
    CoRe = FALSE,
    PrintPlot = TRUE,
    PlotHeight = NULL,
    PlotWidth = NULL,
    PlotUnit = NULL) {

    ################ Save Tables:
    if (!is.null(SaveAs_Table)) {
        ## Excel File: One file with multiple sheets:
        
        for (data_i in names(data)) { ## EDIT: not every entry in list is a data.frame, check and adjust accordingly
                ## make FileName
                
            ##FileName <- paste0(FileName, "_", data_i, "_", Sys.Date(), sep = "") EDIT: is this better?
            if (!CoRe | is.null(CoRe)) {
                FileName_Save <- paste0(FolderPath,"/" , FileName, "_", data_i,"_", Sys.Date(), sep = "")
                ##FileName <- paste0(FolderPath, "/", FileName)
            } else {
                FileName_Save <- paste0(FolderPath,"/CoRe_" , FileName, "_", data_i,"_", Sys.Date(), sep = "")
                ##FileName <- paste0(FolderPath, "/CoRe_", FileName)
            }
                
            if (!is(data[[data_i]], "SummarizedExperiment")) {
                
                ## unlist data_i columns if needed
                data[[data_i]] <- data[[data_i]] %>%
                    mutate(
                        across(
                            where(is.list),
                            ~map_chr(.x, ~ paste(sort(unique(.x)), collapse = "; "))
                        )
                    )
                ## Save table
                ## save Excel
                if (SaveAs_Table == "xlsx") {
                    writexl::write_xlsx(data[[data_i]], paste0(FileName_Save, ".xlsx"),
                                        col_names = TRUE)    
                }
                ## save csv
                if (SaveAs_Table == "csv") {
                    data[[data_i]] %>%
                        readr::write_csv(paste0(FileName_Save, ".csv"))
                } 
                if (SaveAs_Table == "txt") {
                    data[[data_i]] %>%
                        readr::write_delim(paste0(FileName_Save, ".csv"))
                }
            } else {
                saveRDS(data[[data_i]], file = paste0(FileName_Save, ".RDS"))
            }
        }
    }

    ################ Save Plots:
    if(!is.null(SaveAs_Plot)){
        for (plot_i in names(plot)) {
         ## make FileName
            ##FileName <- paste0(FileName, "_", Sys.Date(), sep = "") EDIT: is this better?
        
            if (!CoRe | is.null(CoRe)) {
                FileName_Save <- paste0(FolderPath,"/" , FileName,"_", plot_i, "_",Sys.Date(), sep = "")
                ##FileName <- paste0(FolderPath, "/", FileName)
            } else {
                FileName_Save <- paste0(FolderPath,"/CoRe_" , FileName,"_", plot_i,"_",Sys.Date(), sep = "")
                ##FileName <- paste0(FolderPath, "/CoRe_", FileName)
            }

            ## save
            if (is.null(PlotHeight)) {
                PlotHeight <- 12
            }
            if (is.null(PlotWidth)) {
                PlotWidth <- 16
            }
            if (is.null(PlotUnit)) {
                PlotUnit <- "cm"
            }

            ggplot2::ggsave(filename = paste0(FileName_Save, ".", SaveAs_Plot), 
                plot = plot[[plot_i]], width = PlotWidth, 
                height = PlotHeight, unit = PlotUnit)

            if (PrintPlot) {
                suppressMessages(suppressWarnings(
                    plot(plot[[plot_i]])))
            }
        }
    }
}
