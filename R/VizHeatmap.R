## ---------------------------
##
## Script name: Visualization
##
## Purpose of script: Data Visualisation of the MetaProViz analysis to aid biological interpretation
##
## Author: Dimitrios Prymidis and Christina Schmidt
##
## Date Created: 2022-10-28
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
#'
#' This script allows you to perform different data visualizations using the results of the MetaProViz analysis

###############################
### ### ### Heatmap ### ### ###
###############################

#' @param Input_data DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Includes experimental design and outlier column.
#' @param Plot_SettingsInfo  \emph{Optional: } NULL or Named vector  where you can include vectors or lists for annotation c(individual_Metab= "ColumnName_Plot_SettingsFile_Metab",individual_Sample= "ColumnName_Plot_SettingsFile_Sample", color_Metab="ColumnName_Plot_SettingsFile_Metab", color_Sample= list("ColumnName_Plot_SettingsFile_Sample", "ColumnName_Plot_SettingsFile_Sample",...)).\strong{Default = NULL}
#' @param Plot_SettingsFile_Sample DF which contains information about the samples, which will be combined with your input data based on the unique sample identifiers. and other columns with required PlotSettingInfo.\strong{Default = NULL}
#' @param Plot_SettingsFile_Metab  \emph{Optional: } DF with column "Metabolite" including the Metabolite names (needs to match Metabolite names of Input_data) and other columns with required PlotSettingInfo. \strong{Default = NULL}
#' @param OutputPlotName \emph{Optional: } String which is added to the output files of the plot
#' @param SCALE \emph{Optional: } String with the information for scale row or column. \strong{Default = row}
#' @param Save_as_Plot \emph{Optional: } Select the file type of output plots. Options are svg, pdf or png. \strong{Default = "svg"}
#' @param enforce_FeatureNames \emph{Optional: } If there are more than 100 features no rownames will be shown, which is due to readability. You can enforce this by setting this parameter to TRUE. \strong{Default = FALSE}
#' @param enforce_SampleNames \emph{Optional: } If there are more than 50 sampless no colnames will be shown, which is due to readability. You can enforce this by setting this parameter to TRUE. \strong{Default = FALSE}
#'
#' @keywords Volcano plot, pathways
#' @export
#'
#'


VizHeatmap <- function(Input_data,
                       Plot_SettingsInfo= NULL,
                       Plot_SettingsFile_Sample=NULL,
                       Plot_SettingsFile_Metab= NULL,
                       OutputPlotName= "",
                       SCALE = "row",
                       Save_as_Plot = "svg",
                       enforce_FeatureNames= FALSE,
                       enforce_SampleNames=FALSE
){

  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse", "writexl","pheatmap")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  suppressMessages(library(tidyverse))

  ## ------------ Check Input files ----------- ##
  # 1. The input data:
  if(any(duplicated(row.names(Input_data))) ==  TRUE){# Is the "Input_data" has unique IDs as row names and numeric values in columns?
    stop("Duplicated row.names of Input_data, whilst row.names must be unique")
  } else{
    Test_num <- apply(Input_data, 2, function(x) is.numeric(x))
    if((any(Test_num) ==  FALSE) ==  TRUE){
      stop("Input_data needs to be of class numeric")
    } else(
      Input_data <- Input_data
    )

    # Check if the next lines work correctly in case of duplicated metabolites (colnames)
    if(sum( duplicated(colnames(Input_data))) > 0){
      doublons <- as.character(colnames(Input_data)[duplicated(colnames(Input_data))])#number of duplications
      data <-data[!duplicated(colnames(Input_data)),]#remove duplications
      warning("Input_data contained duplicates based on Metabolite! Dropping duplicate IDs and kept only the first entry. You had ", length(doublons), " duplicates. Note that you should do this before running VizLolipop.")
    }
  }

  #2. Plot_settings
  if(is.null(Plot_SettingsInfo) == FALSE){
    #--- individual_Metab
    if("individual_Metab" %in% names(Plot_SettingsInfo)==TRUE){
      if(Plot_SettingsInfo[["individual_Metab"]] %in% colnames(Plot_SettingsFile_Metab)== FALSE){
        stop("You have chosen individual = ",paste(Plot_SettingsInfo$individual_Metab), ", ", paste(Plot_SettingsInfo$individual_Metab)," does not exist in the Plot_SettingsFile_Metab")
      }
      #Ensure that groups that are assigned NAs do not cause problems:
      Plot_SettingsFile_Metab[[Plot_SettingsInfo[["individual_Metab"]]]] <-ifelse(is.na(Plot_SettingsFile_Metab[[Plot_SettingsInfo[["individual_Metab"]]]]), "NA", Plot_SettingsFile_Metab[[Plot_SettingsInfo[["individual_Metab"]]]])

      unique_paths <- unique(Plot_SettingsFile_Metab[[Plot_SettingsInfo[["individual_Metab"]]]])

      for (i in unique_paths){# Check pathways with 1 metabolite
        selected_path <- Plot_SettingsFile_Metab %>% filter(get(Plot_SettingsInfo[["individual_Metab"]]) == i)
        selected_path_metabs <-  colnames(data) [colnames(data) %in% selected_path$Metabolite]
        if(length(selected_path_metabs)==1 ){
          warning("The metadata group ", paste(i), " includes only 1 metabolite. Heatmap cannot be made for 1 metabolite, thus it will be ignored.")
          unique_paths <- unique_paths[!unique_paths %in% i] # Remove the pathway
        }
      }
    }
    #--- individual_Sample
    if("individual_Sample" %in% names(Plot_SettingsInfo)==TRUE){
        if(Plot_SettingsInfo[["individual_Sample"]] %in% colnames(Plot_SettingsFile_Sample)== FALSE){
          stop("You have chosen individual = ",paste(Plot_SettingsInfo$individual_Sample), ", ", paste(Plot_SettingsInfo$individual_Sample)," does not exist in the Plot_SettingsFile_Metab"   )
        }
        #Ensure that groups that are assigned NAs do not cause problems:
        Plot_SettingsFile_Sample[[Plot_SettingsInfo[["individual_Sample"]]]] <-ifelse(is.na(Plot_SettingsFile_Sample[[Plot_SettingsInfo[["individual_Sample"]]]]), "NA", Plot_SettingsFile_Sample[[Plot_SettingsInfo[["individual_Sample"]]]])

        unique_paths_Sample <- unique(Plot_SettingsFile_Sample[[Plot_SettingsInfo[["individual_Sample"]]]])

        for (i in unique_paths_Sample){# Check pathways with 1 metabolite
          selected_path <- Plot_SettingsFile_Sample %>% filter(get(Plot_SettingsInfo[["individual_Sample"]]) == i)
          selected_path_metabs <-  colnames(data) [colnames(data) %in% selected_path$Metabolite]
          if(length(selected_path_metabs)==1 ){
            warning("The metadata group ", paste(i), " includes only 1 metabolite. Heatmap cannot be made for 1 metabolite, thus it will be ignored.")
            unique_paths_Sample <- unique_paths_Sample[!unique_paths_Sample %in% i] # Remove the pathway
          }
        }
    }
    #--- color_Metab
    if(sum(grepl("color_Metab", names(Plot_SettingsInfo)))>0){ # If color metab exists
      if(is.null(Plot_SettingsFile_Metab)==TRUE){ # Check if Plot_SettingsFile_Metab also exists
        stop("You have chosen Plot_SettingsInfo option that requires you to provide a DF Plot_SettingsFile_Metab.")
      }else{
        if("Metabolite" %in% colnames(Plot_SettingsFile_Metab) == FALSE){  # Check if Metabolite column exists in  Plot_SettingsFile_Metab
          stop("Check input. Plot_SettingsFile_Metab must contain a column named `Metabolite` including the metabolite names.")
        }
        if(sum(colnames(Input_data) %in% Plot_SettingsFile_Metab$Metabolite) != length(Input_data)  ){
          warning("The Input data contains metabolites not found in Plot_SettingsFile_Metab.")
        }
        # Check if color_metab exists in file
        for (metab_color in Plot_SettingsInfo[ grepl("color_Metab", names(Plot_SettingsInfo))]){
          if(metab_color %in% colnames(Plot_SettingsFile_Metab)== FALSE){
            stop("You have chosen color_Metab = ",paste(metab_color), ", ", paste(metab_color)," does not exist in the Plot_SettingsFile_Metab"   )
          }
        }
      }
    }
    #--- color_Sample
    if(sum(grepl("color_Sample", names(Plot_SettingsInfo)))>0){
      if(is.null(Plot_SettingsFile_Sample)==TRUE){ # Check if Plot_SettingsFile_Metab also exists
        stop("You have chosen Plot_SettingsInfo option that requires you to provide a DF Plot_SettingsFile_Sample.")
      }else{
        if(sum(rownames(Input_data) %in% rownames(Plot_SettingsFile_Sample)) != length(rownames(Input_data))){  # Check if all samples exists in the Info_Sample
          stop("Check input. Plot_SettingsFile_Sample must contain the same rownames representing the samples as the Input data.")
        }
        for (samp_color in Plot_SettingsInfo[ grepl("color_Sample", names(Plot_SettingsInfo))]){
          if(samp_color %in% colnames(Plot_SettingsFile_Sample)== FALSE){
            stop("You have chosen color_Sample = ",paste(samp_color), ", ", paste(samp_color)," does not exist in the Plot_SettingsFile_Sample"   )
          }
        }
      }
    }
  }

  # 4. Check other plot-specific parameters:
  Save_as_Plot_options <- c("svg","pdf","png")
  if(Save_as_Plot %in% Save_as_Plot_options == FALSE){
    stop("Check input. The selected Save_as_Plot option is not valid. Please select one of the following: ",paste(Save_as_Plot_options,collapse = ", "),"." )
  }

  ## ------------ Create Output folders ----------- ##
  name <- paste0("MetaProViz_Results_",Sys.Date())
  WorkD <- getwd()
  Results_folder <- file.path(WorkD, name)
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)} # Make Results folder
  Results_folder_plots_Heatmaps_folder = file.path(Results_folder, "Heatmap")
  if (!dir.exists(Results_folder_plots_Heatmaps_folder)) {dir.create(Results_folder_plots_Heatmaps_folder)}  # check and create folder

  #####################################################
  ## -------------- Plot --------------- ##
  data <- Input_data

  if("individual_Metab" %in% names(Plot_SettingsInfo)==TRUE & "individual_Sample" %in% names(Plot_SettingsInfo)==FALSE){
    IndividualPlots <-unique_paths
    PlotList <- list()#Empty list to store all the plots

    for (i in IndividualPlots){
      selected_path <- Plot_SettingsFile_Metab %>% filter(get(Plot_SettingsInfo[["individual_Metab"]]) == i)
      selected_path_metabs <-  colnames(data) [colnames(data) %in% selected_path$Metabolite]
      data_path <- data %>% select(all_of(selected_path_metabs))

      # Column annotation
      col_annot_vars <- Plot_SettingsInfo[grepl("color_Sample", names(Plot_SettingsInfo))]
      col_annot<- NULL
      if(length(col_annot_vars)>0){
        for (x in 1:length(col_annot_vars)){
          annot_sel <- col_annot_vars[[x]]
          col_annot[x] <- Plot_SettingsFile_Sample %>% select(annot_sel) %>% as.data.frame()
          names(col_annot)[x] <- annot_sel
        }
        col_annot<- as.data.frame(col_annot)
        rownames(col_annot) <- rownames(data_path)
      }

      # Row annotation
      row_annot_vars <- Plot_SettingsInfo[grepl("color_Metab", names(Plot_SettingsInfo))]
      row_annot<- NULL
      if(length(row_annot_vars)>0){
        for (y in 1:length(row_annot_vars)){
          annot_sel <- row_annot_vars[[y]]
          row_annot[y] <- Plot_SettingsFile_Metab %>% select(all_of(annot_sel))
          row_annot <- row_annot %>% as.data.frame()
          names(row_annot)[y] <- annot_sel
        }
        rownames(row_annot) <- Plot_SettingsFile_Metab[["Metabolite"]]
      }

      #Check number of features:
      Features <- as.data.frame(t(data_path))
      if(enforce_FeatureNames==TRUE){
        show_rownames <- TRUE
        cellheight_Feature <- 9
      }else if(nrow(Features)>100){
        show_rownames <- FALSE
        cellheight_Feature <- 1
      }else{
        show_rownames <- TRUE
        cellheight_Feature <- 9
      }

      #Check number of samples
      if(enforce_SampleNames==TRUE){
        show_colnames <- TRUE
        cellwidth_Sample <- 9
      }else if(nrow(data_path)>50){
        show_colnames <- FALSE
        cellwidth_Sample <- 1
      }else{
        show_colnames <- TRUE
        cellwidth_Sample <- 9
      }

      # Make the plot
      set.seed(1234)
      heatmap <- pheatmap::pheatmap(t(data_path),
                                     show_rownames = as.logical(show_rownames),
                                     show_colnames = as.logical(show_colnames),
                                     clustering_method =  "complete",
                                     scale = SCALE,
                                     clustering_distance_rows = "correlation",
                                     annotation_col = col_annot,
                                     annotation_row = row_annot,
                                     legend = T,
                                     cellwidth = cellwidth_Sample,
                                     cellheight = cellheight_Feature,
                                     fontsize_row= 10,
                                     fontsize_col = 10,
                                     fontsize=9,
                                     main = paste(OutputPlotName, " Metabolites: ", i, sep=" " ))

      #Width and height according to Sample and metabolite number
      Plot_Sized <- plotGrob_Heatmap(Input=heatmap, Plot_SettingsInfo=Plot_SettingsInfo, data_path=data_path, show_rownames=show_rownames, show_colnames=show_colnames,  Plot_SettingsFile_Sample= Plot_SettingsFile_Sample,  Plot_SettingsFile_Metab= Plot_SettingsFile_Metab)
      heatmap <-Plot_Sized[[3]]

      #----- Save
      cleaned_i <- gsub("[[:space:],/\\\\]", "-", i)#removes empty spaces and replaces /,\ with -
      ggsave(file=paste(Results_folder_plots_Heatmaps_folder,"/", "Heatmap_",cleaned_i,"_",OutputPlotName, ".",Save_as_Plot, sep=""), plot=heatmap, width=Plot_Sized[[2]], height=Plot_Sized[[1]], unit="cm")
      #plot(heatmap)
      }
    }else if("individual_Metab" %in% names(Plot_SettingsInfo)==FALSE & "individual_Sample" %in% names(Plot_SettingsInfo)==TRUE){
      IndividualPlots <-unique_paths_Sample
      PlotList <- list()#Empty list to store all the plots

      for (i in IndividualPlots){
        #Select the data:
        selected_path <- Plot_SettingsFile_Sample %>% filter(get(Plot_SettingsInfo[["individual_Sample"]]) == i)%>%
          rownames_to_column("UniqueID")
        selected_path <- as.data.frame(selected_path[,1])%>%
          dplyr::rename("UniqueID"=1)
        data_path <- merge(selected_path, data%>%rownames_to_column("UniqueID"), by="UniqueID", all.x=TRUE)
        data_path <- data_path%>%
          column_to_rownames("UniqueID")

        # Column annotation
        selected_Plot_SettingsFile_Sample <- merge(selected_path, Plot_SettingsFile_Sample%>%rownames_to_column("UniqueID"), by="UniqueID", all.x=TRUE)

        col_annot_vars <- Plot_SettingsInfo[grepl("color_Sample", names(Plot_SettingsInfo))]
        col_annot<- NULL
        if(length(col_annot_vars)>0){
          for (x in 1:length(col_annot_vars)){
            annot_sel <- col_annot_vars[[x]]
            col_annot[x] <- selected_Plot_SettingsFile_Sample %>% select(annot_sel) %>% as.data.frame()
            names(col_annot)[x] <- annot_sel
          }
          col_annot<- as.data.frame(col_annot)
          rownames(col_annot) <- rownames(data_path)
        }

        # Row annotation
        row_annot_vars <- Plot_SettingsInfo[grepl("color_Metab", names(Plot_SettingsInfo))]
        row_annot<- NULL
        if(length(row_annot_vars)>0){
          for (y in 1:length(row_annot_vars)){
            annot_sel <- row_annot_vars[[y]]
            row_annot[y] <- Plot_SettingsFile_Metab %>% select(all_of(annot_sel))
            row_annot <- row_annot %>% as.data.frame()
            names(row_annot)[y] <- annot_sel
          }
          rownames(row_annot) <- Plot_SettingsFile_Metab[["Metabolite"]]
        }

        #Check number of features:
        Features <- as.data.frame(t(data_path))
        if(enforce_FeatureNames==TRUE){
          show_rownames <- TRUE
          cellheight_Feature <- 9
        }else if(nrow(Features)>100){
          show_rownames <- FALSE
          cellheight_Feature <- 1
        }else{
          show_rownames <- TRUE
          cellheight_Feature <- 9
        }

        #Check number of samples
        if(enforce_SampleNames==TRUE){
          show_colnames <- TRUE
          cellwidth_Sample <- 9
        }else if(nrow(data_path)>50){
          show_colnames <- FALSE
          cellwidth_Sample <- 1
        }else{
          show_colnames <- TRUE
          cellwidth_Sample <- 9
        }

        # Make the plot
        set.seed(1234)
        heatmap <- pheatmap::pheatmap(t(data_path),
                                      show_rownames = as.logical(show_rownames),
                                      show_colnames = as.logical(show_colnames),
                                      clustering_method =  "complete",
                                      scale = SCALE,
                                      clustering_distance_rows = "correlation",
                                      annotation_col = col_annot,
                                      annotation_row = row_annot,
                                      legend = T,
                                      cellwidth = cellwidth_Sample,
                                      cellheight = cellheight_Feature,
                                      fontsize_row= 10,
                                      fontsize_col = 10,
                                      fontsize=9,
                                      main = paste(OutputPlotName," Samples: ", i, sep=" " ))

        #Width and height according to Sample and metabolite number
        Plot_Sized <- plotGrob_Heatmap(Input=heatmap, Plot_SettingsInfo=Plot_SettingsInfo, data_path=data_path, show_rownames=show_rownames, show_colnames=show_colnames,  Plot_SettingsFile_Sample= Plot_SettingsFile_Sample,  Plot_SettingsFile_Metab= Plot_SettingsFile_Metab)
        heatmap <-Plot_Sized[[3]]

        #----- Save
        cleaned_i <- gsub("[[:space:],/\\\\]", "-", i)#removes empty spaces and replaces /,\ with -
        ggsave(file=paste(Results_folder_plots_Heatmaps_folder,"/", "Heatmap_",cleaned_i,"_",OutputPlotName, ".",Save_as_Plot, sep=""), plot=heatmap, width=Plot_Sized[[2]], height=Plot_Sized[[1]], unit="cm")
        #plot(heatmap)
        }
      }else if("individual_Metab" %in% names(Plot_SettingsInfo)==TRUE & "individual_Sample" %in% names(Plot_SettingsInfo)==TRUE){
        IndividualPlots_Metab <-unique_paths
        IndividualPlots_Sample <-unique_paths_Sample

        PlotList <- list()#Empty list to store all the plots
        for (i in IndividualPlots_Metab){
          selected_path <- Plot_SettingsFile_Metab %>% filter(get(Plot_SettingsInfo[["individual_Metab"]]) == i)
          selected_path_metabs <-  colnames(data) [colnames(data) %in% selected_path$Metabolite]
          data_path_metab <- data %>% select(all_of(selected_path_metabs))

          # Row annotation
          row_annot_vars <- Plot_SettingsInfo[grepl("color_Metab", names(Plot_SettingsInfo))]
          row_annot<- NULL
          if(length(row_annot_vars)>0){
            for (y in 1:length(row_annot_vars)){
              annot_sel <- row_annot_vars[[y]]
              row_annot[y] <- Plot_SettingsFile_Metab %>% select(all_of(annot_sel))
              row_annot <- row_annot %>% as.data.frame()
              names(row_annot)[y] <- annot_sel
            }
            rownames(row_annot) <- Plot_SettingsFile_Metab[["Metabolite"]]
          }

          #Col annotation:
          for (s in IndividualPlots_Sample){
            #Select the data:
            selected_path <- Plot_SettingsFile_Sample %>% filter(get(Plot_SettingsInfo[["individual_Sample"]]) == s)%>%
              rownames_to_column("UniqueID")
            selected_path <- as.data.frame(selected_path[,1])%>%
              dplyr::rename("UniqueID"=1)
            data_path <- merge(selected_path, data_path_metab%>%rownames_to_column("UniqueID"), by="UniqueID", all.x=TRUE)
            data_path <- data_path%>%
              column_to_rownames("UniqueID")

            # Column annotation
            selected_Plot_SettingsFile_Sample <- merge(selected_path, Plot_SettingsFile_Sample%>%rownames_to_column("UniqueID"), by="UniqueID", all.x=TRUE)

            col_annot_vars <- Plot_SettingsInfo[grepl("color_Sample", names(Plot_SettingsInfo))]
            col_annot<- NULL
            if(length(col_annot_vars)>0){
              for (x in 1:length(col_annot_vars)){
                annot_sel <- col_annot_vars[[x]]
                col_annot[x] <- selected_Plot_SettingsFile_Sample %>% select(annot_sel) %>% as.data.frame()
                names(col_annot)[x] <- annot_sel
              }
              col_annot<- as.data.frame(col_annot)
              rownames(col_annot) <- rownames(data_path)
            }

            #Check number of features:
            Features <- as.data.frame(t(data_path))
            if(enforce_FeatureNames==TRUE){
              show_rownames <- TRUE
              cellheight_Feature <- 9
            }else if(nrow(Features)>100){
              show_rownames <- FALSE
              cellheight_Feature <- 1
            }else{
              show_rownames <- TRUE
              cellheight_Feature <- 9
            }

            #Check number of samples
            if(enforce_SampleNames==TRUE){
              show_colnames <- TRUE
              cellwidth_Sample <- 9
            }else if(nrow(data_path)>50){
              show_colnames <- FALSE
              cellwidth_Sample <- 1
            }else{
              show_colnames <- TRUE
              cellwidth_Sample <- 9
            }

            # Make the plot
            set.seed(1234)
            heatmap <- pheatmap::pheatmap(t(data_path),
                                          show_rownames = as.logical(show_rownames),
                                          show_colnames = as.logical(show_colnames),
                                          clustering_method =  "complete",
                                          scale = SCALE,
                                          clustering_distance_rows = "correlation",
                                          annotation_col = col_annot,
                                          annotation_row = row_annot,
                                          legend = T,
                                          cellwidth = cellwidth_Sample,
                                          cellheight = cellheight_Feature,
                                          fontsize_row= 10,
                                          fontsize_col = 10,
                                          fontsize=9,
                                          main = paste(OutputPlotName," Metabolites: ", i, " Sample:", s, sep="" ))

            #-------- Plot width and heights
            #Width and height according to Sample and metabolite number
            Plot_Sized <- plotGrob_Heatmap(Input=heatmap, Plot_SettingsInfo=Plot_SettingsInfo, data_path=data_path, show_rownames=show_rownames, show_colnames=show_colnames,  Plot_SettingsFile_Sample= Plot_SettingsFile_Sample,  Plot_SettingsFile_Metab= Plot_SettingsFile_Metab)
            heatmap <-Plot_Sized[[3]]

            #----- Save
            cleaned_i <- gsub("[[:space:],/\\\\]", "-", i)#removes empty spaces and replaces /,\ with -
            cleaned_s <- gsub("[[:space:],/\\\\]", "-", s)#removes empty spaces and replaces /,\ with -
            ggsave(file=paste(Results_folder_plots_Heatmaps_folder,"/", "Heatmap_",cleaned_i,"_",cleaned_s, "_",OutputPlotName, ".",Save_as_Plot, sep=""), plot=heatmap, width=Plot_Sized[[2]], height=Plot_Sized[[1]], units = "cm")

            #plot(heatmap)
            }
        }
    } else if("individual_Metab" %in% names(Plot_SettingsInfo)==FALSE & "individual_Sample" %in% names(Plot_SettingsInfo)==FALSE){

    # Column annotation
    col_annot_vars <- Plot_SettingsInfo[grepl("color_Sample", names(Plot_SettingsInfo))]
    col_annot<- NULL
    if(length(col_annot_vars)>0){
      for (i in 1:length(col_annot_vars)){
        annot_sel <- col_annot_vars[[i]]
        col_annot[i] <- Plot_SettingsFile_Sample %>% select(annot_sel) %>% as.data.frame()
        names(col_annot)[i] <- annot_sel
      }
      col_annot<- as.data.frame(col_annot)
      rownames(col_annot) <- rownames(data)
    }

    # Row annotation
    row_annot_vars <- Plot_SettingsInfo[grepl("color_Metab", names(Plot_SettingsInfo))]
    row_annot<- NULL
    if(length(row_annot_vars)>0){
      for (i in 1:length(row_annot_vars)){
        annot_sel <- row_annot_vars[[i]]
        row_annot[i] <- Plot_SettingsFile_Metab %>% select(all_of(annot_sel))
        row_annot <- row_annot %>% as.data.frame()
        names(row_annot)[i] <- annot_sel
      }
      rownames(row_annot) <- Plot_SettingsFile_Metab[["Metabolite"]]
    }

    #Check number of features:
    Features <- as.data.frame(t(data))
    if(enforce_FeatureNames==TRUE){
      show_rownames <- TRUE
      cellheight_Feature <- 9
    }else if(nrow(Features)>100){
      show_rownames <- FALSE
      cellheight_Feature <- 1
    }else{
      show_rownames <- TRUE
      cellheight_Feature <- 9
    }

    #Check number of samples
    if(enforce_SampleNames==TRUE){
      show_colnames <- TRUE
      cellwidth_Sample <- 9
    }else if(nrow(data)>50){
      show_colnames <- FALSE
      cellwidth_Sample <- 1
    }else{
      show_colnames <- TRUE
      cellwidth_Sample <- 9
    }

    #Make the plot:
    set.seed(1234)
    heatmap <- pheatmap::pheatmap(t(data),
                                  show_rownames = as.logical(show_rownames),
                                  show_colnames = as.logical(show_colnames),
                                  clustering_method =  "complete",
                                  scale = SCALE,
                                  clustering_distance_rows = "correlation",
                                  annotation_col = col_annot,
                                  annotation_row = row_annot,
                                  legend = T,
                                  cellwidth = cellwidth_Sample,
                                  cellheight = cellheight_Feature,
                                  fontsize_row= 10,
                                  fontsize_col = 10,
                                  fontsize=9,
                                  main = OutputPlotName)

    #-------- Plot width and heights
    #Width and height according to Sample and metabolite number
    Plot_Sized <- plotGrob_Heatmap(Input=heatmap, Plot_SettingsInfo=Plot_SettingsInfo, data_path=data, show_rownames=show_rownames, show_colnames=show_colnames,  Plot_SettingsFile_Sample= Plot_SettingsFile_Sample,  Plot_SettingsFile_Metab= Plot_SettingsFile_Metab)
    heatmap <-Plot_Sized[[3]]

    #----- Save
    ggsave(file=paste(Results_folder_plots_Heatmaps_folder,"/", "Heatmap",OutputPlotName, ".", Save_as_Plot ,sep=""), plot=heatmap, width=Plot_Sized[[2]], height=Plot_Sized[[1]], units = "cm")

    #plot(heatmap)
  }
}



##############################################################
### ### ### Heatmap helper function: Internal Function ### ### ###
##############################################################

#' @param Input This is the ggplot object generated within the VizHeatmap function.
#' @param Plot_SettingsInfo Passed to VizHeatmap
#' @param data Generated in VizVolcano from Input_data and used to plot heatmap
#' @param show_rownames Generated in VizVolcano, is logical TRUE or FALSE
#' @param show_colnames Generated in VizVolcano, is logical TRUE or FALSE
#' @param Plot_SettingsFile_Sample Passed to VizHeatmap
#' @param Plot_SettingsFile_Metab Passed to VizHeatmap
#'
#' @keywords Heatmap helper function
#' @noRd

plotGrob_Heatmap <- function(Input=heatmap, Plot_SettingsInfo, data_path, show_rownames, show_colnames, Plot_SettingsFile_Sample, Plot_SettingsFile_Metab){
  #------- Set the total heights and widths
  #we need ggplot_grob to edit the gtable of the ggplot object. Using this we can manipulate the gtable arguments directly.
  plottable<- Input$gtable#get a gtable --> gtable::gtable_show_layout(plottable)

  #----- heights
  plottable$heights[1] <- unit(1, "cm")#space for headline
  plottable$heights[2] <- unit((grid::convertX(unit(as.numeric(plottable$heights[2]), "bigpts"), "cm", valueOnly = TRUE)), "cm")#heights x-axis clustering to heatmap, do not change
  plottable$heights[3] <- unit((grid::convertX(unit(as.numeric(plottable$heights[3]), "bigpts"), "cm", valueOnly = TRUE)), "cm")#distance x-axis clustering to heatmap, do not change
  plottable$heights[4] <- unit((grid::convertX(unit(as.numeric(plottable$heights[4]), "bigpts"), "cm", valueOnly = TRUE)), "cm")#heatmap itself, do not change

  if(sum(grepl("color_Metab", names(Plot_SettingsInfo)))>0){#We have a legend for the color on the y-axis
    # Find the longest column name/color name and count its characters:
    names <- Plot_SettingsInfo[grepl("color_Metab", names(Plot_SettingsInfo))]
    colour_names <- NULL
    for (x in 1:length(names)){
      names_sel <- names[[x]]
      colour_names[x] <- names_sel
    }

    if(show_colnames==TRUE){#check if we also display sample names
      column_names <- c(rownames(data_path),colour_names)
    }else{
      column_names <- colour_names
    }

    longest_column_name <- column_names[which.max(nchar(column_names))]
    character_count <- nchar(longest_column_name)

    plottable$heights[5] <- unit(character_count*0.3, "cm")#sample labels on heatmap
    plot_heights <- 1+(character_count*0.3) + (as.numeric(gsub("^(\\d+\\.?\\d*).*", "\\1", plottable$heights[2])))+ (as.numeric(gsub("^(\\d+\\.?\\d*).*", "\\1", plottable$heights[3])))+ (as.numeric(gsub("^(\\d+\\.?\\d*).*", "\\1", plottable$heights[3])))+ (as.numeric(gsub("^(\\d+\\.?\\d*).*", "\\1", plottable$heights[4])))
  }else if(sum(grepl("color_Metab", names(Plot_SettingsInfo)))==0){
    # Find the longest column name and count its characters:
    if(show_colnames==TRUE){#check if we also display sample names
      column_names <- rownames(data_path)
      longest_column_name <- column_names[which.max(nchar(column_names))]
      character_count <- nchar(longest_column_name)
    }else{
      character_count <- 2
    }

    plottable$heights[5] <- unit(character_count*0.3, "cm")#sample labels on heatmap

    #Summary:
    plot_heights <- 1+(character_count*0.3) + (as.numeric(gsub("^(\\d+\\.?\\d*).*", "\\1", plottable$heights[2])))+ (as.numeric(gsub("^(\\d+\\.?\\d*).*", "\\1", plottable$heights[3])))+ (as.numeric(gsub("^(\\d+\\.?\\d*).*", "\\1", plottable$heights[3])))+ (as.numeric(gsub("^(\\d+\\.?\\d*).*", "\\1", plottable$heights[4])))
    }

  #--- widths:
  plottable$widths[1] <- unit((grid::convertX(unit(as.numeric(plottable$widths[1]), "bigpts"), "cm", valueOnly = TRUE)), "cm")#clustering, do not change
  plottable$widths[2] <- unit((grid::convertX(unit(as.numeric(plottable$widths[2]), "bigpts"), "cm", valueOnly = TRUE)), "cm")##distance clustering to plot, do not change
  plottable$widths[3] <- unit((grid::convertX(unit(as.numeric(plottable$widths[3]), "bigpts"), "cm", valueOnly = TRUE)), "cm")#heatmap itself, do not change
  plottable$widths[5] <- unit(2, "cm")#space between heatmap colour legend and legend

  #legend width
  if((sum(grepl("color_Sample", names(Plot_SettingsInfo)))>0) | (sum(grepl("color_Metab", names(Plot_SettingsInfo)))>0)){
    if(sum(grepl("color_Sample", names(Plot_SettingsInfo)))>0){
      names <- Plot_SettingsInfo[grepl("color_Sample", names(Plot_SettingsInfo))]
      colour_names <- NULL
      legend_names <- NULL
      for (x in 1:length(names)){
        names_sel <- names[[x]]
        legend_names[x] <- names_sel
        colour_names[x] <- Plot_SettingsFile_Sample[names[[x]]]
      }
      }else{
        colour_names <- NULL
        legend_names <- NULL
        }

    if(sum(grepl("color_Metab", names(Plot_SettingsInfo)))>0){
      names <- Plot_SettingsInfo[grepl("color_Metab", names(Plot_SettingsInfo))]
      colour_names_M <- NULL
      legend_names_M <- NULL
      for (x in 1:length(names)){
        names_sel <- names[[x]]
        legend_names_M[x] <- names_sel
        colour_names_M[x] <- Plot_SettingsFile_Metab[names[[x]]]
      }
      }else{
        colour_names_M <- NULL
        legend_names_M <- NULL
        }

    legend_head <- c(legend_names, legend_names_M)
    longest_name <- legend_head[which.max(nchar(legend_head[[1]]))]
    character_count_head <- nchar(longest_name)+4#This is the length of the legend title name

    legend_names <- c(unlist(colour_names), unlist(colour_names_M))
    longest_name <- legend_names[which.max(nchar(legend_names[[1]]))]
    character_count <- nchar(longest_name)#This is the length of the legend colour names


    plottable$widths[6] <- unit(((max(character_count_head, character_count))*0.2)+0.7, "cm")#legend space
  }else{
    plottable$widths[6] <- unit(0, "cm")#legend space
  }

  #heatmap labels
  if(sum(grepl("color_Sample", names(Plot_SettingsInfo)))>0 ){#We have a legend for the color on the x-axis
    # Find the longest column name/color name and count its characters:
    names <- Plot_SettingsInfo[grepl("color_Sample", names(Plot_SettingsInfo))]
    colour_names <- NULL
    for (x in 1:length(names)){
        names_sel <- names[[x]]
        colour_names[x] <- names_sel
    }

    if(show_rownames==TRUE){#check if we also display feature names
      column_names <- c(colnames(data_path),colour_names)
      longest_column_name <- column_names[which.max(nchar(column_names))]
      character_count <- nchar(longest_column_name)
    }else{
      character_count <- 2
    }

    plottable$widths[4] <- unit(character_count*0.2, "cm")#feature labels on heatmap --> unit(sum(grid::convertX(unit(as.numeric(1), "grobwidth", data=plottable), "cm", valueOnly = TRUE), grid::convertX(unit(as.numeric(10), "bigpts", data=plottable), "cm", valueOnly = TRUE)), "cm")

    #Summary:
    plot_widths <- 2+(as.numeric(gsub("^(\\d+\\.?\\d*).*", "\\1", plottable$widths[4])))+ (as.numeric(gsub("^(\\d+\\.?\\d*).*", "\\1", plottable$widths[1])))+ (as.numeric(gsub("^(\\d+\\.?\\d*).*", "\\1", plottable$widths[2])))+ (as.numeric(gsub("^(\\d+\\.?\\d*).*", "\\1", plottable$widths[3])))+ (as.numeric(gsub("^(\\d+\\.?\\d*).*", "\\1", plottable$widths[6])))

  }else if(sum(grepl("color_Sample", names(Plot_SettingsInfo)))==0){
    # Find the longest column name and count its characters:
    if(show_rownames==TRUE){#check if we also display feature names
      column_names <- colnames(data_path)
      longest_column_name <- column_names[which.max(nchar(column_names))]
      character_count <- nchar(longest_column_name)
    }else{
      character_count <- 2
    }

    plottable$widths[4] <- unit(character_count*0.2, "cm")#feature labels on heatmap --> unit(sum(grid::convertX(unit(as.numeric(1), "grobwidth", data=plottable), "cm", valueOnly = TRUE), grid::convertX(unit(as.numeric(10), "bigpts", data=plottable), "cm", valueOnly = TRUE)), "cm")

    #Summary:
    plot_widths <- 2+(as.numeric(gsub("^(\\d+\\.?\\d*).*", "\\1", plottable$widths[4])))+ (as.numeric(gsub("^(\\d+\\.?\\d*).*", "\\1", plottable$widths[1])))+ (as.numeric(gsub("^(\\d+\\.?\\d*).*", "\\1", plottable$widths[2])))+ (as.numeric(gsub("^(\\d+\\.?\\d*).*", "\\1", plottable$widths[3])))+ (as.numeric(gsub("^(\\d+\\.?\\d*).*", "\\1", plottable$widths[6])))
    }

  #plot_param <-c(plot_heights=plot_heights, plot_widths=plot_widths)
  Output<- list(plot_heights, plot_widths, plottable)
}



