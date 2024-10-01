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

#' @param InputData DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Includes experimental design and outlier column.
#' @param SettingsInfo  \emph{Optional: } NULL or Named vector  where you can include vectors or lists for annotation c(individual_Metab= "ColumnName_SettingsFile_Metab",individual_Sample= "ColumnName_SettingsFile_Sample", color_Metab="ColumnName_SettingsFile_Metab", color_Sample= list("ColumnName_SettingsFile_Sample", "ColumnName_SettingsFile_Sample",...)).\strong{Default = NULL}
#' @param SettingsFile_Sample DF which contains information about the samples, which will be combined with your input data based on the unique sample identifiers. and other columns with required PlotSettingInfo.\strong{Default = NULL}
#' @param SettingsFile_Metab  \emph{Optional: } DF with column "Metabolite" including the Metabolite names (needs to match Metabolite names of Input_data) and other columns with required PlotSettingInfo. \strong{Default = NULL}
#' @param PlotName \emph{Optional: } String which is added to the output files of the plot
#' @param Scale \emph{Optional: } String with the information for Scale row, column or none. \strong{Default = row}
#' @param SaveAs_Plot \emph{Optional: } Select the file type of output plots. Options are svg, pdf, png or NULL. \strong{Default = "svg"}
#' @param Enforce_FeatureNames \emph{Optional: } If there are more than 100 features no rownames will be shown, which is due to readability. You can Enforce this by setting this parameter to TRUE. \strong{Default = FALSE}
#' @param Enforce_SampleNames \emph{Optional: } If there are more than 50 sampless no colnames will be shown, which is due to readability. You can Enforce this by setting this parameter to TRUE. \strong{Default = FALSE}
#' @param Folder_Name {Optional:} String which is added to the resulting folder name \strong{default: NULL}
#'
#' @keywords Volcano plot, pathways
#' @export
#'
#'


VizHeatmap <- function(InputData,
                       SettingsInfo= NULL,
                       SettingsFile_Sample=NULL,
                       SettingsFile_Metab= NULL,
                       PlotName= "",
                       Scale = "row",
                       SaveAs_Plot = "svg",
                       Enforce_FeatureNames= FALSE,
                       Enforce_SampleNames= FALSE,
                       PrintPlot=TRUE,
                       FolderPath = NULL
){

  ## ------------ Check Input files ----------- ##
  # HelperFunction `CheckInput`
  MetaProViz:::CheckInput(InputData=InputData,
                          SettingsFile_Sample=SettingsFile_Sample,
                          SettingsFile_Metab=SettingsFile_Metab,
                          SettingsInfo=SettingsInfo,
                          SaveAs_Plot=SaveAs_Plot,
                          SaveAs_Table=NULL,
                          CoRe=FALSE,
                          PrintPlot= PrintPlot)

  # CheckInput` Specific
  if(is.logical(Enforce_FeatureNames) == FALSE | is.logical(Enforce_SampleNames) == FALSE){
    stop("Check input. The Enforce_FeatureNames and Enforce_SampleNames value should be either =TRUE or = FALSE.")
  }

  Scale_options <- c("row","column", "none")
  if(Scale %in% Scale_options == FALSE){
      stop("Check input. The selected Scale option is not valid. Please select one of the folowwing: ",paste(Scale_options,collapse = ", "),"." )
    }


  ## ------------ Create Results output folder ----------- ##
  if(is.null(SaveAs_Plot)==FALSE){
    Folder <- MetaProViz:::SavePath(FolderName= "Heatmap",
                                    FolderPath=FolderPath)
  }

  #####################################################
  ## -------------- Load Data --------------- ##
  data <- InputData

  if(is.null(SettingsFile_Metab)==FALSE){#removes information about metabolites that are not included in the InputData
    SettingsFile_Metab <- merge(x=SettingsFile_Metab, y=as.data.frame(t(InputData)), by=0, all.y=TRUE)%>%
      column_to_rownames("Row.names")
    SettingsFile_Metab <- SettingsFile_Metab[,-c((ncol(SettingsFile_Metab)-nrow(InputData)+1):ncol(SettingsFile_Metab))]
  }

    if(is.null(SettingsFile_Sample)==FALSE){#removes information about samples that are not included in the InputData
      SettingsFile_Sample <- merge(x=SettingsFile_Sample, y=InputData, by=0, all.y=TRUE)%>%
        column_to_rownames("Row.names")
      SettingsFile_Sample <- SettingsFile_Sample[,-c((ncol(SettingsFile_Sample)-ncol(InputData)+1):ncol(SettingsFile_Sample))]
    }


  ## -------------- Plot --------------- ##
  if("individual_Metab" %in% names(SettingsInfo)==TRUE & "individual_Sample" %in% names(SettingsInfo)==FALSE){
    #Ensure that groups that are assigned NAs do not cause problems:
    SettingsFile_Metab[[SettingsInfo[["individual_Metab"]]]] <-ifelse(is.na(SettingsFile_Metab[[SettingsInfo[["individual_Metab"]]]]), "NA", SettingsFile_Metab[[SettingsInfo[["individual_Metab"]]]])
    unique_paths <- unique(SettingsFile_Metab[[SettingsInfo[["individual_Metab"]]]])

    for (i in unique_paths){# Check pathways with 1 metabolite
      selected_path <- SettingsFile_Metab %>% filter(get(SettingsInfo[["individual_Metab"]]) == i)
      selected_path_metabs <-  colnames(data) [colnames(data) %in% row.names(selected_path)]
      if(length(selected_path_metabs)==1 ){
        warning("The metadata group ", paste(i), " includes only 1 metabolite. Heatmap cannot be made for 1 metabolite, thus it will be ignored.")
        unique_paths <- unique_paths[!unique_paths %in% i] # Remove the pathway
      }
    }

    IndividualPlots <-unique_paths

    PlotList <- list()#Empty list to store all the plots
    PlotList_adaptedGrid <- list()#Empty list to store all the plots

    for (i in IndividualPlots){
      selected_path <- SettingsFile_Metab %>% filter(get(SettingsInfo[["individual_Metab"]]) == i)
      selected_path_metabs <-  colnames(data) [colnames(data) %in% row.names(selected_path)]
      data_path <- data %>% select(all_of(selected_path_metabs))

      # Column annotation
      col_annot_vars <- SettingsInfo[grepl("color_Sample", names(SettingsInfo))]
      col_annot<- NULL
      if(length(col_annot_vars)>0){
        for (x in 1:length(col_annot_vars)){
          annot_sel <- col_annot_vars[[x]]
          col_annot[x] <- SettingsFile_Sample %>% select(annot_sel) %>% as.data.frame()
          names(col_annot)[x] <- annot_sel
        }
        col_annot<- as.data.frame(col_annot)
        rownames(col_annot) <- rownames(data_path)
      }

      # Row annotation
      row_annot_vars <- SettingsInfo[grepl("color_Metab", names(SettingsInfo))]
      row_annot<- NULL
      if(length(row_annot_vars)>0){
        for (y in 1:length(row_annot_vars)){
          annot_sel <- row_annot_vars[[y]]
          row_annot[y] <- SettingsFile_Metab %>% select(all_of(annot_sel))
          row_annot <- row_annot %>% as.data.frame()
          names(row_annot)[y] <- annot_sel
        }
        row_annot<- as.data.frame(row_annot)
        rownames(row_annot) <- rownames(SettingsFile_Metab)
      }

      #Check number of features:
      Features <- as.data.frame(t(data_path))
      if(Enforce_FeatureNames==TRUE){
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
      if(Enforce_SampleNames==TRUE){
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
      if(nrow(t(data_path))>= 2){
        set.seed(1234)
        heatmap <- pheatmap::pheatmap(t(data_path),
                                     show_rownames = as.logical(show_rownames),
                                     show_colnames = as.logical(show_colnames),
                                     clustering_method =  "complete",
                                     scale = Scale,
                                     clustering_distance_rows = "correlation",
                                     annotation_col = col_annot,
                                     annotation_row = row_annot,
                                     legend = T,
                                     cellwidth = cellwidth_Sample,
                                     cellheight = cellheight_Feature,
                                     fontsize_row= 10,
                                     fontsize_col = 10,
                                     fontsize=9,
                                     main = paste(PlotName, " Metabolites: ", i, sep=" " ))

        ## Store the plot in the 'plots' list
        cleaned_i <- gsub("[[:space:],/\\\\]", "-", i)#removes empty spaces and replaces /,\ with -
        PlotList[[cleaned_i]] <- heatmap

        #Width and height according to Sample and metabolite number
        Plot_Sized <- PlotGrob_Heatmap(InputPlot=heatmap, SettingsInfo=SettingsInfo, PlotName= cleaned_i)
        PlotHeight <- as.numeric(Plot_Sized$height)
        PlotWidth <- as.numeric(Plot_Sized$width)
        Plot_Sized %<>%
          {ggplot2::ggplot() + annotation_custom(.)} %>%
          add(theme(panel.background = element_rect(fill = "transparent")))

        PlotList_adaptedGrid[[cleaned_i]] <- Plot_Sized

        #----- Save
        suppressMessages(suppressWarnings(
          MetaProViz:::SaveRes(InputList_DF=NULL,
                               InputList_Plot= PlotList_adaptedGrid,
                               SaveAs_Table=NULL,
                               SaveAs_Plot=SaveAs_Plot,
                               FolderPath= Folder,
                               FileName=paste("Heatmap_",PlotName, sep=""),
                               CoRe=FALSE,
                               PrintPlot=PrintPlot,
                               PlotHeight=PlotHeight,
                               PlotWidth=PlotWidth,
                               PlotUnit="cm")))

      }else{
          message(i , " includes <= 2 objects and is hence not plotted.")
        }
    }
    #Return if assigned:
    return(invisible(list("Plot"=PlotList,"Plot_Sized" = PlotList_adaptedGrid)))

    }else if("individual_Metab" %in% names(SettingsInfo)==FALSE & "individual_Sample" %in% names(SettingsInfo)==TRUE){
      #Ensure that groups that are assigned NAs do not cause problems:
      SettingsFile_Sample[[SettingsInfo[["individual_Sample"]]]] <-ifelse(is.na(SettingsFile_Sample[[SettingsInfo[["individual_Sample"]]]]), "NA", SettingsFile_Sample[[SettingsInfo[["individual_Sample"]]]])

      unique_paths_Sample <- unique(SettingsFile_Sample[[SettingsInfo[["individual_Sample"]]]])

      for (i in unique_paths_Sample){# Check pathways with 1 metabolite
        selected_path <- SettingsFile_Sample %>% filter(get(SettingsInfo[["individual_Sample"]]) == i)
        selected_path_metabs <-  colnames(data) [colnames(data) %in% row.names(selected_path)]
        if(length(selected_path_metabs)==1 ){
          warning("The metadata group ", paste(i), " includes only 1 metabolite. Heatmap cannot be made for 1 metabolite, thus it will be ignored.")
          unique_paths_Sample <- unique_paths_Sample[!unique_paths_Sample %in% i] # Remove the pathway
        }
      }

      IndividualPlots <-unique_paths_Sample
      PlotList <- list()#Empty list to store all the plots
      PlotList_adaptedGrid <- list()#Empty list to store all the plots

      for (i in IndividualPlots){
        #Select the data:
        selected_path <- SettingsFile_Sample %>% filter(get(SettingsInfo[["individual_Sample"]]) == i)%>%
          rownames_to_column("UniqueID")
        selected_path <- as.data.frame(selected_path[,1])%>%
          dplyr::rename("UniqueID"=1)
        data_path <- merge(selected_path, data%>%rownames_to_column("UniqueID"), by="UniqueID", all.x=TRUE)
        data_path <- data_path%>%
          column_to_rownames("UniqueID")

        # Column annotation
        selected_SettingsFile_Sample <- merge(selected_path, SettingsFile_Sample%>%rownames_to_column("UniqueID"), by="UniqueID", all.x=TRUE)

        col_annot_vars <- SettingsInfo[grepl("color_Sample", names(SettingsInfo))]
        col_annot<- NULL
        if(length(col_annot_vars)>0){
          for (x in 1:length(col_annot_vars)){
            annot_sel <- col_annot_vars[[x]]
            col_annot[x] <- selected_SettingsFile_Sample %>% select(annot_sel) %>% as.data.frame()
            names(col_annot)[x] <- annot_sel
          }
          col_annot<- as.data.frame(col_annot)
          rownames(col_annot) <- rownames(data_path)
        }

        # Row annotation
        row_annot_vars <- SettingsInfo[grepl("color_Metab", names(SettingsInfo))]
        row_annot<- NULL
        if(length(row_annot_vars)>0){
          for (y in 1:length(row_annot_vars)){
            annot_sel <- row_annot_vars[[y]]
            row_annot[y] <- SettingsFile_Metab %>% select(all_of(annot_sel))
            row_annot <- row_annot %>% as.data.frame()
            names(row_annot)[y] <- annot_sel
          }
          row_annot<- as.data.frame(row_annot)
          rownames(row_annot) <- rownames(SettingsFile_Metab)
        }

        #Check number of features:
        Features <- as.data.frame(t(data_path))
        if(Enforce_FeatureNames==TRUE){
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
        if(Enforce_SampleNames==TRUE){
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
        if(nrow(t(data_path))>= 2){
        set.seed(1234)
        heatmap <- pheatmap::pheatmap(t(data_path),
                                      show_rownames = as.logical(show_rownames),
                                      show_colnames = as.logical(show_colnames),
                                      clustering_method =  "complete",
                                      scale = Scale,
                                      clustering_distance_rows = "correlation",
                                      annotation_col = col_annot,
                                      annotation_row = row_annot,
                                      legend = T,
                                      cellwidth = cellwidth_Sample,
                                      cellheight = cellheight_Feature,
                                      fontsize_row= 10,
                                      fontsize_col = 10,
                                      fontsize=9,
                                      main = paste(PlotName," Samples: ", i, sep=" " ))

        #----- Save
        cleaned_i <- gsub("[[:space:],/\\\\]", "-", i)#removes empty spaces and replaces /,\ with -
        PlotList[[cleaned_i]] <- heatmap

        #Width and height according to Sample and metabolite number
        Plot_Sized <- PlotGrob_Heatmap(InputPlot=heatmap, SettingsInfo=SettingsInfo, PlotName= cleaned_i)
        PlotHeight <- as.numeric(Plot_Sized$height)
        PlotWidth <- as.numeric(Plot_Sized$width)
        Plot_Sized %<>%
          {ggplot2::ggplot() + annotation_custom(.)} %>%
          add(theme(panel.background = element_rect(fill = "transparent")))

        PlotList_adaptedGrid[[cleaned_i]] <- Plot_Sized

        #----- Save
        suppressMessages(suppressWarnings(
          MetaProViz:::SaveRes(InputList_DF=NULL,
                               InputList_Plot= PlotList_adaptedGrid,
                               SaveAs_Table=NULL,
                               SaveAs_Plot=SaveAs_Plot,
                               FolderPath= Folder,
                               FileName= paste("Heatmap_",PlotName, sep=""),
                               CoRe=FALSE,
                               PrintPlot=PrintPlot,
                               PlotHeight=PlotHeight,
                               PlotWidth=PlotWidth,
                               PlotUnit="cm")))
        }else{
          message(i , " includes <= 2 objects and is hence not plotted.")
        }
        }
      #Return if assigned:
      return(invisible(list("Plot"=PlotList,"Plot_Sized" = PlotList_adaptedGrid)))

      }else if("individual_Metab" %in% names(SettingsInfo)==TRUE & "individual_Sample" %in% names(SettingsInfo)==TRUE){
        #Ensure that groups that are assigned NAs do not cause problems:
        SettingsFile_Metab[[SettingsInfo[["individual_Metab"]]]] <-ifelse(is.na(SettingsFile_Metab[[SettingsInfo[["individual_Metab"]]]]), "NA", SettingsFile_Metab[[SettingsInfo[["individual_Metab"]]]])

        unique_paths <- unique(SettingsFile_Metab[[SettingsInfo[["individual_Metab"]]]])

        for (i in unique_paths){# Check pathways with 1 metabolite
          selected_path <- SettingsFile_Metab %>% filter(get(SettingsInfo[["individual_Metab"]]) == i)
          selected_path_metabs <-  colnames(data) [colnames(data) %in% row.names(selected_path)]
          if(length(selected_path_metabs)==1 ){
            warning("The metadata group ", paste(i), " includes only 1 metabolite. Heatmap cannot be made for 1 metabolite, thus it will be ignored.")
            unique_paths <- unique_paths[!unique_paths %in% i] # Remove the pathway
          }
        }

        #Ensure that groups that are assigned NAs do not cause problems:
        SettingsFile_Sample[[SettingsInfo[["individual_Sample"]]]] <-ifelse(is.na(SettingsFile_Sample[[SettingsInfo[["individual_Sample"]]]]), "NA", SettingsFile_Sample[[SettingsInfo[["individual_Sample"]]]])

        unique_paths_Sample <- unique(SettingsFile_Sample[[SettingsInfo[["individual_Sample"]]]])

        for (i in unique_paths_Sample){# Check pathways with 1 metabolite
          selected_path <- SettingsFile_Sample %>% filter(get(SettingsInfo[["individual_Sample"]]) == i)
          selected_path_metabs <-  colnames(data) [colnames(data) %in% row.names(selected_path)]
          if(length(selected_path_metabs)==1 ){
            warning("The metadata group ", paste(i), " includes only 1 metabolite. Heatmap cannot be made for 1 metabolite, thus it will be ignored.")
            unique_paths_Sample <- unique_paths_Sample[!unique_paths_Sample %in% i] # Remove the pathway
          }
        }

        IndividualPlots_Metab <-unique_paths
        IndividualPlots_Sample <-unique_paths_Sample

        PlotList <- list()#Empty list to store all the plots
        PlotList_adaptedGrid <- list()#Empty list to store all the plots

        for (i in IndividualPlots_Metab){
          selected_path <- SettingsFile_Metab %>% filter(get(SettingsInfo[["individual_Metab"]]) == i)
          selected_path_metabs <-  colnames(data) [colnames(data) %in% row.names(selected_path)]
          data_path_metab <- data %>% select(all_of(selected_path_metabs))

          # Row annotation
          row_annot_vars <- SettingsInfo[grepl("color_Metab", names(SettingsInfo))]
          row_annot<- NULL
          if(length(row_annot_vars)>0){
            for (y in 1:length(row_annot_vars)){
              annot_sel <- row_annot_vars[[y]]
              row_annot[y] <- SettingsFile_Metab %>% select(all_of(annot_sel))
              row_annot <- row_annot %>% as.data.frame()
              names(row_annot)[y] <- annot_sel
            }
            row_annot<- as.data.frame(row_annot)
            rownames(row_annot) <- rownames(SettingsFile_Metab)
          }

          #Col annotation:
          for (s in IndividualPlots_Sample){
            #Select the data:
            selected_path <- SettingsFile_Sample %>% filter(get(SettingsInfo[["individual_Sample"]]) == s)%>%
              rownames_to_column("UniqueID")
            selected_path <- as.data.frame(selected_path[,1])%>%
              dplyr::rename("UniqueID"=1)
            data_path <- merge(selected_path, data_path_metab%>%rownames_to_column("UniqueID"), by="UniqueID", all.x=TRUE)
            data_path <- data_path%>%
              column_to_rownames("UniqueID")

            # Column annotation
            selected_SettingsFile_Sample <- merge(selected_path, SettingsFile_Sample%>%rownames_to_column("UniqueID"), by="UniqueID", all.x=TRUE)

            col_annot_vars <- SettingsInfo[grepl("color_Sample", names(SettingsInfo))]
            col_annot<- NULL
            if(length(col_annot_vars)>0){
              for (x in 1:length(col_annot_vars)){
                annot_sel <- col_annot_vars[[x]]
                col_annot[x] <- selected_SettingsFile_Sample %>% select(annot_sel) %>% as.data.frame()
                names(col_annot)[x] <- annot_sel
              }
              col_annot<- as.data.frame(col_annot)
              rownames(col_annot) <- rownames(data_path)
            }

            #Check number of features:
            Features <- as.data.frame(t(data_path))
            if(Enforce_FeatureNames==TRUE){
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
            if(Enforce_SampleNames==TRUE){
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
            if(nrow(t(data_path))>= 2){
            set.seed(1234)
            heatmap <- pheatmap::pheatmap(t(data_path),
                                          show_rownames = as.logical(show_rownames),
                                          show_colnames = as.logical(show_colnames),
                                          clustering_method =  "complete",
                                          scale = Scale,
                                          clustering_distance_rows = "correlation",
                                          annotation_col = col_annot,
                                          annotation_row = row_annot,
                                          legend = T,
                                          cellwidth = cellwidth_Sample,
                                          cellheight = cellheight_Feature,
                                          fontsize_row= 10,
                                          fontsize_col = 10,
                                          fontsize=9,
                                          main = paste(PlotName," Metabolites: ", i, " Sample:", s, sep="" ))

            ## Store the plot in the 'plots' list
            cleaned_i <- gsub("[[:space:],/\\\\]", "-", i)#removes empty spaces and replaces /,\ with -
            cleaned_s <- gsub("[[:space:],/\\\\]", "-", s)#removes empty spaces and replaces /,\ with -
            PlotList[[paste(cleaned_i,cleaned_s, sep="_")]] <- heatmap

            #-------- Plot width and heights
            #Width and height according to Sample and metabolite number
            PlotName <- paste(cleaned_i,cleaned_s, sep="_")
            Plot_Sized <- PlotGrob_Heatmap(InputPlot=heatmap, SettingsInfo=SettingsInfo, PlotName= PlotName)
            PlotHeight <- as.numeric(Plot_Sized$height)
            PlotWidth <- as.numeric(Plot_Sized$width)
            Plot_Sized %<>%
              {ggplot2::ggplot() + annotation_custom(.)} %>%
              add(theme(panel.background = element_rect(fill = "transparent")))

            PlotList_adaptedGrid[[paste(cleaned_i,cleaned_s, sep="_")]] <- Plot_Sized

            #----- Save
            suppressMessages(suppressWarnings(
              MetaProViz:::SaveRes(InputList_DF=NULL,
                                   InputList_Plot= PlotList_adaptedGrid,
                                   SaveAs_Table=NULL,
                                   SaveAs_Plot=SaveAs_Plot,
                                   FolderPath= Folder,
                                   FileName=paste("Heatmap_",PlotName, sep=""),
                                   CoRe=FALSE,
                                   PrintPlot=PrintPlot,
                                   PlotHeight=PlotHeight,
                                   PlotWidth= PlotWidth,
                                   PlotUnit="cm")))


            }
            else{
              message(i , " includes <= 2 objects and is hence not plotted.")
            }
          }
        }
        return(invisible(list("Plot"=PlotList,"Plot_Sized" = PlotList_adaptedGrid)))
    } else if("individual_Metab" %in% names(SettingsInfo)==FALSE & "individual_Sample" %in% names(SettingsInfo)==FALSE){

    PlotList <- list()#Empty list to store all the plots
    PlotList_adaptedGrid <- list()#Empty list to store all the plots

    # Column annotation
    col_annot_vars <- SettingsInfo[grepl("color_Sample", names(SettingsInfo))]
    col_annot<- NULL
    if(length(col_annot_vars)>0){
      for (i in 1:length(col_annot_vars)){
        annot_sel <- col_annot_vars[[i]]
        col_annot[i] <- SettingsFile_Sample %>% select(annot_sel) %>% as.data.frame()
        names(col_annot)[i] <- annot_sel
      }
      col_annot<- as.data.frame(col_annot)
      rownames(col_annot) <- rownames(data)
    }

    # Row annotation
    row_annot_vars <- SettingsInfo[grepl("color_Metab", names(SettingsInfo))]
    row_annot<- NULL
    if(length(row_annot_vars)>0){
      for (i in 1:length(row_annot_vars)){
        annot_sel <- row_annot_vars[[i]]
        row_annot[i] <- SettingsFile_Metab %>% select(all_of(annot_sel))
        row_annot <- row_annot %>% as.data.frame()
        names(row_annot)[i] <- annot_sel
      }
      row_annot<- as.data.frame(row_annot)
      rownames(row_annot) <- rownames(SettingsFile_Metab)
    }

    #Check number of features:
    Features <- as.data.frame(t(data))
    if(Enforce_FeatureNames==TRUE){
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
    if(Enforce_SampleNames==TRUE){
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
    if(nrow(t(data))>= 2){
    set.seed(1234)
    heatmap <- pheatmap::pheatmap(t(data),
                                  show_rownames = as.logical(show_rownames),
                                  show_colnames = as.logical(show_colnames),
                                  clustering_method =  "complete",
                                  scale = Scale,
                                  clustering_distance_rows = "correlation",
                                  annotation_col = col_annot,
                                  annotation_row = row_annot,
                                  legend = T,
                                  cellwidth = cellwidth_Sample,
                                  cellheight = cellheight_Feature,
                                  fontsize_row= 10,
                                  fontsize_col = 10,
                                  fontsize=9,
                                  main = PlotName)

    ## Store the plot in the 'plots' list
    PlotList[[PlotName]] <- heatmap

    #-------- Plot width and heights
    #Width and height according to Sample and metabolite number
    Plot_Sized <- PlotGrob_Heatmap(InputPlot=heatmap, SettingsInfo=SettingsInfo, PlotName= PlotName)
    PlotHeight <- as.numeric(Plot_Sized$height)
    PlotWidth <- as.numeric(Plot_Sized$width)
    Plot_Sized %<>%
      {ggplot2::ggplot() + annotation_custom(.)} %>%
      add(theme(panel.background = element_rect(fill = "transparent")))

    PlotList_adaptedGrid[[PlotName]] <- Plot_Sized

    #----- Save
    suppressMessages(suppressWarnings(
      MetaProViz:::SaveRes(InputList_DF=NULL,
                           InputList_Plot= PlotList_adaptedGrid,
                           SaveAs_Table=NULL,
                           SaveAs_Plot=SaveAs_Plot,
                           FolderPath= Folder,
                           FileName= paste("Heatmap_",PlotName, sep=""),
                           CoRe=FALSE,
                           PrintPlot=PrintPlot,
                           PlotHeight=PlotHeight,
                           PlotWidth=PlotHeight,
                           PlotUnit="cm")))



    }else{
      message(i , " includes <= 2 objects and is hence not plotted.")
    }
    }
   return(invisible(list("Plot"=PlotList,"Plot_Sized" = PlotList_adaptedGrid)))
}


##############################################################
### ### ### Heatmap helper function: Internal Function ### ### ###
##############################################################

#' @param InputPlot This is the ggplot object generated within the VizHeatmap function.
#' @param SettingsInfo Passed to VizHeatmap
#' @param PlotName Passed to VizHeatmap
#'
#' @keywords Heatmap helper function
#' @noRd

PlotGrob_Heatmap <- function(InputPlot, SettingsInfo, PlotName){

  # Set the parameters for the plot we would like to use as a basis, before we start adjusting it:
  HEAT_PARAM <- list(
    widths = list(
      list("legend", "2cm")
    ),
    heights = list(
      list("main", "1cm")
    )
  )


  #Adjust the parameters:
  Input <- InputPlot$gtable

  Plot_Sized <- Input  %>%
    withCanvasSize(width = 12, height = 11) %>%
    adjust_layout(HEAT_PARAM) #%>%
    #adjust_title(c(PlotName)) %>%
    #adjust_legend(
    #  Input,
    #  sections = c("color_Sample", "color_Metab"),#How can we do this as here we have multiple once!
    #  SettingsInfo = SettingsInfo
    #)

  Output <- Plot_Sized

}



