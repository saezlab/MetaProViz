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
#' @param Plot_SettingsInfo  \emph{Optional: } NULL or Named vector  where you can include vectors or lists for annotation c(individual= "", color_Metab="ColumnName_Plot_SettingsFile_Metab", color_Sample= list("ColumnName_Plot_SettingsFile_Sample", "ColumnName_Plot_SettingsFile_Sample",...)).\strong{Default = NULL}
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
    if("individual" %in% names(Plot_SettingsInfo)==TRUE){
      if(Plot_SettingsInfo[["individual"]] %in% colnames(Plot_SettingsFile_Metab)== FALSE){
        stop("You have chosen individual = ",paste(Plot_SettingsInfo$individual), ", ", paste(Plot_SettingsInfo$individual)," does not exist in the Plot_SettingsFile_Metab"   )
      }
      # Check pathways with 1 metabolite
      unique_paths <- unique(Plot_SettingsFile_Metab[[Plot_SettingsInfo[["individual"]]]])

      for (i in unique_paths){
        selected_path <- Plot_SettingsFile_Metab %>% filter(get(Plot_SettingsInfo[["individual"]]) == i)
        selected_path_metabs <-  colnames(data) [colnames(data) %in% selected_path$Metabolite]

        if(length(selected_path_metabs)==1 ){
          warning("The pathway ", paste(i), " includes only 1 metabolite. Heatmap cannot be made for 1 metabolite, thus it will be ignored.")
          unique_paths <- unique_paths[!unique_paths %in% i] # Remove the pathway
        }
      }
    }

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

  data <- Input_data

  ## ------------ Create Output folders ----------- ##
  name <- paste0("MetaProViz_Results_",Sys.Date())
  WorkD <- getwd()
  Results_folder <- file.path(WorkD, name)
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)} # Make Results folder
  Results_folder_plots_Heatmaps_folder = file.path(Results_folder, "Heatmap")
  if (!dir.exists(Results_folder_plots_Heatmaps_folder)) {dir.create(Results_folder_plots_Heatmaps_folder)}  # check and create folder

  #####################################################
  ## -------------- Plot --------------- ##

  if("individual" %in% names(Plot_SettingsInfo)==TRUE){
    #individual_selection <-  Plot_SettingsInfo[["individual"]]
    # Create the list of individual plots that should be made:
    #IndividualPlots <- Plot_SettingsFile_Metab[!duplicated(Plot_SettingsFile_Metab[individual_selection]),]
    #IndividualPlots <- IndividualPlots[individual_selection]%>% unlist() %>% as.vector
    IndividualPlots <-unique_paths
    PlotList <- list()#Empty list to store all the plots

    for (i in IndividualPlots){
      # Select the data
      selected_path <- Plot_SettingsFile_Metab %>% filter(get(Plot_SettingsInfo[["individual"]]) == i)
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
                                    main = paste(OutputPlotName, i, sep=" " ))

      #Width and height according to Sample and metabolite number
      #Helper function to understand our plots parameters: gtable::gtable_show_layout(heatmap$gtable)
      #-------- Height
      heights <- heatmap[["gtable"]][["heights"]]

      mylist <- list() #create an empty list
      for (n in 1:length(heights)) {
        x <-as.character(heights[[n]])
        y <-regmatches(x, gregexpr("[0-9.]+",x))
        y <- paste(y[[1]], collapse = ", ")
        mylist[[n]] <-c(x,y)
      }

      height_DF <- as.data.frame(do.call("rbind",mylist)) #combine all vectors
      height_DF$Original <- height_DF$V1
      height_DF <- height_DF%>%
        mutate(V2 = strsplit(as.character(V2), ", ")) %>%
        unnest(V2) %>%
        mutate(V1 = strsplit(as.character(V1), ",")) %>%
        unnest(V1)
      height_DF$keep <- apply(height_DF, 1, function(row){
        grepl(row["V2"], row["V1"])
      })
      height_DF <- height_DF%>%
        mutate(Unit= case_when(grepl("bigpts", height_DF$V1)  ~ 'bigpts',
                               grepl("grobheight", height_DF$V1)  ~ 'grobheight',
                               TRUE ~ 'NA'))%>%
        subset(keep==TRUE, select = c(2,5))
      height_DF$Value <- as.numeric(height_DF$V2)
      height_DF<- height_DF[,c(3,2)] %>%
        group_by(Unit) %>%
        summarize(total_value = sum(Value))

      as.numeric(height_DF%>%subset(Unit=="grobheight", select=c(2)))

      plot_height <- (grid::convertX(unit(as.numeric(height_DF%>%subset(Unit=="bigpts", select=c(2))), "bigpts"), "cm", valueOnly = TRUE))+((grid::convertX(unit(as.numeric(height_DF%>%subset(Unit=="grobheight", select=c(2))), "npc"), "cm", valueOnly = TRUE))/9.2) #grobheight*3 (for annotation_row=TRUE), otherwise 1.2

      #-------- Width
      widths <- heatmap[["gtable"]][["widths"]]

      mylist <- list() #create an empty list
      for (m in 1:length(widths)) {
        x <-as.character(widths[[m]])
        y <-regmatches(x, gregexpr("[0-9.]+",x))
        y <- paste(y[[1]], collapse = ", ")
        mylist[[m]] <-c(x,y)
      }

      width_DF <- as.data.frame(do.call("rbind",mylist)) #combine all vectors
      width_DF$Original <- width_DF$V1
      width_DF <- width_DF%>%
        mutate(V2 = strsplit(as.character(V2), ", ")) %>%
        unnest(V2) %>%
        mutate(V1 = strsplit(as.character(V1), ",")) %>%
        unnest(V1)
      width_DF$keep <- apply(width_DF, 1, function(row){
        grepl(row["V2"], row["V1"])
      })
      width_DF <- width_DF%>%
        mutate(Unit= case_when(grepl("bigpts", width_DF$V1)  ~ 'bigpts',
                               grepl("grobwidth", width_DF$V1)  ~ 'grobwidth',
                               TRUE ~ 'NA'))%>%
        subset(keep==TRUE, select = c(2,5))
      width_DF$Value <- as.numeric(width_DF$V2)
      width_DF<- width_DF[,c(3,2)] %>%
        group_by(Unit) %>%
        summarize(total_value = sum(Value))

      as.numeric(width_DF%>%subset(Unit=="grobwidth", select=c(2)))

      plot_width <- (grid::convertX(unit(as.numeric(width_DF%>%subset(Unit=="bigpts", select=c(2))), "bigpts"), "cm", valueOnly = TRUE))+((grid::convertX(unit(as.numeric(width_DF%>%subset(Unit=="grobwidth", select=c(2))), "npc"), "cm", valueOnly = TRUE))/9.2)

      #--------Save
      cleaned_i <- gsub("[[:space:],/\\\\]", "-", i)#removes empty spaces and replaces /,\ with -
      ggsave(file=paste(Results_folder_plots_Heatmaps_folder,"/", "Heatmap_",cleaned_i,"_",OutputPlotName, ".",Save_as_Plot, sep=""), plot=heatmap, width=plot_width, height= plot_height, units = "cm")

    }

  }else if("individual" %in% names(Plot_SettingsInfo)==FALSE){

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

    #Width and height according to Sample and metabolite number
    #Helper function to understand our plots parameters: gtable::gtable_show_layout(heatmap$gtable)
    #-------- Height
    heights <- heatmap[["gtable"]][["heights"]]

    mylist <- list() #create an empty list
    for (i in 1:length(heights)) {
      x <-as.character(heights[[i]])
      y <-regmatches(x, gregexpr("[0-9.]+",x))
      y <- paste(y[[1]], collapse = ", ")
      mylist[[i]] <-c(x,y)
    }

    height_DF <- as.data.frame(do.call("rbind",mylist)) #combine all vectors
    height_DF$Original <- height_DF$V1
    height_DF <- height_DF%>%
      mutate(V2 = strsplit(as.character(V2), ", ")) %>%
      unnest(V2) %>%
      mutate(V1 = strsplit(as.character(V1), ",")) %>%
      unnest(V1)
    height_DF$keep <- apply(height_DF, 1, function(row){
      grepl(row["V2"], row["V1"])
    })
    height_DF <- height_DF%>%
      mutate(Unit= case_when(grepl("bigpts", height_DF$V1)  ~ 'bigpts',
                             grepl("grobheight", height_DF$V1)  ~ 'grobheight',
                             TRUE ~ 'NA'))%>%
      subset(keep==TRUE, select = c(2,5))
    height_DF$Value <- as.numeric(height_DF$V2)
    height_DF<- height_DF[,c(3,2)] %>%
      group_by(Unit) %>%
      summarize(total_value = sum(Value))

    as.numeric(height_DF%>%subset(Unit=="grobheight", select=c(2)))

    plot_height <- (grid::convertX(unit(as.numeric(height_DF%>%subset(Unit=="bigpts", select=c(2))), "bigpts"), "cm", valueOnly = TRUE))+((grid::convertX(unit(as.numeric(height_DF%>%subset(Unit=="grobheight", select=c(2))), "npc"), "cm", valueOnly = TRUE))/9.2) #grobheight*3 (for annotation_row=TRUE), otherwise 1.2

    #-------- Width
    widths <- heatmap[["gtable"]][["widths"]]

    mylist <- list() #create an empty list
    for (i in 1:length(widths)) {
      x <-as.character(widths[[i]])
      y <-regmatches(x, gregexpr("[0-9.]+",x))
      y <- paste(y[[1]], collapse = ", ")
      mylist[[i]] <-c(x,y)
    }

    width_DF <- as.data.frame(do.call("rbind",mylist)) #combine all vectors
    width_DF$Original <- width_DF$V1
    width_DF <- width_DF%>%
      mutate(V2 = strsplit(as.character(V2), ", ")) %>%
      unnest(V2) %>%
      mutate(V1 = strsplit(as.character(V1), ",")) %>%
      unnest(V1)
    width_DF$keep <- apply(width_DF, 1, function(row){
      grepl(row["V2"], row["V1"])
    })
    width_DF <- width_DF%>%
      mutate(Unit= case_when(grepl("bigpts", width_DF$V1)  ~ 'bigpts',
                             grepl("grobwidth", width_DF$V1)  ~ 'grobwidth',
                             TRUE ~ 'NA'))%>%
      subset(keep==TRUE, select = c(2,5))
    width_DF$Value <- as.numeric(width_DF$V2)
    width_DF<- width_DF[,c(3,2)] %>%
      group_by(Unit) %>%
      summarize(total_value = sum(Value))

    as.numeric(width_DF%>%subset(Unit=="grobwidth", select=c(2)))

    plot_width <- (grid::convertX(unit(as.numeric(width_DF%>%subset(Unit=="bigpts", select=c(2))), "bigpts"), "cm", valueOnly = TRUE))+((grid::convertX(unit(as.numeric(width_DF%>%subset(Unit=="grobwidth", select=c(2))), "npc"), "cm", valueOnly = TRUE))/9.2)

    #--------Save
    ggsave(file=paste(Results_folder_plots_Heatmaps_folder,"/", "Heatmap",OutputPlotName, ".", Save_as_Plot ,sep=""), plot=heatmap, width=plot_width, height= plot_height, units = "cm")
  }
}
