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


#################################
### ### ### PCA Plots ### ### ###
#################################

#' @param Input_data DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected. includes experimental design and outlier column.
#' @param Experimental_design DF which contains information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames. Column "Conditions" with information about the sample conditions (e.g. "N" and "T" or "Normal" and "Tumor"), can be used for feature filtering and colour coding in the PCA. Column "AnalyticalReplicate" including numerical values, defines technical repetitions of measurements, which will be summarised. Column "BiologicalReplicates" including numerical values. Please use the following names: "Conditions", "Biological_Replicates", "Analytical_Replicates".
#' @param Color \emph{Optional: }String which contains the name of the output file of the Metabolic Clusters
#' @param Shape \emph{Optional: }String which contains the name of the output file of the Metabolic Clusters
#' @param Show_Loadings  \emph{Optional: } TRUE or FALSE for whether PCA loadings are also plotted on the PCA (biplot) \strong{Default = FALSE}
#' @param Scaling  \emph{Optional: } TRUE or FALSE for whether a data scaling is used \strong{Default = TRUE}
# @param  Palette \emph{Optional: } ??
#' @param Theme \emph{Optional: } Selection of theme for plots from ggplot2. \strong{Default = theme_classic} ??
#' @param OutputPlotName \emph{Optional: } String which is added to the output files of the PCA
#' @param Save_as \emph{Optional: } Select the file type of output plots. Options are svg, png, pdf, jpeg, tiff, bmp. \strong{Default = svg}
#'
#' @keywords PCA
#' @export

VizPCA <- function(Input_data,
                   Experimental_design,
                   Color = FALSE,
                   Shape = FALSE,
                   Show_Loadings = FALSE,
                   Scaling = TRUE,
                   #Palette= "Set2".
                   Theme=theme_classic(),
                   OutputPlotName= '',
                   Save_as = "svg"
                  ){

  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse","ggfortify", "ggplot2")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  suppressMessages(library(tidyverse))
  suppressMessages(library("ggfortify"))

  ## ------------ Check Input files ----------- ##
  #1. Input_data

  if(any(duplicated(row.names(Input_data)))==TRUE){
    stop("Duplicated row.names of Input_data, whilst row.names must be unique")
  } else if("Conditions" %in% colnames(Experimental_design)==FALSE){
    stop("There is no column named `Conditions` in Experimental_design to obtain Condition1 and Condition2.")
  } else{
    Test_num <- apply(Input_data, 2, function(x) is.numeric(x))
    if((any(Test_num) ==  FALSE) ==  TRUE){
      stop("Input_data needs to be of class numeric")
    } else{
      Test_match <- merge(Experimental_design, Input_data, by.x = "row.names",by.y = "row.names", all =  FALSE) # Do the unique IDs of the "Input_data" match the row names of the "Experimental_design"?
      if(nrow(Test_match) ==  0){
        stop("row.names Input_data need to match row.names Experimental_design.")
      } else(
        Input_data <- Input_data
      )
    }
    Design <- Experimental_design
  }

  #2. Parameters
  if(Color != FALSE & Color %in% colnames(Design)==FALSE){
    stop(paste("PCA with Color was selected. However, there is no column named: " ,Color," in Input_data.",sep = "") )
  }
  if(Shape != FALSE & Shape %in% colnames(Design)==FALSE){
    stop(paste("PCA with Shapes was selected. However, there is no column named: " ,Shape," in Input_data.",sep = "") )
  }
  if(length(unique(Design[,Shape]))>6){
    stop("Error. You tried to plot more than 6 shapes. It would be preferable to use color instead of shape" )
  }
  if(is.logical(Show_Loadings) == FALSE){
    stop("Check input. The Show_Loadings value should be either =TRUE if loadings are to also be shown in the PCA plot or = FALSE if not.")
  }
  if(is.logical(Scaling) == FALSE){
    stop("Check input. The Scaling value should be either =TRUE if data scaling is to be performed prior to the PCA or = FALSE if not.")
  }
  # Theme ???
  Save_as_options <- c("svg","pdf", "jpeg", "tiff", "png", "bmp", "wmf","eps", "ps", "tex" )
  if(Save_as %in% Save_as_options == FALSE){
    stop("Check input. The selected Save_as option is not valid. Please select one of the folowwing: ",paste(Save_as_options,collapse = ", "),"." )
  }

  ## ------------ Create Output folders ----------- ##
  name <- paste0("MetaProViz_Results_",Sys.Date())
  WorkD <- getwd()
  Results_folder <- file.path(WorkD, name)
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)} # Make Results folder
  Results_folder_plots_PCA_folder = file.path(Results_folder, "PCA")  # This searches for a folder called "Preprocessing" within the "Results" folder in the current working directory and if its not found it creates one
  if (!dir.exists(Results_folder_plots_PCA_folder)) {dir.create(Results_folder_plots_PCA_folder)}  # check and create folder

  ### select plot based on arguments
  #1
  if (Color != FALSE & Shape != FALSE & Show_Loadings == TRUE & Scaling == TRUE){
    message(paste("Selected option for PCA: Color = ",Color,", Shape = ",Shape,", Show_Loadings = TRUE, Scaling = TRUE",sep = ""))

    mdata <- merge(Input_data, Design, by=0) # merge the data with the design to get only the kept samples
    mdata <- column_to_rownames(mdata, "Row.names")

    data <- merge(Input_data, mdata %>% select(all_of(Color)), by = 0) # Add the selected columns(the ones wanted to plot) and add them to the data
    data <- column_to_rownames(data, "Row.names")

    data <- merge(data, mdata %>% select(all_of(Shape)), by = 0) # merge the data with the design to get only the kept samples
    data <- column_to_rownames(data, "Row.names")
    data[,Shape] <- as.factor(data[,Shape]) # make the shape into a factor to be discrete

    # For color if we have character we go with discrete. If we have numeric we to discrete until 4groups. if e have more we go for continuous
    if(is.numeric(data[,Color]) == TRUE | is.integer(data[,Color]) == TRUE){
      if(length(unique(data[,Color])) > 4){ # change this to change the number after which color becomes from distinct to continuous
        data[,Color] <- as.numeric(data[,Color])
      }else(data[,Color] <- as.factor(data[,Color]))
    }

    PCA <- autoplot(prcomp(as.matrix(data %>% select(-all_of(c(Shape, Color)))),scale. = TRUE),   # Run and plot PCA
                   data= data,
                   colour = Color,
                   fill = Color,
                   shape = Shape,
                   size = 3,
                   alpha = 0.8,
                   label=T,
                   label.size=2.5,
                   label.repel = TRUE,
                   loadings=T, #draws Eigenvectors
                   loadings.label = TRUE,
                   loadings.label.vjust = 1.2,
                   loadings.label.size=2.5,
                   loadings.colour="grey10",
                   loadings.label.colour="grey10") +
      scale_shape_manual(values=c(22,21,24,23,25,7,8,11,12))+ #needed if more than 6 shapes are in place
      ggtitle(paste(OutputPlotName)) +
      Theme +
      geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
      geom_vline(xintercept=0,  color = "black", linewidth=0.1)

    loading_data <- prcomp(as.matrix(data %>% select(-all_of(c(Shape, Color)))),scale. = TRUE)

    #2
  }
  else if(Color != FALSE & Shape != FALSE & Show_Loadings == TRUE & Scaling == FALSE){
    message(paste("Selected option for PCA: Color = ",Color,", Shape = ",Shape,", Show_Loadings = TRUE, Scaling = FALSE",sep = ""))

    mdata <- merge(Input_data, Design, by=0) # merge the data with the design to get only the kept samples
    mdata <- column_to_rownames(mdata, "Row.names")

    data <- merge(Input_data, mdata %>% select(all_of(Color)), by = 0) # Add the selected columns(the ones wanted to plot) and add them to the dat
    data <- column_to_rownames(data, "Row.names")

    data <- merge(data, mdata %>% select(all_of(Shape)), by = 0) # merge the data with the design to get only the kept samples
    data <- column_to_rownames(data, "Row.names")
    data[,Shape] <- as.factor(data[,Shape]) # make the shape into a factor to be discrete

    # For color if we have character we go with discrete. If we have numeric we to discrete until 4groups. if e have more we go for continuous
    if(is.numeric(data[,Color]) == TRUE | is.integer(data[,Color]) == TRUE){
      if(length(unique(data[,Color])) > 4){ # change this to change the number after which color becomes from distinct to continuous
        data[,Color] <- as.numeric(data[,Color])
      }else(data[,Color] <- as.factor(data[,Color]))
    }


    PCA <- autoplot(prcomp(as.matrix(data %>% select(-all_of(c(Shape, Color)))),scale. = FALSE),   # Run and plot PCA
                   data= data,
                   colour = Color,
                   fill = Color,
                   shape = Shape,
                   size = 3,
                   alpha = 0.8,
                   label = T,
                   label.size =2.5,
                   label.repel = TRUE,
                   loadings = T, #draws Eigenvectors
                   loadings.label = TRUE,
                   loadings.label.vjust = 1.2,
                   loadings.label.size = 2.5,
                   loadings.colour = "grey10",
                   loadings.label.colour = "grey10") +
      scale_shape_manual(values = c(22,21,24,23,25,7,8,11,12)) + #needed if more than 6 shapes are in place
      ggtitle(paste(OutputPlotName)) +
      Theme +
      geom_hline(yintercept = 0,  color = "black", linewidth = 0.1) +
      geom_vline(xintercept = 0,  color = "black", linewidth = 0.1)

    loading_data <- prcomp(as.matrix(data%>% select(-all_of(c(Shape, Color)))),scale. = FALSE)

    #3 Done
  }
  else if(Color != FALSE & Shape != FALSE & Show_Loadings == FALSE & Scaling== TRUE){
    message(paste("Selected option for PCA: Color = ",Color,", Shape = ",Shape,", Show_Loadings = FALSE, Scaling = TRUE",sep=""))

    mdata<- merge(Input_data, Design, by=0) # merge the data with the design to get only the kept samples
    mdata <- column_to_rownames(mdata, "Row.names")

    data<- merge(Input_data, mdata %>%select(all_of(Color)), by=0) # Add the selected columns(the ones wanted to plot) and add them to the dat
    data <- column_to_rownames(data, "Row.names")

    data<- merge(data, mdata %>%select(all_of(Shape)), by=0) # merge the data with the design to get only the kept samples
    data <- column_to_rownames(data, "Row.names")
    data[,Shape] <- as.factor(data[,Shape]) # make the shape into a factor to be discrete

    # For color if we have character we go with disctere. If we have numeriic we to discrete untill 4groupd. if e have more we go for continouus
    if(is.numeric(data[,Color])==TRUE | is.integer(data[,Color])==TRUE){
      if(length(unique(data[,Color]))>4){ # change this to change the number after which color becomes from distinct to continuous
        data[,Color] <- as.numeric(data[,Color])
      }else(data[,Color] <- as.factor(data[,Color]))
    }

    PCA<- autoplot(prcomp(as.matrix(data%>% select(-all_of(c(Shape, Color)))),scale. = TRUE),   # Run and plot PCA
                   data= data,
                   colour = Color,
                   fill = Color,
                   shape = Shape,
                   size = 3,
                   alpha = 0.8,
                   label=T,
                   label.size=2.5,
                   label.repel = TRUE) +
      scale_shape_manual(values=c(22,21,24,23,25,7,8,11,12))+ #needed if more than 6 shapes are in place
      ggtitle(paste(OutputPlotName)) +
      Theme +
      geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
      geom_vline(xintercept=0,  color = "black", linewidth=0.1)

    loading_data <- prcomp(as.matrix(data%>% select(-all_of(c(Shape, Color)))),scale. = TRUE)

    #4 Done
  }
  else if(Color != FALSE & Shape != FALSE & Show_Loadings == FALSE & Scaling== FALSE){
    message(paste("Selected option for PCA: Color = ",Color,", Shape = ",Shape,", Show_Loadings = FALSE, Scaling = FALSE",sep=""))

    mdata<- merge(Input_data, Design, by=0) # merge the data with the design to get only the kept samples
    mdata <- column_to_rownames(mdata, "Row.names")

    data<- merge(Input_data, mdata %>%select(all_of(Color)), by=0) # Add the selected columns(the ones wanted to plot) and add them to the dat
    data <- column_to_rownames(data, "Row.names")

    data<- merge(data, mdata %>%select(all_of(Shape)), by=0) # merge the data with the design to get only the kept samples
    data <- column_to_rownames(data, "Row.names")
    data[,Shape] <- as.factor(data[,Shape]) # make the shape into a factor to be discrete

    # For color if we have character we go with disctere. If we have numeriic we to discrete untill 4 groups. if e have more we go for continouus
    if(is.numeric(data[,Color])==TRUE | is.integer(data[,Color])==TRUE){
      if(length(unique(data[,Color]))>4){ # change this to change the number after which color becomes from distinct to continuous
        data[,Color] <- as.numeric(data[,Color])
      }else(data[,Color] <- as.factor(data[,Color]))
    }


    PCA<- autoplot(prcomp(as.matrix(data%>% select(-all_of(c(Shape, Color)))),scale. = FALSE),   # Run and plot PCA
                   data= data,
                   colour = Color,
                   fill = Color,
                   shape = Shape,
                   size = 3,
                   alpha = 0.8,
                   label=T,
                   label.size=2.5,
                   label.repel = TRUE) +
      scale_shape_manual(values=c(22,21,24,23,25,7,8,11,12))+ #needed if more than 6 shapes are in place
      ggtitle(paste(OutputPlotName)) +
      Theme +
      geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
      geom_vline(xintercept=0,  color = "black", linewidth=0.1)

    loading_data <- prcomp(as.matrix(data%>% select(-all_of(c(Shape, Color)))),scale. = FALSE)

    #5 Done
  }
  else if (Color == FALSE & Shape != FALSE & Show_Loadings == TRUE & Scaling== TRUE){
    message(paste("Selected option for PCA: Color = FALSE, Shape = ",Shape,", Show_Loadings = TRUE, Scaling = TRUE",sep=""))

    mdata<- merge(Input_data, Design, by=0) # merge the data with the design to get only the kept samples
    mdata <- column_to_rownames(mdata, "Row.names")

    data<- merge(Input_data, mdata %>%select(all_of(Shape)), by=0) # merge the data with the design to get only the kept samples
    data <- column_to_rownames(data, "Row.names")
    data[,Shape] <- as.factor(data[,Shape]) # make the shape into a factor to be discrete

    PCA<- autoplot(prcomp(as.matrix(data%>% select(-all_of(Shape))),scale. = TRUE),   # Run and plot PCA
                   data= data,
                   shape = Shape,
                   size = 3,
                   alpha = 0.8,
                   label=T,
                   label.size=2.5,
                   label.repel = TRUE,
                   loadings=T, #draws Eigenvectors
                   loadings.label = TRUE,
                   loadings.label.vjust = 1.2,
                   loadings.label.size=2.5,
                   loadings.colour="grey10",
                   loadings.label.colour="grey10") +
      scale_shape_manual(values=c(22,21,24,23,25,7,8,11,12))+ #needed if more than 6 shapes are in place
      ggtitle(paste(OutputPlotName)) +
      Theme +
      geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
      geom_vline(xintercept=0,  color = "black", linewidth=0.1)

    loading_data <- prcomp(as.matrix(data%>% select(-all_of( Shape))),scale. = TRUE)

    #6 Done
  }
  else if(Color == FALSE & Shape != FALSE & Show_Loadings == TRUE & Scaling== FALSE){
    message(paste("Selected option for PCA: Color = FALSE, Shape = ",Shape,", Show_Loadings = TRUE, Scaling = FALSE",sep=""))

    mdata<- merge(Input_data, Design, by=0) # merge the data with the design to get only the kept samples
    mdata <- column_to_rownames(mdata, "Row.names")

    data<- merge(Input_data, mdata %>%select(all_of(Shape)), by=0) # merge the data with the design to get only the kept samples
    data <- column_to_rownames(data, "Row.names")
    data[,Shape] <- as.factor(data[,Shape]) # make the shape into a factor to be discrete

    PCA<- autoplot(prcomp(as.matrix(data%>% select(-all_of(Shape))),scale. = FALSE),   # Run and plot PCA
                   data= data,
                   shape = Shape,
                   size = 3,
                   alpha = 0.8,
                   label=T,
                   label.size=2.5,
                   label.repel = TRUE,
                   loadings=T, #draws Eigenvectors
                   loadings.label = TRUE,
                   loadings.label.vjust = 1.2,
                   loadings.label.size=2.5,
                   loadings.colour="grey10",
                   loadings.label.colour="grey10") +
      scale_shape_manual(values=c(22,21,24,23,25,7,8,11,12))+ #needed if more than 6 shapes are in place
      ggtitle(paste(OutputPlotName)) +
      Theme +
      geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
      geom_vline(xintercept=0,  color = "black", linewidth=0.1)

    loading_data <- prcomp(as.matrix(data%>% select(-all_of( Shape))),scale. = FALSE)

    #7 Done
  }
  else if(Color == FALSE & Shape != FALSE & Show_Loadings == FALSE & Scaling== TRUE){
    message(paste("Selected option for PCA: Color = FALSE, Shape = ",Shape,", Show_Loadings = FALSE, Scaling = TRUE",sep=""))

    mdata<- merge(Input_data, Design, by=0) # merge the data with the design to get only the kept samples
    mdata <- column_to_rownames(mdata, "Row.names")

    data<- merge(Input_data, mdata %>%select(all_of(Shape)), by=0) # merge the data with the design to get only the kept samples
    data <- column_to_rownames(data, "Row.names")
    data[,Shape] <- as.factor(data[,Shape]) # make the shape into a factor to be discrete

    PCA<- autoplot(prcomp(as.matrix(data%>% select(-all_of(Shape))),scale. = TRUE),   # Run and plot PCA
                   data= data,
                   shape = Shape,
                   size = 3,
                   alpha = 0.8,
                   label=T,
                   label.size=2.5,
                   label.repel = TRUE) +
      scale_shape_manual(values=c(22,21,24,23,25,7,8,11,12))+ #needed if more than 6 shapes are in place
      ggtitle(paste(OutputPlotName)) +
      Theme +
      geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
      geom_vline(xintercept=0,  color = "black", linewidth=0.1)

    loading_data <- prcomp(as.matrix(data%>% select(-all_of( Shape))),scale. = TRUE)

    #8 Done
  }
  else if(Color == FALSE & Shape != FALSE & Show_Loadings == FALSE & Scaling== FALSE){
    message(paste("Selected option for PCA: Color = FALSE, Shape = ",Shape,", Show_Loadings = FALSE, Scaling = FALSE",sep=""))

    mdata<- merge(Input_data, Design, by=0) # merge the data with the design to get only the kept samples
    mdata <- column_to_rownames(mdata, "Row.names")

    data<- merge(Input_data, mdata %>%select(all_of(Shape)), by=0) # merge the data with the design to get only the kept samples
    data <- column_to_rownames(data, "Row.names")
    data[,Shape] <- as.factor(data[,Shape]) # make the shape into a factor to be discrete

    PCA<- autoplot(prcomp(as.matrix(data%>% select(-all_of(Shape))),scale. = FALSE),   # Run and plot PCA
                   data= data,
                   shape = Shape,
                   size = 3,
                   alpha = 0.8,
                   label=T,
                   label.size=2.5,
                   label.repel = TRUE) +
      scale_shape_manual(values=c(22,21,24,23,25,7,8,11,12))+ #needed if more than 6 shapes are in place
      ggtitle(paste(OutputPlotName)) +
      Theme +
      geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
      geom_vline(xintercept=0,  color = "black", linewidth=0.1)

    loading_data <- prcomp(as.matrix(data%>% select(-all_of( Shape))),scale. = FALSE)

    #9 Done
  }
  else if (Color != FALSE & Shape == FALSE & Show_Loadings == TRUE & Scaling== TRUE){
    message(paste("Selected option for PCA: Color = ",Color,", Shape = FALSE, Show_Loadings = TRUE, Scaling = TRUE",sep=""))

    mdata<- merge(Input_data, Design, by=0) # merge the data with the design to get only the kept samples
    mdata <- column_to_rownames(mdata, "Row.names")

    data<- merge(Input_data, mdata %>%select(all_of(Color)), by=0) # Add the selected columns(the ones wanted to plot) and add them to the dat
    data <- column_to_rownames(data, "Row.names")

    # For color if we have character we go with disctere. If we have numeriic we to discrete untill 4groupd. if e have more we go for continouus
    if(is.numeric(data[,Color])==TRUE | is.integer(data[,Color])==TRUE){
      if(length(unique(data[,Color]))>4){ # change this to change the number after which color becomes from distinct to continuous
        data[,Color] <- as.numeric(data[,Color])
      }else(data[,Color] <- as.factor(data[,Color]))
    }

    PCA<- autoplot(prcomp(as.matrix(data%>% select(-all_of(Color))),scale. = TRUE),   # Run and plot PCA
                   data= data,
                   colour = Color,
                   size = 3,
                   alpha = 0.8,
                   label=T,
                   label.size=2.5,
                   label.repel = TRUE,
                   loadings=T, #draws Eigenvectors
                   loadings.label = TRUE,
                   loadings.label.vjust = 1.2,
                   loadings.label.size=2.5,
                   loadings.colour="grey10",
                   loadings.label.colour="grey10") +
      ggtitle(paste(OutputPlotName)) +
      Theme +
      geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
      geom_vline(xintercept=0,  color = "black", linewidth=0.1)

    loading_data <- prcomp(as.matrix(data%>% select(-all_of( Color))),scale. = TRUE)

    #10 Done
  }
  else if(Color != FALSE & Shape == FALSE & Show_Loadings == TRUE & Scaling== FALSE){
    message(paste("Selected option for PCA: Color = ",Color,", Shape = FALSE, Show_Loadings = TRUE, Scaling = FALSE",sep=""))

    mdata<- merge(Input_data, Design, by=0) # merge the data with the design to get only the kept samples
    mdata <- column_to_rownames(mdata, "Row.names")

    data<- merge(Input_data, mdata %>%select(all_of(Color)), by=0) # Add the selected columns(the ones wanted to plot) and add them to the dat
    data <- column_to_rownames(data, "Row.names")

    # For color if we have character we go with disctere. If we have numeriic we to discrete untill 4groupd. if e have more we go for continouus
    if(is.numeric(data[,Color])==TRUE | is.integer(data[,Color])==TRUE){
      if(length(unique(data[,Color]))>4){ # change this to change the number after which color becomes from distinct to continuous
        data[,Color] <- as.numeric(data[,Color])
      }else(data[,Color] <- as.factor(data[,Color]))
    }

    PCA<- autoplot(prcomp(as.matrix(data%>% select(-all_of(Color))),scale. = FALSE),   # Run and plot PCA
                   data= data,
                   colour = Color,
                   size = 3,
                   alpha = 0.8,
                   label=T,
                   label.size=2.5,
                   label.repel = TRUE,
                   loadings=T, #draws Eigenvectors
                   loadings.label = TRUE,
                   loadings.label.vjust = 1.2,
                   loadings.label.size=2.5,
                   loadings.colour="grey10",
                   loadings.label.colour="grey10") +
      ggtitle(paste(OutputPlotName)) +
      Theme +
      geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
      geom_vline(xintercept=0,  color = "black", linewidth=0.1)

    loading_data <- prcomp(as.matrix(data%>% select(-all_of( Color))),scale. = FALSE)

    #11 Done
  }
  else if(Color != FALSE & Shape == FALSE & Show_Loadings == FALSE & Scaling== TRUE){
    message(paste("Selected option for PCA: Color = ",Color,", Shape = FALSE, Show_Loadings = TRUE, Scaling = TRUE",sep=""))

    mdata<- merge(Input_data, Design, by=0) # merge the data with the design to get only the kept samples
    mdata <- column_to_rownames(mdata, "Row.names")

    data<- merge(Input_data, mdata %>%select(all_of(Color)), by=0) # Add the selected columns(the ones wanted to plot) and add them to the dat
    data <- column_to_rownames(data, "Row.names")

    # For color if we have character we go with disctere. If we have numeriic we to discrete untill 4groupd. if e have more we go for continouus
    if(is.numeric(data[,Color])==TRUE | is.integer(data[,Color])==TRUE){
      if(length(unique(data[,Color]))>4){ # change this to change the number after which color becomes from distinct to continuous
        data[,Color] <- as.numeric(data[,Color])
      }else(data[,Color] <- as.factor(data[,Color]))
    }

    PCA<- autoplot(prcomp(as.matrix(data%>% select(-all_of(Color))),scale. = TRUE),   # Run and plot PCA
                   data= data,
                   colour = Color,
                   size = 3,
                   alpha = 0.8,
                   label=T,
                   label.size=2.5,
                   label.repel = TRUE) +
      ggtitle(paste(OutputPlotName)) +
      Theme +
      geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
      geom_vline(xintercept=0,  color = "black", linewidth=0.1)

    loading_data <- prcomp(as.matrix(data%>% select(-all_of( Color))),scale. = TRUE)

    #12 Done
  }
  else if(Color != FALSE & Shape == FALSE & Show_Loadings == FALSE & Scaling== FALSE){
    message(paste("Selected option for PCA: Color = ",Color,", Shape = FALSE, Show_Loadings = FALSE, Scaling = FALSE",sep=""))

    mdata<- merge(Input_data, Design, by=0) # merge the data with the design to get only the kept samples
    mdata <- column_to_rownames(mdata, "Row.names")
    data<- merge(Input_data, mdata %>%select(all_of(Color)), by=0) # Add the selected columns(the ones wanted to plot) and add them to the dat
    data <- column_to_rownames(data, "Row.names")

    # For color if we have character we go with disctere. If we have numeriic we to discrete untill 4groupd. if e have more we go for continouus
    if(is.numeric(data[,Color])==TRUE | is.integer(data[,Color])==TRUE){
      if(length(unique(data[,Color]))>4){ # change this to change the number after which color becomes from distinct to continuous
        data[,Color] <- as.numeric(data[,Color])
      }else(data[,Color] <- as.factor(data[,Color]))
    }

    PCA<- autoplot(prcomp(as.matrix(data%>% select(-all_of(Color))),scale. = FALSE),   # Run and plot PCA
                   data= data,
                   colour = Color,
                   size = 3,
                   alpha = 0.8,
                   label=T,
                   label.size=2.5,
                   label.repel = TRUE) +
      ggtitle(paste(OutputPlotName)) +
      Theme +
      geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
      geom_vline(xintercept=0,  color = "black", linewidth=0.1)

    loading_data <- prcomp(as.matrix(data%>% select(-all_of( Color))),scale. = FALSE)

    #13 Done
  }
  else if (Color == FALSE & Shape == FALSE & Show_Loadings == TRUE & Scaling== TRUE){
    message(paste("Selected option for PCA: Color = FALSE, Shape = FALSE, Show_Loadings = TRUE, Scaling = TRUE"))

    data<- Input_data
    PCA<- autoplot(prcomp(as.matrix(data),scale. = TRUE),   # Run and plot PCA
                   data= data,
                   size = 3,
                   alpha = 0.8,
                   label=T,
                   label.size=2.5,
                   label.repel = TRUE,
                   loadings=T, #draws Eigenvectors
                   loadings.label = TRUE,
                   loadings.label.vjust = 1.2,
                   loadings.label.size=2.5,
                   loadings.colour="grey10",
                   loadings.label.colour="grey10") +
      ggtitle(paste(OutputPlotName)) +
      Theme +
      geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
      geom_vline(xintercept=0,  color = "black", linewidth=0.1)

    loading_data <- prcomp(as.matrix(data),scale. = TRUE)

    #14 Done
  }
  else if(Color == FALSE & Shape == FALSE & Show_Loadings == TRUE & Scaling== FALSE){
    message(paste("Selected option for PCA: Color = FALSE, Shape = FALSE, Show_Loadings = TRUE, Scaling = FALSE"))

    data<- Input_data
    PCA<- autoplot(prcomp(as.matrix(data),scale. = FALSE),   # Run and plot PCA
                   data= data,
                   size = 3,
                   alpha = 0.8,
                   label=T,
                   label.size=2.5,
                   label.repel = TRUE,
                   loadings=T, #draws Eigenvectors
                   loadings.label = TRUE,
                   loadings.label.vjust = 1.2,
                   loadings.label.size=2.5,
                   loadings.colour="grey10",
                   loadings.label.colour="grey10") +
      ggtitle(paste(OutputPlotName)) +
      Theme +
      geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
      geom_vline(xintercept=0,  color = "black", linewidth=0.1)

    loading_data <- prcomp(as.matrix(data),scale. = FALSE)

    #15 Done
  }
  else if(Color == FALSE &  Shape == FALSE & Show_Loadings == FALSE & Scaling== TRUE){
    message(paste("Selected option for PCA: Color = FALSE, Shape = FALSE, Show_Loadings = FALSE, Scaling = TRUE"))

    data<- Input_data
    PCA<- autoplot(prcomp(as.matrix(data),scale. = TRUE),   # Run and plot PCA
                   data= data,
                   size = 3,
                   alpha = 0.8,
                   label=T,
                   label.size=2.5,
                   label.repel = TRUE) +
      ggtitle(paste(OutputPlotName)) +
      Theme +
      geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
      geom_vline(xintercept=0,  color = "black", linewidth=0.1)

    loading_data <- prcomp(as.matrix(data),scale. = TRUE)

    #16 Done
  }
  else if(Color == FALSE &  Shape == FALSE & Show_Loadings == FALSE & Scaling== FALSE){
    message(paste("Selected option for PCA: Color = FALSE, Shape = FALSE, Show_Loadings = FALSE, Scaling = FALSE"))

    data<- Input_data
    PCA<- autoplot(prcomp(as.matrix(data),scale. = FALSE),   # Run and plot PCA
                   data= data,
                   size = 3,
                   alpha = 0.8,
                   label=T,
                   label.size=2.5,
                   label.repel = TRUE) +
      ggtitle(paste(OutputPlotName)) +
      Theme +
      geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
      geom_vline(xintercept=0,  color = "black", linewidth=0.1)

    loading_data <- prcomp(as.matrix(data),scale. = FALSE)

    # For the selected few
  }
  else {
    print("What are you doing? You've got it all wrong!")
  }

  loading_data_table <-as.data.frame(loading_data$rotation)
  loading_data_table <- loading_data_table[,1:2]
  loading_data_table <- tibble::rownames_to_column(loading_data_table, "Metabolite")

  # Save output
  writexl::write_xlsx(loading_data_table, paste(Results_folder_plots_PCA_folder,"/", OutputPlotName, "_Loadings.xlsx", sep=""), col_names = TRUE)

  if(OutputPlotName ==""){
    ggsave(file=paste(Results_folder_plots_PCA_folder,"/", "PCA", OutputPlotName, ".",Save_as, sep=""), plot=PCA, width=10, height=10)
  }else{
    ggsave(file=paste(Results_folder_plots_PCA_folder,"/", "PCA_", OutputPlotName, ".",Save_as, sep=""), plot=PCA, width=10, height=10)
  }

}


##------- How to use -------##
# PCA(Input_data = preprocessing_output_Intra$Processed_data)
# PCA(Input_data = preprocessing_output_Intra$Processed_data, Color = "Conditions",  Shape = "Symbols_for_PCA")

# Notes
# The x=0 and y=0 black lines in PCA will always be there regardless the theme change. I cannot yet make it to be there by default and not be there when you change theme.
# Well It can be done but it requires a lot additional of work. So for now the lines will be there
# Palette changing is still missing
# To do: select a better palette and add option to the user to change the palette to whatever they like




#####################################
### ### ### Volcano Plots ### ### ###
#####################################
#' @param Plot_Settings \emph{Optional: } Choose between "Standard" (Input_data), "Compare" (plot two comparisons together Input_data and Input_data2) or "GSE" (Gene Set Enrichment results) \strong{Default = "Standard"}
#' @param Plot_SettingsInfo \emph{Optional: } Named vector including at least one of those three information: c(color="File$ColumnName", shape= "ColumnName", individual="ColumnName") \strong{Default = NULL}
#' @param Plot_SettingsFile \emph{Optional: } DF with column "Metabolite" including the Metabolite names (needs to match Metabolite names of Input_data) and other columns with required PlotSettingInfo. \strong{Default = NULL}
#' @param Input_data DF with column "Metabolite" including the Metabolite names, Log2FC, pvalue/padjusted values. Can also include additional columns with metadata usable for Plot_Setting_Info.
#' @param y \emph{Optional: } Column name including the values that should be used for y-axis. Usually this would include the p.adjusted value. \strong{Default = "p.adj"}
#' @param x \emph{Optional: } Column name including the values that should be used for x-axis. Usually this would include the Log2FC value. \strong{Default = "Log2FC"}
#' @param AdditionalInput_data \emph{Optional: } DF to compare to main Input_data with the same column names x and y (Plot_Settings="Compare") or Gene set enrichment analysis results (Plot_Settings="GSE"). \strong{Default = NULL}
#' @param Output_Name \emph{Optional: } String which is added to the output files of the plot. \strong{Default = ""}
#'
#' @param Comparison_name \emph{Optional: } amed vector including those information about the two datasets that are compared on the plots when choosing Plot_Settings= "Compare". \strong{Default = c(Input_data="Cond1", AdditionalInput_data= "Cond2")}
#' @param xlab \emph{Optional: } String to replace x-axis label in plot. \strong{Default = NULL}
#' @param ylab \emph{Optional: } String to replace y-axis label in plot. \strong{Default = NULL}
#' @param pCutoff \emph{Optional: } Number of the desired p value cutoff for assessing significance. \strong{Default = 0.05}
#' @param FCcutoff \emph{Optional: } Number of the desired log fold change cutoff for assessing significance. \strong{Default = 0.5}
#' @param Connectors \emph{Optional: } TRUE or FALSE for whether Connectors from names to points are to be added to the plot. \strong{Default =  FALSE}
#' @param Subtitle \emph{Optional: } \strong{Default = NULL}
#' @param Theme \emph{Optional: } Selection of theme for plot. \strong{Default = theme_classic}
#'
#' @param Save_as \emph{Optional: } Select the file type of output plots. Options are svg or pdf. \strong{Default = "svg"}
#'
#' @keywords Volcano plot, pathways
#' @export

# Helper function needed for adding column to pathway file defining if this metabolite is unique/multiple pathways
VizVolcano <- function(Plot_Settings="Standard",
                       Plot_SettingsInfo= NULL,
                       Plot_SettingsFile= NULL,
                       Input_data,
                       y= "p.adj",
                       x= "Log2FC",
                       AdditionalInput_data= NULL,
                       OutputPlotName= "",
                       Comparison_name= c(Input_data="Cond1", AdditionalInput_data= "Cond2"),
                       xlab= NULL,
                       ylab= NULL,
                       pCutoff= 0.05,
                       FCcutoff= 0.5,
                       Connectors=  FALSE,
                       Subtitle= "",
                       Theme= NULL,
                       Save_as= "svg"
                       #add assign=TRUE/FALSE if the user wants the plotlist returned
                       #add parameter for margins that the plot should be saved in
                       ){

  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse", "EnhancedVolcano")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  suppressMessages(library(tidyverse))

  ## ------------ Check Input files ----------- ##
  # 1. The input data:
  if(any(duplicated(Input_data$Metabolite))==TRUE){
    stop("Duplicated Metabolites in Input_data.")#remove duplications, but tell user that this has hapened and should have happened before
  }
  if( is.numeric(pCutoff)== FALSE |pCutoff > 1 | pCutoff < 0){
      stop("Check input. The selected pCutoff value should be numeric and between 0 and 1.")
    }
  if( is.numeric(FCcutoff)== FALSE  | FCcutoff < 0){
      stop("Check input. The selected pCutoff value should be numeric and between 0 and +oo.")
    }
  if(paste(x) %in% colnames(Input_data)==TRUE & paste(y) %in% colnames(Input_data)==TRUE){
      Input1 <- as.data.frame(Input_data)%>%
        dplyr::rename("Log2FC"=paste(x),
                      "p.adj"=paste(y))
    } else{
      stop("Check your input. The column name of x and/ore y does not exist in Input_data.")
    }

    # 2. The Plot_settings
    Plot_options <- c("Standard", "Compare", "GSE")
    if (Plot_Settings %in% Plot_options == FALSE){
      stop("Plot_Settings option is incorrect. The allowed options are the following: ",paste(Plot_options, collapse = ", "),"." )
    }
    if(is.vector(Plot_SettingsInfo)==TRUE & is.null(Plot_SettingsFile)==TRUE){
    stop("You have chosen Plot_SettingsInfo option that requires you to provide a DF Plot_SettingsFile.")
      }
    if(is.vector(Plot_SettingsInfo)==TRUE){
      if("color" %in% names(Plot_SettingsInfo)==TRUE & "shape" %in% names(Plot_SettingsInfo)==TRUE){
        if((Plot_SettingsInfo[["shape"]] == Plot_SettingsInfo[["color"]])==TRUE){
          Plot_SettingsFile$shape <- Plot_SettingsFile[,paste(Plot_SettingsInfo[["color"]])]
          Plot_SettingsFile<- Plot_SettingsFile%>%
            dplyr::rename("color"=paste(Plot_SettingsInfo[["color"]]))
          }
        if((Plot_SettingsInfo[["shape"]] == Plot_SettingsInfo[["color"]])==FALSE & "color" %in% names(Plot_SettingsInfo)==TRUE){
            Plot_SettingsFile <- Plot_SettingsFile%>%
              dplyr::rename("color"=paste(Plot_SettingsInfo[["color"]]))
            }
        if((Plot_SettingsInfo[["shape"]] == Plot_SettingsInfo[["color"]])==FALSE & "shape" %in% names(Plot_SettingsInfo)==TRUE){
              Plot_SettingsFile <- Plot_SettingsFile%>%
                dplyr::rename("shape"=paste(Plot_SettingsInfo[["shape"]]))
            }
      } else if("color" %in% names(Plot_SettingsInfo)==TRUE & "shape" %in% names(Plot_SettingsInfo)==FALSE){
        Plot_SettingsFile <- Plot_SettingsFile%>%
          dplyr::rename("color"=paste(Plot_SettingsInfo[["color"]]))
      } else if("color" %in% names(Plot_SettingsInfo)==FALSE & "shape" %in% names(Plot_SettingsInfo)==TRUE){
        Plot_SettingsFile <- Plot_SettingsFile%>%
          dplyr::rename("shape"=paste(Plot_SettingsInfo[["shape"]]))
      }

      if("individual" %in% names(Plot_SettingsInfo)==TRUE){
      Plot_SettingsFile <- Plot_SettingsFile%>%
        dplyr::rename("individual"=paste(Plot_SettingsInfo[["individual"]]))
      }

    }else if(is.null(Plot_SettingsInfo)==TRUE){
      #No vector provided
    }else{
      stop("Plot_SettingsInfo must be named vector or NULL.")
    }


 #Check that Theme choosen exists for ggplot!

   #3. AdditionalInput_data
  if(Plot_Settings=="Comparison" & is.data.frame(AdditionalInput_data)==TRUE){
    if(paste(x) %in% colnames(Input_data2)==TRUE & paste(y) %in% colnames(Input_data2)==TRUE){
      Input2 <- Plot_SettingsInfo%>%
        dplyr::rename("Log2FC"=paste(x),
                      "p.adj"=paste(y))%>%
        rownames_to_column("Metabolite")
    } else{
      stop("Check your input. The column name of x and/ore y does not exist in Plot_SettingsInfo.")
    }
    #check that colnames match Input_data (or at least partly match)
    #merge the two DFs keeping row.names of Input_data and also adding the names for cond1 and cond2

    if(Plot_Settings=="Individual_Comparison" & is.data.frame(Input_data2)==TRUE){
      #check additionally that individual matches rownames or that column passed by individual exists in Input_data
      #rename row.names if col name passed by individual!
    }

  } #else if(Plot_Settings=="Comparison" & Input_data2==NULL | Plot_Settings=="Individual_Comparison" & Input_data2==NULL){
  # stop()
  #}

  # 4. Check other parameters:
  #save.as
  #outputplotname
  #if (length(unique(Input_pathways$Pathway)) >  17)){#adapt the message
  #    stop(" The maximum number of pathways in the Input_pathways must be less than ",length(safe_colorblind_palette),". Please summarize sub-pathways together where possible and repeat.")
  #}
#test also lth shape!


  if(is.logical(Connectors) == FALSE){
      stop("Check input. The Connectors value should be either = TRUE if connectors from names to points are to be added to the plot or =FALSE if not.")
    }
    Save_as_options <- c("svg","pdf")
    if(Save_as %in% Save_as_options == FALSE){
      stop("Check input. The selected Save_as option is not valid. Please select one of the following: ",paste(Save_as_options,collapse = ", "),"." )
    }

  ## ------------ Create Output folders ----------- ##
  name <- paste0("MetaProViz_Results_",Sys.Date())
  WorkD <- getwd()
  Results_folder <- file.path(WorkD, name)
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)} # Make Results folder
  Results_folder_plots_Volcano_folder = file.path(Results_folder, "Volcano")  # This searches for a folder called "Preprocessing" within the "Results" folder in the current working directory and if its not found it creates one
  if (!dir.exists(Results_folder_plots_Volcano_folder)) {dir.create(Results_folder_plots_Volcano_folder)}  # check and create folder

  ############################################################################################################
  ## ------------ Prepare setting/data for plot ----------- ##

  # Rename the x and y lab if the information has been passed:
  xlab=bquote(~Log[2]~FC)
  ylab=bquote(~-Log[10]~p.adj)


  ############################################################################################################
  ## ----------- Make the  plot based on the choosen parameters ------------ ##
  #Check the plot type: Comparison/Standard/GSE

  #####--- 1. Standard
  if(Plot_Settings=="Standard"){
    if("individual" %in% names(Plot_SettingsInfo)==TRUE){
      # Create the list of individual plots that should be made:
      IndividualPlots <- Plot_SettingsFile[!duplicated(Plot_SettingsFile$individual),]
      IndividualPlots <- IndividualPlots$individual

      PlotList <- list()#Empty list to store all the plots

      for (i in IndividualPlots){
        Plot_SettingsFile_Select <- subset(Plot_SettingsFile, individual == paste(i))
        InputVolcano  <- merge(x=Plot_SettingsFile_Select,y=Input_data, by="Metabolite", all.x=TRUE)%>%
        na.omit()

        if(nrow(InputVolcano)>=1){
          #Prepare the colour scheme:
          if("color" %in% names(Plot_SettingsInfo)==TRUE){
            safe_colorblind_palette <- c("#88CCEE",  "#DDCC77","#661100",  "#332288", "#AA4499","#999933",  "#44AA99", "#882215",  "#6699CC", "#117733", "#888888","#CC6677", "black","gold1","darkorchid4","red","orange")
            color_select <- safe_colorblind_palette[1:length(unique(InputVolcano$color))]

            keyvals <- c()
            for(row in 1:nrow(InputVolcano)){
              col <- color_select[unique(InputVolcano$color) %in% InputVolcano[row, "color"]]
              names(col) <- InputVolcano$color[row]
              keyvals <- c(keyvals, col)
            }
          } else{
            keyvals <-NULL
          }
          #Prepare the shape scheme:
          if("shape" %in% names(Plot_SettingsInfo)==TRUE){
            safe_shape_palette <- c(22,24,21,23,25,7,8,11,12)
            shape_select <- safe_shape_palette[1:length(unique(InputVolcano$shape))]

            keyvalsshape <- c()
            for(row in 1:nrow(InputVolcano)){
              sha <- shape_select[unique(InputVolcano$shape) %in% InputVolcano[row, "shape"]]
              names(sha) <- InputVolcano$shape[row]
              keyvalsshape <- c(keyvalsshape, sha)
            }
          } else{
            keyvalsshape <-NULL
          }

          #Prepare the Plot:
          Plot<- EnhancedVolcano::EnhancedVolcano(InputVolcano,
                                                  lab = InputVolcano$Metabolite,#Metabolite name
                                                  x = paste(x),
                                                  y = paste(y),
                                                  xlab  =xlab,
                                                  ylab =ylab,
                                                  pCutoff = pCutoff,
                                                  FCcutoff = FCcutoff,#Cut off Log2FC, automatically 2
                                                  pointSize = 3,
                                                  labSize = 3,
                                                  titleLabSize = 16,
                                                  colCustom = keyvals,
                                                  shapeCustom = keyvalsshape,
                                                  colAlpha = 1,
                                                  title= paste(OutputPlotName),
                                                  subtitle = i,
                                                  caption = paste0("total = ", nrow(InputVolcano), " Metabolites"),
                                                  xlim =  c(min(Input_data$Log2FC[is.finite(InputVolcano$Log2FC )])-0.2,max(InputVolcano$Log2FC[is.finite(InputVolcano$Log2FC )])+0.2),
                                                  ylim = c(0,(ceiling(-log10(Reduce(min,InputVolcano$p.adj))))),
                                                  cutoffLineType = "dashed",
                                                  cutoffLineCol = "black",
                                                  cutoffLineWidth = 0.5,
                                                  legendLabels=c(paste(x," < |", FCcutoff, "|"), paste(x," > |", FCcutoff, "|"), paste(y, ' < ', pCutoff) , paste(y, ' < ', pCutoff,' & ',x," < |", FCcutoff, "|")),
                                                  legendPosition = 'right',
                                                  legendLabSize = 7,
                                                  legendIconSize =4,
                                                  gridlines.major = FALSE,
                                                  gridlines.minor = FALSE,
                                                  drawConnectors = Connectors)
          #Add the theme
          if(is.null(Theme)==FALSE){
             Plot <- Plot+Theme
          }

          #save plot and get rid of extra signs before saving i
          cleaned_i <- gsub("[[:space:],/\\\\]", "-", i)#removes empty spaces and replaces /,\ with -

          if(OutputPlotName ==""){
            ggsave(file=paste(Results_folder_plots_Volcano_folder,"/", "Volcano_",cleaned_i, ".",Save_as, sep=""), plot=Plot, width=8, height=6)
          }else{
            ggsave(file=paste(Results_folder_plots_Volcano_folder,"/", "Volcano_", OutputPlotName, "_",cleaned_i, ".",Save_as, sep=""), plot=Plot, width=8, height=6)
          }

          ## Store the plot in the 'plots' list
          PlotList[[cleaned_i]] <- Plot
          plot(Plot)
        }
      }
      # Return PlotList into the environment to enable the user to view the plots directly
      assign("VolcanoPlots", PlotList, envir=.GlobalEnv)
      # Combine plots into a single plot using facet_grid or patchwork::wrap_plots

      } else if("individual" %in% names(Plot_SettingsInfo)==FALSE){
        if(is.null(Plot_SettingsFile)==FALSE){
          InputVolcano  <- merge(x=Plot_SettingsFile,y=Input_data, by="Metabolite", all.x=TRUE)%>%
            na.omit()
          }else{
            InputVolcano  <- Input_data
            }

        if(nrow(InputVolcano)>=1){
          #Prepare the colour scheme:
          if("color" %in% names(Plot_SettingsInfo)==TRUE){
            safe_colorblind_palette <- c("#88CCEE",  "#DDCC77","#661100",  "#332288", "#AA4499","#999933",  "#44AA99", "#882215",  "#6699CC", "#117733", "#888888","#CC6677", "black","gold1","darkorchid4","red","orange")
            color_select <- safe_colorblind_palette[1:length(unique(InputVolcano$color))]

            keyvals <- c()
            for(row in 1:nrow(InputVolcano)){
              col <- color_select[unique(InputVolcano$color) %in% InputVolcano[row, "color"]]
              names(col) <- InputVolcano$color[row]
              keyvals <- c(keyvals, col)
            }
          } else{
            keyvals <-NULL
          }
          #Prepare the shape scheme:
          if("shape" %in% names(Plot_SettingsInfo)==TRUE){
            safe_shape_palette <- c(22,24,21,23,25,7,8,11,12)
            shape_select <- safe_shape_palette[1:length(unique(InputVolcano$shape))]

            keyvalsshape <- c()
            for(row in 1:nrow(InputVolcano)){
              sha <- shape_select[unique(InputVolcano$shape) %in% InputVolcano[row, "shape"]]
              names(sha) <- InputVolcano$shape[row]
              keyvalsshape <- c(keyvalsshape, sha)
            }
          } else{
            keyvalsshape <-NULL
          }

          #Prepare the Plot:
          Plot<- EnhancedVolcano::EnhancedVolcano(InputVolcano,
                                                  lab = InputVolcano$Metabolite,#Metabolite name
                                                  x = paste(x),
                                                  y = paste(y),
                                                  xlab  =xlab,
                                                  ylab =ylab,
                                                  pCutoff = pCutoff,
                                                  FCcutoff = FCcutoff,#Cut off Log2FC, automatically 2
                                                  pointSize = 3,
                                                  labSize = 3,
                                                  titleLabSize = 16,
                                                  colCustom = keyvals,
                                                  shapeCustom = keyvalsshape,
                                                  colAlpha = 1,
                                                  title= paste(OutputPlotName),
                                                  subtitle = Subtitle,
                                                  caption = paste0("total = ", nrow(InputVolcano), " Metabolites"),
                                                  xlim =  c(min(Input_data$Log2FC[is.finite(InputVolcano$Log2FC )])-0.2,max(InputVolcano$Log2FC[is.finite(InputVolcano$Log2FC )])+0.2),
                                                  ylim = c(0,(ceiling(-log10(Reduce(min,InputVolcano$p.adj))))),
                                                  cutoffLineType = "dashed",
                                                  cutoffLineCol = "black",
                                                  cutoffLineWidth = 0.5,
                                                  legendLabels=c(paste(x," < |", FCcutoff, "|"), paste(x," > |", FCcutoff, "|"), paste(y, ' < ', pCutoff) , paste(y, ' < ', pCutoff,' & ',x," < |", FCcutoff, "|")),
                                                  legendPosition = 'right',
                                                  legendLabSize = 7,
                                                  legendIconSize =4,
                                                  gridlines.major = FALSE,
                                                  gridlines.minor = FALSE,
                                                  drawConnectors = Connectors)
          #Add the theme
          if(is.null(Theme)==FALSE){
            Plot <- Plot+Theme
          }

          #save plot and get rid of extra signs before saving i
          if(OutputPlotName ==""){
            ggsave(file=paste(Results_folder_plots_Volcano_folder,"/", "Volcano." ,Save_as, sep=""), plot=Plot, width=8, height=6)
          }else{
            ggsave(file=paste(Results_folder_plots_Volcano_folder,"/", "Volcano_", OutputPlotName, ".",Save_as, sep=""), plot=Plot, width=8, height=6)
          }

          #Plot
          plot(Plot)
        }
      }
  } else if(Plot_Settings=="Condition"){

  } else if(Plot_Settings=="GSE"){

  }
}




######################################
### ### ### Alluvial Plots ### ### ###
######################################=
#' Notes; Check the alluvial in the Metabolic Clusters. This one is a bit old but kept here to return to if nedded

VizAlluvial <- function(Input_data1,
                     Input_data2,
                     Output_Name = "Metabolic_Clusters_Output_Condition1-versus-Condition2",
                     Condition1,
                     Condition2,
                     pCutoff= 0.05 ,
                     FCcutoff=0.5,
                     test = "p.adj",
                     plot_column_names= c("class", "MetaboliteChange_Significant", "Overall_Change", "Metabolite"),
                     safe_colorblind_palette = c("#88CCEE",  "#DDCC77","#661100",  "#332288", "#AA4499","#999933",  "#44AA99", "#882255",  "#6699CC", "#117733", "#888888","#CC6677", "#FFF", "#000"), # https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
                     plot_color_variable = "Overall_Change",
                     plot_color_remove_variable = "SameDirection_NoChange",
                     Save_as = pdf
                     ){


  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse", "alluvial")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0){
    install.packages(new.packages)
  }
  suppressMessages(library(tidyverse))

  ####################################################
  # This searches for a Results directory within the current working directory and if its not found it creates a new one
  Results_folder = paste(getwd(), "/MetaProViz_Results_",Sys.Date(),  sep="")
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)}
  ### Create Volcano plots folder in  result directory ###
  Results_folder_plots_MetabolicCluster_folder = paste(Results_folder,"/Alluvial",  sep="")
  if (!dir.exists(Results_folder_plots_MetabolicCluster_folder)) {dir.create(Results_folder_plots_MetabolicCluster_folder)}



  #####################################################
  ### ### ### make output plot save_as name ### ### ###
  Save_as_var <- Save_as
  Save_as= deparse(substitute(Save_as))


  C1 <- Input_data1
  C1 <- na.omit(C1)
  C1$class <- paste (Condition1)
  C2 <- Input_data2
  C2 <- na.omit(C2)
  C2$class <- paste (Condition2)

  # 0. Overall Regulation: label the selection of metabolites that change or where at least we have a change in one of the two conditions
  if(test== "p.val"){
    C1 <- C1 %>%
      mutate(MetaboliteChange_Significant = case_when(Log2FC >= FCcutoff & p.val < pCutoff ~ 'UP',
                                                      Log2FC <= -FCcutoff & p.val < pCutoff ~ 'DOWN',
                                                      TRUE ~ 'No_Change'))
    C2 <- C2%>%
      mutate(MetaboliteChange_Significant1 = case_when(Log2FC >= FCcutoff & p.val < pCutoff ~ 'UP',
                                                       Log2FC <= -FCcutoff & p.val < pCutoff ~ 'DOWN',
                                                       TRUE ~ 'No_Change'))%>%
      rename("class1"="class",
             "Log2FC1"="Log2FC",
             "p.val1"="p.val",
             "p.adj1"="p.adj")

  }else{
    C1 <- C1 %>%
      mutate(MetaboliteChange_Significant = case_when(Log2FC >= FCcutoff & p.adj < pCutoff ~ 'UP',
                                                      Log2FC <= -FCcutoff & p.adj < pCutoff ~ 'DOWN',
                                                      TRUE ~ 'No_Change'))
    C2 <- C2%>%
      mutate(MetaboliteChange_Significant1 = case_when(Log2FC >= FCcutoff & p.adj < pCutoff ~ 'UP',
                                                       Log2FC <= -FCcutoff & p.adj < pCutoff ~ 'DOWN',
                                                       TRUE ~ 'No_Change'))%>%
      rename("class1"="class",
             "Log2FC1"="Log2FC",
             "p.val1"="p.val",
             "p.adj1"="p.adj")
  }



  MergeDF<- merge(C1, C2[,c("Metabolite","MetaboliteChange_Significant1", "class1","Log2FC1","p.val1","p.adj1")], by="Metabolite")%>%
    mutate(Change_Specific = case_when((class== paste(Condition1)& MetaboliteChange_Significant == "UP")       & (class1== paste(Condition2)& MetaboliteChange_Significant1 == "DOWN") ~ 'OppositeChange',
                                       (class== paste(Condition1)& MetaboliteChange_Significant == "DOWN")     & (class1== paste(Condition2)& MetaboliteChange_Significant1 == "UP") ~ 'OppositeChange',
                                       (class== paste(Condition1)& MetaboliteChange_Significant == "UP")       & (class1== paste(Condition2)& MetaboliteChange_Significant1 == "UP") ~ 'SameDirection_UP',
                                       (class== paste(Condition1)& MetaboliteChange_Significant == "DOWN")     & (class1== paste(Condition2)& MetaboliteChange_Significant1 == "DOWN") ~ 'SameDirection_DOWN',
                                       (class== paste(Condition1)& MetaboliteChange_Significant == "No_Change")& (class1== paste(Condition2)& MetaboliteChange_Significant1 == "DOWN") ~ paste("ChangeOnly", Condition2, "DOWN", sep="_"),
                                       (class== paste(Condition1)& MetaboliteChange_Significant == "No_Change")& (class1== paste(Condition2)& MetaboliteChange_Significant1 == "UP") ~ paste("ChangeOnly", Condition2, "UP", sep="_"),
                                       (class== paste(Condition1)& MetaboliteChange_Significant == "UP")       & (class1== paste(Condition2)& MetaboliteChange_Significant1 == "No_Change") ~ paste("ChangeOnly", Condition1, "UP", sep="_"),
                                       (class== paste(Condition1)& MetaboliteChange_Significant == "DOWN")     & (class1== paste(Condition2)& MetaboliteChange_Significant1 == "No_Change") ~ paste("ChangeOnly", Condition1, "DOWN", sep="_"),
                                       TRUE ~ 'SameDirection_NoChange'))


  MergeDF_C1 <-MergeDF %>% select(-c( "MetaboliteChange_Significant1", "class1", "Log2FC1", "p.val1", "p.adj1")) %>%
    unite(col=UniqueID, c(Metabolite, MetaboliteChange_Significant, class), sep = "_", remove = FALSE, na.rm = FALSE)%>%
    mutate(Change_Specific = case_when(Change_Specific== 'OppositeChange' ~ 'OppositeChange',
                                       Change_Specific== 'SameDirection_UP'~ 'SameDirection_UP',
                                       Change_Specific== 'SameDirection_DOWN' ~ 'SameDirection_DOWN',
                                       Change_Specific== 'SameDirection_NoChange' ~ 'SameDirection_NoChange',
                                       Change_Specific== paste("ChangeOnly", Condition2, "DOWN", sep="_") ~paste("ChangeOnly", Condition2, "DOWN", sep="_"),
                                       Change_Specific== paste("ChangeOnly", Condition2, "UP", sep="_") ~paste("ChangeOnly", Condition2, "UP", sep="_"),
                                       Change_Specific== paste("ChangeOnly", Condition1, "DOWN", sep="_") ~paste("Unique", Condition1,"DOWN", sep="_"),
                                       Change_Specific== paste("ChangeOnly", Condition1, "UP", sep="_") ~paste("Unique", Condition1, "UP", sep="_"),
                                       TRUE ~ paste("Unique", Condition1, sep="_")))

  MergeDF_C2 <- MergeDF %>% select(-c( "MetaboliteChange_Significant", "class", "Log2FC", "p.val", "p.adj")) %>%
    unite(col=UniqueID, c(Metabolite, MetaboliteChange_Significant1, class1), sep = "_", remove = FALSE, na.rm = FALSE)%>%
    rename("class"="class1",
           "MetaboliteChange_Significant"="MetaboliteChange_Significant1",
           "Log2FC"="Log2FC1",
           "p.val"="p.val1",
           "p.adj"="p.adj1")%>%
    mutate(Change_Specific = case_when(Change_Specific== 'OppositeChange' ~ 'OppositeChange',
                                       Change_Specific== 'SameDirection_UP'~ 'SameDirection_UP',
                                       Change_Specific== 'SameDirection_DOWN' ~ 'SameDirection_DOWN',
                                       Change_Specific== 'SameDirection_NoChange' ~ 'SameDirection_NoChange',
                                       Change_Specific== paste("ChangeOnly", Condition1, "DOWN", sep="_") ~paste("ChangeOnly", Condition1, "DOWN", sep="_"),
                                       Change_Specific== paste("ChangeOnly", Condition1, "UP", sep="_") ~paste("ChangeOnly", Condition1, "UP", sep="_"),
                                       Change_Specific== paste("ChangeOnly", Condition2, "DOWN", sep="_") ~paste("Unique", Condition2,"DOWN", sep="_"),
                                       Change_Specific== paste("ChangeOnly", Condition2, "UP", sep="_") ~paste("Unique", Condition2, "UP", sep="_"),
                                       TRUE ~ paste("Unique", Condition1, sep="_")))


  Alluvial_DF <- rbind(MergeDF_C1, MergeDF_C2)
  Alluvial_DF<- Alluvial_DF%>%
    mutate(Amount_Change_Specific = case_when(Change_Specific== 'OppositeChange' ~ paste((sum(Alluvial_DF$Change_Specific=="OppositeChange", na.rm=T))/2),
                                              Change_Specific== 'SameDirection_UP' ~ paste((sum(Alluvial_DF$Change_Specific=="SameDirection_UP", na.rm=T))/2),
                                              Change_Specific== 'SameDirection_DOWN' ~ paste((sum(Alluvial_DF$Change_Specific=="SameDirection_DOWN", na.rm=T))/2),
                                              Change_Specific== paste("ChangeOnly", Condition1, "UP", sep="_") ~ paste(sum(Alluvial_DF$Change_Specific==paste("ChangeOnly", Condition1, "UP", sep="_"), na.rm=T)),
                                              Change_Specific== paste("ChangeOnly", Condition1, "DOWN", sep="_") ~ paste(sum(Alluvial_DF$Change_Specific==paste("ChangeOnly", Condition1, "DOWN", sep="_"), na.rm=T)),
                                              Change_Specific== paste("ChangeOnly", Condition2, "UP" ,sep="_") ~ paste(sum(Alluvial_DF$Change_Specific==paste("ChangeOnly", Condition2, "UP", sep="_"), na.rm=T)),
                                              Change_Specific== paste("ChangeOnly", Condition2, "DOWN", sep="_") ~ paste(sum(Alluvial_DF$Change_Specific==paste("ChangeOnly", Condition2, "DOWN", sep="_"), na.rm=T)),
                                              Change_Specific== 'SameDirection_NoChange' ~ paste((sum(Alluvial_DF$Change_Specific=="SameDirection_NoChange", na.rm=T))/2),
                                              Change_Specific== paste("Unique", Condition1,"DOWN", sep="_") ~paste(sum(Alluvial_DF$Change_Specific==paste("Unique", Condition1,"DOWN", sep="_"), na.rm=T)),
                                              Change_Specific== paste("Unique", Condition1,"UP", sep="_") ~paste(sum(Alluvial_DF$Change_Specific==paste("Unique", Condition1,"UP", sep="_"), na.rm=T)),
                                              Change_Specific== paste("Unique", Condition2,"DOWN", sep="_") ~paste(sum(Alluvial_DF$Change_Specific==paste("Unique", Condition2,"DOWN", sep="_"), na.rm=T)),
                                              Change_Specific== paste("Unique", Condition2,"UP", sep="_") ~paste(sum(Alluvial_DF$Change_Specific==paste("Unique", Condition2,"UP", sep="_"), na.rm=T)),
                                              TRUE ~ 'FALSE'))
  Alluvial_DF<- Alluvial_DF%>%
    mutate(Overall_Change = case_when(Change_Specific== 'OppositeChange' ~ 'OppositeChange',
                                      Change_Specific== 'SameDirection_UP'~ 'SameDirection_UP_or_DOWN',
                                      Change_Specific== 'SameDirection_DOWN' ~ 'SameDirection_UP_or_DOWN',
                                      Change_Specific== 'SameDirection_NoChange' ~ 'SameDirection_NoChange',
                                      Change_Specific== paste("ChangeOnly", Condition2, "DOWN", sep="_") & MetaboliteChange_Significant == "No_Change" ~ paste("ChangeOnly", Condition2, sep="_"),
                                      Change_Specific== paste("ChangeOnly", Condition2, "UP", sep="_") & MetaboliteChange_Significant == "No_Change" ~ paste("ChangeOnly", Condition2, sep="_"),
                                      Change_Specific== paste("Unique", Condition1, "DOWN",sep="_") ~paste("Unique", Condition1, sep="_"),
                                      Change_Specific== paste("Unique", Condition1, "UP",sep="_") ~paste("Unique", Condition1, sep="_"),
                                      Change_Specific== paste("ChangeOnly", Condition1, "DOWN", sep="_") & MetaboliteChange_Significant == "No_Change" ~ paste("ChangeOnly", Condition1, sep="_"),
                                      Change_Specific== paste("ChangeOnly", Condition1, "UP", sep="_") & MetaboliteChange_Significant == "No_Change" ~ paste("ChangeOnly", Condition1, sep="_"),
                                      Change_Specific== paste("Unique", Condition2,"DOWN", sep="_") ~paste("Unique", Condition2, sep="_"),
                                      Change_Specific== paste("Unique", Condition2,"UP", sep="_") ~paste("Unique", Condition2, sep="_"),
                                      TRUE ~ "FALSE"))
  Alluvial_DF <- Alluvial_DF %>%
    mutate(Amount_Overall_Change = case_when(Overall_Change== 'OppositeChange' ~ paste((sum(Alluvial_DF$Overall_Change=="OppositeChange", na.rm=T))/2),
                                             Overall_Change== 'SameDirection_UP_or_DOWN' ~ paste((sum(Alluvial_DF$Overall_Change=="SameDirection_UP_or_DOWN", na.rm=T))/2),
                                             Overall_Change== paste("Unique", Condition1, sep="_") ~ paste(sum(Alluvial_DF$Overall_Change==paste("Unique", Condition1, sep="_"), na.rm=T)),
                                             Overall_Change== paste("Unique", Condition2, sep="_") ~ paste(sum(Alluvial_DF$Overall_Change==paste("Unique", Condition2, sep="_"), na.rm=T)),
                                             Overall_Change== paste("ChangeOnly", Condition1, sep="_") ~ paste(sum(Alluvial_DF$Overall_Change==paste("ChangeOnly", Condition1, sep="_"), na.rm=T)),
                                             Overall_Change== paste("ChangeOnly", Condition2, sep="_") ~ paste(sum(Alluvial_DF$Overall_Change==paste("ChangeOnly", Condition2, sep="_"), na.rm=T)),
                                             Overall_Change== 'SameDirection_NoChange' ~ paste((sum(Alluvial_DF$Overall_Change=="SameDirection_NoChange", na.rm=T))/2),
                                             TRUE ~ 'FALSE'))

  # Create Alluvial final output df
  C1.final<- Alluvial_DF %>% filter(class == Condition1)
  C1.final <- C1.final %>% select(-c("UniqueID", "class", "Change_Specific","Amount_Change_Specific", "Overall_Change", "Amount_Overall_Change" ))
  C2.final<- Alluvial_DF %>% filter(class == Condition2)
  C2.final <- C2.final %>% select(-c("UniqueID", "class" ))
  Alluvial_DF.final <- merge(C1.final, C2.final, by = "Metabolite")
  names(Alluvial_DF.final) <- gsub(".x",paste(".",substr(Condition1, 1, 3), sep = ""),names(Alluvial_DF.final))
  names(Alluvial_DF.final) <- gsub(".y",paste(".",substr(Condition2, 1, 3), sep = ""),names(Alluvial_DF.final))
  names(Alluvial_DF.final) <- gsub(x = names(Alluvial_DF.final), pattern = "MetaboliteChange_Significant", replacement =  paste("MetaboliteChange_Significant_",test,pCutoff,"logFC",FCcutoff, sep = ""))

  ##Write to file
  # This is not needed fot the plots
  #  writexl::write_xlsx(Alluvial_DF.final, paste("Results_", Sys.Date(), "/MetabolicCluster_plots/","Metabolic_Clusters_Output_",Condition1,"-versus-",Condition2,Output_Name,".xlsx", sep = ""))
  # write.csv(Alluvial_DF2, paste("AlluvianDF", Output, ".csv", sep="_"), row.names= TRUE)


  # 1. Regulation:
  Alluvial_DF2  <- Alluvial_DF  %>%
    mutate(MetaboliteChange = case_when(Log2FC >= FCcutoff  ~ 'UP',
                                        Log2FC <= -FCcutoff ~ 'DOWN',
                                        TRUE ~ 'No_Change'))

  if (test == "p.val"){
    # 2. Excluded according to p-value:
    Alluvial_DF2  <- Alluvial_DF2  %>%
      mutate(Excluded_by_pval = case_when(p.val <= pCutoff ~ 'NO',
                                          p.val > pCutoff ~ 'YES'))
    #3. Excluded according to P-value?
    Alluvial_DF2  <- Alluvial_DF2  %>%
      mutate(MetaboliteChange_Excluded = case_when(Log2FC >= FCcutoff & p.val < pCutoff ~ 'UP',
                                                   Log2FC <= -FCcutoff & p.val < pCutoff ~ 'DOWN',
                                                   Log2FC >= FCcutoff & p.val >= pCutoff ~ 'UP_Excluded',
                                                   Log2FC <= -FCcutoff & p.val >= pCutoff ~ 'DOWN_Excluded',
                                                   TRUE ~ 'No_Change'))
    #4. After exlusion by p-value
    Alluvial_DF2  <- Alluvial_DF2  %>%
      mutate(MetaboliteChange_Significant = case_when(Log2FC >= FCcutoff & p.val < pCutoff ~ 'UP',
                                                      Log2FC <= -FCcutoff & p.val < pCutoff ~ 'DOWN',
                                                      TRUE ~ 'No_Change'))
  }else if(test=="p.adj"){
    # 2. Excluded according to p-value:
    Alluvial_DF2  <- Alluvial_DF2  %>%
      mutate(Excluded_by_pval = case_when(p.adj <= pCutoff ~ 'NO',
                                          p.adj > pCutoff ~ 'YES'))
    #3. Excluded according to P-value?
    Alluvial_DF2  <- Alluvial_DF2  %>%
      mutate(MetaboliteChange_Excluded = case_when(Log2FC >= FCcutoff & p.adj < pCutoff ~ 'UP',
                                                   Log2FC <= -FCcutoff & p.adj < pCutoff ~ 'DOWN',
                                                   Log2FC >= FCcutoff & p.adj >= pCutoff ~ 'UP_Excluded',
                                                   Log2FC <= -FCcutoff & p.adj >= pCutoff ~ 'DOWN_Excluded',
                                                   TRUE ~ 'No_Change'))
    #4. After exlusion by p-value
    Alluvial_DF2  <- Alluvial_DF2  %>%
      mutate(MetaboliteChange_Significant = case_when(Log2FC >= FCcutoff & p.adj < pCutoff ~ 'UP',
                                                      Log2FC <= -FCcutoff & p.adj < pCutoff ~ 'DOWN',
                                                      TRUE ~ 'No_Change'))
  }


  #5. Frequency:
  Alluvial_DF2[,"Frequency"]  <- as.numeric("1")
  #6. Safe DF
  # Alluvial_DF3 <- Alluvial_DF2[, c(1:6,12:14,7:8,10,9,11,15)]
  # Alluvial_Plot <- Alluvial_DF3[,c(6:14,2,15)]

  Alluvial_Plot <- Alluvial_DF2

  #Make the SelectionPlot:
  if(plot_color_remove_variable %in% Alluvial_Plot[,plot_color_variable] ){
    Alluvial_Plot<- Alluvial_Plot[-which(Alluvial_Plot[,plot_color_variable]==plot_color_remove_variable),]#remove the metabolites that do not change in either of the conditions
  }


  Save_as_var(paste(Results_folder_plots_MetabolicCluster_folder,"/Metabolic_Clusters_",Condition1,"-versus-",Condition2,"_", Output_Name,  ".",Save_as, sep=""), width=12, height=9)
  par(oma=c(2,2,8,2), mar = c(2, 2, 0.1, 2)+0.1)#https://www.r-graph-gallery.com/74-margin-and-oma-cheatsheet.html
  alluvial::alluvial( Alluvial_Plot %>% select(all_of(plot_column_names)), freq=Alluvial_Plot$Frequency,
            col = case_when(Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[1] ~ safe_colorblind_palette[1],
                            Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[2] ~ safe_colorblind_palette[2],
                            Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[3] ~ safe_colorblind_palette[3],
                            Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[4] ~ safe_colorblind_palette[4],
                            Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[5] ~ safe_colorblind_palette[5],
                            Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[6] ~ safe_colorblind_palette[6],
                            Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[7] ~ safe_colorblind_palette[7],
                            Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[8] ~ safe_colorblind_palette[8],
                            Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[9] ~ safe_colorblind_palette[9],
                            Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[10] ~ safe_colorblind_palette[10],
                            TRUE ~ 'black'),
            border = case_when(Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[1] ~ safe_colorblind_palette[1],
                               Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[2] ~ safe_colorblind_palette[2],
                               Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[3] ~ safe_colorblind_palette[3],
                               Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[4] ~ safe_colorblind_palette[4],
                               Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[5] ~ safe_colorblind_palette[5],
                               Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[6] ~ safe_colorblind_palette[6],
                               Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[7] ~ safe_colorblind_palette[7],
                               Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[8] ~ safe_colorblind_palette[8],
                               Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[9] ~ safe_colorblind_palette[9],
                               Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[10] ~ safe_colorblind_palette[10],

                               TRUE ~ 'black'),
            hide = Alluvial_Plot$Frequency == 0,
            cex = 0.3,
            cex.axis=0.5)
  mtext("Selection of metabolites that change in at least one of the two conditions", side=3, line=6, cex=1.2, col="black", outer=TRUE) #https://www.r-graph-gallery.com/74-margin-and-oma-cheatsheet.html
  mtext(paste("",Output_Name), side=3, line=5, cex=0.8, col="black", outer=TRUE)
  mtext(paste("Underlying comparison: ",Condition1,"-versus-",Condition2), side=2, line=0, cex=0.8, col="black", outer=TRUE)
  #mtext("Legend", side=3, line=5, adj=1.0, cex=1, col="black", outer=TRUE)
  if( is.na(unique(Alluvial_Plot[,plot_color_variable])[1])==FALSE)  {
    mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[1]), side=3, line=6, adj=0, cex=0.6, col=safe_colorblind_palette[1], outer=TRUE)
  }
  if( is.na(unique(Alluvial_Plot[,plot_color_variable])[2])==FALSE)  {
    mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[2]), side=3, line=5, adj=0, cex=0.6, col=safe_colorblind_palette[2], outer=TRUE)
  }
  if( is.na(unique(Alluvial_Plot[,plot_color_variable])[3])==FALSE)  {
    mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[3]), side=3, line=4, adj=0, cex=0.6, col=safe_colorblind_palette[3], outer=TRUE)
  }
  if( is.na(unique(Alluvial_Plot[,plot_color_variable])[4])==FALSE)  {
    mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[4]), side=3, line=3, adj=0, cex=0.6, col=safe_colorblind_palette[4], outer=TRUE)
  }
  if( is.na(unique(Alluvial_Plot[,plot_color_variable])[5])==FALSE)  {
    mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[5]), side=3, line=2, adj=0, cex=0.6, col=safe_colorblind_palette[5], outer=TRUE)
  }
  if( is.na(unique(Alluvial_Plot[,plot_color_variable])[6])==FALSE)  {
    mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[6]), side=3, line=1, adj=0, cex=0.6, col=safe_colorblind_palette[6], outer=TRUE)
  }
  if( is.na(unique(Alluvial_Plot[,plot_color_variable])[7])==FALSE)  {
    mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[7]), side=3, line=0, adj=0, cex=0.6, col=safe_colorblind_palette[7], outer=TRUE)
  }
  if( is.na(unique(Alluvial_Plot[,plot_color_variable])[8])==FALSE)  {
    mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[8]), side=3, line=7, adj=1, cex=0.6, col=safe_colorblind_palette[8], outer=TRUE)
  }
  if( is.na(unique(Alluvial_Plot[,plot_color_variable])[9])==FALSE)  {
    mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[9]), side=3, line=6, adj=1, cex=0.6, col=safe_colorblind_palette[9], outer=TRUE)
  }
  if( is.na(unique(Alluvial_Plot[,plot_color_variable])[10])==FALSE)  {
    mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[10]), side=3, line=5, adj=1, cex=0.6, col=safe_colorblind_palette[10], outer=TRUE)
  }
  if( is.na(unique(Alluvial_Plot[,plot_color_variable])[11])==FALSE)  {
    mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[11]), side=3, line=4, adj=1, cex=0.6, col=safe_colorblind_palette[11], outer=TRUE)
  }
  if( is.na(unique(Alluvial_Plot[,plot_color_variable])[12])==FALSE)  {
    mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[12]), side=3, line=3, adj=1, cex=0.6, col=safe_colorblind_palette[12], outer=TRUE)
  }
  if( is.na(unique(Alluvial_Plot[,plot_color_variable])[13])==FALSE)  {
    mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[13]), side=3, line=2, adj=1, cex=0.6, col=safe_colorblind_palette[13], outer=TRUE)
  }
  if( is.na(unique(Alluvial_Plot[,plot_color_variable])[14])==FALSE)  {
    mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[14]), side=3, line=1, adj=1, cex=0.6, col=safe_colorblind_palette[14], outer=TRUE)
  }
  if( is.na(unique(Alluvial_Plot[,plot_color_variable])[15])==FALSE)  {
    mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[15]), side=3, line=7, adj=1, cex=0.6, col=safe_colorblind_palette[15], outer=TRUE)
  }
  dev.off()# Close the pdf file

}


#####################################
### ### ### Lolipop Plots ### ### ###
#####################################

#' @param Input_data Dataframe which contains metabolites in rows and Log fold changes, pvalues and padjusted values in columns. Multiple Input data frames can be added using a list ie. list(df1, df2, df3)
#' @param test \emph{Optional: } String which selects pvalue or padj for significance. \strong{Default = padj}
#' @param pCutoff \emph{Optional: } Number of the desired p value cutoff for assessing significance. \strong{Default = 0.05}
#' @param FCcutoff \emph{Optional: } Number of the desired log fold change cutoff for assessing significance. \strong{Default = 0.5}
#' @param Output_Name \emph{Optional: } String which is added to the output files of the plot
#' @param Input_pathways \emph{Optional: } DF which contains a 'Metabolite' and a 'Pathway' columns, with pathway information for each metabolite. It can be the same as Input_data or another data frame. \strong{Default = NULL}
#' @param Plot_pathways \emph{Optional: } String with plotting information about Metabolite pathways. Available only when the Input_pathways parameter has a file. Options are "None" if no pathways are to be plotted, "Individual" for plots of each Individual pathway and "Together" for metabolite pathways color-coded on a single volcano plot \strong{Default = "Together"}
#' @param Theme \emph{Optional: } Selection of theme for plot. \strong{Default = theme_classic} ??
#' @param CondNames \emph{Optional: } list of Input dataset names as strings when multiple input datasets are used. \strong{Default = Comparisson 1}
#' @param Comparison \emph{Optional: } String that is placed as the plot title which multiple datasets are used. \strong{Default = Comparisson 1} ## add this to multiple = false
#' @param Connectors \emph{Optional: } TRUE or FALSE for whether Connectors from names to points are to be added to the plot. \strong{Default =  FALSE}
#' @param Save_as \emph{Optional: } Select the file type of output plots. Options are svg or pdf. \strong{Default = svg}
#'
#' @keywords Volcano plot, pathways
#' @export
#'

VizLolipop <- function(Input_data, # a dataframe of list of dataframes
                    test = "p.adj",
                    pCutoff= 0.05 ,
                    FCcutoff=0.5,
                    OutputPlotName= "",
                    Input_pathways = NULL,
                    Plot_pathways = "None",# or "Individual" or "Together
                    Theme=theme_classic(),
                    CondNames = NULL, # list of dataframe names
                    Comparison = "Plot Title",
                    Save_as = "svg"
                    ){


  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse", "showtext", "cowplot")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  suppressMessages(library(tidyverse))

  ## ------------ Check Input files ----------- ##
  if(inherits(Input_data,"list") ==FALSE){
    temp <- list(Input_data)
    Input_data <- temp
  }
  for(data in Input_data){
    if(any(duplicated(row.names(data)))==TRUE){
      stop("Duplicated row.names of Input_data, whilst row.names must be unique")
    }
    if("Metabolite" %in% colnames(data) == FALSE){
      stop("Check Input. Metabolite column is missing from Input_data")
    }
    if(any(duplicated(data$Metabolite))==TRUE){
      stop("Duplicated Metabolite names in Input_data, Metabolites must be unique")
    }
  }

  if( is.numeric(pCutoff)== FALSE |pCutoff > 1 | pCutoff < 0){
    stop("Check input. The selected pCutoff value should be numeric and between 0 and 1.")
  }
  if( is.numeric(FCcutoff)== FALSE  | FCcutoff < 0){
    stop("Check input. The selected pCutoff value should be numeric and between 0 and +oo.")
  }
  if(test != "p.val" & test != "p.adj"){
    stop("Check input. The selected test option for assessing significance is not valid. Please select one of the following: p.adj, p.val.")
  }
  Save_as_options <- c("svg","pdf", "jpeg", "tiff", "png", "bmp", "wmf","eps", "ps", "tex" )
  if(Save_as %in% Save_as_options == FALSE){
    stop("Check input. The selected Save_as option is not valid. Please select one of the folowwing: ",paste(Save_as_options,collapse = ", "),"." )
  }
  Plot_pathways_options <- c("None", "Individual", "Together")
  if (Plot_pathways %in% Plot_pathways_options == FALSE){
    stop("Check Input the Plot_pathways option is incorrect. The Allowed options are the following: ",paste(Plot_pathways_options,collapse = ", "),"." )
  }
  if(is.null(Input_pathways) == TRUE){
    if (Plot_pathways != "None"){
      warning("No Input_pathways were added. Yet the Plot_pathways option was changed. This will have no effect on the plot")
      Plot_pathways = "None"
    }
  }
  if(is.null(Input_pathways) == FALSE){
    if("Metabolite" %in% colnames(Input_pathways) == FALSE){
      stop("Check Input. Metabolite column is missing from Input_pathways")
    }
    if("Pathway" %in% colnames(Input_pathways) == FALSE){
      stop("Check Input. Pathway column is missing from Input_pathways")
    }
    if (sum(duplicated(Input_pathways$Metabolite)) > 0){
      stop("Duplicated Metabolites found in the Input_Pathways. The Metabolites must be unique.")
    }
    for(i in 1:length(Input_data)){ # Here we check every data in the input data with the Input pathways. We dont check common metablites between each input dataset. Should we also add this?
      data <- Input_data[[i]]
      if(identical(sort(data$Metabolite), sort(Input_pathways$Metabolite)) == FALSE){
        warning("The Metabolite column in the Input_data is not the same as the Metabolite column in the Input_pathways. We will take into consideration only the common Metabolites.")
        # find common metabolites
        common_metabolites <- data[data$Metabolite %in% Input_pathways$Metabolite, "Metabolite"]
        # Take the data that have both pval reslults and pathwayss
        Input_data <- data %>% filter(Metabolite %in% common_metabolites)
        Input_pathways <- Input_pathways %>% filter(Metabolite %in% common_metabolites)
      }
    }
  }


  for(i in 1:length(Input_data)){
    data <- Input_data[[i]]
    if("Pathway" %in% colnames(data)){
      data$Pathway <- NULL
      Input_data[[i]] <- data
    }
  }



  ## ------------ Create Output folders ----------- ##
  name <- paste0("MetaProViz_Results_",Sys.Date())
  WorkD <- getwd()
  Results_folder <- file.path(WorkD, name)
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)} # Make Results folder
  Results_folder_plots_Lolipop_folder = file.path(Results_folder, "Lolipop")  # This searches for a folder called "Preprocessing" within the "Results" folder in the current working directory and if its not found it creates one
  if (!dir.exists(Results_folder_plots_Lolipop_folder)) {dir.create(Results_folder_plots_Lolipop_folder)}  # check and create folder


  Multiple = FALSE
  if(length(Input_data) > 1){
    Multiple = TRUE
  }
  if (Plot_pathways == "None"){
    if(is.null(Input_pathways)==FALSE){
      warning("An Input_pathways file has been added. However, Plot_pathways = none ,thus no pathway information will be incorporated in the Volcano plot. To use Pathway information please change Plot_pathways to Individual or Together.")
    }
    Input_pathways = NULL
  }

  if(Multiple == FALSE){
    # Remove rows with NAs
    Input_data <- Input_data[[1]]
    loli.data<- Input_data %>% drop_na()

    if (is.null(Input_pathways)== FALSE){
      Input_pathways <- Input_pathways %>% select(all_of(c("Metabolite", "Pathway")))
      loli.data <- merge(loli.data,Input_pathways, by = "Metabolite" )
    }
    #Select metabolites for the cut offs selected
    loli.data <- loli.data %>% mutate(names=Metabolite) %>% filter( abs(Log2FC) >=FCcutoff & test > pCutoff)

    if(Plot_pathways == "Individual"){

      Pathway_Names <- unique(loli.data$Pathway)
      for (i in Pathway_Names){
        print(i)
        loli.data_path_indi <- loli.data %>% filter(Pathway==i)

        lolipop_plot <- ggplot(loli.data_path_indi , aes(x = Log2FC, y = names)) +
          geom_segment(aes(x = 0, xend = Log2FC, y = names, yend = names)) +
          geom_point(aes(colour = p.adj, size = p.adj ))   +
          scale_size_continuous(range = c(1,5))+# , trans = 'reverse') +
          scale_colour_gradient(low = "red", high = "blue", limits = c(0, max(loli.data[,test]))) +
          ggtitle(label = paste(i," Pathway"), subtitle = paste("Metabolites with > |",FCcutoff,"| logfold change")) + theme(plot.title = element_text(hjust = 0.5)) + ylab("Metabolites")+Theme
        #ggsave(filename = "Loli_plot2.pdf", plot = last_plot(), width=10, height=8)

        if(OutputPlotName ==""){
          ggsave(file=paste(Results_folder_plots_Lolipop_folder, "/",i, ".",Save_as, sep=""), plot=lolipop_plot, width=10, height=10)
        }else{
          ggsave(file=paste(Results_folder_plots_Lolipop_folder,"/",OutputPlotName,"_",i, ".",Save_as, sep=""), plot=lolipop_plot, width=10, height=10)
        }

      }}
    else if(Plot_pathways == "Together"){

      loli.data <- loli.data %>%
        arrange(Pathway, Metabolite)

      loli.data_avg <- loli.data %>%
        arrange(Pathway, Metabolite) %>%
        mutate(Metab_name = row_number()) %>%
        group_by(Pathway) %>%
        mutate(
          avg = mean(Log2FC)
        ) %>%
        ungroup() %>%
        mutate(Pathway = factor(Pathway))


      loli_lines <-   loli.data_avg %>%
        arrange(Pathway, Metabolite) %>%
        group_by(Pathway) %>%
        summarize(
          start_x = min(Metab_name) -0.5,
          end_x = max(Metab_name) + 0.5,
          y = 0#unique(avg)
        ) %>%
        pivot_longer(
          cols = c(start_x, end_x),
          names_to = "type",
          values_to = "x"
        ) %>%
        mutate(
          x_group = if_else(type == "start_x", x + .1, x - .1),
          x_group = if_else(type == "start_x" & x == min(x), x_group - .1, x_group),
          x_group = if_else(type == "end_x" & x == max(x), x_group + .1, x_group) )

      #rm(p2)
      p2 <- loli.data_avg %>%
        ggplot(aes(Metab_name, Log2FC)) + # names in aes ro Metab_name
        geom_hline(
          data = tibble(y = -5:5),
          aes(yintercept = y),
          color = "grey82",
          size = .5 )

      p2 <- p2 + geom_segment(
        aes(
          xend = Metab_name,          # names
          yend = 0,#avg,
          color = Pathway,
          #color = after_scale(colorspace::lighten(color, .2))
        ))

      p2 <- p2 + # geom_line( data = loli_lines, aes(x, y),  color = "grey40"  ) +
        geom_line(
          data = loli_lines,
          aes( x_group, y,
               color = Pathway,
               #  color = after_scale(colorspace::darken(color, .2))
          ), size = 2.5) +  geom_point(aes(size = p.adj, color = Pathway)
          )

      p2<- p2 + coord_flip()
      p2<-p2+ theme(axis.text.x=element_text())

      lab_pos_metab <- loli.data_avg %>% filter(Log2FC>0) %>% select(Metabolite, Metab_name, Log2FC)
      p2<- p2+ annotate("text", x = lab_pos_metab$Metab_name, y = lab_pos_metab$Log2FC+1.5, label = lab_pos_metab$Metabolite, size = 3)

      lab_neg_metab <- loli.data_avg %>% filter(Log2FC<0) %>% select(Metabolite, Metab_name, Log2FC)
      p2<- p2+ annotate("text", x = lab_neg_metab$Metab_name, y = lab_neg_metab$Log2FC-1.5, label = lab_neg_metab$Metabolite, size = 3)

      p2 <- p2+ annotate("text", x = max(lab_neg_metab$Metab_name)+ 3, y = 0, label = "Significantly changed metabolites and their pathways", size = 8)

      p2 <- p2+Theme


      if(OutputPlotName ==""){
        ggsave(file=paste(Results_folder_plots_Lolipop_folder,"/", "Lolipop_Together", OutputPlotName, ".",Save_as, sep=""), plot=p2, width=20, height=20)
      }else{
        ggsave(file=paste(Results_folder_plots_Lolipop_folder,"/", "Lolipop_Together_", OutputPlotName, ".",Save_as, sep=""), plot=p2,  width=20, height=20)
      }

    }
    else if(Plot_pathways == "None"){
      lolipop_plot <- ggplot(loli.data , aes(x = Log2FC, y = names)) +
        geom_segment(aes(x = 0, xend = Log2FC, y = names, yend = names)) +
        geom_point(aes(colour = p.adj, size = p.adj ))   +
        scale_size_continuous(range = c(1,5))+# , trans = 'reverse') +
        scale_colour_gradient(low = "red", high = "blue", limits = c(0, max(loli.data[,test]))) +
        ggtitle(paste("Metabolites with > |",FCcutoff,"| logfold change")) + theme(plot.title = element_text(hjust = 0.5)) + ylab("Metabolites")+Theme

      if(OutputPlotName ==""){
        ggsave(file=paste(Results_folder_plots_Lolipop_folder,"/", "Lolipop", OutputPlotName, ".",Save_as, sep=""), plot=lolipop_plot, width=10, height=10)
      }else{
        ggsave(file=paste(Results_folder_plots_Lolipop_folder,"/", "Lolipop_", OutputPlotName, ".",Save_as, sep=""), plot=lolipop_plot, width=10, height=10)
      }

    }
  }
  else{

    if(Plot_pathways == "Individual"){

      Input_pathways <- Input_pathways %>% select(all_of(c("Metabolite", "Pathway")))
      Combined_Input <- data.frame(matrix(ncol = 5, nrow = 0))
      comb.colnames <- c("Metabolite","Log2FC",test, "Condition","Pathway")
      colnames(Combined_Input) <- comb.colnames

      for ( i in 1:length(Input_data)){
        Input_data[[i]]$Condition <- CondNames[[i]]
        Input_data_pathways <- merge( Input_data[[i]],Input_pathways, by = "Metabolite" )
        Combined_Input <- rbind(Combined_Input, Input_data_pathways %>% select(all_of(c("Metabolite","Log2FC",test, "Condition", "Pathway"))))
      }

      Combined_Input["Pathway"][Combined_Input["Pathway"] == "unknown"] <- "Other"
      Combined_Input[test] <- round(Combined_Input[test], digits = 6)


      for (pathway in unique(Combined_Input$Pathway)){
        Combined_Input_pathway <- Combined_Input %>% filter(Pathway == pathway)

        Dotplot1 <-ggplot(Combined_Input_pathway, aes(x=reorder(Metabolite, + `Log2FC`), y=`Log2FC`, label=`p.adj`)) +
          geom_point(stat='identity', aes(size = abs(`Log2FC`), col=Condition))  +
          geom_segment(aes(y = 0,
                           x = Metabolite,
                           yend = `Log2FC`,
                           xend = Metabolite),
                       color = "black") +
          scale_size(name="abs(Log2FC)",range = c(6,16))+
          geom_text(color="black", size=2) +
          labs(title=paste(Comparison)) +
          ylim(((Reduce(min,Combined_Input_pathway$`Log2FC`))-0.5),((Reduce(max,Combined_Input_pathway$`Log2FC`))+0.5)) +
          theme_minimal() +
          coord_flip()+
          theme(plot.title = element_text(color = "black", size = 12, face = "bold"),
                plot.subtitle = element_text(color = "black", size=10),
                plot.caption = element_text(color = "black",size=9, face = "italic", hjust = 2.5))+
          labs(y="Log2FC", x="")+
          geom_hline(yintercept=0,  color = "black", linewidth=0.1)



        if(OutputPlotName ==""){
          ggsave(file=paste(Results_folder_plots_Lolipop_folder,"/", "Lolipop",pathway, ".",Save_as, sep=""), plot=Dotplot1, width=12, height=14)
        }else{
          ggsave(file=paste(Results_folder_plots_Lolipop_folder,"/", "Lolipop_",pathway,"_", OutputPlotName, ".",Save_as, sep=""), plot=Dotplot1, width=12, height=14)
        }
      }
    }
    else if(Plot_pathways == "Together"){

      Input_pathways <- Input_pathways %>% select(all_of(c("Metabolite", "Pathway")))
      Combined_Input <- data.frame(matrix(ncol = 5, nrow = 0))
      comb.colnames <- c("Metabolite","Log2FC",test, "Condition","Pathway")
      colnames(Combined_Input) <- comb.colnames

      for ( i in 1:length(Input_data)){
        Input_data[[i]]$Condition <- CondNames[[i]]
        Input_data_pathways <- merge( Input_data[[i]],Input_pathways, by = "Metabolite" )
        Combined_Input <- rbind(Combined_Input, Input_data_pathways %>% select(all_of(c("Metabolite","Log2FC",test, "Condition", "Pathway"))))
      }

      Combined_Input["Pathway"][Combined_Input["Pathway"] == "unknown"] <- "Other"
      Combined_Input[test] <- round(Combined_Input[test], digits = 6)


      ######
      loli.data <- Combined_Input

      loli.data <- loli.data %>%
        arrange(Pathway, Metabolite)

      loli.data_avg <- loli.data %>%
        arrange(Pathway, Metabolite) %>%
        mutate(Metab_name = row_number()) %>%
        group_by(Condition, Pathway) %>%
        mutate(
          avg = mean(Log2FC)
        ) %>%
        ungroup() %>%
        mutate(Pathway = factor(Pathway))

      loli.data_avg$Metab_name <- rep(seq.int(length(loli.data_avg$Metab_name)/length(Input_data)), each= length(Input_data) )


      loli_lines <-   loli.data_avg %>%
        arrange(Pathway, Metabolite) %>%
        group_by(Pathway) %>%
        summarize(
          start_x = min(Metab_name) -0.5,
          end_x = max(Metab_name) + 0.5,
          y = 0#unique(avg)
        ) %>%
        pivot_longer(
          cols = c(start_x, end_x),
          names_to = "type",
          values_to = "x"
        ) %>%
        mutate(
          x_group = if_else(type == "start_x", x + .1, x - .1),
          x_group = if_else(type == "start_x" & x == min(x), x_group - .1, x_group),
          x_group = if_else(type == "end_x" & x == max(x), x_group + .1, x_group) )

      p2 <- loli.data_avg %>%
        ggplot(aes(Metab_name, Log2FC)) + # names in aes ro Metab_name
        geom_hline(
          data = tibble(y = -5:5),
          aes(yintercept = y),
          color = "grey82",
          size = .5 )

      p2 <- p2 + geom_segment(
        aes(
          xend = Metab_name,          # names
          yend = 0,#avg,
          color = Pathway,
          #color = after_scale(colorspace::lighten(color, .2))
        ))

      p2 <- p2 + # geom_line( data = loli_lines, aes(x, y),  color = "grey40"  ) +
        geom_line(
          data = loli_lines,
          aes( x_group, y,
               color = Pathway,
               #  color = after_scale(colorspace::darken(color, .2))
          ), size = 2.5) +  geom_point(aes(size = p.adj, color = Pathway, shape = Condition  )
          )+ scale_size(trans = 'reverse')


      p2<- p2 + coord_flip()
      p2<-p2+ theme(axis.text.x=element_text())

      lab_metab <- loli.data_avg  %>% group_by(Metabolite,Metab_name) %>% summarise(max= Log2FC[which.max(abs(Log2FC))])


      lab_pos_metab <- lab_metab[rep(seq_len(nrow(lab_metab)), each = length(Input_data)), ] %>% filter(max>0)
      lab_neg_metab <- lab_metab[rep(seq_len(nrow(lab_metab)), each = length(Input_data)), ] %>% filter(max<0)

      p2<- p2+ annotate("text", x = lab_pos_metab$Metab_name, y = lab_pos_metab$max+1, label = lab_pos_metab$Metabolite, size = 4)
      p2<- p2+ annotate("text", x = lab_neg_metab$Metab_name, y = lab_neg_metab$max-1, label = lab_neg_metab$Metabolite, size = 4)

      p2

      p2 <- p2+ annotate("text", x = max(lab_neg_metab$Metab_name)+10, y = 0, label = "Significantly changed metabolites and their pathways", size = 8)

      p2 <- p2+Theme

      if(OutputPlotName ==""){
        ggsave(file=paste(Results_folder_plots_Lolipop_folder,"/", "Lolipop_Together", OutputPlotName, ".",Save_as, sep=""), plot=p2, width=20, height=20)
      }else{
        ggsave(file=paste(Results_folder_plots_Lolipop_folder,"/", "Lolipop_Together_", OutputPlotName, ".",Save_as, sep=""), plot=p2,  width=20, height=20)
      }
    }
    else if(Plot_pathways == "None"){

      Combined_Input <- data.frame(matrix(ncol = 4, nrow = 0))
      comb.colnames <- c("Metabolite","Log2FC",test, "Condition")
      colnames(Combined_Input) <- comb.colnames

      for (i in 1:length(Input_data)){
        Input_data[[i]]$Condition <- CondNames[[i]]
        Combined_Input <- rbind(Combined_Input, Input_data[[i]] %>% select(all_of(c("Metabolite","Log2FC",test, "Condition"))))

      }

      # Combined_Input <- Combined_Input %>% filter(abs(Log2FC) >=FCcutoff)
      # Combined_Input <- Combined_Input[Combined_Input[test] <= pCutoff,]
      # Combined_Input<- Combined_Input %>% drop_na()
      Combined_Input[test] <- round(Combined_Input[test], digits = 6)

      # Remove the metabolite with inf in logFC because it messes the plot
      Combined_Input <- Combined_Input[is.finite(Combined_Input$Log2FC),]

      Dotplot1 <- ggplot(Combined_Input, aes(x=reorder(Metabolite, + `Log2FC`), y=`Log2FC`, label=`p.adj`)) +
        geom_point(stat = 'identity', aes(size = abs(`Log2FC`), col = Condition))  +
        geom_segment(aes(y =0,
                         x = Metabolite,
                         yend = `Log2FC`,
                         xend = Metabolite),
                     color = "black") +
        scale_size(name="abs(Log2FC)",range = c(6,16))+
        geom_text(color="black", size=2) +
        labs(title=paste(Comparison)) +
        ylim(((Reduce(min,Combined_Input$`Log2FC`))-0.5),((Reduce(max,Combined_Input$`Log2FC`))+0.5)) +
        theme_minimal() +
        coord_flip()+
        theme(plot.title = element_text(color = "black", size = 12, face = "bold"),
              plot.subtitle = element_text(color = "black", size=10),
              plot.caption = element_text(color = "black",size=9, face = "italic", hjust = 2.5))+
        labs(y="Log2FC", x="")

      Dotplot1

      if(OutputPlotName ==""){
        ggsave(file=paste(Results_folder_plots_Lolipop_folder,"/", "Lolipop", ".",Save_as, sep=""), plot=Dotplot1, width=12, height=14)
      }else{
        ggsave(file=paste(Results_folder_plots_Lolipop_folder,"/", "Lolipop_", OutputPlotName, ".",Save_as, sep=""), plot=Dotplot1, width=12, height=14)
      }

    }
  }
}


###########----------------------
# Use function
#Lolipop(Input_data = list(DMA_output_Intra,DMA_output_Intra2,DMA_output_Intra3), CondNames = list("CvsRot", "Cvs3NPA", "CvsRot24h"), Input_pathways = DMA_output_Intra, Plot_pathways = "Individual")
######---------------------------



###############################
### ### ### Heatmap ### ### ###
###############################

#' @param Input_data DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Includes experimental design and outlier column.
#' @param Experimental_design DF which contains information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames. Column "Conditions" with information about the sample conditions (e.g. "N" and "T" or "Normal" and "Tumor"), can be used for feature filtering and colour coding in the PCA. Column "AnalyticalReplicate" including numerical values, defines technical repetitions of measurements, which will be summarised. Column "BiologicalReplicates" including numerical values. Please use the following names: "Conditions", "Biological_Replicates", "Analytical_Replicates".
#' @param Clustering_Condition String vector with the name/s of the Clustering conditions to be used from the Experimental Design.  \strong{Default = Conditions}
#' @param Input_pathways \emph{Optional: } DF which contains a 'Metabolite' and a 'Pathway' columns, with pathway information for each metabolite. \strong{Default = NULL}
#' @param Plot_pathways \emph{Optional: } String with plotting information about Metabolite pathways. Available only when the Input_pathways parameter has a file. Options are "None" if no pathways are to be plotted, "Individual" for plots of each Individual pathway and "Together" for metabolite pathways color-coded on a single volcano plot \strong{Default = "Together"}
#' @param Output_Name \emph{Optional: } String which is added to the output files of the plot
#' @param kMEAN \emph{Optional: } Vector of values for the values of k for k-means clustering of rows (Metabolites). \strong{Default = NA}
#' @param SCALE \emph{Optional: } String with the ??? \strong{Default = row}
#' @param Save_as \emph{Optional: } Select the file type of output plots. Options are svg or pdf. \strong{Default = svg}
#'
#' @keywords Volcano plot, pathways
#' @export
#'
#'

VizHeatmap <- function(Input_data,
                       Experimental_design,
                    Clustering_Condition = "Conditions",
                    Input_pathways = NULL,
                    Plot_pathways = "None",# or "Individual" or "Together=
                    Clustering_method = "single",
                    OutputPlotName= "",
                    kMEAN = NA,
                    SCALE = "row",
                    Save_as = "svg"
                    ){


  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse", "writexl","pheatmap")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  suppressMessages(library(tidyverse))


  ## ------------ Check Input files ----------- ##
  #1. Input_data and Conditions
  if(any(duplicated(row.names(Input_data)))==TRUE){
    stop("Duplicated row.names of Input_data, whilst row.names must be unique")
  } else if("Conditions" %in% colnames(Experimental_design)==FALSE){
    stop("There is no column named `Conditions` in Experimental_design to obtain Condition1 and Condition2.")
  } else{
    Test_num <- apply(Input_data, 2, function(x) is.numeric(x))
    if((any(Test_num) ==  FALSE) ==  TRUE){
      stop("Input_data needs to be of class numeric")
    } else{
      Test_match <- merge(Experimental_design, Input_data, by.x = "row.names",by.y = "row.names", all =  FALSE) # Do the unique IDs of the "Input_data" match the row names of the "Experimental_design"?
      if(nrow(Test_match) ==  0){
        stop("row.names Input_data need to match row.names Experimental_design.")
      } else(
        data <- Input_data
      )
    }
    Design <- Experimental_design
  }



  if(Clustering_Condition %in% colnames(Experimental_design)==FALSE){
    stop("Check Inpit. The clustering congitions does not exist in the column names on the Input data.")
  }
  Save_as_options <- c("svg","pdf", "jpeg", "tiff", "png", "bmp", "wmf","eps", "ps", "tex" )
  if(Save_as %in% Save_as_options == FALSE){
    stop("Check input. The selected Save_as option is not valid. Please select one of the folowwing: ",paste(Save_as_options,collapse = ", "),"." )
  }
  Plot_pathways_options <- c("None", "Individual", "Together")
  if (Plot_pathways %in% Plot_pathways_options == FALSE){
    stop("Check Input the Plot_pathways option is incorrect. The Allowed options are the following: ",paste(Plot_pathways_options,collapse = ", "),"." )
  }
  if(is.null(Input_pathways) == TRUE){
    if (Plot_pathways != "None"){
      warning("No Input_pathways were added. Yet the Plot_pathways option was changed. This will have no effect on the plot")
      Plot_pathways = "None"
    }
  }
  if(is.null(Input_pathways) == FALSE){
    if("Metabolite" %in% colnames(Input_pathways) == FALSE){
      stop("Check Input. Metabolite column is missing from Input_pathways")
    }
    if("Pathway" %in% colnames(Input_pathways) == FALSE){
      stop("Check Input. Pathway column is missing from Input_pathways")
    }
    if (sum(duplicated(Input_pathways$Metabolite)) > 0){
      stop("Duplicated Metabolites found in the Input_Pathways. The Metabolites must be unique.")
    }
    if(identical(sort(colnames(data)), sort(Input_pathways$Metabolite)) == FALSE){
      warning("The Metabolite column in the Input_data is not the same as the Metabolite column in the Input_pathways. We will take into consideration only the common Metabolites.")
      # find common metabolites
      common_metabolites <- colnames(data)[ colnames(data) %in% Input_pathways$Metabolite]
      # Take the data that have both pval reslults and pathwayss
      data <- data %>% select(common_metabolites)
      Input_pathways <- Input_pathways %>% filter(Metabolite %in% common_metabolites)
    }
    Input_pathways <- Input_pathways %>% select(all_of(c("Metabolite", "Pathway")))
  }


  ## ------------ Create Output folders ----------- ##
  name <- paste0("MetaProViz_Results_",Sys.Date())
  WorkD <- getwd()
  Results_folder <- file.path(WorkD, name)
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)} # Make Results folder
  Results_folder_plots_Heatmaps_folder = file.path(Results_folder, "Heatmap")
  if (!dir.exists(Results_folder_plots_Heatmaps_folder)) {dir.create(Results_folder_plots_Heatmaps_folder)}  # check and create folder


  if(Plot_pathways == "Individual"){


    my_annot<- NULL
    for (i in Clustering_Condition){
      my_annot[i] <- Experimental_design %>% select(i) %>% as.data.frame()
    }
    my_annot<- as.data.frame(my_annot)
    rownames(my_annot) <- rownames(data)

    # my_paths<- Input_pathways
    # my_paths <- column_to_rownames(my_paths,"Metabolite" )


    for (path in unique(Input_pathways$Pathway)){

      #  path = unique(Input_pathways$Pathway)[1]
      selected_path <- Input_pathways %>% filter(Pathway == path)
      selected_path_metabs <-  colnames(data) [colnames(data) %in% selected_path$Metabolite]
      data_path <- data %>% select(all_of(selected_path_metabs))


      for (k in kMEAN){
        set.seed(1234)
        out <-pheatmap::pheatmap(t(data_path),
                       clustering_method =  "complete",
                       scale = SCALE,
                       kmeans_k = k,
                       clustering_distance_rows = "correlation",
                       annotation_col = my_annot,
                     #annotation_row = my_paths,
                       main = paste(path, " pathway", sep = ""))

        if(is.na(k)==FALSE){
          Metabolite_clusters <- out[["kmeans"]][["cluster"]] %>% as.data.frame()
          names(Metabolite_clusters) <- "Clusters"
          Mouse_Cluster_Analysis <- merge(t(data), Metabolite_clusters, by = 'row.names' )
          names(Mouse_Cluster_Analysis)[1] <- "Metabolite"
          Mouse_Cluster_Analysis <- column_to_rownames(Mouse_Cluster_Analysis,"Metabolite" )
          Mouse_Cluster_Analysis_selectec <- Mouse_Cluster_Analysis %>% t() %>% as.data.frame()
          Mouse_Cluster_Analysis_selectec <- rownames_to_column(Mouse_Cluster_Analysis_selectec, "Sample")
        }


        if(is.na(k)==FALSE){

          if(OutputPlotName ==""){
            ggsave(file=paste(Results_folder_plots_Heatmaps_folder,"/", "Heatmap_",path,"_kmeans=",k, ".",Save_as, sep=""), plot=out, width=10, height=12)
          }else{
            ggsave(file=paste(Results_folder_plots_Heatmaps_folder,"/", "Heatmap_",path, "_kmeans=",k,"_",OutputPlotName, ".",Save_as, sep=""), plot=out, width=10, height=12)
          }
          writexl::write_xlsx(Mouse_Cluster_Analysis_selectec, paste(Results_folder_plots_Heatmaps_folder,"/Heatmap_",path,"Clustering_k=",k,"_data.xlsx", sep=""),col_names = TRUE)
        }else{

          if(OutputPlotName ==""){
            ggsave(file=paste(Results_folder_plots_Heatmaps_folder,"/", "Heatmap_",path, ".",Save_as, sep=""), plot=out, width=10, height=12)
          }else{
            ggsave(file=paste(Results_folder_plots_Heatmaps_folder,"/", "Heatmap_",path,"_",OutputPlotName, ".",Save_as, sep=""), plot=out, width=10, height=12)
          }
        }
      }
    }


  }
  else if(Plot_pathways == "Together"){



    my_annot<- NULL
    for (i in Clustering_Condition){
      my_annot[i] <- Experimental_design %>% select(i) %>% as.data.frame()
    }
    my_annot<- as.data.frame(my_annot)
    rownames(my_annot) <- rownames(data)

    my_paths<- Input_pathways
    my_paths <- column_to_rownames(my_paths,"Metabolite" )


    for (k in kMEAN){
      set.seed(1234)
      out <-pheatmap::pheatmap(t(data),
                     clustering_method =  "complete",
                     scale = SCALE,
                     kmeans_k = k,
                     clustering_distance_rows = "correlation",
                     annotation_col = my_annot,
                     annotation_row = my_paths)

      if(is.na(k)==FALSE){
        Metabolite_clusters <- out[["kmeans"]][["cluster"]] %>% as.data.frame()
        names(Metabolite_clusters) <- "Clusters"
        Mouse_Cluster_Analysis <- merge(t(data), Metabolite_clusters, by = 'row.names' )
        names(Mouse_Cluster_Analysis)[1] <- "Metabolite"
        Mouse_Cluster_Analysis <- column_to_rownames(Mouse_Cluster_Analysis,"Metabolite" )
        Mouse_Cluster_Analysis_selectec <- Mouse_Cluster_Analysis %>% t() %>% as.data.frame()
        Mouse_Cluster_Analysis_selectec <- rownames_to_column(Mouse_Cluster_Analysis_selectec, "Sample")
      }


      if(is.na(k)==FALSE){

        if(OutputPlotName ==""){
          ggsave(file=paste(Results_folder_plots_Heatmaps_folder,"/", "Heatmap","_kmeans=",k, ".",Save_as, sep=""), plot=out, width=10, height=12)
        }else{
          ggsave(file=paste(Results_folder_plots_Heatmaps_folder,"/", "Heatmap_", "kmeans=",k,"_",OutputPlotName, ".",Save_as, sep=""), plot=out, width=10, height=12)
        }
        writexl::write_xlsx(Mouse_Cluster_Analysis_selectec, paste(Results_folder_plots_Heatmaps_folder,"/","Clustering_k=",k,"_t(data).xlsx", sep=""),col_names = TRUE)
      }else{

        if(OutputPlotName ==""){
          ggsave(file=paste(Results_folder_plots_Heatmaps_folder,"/", "Heatmap", ".",Save_as, sep=""), plot=out, width=10, height=12)
        }else{
          ggsave(file=paste(Results_folder_plots_Heatmaps_folder,"/", "Heatmap_",OutputPlotName, ".",Save_as, sep=""), plot=out, width=10, height=12)
        }
      }
    }

  }
  else if(Plot_pathways == "None"){


    my_annot<- NULL
    for (i in Clustering_Condition){
      my_annot[i] <- Experimental_design %>% select(i) %>% as.data.frame()
    }

    my_annot<- as.data.frame(my_annot)
    rownames(my_annot) <- rownames(data)

    for (k in kMEAN){
      set.seed(1234)
      out <-pheatmap::pheatmap(t(data),
                     clustering_method =  "complete",
                     scale = SCALE,
                     kmeans_k = k,
                     clustering_distance_rows = "correlation",
                     annotation_col = my_annot)

      if(is.na(k)==FALSE){
        Metabolite_clusters <- out[["kmeans"]][["cluster"]] %>% as.data.frame()
        names(Metabolite_clusters) <- "Clusters"
        Mouse_Cluster_Analysis <- merge(t(data), Metabolite_clusters, by = 'row.names' )
        names(Mouse_Cluster_Analysis)[1] <- "Metabolite"
        Mouse_Cluster_Analysis <- column_to_rownames(Mouse_Cluster_Analysis,"Metabolite" )
        Mouse_Cluster_Analysis_selectec <- Mouse_Cluster_Analysis %>% t() %>% as.data.frame()
        Mouse_Cluster_Analysis_selectec <- rownames_to_column(Mouse_Cluster_Analysis_selectec, "Sample")
      }


      if(is.na(k)==FALSE){

        if(OutputPlotName ==""){
          ggsave(file=paste(Results_folder_plots_Heatmaps_folder,"/", "Heatmap","_kmeans=",k, ".",Save_as, sep=""), plot=out, width=10, height=12)
        }else{
          ggsave(file=paste(Results_folder_plots_Heatmaps_folder,"/", "Heatmap_", "kmeans=",k,"_",OutputPlotName, ".",Save_as, sep=""), plot=out, width=10, height=12)
        }
        writexl::write_xlsx(Mouse_Cluster_Analysis_selectec, paste(Results_folder_plots_Heatmaps_folder,"/","Clustering_k=",k,"_t(data).xlsx", sep=""),col_names = TRUE)
      }else{

        if(OutputPlotName ==""){
          ggsave(file=paste(Results_folder_plots_Heatmaps_folder,"/", "Heatmap", ".",Save_as, sep=""), plot=out, width=10, height=12)
        }else{
          ggsave(file=paste(Results_folder_plots_Heatmaps_folder,"/", "Heatmap_",OutputPlotName, ".",Save_as, sep=""), plot=out, width=10, height=12)
        }
      }
    }
  }

}

###########----------------------
# Use function
# Heatmap(Input_data = preprocessing_output_Intra$Processed_data,  Input_pathways = pathwaysdf, Plot_pathways = "Individual")
######---------------------------


################################
### ### ### Barplot  ### ### ###
################################

#' @param Input_data DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Includes experimental design and outlier column.
#' @param Experimental_design DF which contains information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames. Column "Conditions" with information about the sample conditions (e.g. "N" and "T" or "Normal" and "Tumor"), can be used for feature filtering and colour coding in the PCA. Column "AnalyticalReplicate" including numerical values, defines technical repetitions of measurements, which will be summarised. Column "BiologicalReplicates" including numerical values. Please use the following names: "Conditions", "Biological_Replicates", "Analytical_Replicates".
#' @param Output_Name \emph{Optional: } String which is added to the output files of the plot.
#' @param Output_plots String with plot save information. Available options are "Individual" for plots of each Individual metabolite and "Together" for a pdf containing all the plots. \strong{Default = Together}
#' @param Selected_Conditions Vector with names of selected Conditions for the plot. \strong{Default = NULL}
#' @param Selected_Comparisons Logical, TRUE to use t.test between the Selected_Conditions or FALSE. \strong{Default = NULL}
#' @param Theme \emph{Optional: } Selection of theme for plot. \strong{Default = theme_classic} ??
#' @param Save_as \emph{Optional: } Select the file type of output plots. Options are svg or pdf. \strong{Default = svg}
#'
#' @keywords Barplot
#' @export

VizBarplot <- function(Input_data,
                       Experimental_design,
                    OutputPlotName = "",
                    Output_plots = "Together",
                    Selected_Conditions = NULL,
                    Selected_Comparisons = NULL,
                    Theme = theme_classic(),
                    Save_as = "svg"
){

  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse", "ggplot2", "ggpubr")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  suppressMessages(library(tidyverse))


  ## ------------ Check Input files ----------- ##
  #1. Input_data and Conditions
  if(any(duplicated(row.names(Input_data)))==TRUE){
    stop("Duplicated row.names of Input_data, whilst row.names must be unique")
  } else if("Conditions" %in% colnames(Experimental_design)==FALSE){
    stop("There is no column named `Conditions` in Experimental_design to obtain Condition1 and Condition2.")
  } else{
    Test_num <- apply(Input_data, 2, function(x) is.numeric(x))
    if((any(Test_num) ==  FALSE) ==  TRUE){
      stop("Input_data needs to be of class numeric")
    } else{
      Test_match <- merge(Experimental_design, Input_data, by.x = "row.names",by.y = "row.names", all =  FALSE) # Do the unique IDs of the "Input_data" match the row names of the "Experimental_design"?
      if(nrow(Test_match) ==  0){
        stop("row.names Input_data need to match row.names Experimental_design.")
      } else(
        data <- Input_data
      )
    }
    Experimental_design <- Experimental_design
  }




  Output_plots_options <- c("Individual", "Together")
  if (Output_plots %in% Output_plots_options == FALSE){
    stop("Check Input the Plot_pathways option is incorrect. The Allowed options are the following: ",paste(Output_plots_options,collapse = ", "),"." )
  }
  if("Conditions" %in% colnames(Experimental_design)==FALSE){
    stop("There is no column named `Conditions` in Input_data.")
  }
  if(is.null(Selected_Conditions)==FALSE){
    for (Conditions in Selected_Conditions){
      if(Conditions %in% Experimental_design$Conditions==FALSE){
        stop("Check Input. The Selected_Conditions were not found in the Conditions Column.")
      }
    }
  }
  Save_as_options <- c("svg","pdf", "jpeg", "tiff", "png", "bmp", "wmf","eps", "ps", "tex" )
  if(Save_as %in% Save_as_options == FALSE){
    stop("Check input. The selected Save_as option is not valid. Please select one of the folowwing: ",paste(Save_as_options,collapse = ", "),"." )
  }

  ## ------------ Create Output folders ----------- ##
  name <- paste0("MetaProViz_Results_",Sys.Date())
  WorkD <- getwd()
  Results_folder <- file.path(WorkD, name)
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)} # Make Results folder
  Results_folder_plots_Barplots_folder = file.path(Results_folder, "Barplot")
  if (!dir.exists(Results_folder_plots_Barplots_folder)) {dir.create(Results_folder_plots_Barplots_folder)}  # check and create folder


  Metabolite_Names <- colnames(data)

  # make a list for plotting all plots together
  outlier_plot_list <- list()
  k=1

  for (i in Metabolite_Names){
    # i = Metabolite_Names[1]

    barplotdataMeans <- data %>%  select(all_of(i)) %>%  # Get mean & standard deviation by group
      group_by(Conditions=Experimental_design$Conditions) %>%
      summarise_at(vars(i), list(mean = mean, sd = sd)) %>%
      as.data.frame()

    barplotdata <- data %>%  select(i) %>%  group_by(Conditions=Experimental_design$Conditions)  %>%
      as.data.frame()
    names(barplotdata) <- c("Intensity", "Conditions")

    if (is.null(Selected_Conditions) == "FALSE"){
      barplotdataMeans <- barplotdataMeans %>% filter(Conditions %in% Selected_Conditions)
      barplotdata <- barplotdata %>% filter(Conditions %in% Selected_Conditions)
    }
    names(barplotdataMeans)[2] <- "Intensity"


    if(is.null(Selected_Comparisons)== TRUE){
      # names(barplotdataMeans)[2] <- "Intensity"
      # a <- max(barplotdataMeans$Intensity)
      barplot <- ggplot(barplotdata, aes(x = factor(Conditions), y = Intensity)) +
        geom_bar(stat = "summary", fun = "mean", fill = "skyblue") +
        geom_errorbar(data = barplotdataMeans, aes(x=Conditions, ymin=Intensity-sd, ymax=Intensity+sd), width=0.4, colour="black", alpha=0.9, size=0.5)+
        theme(legend.position = "right")+xlab("Conditions")+ ylab("Mean Intensity")

    }else{

      # names(barplotdataMeans)[2] <- "Intensity"
      # a <- max(barplotdataMeans$Intensity)
      barplot <- ggplot(barplotdata, aes(x = factor(Conditions), y = Intensity)) +
        geom_bar(stat = "summary", fun = "mean", fill = "skyblue") +
        geom_errorbar(data = barplotdataMeans, aes(x=Conditions, ymin=Intensity-sd, ymax=Intensity+sd), width=0.4, colour="black", alpha=0.9, size=0.5)+
        ggpubr::stat_compare_means(comparisons = Selected_Comparisons,
                                   label = "p.format", method = "t.test", hide.ns = TRUE, position = position_dodge(0.9), vjust = 0.25, show.legend = FALSE) +
        theme(legend.position = "right")+xlab("Conditions")+ ylab("Mean Intensity")
    }


    barplot <- barplot + Theme
    barplot <- barplot + theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1))
    barplot <- barplot + ggtitle(paste(i))



    if(Output_plots=="Individual"){

      i <- (gsub("/","_",i))#remove "/" cause this can not be safed in a PDF name
      i <- (gsub(":","_",i))

      if(OutputPlotName ==""){
        ggsave(file=paste(Results_folder_plots_Barplots_folder, "/",i, ".",Save_as, sep=""), plot=barplot, width=10, height=8)
      }else{
        ggsave(file=paste(Results_folder_plots_Barplots_folder, "/",OutputPlotName,"_",i, ".",Save_as, sep=""), plot=barplot, width=10, height=8)
      }


    } else if(Output_plots=="Together"){

      plot(barplot)
      outlier_plot_list[[k]] <- recordPlot()
      dev.off()
      k=k+1
    }
  }

  if(Output_plots=="Together"){
    if(OutputPlotName ==""){
      pdf(file= paste(Results_folder_plots_Barplots_folder,"/Barplots", OutputPlotName,".pdf", sep = ""), onefile = TRUE ) # or multivariate quality control chart
    }else{
      pdf(file= paste(Results_folder_plots_Barplots_folder,"/Barplots_", OutputPlotName,".pdf", sep = ""), onefile = TRUE ) # or multivariate quality control chart
    }
    for (plot in outlier_plot_list){
      replayPlot(plot)
    }
    dev.off()
  }
}

###--------------------------------
# Use function
#Barplot(Input_data = preprocessing_output_Intra$Processed_data)
#Barplot(Input_data = preprocessing_output_Intra$Processed_data, Selected_Conditions = c("Control", "Rot", "3NPA"),  Output_plots = "Individual")
#Barplot(Input_data = preprocessing_output_Intra$Processed_data, Selected_Conditions = c("Control", "Rot", "3NPA"), Selected_Comparisons = list(c(1,2), c(1,3), c(2,3)),  Output_plots = "Individual")
###-----------------------------------

################################
### ### ### Boxplots ### ### ###
################################

#' @param Input_data DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Includes experimental design and outlier column.
#' @param Experimental_design DF which contains information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames. Column "Conditions" with information about the sample conditions (e.g. "N" and "T" or "Normal" and "Tumor"), can be used for feature filtering and colour coding in the PCA. Column "AnalyticalReplicate" including numerical values, defines technical repetitions of measurements, which will be summarised. Column "BiologicalReplicates" including numerical values. Please use the following names: "Conditions", "Biological_Replicates", "Analytical_Replicates".
#' @param Output_Name \emph{Optional: } String which is added to the output files of the plot.
#' @param Output_plots String with plot save information. Available options are "Individual" for plots of each Individual metabolite and "Together" for a pdf containing all the plots. \strong{Default = Together}
#' @param Selected_Conditions Vector with names of selected Conditions for the plot. \strong{Default = NULL}
#' @param Selected_Comparisons Logical, TRUE to use t.test between the Selected_Conditions or FALSE. \strong{Default = NULL}
#' @param Theme \emph{Optional: } Selection of theme for plot. \strong{Default = theme_classic} ??
#' @param Save_as \emph{Optional: } Select the file type of output plots. Options are svg or pdf. \strong{Default = svg}
#'
#' @keywords Boxplot
#' @export
#'
VizBoxplot <- function(Input_data,
                       Experimental_design,
                    OutputPlotName = "",
                    Output_plots = "Together",
                    Selected_Conditions = NULL,
                    Selected_Comparisons = NULL,
                    Theme = theme_classic(),
                    Save_as = "svg"
){

  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse", "ggplot2", "ggpubr")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  suppressMessages(library(tidyverse))


  ## ------------ Check Input files ----------- ##
  #1. Input_data and Conditions
  RequiredPackages <- c("tidyverse", "ggplot2", "ggpubr")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  suppressMessages(library(tidyverse))


  ## ------------ Check Input files ----------- ##
  #1. Input_data and Conditions
  if(any(duplicated(row.names(Input_data)))==TRUE){
    stop("Duplicated row.names of Input_data, whilst row.names must be unique")
  } else if("Conditions" %in% colnames(Experimental_design)==FALSE){
    stop("There is no column named `Conditions` in Experimental_design to obtain Condition1 and Condition2.")
  } else{
    Test_num <- apply(Input_data, 2, function(x) is.numeric(x))
    if((any(Test_num) ==  FALSE) ==  TRUE){
      stop("Input_data needs to be of class numeric")
    } else{
      Test_match <- merge(Experimental_design, Input_data, by.x = "row.names",by.y = "row.names", all =  FALSE) # Do the unique IDs of the "Input_data" match the row names of the "Experimental_design"?
      if(nrow(Test_match) ==  0){
        stop("row.names Input_data need to match row.names Experimental_design.")
      } else(
        data <- Input_data
      )
    }
    Experimental_design <- Experimental_design
  }

  Output_plots_options <- c("Individual", "Together")
  if (Output_plots %in% Output_plots_options == FALSE){
    stop("Check Input the Plot_pathways option is incorrect. The Allowed options are the following: ",paste(Output_plots_options,collapse = ", "),"." )
  }
  if("Conditions" %in% colnames(Experimental_design)==FALSE){
    stop("There is no column named `Conditions` in Input_data.")
  }
  if(is.null(Selected_Conditions)==FALSE){
    for (Conditions in Selected_Conditions){
      if(Conditions %in% Experimental_design$Conditions==FALSE){
        stop("Check Input. The Selected_Conditions were not found in the Conditions Column.")
      }
    }
  }
  Save_as_options <- c("svg","pdf", "jpeg", "tiff", "png", "bmp", "wmf","eps", "ps", "tex" )
  if(Save_as %in% Save_as_options == FALSE){
    stop("Check input. The selected Save_as option is not valid. Please select one of the folowwing: ",paste(Save_as_options,collapse = ", "),"." )
  }

  ## ------------ Create Output folders ----------- ##
  name <- paste0("MetaProViz_Results_",Sys.Date())
  WorkD <- getwd()
  Results_folder <- file.path(WorkD, name)
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)} # Make Results folder
  Results_folder_plots_Boxplots_folder = file.path(Results_folder, "Boxplot")
  if (!dir.exists(Results_folder_plots_Boxplots_folder)) {dir.create(Results_folder_plots_Boxplots_folder)}  # check and create folder


  Metabolite_Names <- colnames(data)

  # make a list for plotting all plots together
  box_plot_list <- list()
  k=1

  for (i in Metabolite_Names){
    # i = Metabolite_Names[2]

    boxplotdata <- data %>%  select(all_of(i)) %>%                        # Get mean & standard deviation by group
      group_by(Conditions=Experimental_design$Conditions)
    names(boxplotdata) <- c("Intensity", "Conditions")

    if (is.null(Selected_Conditions) == "FALSE"){
      boxplotdata <- boxplotdata %>% filter(Conditions %in% Selected_Conditions)
    }




    #names(barplotdataMeans)[2] <- "Intensity"
    if(is.null(Selected_Comparisons)== TRUE){
      # names(barplotdataMeans)[2] <- "Intensity"
      # a <- max(barplotdataMeans$Intensity)
      boxplot <- ggplot(boxplotdata, aes(x=Conditions, y=Intensity)) +
        geom_boxplot(fill="skyblue") +
        geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.7)+xlab("Conditions")+ ylab("Mean Intensity")

    }else{

      # names(barplotdataMeans)[2] <- "Intensity"
      # a <- max(barplotdataMeans$Intensity)
      boxplot <- ggplot(boxplotdata, aes(x=Conditions, y=Intensity)) +
        geom_boxplot(fill="skyblue") +
        geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.7)+xlab("Conditions")+ ylab("Mean Intensity")+
        ggpubr::stat_compare_means(comparisons = Selected_Comparisons,
                                   label = "p.format", method = "t.test", hide.ns = TRUE, position = position_dodge(0.9), vjust = 0.25, show.legend = FALSE) +
        theme(legend.position = "right")+xlab("Conditions")+ ylab("Mean Intensity")
    }




    boxplot <- boxplot + Theme
    boxplot <- boxplot + theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1))
    boxplot <- boxplot + ggtitle(paste(i))


    if(Output_plots=="Individual"){

      i <- (gsub("/","_",i))#remove "/" cause this can not be safed in a PDF name
      i <- (gsub(":","_",i))

      if(OutputPlotName ==""){
        ggsave(file=paste(Results_folder_plots_Boxplots_folder, "/",i, ".",Save_as, sep=""), plot=boxplot, width=10, height=8)
      }else{
        ggsave(file=paste(Results_folder_plots_Boxplots_folder, "/",OutputPlotName,"_",i, ".",Save_as, sep=""), plot=boxplot, width=10, height=8)
      }

    } else if(Output_plots=="Together"){

      plot(boxplot)
      box_plot_list[[k]] <- recordPlot()
      dev.off()
      k=k+1
    }
  }


  if(Output_plots=="Together"){
    if(OutputPlotName ==""){
      pdf(file= paste(Results_folder_plots_Boxplots_folder,"/Boxplots", OutputPlotName,".pdf", sep = ""), onefile = TRUE )
    }else{
      pdf(file= paste(Results_folder_plots_Boxplots_folder,"/Boxplots", OutputPlotName,".pdf", sep = ""), onefile = TRUE )
    }
    for (plot in box_plot_list){
      replayPlot(plot)
    }
    dev.off()
  }

}

###--------------------------------
# Use function
#Boxplot(Input_data = preprocessing_output_Intra$Processed_data)
#Boxplot(Input_data = preprocessing_output_Intra$Processed_data, Selected_Conditions = c("Control", "Rot", "3NPA"),  Output_plots = "Individual")
#Boxplot(Input_data = preprocessing_output_Intra$Processed_data, Selected_Conditions = c("Control", "Rot", "3NPA"), Selected_Comparisons = list(c(1,2), c(1,3), c(2,3)),  Output_plots = "Together")
###-----------------------------------d_Conditions = c("Control", "Rot", "3NPA"), Selected_Comparisons = list(c(1,2), c(1,3)) )
#####----

####################################
### ### ### Violin Plots ### ### ###
#####################################

################################
### ### ### Violin ### ### ###
################################

#' @param Input_data DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Includes experimental design and outlier column.
#' @param Experimental_design DF which contains information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames. Column "Conditions" with information about the sample conditions (e.g. "N" and "T" or "Normal" and "Tumor"), can be used for feature filtering and colour coding in the PCA. Column "AnalyticalReplicate" including numerical values, defines technical repetitions of measurements, which will be summarised. Column "BiologicalReplicates" including numerical values. Please use the following names: "Conditions", "Biological_Replicates", "Analytical_Replicates".
#' @param Output_Name \emph{Optional: } String which is added to the output files of the plot.
#' @param Output_plots String with plot save information. Available options are "Individual" for plots of each Individual metabolite and "Together" for a pdf containing all the plots. \strong{Default = Together}
#' @param Selected_Conditions Vector with names of selected Conditions for the plot. \strong{Default = NULL}
#' @param Selected_Comparisons Logical, TRUE to use t.test between the Selected_Conditions or FALSE. \strong{Default = NULL}
#' @param Theme \emph{Optional: } Selection of theme for plot. \strong{Default = theme_classic} ??
#' @param Save_as \emph{Optional: } Select the file type of output plots. Options are svg or pdf. \strong{Default = svg}
#'
#' @keywords Violinplot
#' @export
#'
VizViolinplot <- function(Input_data,
                          Experimental_design,
                       OutputPlotName = "",
                       Output_plots = "Together",
                       Selected_Conditions = NULL,
                       Selected_Comparisons = NULL,
                       Theme = theme_classic(),
                       Save_as = "svg"
                       ){

  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse", "ggplot2", "ggpubr")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  suppressMessages(library(tidyverse))


  ## ------------ Check Input files ----------- ##
  #1. Input_data and Conditions
  RequiredPackages <- c("tidyverse", "ggplot2", "ggpubr")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  suppressMessages(library(tidyverse))


  ## ------------ Check Input files ----------- ##
  #1. Input_data and Conditions
  if(any(duplicated(row.names(Input_data)))==TRUE){
    stop("Duplicated row.names of Input_data, whilst row.names must be unique")
  } else if("Conditions" %in% colnames(Experimental_design)==FALSE){
    stop("There is no column named `Conditions` in Experimental_design to obtain Condition1 and Condition2.")
  } else{
    Test_num <- apply(Input_data, 2, function(x) is.numeric(x))
    if((any(Test_num) ==  FALSE) ==  TRUE){
      stop("Input_data needs to be of class numeric")
    } else{
      Test_match <- merge(Experimental_design, Input_data, by.x = "row.names",by.y = "row.names", all =  FALSE) # Do the unique IDs of the "Input_data" match the row names of the "Experimental_design"?
      if(nrow(Test_match) ==  0){
        stop("row.names Input_data need to match row.names Experimental_design.")
      } else(
        data <- Input_data
      )
    }
    Experimental_design <- Experimental_design
  }

  Output_plots_options <- c("Individual", "Together")
  if (Output_plots %in% Output_plots_options == FALSE){
    stop("Check Input the Plot_pathways option is incorrect. The Allowed options are the following: ",paste(Output_plots_options,collapse = ", "),"." )
  }
  if("Conditions" %in% colnames(Experimental_design)==FALSE){
    stop("There is no column named `Conditions` in Input_data.")
  }
  if(is.null(Selected_Conditions)==FALSE){
    for (Conditions in Selected_Conditions){
      if(Conditions %in% Experimental_design$Conditions==FALSE){
        stop("Check Input. The Selected_Conditions were not found in the Conditions Column.")
      }
    }
  }
  Save_as_options <- c("svg","pdf", "jpeg", "tiff", "png", "bmp", "wmf","eps", "ps", "tex" )
  if(Save_as %in% Save_as_options == FALSE){
    stop("Check input. The selected Save_as option is not valid. Please select one of the folowwing: ",paste(Save_as_options,collapse = ", "),"." )
  }

  ## ------------ Create Output folders ----------- ##
  name <- paste0("MetaProViz_Results_",Sys.Date())
  WorkD <- getwd()
  Results_folder <- file.path(WorkD, name)
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)} # Make Results folder
  Results_folder_plots_Violinplots_folder = file.path(Results_folder, "Violinplot")
  if (!dir.exists(Results_folder_plots_Violinplots_folder)) {dir.create(Results_folder_plots_Violinplots_folder)}  # check and create folder


  Metabolite_Names <- colnames(data)

  # make a list for plotting all plots together
  violin_plot_list <- list()
  k=1

  for (i in Metabolite_Names){
    # i = Metabolite_Names[2]

    violinplotdata <- data %>%  select(all_of(i)) %>%                        # Get mean & standard deviation by group
      group_by(Conditions=Experimental_design$Conditions)
    names(violinplotdata) <- c("Intensity", "Conditions")

    if (is.null(Selected_Conditions) == "FALSE"){
      violinplotdata <- violinplotdata %>% filter(Conditions %in% Selected_Conditions)
    }

    if(is.null(Selected_Comparisons)== TRUE){
      violinplot <- ggplot(violinplotdata, aes(x=Conditions, y=Intensity)) +
        geom_violin(fill="skyblue",width = 1)  +
        geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.7)+xlab("Conditions")+ ylab("Mean Intensity")

    }else{

      # names(barplotdataMeans)[2] <- "Intensity"
      # a <- max(barplotdataMeans$Intensity)
      violinplot <- ggplot(violinplotdata, aes(x=Conditions, y=Intensity)) +
        geom_violin(fill="skyblue") +
        geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.7)+xlab("Conditions")+ ylab("Mean Intensity")+
        ggpubr::stat_compare_means(comparisons = Selected_Comparisons,
                                   label = "p.format", method = "t.test", hide.ns = TRUE, position = position_dodge(0.9), vjust = 0.25, show.legend = FALSE) +
        theme(legend.position = "right")+xlab("Conditions")+ ylab("Mean Intensity")
    }




    violinplot <- violinplot + Theme
    violinplot <- violinplot + theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1))
    violinplot <- violinplot + ggtitle(paste(i))


    if(Output_plots=="Individual"){

      i <- (gsub("/","_",i))#remove "/" cause this can not be safed in a PDF name
      i <- (gsub(":","_",i))

      if(OutputPlotName ==""){
        ggsave(file=paste(Results_folder_plots_Violinplots_folder, "/",i, ".",Save_as, sep=""), plot=violinplot, width=10, height=8)
      }else{
        ggsave(file=paste(Results_folder_plots_Violinplots_folder, "/",OutputPlotName,"_",i, ".",Save_as, sep=""), plot=violinplot, width=10, height=8)
      }

    } else if(Output_plots=="Together"){

      plot(violinplot)
      violin_plot_list[[k]] <- recordPlot()
      dev.off()
      k=k+1
    }
  }


  if(Output_plots=="Together"){
    if(OutputPlotName ==""){
      pdf(file= paste(Results_folder_plots_violinplots_folder,"/Violinplots", OutputPlotName,".pdf", sep = ""), onefile = TRUE )
    }else{
      pdf(file= paste(Results_folder_plots_violinplots_folder,"/Violinplots", OutputPlotName,".pdf", sep = ""), onefile = TRUE )
    }
    for (plot in violin_plot_list){
      replayPlot(plot)
    }
    dev.off()
  }

}

###--------------------------------
# Use function
#Violinplot(Input_data = preprocessing_output_Intra$Processed_data)
#Violinplot(Input_data = preprocessing_output_Intra$Processed_data, Selected_Conditions = c("Control", "Rot", "3NPA"),  Output_plots = "Individual")
#Violinplot(Input_data = preprocessing_output_Intra$Processed_data, Selected_Conditions = c("Control", "Rot", "3NPA"), Selected_Comparisons = list(c(1,2), c(1,3), c(2,3)),  Output_plots = "Together")
###-----------------------------------d_Conditions = c("Control", "Rot", "3NPA"), Selected_Comparisons = list(c(1,2), c(1,3)) )
#####----



###################################
### ### ### super plots ### ### ###
###################################

#' @param Input_data DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Includes experimental design and outlier column.
#' @param Experimental_design DF which contains information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames. Column "Conditions" with information about the sample conditions (e.g. "N" and "T" or "Normal" and "Tumor"), can be used for feature filtering and colour coding in the PCA. Column "AnalyticalReplicate" including numerical values, defines technical repetitions of measurements, which will be summarised. Column "BiologicalReplicates" including numerical values. Please use the following names: "Conditions", "Biological_Replicates", "Analytical_Replicates".
#' @param Output_Name \emph{Optional: } String which is added to the output files of the plot.
#' @param Output_plots String with plot save information. Available options are "Individual" for plots of each Individual metabolite and "Together" for a pdf containing all the plots. \strong{Default = Together}
#' @param Selected_Conditions Vector with names of selected Conditions for the plot. \strong{Default = NULL}
#' @param Selected_Comparisons Logical, TRUE to use t.test between the Selected_Conditions or FALSE. \strong{Default = NULL}
#' @param Theme \emph{Optional: } Selection of theme for plot. \strong{Default = theme_classic} ??
#' @param Save_as \emph{Optional: } Select the file type of output plots. Options are svg or pdf. \strong{Default = svg}
#'
#' @keywords Barplot
#' @export
#'
VizSuperplot <- function(Input_data,
                         Experimental_design,
                      OutputPlotName = "",
                      Output_plots = "Together",
                      Selected_Conditions = NULL,
                      Selected_Comparisons = NULL,
                      Theme = theme_classic(),
                      Save_as = "svg"
){

  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse", "ggplot2", "ggpubr", "ggbeeswarm")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  suppressMessages(library(tidyverse))


  ## ------------ Check Input files ----------- ##
  #1. Input_data and Conditions
  RequiredPackages <- c("tidyverse", "ggplot2", "ggpubr")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  suppressMessages(library(tidyverse))


  ## ------------ Check Input files ----------- ##
  #1. Input_data and Conditions
  if(any(duplicated(row.names(Input_data)))==TRUE){
    stop("Duplicated row.names of Input_data, whilst row.names must be unique")
  } else if("Conditions" %in% colnames(Experimental_design)==FALSE){
    stop("There is no column named `Conditions` in Experimental_design to obtain Condition1 and Condition2.")
  } else{
    Test_num <- apply(Input_data, 2, function(x) is.numeric(x))
    if((any(Test_num) ==  FALSE) ==  TRUE){
      stop("Input_data needs to be of class numeric")
    } else{
      Test_match <- merge(Experimental_design, Input_data, by.x = "row.names",by.y = "row.names", all =  FALSE) # Do the unique IDs of the "Input_data" match the row names of the "Experimental_design"?
      if(nrow(Test_match) ==  0){
        stop("row.names Input_data need to match row.names Experimental_design.")
      } else(
        data <- Input_data
      )
    }
    Experimental_design <- Experimental_design
  }

  Output_plots_options <- c("Individual", "Together")
  if (Output_plots %in% Output_plots_options == FALSE){
    stop("Check Input the Plot_pathways option is incorrect. The Allowed options are the following: ",paste(Output_plots_options,collapse = ", "),"." )
  }
  if("Conditions" %in% colnames(Experimental_design)==FALSE){
    stop("There is no column named `Conditions` in Input_data.")
  }
  if(is.null(Selected_Conditions)==FALSE){
    for (Conditions in Selected_Conditions){
      if(Conditions %in% Experimental_design$Conditions==FALSE){
        stop("Check Input. The Selected_Conditions were not found in the Conditions Column.")
      }
    }
  }
  Save_as_options <- c("svg","pdf", "jpeg", "tiff", "png", "bmp", "wmf","eps", "ps", "tex" )
  if(Save_as %in% Save_as_options == FALSE){
    stop("Check input. The selected Save_as option is not valid. Please select one of the folowwing: ",paste(Save_as_options,collapse = ", "),"." )
  }

  ## ------------ Create Output folders ----------- ##
  name <- paste0("MetaProViz_Results_",Sys.Date())
  WorkD <- getwd()
  Results_folder <- file.path(WorkD, name)
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)} # Make Results folder
  Results_folder_plots_Superplots_folder = file.path(Results_folder, "superplot")
  if (!dir.exists(Results_folder_plots_Superplots_folder)) {dir.create(Results_folder_plots_Superplots_folder)}  # check and create folder


  Metabolite_Names <- colnames(data)

  # make a list for plotting all plots together
  super_plot_list <- list()
  k=1

  for (i in Metabolite_Names){
    # i = Metabolite_Names[2]
    # superplotdata <- data %>%  select(all_of(i)) %>%                        # Get mean & standard deviation by group
    #   group_by(Conditions=Experimental_design$Conditions)
    # names(superplotdata) <- c("Intensity", "Conditions")

    # if (is.null(Selected_Conditions) == "FALSE"){
    #   superplotdata <- superplotdata %>% filter(Conditions %in% Selected_Conditions)
    # }


    cond_selected <- data %>% select(i)
    names(cond_selected) <- "metabolite"
    cond_selected$Conditions <-  Experimental_design$Conditions
    cond_selected$Biological_Replicates <- Experimental_design$Biological_Replicates

    if (is.null(Selected_Conditions) == "FALSE"){
      cond_selected <- cond_selected %>% filter(Conditions %in% Selected_Conditions)
    }


    ReplicateAverages <- cond_selected %>%
      group_by(Conditions, Biological_Replicates) %>% summarise_each(list(mean))

    if(is.null(Selected_Comparisons)== TRUE){
      superplot<- ggplot(cond_selected, aes(x=Conditions,y=metabolite,color=factor(Biological_Replicates))) +
        ggbeeswarm::geom_beeswarm(cex=1) + scale_colour_brewer(palette = "Set1") +
        ggbeeswarm::geom_beeswarm(data=ReplicateAverages, size=4)

    }else{
      superplot<- ggplot(cond_selected, aes(x=Conditions,y=metabolite,color=factor(Biological_Replicates))) +
        ggbeeswarm::geom_beeswarm(cex=1) + scale_colour_brewer(palette = "Set1") +
        ggbeeswarm::geom_beeswarm(data=ReplicateAverages, size=4)+
        ggpubr::stat_compare_means(comparisons = Selected_Comparisons,
                                   label = "p.format", method = "t.test", hide.ns = TRUE, position = position_dodge(0.9), vjust = 0.25, show.legend = FALSE) +
        theme(legend.position = "right")+xlab("Conditions")+ ylab("Mean Intensity")
    }

    superplot<-  superplot+ theme(legend.position="right")
    superplot <- superplot + Theme
    superplot <- superplot + theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1))
    superplot <- superplot + ggtitle(paste(i)) +ylab(NULL)


    if(Output_plots=="Individual"){

      i <- (gsub("/","_",i))#remove "/" cause this can not be safed in a PDF name
      i <- (gsub(":","_",i))

      if(OutputPlotName ==""){
        ggsave(file=paste(Results_folder_plots_Superplots_folder, "/",i, ".",Save_as, sep=""), plot=superplot, width=10, height=8)
      }else{
        ggsave(file=paste(Results_folder_plots_Superplots_folder, "/",OutputPlotName,"_",i, ".",Save_as, sep=""), plot=superplot, width=10, height=8)
      }

    } else if(Output_plots=="Together"){

      plot(superplot)
      super_plot_list[[k]] <- recordPlot()
      dev.off()
      k=k+1
    }
  }


  if(Output_plots=="Together"){
    if(OutputPlotName ==""){
      pdf(file= paste(Results_folder_plots_Superplots_folder,"/Superplots", OutputPlotName,".pdf", sep = ""), onefile = TRUE )
    }else{
      pdf(file= paste(Results_folder_plots_Superplots_folder,"/Superplots", OutputPlotName,".pdf", sep = ""), onefile = TRUE )
    }
    for (plot in super_plot_list){
      replayPlot(plot)
    }
    dev.off()
  }

}

###--------------------------------
# Use function
#superplot(Input_data = preprocessing_output_Intra$Processed_data)
#superplot(Input_data = preprocessing_output_Intra$Processed_data, Selected_Conditions = c("Control", "Rot", "3NPA"),  Output_plots = "Individual")
#superplot(Input_data = preprocessing_output_Intra$Processed_data, Selected_Conditions = c("Control", "Rot", "3NPA"), Selected_Comparisons = list(c(1,2), c(1,3), c(2,3)),  Output_plots = "Individual")
###-----------------------------------
