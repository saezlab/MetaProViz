## ---------------------------
##
## Script name: Visualization  PCA plot
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
#' This script allows you to perform PCA plot visualization using the results of the MetaProViz analysis


#################################
### ### ### PCA Plots ### ### ###
#################################

#' @param Plot_SettingsInfo \emph{Optional: } NULL or Named vector including at least one of those three information : c(color="ColumnName_Plot_SettingsFile", shape= "ColumnName_Plot_SettingsFile"). \strong{Default = NULL}
#' @param Plot_SettingsFile \emph{Optional: } DF which contains information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames. Column "Conditions" with information about the sample conditions (e.g. "N" and "T" or "Normal" and "Tumor"), can be used for feature filtering and colour coding in the PCA. Column "AnalyticalReplicate" including numerical values, defines technical repetitions of measurements, which will be summarised. Column "BiologicalReplicates" including numerical values. Please use the following names: "Conditions", "Biological_Replicates", "Analytical_Replicates".\strong{Default = NULL}
#' @param Input_data DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected. includes experimental design and outlier column.
#' @param color_palette \emph{Optional: } Provide customiced color-palette in vector format. \strong{Default = NULL}
#' @param shape_palette \emph{Optional: } Provide customiced shape-palette in vector format. \strong{Default = NULL}
#' @param Show_Loadings \emph{Optional: } TRUE or FALSE for whether PCA loadings are also plotted on the PCA (biplot) \strong{Default = FALSE}
#' @param Scaling \emph{Optional: } TRUE or FALSE for whether a data scaling is used \strong{Default = TRUE}
#' @param Theme \emph{Optional: } Selection of theme for plot, e.g. theme_grey(). You can check for complete themes here: https://ggplot2.tidyverse.org/reference/ggtheme.html. If default=NULL we use theme_classic(). \strong{Default = NULL}
#' @param color_scale \emph{Optional: } Either "continuous" or "discrete" colour scale. For numeric or integer you can choose either, for character you have to choose discrete. If set to NULL this is done automatically. \strong{Default = NULL}
#' @param OutputPlotName \emph{Optional: } String which is added to the output files of the PCA \strong{Default = ""}
#' @param Save_as_Plot \emph{Optional: } Select the file type of output plots. Options are svg, png, pdf or NULL. \strong{Default = svg}
#'
#' @keywords PCA
#' @export

VizPCA <- function(Plot_SettingsInfo= NULL,
                   Plot_SettingsFile= NULL,#Design
                   Input_data,
                   color_palette= NULL,
                   shape_palette=NULL,
                   Show_Loadings = FALSE,
                   Scaling = TRUE,
                   Theme=NULL,#theme_classic()
                   color_scale=NULL,
                   OutputPlotName= '',
                   Save_as_Plot = "svg"
){

  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse","ggfortify", "ggplot2")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  suppressMessages(library(tidyverse))
  suppressMessages(library("ggfortify"))
  suppressMessages(library("ggplot2"))

  ## ------------ Check Input files ----------- ##
  # 1. The input data:
  if(any(duplicated(row.names(Input_data))) ==  TRUE){# Is the "Input_data" has unique IDs as row names and numeric values in columns?
    stop("Duplicated row.names of Input_data, whilst row.names must be unique")
  }else{
    Input_data_m<- Input_data
    Test_num <- apply(Input_data_m, 2, function(x) is.numeric(x))
    if((any(Test_num) ==  FALSE) ==  TRUE){
      stop("Input_data needs to be of class numeric in all columns.")
    }
  }

  # 2. The Plot_settings: Plot_Settings, Plot_SettingInfo and Plot_SettingFile
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
  }

  if(is.vector(Plot_SettingsInfo)==FALSE & is.null(Plot_SettingsFile)==FALSE){
    stop("Plot_SettingsInfo must be named vector or NULL.")
  }

  #3. Check other plot-specific parameters:
  if(is.null(color_palette)){
    safe_colorblind_palette <- c("#88CCEE",  "#DDCC77","#661100",  "#332288", "#AA4499","#999933",  "#44AA99", "#882215",  "#6699CC", "#117733", "#888888","#CC6677", "black","gold1","darkorchid4","red","orange")
    #check that length is enough for what the user wants to colour
    #stop(" The maximum number of pathways in the Input_pathways must be less than ",length(safe_colorblind_palette),". Please summarize sub-pathways together where possible and repeat.")
  } else{
    safe_colorblind_palette <-color_palette
    #check that length is enough for what the user wants to colour
  }
  if(is.null(shape_palette)){
    safe_shape_palette <- c(15,17,16,18,25,7,8,11,12)
    #check that length is enough for what the user wants to shape
  } else{
    safe_shape_palette <-shape_palette
    #check that length is enough for what the user wants to shape
  }

  if(is.logical(Show_Loadings) == FALSE){
    stop("Check input. The Show_Loadings value should be either =TRUE if loadings are to also be shown in the PCA plot or = FALSE if not.")
  }
  if(is.logical(Scaling) == FALSE){
    stop("Check input. The Scaling value should be either =TRUE if data scaling is to be performed prior to the PCA or = FALSE if not.")
  }

  # Theme check
  if(is.null(Theme)==FALSE){
    Theme_options <- c("theme_grey()", "theme_gray()", "theme_bw()", "theme_linedraw()", "theme_light()", "theme_dark()", "theme_minimal()", "theme_classic()", "theme_void()", "theme_test()")
    if (Theme %in% Theme_options == FALSE){
      stop("Theme option is incorrect. You can check for complete themes here: https://ggplot2.tidyverse.org/reference/ggtheme.html. Options are the following: ",paste(Theme_options, collapse = ", "),"." )
    }
  }

  # color_scale check

  if (!is.null(Save_as_Plot)) {
    Save_as_Plot_options <- c("svg","pdf", "png")
    if(Save_as_Plot %in% Save_as_Plot_options == FALSE){
      stop("Check input. The selected Save_as_Plot option is not valid. Please select one of the folowing: ",paste(Save_as_Plot_options,collapse = ", "),"." )
    }
  }

  ## ------------ Create Output folders ----------- ##
  if (!is.null(Save_as_Plot)) {
    name <- paste0("MetaProViz_Results_",Sys.Date())
    WorkD <- getwd()
    Results_folder <- file.path(WorkD, name)
    if (!dir.exists(Results_folder)) {dir.create(Results_folder)} # Make Results folder
    Results_folder_plots_PCA_folder = file.path(Results_folder, "PCA")  # This searches for a folder called "Preprocessing" within the "Results" folder in the current working directory and if its not found it creates one
    if (!dir.exists(Results_folder_plots_PCA_folder)) {dir.create(Results_folder_plots_PCA_folder)}  # check and create folder
  }

  ############################################################################################################
  ## ----------- Set the plot parameters: ------------ ##
  if(is.null(Plot_SettingsFile)==FALSE){
    InputPCA  <- merge(x=Plot_SettingsFile%>%rownames_to_column("UniqueID") , y=Input_data%>%rownames_to_column("UniqueID"), by="UniqueID", all.y=TRUE)%>%
      column_to_rownames("UniqueID")
  }else{
    InputPCA  <- Input_data
  }

  #Prepare the color scheme:
  if("color" %in% names(Plot_SettingsInfo)==TRUE){
    #color that will be used
    color_select <- safe_colorblind_palette[1:length(unique(InputPCA$color))]

    # numeric scale or continuous
    if(is.null(color_scale)==TRUE){#Will be set automatically:
      if(is.numeric(InputPCA$color) == TRUE | is.integer(InputPCA$color) == TRUE){
        if(length(unique(InputPCA$color)) > 4){ # change this to change the number after which color is not distinct but continuous
          InputPCA$color <- as.numeric(InputPCA$color)
        }else{
          InputPCA$color <- as.factor(InputPCA$color)
        }
      }
    }else if(color_scale=="discrete"){
      InputPCA$color <- as.factor(InputPCA$color)
    }else if(color_scale=="continuous"){
      if(is.numeric(InputPCA$color) == TRUE | is.integer(InputPCA$color) == TRUE){
        InputPCA$color <- as.numeric(InputPCA$color)
      }else{
        InputPCA$color <- as.factor(InputPCA$color)
        warning("color_scale=continuous, but is.numeric or is.integer is FALSE, hence colour scale is set to discrete.")
      }
    }


    #assign column and legend name
    InputPCA  <- InputPCA%>%
      dplyr::rename(!!paste(Plot_SettingsInfo[["color"]]) :="color")
    Param_Col <-paste(Plot_SettingsInfo[["color"]])
  } else{
    color_select <- NULL
    Param_Col <- NULL
  }
  #Prepare the shape scheme:
  if("shape" %in% names(Plot_SettingsInfo)==TRUE){
    #shapes that will be used
    shape_select <- safe_shape_palette[1:length(unique(InputPCA$shape))]

    #character
    if (!is.character(InputPCA$shape)) {
      # Convert the column to character
      InputPCA$shape <- as.character(InputPCA$shape)
    }

    #assign column and legend name
    InputPCA  <- InputPCA%>%
      dplyr::rename(!!paste(Plot_SettingsInfo[["shape"]]) :="shape")
    Param_Sha <-paste(Plot_SettingsInfo[["shape"]])
  } else{
    shape_select <-NULL
    Param_Sha <-NULL
  }

  ## ----------- Make the  plot based on the choosen parameters ------------ ##
  #Make the plot:
  PCA <- autoplot(prcomp(as.matrix(Input_data_m), scale. = as.logical(Scaling)),
                  data= InputPCA,
                  colour = Param_Col,
                  fill =  Param_Col,
                  shape = Param_Sha,
                  size = 3,
                  alpha = 0.8,
                  label=T,
                  label.size=2.5,
                  label.repel = TRUE,
                  loadings= as.logical(Show_Loadings), #draws Eigenvectors
                  loadings.label = as.logical(Show_Loadings),
                  loadings.label.vjust = 1.2,
                  loadings.label.size=2.5,
                  loadings.colour="grey10",
                  loadings.label.colour="grey10") +
    scale_shape_manual(values=shape_select)+
    scale_color_manual(values=color_select)+
    ggtitle(paste(OutputPlotName)) +
    geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
    geom_vline(xintercept=0,  color = "black", linewidth=0.1)

  #Add the theme
  if(is.null(Theme)==FALSE){
    PCA <- PCA+Theme
  }else{
    PCA <- PCA+theme_classic()
  }

  #Set the total heights and widths
  Plot_Sized <- MetaProViz:::plotGrob_PCA(Input=PCA, Plot_SettingsInfo=Plot_SettingsInfo, OutputPlotName=OutputPlotName)
  PCA <-Plot_Sized[[3]]

  if (!is.null(Save_as_Plot)) {
    if(OutputPlotName ==""){
      ggsave(file=paste(Results_folder_plots_PCA_folder,"/", "PCA", OutputPlotName, ".",Save_as_Plot, sep=""), plot=PCA, width=Plot_Sized[[2]], height=Plot_Sized[[1]], unit="cm")
    }else{
      ggsave(file=paste(Results_folder_plots_PCA_folder,"/", "PCA_", OutputPlotName, ".",Save_as_Plot, sep=""), plot=PCA, width=Plot_Sized[[2]], height=Plot_Sized[[1]], unit="cm")
    }
  }
  plot(PCA)

  # Return DF which can be assigned, but which does not print when they are not assigned.
  invisible(PCA)
}


##############################################################
### ### ### PCA helper function: Internal Function ### ### ###
##############################################################

#' @param Input This is the ggplot object generated within the VizPCA function.
#' @param Plot_SettingsInfo Passed to VizPCA
#' @param Plot_SettingsInfo Passed to VizPCA
#'
#' @keywords PCA helper function
#' @noRd

plotGrob_PCA <- function(Input, Plot_SettingsInfo, OutputPlotName){
  #------- Set the total heights and widths
  #we need ggplot_grob to edit the gtable of the ggplot object. Using this we can manipulate the gtable arguments directly.
  plottable<- ggplot2::ggplotGrob(Input) # Convert the plot to a gtable
  if(is.null(Plot_SettingsInfo)==TRUE){
    #-----widths
    plottable$widths[5] <- unit(8, "cm")#controls x-axis
    plottable$widths[c(3)] <- unit(2,"cm")#controls margins --> y-axis label is there
    plottable$widths[c(1,2,4)] <- unit(0,"cm")#controls margins --> not needed
    plottable$widths[c(6)] <- unit(1,"cm")#controls margins --> start Figure legend
    plottable$widths[c(10)] <- unit(0,"cm")#controls margins --> Figure legend
    plottable$widths[c(7,8,9,11)] <- unit(0,"cm")#controls margins --> not needed
    plot_widths <- 11

    if((OutputPlotName=="")==FALSE){#Check how much width is needed for the figure title/subtitle
      character_count <- nchar(OutputPlotName)
      Titles_width <- (character_count*0.25)+0.8
      if(Titles_width>plot_widths){#If the title needs more space than the plot offers:
        plottable$widths[11] <- unit(Titles_width-plot_widths,"cm")#controls margins --> start Figure legend
        plot_widths <- Titles_width
      }
    }

    #-----heigths
    plottable$heights[7] <- unit(8, "cm")#controls x-axis
    plottable$heights[c(8)] <- unit(1,"cm")#controls margins --> x-axis label
    plottable$heights[c(10)] <- unit(1,"cm")#controls margins --> Figure caption
    plottable$heights[c(9,11,12)] <- unit(0,"cm")#controls margins --> not needed

    if(OutputPlotName==""){
      plottable$heights[c(6)] <- unit(0.5,"cm")#controls margins --> Some space above the plot
      plottable$heights[c(1,2,3,4,5)] <- unit(0,"cm")#controls margins --> not needed
      plot_heights <- 10.5
    } else{
      plottable$heights[c(3)] <- unit(1,"cm")#controls margins --> OutputPlotName and subtitle
      plottable$heights[c(1,2,4,5,6)] <- unit(0,"cm")#controls margins --> not needed
      plot_heights <-11
    }
  }else if("color" %in% names(Plot_SettingsInfo)==TRUE & "shape" %in% names(Plot_SettingsInfo)==TRUE){
    #------- Legend heights
    Legend <- ggpubr::get_legend(Input) # Extract legend to adjust separately
    Legend_heights <- (round(as.numeric(Legend$heights[3]),1))+(round(as.numeric(Legend$heights[5]),1))

    #-----Plot widths
    plottable$widths[5] <- unit(8, "cm")#controls x-axis
    plottable$widths[c(3)] <- unit(2,"cm")#controls margins --> y-axis label is there
    plottable$widths[c(1,2,4)] <- unit(0,"cm")#controls margins --> not needed
    plottable$widths[c(6)] <- unit(1,"cm")#controls margins --> start Figure legend
    plottable$widths[c(7,8,10,11)] <- unit(0,"cm")#controls margins --> not needed

    Value <- round(as.numeric(plottable$widths[9]),1) #plottable$widths[9] is a <unit/unit_v2> object and we can extract the extract the numeric part
    plot_widths <- 11+Value

    if((OutputPlotName=="")==FALSE){#Check how much width is needed for the figure title/subtitle
      character_count <- nchar(OutputPlotName)
      Titles_width <- (character_count*0.25)+0.8
      if(Titles_width>plot_widths){#If the title needs more space than the plot offers:
        plottable$widths[11] <- unit(Titles_width-plot_widths,"cm")#controls margins --> start Figure legend
        plot_widths <- Titles_width
      }
    }

    #-----Plot heigths
    plottable$heights[7] <- unit(8, "cm")#controls x-axis
    plottable$heights[c(8)] <- unit(1,"cm")#controls margins --> x-axis label
    plottable$heights[c(10)] <- unit(1,"cm")#controls margins --> Figure caption
    plottable$heights[c(9,11)] <- unit(0,"cm")#controls margins --> not needed

    if(OutputPlotName==""){
      plottable$heights[c(6)] <- unit(0.5,"cm")#controls margins --> Some space above the plot
      plottable$heights[c(2,3,4,5)] <- unit(0,"cm")#controls margins --> not needed

      if(Legend_heights>10.5){#If the legend requires more heights than the Plot
        Add <- (Legend_heights-10.5)/2
        plottable$heights[1] <- unit(Add,"cm")#controls margins --> Can be increased if Figure legend needs more space on the top
        plottable$heights[12] <- unit(Add,"cm")#controls margins --> Can be increased if Figure legend needs more space on the bottom
        plot_heights <- Legend_heights
      }else{
        plottable$heights[1] <- unit(0,"cm")#controls margins --> Can be increased if Figure legend needs more space on the top
        plottable$heights[12] <- unit(0,"cm")#controls margins --> Can be increased if Figure legend needs more space on the bottom
        plot_heights <- 10.5
      }
    } else{#If we do have Title and or subtitle
      plottable$heights[c(3)] <- unit(1,"cm")#controls margins --> OutputPlotName and subtitle
      plottable$heights[c(2,4,5,6)] <- unit(0,"cm")#controls margins --> not needed
      if(Legend_heights>11){#If the legend requires more heights than the Plot
        Add <- (Legend_heights-11)/2
        plottable$heights[1] <- unit(Add,"cm")#controls margins --> Can be increased if Figure legend needs more space on the top
        plottable$heights[12] <- unit(Add,"cm")#controls margins --> Can be increased if Figure legend needs more space on the bottom
        plot_heights <- Legend_heights
      }else{
        plottable$heights[1] <- unit(0,"cm")#controls margins --> Can be increased if Figure legend needs more space on the top
        plottable$heights[12] <- unit(0,"cm")#controls margins --> Can be increased if Figure legend needs more space on the bottom
        plot_heights <- 11
      }
    }
  }else if("color" %in% names(Plot_SettingsInfo)==TRUE | "shape" %in% names(Plot_SettingsInfo)==TRUE){
    #------- Legend heights
    Legend <- ggpubr::get_legend(Input) # Extract legend to adjust separately
    Legend_heights <- (round(as.numeric(Legend$heights[3]),1))

    #----- Plot widths
    plottable$widths[5] <- unit(8, "cm")#controls x-axis
    plottable$widths[c(3)] <- unit(2,"cm")#controls margins --> y-axis label is there
    plottable$widths[c(1,2,4)] <- unit(0,"cm")#controls margins --> not needed
    plottable$widths[c(6)] <- unit(1,"cm")#controls margins --> start Figure legend
    plottable$widths[c(7,8,10,11)] <- unit(0,"cm")#controls margins --> not needed

    Value <- round(as.numeric(plottable$widths[9]),1) #plottable$widths[9] is a <unit/unit_v2> object and we can extract the extract the numeric part
    plot_widths <- 11+Value

    if((OutputPlotName=="")==FALSE){#Check how much width is needed for the figure title/subtitle
      character_count <- nchar(OutputPlotName)
      Titles_width <- (character_count*0.25)+0.8
      if(Titles_width>plot_widths){#If the title needs more space than the plot offers:
        plottable$widths[11] <- unit(Titles_width-plot_widths,"cm")#controls margins --> start Figure legend
        plot_widths <- Titles_width
      }
    }

    #-----Plot heigths
    plottable$heights[7] <- unit(8, "cm")#controls x-axis
    plottable$heights[c(8)] <- unit(1,"cm")#controls margins --> x-axis label
    plottable$heights[c(10)] <- unit(1,"cm")#controls margins --> Figure caption
    plottable$heights[c(9,11)] <- unit(0,"cm")#controls margins --> not needed

    if(OutputPlotName==""){
      plottable$heights[c(6)] <- unit(0.5,"cm")#controls margins --> Some space above the plot
      plottable$heights[c(2,3,4,5)] <- unit(0,"cm")#controls margins --> not needed

      if(Legend_heights>10.5){#If the legend requires more heights than the Plot
        Add <- (Legend_heights-10.5)/2
        plottable$heights[1] <- unit(Add,"cm")#controls margins --> Can be increased if Figure legend needs more space on the top
        plottable$heights[12] <- unit(Add,"cm")#controls margins --> Can be increased if Figure legend needs more space on the bottom
        plot_heights <- Legend_heights
      }else{
        plottable$heights[1] <- unit(0,"cm")#controls margins --> Can be increased if Figure legend needs more space on the top
        plottable$heights[12] <- unit(0,"cm")#controls margins --> Can be increased if Figure legend needs more space on the bottom
        plot_heights <- 10.5
      }
    }else{#If we do have Title and or subtitle
      plottable$heights[c(3)] <- unit(1,"cm")#controls margins --> OutputPlotName and subtitle
      plottable$heights[c(2,4,5,6)] <- unit(0,"cm")#controls margins --> not needed
      if(Legend_heights>11){#If the legend requires more heights than the Plot
        Add <- (Legend_heights-11)/2
        plottable$heights[1] <- unit(Add,"cm")#controls margins --> Can be increased if Figure legend needs more space on the top
        plottable$heights[12] <- unit(Add,"cm")#controls margins --> Can be increased if Figure legend needs more space on the bottom
        plot_heights <- Legend_heights
      }else{
        plottable$heights[1] <- unit(0,"cm")#controls margins --> Can be increased if Figure legend needs more space on the top
        plottable$heights[12] <- unit(0,"cm")#controls margins --> Can be increased if Figure legend needs more space on the bottom
        plot_heights <- 11
      }
    }
  }
  #plot_param <-c(plot_heights=plot_heights, plot_widths=plot_widths)
  Output<- list(plot_heights, plot_widths, plottable)
}



