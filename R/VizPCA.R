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

#' @param Input_data DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected. includes experimental design and outlier column.
#' @param SettingsInfo \emph{Optional: } NULL or Named vector including at least one of those three information : c(color="ColumnName_Plot_SettingsFile", shape= "ColumnName_Plot_SettingsFile"). \strong{Default = NULL}
#' @param SettingsFile_Sample \emph{Optional: } DF which contains information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames. Column "Conditions" with information about the sample conditions (e.g. "N" and "T" or "Normal" and "Tumor"), can be used for feature filtering and colour coding in the PCA. Column "AnalyticalReplicate" including numerical values, defines technical repetitions of measurements, which will be summarised. Column "BiologicalReplicates" including numerical values. Please use the following names: "Conditions", "Biological_Replicates", "Analytical_Replicates".\strong{Default = NULL}
#' @param ColorPalette \emph{Optional: } Provide customiced color-palette in vector format. For continuous scale use e.g. scale_color_gradient(low = "#88CCEE", high = "red") and for discrete scale c("#88CCEE",  "#DDCC77","#661100",  "#332288")\strong{Default = NULL}
#' @param ShapePalette \emph{Optional: } Provide customiced shape-palette in vector format. \strong{Default = NULL}
#' @param ShowLoadings \emph{Optional: } TRUE or FALSE for whether PCA loadings are also plotted on the PCA (biplot) \strong{Default = FALSE}
#' @param Scaling \emph{Optional: } TRUE or FALSE for whether a data scaling is used \strong{Default = TRUE}
#' @param Theme \emph{Optional: } Selection of theme for plot, e.g. theme_grey(). You can check for complete themes here: https://ggplot2.tidyverse.org/reference/ggtheme.html. If default=NULL we use theme_classic(). \strong{Default = "discrete"}
#' @param ColorScale \emph{Optional: } Either "continuous" or "discrete" colour scale. For numeric or integer you can choose either, for character you have to choose discrete. \strong{Default = NULL}
#' @param PlotName \emph{Optional: } String which is added to the output files of the PCA \strong{Default = ""}
#' @param SaveAs_Plot \emph{Optional: } Select the file type of output plots. Options are svg, png, pdf or NULL. \strong{Default = svg}
#' @param PrintPlot \emph{Optional: } TRUE or FALSE, if TRUE Volcano plot is saved as an overview of the results. \strong{Default = TRUE}
#' @param FolderPath \emph{Optional:} Path to the folder the results should be saved at. \strong(Default = NULL)
#'
#' @keywords PCA
#' @export

VizPCA <- function(InputData,
                   SettingsInfo= NULL,
                   SettingsFile_Sample = NULL,
                   ColorPalette= NULL,
                   ColorScale="discrete",
                   ShapePalette=NULL,
                   ShowLoadings = FALSE,
                   Scaling = TRUE,
                   Theme=NULL,#theme_classic()
                   PlotName= '',
                   SaveAs_Plot = "svg",
                   PrintPlot=TRUE,
                   FolderPath = NULL
){

  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse","ggfortify", "ggplot2")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  suppressMessages(library(tidyverse))
  suppressMessages(library("ggfortify"))

  ################################################################################################################################################################################################
  ## ------------ Check Input files ----------- ##
  # HelperFunction `CheckInput`
  MetaProViz:::CheckInput(InputData=InputData,
                          SettingsFile_Sample=SettingsFile_Sample,
                          SettingsFile_Metab=NULL,
                          SettingsInfo=SettingsInfo,
                          SaveAs_Plot=SaveAs_Plot,
                          SaveAs_Table=NULL,
                          CoRe=FALSE,
                          PrintPlot= PrintPlot)

  # CheckInput` Specific
  if(is.logical(ShowLoadings) == FALSE){
    stop("Check input. The Show_Loadings value should be either =TRUE if loadings are to be shown on the PCA plot or = FALSE if not.")
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


  ## ------------ Create Results output folder ----------- ##
  if(is.null(SaveAs_Plot)==FALSE){
    Folder <- MetaProViz:::SavePath(FolderName= "PCA_Plots",
                                    FolderPath=FolderPath)
  }


  ############################################################################################################
  ## ----------- Set the plot parameters: ------------ ##
  ##--- Prepare colour and shape palette
  if(is.null(ColorPalette)){
    if((ColorScale=="discrete")==TRUE){
      safe_colorblind_palette <- c("#88CCEE",  "#DDCC77","#661100",  "#332288", "#AA4499","#999933",  "#44AA99", "#882215",  "#6699CC", "#117733", "#888888","#CC6677", "black","gold1","darkorchid4","red","orange")
    }else if(ColorScale=="continuous"){
      safe_colorblind_palette <- NULL
    }
  } else{
    safe_colorblind_palette <-ColorPalette
    #check that length is enough for what the user wants to colour
  }
  if(is.null(ShapePalette)){
    safe_shape_palette <- c(15,17,16,18,6,7,8,11,12)
    #check that length is enough for what the user wants to shape
  } else{
    safe_shape_palette <-ShapePalette
    #check that length is enough for what the user wants to shape
    #stop(" The maximum number of pathways in the Input_pathways must be less than ",length(safe_colorblind_palette),". Please summarize sub-pathways together where possible and repeat.")
  }

  ##--- Prepare the color scheme:
  if("color" %in% names(SettingsInfo)==TRUE & "shape" %in% names(SettingsInfo)==TRUE){
    if((SettingsInfo[["shape"]] == SettingsInfo[["color"]])==TRUE){
      SettingsFile_Sample$shape <- SettingsFile_Sample[,paste(SettingsInfo[["color"]])]
      SettingsFile_Sample<- SettingsFile_Sample%>%
        dplyr::rename("color"=paste(SettingsInfo[["color"]]))
    }else{
      SettingsFile_Sample <- SettingsFile_Sample%>%
          dplyr::rename("color"=paste(SettingsInfo[["color"]]),
                        "shape"=paste(SettingsInfo[["shape"]]))
      }
 }else if("color" %in% names(SettingsInfo)==TRUE & "shape" %in% names(SettingsInfo)==FALSE){
   if("color" %in% names(SettingsInfo)==TRUE){
     SettingsFile_Sample <- SettingsFile_Sample%>%
        dplyr::rename("color"=paste(SettingsInfo[["color"]]))
   }
   if("shape" %in% names(SettingsInfo)==TRUE){
      SettingsFile_Sample <- SettingsFile_Sample%>%
        dplyr::rename("shape"=paste(SettingsInfo[["shape"]]))
   }
 }


  ##--- Prepare Input Data:
  if(is.null(SettingsFile_Sample)==FALSE){
    InputPCA  <- merge(x=SettingsFile_Sample%>%rownames_to_column("UniqueID") , y=InputData%>%rownames_to_column("UniqueID"), by="UniqueID", all.y=TRUE)%>%
      column_to_rownames("UniqueID")
  }else{
    InputPCA  <- InputData
  }

 ##--- Prepare the color and shape settings:
  if("color" %in% names(SettingsFile_Sample)==TRUE){
    if(ColorScale=="discrete"){
      InputPCA$color <- as.factor(InputPCA$color)
      color_select <- safe_colorblind_palette[1:length(unique(InputPCA$color))]
    }else if(ColorScale=="continuous"){
      if(is.numeric(InputPCA$color) == TRUE | is.integer(InputPCA$color) == TRUE){
        InputPCA$color <- as.numeric(InputPCA$color)
        color_select <- safe_colorblind_palette
      }else{
        InputPCA$color <- as.factor(InputPCA$color)
        #Overwrite color pallette
        safe_colorblind_palette <- c("#88CCEE",  "#DDCC77","#661100",  "#332288", "#AA4499","#999933",  "#44AA99", "#882215",  "#6699CC", "#117733", "#888888","#CC6677", "black","gold1","darkorchid4","red","orange")
        #color that will be used for distinct
        color_select <- safe_colorblind_palette[1:length(unique(InputPCA$color))]
        #Overwrite color_scale
        color_scale<- "discrete"
        warning("ColorScale=continuous, but is.numeric or is.integer is FALSE, hence colour scale is set to discrete.")
      }
    }
  }

  if("shape" %in% names(SettingsFile_Sample)==TRUE){
    shape_select <- safe_shape_palette[1:length(unique(InputPCA$shape))]

    if (!is.character(InputPCA$shape)) {
        # Convert the column to character
        InputPCA$shape <- as.character(InputPCA$shape)
    }
  }

  ##---  #assign column and legend name
  if("color" %in% names(SettingsFile_Sample)==TRUE){
    InputPCA  <- InputPCA%>%
      dplyr::rename(!!paste(SettingsInfo[["color"]]) :="color")
    Param_Col <-paste(SettingsInfo[["color"]])
  } else{
    color_select <- NULL
    Param_Col <- NULL
  }

  if("shape" %in% names(SettingsFile_Sample)==TRUE){
    InputPCA  <- InputPCA%>%
      dplyr::rename(!!paste(SettingsInfo[["shape"]]) :="shape")
    Param_Sha <-paste(SettingsInfo[["shape"]])
  } else{
    shape_select <-NULL
    Param_Sha <-NULL
  }

  ## ----------- Make the  plot based on the choosen parameters ------------ ##
  PlotList <- list()#Empty list to store all the plots

  #Make the plot:
  PCA <- autoplot(prcomp(as.matrix(InputData), scale. = as.logical(Scaling)),
                  data= InputPCA,
                  colour = Param_Col,
                  fill =  Param_Col,
                  shape = Param_Sha,
                  size = 3,
                  alpha = 0.8,
                  label=T,
                  label.size=2.5,
                  label.repel = TRUE,
                  loadings= as.logical(ShowLoadings), #draws Eigenvectors
                  loadings.label = as.logical(ShowLoadings),
                  loadings.label.vjust = 1.2,
                  loadings.label.size=2.5,
                  loadings.colour="grey10",
                  loadings.label.colour="grey10" ) +
    scale_shape_manual(values=shape_select)+
    ggtitle(paste(PlotName)) +
    geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
    geom_vline(xintercept=0,  color = "black", linewidth=0.1)

    if(ColorScale=="discrete"){
      PCA <-PCA + scale_color_manual(values=color_select)
    }else if(ColorScale=="continuous" & is.null(color_palette)){
      PCA <-PCA + color_select
    }

  #Add the theme
  if(is.null(Theme)==FALSE){
    PCA <- PCA+Theme
  }else{
    PCA <- PCA+theme_classic()
  }

  ## Store the plot in the 'plots' list
  PlotList[[Plot]] <- PCA

  #Set the total heights and widths
  Plot_Sized <- MetaProViz:::PlotGrob_PCA(InputPlot=PCA, SettingsInfo=SettingsInfo, PlotName=PlotName)
  PCA <- Plot_Sized[[3]]
  PCA <- ggplot2::ggplot() +
    annotation_custom(PCA)
  PCA  <-PCA  + theme(panel.background = element_rect(fill = "transparent"))
  PlotList[[Plot_Sized]] <- PCA


  ######################################################################################################################################################################
  ##----- Save and Return
  #Here we make a list in which we will save the outputs:
  suppressMessages(suppressWarnings(
      MetaProViz:::SaveRes(InputList_DF=NULL,
                           InputList_Plot= PlotList[[Plot_Sized]],
                           SaveAs_Table=NULL,
                           SaveAs_Plot=SaveAs_Plot,
                           FolderPath= Folder,
                           FileName= paste("PCA_", PlotName),
                           CoRe=FALSE,
                           PrintPlot=PrintPlot,
                           PlotHeight=Plot_Sized[[1]],
                           PlotWidth=Plot_Sized[[2]],
                           PlotUnit="cm")))

  return(invisible(list("Plot"=PlotList,"Plot_Sized" = PlotList_adaptedGrid)))
}


##############################################################
### ### ### PCA helper function: Internal Function ### ### ###
##############################################################

#' @param InputPlot This is the ggplot object generated within the VizPCA function.
#' @param SettingsInfo Passed to VizPCA
#' @param PlotName Passed to VizPCA
#'
#' @keywords PCA helper function
#' @noRd

PlotGrob_PCA <- function(InputPlot, SettingsInfo, PlotName){
  #------- Set the total heights and widths
  #we need ggplot_grob to edit the gtable of the ggplot object. Using this we can manipulate the gtable arguments directly.
  plottable<- ggplot2::ggplotGrob(InputPlot) # Convert the plot to a gtable
  if(is.null(SettingsInfo)==TRUE){
    #-----widths
    plottable$widths[5] <- unit(8, "cm")#controls x-axis
    plottable$widths[c(3)] <- unit(2,"cm")#controls margins --> y-axis label is there
    plottable$widths[c(1,2,4)] <- unit(0,"cm")#controls margins --> not needed
    plottable$widths[c(6)] <- unit(1,"cm")#controls margins --> start Figure legend
    plottable$widths[c(10)] <- unit(0,"cm")#controls margins --> Figure legend
    plottable$widths[c(7,8,9,11)] <- unit(0,"cm")#controls margins --> not needed
    plot_widths <- 11

    if((PlotName=="")==FALSE){#Check how much width is needed for the figure title/subtitle
      character_count <- nchar(PlotName)
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

    if(PlotName==""){
      plottable$heights[c(6)] <- unit(0.5,"cm")#controls margins --> Some space above the plot
      plottable$heights[c(1,2,3,4,5)] <- unit(0,"cm")#controls margins --> not needed
      plot_heights <- 10.5
    } else{
      plottable$heights[c(3)] <- unit(1,"cm")#controls margins --> PlotName and subtitle
      plottable$heights[c(1,2,4,5,6)] <- unit(0,"cm")#controls margins --> not needed
      plot_heights <-11
    }
  }else if("color" %in% names(SettingsInfo)==TRUE & "shape" %in% names(SettingsInfo)==TRUE){
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

    if((PlotName=="")==FALSE){#Check how much width is needed for the figure title/subtitle
      character_count <- nchar(PlotName)
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

    if(PlotName==""){
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
      plottable$heights[c(3)] <- unit(1,"cm")#controls margins --> PlotName and subtitle
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
  }else if("color" %in% names(SettingsInfo)==TRUE | "shape" %in% names(SettingsInfo)==TRUE){
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

    if((PlotName=="")==FALSE){#Check how much width is needed for the figure title/subtitle
      character_count <- nchar(PlotName)
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

    if(PlotName==""){
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
      plottable$heights[c(3)] <- unit(1,"cm")#controls margins --> PlotName and subtitle
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



