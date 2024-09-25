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


# REFACT: Each docstring should start with a title that fits a single line,
# i.e. <77 characters long

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
#' @param FolderPath \emph{Optional:} Path to the folder the results should be saved at. \strong{default: NULL}
# REFACT: use the @return directive to describe the return value of the function
#'
#' @keywords PCA
#' @importFrom ggplot2 ggplot theme element_rect
#' @importFrom dplyr rename
#' @importFrom magrittr %>% %<>%
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom rlang !! :=
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


  ###########################################################################
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

  # REFACT: What is the likelihood that these simple logical parameters have a
  # wrong type of value? Even if that happens, we don't even tell the user here
  # in which function the error happened, and which other function called it with
  # the wrong arguments.
  # CheckInput` Specific
  if(is.logical(ShowLoadings) == FALSE){
    stop("Check input. The Show_Loadings value should be either =TRUE if loadings are to be shown on the PCA plot or = FALSE if not.")
  }
  if(is.logical(Scaling) == FALSE){
    stop("Check input. The Scaling value should be either =TRUE if data scaling is to be performed prior to the PCA or = FALSE if not.")
  }

  if(any(is.na(InputData))==TRUE){
     InputData[is.na(InputData)] <- 0#replace NA with 0
     message("NA values are included in InputData that were set to 0 prior to performing PCA.")
  }


  ## ------------ Create Results output folder ----------- ##
  if(is.null(SaveAs_Plot)==FALSE){
    Folder <- MetaProViz:::SavePath(FolderName= "PCAPlots",
                                    FolderPath=FolderPath)
  }


  ###########################################################################
  ## ----------- Set the plot parameters: ------------ ##
  ##--- Prepare colour and shape palette
  if(is.null(ColorPalette)){
    if((ColorScale=="discrete")==TRUE){
      safe_colorblind_palette <- metaproviz_palette()
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
        safe_colorblind_palette <- metaproviz_palette()
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
  PlotList_adaptedGrid <- list()

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
    }else if(ColorScale=="continuous" & is.null(ColorPalette)){
      PCA <-PCA + color_select
    }

  #Add the theme
  if(is.null(Theme)==FALSE){
    PCA <- PCA+Theme
  }else{
    PCA <- PCA+theme_classic()
  }

  ## Store the plot in the 'plots' list
  PlotList[["Plot"]] <- PCA

  #Set the total heights and widths
  PCA %<>% PlotGrob_PCA(SettingsInfo=SettingsInfo, PlotName=PlotName)
  PlotHeight <- as.numeric(PCA$height)
  PlotWidth <- as.numeric(PCA$width)
  PCA %<>%
    {ggplot2::ggplot() + annotation_custom(.)} %>%
    add(theme(panel.background = element_rect(fill = "transparent")))

  PlotList_adaptedGrid[["Plot_Sized"]] <- PCA


  ###########################################################################
  ##----- Save and Return
  #Here we make a list in which we will save the outputs:
  FileName <- PlotName %>% {`if`(nchar(.), sprintf('PCA_%s', .), 'PCA')}

  suppressMessages(suppressWarnings(
    SaveRes(
      InputList_DF=NULL,
      InputList_Plot= PlotList_adaptedGrid,
      SaveAs_Table=NULL,
      SaveAs_Plot=SaveAs_Plot,
      FolderPath= Folder,
      FileName= FileName,
      CoRe=FALSE,
      PrintPlot=PrintPlot,
      PlotHeight=PlotHeight,
      PlotWidth=PlotWidth,
      PlotUnit="cm"
    )
  ))

  invisible(list(Plot = PlotList, Plot_Sized = PlotList_adaptedGrid))

}


##############################################################
### ### ### PCA helper function: Internal Function ### ### ###
##############################################################

#' @param InputPlot This is the ggplot object generated within the VizPCA function.
#' @param SettingsInfo Passed to VizPCA
#' @param PlotName Passed to VizPCA
#'
#' @keywords PCA helper function
#' @importFrom ggplot2 ggplotGrob
#' @importFrom magrittr %<>% %>% add multiply_by
#' @importFrom grid convertUnit unit
#' @importFrom ggpubr get_legend
#' @importFrom logger log_trace
#' @importFrom purrr map
#' @noRd
PlotGrob_PCA <- function(InputPlot, SettingsInfo, PlotName){

  #------- Set the total heights and widths
  #we need ggplot_grob to edit the gtable of the ggplot object. Using this we can manipulate the gtable arguments directly.
  plottable <-
    ggplotGrob(InputPlot) %>%
    withCanvasSize(width = 12, height = 11) # Convert the plot to a gtable
  ptb <<- plottable

  plottable %<>%
    ##############################################
    #-----widths general
    set_width("axis-b", "8cm") %>%
    set_width("ylab-l", "0cm", offset = -4L, ifempty = FALSE) %>%
    set_width("axis-l", "1cm") %>%
    set_width("ylab-l", "1cm") %>%
    set_width("guide-box-left", "0cm") %>%
    set_width("axis-r", "0cm") %>%
    set_width("ylab-r", "0cm") %>%
    set_width("ylab-l", "1cm", offset = -1L) %>%
    set_width("guide-box-right", "1cm") %>%
    #############################################
    #-----heigths general
    set_height("axis-l", "8cm") %>%
    set_height("axis-b", "1cm") %>%
    set_height("xlab-b", ".5cm") %>%
    set_height("xlab-b", "1cm", offset = 1L) %>%
    set_height("title", "0cm", offset = -2L, ifempty = FALSE) %>%
    set_height("title", "0cm", offset = -1L) %>%
    set_height("title", "1cm") %>%
    set_height("subtitle", "0cm") %>%
    set_height("guide-box-top", "0cm") %>%
    set_height("xlab-t", "0cm", offset = -1L)

  #############################################
  #----- Plot Name
    #------ Height
  if (nchar(PlotName) > 0L) {
    plottable %<>% set_height("title", "1.5cm")#controls margins --> PlotName and subtitle
    # Sum up total heights:
    plottable$height %<>% add(unit(.5, "cm"))

    #------- Width: Check how much width is needed for the figure title/subtitle
    title_width <- unit(nchar(PlotName) * .25 + .8, "cm")
    plottable %<>% set_width(
      "guide-box-right",
      sprintf('%.02fcm', title_width - plottable$width),
      callback = max
    )
    plottable$width %<>% max(title_width)

  }

  #############################################
  #-------- Figure legend
  legend_sections <- c("color", "shape")
  if(any(legend_sections %in% names(SettingsInfo))) {

    log_trace('color and/or shape in SettingsInfo')
    Legend <- get_legend(InputPlot) # Extract legend to adjust separately
    leg <<- Legend
    #gtable::gtable_show_layout(leg)
    #plot(leg)

    #------- Legend widths
    ## Legend titles:
    legend_nchar <-
       legend_sections %>%
       map(nchar) %>%
       unlist %>%
       max %>%
       multiply_by(.25)

    legend_width <-
      Legend$width[3L] %>%
      as.numeric %>%
      round(1L) %>%
      max(legend_nchar)

    ## Legend space:
    plottable %<>% set_width(
      "guide-box-right",
      sprintf('%.02fcm', legend_width),
      callback = max,
      grow = TRUE
    )

    #------- Legend heights
    Legend_heights <- (round(as.numeric(Legend$heights[3L]), 1L))+(round(as.numeric(Legend$heights[5L]), 1L)) + 2 #+2 to secure space above and below plot
    if(as.numeric(plottable$height) < Legend_heights){
      Add <- (Legend_heights-plottable$height) / 2
      plottable$heights[1] <- unit(Add,"cm")#controls margins --> Can be increased if Figure legend needs more space on the top
      plottable$heights[12] <- unit(Add,"cm")#controls margins --> Can be increased if Figure legend needs more space on the bottom
      plottable$height <- unit(Legend_heights, "cm")
    }

  }else{
    log_trace('No SettingsInfo')# = No figure legend
  }
  ptb1 <<- plottable

  plottable
}



