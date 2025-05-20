## ---------------------------
##
## Script name: Visualization  PCA plot
##
## Purpose of script: data Visualisation of the MetaProViz analysis to aid biological interpretation
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

#' PCA plot visualization
#'
#' @param data DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected. includes experimental design and outlier column.
#' @param metadata_info \emph{Optional: } NULL or Named vector including at least one of those three information : c(color="ColumnName_Plot_SettingsFile", shape= "ColumnName_Plot_SettingsFile"). \strong{Default = NULL}
#' @param metadata_sample \emph{Optional: } DF which contains information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames. Column "Conditions" with information about the sample conditions (e.g. "N" and "T" or "Normal" and "Tumor"), can be used for feature filtering and colour coding in the PCA. Column "AnalyticalReplicate" including numerical values, defines technical repetitions of measurements, which will be summarised. Column "BiologicalReplicates" including numerical values. Please use the following names: "Conditions", "Biological_Replicates", "Analytical_Replicates".\strong{Default = NULL}
#' @param color_palette \emph{Optional: } Provide customiced color-palette in vector format. For continuous scale use e.g. scale_color_gradient(low = "#88CCEE", high = "red") and for discrete scale c("#88CCEE",  "#DDCC77","#661100",  "#332288")\strong{Default = NULL}
#' @param scale_color \emph{Optional: } Either "continuous" or "discrete" colour scale. For numeric or integer you can choose either, for character you have to choose discrete. \strong{Default = NULL}
#' @param shape_palette \emph{Optional: } Provide customiced shape-palette in vector format. \strong{Default = NULL}
#' @param show_loadings \emph{Optional: } TRUE or FALSE for whether PCA loadings are also plotted on the PCA (biplot) \strong{Default = FALSE}
#' @param scaling \emph{Optional: } TRUE or FALSE for whether a data scaling is used \strong{Default = TRUE}
#' @param pcx \emph{Optional: } Numeric value of the PC that should be plotted on the x-axis \strong{Default = 1}
#' @param pcy \emph{Optional: } Numeric value of the PC that should be plotted on the y-axis \strong{Default = 2}
#' @param theme \emph{Optional: } Selection of theme for plot, e.g. theme_grey(). You can check for complete themes here: https://ggplot2.tidyverse.org/reference/ggtheme.html. If default=NULL we use theme_classic(). \strong{Default = "discrete"}
#' @param plot_name \emph{Optional: } String which is added to the output files of the PCA \strong{Default = ""}
#' @param save_plot \emph{Optional: } Select the file type of output plots. Options are svg, png, pdf or NULL. \strong{Default = svg}
#' @param print_plot \emph{Optional: } TRUE or FALSE, if TRUE Volcano plot is saved as an overview of the results. \strong{Default = TRUE}
#' @param path \emph{Optional:} Path to the folder the results should be saved at. \strong{default: NULL}
#'
#' @return List with two elements: Plot and Plot_Sized
#'
#' @examples
#' Intra <- intracell_raw[,-c(1:3)]%>%tibble::column_to_rownames("Code")
#' Res <- viz_pca(Intra)
#'
#' @keywords PCA
#'
#' @importFrom ggplot2 ggplot theme element_rect autoplot scale_shape_manual geom_hline geom_vline ggtitle
#  required for autoplot.prcomp:
#' @importFrom ggfortify fortify_map
#' @importFrom dplyr rename
#' @importFrom magrittr %>% %<>%
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom rlang !! :=
#' @importFrom logger log_info log_trace
#'
#' @export
#'
viz_pca <- function(data,
                   metadata_info= NULL,
                   metadata_sample = NULL,
                   color_palette= NULL,
                   scale_color="discrete",
                   shape_palette=NULL,
                   show_loadings = FALSE,
                   scaling = TRUE,
                   pcx=1,
                   pcy=2,
                   theme=NULL,#theme_classic()
                   plot_name= '',
                   save_plot = "svg",
                   print_plot=TRUE,
                   path = NULL
){

  ###########################################################################
  ## ------------ Create log file ----------- ##
  metaproviz_init()

  logger::log_info("viz_pca: PCA plot visualization")
  ## ------------ Check Input files ----------- ##
  # HelperFunction `check_param`
  check_param(data=data,
                          metadata_sample=metadata_sample,
                          metadata_feature=NULL,
                          metadata_info=metadata_info,
                          save_plot=save_plot,
                          save_table=NULL,
                          core=FALSE,
                          print_plot= print_plot)

  # check_param` Specific
  if(is.logical(show_loadings) == FALSE){
    message <- paste("The Show_Loadings value should be either =TRUE if loadings are to be shown on the PCA plot or = FALSE if not.")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }
  if(is.logical(scaling) == FALSE){
    message <- paste("The scaling value should be either =TRUE if data scaling is to be performed prior to the PCA or = FALSE if not.")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }

  if(any(is.na(data))) {
     data[is.na(data)] <- 0#replace NA with 0
     message <- paste("NA values are included in data that were set to 0 prior to performing PCA.")
     logger::log_info(message)
     message(message)
  }

  ## ------------ Create Results output folder ----------- ##
  folder <- NULL
  if(is.null(save_plot)==FALSE){
    folder <- save_path(folder_name= "PCAPlots",
                                    path=path)
  }
  logger::log_info("viz_pca results saved at ", folder)

  ###########################################################################
  ## ----------- Set the plot parameters: ------------ ##
  ##--- Prepare colour and shape palette
  if(is.null(color_palette)){
    if(scale_color=="discrete"){
      safe_colorblind_palette <- c("#88CCEE",  "#DDCC77","#661100",  "#332288", "#AA4499","#999933",  "#44AA99", "#882215",  "#6699CC", "#117733", "#888888","#CC6677", "black","gold1","darkorchid4","red","orange", "blue")
    }else if(scale_color=="continuous"){
      safe_colorblind_palette <- NULL
    }
  } else{
    safe_colorblind_palette <-color_palette
  }
  if(is.null(shape_palette)){
    safe_shape_palette <- c(15,17,16,18,6,7,8,11,12)
  } else{
    safe_shape_palette <-shape_palette
  }

  logger::log_info(paste("viz_pca colour:", paste(safe_colorblind_palette, collapse = ", ")))
  logger::log_info(paste("viz_pca shape:", paste(safe_shape_palette, collapse = ", ")))

  ##--- Prepare the color scheme:
  if("color" %in% names(metadata_info)==TRUE & "shape" %in% names(metadata_info)==TRUE){
    if((metadata_info[["shape"]] == metadata_info[["color"]])==TRUE){
      metadata_sample$shape <- metadata_sample[,paste(metadata_info[["color"]])]
      metadata_sample<- metadata_sample%>%
        dplyr::rename("color"=paste(metadata_info[["color"]]))
    }else{
      metadata_sample <- metadata_sample%>%
          dplyr::rename("color"=paste(metadata_info[["color"]]),
                        "shape"=paste(metadata_info[["shape"]]))
      }
 }else if("color" %in% names(metadata_info)==TRUE & "shape" %in% names(metadata_info)==FALSE){
   if("color" %in% names(metadata_info)==TRUE){
     metadata_sample <- metadata_sample%>%
        dplyr::rename("color"=paste(metadata_info[["color"]]))
   }
   if("shape" %in% names(metadata_info)==TRUE){
      metadata_sample <- metadata_sample%>%
        dplyr::rename("shape"=paste(metadata_info[["shape"]]))
   }
 }

  ##--- Prepare Input data:
  if(is.null(metadata_sample)==FALSE){
    InputPCA  <- merge(x=metadata_sample%>%tibble::rownames_to_column("UniqueID") , y=data%>%tibble::rownames_to_column("UniqueID"), by="UniqueID", all.y=TRUE)%>%
      tibble::column_to_rownames("UniqueID")
  }else{
    InputPCA  <- data
  }

 ##--- Prepare the color and shape settings:
  if("color" %in% names(metadata_sample)==TRUE){
    if(scale_color=="discrete"){
      InputPCA$color <- as.factor(InputPCA$color)
      color_select <- safe_colorblind_palette[1:length(unique(InputPCA$color))]
    }else if(scale_color=="continuous"){
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
        scale_color <- "discrete"
        logger::log_info("Warning: scale_color=continuous, but is.numeric or is.integer is FALSE, hence colour scale is set to discrete.")
        warning("scale_color=continuous, but is.numeric or is.integer is FALSE, hence colour scale is set to discrete.")
      }
    }
  }

  logger::log_info("viz_pca scale_color: ", scale_color)

  if("shape" %in% names(metadata_sample)==TRUE){
    shape_select <- safe_shape_palette[1:length(unique(InputPCA$shape))]

    if (!is.character(InputPCA$shape)) {
        # Convert the column to character
        InputPCA$shape <- as.character(InputPCA$shape)
    }
  }

  ##---  #assign column and legend name
  if("color" %in% names(metadata_sample)==TRUE){
    InputPCA  <- InputPCA%>%
      dplyr::rename(!!paste(metadata_info[["color"]]) :="color")
    Param_Col <-paste(metadata_info[["color"]])
  } else{
    color_select <- NULL
    Param_Col <- NULL
  }

  if("shape" %in% names(metadata_sample)==TRUE){
    InputPCA  <- InputPCA%>%
      dplyr::rename(!!paste(metadata_info[["shape"]]) :="shape")
    Param_Sha <-paste(metadata_info[["shape"]])
  } else{
    shape_select <-NULL
    Param_Sha <-NULL
  }

  ## ----------- Make the  plot based on the choosen parameters ------------ ##
  PlotList <- list()#Empty list to store all the plots
  PlotList_adaptedGrid <- list()

  #Make the plot:
  PCA <- ggplot2::autoplot(stats::prcomp(as.matrix(data), scale. = as.logical(scaling)),
                  data= InputPCA,
                  x= pcx ,
                  y= pcy,
                  colour = Param_Col,
                  fill =  Param_Col,
                  shape = Param_Sha,
                  size = 3,
                  alpha = 0.8,
                  label=T,
                  label.size=2.5,
                  label.repel = TRUE,
                  loadings= as.logical(show_loadings), #draws Eigenvectors
                  loadings.label = as.logical(show_loadings),
                  loadings.label.vjust = 1.2,
                  loadings.label.size=2.5,
                  loadings.colour="grey10",
                  loadings.label.colour="grey10" ) +
    ggplot2::scale_shape_manual(values=shape_select)+
    ggplot2::ggtitle(paste(plot_name)) +
    ggplot2::geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
    ggplot2::geom_vline(xintercept=0,  color = "black", linewidth=0.1)

    if(scale_color=="discrete"){
      PCA <-PCA + ggplot2::scale_color_manual(values=color_select)
    }else if(scale_color=="continuous" & is.null(color_palette)){
      PCA <-PCA + color_select
    }

  #Add the theme
  if(is.null(theme)==FALSE){
    PCA <- PCA+theme
  }else{
    PCA <- PCA+ggplot2::theme_classic()
  }

  ## Store the plot in the 'plots' list
  PlotList[["Plot"]] <- PCA

  #Set the total heights and widths
  PCA %<>% plot_grob_pca(metadata_info=metadata_info,plot_name=plot_name)
  plot_height <- grid::convertUnit(PCA$height, 'cm', valueOnly = TRUE)
  plot_width <- grid::convertUnit(PCA$width, 'cm', valueOnly = TRUE)
  PCA %<>%
    {ggplot2::ggplot() + annotation_custom(.)} %>%
    add(theme(panel.background = ggplot2::element_rect(fill = "transparent")))

  PlotList_adaptedGrid[["Plot_Sized"]] <- PCA

  ###########################################################################
  ##----- Save and Return
  #Here we make a list in which we will save the outputs:
  file_name <-plot_name %>% {`if`(nchar(.), sprintf('PCA_%s', .), 'PCA')}

  suppressMessages(suppressWarnings(
    save_res(
      inputlist_df=NULL,
      inputlist_plot= PlotList_adaptedGrid,
      save_table=NULL,
      save_plot=save_plot,
      path= folder,
      file_name= file_name,
      core=FALSE,
      print_plot=print_plot,
      plot_height=plot_height,
      plot_width=plot_width,
      plot_unit="cm")
  ))

  invisible(list(Plot = PlotList, Plot_Sized = PlotList_adaptedGrid))
}
