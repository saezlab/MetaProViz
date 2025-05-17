## ---------------------------
##
## Script name: Visualization
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
#' This script allows you to perform different data visualizations using the results of the MetaProViz analysis

#####################################
### ### ### Volcano Plots ### ### ###
#####################################

#' Volcano plot visualization
#'
#' @param plot_types \emph{Optional: } Choose between "Standard" (data), "Compare" (plot two comparisons together data and data2) or "PEA" (Pathway Enrichment Analysis) \strong{Default = "Standard"}
#' @param metadata_info \emph{Optional: } NULL or Named vector including at least one of those three information for Settings="Standard" or "Compare": c(color ="ColumnName_metadata_feature", shape = "ColumnName_metadata_feature", individual="ColumnName_metadata_feature"). For Settings="PEA" a named vector with: PEA_Pathway="ColumnName_data2"=each pathway will be plotted, PEA_score="ColumnName_data2", PEA_stat= "ColumnName_data2"= usually p.adj column, "PEA_Feature="ColumnName_data2"= usually Metabolites), optionally you can additionally include c(color_Metab="ColumnName_metadata_feature", shape= "ColumnName_metadata_feature").\strong{Default = NULL}
#' @param metadata_feature \emph{Optional: } DF with column including the Metabolite names (needs to match Metabolite names and Metabolite column name of data) and other columns with required plot_typeInfo. \strong{Default = NULL}
#' @param data DF with metabolites as row names and columns including Log2FC and stat (p-value, p.adjusted) value columns.
#' @param data2 \emph{Optional: } DF to compare to main Input_data with the same column names x and y (Settings="Compare") and metabolites as row names or Pathway enrichment analysis results (Settings="PEA"). \strong{Default = NULL}
#' @param y \emph{Optional: } Column name including the values that should be used for y-axis. Usually this would include the p.adjusted value. \strong{Default = "p.adj"}
#' @param x \emph{Optional: } Column name including the values that should be used for x-axis. Usually this would include the Log2FC value. \strong{Default = "Log2FC"}
#' @param plot_name \emph{Optional: } String which is added to the output files of the plot. \strong{Default = ""}
#' @param name_comparison \emph{Optional: } Named vector including those information about the two datasets that are compared on the plots when choosing Settings= "Compare". \strong{Default = c(data="Cond1", data2= "Cond2")}
#' @param xlab \emph{Optional: } String to replace x-axis label in plot. \strong{Default = NULL}
#' @param ylab \emph{Optional: } String to replace y-axis label in plot. \strong{Default = NULL}
#' @param cutoff_x \emph{Optional: } Number of the desired log fold change cutoff for assessing significance. \strong{Default = 0.5}
#' @param cutoff_y \emph{Optional: } Number of the desired p value cutoff for assessing significance. \strong{Default = 0.05}
#' @param color_palette \emph{Optional: } Provide customiced color-palette in vector format. \strong{Default = NULL}
#' @param shape_palette \emph{Optional: } Provide customiced shape-palette in vector format. \strong{Default = NULL}
#' @param select_label \emph{Optional: } If set to NULL, feature labels will be plotted randomly. If vector is provided, e.g. c("MetaboliteName1", "MetaboliteName2"), selected names will be plotted. If set to default "", no feature names will be plotted. \strong{Default = ""}
#' @param connectors \emph{Optional: } TRUE or FALSE for whether connectors from names to points are to be added to the plot. \strong{Default =  FALSE}
#' @param subtitle \emph{Optional: } \strong{Default = ""}
#' @param theme \emph{Optional: } Selection of theme for plot, e.g. theme_grey(). You can check for complete themes here: https://ggplot2.tidyverse.org/reference/ggtheme.html. \strong{Default = NULL}
#' @param path {Optional:} Path to the folder the results should be saved at. \strong{default: NULL}
#' @param feature \emph{Optional: } Name of the feature that are plotted, e.g. "Metabolites", "RNA", "Proteins", "Genes", etc. \strong{Default = "metabolites"}
#' @param save_plot \emph{Optional: } Select the file type of output plots. Options are svg, pdf, png or NULL. \strong{Default = "svg"}
#' @param print_plot \emph{Optional: } print the plots to the active graphic
#' device.
#'
#' @return List with two elements: Plot and Plot_Sized
#'
#' @examples
#' Intra <- MetaProViz::ToyData("IntraCells_DMA")
#' Res <- MetaProViz::viz_volcano(data=Intra)
#'
#' @keywords Volcano plot, pathways
#'
#' @importFrom ggplot2 ggplot theme
#' @importFrom dplyr rename filter mutate rename_with
#' @importFrom magrittr %>% %<>%
#' @importFrom tibble rownames_to_column column_to_rownames remove_rownames
#' @importFrom logger log_trace
#' @importFrom tidyselect all_of
#'
#' @export
#'
viz_volcano <- function(plot_types="Standard",
                       data,
                       metadata_info= NULL,
                       metadata_feature=NULL,
                       data2= NULL,
                       y= "p.adj",
                       x= "Log2FC",
                       xlab= NULL,#"~Log[2]~FC"
                       ylab= NULL,#"~-Log[10]~p.adj"
                       cutoff_x= 0.5,
                       cutoff_y= 0.05,
                       connectors=  FALSE,
                       select_label= "",
                       plot_name= "",
                       subtitle= "",
                       name_comparison= c(data="Cond1", data2= "Cond2"),
                       color_palette= NULL,
                       shape_palette=NULL,
                       theme= NULL,
                       save_plot= "svg",
                       path = NULL,
                       feature="Metabolites",
                       print_plot=TRUE){
  ## ------------ Create log file ----------- ##
  metaproviz_init()

  ## ------------ Check Input files ----------- ##
  # HelperFunction `check_param`
  if(plot_types=="PEA"){
    #Those relationships are checked in the viz_volcano_pea() function!
    SettingsFile <- NULL # For PEA the metadata_feature is the prior knowledge file, and hence this will not have feature as row names.
    Info <- NULL # If SettingsFileMetab=NULL, SetingsInfo has to be NULL to, otherwise we will get an error.
  }else{
    SettingsFile <-metadata_feature
    Info <- metadata_info
  }

  check_param(data=as.data.frame(t(data)),
                          data_num=FALSE,
                          metadata_sample=NULL,
                          metadata_feature=SettingsFile,#Set above
                          metadata_info=Info,#Set above
                          save_plot=save_plot,
                          save_table=NULL,
                          core=FALSE,
                          print_plot= print_plot,
                          plot_types="Feature")

  # check_param` Specific:
  if(is.numeric(cutoff_y)== FALSE |cutoff_y > 1 | cutoff_y < 0){
    message<- paste0("Check input. The selected cutoff_y value should be numeric and between 0 and 1.")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
    }
  if(is.numeric(cutoff_x)== FALSE  | cutoff_x < 0){
    message<- paste0("Check input. The selected cutoff_x value should be numeric and between 0 and +oo.")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
    }
  if(paste(x) %in% colnames(data)==FALSE | paste(y) %in% colnames(data)==FALSE){
    message<- paste0("Check your input. The column name of x and/ore y does not exist in Input_data.")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }

  if(is.null(select_label)==FALSE & is.vector(select_label)==FALSE){
    message<- paste0("Check input. select_label must be either NULL or a vector.")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }

  if(is.logical(connectors) == FALSE){
    message<- paste0("Check input. The connectors value should be either = TRUE if connectors from names to points are to be added to the plot or =FALSE if not.")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }

  if(is.null(PlotName)==FALSE & is.vector(PlotName)==FALSE){
    message<- paste0("Check input.plot_name must be either NULL or a vector.")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }

  Plot_options <- c("Standard", "Compare", "PEA")
  if (plot_types %in% Plot_options == FALSE){
    message<- paste0("plot_types option is incorrect. The allowed options are the following: ",paste(Plot_options, collapse = ", "),"." )
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }

  ## ------------ Create Results output folder ----------- ##
  if(is.null(save_plot)==FALSE){
    folder <- save_path(folder_name= "VolcanoPlots",
                                    path=path)
  }

  ############################################################################################################
  ## ----------- Prepare data ------------ ##
  #Extract required columns and merge with SettingsFile
   if(is.null(metadata_feature)==FALSE){
     ##--- Prepare the color scheme:
     if("color" %in% names(metadata_info)==TRUE & "shape" %in% names(metadata_info)==TRUE){
       if((metadata_info[["shape"]] == metadata_info[["color"]])==TRUE){
         metadata_feature$shape <- metadata_feature[,paste(metadata_info[["color"]])]
         metadata_feature<- metadata_feature%>%
           dplyr::rename("color"=paste(metadata_info[["color"]]))
       }else{
         metadata_feature <- metadata_feature%>%
           dplyr::rename("color"=paste(metadata_info[["color"]]),
                         "shape"=paste(metadata_info[["shape"]]))
       }
     }else if("color" %in% names(metadata_info)==TRUE & "shape" %in% names(metadata_info)==FALSE){
       metadata_feature <- metadata_feature%>%
         dplyr::rename("color"=paste(metadata_info[["color"]]))
     }else if("color" %in% names(metadata_info)==FALSE & "shape" %in% names(metadata_info)==TRUE){
       metadata_feature <- metadata_feature%>%
         dplyr::rename("shape"=paste(metadata_info[["shape"]]))
     }
     if("individual" %in% names(metadata_info)==TRUE){
       metadata_feature <- metadata_feature%>%
         dplyr::rename("individual"=paste(metadata_info[["individual"]]))
     }


    ##--- Merge data with SettingsFile:
     common_columns <- character(0)  # Initialize an empty character vector
     for(col_name in colnames(data[, c(x, y)])) {
       if(col_name %in% colnames(metadata_feature)) {
         common_columns <- c(common_columns, col_name)  # Add the common column name to the vector
       }
     }

     if (length(common_columns)) {

       metadata_feature %<>%  # rename those column since they otherwise
                                   # will cause issues when we merge the DFs
                                   # later
        # this should not be handled like this, use suffixes for dplyr::join
        # instead
        dplyr::rename_with(
          ~paste0(.x, "_metadata_feature"),
          tidyselect::all_of(common_columns)
        )
    }

    if(plot_types=="PEA"){
      Volcanodata <- merge(x=metadata_feature ,y=data[, c(x, y)], by.x=metadata_info[["PEA_Feature"]] , by.y=0, all.y=TRUE)%>%
        tibble::remove_rownames()%>%
        dplyr::mutate(FeatureNames = metadata_info[["PEA_Feature"]])%>%
        dplyr::filter(!is.na(x) | !is.na(x))
    }else{
     Volcanodata <- merge(x=metadata_feature ,y=data[, c(x, y)], by=0, all.y=TRUE)%>%
       tibble::remove_rownames()%>%
       tibble::column_to_rownames("Row.names")%>%
       dplyr::mutate(FeatureNames = rownames(data))%>%
       dplyr::filter(!is.na(x) | !is.na(x))
    }

   }else{
     Volcanodata <- data[, c(x, y)]%>%
       dplyr::mutate(FeatureNames = rownames(data))%>%
       dplyr::filter(!is.na(x) | !is.na(x))
  }

  # Rename the x and y lab if the information has been passed:
  if(is.null(xlab)==TRUE){#use column name of x provided by user
    xlab <- bquote(.(as.symbol(x)))
    }else if(is.null(xlab)==FALSE){
    xlab <- bquote(.(as.symbol(xlab)))
    }

  if(is.null(ylab)==TRUE){#use column name of x provided by user
    ylab <- bquote(.(as.symbol(y)))
    }else if(is.null(ylab)==FALSE){
      ylab <- bquote(.(as.symbol(ylab)))
      }

  ## ----------- Set the plot parameters: ------------ ##
  ##--- Prepare colour and shape palette
  if(is.null(color_palette)){
    if("color" %in% names(metadata_info)==TRUE){
      safe_colorblind_palette <- c("#88CCEE",  "#DDCC77","#661100",  "#332288", "#AA4499","#999933",  "#44AA99", "#882215",  "#6699CC", "#117733", "#888888","#CC6677", "black","gold1","darkorchid4","red","orange", "blue")
    }else{
      safe_colorblind_palette <- c("#888888", "#44AA99", "#44AA99","#CC6677")
    }

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

  ############################################################################################################
  ## ----------- Make the  plot based on the chosen parameters ------------ ##

  if(plot_types=="Standard"){#####--- 1. Standard
    VolcanoRes <- viz_volcano_standard(data= Volcanodata,
                                                   metadata_feature=metadata_feature,
                                                   metadata_info=metadata_info,
                                                   y= y,
                                                   x= x,
                                                   xlab= xlab,
                                                   ylab= ylab,
                                                   cutoff_x= cutoff_x,
                                                   cutoff_y= cutoff_y,
                                                   connectors= connectors,
                                                   select_label=select_label,
                                                  plot_name=plot_name,
                                                   subtitle= subtitle,
                                                   color_palette=safe_colorblind_palette,
                                                   shape_palette=safe_shape_palette,
                                                   theme= theme,
                                                   feature=feature,
                                                   save_plot=save_plot,
                                                   print_plot=print_plot,
                                                   folder=folder)

  }else if(plot_types=="Compare"){#####--- 2. Compare
    VolcanoRes <- viz_volcano_compare(data= Volcanodata,
                                                  data2=data2,
                                                  metadata_feature=metadata_feature,
                                                  metadata_info=metadata_info,
                                                  y= y,
                                                  x= x,
                                                  xlab= xlab,
                                                  ylab= ylab,
                                                  cutoff_x= cutoff_x,
                                                  cutoff_y= cutoff_y,
                                                  connectors= connectors,
                                                  select_label=select_label,
                                                 plot_name=plot_name,
                                                  subtitle= subtitle,
                                                  color_palette=safe_colorblind_palette,
                                                  shape_palette=safe_shape_palette,
                                                  theme= theme,
                                                  feature=feature,
                                                  name_comparison=name_comparison,
                                                  save_plot=save_plot,
                                                  print_plot=print_plot,
                                                  folder=folder)

  } else if(plot_types=="PEA"){#####--- 3. PEA
    VolcanoRes <- viz_volcano_pea(data= Volcanodata,
                                              data2=data2,
                                              metadata_feature=metadata_feature,#Problem: we need to know the column name of the feature!
                                              metadata_info=metadata_info,
                                              y= y,
                                              x= x,
                                              xlab= xlab,
                                              ylab= ylab,
                                              cutoff_x= cutoff_x,
                                              cutoff_y= cutoff_y,
                                              connectors= connectors,
                                              select_label=select_label,
                                             plot_name=plot_name,
                                              subtitle= subtitle,
                                              color_palette=safe_colorblind_palette,
                                              shape_palette=safe_shape_palette,
                                              theme= theme,
                                              feature=feature,
                                              save_plot=save_plot,
                                              print_plot=print_plot,
                                              folder=folder)
  }
  return(invisible(VolcanoRes))
}

################################################################################################
### ### ### viz_volcano helper function: Internal Function for plot_types Standard ### ### ###
################################################################################################

#' viz_volcano_standard
#'
#' @param data Passed to main function viz_volcano()
#' @param metadata_feature Passed to main function viz_volcano()
#' @param metadata_info Passed to main function viz_volcano()
#' @param y \emph{Optional: } Passed to main function viz_volcano() \strong{Default = "p.adj"}
#' @param x \emph{Optional: } Passed to main function viz_volcano() \strong{Default = "Log2FC"}
#' @param plot_name \emph{Optional: } Passed to main function viz_volcano() \strong{Default = ""}
#' @param xlab \emph{Optional: } Passed to main function viz_volcano()  \strong{Default = NULL}
#' @param ylab \emph{Optional: } Passed to main function viz_volcano() \strong{Default = NULL}
#' @param cutoff_x \emph{Optional: } Passed to main function viz_volcano() \strong{Default = 0.5}
#' @param cutoff_y \emph{Optional: } Passed to main function viz_volcano() \strong{Default = 0.05}
#' @param select_label \emph{Optional: } Passed to main function viz_volcano() \strong{Default = ""}
#' @param connectors \emph{Optional: } Passed to main function viz_volcano() \strong{Default =  FALSE}
#' @param subtitle \emph{Optional: } Passed to main function viz_volcano() \strong{Default = ""}
#' @param color_palette Created in viz_volcano() based on color_palette passed to main function viz_volcano()
#' @param shape_palette Created in viz_volcano() based on shape_palette passed to main function viz_volcano()
#' @param theme \emph{Optional: } Selection of theme for plot, e.g. theme_grey(). You can check for complete themes here: https://ggplot2.tidyverse.org/reference/ggtheme.html. \strong{Default = NULL}
#' @param feature \emph{Optional: } Name of the feature that are plotted, e.g. "Metabolites", "RNA", "Proteins", "Genes", etc. \strong{Default = "Metabolites"}
#' @param save_plot Passed to main function viz_volcano()
#' @param print_plot Passed to main function viz_volcano()
#' @param folder Created in viz_volcano(). Path to the folder where files are saved.
#'
#' @return List with two elements: Plot and Plot_Sized
#'
#' @keywords Standard volcano plots
#'
#' @importFrom ggplot2 ggplot theme
#' @importFrom dplyr rename filter mutate
#' @importFrom magrittr %>% %<>%
#' @importFrom tibble rownames_to_column column_to_rownames remove_rownames
#'
#' @noRd
#'
viz_volcano_standard <- function(data,
                                metadata_feature,
                                metadata_info,
                                y= "p.adj",
                                x= "Log2FC",
                                xlab= NULL,#"~Log[2]~FC"
                                ylab= NULL,#"~-Log[10]~p.adj"
                                cutoff_x= 0.5,
                                cutoff_y= 0.05,
                                connectors=  FALSE,
                                select_label= "",
                               plot_name= "",
                                subtitle= "",
                                color_palette,
                                shape_palette,
                                theme= NULL,
                                feature="Metabolites",
                                save_plot,
                                print_plot,
                                folder){

  #Pass colours/shapes
  safe_colorblind_palette <- color_palette
  safe_shape_palette <- shape_palette

  #Plots
  if("individual" %in% names(metadata_info)==TRUE){
    # Create the list of individual plots that should be made:
    IndividualPlots <- unique(data$individual)

    PlotList <- list()#Empty list to store all the plots
    PlotList_adaptedGrid <- list()#Empty list to store all the plots

    for (i in IndividualPlots){
      InputVolcano <- subset(data, individual == paste(i))

      if(nrow(InputVolcano)>=1){
        if("color" %in% names(metadata_info)==TRUE ){
          color_select <- safe_colorblind_palette[1:length(unique(InputVolcano$color))]

          keyvals <- c()
          for(row in 1:nrow(InputVolcano)){
            col <- color_select[unique(InputVolcano$color) %in% InputVolcano[row, "color"]]
            names(col) <- InputVolcano$color[row]
            keyvals <- c(keyvals, col)
          }

          LegendPos<- "right"
        } else{
          keyvals <-NULL
        }
        #Prepare the shape scheme:
        if("shape" %in% names(metadata_info)==TRUE){
          shape_select <- safe_shape_palette[1:length(unique(InputVolcano$shape))]

          keyvalsshape <- c()
          for(row in 1:nrow(InputVolcano)){
            sha <- shape_select[unique(InputVolcano$shape) %in% InputVolcano[row, "shape"]]
            names(sha) <- InputVolcano$shape[row]
            keyvalsshape <- c(keyvalsshape, sha)
          }

          LegendPos<- "right"
        } else{
          keyvalsshape <-NULL
        }

        if("color" %in% names(metadata_info)==FALSE & "shape" %in% names(metadata_info)==FALSE){
          LegendPos<- "none"
        }

        #Prepare the Plot:
        Plot<- EnhancedVolcano::EnhancedVolcano(InputVolcano,
                                                lab = InputVolcano$FeatureNames,#Metabolite name
                                                select_label = select_label,
                                                x = paste(x),
                                                y = paste(y),
                                                xlab  =xlab,
                                                ylab =ylab,
                                                cutoff_stat = cutoff_y,
                                                FCcutoff = cutoff_x,#Cut off Log2FC, automatically 2
                                                pointSize = 3,
                                                labSize = 3,
                                                axisLabSize = 10,
                                                titleLabSize = 12,
                                                subtitleLabSize = 11,
                                                captionLabSize = 10,
                                                col=safe_colorblind_palette,
                                                colCustom = keyvals,
                                                shapeCustom = keyvalsshape,
                                                colAlpha = 1,
                                                title= paste(PlotName, ": ", i, sep=""),
                                                subtitle = subtitle,
                                                caption = paste0("total = ", nrow(InputVolcano), " ", feature),
                                                xlim =  c(min(InputVolcano[[x]][is.finite(InputVolcano[[x]] )])-0.2, max(InputVolcano[[x]][is.finite(InputVolcano[[x]])])+1.2),
                                                ylim = c(0,(ceiling(-log10(Reduce(min,InputVolcano[[y]]))))),
                                                cutoffLineType = "dashed",
                                                cutoffLineCol = "black",
                                                cutoffLineWidth = 0.5,
                                                legendLabels=c(paste(x," < |", cutoff_x, "|"), paste(x," > |", cutoff_x, "|"), paste(y, ' < ', cutoff_y) , paste(y, ' < ', cutoff_y,' & ',x," < |", cutoff_x, "|")),
                                                legendPosition = LegendPos,
                                                legendLabSize = 7,
                                                legendIconSize =4,
                                                gridlines.major = FALSE,
                                                gridlines.minor = FALSE,
                                                drawConnectors = connectors)
        #Add the theme
        if(is.null(theme)==FALSE){
          Plot <- Plot+theme
        }

        ## Store the plot in the 'plots' list
        PlotList[[i]] <- Plot

        #Set the total heights and widths
        PlotTitle <- paste(PlotName, ": ", i, sep="")
        Plot_Sized <-  plot_grob_volcano(input_plot=Plot, metadata_info=metadata_info, plot_name = PlotTitle, subtitle = subtitle)
        plot_height <- grid::convertUnit(Plot_Sized$height, 'cm', valueOnly = TRUE)
        plot_width <- grid::convertUnit(Plot_Sized$width, 'cm', valueOnly = TRUE)
        Plot_Sized %<>%
          {ggplot2::ggplot() + annotation_custom(.)} %>%
          add(theme(panel.background = ggplot2::element_rect(fill = "transparent")))

        cleaned_i <- gsub("[[:space:],/\\\\]", "-", i)#removes empty spaces and replaces /,\ with -
        PlotList_adaptedGrid[[cleaned_i]] <- Plot_Sized

        SaveList <- list()
        SaveList[[cleaned_i]] <- Plot_Sized

        #----- Save
        suppressMessages(suppressWarnings(
          save_res(inputlist_df=NULL,
                               inputlist_plot= SaveList,
                               save_table=NULL,
                               save_plot=save_plot,
                               path= folder,
                               file_name= paste("Volcano_",PlotName, sep=""),
                               core=FALSE,
                               print_plot=print_plot,
                               plot_height= plot_height,
                               plot_width=plot_width,
                               plot_unit="cm")))
      }
    }
  }else if("individual" %in% names(metadata_info)==FALSE){
    PlotList <- list()#Empty list to store all the plots
    PlotList_adaptedGrid <- list()#Empty list to store all the plots

    InputVolcano <- data
    if(nrow(InputVolcano)>=1){
      if("color" %in% names(metadata_info)==TRUE ){
        color_select <- safe_colorblind_palette[1:length(unique(InputVolcano$color))]

        keyvals <- c()
        for(row in 1:nrow(InputVolcano)){
          col <- color_select[unique(InputVolcano$color) %in% InputVolcano[row, "color"]]
          names(col) <- InputVolcano$color[row]
          keyvals <- c(keyvals, col)
        }

        LegendPos<- "right"
      } else{
        keyvals <-NULL
      }
      #Prepare the shape scheme:
      if("shape" %in% names(metadata_info)==TRUE){
        shape_select <- safe_shape_palette[1:length(unique(InputVolcano$shape))]

        keyvalsshape <- c()
        for(row in 1:nrow(InputVolcano)){
          sha <- shape_select[unique(InputVolcano$shape) %in% InputVolcano[row, "shape"]]
          names(sha) <- InputVolcano$shape[row]
          keyvalsshape <- c(keyvalsshape, sha)
        }

        LegendPos<- "right"
      } else{
        keyvalsshape <-NULL
      }

      if("color" %in% names(metadata_info)==FALSE & "shape" %in% names(metadata_info)==FALSE){
        LegendPos<- "none"
      }

      #Prepare the Plot:
      Plot<- EnhancedVolcano::EnhancedVolcano(InputVolcano,
                                              lab = InputVolcano$FeatureNames,#Metabolite name
                                              select_label = select_label,
                                              x = paste(x),
                                              y = paste(y),
                                              xlab  =xlab,
                                              ylab =ylab,
                                              cutoff_stat = cutoff_y,
                                              FCcutoff = cutoff_x,#Cut off Log2FC, automatically 2
                                              pointSize = 3,
                                              labSize = 3,
                                              axisLabSize = 10,
                                              titleLabSize = 12,
                                              subtitleLabSize = 11,
                                              captionLabSize = 10,
                                              col=safe_colorblind_palette,
                                              colCustom = keyvals,
                                              shapeCustom = keyvalsshape,
                                              colAlpha = 1,
                                              title= paste(PlotName),
                                              subtitle = subtitle,
                                              caption = paste0("total = ", nrow(InputVolcano), " ", feature),
                                              xlim =  c(min(InputVolcano[[x]][is.finite(InputVolcano[[x]] )])-0.2, max(InputVolcano[[x]][is.finite(InputVolcano[[x]])])+1.2),
                                              ylim = c(0,(ceiling(-log10(Reduce(min,InputVolcano[[y]]))))),
                                              cutoffLineType = "dashed",
                                              cutoffLineCol = "black",
                                              cutoffLineWidth = 0.5,
                                              legendLabels=c(paste(x," < |", cutoff_x, "|"), paste(x," > |", cutoff_x, "|"), paste(y, ' < ', cutoff_y) , paste(y, ' < ', cutoff_y,' & ',x," < |", cutoff_x, "|")),
                                              legendPosition = LegendPos,
                                              legendLabSize = 9,
                                              legendIconSize =4,
                                              gridlines.major = FALSE,
                                              gridlines.minor = FALSE,
                                              drawConnectors = connectors)
      #Add the theme
      if(is.null(theme)==FALSE){
        Plot <- Plot+theme
      }

      ## Store the plot in the 'plots' list
      PlotList[["Plot"]] <- Plot

      #Set the total heights and widths
      Plot_Sized <-  plot_grob_volcano(input_plot=Plot, metadata_info=metadata_info, plot_name =plot_name, subtitle = subtitle)
      plot_height <- grid::convertUnit(Plot_Sized$height, 'cm', valueOnly = TRUE)
      plot_width <- grid::convertUnit(Plot_Sized$width, 'cm', valueOnly = TRUE)
      Plot_Sized %<>%
        {ggplot2::ggplot() + annotation_custom(.)} %>%
        add(theme(panel.background = ggplot2::element_rect(fill = "transparent")))

      PlotList_adaptedGrid[["Plot_Sized"]] <- Plot_Sized

      #----- Save
      suppressMessages(suppressWarnings(
      save_res(inputlist_df=NULL,
                             inputlist_plot= list("Plot_Sized"= PlotList_adaptedGrid[["Plot_Sized"]]),
                             save_table=NULL,
                             save_plot=save_plot,
                             path= folder,
                             file_name= paste("Volcano_",plot_name, sep=""),
                             core=FALSE,
                             print_plot=print_plot,
                             plot_height=plot_height,
                             plot_width=plot_width,
                             plot_unit="cm")))
    }
  }
  return(invisible(list("Plot"=PlotList,"Plot_Sized" = PlotList_adaptedGrid)))
}




################################################################################################
### ### ### viz_volcano helper function: Internal Function for plot_types Compare ### ### ###
################################################################################################

#' Check input parameters
#'
#' @param data Passed to main function viz_volcano()
#' @param data2 Passed to main function viz_volcano()
#' @param metadata_feature Passed to main function viz_volcano()
#' @param metadata_info Passed to main function viz_volcano()
#' @param y \emph{Optional: } Passed to main function viz_volcano() \strong{Default = "p.adj"}
#' @param x \emph{Optional: } Passed to main function viz_volcano() \strong{Default = "Log2FC"}
#' @param plot_name \emph{Optional: } Passed to main function viz_volcano() \strong{Default = ""}
#' @param xlab \emph{Optional: } Passed to main function viz_volcano()  \strong{Default = NULL}
#' @param ylab \emph{Optional: } Passed to main function viz_volcano() \strong{Default = NULL}
#' @param cutoff_x \emph{Optional: } Passed to main function viz_volcano() \strong{Default = 0.5}
#' @param cutoff_y \emph{Optional: } Passed to main function viz_volcano() \strong{Default = 0.05}
#' @param select_label \emph{Optional: } Passed to main function viz_volcano() \strong{Default = ""}
#' @param connectors \emph{Optional: } Passed to main function viz_volcano() \strong{Default =  FALSE}
#' @param subtitle \emph{Optional: } Passed to main function viz_volcano() \strong{Default = ""}
#' @param color_palette Created in viz_volcano() based on color_palette passed to main function viz_volcano()
#' @param shape_palette Created in viz_volcano() based on shape_palette passed to main function viz_volcano()
#' @param theme \emph{Optional: } Selection of theme for plot, e.g. theme_grey(). You can check for complete themes here: https://ggplot2.tidyverse.org/reference/ggtheme.html. \strong{Default = NULL}
#' @param feature \emph{Optional: } Name of the feature that are plotted, e.g. "Metabolites", "RNA", "Proteins", "Genes", etc. \strong{Default = "Metabolites"}
#' @param name_comparison Passed to main function viz_volcano()
#' @param save_plot Passed to main function viz_volcano()
#' @param print_plot Passed to main function viz_volcano()
#' @param folder Created in viz_volcano(). Path to the folder where files are saved.
#'
#' @return List with two elements: Plot and Plot_Sized
#'
#' @keywords Compare volcano plots
#'
#' @importFrom ggplot2 ggplot theme
#' @importFrom dplyr rename filter mutate
#' @importFrom magrittr %>% %<>%
#' @importFrom tibble rownames_to_column column_to_rownames remove_rownames
#' @importFrom logger log_trace
#'
#' @noRd
#'
viz_volcano_compare <- function(data,
                               data2,
                               metadata_feature,
                               metadata_info,
                               y= "p.adj",
                               x= "Log2FC",
                               xlab= NULL,#"~Log[2]~FC"
                               ylab= NULL,#"~-Log[10]~p.adj"
                               cutoff_x= 0.5,
                               cutoff_y= 0.05,
                               connectors=  FALSE,
                               select_label= "",
                              plot_name= "",
                               subtitle= "",
                               color_palette,
                               shape_palette,
                               theme= NULL,
                               feature="Metabolites",
                               name_comparison,
                               save_plot,
                               print_plot,
                               folder){

  #####################
  ##--- Check data
  if(is.data.frame(data2)==FALSE){
    if(paste(x) %in% colnames(data2)==FALSE | paste(y) %in% colnames(data2)==FALSE){
      message <- paste("Check your data2. The column name of ", x, " and/or ", y, " does not exist in data2.")
      logger::log_trace(paste("Error ", message, sep=""))
      stop(message)
    }
    }

  if(any(duplicated(row.names(data2)))==TRUE){
    message <- paste("Duplicated row.names of data2, whilst row.names must be unique")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }

  #Pass colours/shapes
  safe_colorblind_palette <- color_palette
  safe_shape_palette <- shape_palette

  ##--- Prepare Input data
  if(is.null(metadata_feature)==FALSE){
  data2 <- merge(x=metadata_feature%>%tibble::rownames_to_column("FeatureNames") , y=data2[, c(x, y)]%>%tibble::rownames_to_column("FeatureNames") , by="FeatureNames", all.y=TRUE)%>%
    filter(!is.na(x) | !is.na(x))
  data[,"comparison"]  <- as.character(paste(name_comparison[["data"]]))
  data2[,"comparison"]  <- as.character(paste(name_comparison[["data2"]]))
  InputCompare  <- rbind(data,data2)

  }else{
   data2 <-  data2[, c(x, y)]%>%
    mutate(FeatureNames = rownames(data2))%>%
    na.omit()

   #Combine DFs and add appropriate column names
   data[,"comparison"]  <- as.character(paste(name_comparison[["data"]]))
   data2[,"comparison"]  <- as.character(paste(name_comparison[["data2"]]))
   InputCompare  <- rbind(data[,c("FeatureNames", x, y, "comparison")],data2[,c("FeatureNames", x, y, "comparison")])
  }



  #####################
  ##--- Plots
  if("individual" %in% names(metadata_info)==TRUE){
    # Create the list of individual plots that should be made:
    IndividualPlots <- unique(InputCompare$individual)

    PlotList <- list()#Empty list to store all the plots
    PlotList_adaptedGrid <- list()#Empty list to store all the plots

    for (i in IndividualPlots){
      InputVolcano <- subset(InputCompare, individual == paste(i))

      if(nrow(InputVolcano)>=1){
        #Prepare the colour scheme:
        if("color" %in% names(metadata_info)==TRUE){
          color_select <- safe_colorblind_palette[1:length(unique(InputVolcano$color))]

          keyvals <- c()
          for(row in 1:nrow(InputVolcano)){
            col <- color_select[unique(InputVolcano$color) %in% InputVolcano[row, "color"]]
            names(col) <- InputVolcano$color[row]
            keyvals <- c(keyvals, col)
          }
        } else{#here we will use the conditions if no other color is provided!
          color_select <- safe_colorblind_palette[1:length(unique(InputVolcano$comparison))]

          keyvals <- c()
          for(row in 1:nrow(InputVolcano)){
            col <- color_select[unique(InputVolcano$comparison) %in% InputVolcano[row, "comparison"]]
            names(col) <- InputVolcano$comparison[row]
            keyvals <- c(keyvals, col)
          }
        }
        #Prepare the shape scheme:
        if("shape" %in% names(metadata_info)==TRUE & "color" %in% names(metadata_info)==FALSE){
          shape_select <- safe_shape_palette[1:length(unique(InputVolcano$shape))]

          keyvalsshape <- c()
          for(row in 1:nrow(InputVolcano)){
            sha <- shape_select[unique(InputVolcano$shape) %in% InputVolcano[row, "shape"]]
            names(sha) <- InputVolcano$shape[row]
            keyvalsshape <- c(keyvalsshape, sha)
          }
        } else if("shape" %in% names(metadata_info)==TRUE & "color" %in% names(metadata_info)==TRUE){
          #Here we have already used color from metadata_info and we need to use shape for the conditions
          message("For Plot_setting= `Consitions`we can only use colour or shape from metadata_feature. We ignore shape and use it to label the Comparison_name.")
          shape_select <- safe_shape_palette[1:length(unique(InputVolcano$comparison))]

          keyvalsshape <- c()
          for(row in 1:nrow(InputVolcano)){
            sha <- shape_select[unique(InputVolcano$comparison) %in% InputVolcano[row, "comparison"]]
            names(sha) <- InputVolcano$comparison[row]
            keyvalsshape <- c(keyvalsshape, sha)
          }
        } else if("shape" %in% names(metadata_info)==FALSE & "color" %in% names(metadata_info)==FALSE | "shape" %in% names(metadata_info)==FALSE & "color" %in% names(metadata_info)==TRUE){
          shape_select <- safe_shape_palette[1:length(unique(InputVolcano$comparison))]

          keyvalsshape <- c()
          for(row in 1:nrow(InputVolcano)){
            sha <- shape_select[unique(InputVolcano$comparison) %in% InputVolcano[row, "comparison"]]
            names(sha) <- InputVolcano$comparison[row]
            keyvalsshape <- c(keyvalsshape, sha)
          }
        }
        #Prepare the Plot:
        Plot<- EnhancedVolcano::EnhancedVolcano(InputVolcano,
                                                lab = InputVolcano$FeatureNames,#Metabolite name
                                                select_label = select_label,
                                                x = paste(x),
                                                y = paste(y),
                                                xlab  =xlab,
                                                ylab =ylab,
                                                cutoff_stat = cutoff_y,
                                                FCcutoff = cutoff_x,#Cut off Log2FC, automatically 2
                                                pointSize = 3,
                                                labSize = 3,
                                                axisLabSize = 10,
                                                titleLabSize = 12,
                                                subtitleLabSize = 11,
                                                captionLabSize = 10,
                                                col=safe_colorblind_palette,
                                                colCustom = keyvals,
                                                shapeCustom = keyvalsshape,
                                                colAlpha = 1,
                                                title= paste(PlotName, ": ", i, sep=""),
                                                subtitle = subtitle,
                                                caption = paste0("total = ", (nrow(InputVolcano)/2), " ", feature),
                                                xlim =  c(min(InputVolcano[[x]][is.finite(InputVolcano[[x]] )])-0.2, max(InputVolcano[[x]][is.finite(InputVolcano[[x]])])+1.2),
                                                ylim = c(0,(ceiling(-log10(Reduce(min,InputVolcano[[y]]))))),
                                                cutoffLineType = "dashed",
                                                cutoffLineCol = "black",
                                                cutoffLineWidth = 0.5,
                                                legendLabels=c(paste(x," < |", cutoff_x, "|"), paste(x," > |", cutoff_x, "|"), paste(y, ' < ', cutoff_y) , paste(y, ' < ', cutoff_y,' & ',x," < |", cutoff_x, "|")),
                                                legendPosition = 'right',
                                                legendLabSize = 7,
                                                legendIconSize =4,
                                                gridlines.major = FALSE,
                                                gridlines.minor = FALSE,
                                                drawConnectors = connectors)
        #Add the theme
        if(is.null(theme)==FALSE){
          Plot <- Plot+theme
        }

        ## Store the plot in the 'plots' list
        PlotList[[i]] <- Plot

        #Set the total heights and widths
        PlotTitle <- paste(PlotName, ": ", i, sep="")
        Plot_Sized <-  plot_grob_volcano(input_plot=Plot, metadata_info=metadata_info, plot_name = PlotTitle, subtitle = subtitle)
        plot_height <- grid::convertUnit(Plot_Sized$height, 'cm', valueOnly = TRUE)
        plot_width <- grid::convertUnit(Plot_Sized$width, 'cm', valueOnly = TRUE)
        Plot_Sized %<>%
          {ggplot2::ggplot() + annotation_custom(.)} %>%
          add(theme(panel.background = ggplot2::element_rect(fill = "transparent")))

        cleaned_i <- gsub("[[:space:],/\\\\]", "-", i)#removes empty spaces and replaces /,\ with -
        PlotList_adaptedGrid[[cleaned_i]] <- Plot_Sized

        SaveList <- list()
        SaveList[[cleaned_i]] <- Plot_Sized

        #----- Save
        suppressMessages(suppressWarnings(
        save_res(inputlist_df=NULL,
                           inputlist_plot= SaveList,
                           save_table=NULL,
                           save_plot=save_plot,
                           path= folder,
                           file_name= paste("Volcano_",PlotName, sep=""),
                           core=FALSE,
                           print_plot=print_plot,
                           plot_height=plot_height,
                           plot_width=plot_width,
                           plot_unit="cm")))
      }
    }


  } else if("individual" %in% names(metadata_info)==FALSE){
    PlotList <- list()#Empty list to store all the plots
    PlotList_adaptedGrid <- list()#Empty list to store all the plots

    if(nrow(InputCompare)>=1){
      InputVolcano <- InputCompare
      #Prepare the colour scheme:
      if("color" %in% names(metadata_info)==TRUE){
        color_select <- safe_colorblind_palette[1:length(unique(InputVolcano$color))]

        keyvals <- c()
        for(row in 1:nrow(InputVolcano)){
          col <- color_select[unique(InputVolcano$color) %in% InputVolcano[row, "color"]]
          names(col) <- InputVolcano$color[row]
          keyvals <- c(keyvals, col)
        }
      } else{#here we will use the conditions if no other color is provided!
        color_select <- safe_colorblind_palette[1:length(unique(InputVolcano$comparison))]

        keyvals <- c()
        for(row in 1:nrow(InputVolcano)){
          col <- color_select[unique(InputVolcano$comparison) %in% InputVolcano[row, "comparison"]]
          names(col) <- InputVolcano$comparison[row]
          keyvals <- c(keyvals, col)
        }
      }
      #Prepare the shape scheme:
      if("shape" %in% names(metadata_info)==TRUE & "color" %in% names(metadata_info)==FALSE){
        shape_select <- safe_shape_palette[1:length(unique(InputVolcano$shape))]

        keyvalsshape <- c()
        for(row in 1:nrow(InputVolcano)){
          sha <- shape_select[unique(InputVolcano$shape) %in% InputVolcano[row, "shape"]]
          names(sha) <- InputVolcano$shape[row]
          keyvalsshape <- c(keyvalsshape, sha)
        }
      } else if("shape" %in% names(metadata_info)==TRUE & "color" %in% names(metadata_info)==TRUE){
        #Here we have already used color from metadata_info and we need to use shape for the conditions
        message("For plot_types Comparison we can only use colour or shape from metadata_feature. Hence, we ignore shape and use it to label the name_comparison.")
        shape_select <- safe_shape_palette[1:length(unique(InputVolcano$comparison))]

        keyvalsshape <- c()
        for(row in 1:nrow(InputVolcano)){
          sha <- shape_select[unique(InputVolcano$comparison) %in% InputVolcano[row, "comparison"]]
          names(sha) <- InputVolcano$comparison[row]
          keyvalsshape <- c(keyvalsshape, sha)
        }
      } else if("shape" %in% names(metadata_info)==FALSE & "color" %in% names(metadata_info)==FALSE | "shape" %in% names(metadata_info)==FALSE & "color" %in% names(metadata_info)==TRUE){
        shape_select <- safe_shape_palette[1:length(unique(InputVolcano$comparison))]

        keyvalsshape <- c()
        for(row in 1:nrow(InputVolcano)){
          sha <- shape_select[unique(InputVolcano$comparison) %in% InputVolcano[row, "comparison"]]
          names(sha) <- InputVolcano$comparison[row]
          keyvalsshape <- c(keyvalsshape, sha)
        }
      }
      #Prepare the Plot:
      Plot<- EnhancedVolcano::EnhancedVolcano(InputVolcano,
                                              lab = InputVolcano$FeatureNames,#Metabolite name
                                              select_label = select_label,
                                              x = paste(x),
                                              y = paste(y),
                                              xlab  =xlab,
                                              ylab =ylab,
                                              cutoff_stat = cutoff_y,
                                              FCcutoff = cutoff_x,#Cut off Log2FC, automatically 2
                                              pointSize = 3,
                                              labSize = 3,
                                              axisLabSize = 10,
                                              titleLabSize = 12,
                                              subtitleLabSize = 11,
                                              captionLabSize = 10,
                                              col=safe_colorblind_palette,
                                              colCustom = keyvals,
                                              shapeCustom = keyvalsshape,
                                              colAlpha = 1,
                                              title= paste(PlotName),
                                              subtitle = subtitle,
                                              caption = paste0("total = ", (nrow(InputVolcano)/2)," ", feature),
                                              xlim =  c(min(InputVolcano[[x]][is.finite(InputVolcano[[x]] )])-0.2, max(InputVolcano[[x]][is.finite(InputVolcano[[x]])])+1.2),
                                              ylim = c(0,(ceiling(-log10(Reduce(min,InputVolcano[[y]]))))),
                                              cutoffLineType = "dashed",
                                              cutoffLineCol = "black",
                                              cutoffLineWidth = 0.5,
                                              legendLabels=c(paste(x," < |", cutoff_x, "|"), paste(x," > |", cutoff_x, "|"), paste(y, ' < ', cutoff_y) , paste(y, ' < ', cutoff_y,' & ',x," < |", cutoff_x, "|")),
                                              legendPosition = 'right',
                                              legendLabSize = 7,
                                              legendIconSize =4,
                                              gridlines.major = FALSE,
                                              gridlines.minor = FALSE,
                                              drawConnectors = connectors)
      #Add the theme
      if(is.null(theme)==FALSE){
        Plot <- Plot+theme
      }

      ## Store the plot in the 'plots' list
      PlotList[["Plot"]] <- Plot

      Plot_Sized <-  plot_grob_volcano(input_plot=Plot, metadata_info=metadata_info, plot_name =plot_name, subtitle = subtitle)
      plot_height <- grid::convertUnit(Plot_Sized$height, 'cm', valueOnly = TRUE)
      plot_width <- grid::convertUnit(Plot_Sized$width, 'cm', valueOnly = TRUE)
      Plot_Sized %<>%
        {ggplot2::ggplot() + annotation_custom(.)} %>%
        add(theme(panel.background = ggplot2::element_rect(fill = "transparent")))

      PlotList_adaptedGrid[["Plot_Sized"]] <- Plot_Sized

       #----- Save
      suppressMessages(suppressWarnings(
      save_res(inputlist_df=NULL,
                             inputlist_plot= list("Plot_Sized"= PlotList_adaptedGrid[["Plot_Sized"]]),
                             save_table=NULL,
                             save_plot=save_plot,
                             path= folder,
                             file_name= paste("Volcano_",plot_name, sep=""),
                             core=FALSE,
                             print_plot=print_plot,
                             plot_height=plot_height,
                             plot_width=plot_width,
                             plot_unit="cm")))

    }
  }
  return(invisible(list("Plot"=PlotList,"Plot_Sized" = PlotList_adaptedGrid)))
}



################################################################################################
### ### ### viz_volcano helper function: Internal Function for plot_types PEA ### ### ###
################################################################################################

#' Check input parameters
#'
#' @param data Passed to main function viz_volcano()
#' @param data2 Passed to main function viz_volcano()
#' @param metadata_feature Passed to main function viz_volcano()
#' @param metadata_info Passed to main function viz_volcano()
#' @param y \emph{Optional: } Passed to main function viz_volcano() \strong{Default = "p.adj"}
#' @param x \emph{Optional: } Passed to main function viz_volcano() \strong{Default = "Log2FC"}
#' @param plot_name \emph{Optional: } Passed to main function viz_volcano() \strong{Default = ""}
#' @param xlab \emph{Optional: } Passed to main function viz_volcano()  \strong{Default = NULL}
#' @param ylab \emph{Optional: } Passed to main function viz_volcano() \strong{Default = NULL}
#' @param cutoff_x \emph{Optional: } Passed to main function viz_volcano() \strong{Default = 0.5}
#' @param cutoff_y \emph{Optional: } Passed to main function viz_volcano() \strong{Default = 0.05}
#' @param select_label \emph{Optional: } Passed to main function viz_volcano() \strong{Default = ""}
#' @param connectors \emph{Optional: } Passed to main function viz_volcano() \strong{Default =  FALSE}
#' @param subtitle \emph{Optional: } Passed to main function viz_volcano() \strong{Default = ""}
#' @param color_palette Created in viz_volcano() based on color_palette passed to main function viz_volcano()
#' @param shape_palette Created in viz_volcano() based on shape_palette passed to main function viz_volcano()
#' @param theme \emph{Optional: } Selection of theme for plot, e.g. theme_grey(). You can check for complete themes here: https://ggplot2.tidyverse.org/reference/ggtheme.html. \strong{Default = NULL}
#' @param feature \emph{Optional: } Name of the feature that are plotted, e.g. "Metabolites", "RNA", "Proteins", "Genes", etc. \strong{Default = "Metabolites"}
#' @param save_plot Passed to main function viz_volcano()
#' @param print_plot Passed to main function viz_volcano()
#' @param folder Created in viz_volcano(). Path to the folder where files are saved.
#'
#' @return List with two elements: Plot and Plot_Sized
#'
#' @keywords Volcano plots of pathway enrichment results
#'
#' @importFrom ggplot2 ggplot theme
#' @importFrom dplyr rename filter mutate
#' @importFrom magrittr %>% %<>%
#' @importFrom tibble rownames_to_column column_to_rownames remove_rownames
#' @importFrom logger log_trace
#'
#' @noRd
#'
viz_volcano_pea <- function(data,
                           data2,
                           metadata_feature,
                           metadata_info,
                           y= "p.adj",
                           x= "Log2FC",
                           xlab= NULL,#"~Log[2]~FC"
                           ylab= NULL,#"~-Log[10]~p.adj"
                           cutoff_x= 0.5,
                           cutoff_y= 0.05,
                           connectors=  FALSE,
                           select_label= "",
                          plot_name= "",
                           subtitle= "",
                           color_palette,
                           shape_palette,
                           theme= NULL,
                           feature="Metabolites",
                           save_plot,
                           print_plot,
                           folder){
  #####################
  ##--- Check PEA settings
  if(is.vector(metadata_info)==FALSE){
    message <- paste0(
      "You have chosen Settings=`PEA` that requires you to ",
      "provide a vector for metadata_info."
    )
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }
  if(is.null(metadata_feature)==TRUE){
    message <- paste0(
      "You have chosen Settings=`PEA` that requires you to provide a DF ",
      "metadata_feature including the pathways used for the enrichment analysis."
    )
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }
  if(is.null(metadata_feature)==FALSE & is.null(metadata_feature)==FALSE){
    if("PEA_Feature" %in% names(metadata_info)==FALSE | "PEA_score" %in% names(metadata_info)==FALSE | "PEA_stat" %in% names(metadata_info)==FALSE | "PEA_Pathway" %in% names(metadata_info)==FALSE){
      message <- paste0("You have chosen Settings=`PEA` that requires you to provide a vector for metadata_info including `PEA_Feature`, `PEA_Pathway`, `PEA_stat` and `PEA_score`.")
      logger::log_trace(paste("Error ", message, sep=""))
      stop(message)
    }
  }

  #Pass colours/shapes
  safe_colorblind_palette <- color_palette
  safe_shape_palette <- shape_palette

  #Prepare data:
  data <- data%>%
    dplyr::rename("PEA_Feature"= !!metadata_info[["PEA_Feature"]])


  data2 <- data2%>%
    dplyr::rename("PEA_score"= !!metadata_info[["PEA_score"]],
                  "PEA_stat"= !!metadata_info[["PEA_stat"]],
                  "PEA_Pathway"= !!metadata_info[["PEA_Pathway"]])

  metadata_feature <- metadata_feature%>%
    dplyr::rename("PEA_Pathway"= !!metadata_info[["PEA_Pathway"]],
                  "PEA_Feature"= !!metadata_info[["PEA_Feature"]])

  #################
  ##--- Plot
  # Create the list of individual plots that should be made:
  IndividualPlots <- unique(data2$PEA_Pathway)

  PlotList <- list()#Empty list to store all the plots
  PlotList_adaptedGrid <- list()#Empty list to store all the plots

  for (i in IndividualPlots){
    data2_Select<- data2%>%
      filter(PEA_Pathway == paste(i)) #Select pathway we plot and use the score and stats

    metadata_feature_Select <- metadata_feature%>%
      filter(PEA_Pathway == paste(i))

    InputVolcano <-merge(metadata_feature_Select, data, by="PEA_Feature", all.x=TRUE)%>%
      distinct(PEA_Feature, .keep_all = TRUE) %>%
      filter(!is.na(!!sym(y)) & !is.na(!!sym(x)))

    if(nrow(InputVolcano)>=1){
      #Prepare the colour scheme:
      if("color" %in% names(metadata_info)==TRUE){
        color_select <- safe_colorblind_palette[1:length(unique(InputVolcano$color))]

        keyvals <- c()
        for(row in 1:nrow(InputVolcano)){
          col <- color_select[unique(InputVolcano$color) %in% InputVolcano[row, "color"]]
          names(col) <- InputVolcano$color[row]
          keyvals <- c(keyvals, col)
        }

        LegendPos<- "right"
      } else{
        keyvals <-NULL
      }
      #Prepare the shape scheme:
      if("shape" %in% names(metadata_info)==TRUE){
        shape_select <- safe_shape_palette[1:length(unique(InputVolcano$shape))]

        keyvalsshape <- c()
        for(row in 1:nrow(InputVolcano)){
          sha <- shape_select[unique(InputVolcano$shape) %in% InputVolcano[row, "shape"]]
          names(sha) <- InputVolcano$shape[row]
          keyvalsshape <- c(keyvalsshape, sha)
        }

        LegendPos<- "right"
      } else{
        keyvalsshape <-NULL
      }

      if("color" %in% names(metadata_info)==FALSE & "shape" %in% names(metadata_info)==FALSE){
        LegendPos<- "none"
      }

      #Prepare the Plot:
      Plot<- EnhancedVolcano::EnhancedVolcano(InputVolcano,
                                              lab = InputVolcano$PEA_Feature,#Metabolite name
                                              select_label = select_label,
                                              x = paste(x),
                                              y = paste(y),
                                              xlab  =xlab,
                                              ylab =ylab,
                                              cutoff_stat = cutoff_y,
                                              FCcutoff = cutoff_x,#Cut off Log2FC, automatically 2
                                              pointSize = 3,
                                              labSize = 3,
                                              axisLabSize = 10,
                                              titleLabSize = 12,
                                              subtitleLabSize = 11,
                                              captionLabSize = 10,
                                              col=safe_colorblind_palette,
                                              colCustom = keyvals,
                                              shapeCustom = keyvalsshape,
                                              colAlpha = 1,
                                              title= paste(PlotName, ": ", i, sep=""),
                                              subtitle = paste(metadata_info[["PEA_score"]],"= ", data2_Select$PEA_score, ", ",metadata_info[["PEA_stat"]] , "= ", data2_Select$PEA_stat, sep=""),
                                              caption = paste0("total = ", nrow(InputVolcano), " of ", nrow(metadata_feature_Select), " ", feature, " in pathway"),
                                              xlim =  c(min(InputVolcano[[x]][is.finite(InputVolcano[[x]] )])-0.2, max(InputVolcano[[x]][is.finite(InputVolcano[[x]])])+1.2),
                                              ylim = c(0,(ceiling(-log10(Reduce(min,InputVolcano[[y]]))))),
                                              cutoffLineType = "dashed",
                                              cutoffLineCol = "black",
                                              cutoffLineWidth = 0.5,
                                              legendLabels=c(paste(x," < |", cutoff_x, "|"), paste(x," > |", cutoff_x, "|"), paste(y, ' < ', cutoff_y) , paste(y, ' < ', cutoff_y,' & ',x," < |", cutoff_x, "|")),
                                              legendPosition = LegendPos,
                                              legendLabSize = 7,
                                              legendIconSize =4,
                                              gridlines.major = FALSE,
                                              gridlines.minor = FALSE,
                                              drawConnectors = connectors)
      #Add the theme
      if(is.null(theme)==FALSE){
        Plot <- Plot+theme
      }

      ## Store the plot in the 'plots' list
      PlotList[[i]] <- Plot

      #Set the total heights and widths
      PlotTitle <- paste(PlotName, ": ", i, sep="")
      Plot_Sized <-  plot_grob_volcano(input_plot=Plot, metadata_info=metadata_info, plot_name = PlotTitle, subtitle = subtitle)
      plot_height <- grid::convertUnit(Plot_Sized$height, 'cm', valueOnly = TRUE)
      plot_width <- grid::convertUnit(Plot_Sized$width, 'cm', valueOnly = TRUE)

      Plot_Sized %<>%
        {ggplot2::ggplot() + annotation_custom(.)} %>%
        add(theme(panel.background = ggplot2::element_rect(fill = "transparent")))

      cleaned_i <- gsub("[[:space:],/\\\\]", "-", i)#removes empty spaces and replaces /,\ with -
      PlotList_adaptedGrid[[cleaned_i]] <- Plot_Sized

      SaveList <- list()
      SaveList[[cleaned_i]] <- Plot_Sized

      #----- Save
      suppressMessages(suppressWarnings(
        save_res(inputlist_df=NULL,
                             inputlist_plot= SaveList,
                             save_table=NULL,
                             save_plot=save_plot,
                             path= folder,
                             file_name= paste("Volcano_",PlotName, sep=""),
                             core=FALSE,
                             print_plot=print_plot,
                             plot_height=plot_height,
                             plot_width=plot_width,
                             plot_unit="cm")))

    }
  }

 return(invisible(list("Plot"=PlotList,"Plot_Sized" = PlotList_adaptedGrid)))
}
