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

#####################################
### ### ### Volcano Plots ### ### ###
#####################################
#' @param Settings \emph{Optional: } Choose between "Standard" (InputData), "Compare" (plot two comparisons together InputData and InputData2) or "PEA" (Pathway Enrichment Analysis) \strong{Default = "Standard"}
#' @param SettingsInfo \emph{Optional: } NULL or Named vector including at least one of those three information for Settings="Standard" or "Compare": c(color ="ColumnName_SettingsFile_Metab", shape = "ColumnName_SettingsFile_Metab", individual="ColumnName_SettingsFile_Metab"). For Settings="PEA" a named vector with: PEA_Pathway="ColumnName_InputData2"=each pathway will be plotted, PEA_score="ColumnName_InputData2", PEA_stat= "ColumnName_InputData2"= usually p.adj column, "PEA_Feature="ColumnName_InputData2"= usually Metabolites), optionally you can additionally include c(color_Metab="ColumnName_SettingsFile_Metab", shape= "ColumnName_SettingsFile_Metab").\strong{Default = NULL}
#' @param SettingsFile_Metab \emph{Optional: } DF with column including the Metabolite names (needs to match Metabolite names and Metabolite column name of InputData) and other columns with required PlotSettingInfo. \strong{Default = NULL}
#' @param Input_data DF with metabolites as row names and columns including Log2FC and stat (p-value, p.adjusted) value columns.
#' @param InputData2 \emph{Optional: } DF to compare to main Input_data with the same column names x and y (Settings="Compare") and metabolites as row names or Pathway enrichment analysis results (Settings="PEA"). \strong{Default = NULL}
#' @param y \emph{Optional: } Column name including the values that should be used for y-axis. Usually this would include the p.adjusted value. \strong{Default = "p.adj"}
#' @param x \emph{Optional: } Column name including the values that should be used for x-axis. Usually this would include the Log2FC value. \strong{Default = "Log2FC"}
#' @param FeatureID {Optional: } Column name including the feature names, e.g. metabolite names. \strong{Default = "Metabolite"}
#' @param PlotName \emph{Optional: } String which is added to the output files of the plot. \strong{Default = ""}
#' @param ComparisonName \emph{Optional: } Named vector including those information about the two datasets that are compared on the plots when choosing Settings= "Compare". \strong{Default = c(InputData="Cond1", InputData2= "Cond2")}
#' @param xlab \emph{Optional: } String to replace x-axis label in plot. \strong{Default = NULL}
#' @param ylab \emph{Optional: } String to replace y-axis label in plot. \strong{Default = NULL}
#' @param xCutoff \emph{Optional: } Number of the desired log fold change cutoff for assessing significance. \strong{Default = 0.5}
#' @param yCutoff \emph{Optional: } Number of the desired p value cutoff for assessing significance. \strong{Default = 0.05}
#' @param ColorPalette \emph{Optional: } Provide customiced color-palette in vector format. \strong{Default = NULL}
#' @param ShapePalette \emph{Optional: } Provide customiced shape-palette in vector format. \strong{Default = NULL}
#' @param SelectLab \emph{Optional: } If set to NULL, feature labels will be plotted randomly. If vector is provided, e.g. c("MetaboliteName1", "MetaboliteName2"), selected names will be plotted. If set to default "", no feature names will be plotted. \strong{Default = ""}
#' @param Connectors \emph{Optional: } TRUE or FALSE for whether Connectors from names to points are to be added to the plot. \strong{Default =  FALSE}
#' @param Subtitle \emph{Optional: } \strong{Default = ""}
#' @param Theme \emph{Optional: } Selection of theme for plot, e.g. theme_grey(). You can check for complete themes here: https://ggplot2.tidyverse.org/reference/ggtheme.html. \strong{Default = NULL}
#' @param FolderPath {Optional:} Path to the folder the results should be saved at. \strong{default: NULL}
#' @param Features \emph{Optional: } Name of the features that are plotted, e.g. "Metabolites", "RNA", "Proteins", "Genes", etc. \strong{Default = "metabolites"}
#' @param SaveAs_Plot \emph{Optional: } Select the file type of output plots. Options are svg, pdf, png or NULL. \strong{Default = "svg"}
#'
#' @keywords Volcano plot, pathways
#' @export

# Helper function needed for adding column to pathway file defining if this metabolite is unique/multiple pathways
VizVolcano <- function(PlotSettings="Standard",
                       InputData,
                       SettingsInfo= NULL,
                       SettingsFile_Metab=NULL,
                       InputData2= NULL,
                       y= "p.adj",
                       x= "Log2FC",
                       xlab= NULL,#"~Log[2]~FC"
                       ylab= NULL,#"~-Log[10]~p.adj"
                       xCutoff= 0.5,
                       yCutoff= 0.05,
                       Connectors=  FALSE,
                       SelectLab= "",
                       PlotName= "",
                       Subtitle= "",
                       ComparisonName= c(InputData="Cond1", InputData2= "Cond2"),
                       ColorPalette= NULL,
                       ShapePalette=NULL,
                       Theme= NULL,
                       SaveAs_Plot= "svg",
                       FolderPath = NULL,
                       Features="Metabolites",
                       PrintPlot=TRUE){

  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse", "EnhancedVolcano")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0){
    install.packages(new.packages)
  }
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0){
    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")

    BiocManager::install(new.packages)
  }
  suppressMessages(library(tidyverse))

  ################################################################################################################################################################################################
  ## ------------ Check Input files ----------- ##
  # HelperFunction `CheckInput`
  if(PlotSettings=="PEA"){
    #Those relationships are checked in the VizVolcano_PEA() function!
    SettingsFile <- NULL # For PEA the SettingsFile_Metab is the prior knowledge file, and hence this will not have features as row names.
    Info <- NULL # If SettingsFileMetab=NULL, SetingsInfo has to be NULL to, otherwise we will get an error.
  }else{
    SettingsFile <-SettingsFile_Metab
    Info <- SettingsInfo
  }

  MetaProViz:::CheckInput(InputData=as.data.frame(t(InputData)),
                          InputData_Num=FALSE,
                          SettingsFile_Sample=NULL,
                          SettingsFile_Metab=SettingsFile,#Set above
                          SettingsInfo=Info,#Set above
                          SaveAs_Plot=SaveAs_Plot,
                          SaveAs_Table=NULL,
                          CoRe=FALSE,
                          PrintPlot= PrintPlot,
                          PlotSettings="Feature")

  # CheckInput` Specific:
  if(is.numeric(yCutoff)== FALSE |yCutoff > 1 | yCutoff < 0){
      stop("Check input. The selected yCutoff value should be numeric and between 0 and 1.")
    }
  if(is.numeric(xCutoff)== FALSE  | xCutoff < 0){
      stop("Check input. The selected xCutoff value should be numeric and between 0 and +oo.")
    }
  if(paste(x) %in% colnames(InputData)==FALSE | paste(y) %in% colnames(InputData)==FALSE){
    stop("Check your input. The column name of x and/ore y does not exist in Input_data.")
  }

  if(is.null(SelectLab)==FALSE & is.vector(SelectLab)==FALSE){
    stop("Check input. SelectedLab must be either NULL or a vector.")
  }

  if(is.logical(Connectors) == FALSE){
    stop("Check input. The Connectors value should be either = TRUE if connectors from names to points are to be added to the plot or =FALSE if not.")
  }

  if(is.null(PlotName)==FALSE & is.vector(PlotName)==FALSE){
    stop("Check input. PlotName must be either NULL or a vector.")
  }

  Plot_options <- c("Standard", "Compare", "PEA")
  if (PlotSettings %in% Plot_options == FALSE){
    stop("PlotSettings option is incorrect. The allowed options are the following: ",paste(Plot_options, collapse = ", "),"." )
  }


  ## ------------ Create Results output folder ----------- ##
  if(is.null(SaveAs_Plot)==FALSE){
    Folder <- MetaProViz:::SavePath(FolderName= "VolcanoPlots",
                                    FolderPath=FolderPath)
  }

  ############################################################################################################
  ## ----------- Prepare InputData ------------ ##
  #Extract required columns and merge with SettingsFile
   if(is.null(SettingsFile_Metab)==FALSE){
     ##--- Prepare the color scheme:
     if("color" %in% names(SettingsInfo)==TRUE & "shape" %in% names(SettingsInfo)==TRUE){
       if((SettingsInfo[["shape"]] == SettingsInfo[["color"]])==TRUE){
         SettingsFile_Metab$shape <- SettingsFile_Metab[,paste(SettingsInfo[["color"]])]
         SettingsFile_Metab<- SettingsFile_Metab%>%
           dplyr::rename("color"=paste(SettingsInfo[["color"]]))
       }else{
         SettingsFile_Metab <- SettingsFile_Metab%>%
           dplyr::rename("color"=paste(SettingsInfo[["color"]]),
                         "shape"=paste(SettingsInfo[["shape"]]))
       }
     }else if("color" %in% names(SettingsInfo)==TRUE & "shape" %in% names(SettingsInfo)==FALSE){
       SettingsFile_Metab <- SettingsFile_Metab%>%
         dplyr::rename("color"=paste(SettingsInfo[["color"]]))
     }else if("color" %in% names(SettingsInfo)==FALSE & "shape" %in% names(SettingsInfo)==TRUE){
       SettingsFile_Metab <- SettingsFile_Metab%>%
         dplyr::rename("shape"=paste(SettingsInfo[["shape"]]))
     }
     if("individual" %in% names(SettingsInfo)==TRUE){
       SettingsFile_Metab <- SettingsFile_Metab%>%
         dplyr::rename("individual"=paste(SettingsInfo[["individual"]]))
     }


    ##--- Merge InputData with SettingsFile:
     common_columns <- character(0)  # Initialize an empty character vector
     for(col_name in colnames(InputData[, c(x, y)])) {
       if(col_name %in% colnames(SettingsFile_Metab)) {
         common_columns <- c(common_columns, col_name)  # Add the common column name to the vector
       }
     }
     SettingsFile_Metab <- SettingsFile_Metab%>%#rename those column since they otherwise will cause issues when we merge the DFs later
      dplyr::rename_at(vars(common_columns), ~ paste0(., "_SettingsFile_Metab"))

    if(PlotSettings=="PEA"){
      VolcanoData <- merge(x=SettingsFile_Metab ,y=InputData[, c(x, y)], by.x=SettingsInfo[["PEA_Feature"]] , by.y=0, all.y=TRUE)%>%
        remove_rownames()%>%
        mutate(FeatureNames = SettingsInfo[["PEA_Feature"]])%>%
        filter(!is.na(x) | !is.na(x))
    }else{
     VolcanoData <- merge(x=SettingsFile_Metab ,y=InputData[, c(x, y)], by=0, all.y=TRUE)%>%
      remove_rownames()%>%
      column_to_rownames("Row.names")%>%
      mutate(FeatureNames = rownames(InputData))%>%
      filter(!is.na(x) | !is.na(x))
    }

   }else{
     VolcanoData <- InputData[, c(x, y)]%>%
       mutate(FeatureNames = rownames(InputData))%>%
       filter(!is.na(x) | !is.na(x))
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
  if(is.null(ColorPalette)){
    if("color" %in% names(SettingsInfo)==TRUE){
      safe_colorblind_palette <- c("#88CCEE",  "#DDCC77","#661100",  "#332288", "#AA4499","#999933",  "#44AA99", "#882215",  "#6699CC", "#117733", "#888888","#CC6677", "black","gold1","darkorchid4","red","orange")
    }else{
      safe_colorblind_palette <- c("#888888", "#44AA99", "#44AA99","#CC6677")
    }

    #check that length is enough for what the user wants to colour
    #stop(" The maximum number of pathways in the Input_pathways must be less than ",length(safe_colorblind_palette),". Please summarize sub-pathways together where possible and repeat.")
  } else{
    safe_colorblind_palette <-ColorPalette
    #check that length is enough for what the user wants to colour
  }
  if(is.null(ShapePalette)){
    safe_shape_palette <- c(15,17,16,18,25,7,8,11,12)
    #check that length is enough for what the user wants to shape
  } else{
    safe_shape_palette <-shape_palette
    #check that length is enough for what the user wants to shape
  }

  ############################################################################################################
  ## ----------- Make the  plot based on the chosen parameters ------------ ##

  if(PlotSettings=="Standard"){#####--- 1. Standard
    VolcanoRes <- MetaProViz:::VizVolcano_Standard(InputData= VolcanoData,
                                                   SettingsFile_Metab=SettingsFile_Metab,
                                                   SettingsInfo=SettingsInfo,
                                                   y= y,
                                                   x= x,
                                                   xlab= xlab,
                                                   ylab= ylab,
                                                   xCutoff= xCutoff,
                                                   yCutoff= yCutoff,
                                                   Connectors= Connectors,
                                                   SelectLab=SelectLab,
                                                   PlotName= PlotName,
                                                   Subtitle= Subtitle,
                                                   ColorPalette=safe_colorblind_palette,
                                                   ShapePalette=safe_shape_palette,
                                                   Theme= Theme,
                                                   Features=Features,
                                                   SaveAs_Plot=SaveAs_Plot,
                                                   PrintPlot=PrintPlot,
                                                   Folder=Folder)

  }else if(PlotSettings=="Compare"){#####--- 2. Compare
    VolcanoRes <- MetaProViz:::VizVolcano_Compare(InputData= VolcanoData,
                                                  InputData2=InputData2,
                                                  SettingsFile_Metab=SettingsFile_Metab,
                                                  SettingsInfo=SettingsInfo,
                                                  y= y,
                                                  x= x,
                                                  xlab= xlab,
                                                  ylab= ylab,
                                                  xCutoff= xCutoff,
                                                  yCutoff= yCutoff,
                                                  Connectors= Connectors,
                                                  SelectLab=SelectLab,
                                                  PlotName= PlotName,
                                                  Subtitle= Subtitle,
                                                  ColorPalette=safe_colorblind_palette,
                                                  ShapePalette=safe_shape_palette,
                                                  Theme= Theme,
                                                  Features=Features,
                                                  ComparisonName=ComparisonName,
                                                  SaveAs_Plot=SaveAs_Plot,
                                                  PrintPlot=PrintPlot,
                                                  Folder=Folder)

  } else if(PlotSettings=="PEA"){#####--- 3. PEA
    VolcanoRes <- MetaProViz:::VizVolcano_PEA(InputData= VolcanoData,
                                              InputData2=InputData2,
                                              SettingsFile_Metab=SettingsFile_Metab,#Problem: we need to know the column name of the features!
                                              SettingsInfo=SettingsInfo,
                                              y= y,
                                              x= x,
                                              xlab= xlab,
                                              ylab= ylab,
                                              xCutoff= xCutoff,
                                              yCutoff= yCutoff,
                                              Connectors= Connectors,
                                              SelectLab=SelectLab,
                                              PlotName= PlotName,
                                              Subtitle= Subtitle,
                                              ColorPalette=safe_colorblind_palette,
                                              ShapePalette=safe_shape_palette,
                                              Theme= Theme,
                                              Features=Features,
                                              SaveAs_Plot=SaveAs_Plot,
                                              PrintPlot=PrintPlot,
                                              Folder=Folder)
  }
  return(invisible(VolcanoRes))
}

################################################################################################
### ### ### VizVolcano helper function: Internal Function for PlotSettings Standard ### ### ###
################################################################################################

#' Check input parameters
#'
#' @param InputData Passed to main function MetaProViz::VizVolcano()
#' @param SettingsFile_Metab Passed to main function MetaProViz::VizVolcano()
#' @param SettingsInfo Passed to main function MetaProViz::VizVolcano()
#' @param y \emph{Optional: } Passed to main function MetaProViz::VizVolcano() \strong{Default = "p.adj"}
#' @param x \emph{Optional: } Passed to main function MetaProViz::VizVolcano() \strong{Default = "Log2FC"}
#' @param PlotName \emph{Optional: } Passed to main function MetaProViz::VizVolcano() \strong{Default = ""}
#' @param xlab \emph{Optional: } Passed to main function MetaProViz::VizVolcano()  \strong{Default = NULL}
#' @param ylab \emph{Optional: } Passed to main function MetaProViz::VizVolcano() \strong{Default = NULL}
#' @param xCutoff \emph{Optional: } Passed to main function MetaProViz::VizVolcano() \strong{Default = 0.5}
#' @param ycutoff \emph{Optional: } Passed to main function MetaProViz::VizVolcano() \strong{Default = 0.05}
#' @param SelectLab \emph{Optional: } Passed to main function MetaProViz::VizVolcano() \strong{Default = ""}
#' @param Connectors \emph{Optional: } Passed to main function MetaProViz::VizVolcano() \strong{Default =  FALSE}
#' @param Subtitle \emph{Optional: } Passed to main function MetaProViz::VizVolcano() \strong{Default = ""}
#' @param ColorPalette Created in MetaProViz::VizVolcano() based on ColorPalette passed to main function MetaProViz::VizVolcano()
#' @param ShapePalette Created in MetaProViz::VizVolcano() based on ShapePalette passed to main function MetaProViz::VizVolcano()
#' @param Theme \emph{Optional: } Selection of theme for plot, e.g. theme_grey(). You can check for complete themes here: https://ggplot2.tidyverse.org/reference/ggtheme.html. \strong{Default = NULL}
#' @param Features \emph{Optional: } Name of the features that are plotted, e.g. "Metabolites", "RNA", "Proteins", "Genes", etc. \strong{Default = "Metabolites"}
#' @param SaveAs_Plot Passed to main function MetaProViz::VizVolcano()
#' @param PrintPlot Passed to main function MetaProViz::VizVolcano()
#' @param Folder Created in MetaProViz::VizVolcano(). Path to the folder where files are saved.
#'

#' @keywords Standard volcano plots
#' @noRd
#'
#'

VizVolcano_Standard <- function(InputData,
                                SettingsFile_Metab,
                                SettingsInfo,
                                y= "p.adj",
                                x= "Log2FC",
                                xlab= NULL,#"~Log[2]~FC"
                                ylab= NULL,#"~-Log[10]~p.adj"
                                xCutoff= 0.5,
                                yCutoff= 0.05,
                                Connectors=  FALSE,
                                SelectLab= "",
                                PlotName= "",
                                Subtitle= "",
                                ColorPalette,
                                ShapePalette,
                                Theme= NULL,
                                Features="Metabolites",
                                SaveAs_Plot,
                                PrintPlot,
                                Folder){

  #Pass colours/shapes
  safe_colorblind_palette <- ColorPalette
  safe_shape_palette <- ShapePalette

  #Plots
  if("individual" %in% names(SettingsInfo)==TRUE){
    # Create the list of individual plots that should be made:
    IndividualPlots <- unique(InputData$individual)

    PlotList <- list()#Empty list to store all the plots
    PlotList_adaptedGrid <- list()#Empty list to store all the plots

    for (i in IndividualPlots){
      InputVolcano <- subset(InputData, individual == paste(i))

      if(nrow(InputVolcano)>=1){
        if("color" %in% names(SettingsInfo)==TRUE ){
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
        if("shape" %in% names(SettingsInfo)==TRUE){
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

        if("color" %in% names(SettingsInfo)==FALSE & "shape" %in% names(SettingsInfo)==FALSE){
          LegendPos<- "none"
        }

        #Prepare the Plot:
        Plot<- EnhancedVolcano::EnhancedVolcano(InputVolcano,
                                                lab = InputVolcano$FeatureNames,#Metabolite name
                                                selectLab = SelectLab,
                                                x = paste(x),
                                                y = paste(y),
                                                xlab  =xlab,
                                                ylab =ylab,
                                                pCutoff = yCutoff,
                                                FCcutoff = xCutoff,#Cut off Log2FC, automatically 2
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
                                                subtitle = Subtitle,
                                                caption = paste0("Total = ", nrow(InputVolcano), " ", Features),
                                                xlim =  c(min(InputVolcano[[x]][is.finite(InputVolcano[[x]] )])-0.2, max(InputVolcano[[x]][is.finite(InputVolcano[[x]])])+1.2),
                                                ylim = c(0,(ceiling(-log10(Reduce(min,InputVolcano[[y]]))))),
                                                cutoffLineType = "dashed",
                                                cutoffLineCol = "black",
                                                cutoffLineWidth = 0.5,
                                                legendLabels=c(paste(x," < |", xCutoff, "|"), paste(x," > |", xCutoff, "|"), paste(y, ' < ', yCutoff) , paste(y, ' < ', yCutoff,' & ',x," < |", xCutoff, "|")),
                                                legendPosition = LegendPos,
                                                legendLabSize = 7,
                                                legendIconSize =4,
                                                gridlines.major = FALSE,
                                                gridlines.minor = FALSE,
                                                drawConnectors = Connectors)
        #Add the theme
        if(is.null(Theme)==FALSE){
          Plot <- Plot+Theme
        }

        ## Store the plot in the 'plots' list
        PlotList[[i]] <- Plot

        #Set the total heights and widths
        PlotTitle <- paste(PlotName, ": ", i, sep="")
        Plot_Sized <-  MetaProViz:::plotGrob_Volcano(Input=Plot, keyvals = keyvals, keyvalsshape = keyvalsshape, PlotName = PlotTitle, Subtitle = Subtitle)
        Plot <-Plot_Sized[[3]]
        Plot <- ggplot2::ggplot() +
          annotation_custom(Plot)
        Plot <-Plot + theme(panel.background = element_rect(fill = "transparent"))

        cleaned_i <- gsub("[[:space:],/\\\\]", "-", i)#removes empty spaces and replaces /,\ with -
        PlotList_adaptedGrid[[cleaned_i]] <- Plot

        SaveList <- list()
        SaveList[[cleaned_i]] <- Plot
        #----- Save
        suppressMessages(suppressWarnings(
          MetaProViz:::SaveRes(InputList_DF=NULL,
                               InputList_Plot= SaveList,
                               SaveAs_Table=NULL,
                               SaveAs_Plot=SaveAs_Plot,
                               FolderPath= Folder,
                               FileName= paste("Volcano_",PlotName, sep=""),
                               CoRe=FALSE,
                               PrintPlot=PrintPlot,
                               PlotHeight=Plot_Sized[[1]],
                               PlotWidth=Plot_Sized[[2]],
                               PlotUnit="cm")))
      }
    }
  }else if("individual" %in% names(SettingsInfo)==FALSE){
    PlotList <- list()#Empty list to store all the plots
    PlotList_adaptedGrid <- list()#Empty list to store all the plots

    InputVolcano <- InputData
    if(nrow(InputVolcano)>=1){
      if("color" %in% names(SettingsInfo)==TRUE ){
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
      if("shape" %in% names(SettingsInfo)==TRUE){
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

      if("color" %in% names(SettingsInfo)==FALSE & "shape" %in% names(SettingsInfo)==FALSE){
        LegendPos<- "none"
      }

      #Prepare the Plot:
      Plot<- EnhancedVolcano::EnhancedVolcano(InputVolcano,
                                              lab = InputVolcano$FeatureNames,#Metabolite name
                                              selectLab = SelectLab,
                                              x = paste(x),
                                              y = paste(y),
                                              xlab  =xlab,
                                              ylab =ylab,
                                              pCutoff = xCutoff,
                                              FCcutoff = yCutoff,#Cut off Log2FC, automatically 2
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
                                              subtitle = Subtitle,
                                              caption = paste0("Total = ", nrow(InputVolcano), " ", Features),
                                              xlim =  c(min(InputVolcano[[x]][is.finite(InputVolcano[[x]] )])-0.2, max(InputVolcano[[x]][is.finite(InputVolcano[[x]])])+1.2),
                                              ylim = c(0,(ceiling(-log10(Reduce(min,InputVolcano[[y]]))))),
                                              cutoffLineType = "dashed",
                                              cutoffLineCol = "black",
                                              cutoffLineWidth = 0.5,
                                              legendLabels=c(paste(x," < |", yCutoff, "|"), paste(x," > |", yCutoff, "|"), paste(y, ' < ', xCutoff) , paste(y, ' < ', xCutoff,' & ',x," < |", yCutoff, "|")),
                                              legendPosition = LegendPos,
                                              legendLabSize = 9,
                                              legendIconSize =4,
                                              gridlines.major = FALSE,
                                              gridlines.minor = FALSE,
                                              drawConnectors = Connectors)
      #Add the theme
      if(is.null(Theme)==FALSE){
        Plot <- Plot+Theme
      }

      ## Store the plot in the 'plots' list
      PlotList[["Plot"]] <- Plot

      #Set the total heights and widths
      Plot_Sized <-  MetaProViz:::plotGrob_Volcano(Input=Plot, keyvals = keyvals, keyvalsshape = keyvalsshape, PlotName = PlotName, Subtitle = Subtitle)
      Plot <-Plot_Sized[[3]]
      Plot <- ggplot2::ggplot() +
        annotation_custom(Plot)
      Plot <-Plot + theme(panel.background = element_rect(fill = "transparent"))
      PlotList_adaptedGrid[["Plot_Sized"]] <- Plot

      #----- Save
      suppressMessages(suppressWarnings(
      MetaProViz:::SaveRes(InputList_DF=NULL,
                             InputList_Plot= list("Plot_Sized"= PlotList_adaptedGrid[["Plot_Sized"]]),
                             SaveAs_Table=NULL,
                             SaveAs_Plot=SaveAs_Plot,
                             FolderPath= Folder,
                             FileName= paste("Volcano_", PlotName, sep=""),
                             CoRe=FALSE,
                             PrintPlot=PrintPlot,
                             PlotHeight=Plot_Sized[[1]],
                             PlotWidth=Plot_Sized[[2]],
                             PlotUnit="cm")))
    }
  }
  return(invisible(list("Plot"=PlotList,"Plot_Sized" = PlotList_adaptedGrid)))
}




################################################################################################
### ### ### VizVolcano helper function: Internal Function for PlotSettings Compare ### ### ###
################################################################################################

#' Check input parameters
#'
#' @param InputData Passed to main function MetaProViz::VizVolcano()
#' @param InputData2 Passed to main function MetaProViz::VizVolcano()
#' @param SettingsFile_Metab Passed to main function MetaProViz::VizVolcano()
#' @param SettingsInfo Passed to main function MetaProViz::VizVolcano()
#' @param y \emph{Optional: } Passed to main function MetaProViz::VizVolcano() \strong{Default = "p.adj"}
#' @param x \emph{Optional: } Passed to main function MetaProViz::VizVolcano() \strong{Default = "Log2FC"}
#' @param PlotName \emph{Optional: } Passed to main function MetaProViz::VizVolcano() \strong{Default = ""}
#' @param xlab \emph{Optional: } Passed to main function MetaProViz::VizVolcano()  \strong{Default = NULL}
#' @param ylab \emph{Optional: } Passed to main function MetaProViz::VizVolcano() \strong{Default = NULL}
#' @param xCutoff \emph{Optional: } Passed to main function MetaProViz::VizVolcano() \strong{Default = 0.5}
#' @param ycutoff \emph{Optional: } Passed to main function MetaProViz::VizVolcano() \strong{Default = 0.05}
#' @param SelectLab \emph{Optional: } Passed to main function MetaProViz::VizVolcano() \strong{Default = ""}
#' @param Connectors \emph{Optional: } Passed to main function MetaProViz::VizVolcano() \strong{Default =  FALSE}
#' @param Subtitle \emph{Optional: } Passed to main function MetaProViz::VizVolcano() \strong{Default = ""}
#' @param ColorPalette Created in MetaProViz::VizVolcano() based on ColorPalette passed to main function MetaProViz::VizVolcano()
#' @param ShapePalette Created in MetaProViz::VizVolcano() based on ShapePalette passed to main function MetaProViz::VizVolcano()
#' @param Theme \emph{Optional: } Selection of theme for plot, e.g. theme_grey(). You can check for complete themes here: https://ggplot2.tidyverse.org/reference/ggtheme.html. \strong{Default = NULL}
#' @param Features \emph{Optional: } Name of the features that are plotted, e.g. "Metabolites", "RNA", "Proteins", "Genes", etc. \strong{Default = "Metabolites"}
#' @param ComparisonName Passed to main function MetaProViz::VizVolcano()
#' @param SaveAs_Plot Passed to main function MetaProViz::VizVolcano()
#' @param PrintPlot Passed to main function MetaProViz::VizVolcano()
#' @param Folder Created in MetaProViz::VizVolcano(). Path to the folder where files are saved.
#'
#'
#' @keywords Compare volcano plots
#' @noRd
#'
#'

VizVolcano_Compare <- function(InputData,
                               InputData2,
                               SettingsFile_Metab,
                               SettingsInfo,
                               y= "p.adj",
                               x= "Log2FC",
                               xlab= NULL,#"~Log[2]~FC"
                               ylab= NULL,#"~-Log[10]~p.adj"
                               xCutoff= 0.5,
                               yCutoff= 0.05,
                               Connectors=  FALSE,
                               SelectLab= "",
                               PlotName= "",
                               Subtitle= "",
                               ColorPalette,
                               ShapePalette,
                               Theme= NULL,
                               Features="Metabolites",
                               ComparisonName,
                               SaveAs_Plot,
                               PrintPlot,
                               Folder){

  #####################
  ##--- Check InputData
  if(is.data.frame(InputData2)==FALSE){
    if(paste(x) %in% colnames(InputData2)==FALSE | paste(y) %in% colnames(InputData2)==FALSE){
      stop("Check your InputData2. The column name of x and/or y does not exist in InputData2.")
    }
    }

  if(any(duplicated(row.names(InputData2)))==TRUE){
    stop("Duplicated row.names of InputData2, whilst row.names must be unique")
  }

  #Pass colours/shapes
  safe_colorblind_palette <- ColorPalette
  safe_shape_palette <- ShapePalette

  ##--- Prepare Input Data
  if(is.null(SettingsFile_Metab)==FALSE){
  InputData2 <- merge(x=SettingsFile_Metab%>%rownames_to_column("FeatureNames") , y=InputData2[, c(x, y)]%>%rownames_to_column("FeatureNames") , by="FeatureNames", all.y=TRUE)%>%
    na.omit()
  InputData[,"comparison"]  <- as.character(paste(ComparisonName[["InputData"]]))
  InputData2[,"comparison"]  <- as.character(paste(ComparisonName[["InputData2"]]))
  InputCompare  <- rbind(InputData,InputData2)

  }else{
   InputData2 <-  InputData2[, c(x, y)]%>%
    mutate(FeatureNames = rownames(InputData2))%>%
    na.omit()

   #Combine DFs and add appropriate column names
   InputData[,"comparison"]  <- as.character(paste(ComparisonName[["InputData"]]))
   InputData2[,"comparison"]  <- as.character(paste(ComparisonName[["InputData2"]]))
   InputCompare  <- rbind(InputData[,c("FeatureNames", x, y, "comparison")],InputData2[,c("FeatureNames", x, y, "comparison")])
  }



  #####################
  ##--- Plots
  if("individual" %in% names(SettingsInfo)==TRUE){
    # Create the list of individual plots that should be made:
    IndividualPlots <- unique(InputCompare$individual)

    PlotList <- list()#Empty list to store all the plots
    PlotList_adaptedGrid <- list()#Empty list to store all the plots

    for (i in IndividualPlots){
      InputVolcano <- subset(InputCompare, individual == paste(i))

      if(nrow(InputVolcano)>=1){
        #Prepare the colour scheme:
        if("color" %in% names(SettingsInfo)==TRUE){
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
        if("shape" %in% names(SettingsInfo)==TRUE & "color" %in% names(SettingsInfo)==FALSE){
          shape_select <- safe_shape_palette[1:length(unique(InputVolcano$shape))]

          keyvalsshape <- c()
          for(row in 1:nrow(InputVolcano)){
            sha <- shape_select[unique(InputVolcano$shape) %in% InputVolcano[row, "shape"]]
            names(sha) <- InputVolcano$shape[row]
            keyvalsshape <- c(keyvalsshape, sha)
          }
        } else if("shape" %in% names(SettingsInfo)==TRUE & "color" %in% names(SettingsInfo)==TRUE){
          #Here we have already used color from SettingsInfo and we need to use shape for the conditions
          message("For Plot_setting= `Consitions`we can only use colour or shape from SettingsFile_Metab. We ignore shape and use it to label the Comparison_name.")
          shape_select <- safe_shape_palette[1:length(unique(InputVolcano$comparison))]

          keyvalsshape <- c()
          for(row in 1:nrow(InputVolcano)){
            sha <- shape_select[unique(InputVolcano$comparison) %in% InputVolcano[row, "comparison"]]
            names(sha) <- InputVolcano$comparison[row]
            keyvalsshape <- c(keyvalsshape, sha)
          }
        } else if("shape" %in% names(SettingsInfo)==FALSE & "color" %in% names(SettingsInfo)==FALSE | "shape" %in% names(SettingsInfo)==FALSE & "color" %in% names(SettingsInfo)==TRUE){
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
                                                selectLab = SelectLab,
                                                x = paste(x),
                                                y = paste(y),
                                                xlab  =xlab,
                                                ylab =ylab,
                                                pCutoff = yCutoff,
                                                FCcutoff = xCutoff,#Cut off Log2FC, automatically 2
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
                                                subtitle = Subtitle,
                                                caption = paste0("Total = ", (nrow(InputVolcano)/2), " ", Features),
                                                xlim =  c(min(InputVolcano[[x]][is.finite(InputVolcano[[x]] )])-0.2, max(InputVolcano[[x]][is.finite(InputVolcano[[x]])])+1.2),
                                                ylim = c(0,(ceiling(-log10(Reduce(min,InputVolcano[[y]]))))),
                                                cutoffLineType = "dashed",
                                                cutoffLineCol = "black",
                                                cutoffLineWidth = 0.5,
                                                legendLabels=c(paste(x," < |", xCutoff, "|"), paste(x," > |", xCutoff, "|"), paste(y, ' < ', yCutoff) , paste(y, ' < ', yCutoff,' & ',x," < |", xCutoff, "|")),
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

        ## Store the plot in the 'plots' list
        PlotList[[i]] <- Plot

        #Set the total heights and widths
        PlotTitle <- paste(PlotName, ": ", i, sep="")
        Plot_Sized <-  MetaProViz:::plotGrob_Volcano(Input=Plot, keyvals = keyvals, keyvalsshape = keyvalsshape, PlotName = PlotTitle, Subtitle = Subtitle)
        Plot <-Plot_Sized[[3]]
        Plot <- ggplot2::ggplot() +
          annotation_custom(Plot)
        Plot <-Plot + theme(panel.background = element_rect(fill = "transparent"))

        #save plot and get rid of extra signs before saving
        cleaned_i <- gsub("[[:space:],/\\\\]", "-", i)#removes empty spaces and replaces /,\ with -
        PlotList_adaptedGrid[[cleaned_i]] <- Plot

        SaveList <- list()
        SaveList[[cleaned_i]] <- Plot
        #----- Save
        suppressMessages(suppressWarnings(
        MetaProViz:::SaveRes(InputList_DF=NULL,
                           InputList_Plot= SaveList,
                           SaveAs_Table=NULL,
                           SaveAs_Plot=SaveAs_Plot,
                           FolderPath= Folder,
                           FileName= paste("Volcano_",PlotName, sep=""),
                           CoRe=FALSE,
                           PrintPlot=PrintPlot,
                           PlotHeight=Plot_Sized[[1]],
                           PlotWidth=Plot_Sized[[2]],
                           PlotUnit="cm")))
      }
    }


  } else if("individual" %in% names(SettingsInfo)==FALSE){
    PlotList <- list()#Empty list to store all the plots
    PlotList_adaptedGrid <- list()#Empty list to store all the plots

    if(nrow(InputCompare)>=1){
      InputVolcano <- InputCompare
      #Prepare the colour scheme:
      if("color" %in% names(SettingsInfo)==TRUE){
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
      if("shape" %in% names(SettingsInfo)==TRUE & "color" %in% names(SettingsInfo)==FALSE){
        shape_select <- safe_shape_palette[1:length(unique(InputVolcano$shape))]

        keyvalsshape <- c()
        for(row in 1:nrow(InputVolcano)){
          sha <- shape_select[unique(InputVolcano$shape) %in% InputVolcano[row, "shape"]]
          names(sha) <- InputVolcano$shape[row]
          keyvalsshape <- c(keyvalsshape, sha)
        }
      } else if("shape" %in% names(SettingsInfo)==TRUE & "color" %in% names(SettingsInfo)==TRUE){
        #Here we have already used color from SettingsInfo and we need to use shape for the conditions
        message("For PlotSettings Comparison we can only use colour or shape from SettingsFile_Metab. Hence, we ignore shape and use it to label the ComparisonName.")
        shape_select <- safe_shape_palette[1:length(unique(InputVolcano$comparison))]

        keyvalsshape <- c()
        for(row in 1:nrow(InputVolcano)){
          sha <- shape_select[unique(InputVolcano$comparison) %in% InputVolcano[row, "comparison"]]
          names(sha) <- InputVolcano$comparison[row]
          keyvalsshape <- c(keyvalsshape, sha)
        }
      } else if("shape" %in% names(SettingsInfo)==FALSE & "color" %in% names(SettingsInfo)==FALSE | "shape" %in% names(SettingsInfo)==FALSE & "color" %in% names(SettingsInfo)==TRUE){
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
                                              selectLab = SelectLab,
                                              x = paste(x),
                                              y = paste(y),
                                              xlab  =xlab,
                                              ylab =ylab,
                                              pCutoff = xCutoff,
                                              FCcutoff = yCutoff,#Cut off Log2FC, automatically 2
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
                                              subtitle = Subtitle,
                                              caption = paste0("Total = ", (nrow(InputVolcano)/2)," ", Features),
                                              xlim =  c(min(InputVolcano[[x]][is.finite(InputVolcano[[x]] )])-0.2, max(InputVolcano[[x]][is.finite(InputVolcano[[x]])])+1.2),
                                              ylim = c(0,(ceiling(-log10(Reduce(min,InputVolcano[[y]]))))),
                                              cutoffLineType = "dashed",
                                              cutoffLineCol = "black",
                                              cutoffLineWidth = 0.5,
                                              legendLabels=c(paste(x," < |", yCutoff, "|"), paste(x," > |", yCutoff, "|"), paste(y, ' < ', xCutoff) , paste(y, ' < ', xCutoff,' & ',x," < |", yCutoff, "|")),
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

      ## Store the plot in the 'plots' list
      PlotList[["Plot"]] <- Plot

      #Set the total heights and widths
      Plot_Sized <-  MetaProViz:::plotGrob_Volcano(Input=Plot, keyvals = keyvals, keyvalsshape = keyvalsshape, PlotName = PlotName, Subtitle = Subtitle)
      Plot <-Plot_Sized[[3]]
      Plot <- ggplot2::ggplot() +
        annotation_custom(Plot)
      Plot <-Plot + theme(panel.background = element_rect(fill = "transparent"))
      PlotList_adaptedGrid[["Plot_Sized"]] <- Plot

      #----- Save
      suppressMessages(suppressWarnings(
      MetaProViz:::SaveRes(InputList_DF=NULL,
                             InputList_Plot= list("Plot_Sized"= PlotList_adaptedGrid[["Plot_Sized"]]),
                             SaveAs_Table=NULL,
                             SaveAs_Plot=SaveAs_Plot,
                             FolderPath= Folder,
                             FileName= paste("Volcano_", PlotName, sep=""),
                             CoRe=FALSE,
                             PrintPlot=PrintPlot,
                             PlotHeight=Plot_Sized[[1]],
                             PlotWidth=Plot_Sized[[2]],
                             PlotUnit="cm")))

    }
  }
  return(invisible(list("Plot"=PlotList,"Plot_Sized" = PlotList_adaptedGrid)))
}



################################################################################################
### ### ### VizVolcano helper function: Internal Function for PlotSettings PEA ### ### ###
################################################################################################

#' Check input parameters
#'
#' @param InputData Passed to main function MetaProViz::VizVolcano()
#' @param InputData2 Passed to main function MetaProViz::VizVolcano()
#' @param SettingsFile_Metab Passed to main function MetaProViz::VizVolcano()
#' @param SettingsInfo Passed to main function MetaProViz::VizVolcano()
#' @param y \emph{Optional: } Passed to main function MetaProViz::VizVolcano() \strong{Default = "p.adj"}
#' @param x \emph{Optional: } Passed to main function MetaProViz::VizVolcano() \strong{Default = "Log2FC"}
#' @param PlotName \emph{Optional: } Passed to main function MetaProViz::VizVolcano() \strong{Default = ""}
#' @param xlab \emph{Optional: } Passed to main function MetaProViz::VizVolcano()  \strong{Default = NULL}
#' @param ylab \emph{Optional: } Passed to main function MetaProViz::VizVolcano() \strong{Default = NULL}
#' @param xCutoff \emph{Optional: } Passed to main function MetaProViz::VizVolcano() \strong{Default = 0.5}
#' @param yCutoff \emph{Optional: } Passed to main function MetaProViz::VizVolcano() \strong{Default = 0.05}
#' @param SelectLab \emph{Optional: } Passed to main function MetaProViz::VizVolcano() \strong{Default = ""}
#' @param Connectors \emph{Optional: } Passed to main function MetaProViz::VizVolcano() \strong{Default =  FALSE}
#' @param Subtitle \emph{Optional: } Passed to main function MetaProViz::VizVolcano() \strong{Default = ""}
#' @param ColorPalette Created in MetaProViz::VizVolcano() based on ColorPalette passed to main function MetaProViz::VizVolcano()
#' @param ShapePalette Created in MetaProViz::VizVolcano() based on ShapePalette passed to main function MetaProViz::VizVolcano()
#' @param Theme \emph{Optional: } Selection of theme for plot, e.g. theme_grey(). You can check for complete themes here: https://ggplot2.tidyverse.org/reference/ggtheme.html. \strong{Default = NULL}
#' @param Features \emph{Optional: } Name of the features that are plotted, e.g. "Metabolites", "RNA", "Proteins", "Genes", etc. \strong{Default = "Metabolites"}
#' @param SaveAs_Plot Passed to main function MetaProViz::VizVolcano()
#' @param PrintPlot Passed to main function MetaProViz::VizVolcano()
#' @param Folder Created in MetaProViz::VizVolcano(). Path to the folder where files are saved.
#'
#' @keywords Volcano plots of pathway enrichment results
#' @noRd
#'
#'

VizVolcano_PEA <- function(InputData,
                           InputData2,
                           SettingsFile_Metab,
                           SettingsInfo,
                           y= "p.adj",
                           x= "Log2FC",
                           xlab= NULL,#"~Log[2]~FC"
                           ylab= NULL,#"~-Log[10]~p.adj"
                           xCutoff= 0.5,
                           yCutoff= 0.05,
                           Connectors=  FALSE,
                           SelectLab= "",
                           PlotName= "",
                           Subtitle= "",
                           ColorPalette,
                           ShapePalette,
                           Theme= NULL,
                           Features="Metabolites",
                           SaveAs_Plot,
                           PrintPlot,
                           Folder){
  #####################
  ##--- Check PEA settings
  if(is.vector(SettingsInfo)==FALSE){
    stop("You have chosen Settings=`PEA` that requires you to provide a vector for SettingsInfo.")
  }
  if(is.null(SettingsFile_Metab)==TRUE){
    stop("You have chosen Settings=`PEA` that requires you to provide a DF SettingsFile_Metab including the pathways used for the enrichment analysis.")
  }
  if(is.null(SettingsFile_Metab)==FALSE & is.null(SettingsFile_Metab)==FALSE){
    if("PEA_Feature" %in% names(SettingsInfo)==FALSE | "PEA_score" %in% names(SettingsInfo)==FALSE | "PEA_stat" %in% names(SettingsInfo)==FALSE | "PEA_Pathway" %in% names(SettingsInfo)==FALSE){
      stop("You have chosen Settings=`PEA` that requires you to provide a vector for SettingsInfo including `PEA_Feature`, `PEA_Pathway`, `PEA_stat` and `PEA_score`.")
    }
  }

  #Pass colours/shapes
  safe_colorblind_palette <- ColorPalette
  safe_shape_palette <- ShapePalette

  #Prepare data:
  InputData <- InputData%>%
    dplyr::rename("PEA_Feature"= !!SettingsInfo[["PEA_Feature"]])


  InputData2 <- InputData2%>%
    dplyr::rename("PEA_score"= !!SettingsInfo[["PEA_score"]],
                  "PEA_stat"= !!SettingsInfo[["PEA_stat"]],
                  "PEA_Pathway"= !!SettingsInfo[["PEA_Pathway"]])

  SettingsFile_Metab <- SettingsFile_Metab%>%
    dplyr::rename("PEA_Pathway"= !!SettingsInfo[["PEA_Pathway"]],
                  "PEA_Feature"= !!SettingsInfo[["PEA_Feature"]])

  #################
  ##--- Plot
  # Create the list of individual plots that should be made:
  IndividualPlots <- unique(InputData2$PEA_Pathway)

  PlotList <- list()#Empty list to store all the plots
  PlotList_adaptedGrid <- list()#Empty list to store all the plots

  for (i in IndividualPlots){
    InputData2_Select<- InputData2%>%
      filter(PEA_Pathway == paste(i)) #Select pathway we plot and use the score and stats

    SettingsFile_Metab_Select <- SettingsFile_Metab%>%
      filter(PEA_Pathway == paste(i))

    InputVolcano <-merge(SettingsFile_Metab_Select, InputData, by="PEA_Feature", all.x=TRUE)%>%
      distinct(PEA_Feature, .keep_all = TRUE) %>%
      filter(!is.na(!!sym(y)) & !is.na(!!sym(x)))

    if(nrow(InputVolcano)>=1){
      #Prepare the colour scheme:
      if("color" %in% names(SettingsInfo)==TRUE){
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
      if("shape" %in% names(SettingsInfo)==TRUE){
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

      if("color" %in% names(SettingsInfo)==FALSE & "shape" %in% names(SettingsInfo)==FALSE){
        LegendPos<- "none"
      }

      #Prepare the Plot:
      Plot<- EnhancedVolcano::EnhancedVolcano(InputVolcano,
                                              lab = InputVolcano$PEA_Feature,#Metabolite name
                                              selectLab = SelectLab,
                                              x = paste(x),
                                              y = paste(y),
                                              xlab  =xlab,
                                              ylab =ylab,
                                              pCutoff = yCutoff,
                                              FCcutoff = xCutoff,#Cut off Log2FC, automatically 2
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
                                              subtitle = paste(SettingsInfo[["PEA_score"]],"= ", InputData2_Select$PEA_score, ", ",SettingsInfo[["PEA_stat"]] , "= ", InputData2_Select$PEA_stat, sep=""),
                                              caption = paste0("Total = ", nrow(InputVolcano), " of ", nrow(SettingsFile_Metab_Select), " ", Features, " in pathway"),
                                              xlim =  c(min(InputVolcano[[x]][is.finite(InputVolcano[[x]] )])-0.2, max(InputVolcano[[x]][is.finite(InputVolcano[[x]])])+1.2),
                                              ylim = c(0,(ceiling(-log10(Reduce(min,InputVolcano[[y]]))))),
                                              cutoffLineType = "dashed",
                                              cutoffLineCol = "black",
                                              cutoffLineWidth = 0.5,
                                              legendLabels=c(paste(x," < |", xCutoff, "|"), paste(x," > |", xCutoff, "|"), paste(y, ' < ', yCutoff) , paste(y, ' < ', yCutoff,' & ',x," < |", xCutoff, "|")),
                                              legendPosition = LegendPos,
                                              legendLabSize = 7,
                                              legendIconSize =4,
                                              gridlines.major = FALSE,
                                              gridlines.minor = FALSE,
                                              drawConnectors = Connectors)
      #Add the theme
      if(is.null(Theme)==FALSE){
        Plot <- Plot+Theme
      }

      ## Store the plot in the 'plots' list
      PlotList[[i]] <- Plot

      #Set the total heights and widths
      PlotTitle <- paste(PlotName, ": ", i, sep="")
      Plot_Sized <-  MetaProViz:::plotGrob_Volcano(Input=Plot, keyvals = keyvals, keyvalsshape = keyvalsshape, PlotName = PlotTitle, Subtitle = Subtitle)
      Plot <-Plot_Sized[[3]]

      # First we want to convert the plot back into a ggplot object:
      Plot <- ggplot2::ggplot() +
        annotation_custom(Plot)
      Plot <-Plot + theme(panel.background = element_rect(fill = "transparent"))

      #save plot and get rid of extra signs before saving
      cleaned_i <- gsub("[[:space:],/\\\\]", "-", i)#removes empty spaces and replaces /,\ with -
      PlotList_adaptedGrid[[cleaned_i]] <- Plot
      SaveList <- list()
      SaveList[[cleaned_i]] <- Plot

      #----- Save
      suppressMessages(suppressWarnings(
        MetaProViz:::SaveRes(InputList_DF=NULL,
                             InputList_Plot= SaveList,
                             SaveAs_Table=NULL,
                             SaveAs_Plot=SaveAs_Plot,
                             FolderPath= Folder,
                             FileName= paste("Volcano_",PlotName, sep=""),
                             CoRe=FALSE,
                             PrintPlot=PrintPlot,
                             PlotHeight=Plot_Sized[[1]],
                             PlotWidth=Plot_Sized[[2]],
                             PlotUnit="cm")))

    }
  }

 return(invisible(list("Plot"=PlotList,"Plot_Sized" = PlotList_adaptedGrid)))
}


##############################################################
### ### ### Volcano helper function: Internal Function ### ### ###
##############################################################

#' @param Input This is the ggplot object generated within the VizVolcano function.
#' @param keyvals Generated in VizVolcano
#' @param keyvalsshape Generated in VizVolcano
#' @param PlotName Passed to VizVolcano
#' @param Subtitle
#'
#' @keywords Volcano helper function
#' @noRd
#'

plotGrob_Volcano <- function(Input, keyvals, keyvalsshape, PlotName, Subtitle){
  #------- Set the total heights and widths
  #we need ggplot_grob to edit the gtable of the ggplot object. Using this we can manipulate the gtable arguments directly.
  plottable<- ggplot2::ggplotGrob(Input) # Convert the plot to a gtable:  gtable::gtable_show_layout(plottable)
  if(is.null(keyvals)==TRUE & is.null(keyvalsshape)==TRUE){
    #-----widths
    plottable$widths[5] <- unit(6, "cm")#controls x-axis
    plottable$widths[c(3)] <- unit(2,"cm")#controls margins --> y-axis label is there
    plottable$widths[c(1,2,4)] <- unit(0,"cm")#controls margins --> not needed
    plottable$widths[c(6)] <- unit(1,"cm")#controls margins --> start Figure legend
    plottable$widths[c(10)] <- unit(5,"cm")#controls margins --> Figure legend
    plottable$widths[c(7,8,9,11)] <- unit(0,"cm")#controls margins --> not needed
    plot_widths <- 14

    if((PlotName=="" | Subtitle=="")==FALSE){#Check how much width is needed for the figure title/subtitle
      Titles <- c(PlotName, Subtitle)
      longest_title <- Titles[which.max(nchar(Titles))]
      character_count <- nchar(longest_title)
      Titles_width <- (character_count*0.25)+0.8
      if(Titles_width>plot_widths){#If the title needs more space than the plot offers:
        plottable$widths[11] <- unit(Titles_width-plot_widths,"cm")#controls margins --> start Figure legend
        plot_widths <- Titles_width
      }
    }

    #-----heigths
    plottable$heights[7] <- unit(8, "cm")#controls x-axis
    plottable$heights[c(8)] <- unit(1,"cm")#controls margins --> x-axis label
    plottable$heights[c(10)] <- unit(1.5,"cm")#controls margins --> Figure caption
    plottable$heights[c(9,11,12)] <- unit(0,"cm")#controls margins --> not needed

    if(PlotName=="" & Subtitle==""){
      plottable$heights[c(6)] <- unit(0.5,"cm")#controls margins --> Some space above the plot
      plottable$heights[c(1,2,3,4,5)] <- unit(0,"cm")#controls margins --> not needed
      plot_heights <- 11
    } else{
      plottable$heights[c(3)] <- unit(1,"cm")#controls margins --> PlotName and subtitle
      plottable$heights[c(2,4,5,6)] <- unit(0,"cm")#controls margins --> not needed
      plot_heights <-11.5
    }
  }else if(is.null(keyvals)==FALSE & is.null(keyvalsshape)==FALSE){
    #------- Legend heights
    Legend <- ggpubr::get_legend(Input) # Extract legend to adjust separately
    Legend_heights <- (round(as.numeric(Legend$heights[3]),1))+(round(as.numeric(Legend$heights[5]),1))

    #-----Plot widths
    plottable$widths[5] <- unit(6, "cm")#controls x-axis
    plottable$widths[c(3)] <- unit(2,"cm")#controls margins --> y-axis label is there
    plottable$widths[c(1,2,4)] <- unit(0,"cm")#controls margins --> not needed
    plottable$widths[c(6)] <- unit(1,"cm")#controls margins --> start Figure legend
    plottable$widths[c(7,8,10,11)] <- unit(0,"cm")#controls margins --> not needed

    Value <- round(as.numeric(plottable$widths[9]),1) #plottable$widths[9] is a <unit/unit_v2> object and we can extract the extract the numeric part
    plot_widths <- 9+Value

    if((PlotName=="" | Subtitle=="")==FALSE){#Check how much width is needed for the figure title/subtitle
      Titles <- c(PlotName, Subtitle)
      longest_title <- Titles[which.max(nchar(Titles))]
      character_count <- nchar(longest_title)
      Titles_width <- (character_count*0.25)+0.8
      if(Titles_width>plot_widths){#If the title needs more space than the plot offers:
        plottable$widths[11] <- unit(Titles_width-plot_widths,"cm")#controls margins --> start Figure legend
        plot_widths <- Titles_width
      }
    }

    #-----Plot heigths
    plottable$heights[7] <- unit(8, "cm")#controls x-axis
    plottable$heights[c(8)] <- unit(1,"cm")#controls margins --> x-axis label
    plottable$heights[c(10)] <- unit(1.5,"cm")#controls margins --> Figure caption
    plottable$heights[c(9,11)] <- unit(0,"cm")#controls margins --> not needed

    if(PlotName=="" & Subtitle==""){
      plottable$heights[c(6)] <- unit(0.5,"cm")#controls margins --> Some space above the plot
      plottable$heights[c(2,3,4,5)] <- unit(0,"cm")#controls margins --> not needed

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
    } else{#If we do have Title and or subtitle
      plottable$heights[c(3)] <- unit(1,"cm")#controls margins --> PlotName and subtitle
      plottable$heights[c(2,4,5,6)] <- unit(0,"cm")#controls margins --> not needed
      if(Legend_heights>11.5){#If the legend requires more heights than the Plot (excluding title space)
          Add <- (Legend_heights-11.5)/2
          plottable$heights[1] <- unit(Add,"cm")#controls margins --> Can be increased if Figure legend needs more space on the top
          plottable$heights[12] <- unit(Add,"cm")#controls margins --> Can be increased if Figure legend needs more space on the bottom
          plot_heights <- Legend_heights
        }else{
        plottable$heights[1] <- unit(0,"cm")#controls margins --> Can be increased if Figure legend needs more space on the top
        plottable$heights[12] <- unit(0,"cm")#controls margins --> Can be increased if Figure legend needs more space on the bottom
        plot_heights <- 11.5
      }
    }
  }else if(is.null(keyvals)==FALSE | is.null(keyvalsshape)==FALSE){
    #------- Legend heights
    Legend <- ggpubr::get_legend(Input) # Extract legend to adjust separately
    Legend_heights <- (round(as.numeric(Legend$heights[3]),1))

    #----- Plot widths
    plottable$widths[5] <- unit(6, "cm")#controls x-axis
    plottable$widths[c(3)] <- unit(2,"cm")#controls margins --> y-axis label is there
    plottable$widths[c(1,2,4)] <- unit(0,"cm")#controls margins --> not needed
    plottable$widths[c(6)] <- unit(1,"cm")#controls margins --> start Figure legend
    plottable$widths[c(7,8,10,11)] <- unit(0,"cm")#controls margins --> not needed

    Value <- round(as.numeric(plottable$widths[9]),1) #plottable$widths[9] is a <unit/unit_v2> object and we can extract the extract the numeric part
    plot_widths <- 9+Value

    if((PlotName=="" | Subtitle=="")==FALSE){#Check how much width is needed for the figure title/subtitle
      Titles <- c(PlotName, Subtitle)
      longest_title <- Titles[which.max(nchar(Titles))]
      character_count <- nchar(longest_title)
      Titles_width <- (character_count*0.25)+0.8
      if(Titles_width>plot_widths){#If the title needs more space than the plot offers:
        plottable$widths[11] <- unit(Titles_width-plot_widths,"cm")#controls margins --> start Figure legend
        plot_widths <- Titles_width
      }
    }

    #-----Plot heigths
    plottable$heights[7] <- unit(8, "cm")#controls x-axis
    plottable$heights[c(8)] <- unit(1,"cm")#controls margins --> x-axis label
    plottable$heights[c(10)] <- unit(1.5,"cm")#controls margins --> Figure caption
    plottable$heights[c(9,11)] <- unit(0,"cm")#controls margins --> not needed

    if(PlotName=="" & Subtitle==""){
      plottable$heights[c(6)] <- unit(0.5,"cm")#controls margins --> Some space above the plot
      plottable$heights[c(2,3,4,5)] <- unit(0,"cm")#controls margins --> not needed

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
    }else{#If we do have Title and or subtitle
      plottable$heights[c(3)] <- unit(1,"cm")#controls margins --> PlotName and subtitle
      plottable$heights[c(2,4,5,6)] <- unit(0,"cm")#controls margins --> not needed
      if(Legend_heights>11){#If the legend requires more heights than the Plot (excluding title space)
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
