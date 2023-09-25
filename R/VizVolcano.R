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
#' @param Plot_Settings \emph{Optional: } Choose between "Standard" (Input_data), "Compare" (plot two comparisons together Input_data and Input_data2) or "PEA" (Pathway Enrichment Analysis) \strong{Default = "Standard"}
#' @param Plot_SettingsInfo \emph{Optional: } NULL or Named vector including at least one of those three information for Plot_Settings="Standard" or "Compare": c(color="ColumnName_Plot_SettingsFile", shape= "ColumnName_Plot_SettingsFile", individual="ColumnName_Plot_SettingsFile"). For Plot_Settings="PEA" a named vector with c(PEA_Pathway="ColumnNameAdditionalInput_data", PEA_score="ColumnNameAdditionalInput_data", PEA_stat= "ColumnNameAdditionalInput_data", individual="Plot_SettingsFile), optionally you can additionally include c(color="ColumnName_Plot_SettingsFile", shape= "ColumnName_Plot_SettingsFile").\strong{Default = NULL}
#' @param Plot_SettingsFile \emph{Optional: } DF with column "Metabolite" including the Metabolite names (needs to match Metabolite names of Input_data) and other columns with required PlotSettingInfo. \strong{Default = NULL}
#' @param Input_data DF with column "Metabolite" including the Metabolite names, Log2FC, pvalue/padjusted values. Can also include additional columns with metadata usable for Plot_Setting_Info.
#' @param y \emph{Optional: } Column name including the values that should be used for y-axis. Usually this would include the p.adjusted value. \strong{Default = "p.adj"}
#' @param x \emph{Optional: } Column name including the values that should be used for x-axis. Usually this would include the Log2FC value. \strong{Default = "Log2FC"}
#' @param AdditionalInput_data \emph{Optional: } DF to compare to main Input_data with the same column names x and y (Plot_Settings="Compare") or Pathway enrichment analysis results (Plot_Settings="PEA"). \strong{Default = NULL}
#' @param OutputPlotName \emph{Optional: } String which is added to the output files of the plot. \strong{Default = ""}
#' @param Comparison_name \emph{Optional: } Named vector including those information about the two datasets that are compared on the plots when choosing Plot_Settings= "Compare". \strong{Default = c(Input_data="Cond1", AdditionalInput_data= "Cond2")}
#' @param xlab \emph{Optional: } String to replace x-axis label in plot. \strong{Default = NULL}
#' @param ylab \emph{Optional: } String to replace y-axis label in plot. \strong{Default = NULL}
#' @param pCutoff \emph{Optional: } Number of the desired p value cutoff for assessing significance. \strong{Default = 0.05}
#' @param FCcutoff \emph{Optional: } Number of the desired log fold change cutoff for assessing significance. \strong{Default = 0.5}
#' @param color_palette \emph{Optional: } Provide customiced color-palette in vector format. \strong{Default = NULL}
#' @param shape_palette \emph{Optional: } Provide customiced shape-palette in vector format. \strong{Default = NULL}
#' @param SelectLab \emph{Optional: } If set to NULL, feature labels will be plotted randomly. If vector is provided, e.g. c("MetaboliteName1", "MetaboliteName2"), selected names will be plotted. If set to default "", no feature names will be plotted. \strong{Default = ""}
#' @param Connectors \emph{Optional: } TRUE or FALSE for whether Connectors from names to points are to be added to the plot. \strong{Default =  FALSE}
#' @param Subtitle \emph{Optional: } \strong{Default = ""}
#' @param Theme \emph{Optional: } Selection of theme for plot, e.g. theme_grey(). You can check for complete themes here: https://ggplot2.tidyverse.org/reference/ggtheme.html. \strong{Default = NULL}
#'
#' @param Save_as_Plot \emph{Optional: } Select the file type of output plots. Options are svg, pdf, png or NULL. \strong{Default = "svg"}
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
                       xlab= NULL,#"~Log[2]~FC"
                       ylab= NULL,#"~-Log[10]~p.adj"
                       pCutoff= 0.05,
                       FCcutoff= 0.5,
                       color_palette= NULL,
                       shape_palette=NULL,
                       SelectLab= "",
                       Connectors=  FALSE,
                       Subtitle= "",
                       Theme= NULL,
                       Save_as_Plot= "svg"
                       ){

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

  ## ------------ Check Input files ----------- ##
  # 1. The input data:
  if("Metabolite" %in% names(Input_data)==FALSE){
    stop("Check input. Input_data must contain a column named `Metabolite` including the metabolite names.")
  }
  if(length(Input_data[duplicated(Input_data$Metabolite), "Metabolite"]) > 0){
    doublons <- as.character(Input_data[duplicated(Input_data$Metabolite), "Metabolite"])#number of duplications
    Input_data <-Input_data[!duplicated(Input_data$Metabolite),]#remove duplications
    warning("Input_data contained duplicates based on Metabolite! Dropping duplicate IDs and kept only the first entry. You had ", length(doublons), " duplicates. Note that you should do this before running VizVolcano.")
  }
  if( is.numeric(pCutoff)== FALSE |pCutoff > 1 | pCutoff < 0){
      stop("Check input. The selected pCutoff value should be numeric and between 0 and 1.")
    }
  if( is.numeric(FCcutoff)== FALSE  | FCcutoff < 0){
      stop("Check input. The selected pCutoff value should be numeric and between 0 and +oo.")
    }
  if(paste(x) %in% colnames(Input_data)==FALSE | paste(y) %in% colnames(Input_data)==FALSE){
    stop("Check your input. The column name of x and/ore y does not exist in Input_data.")
  }

  # 2. The Plot_settings: Plot_Settings, Plot_SettingInfo and Plot_SettingFile
  Plot_options <- c("Standard", "Compare", "PEA")
  if (Plot_Settings %in% Plot_options == FALSE){
    stop("Plot_Settings option is incorrect. The allowed options are the following: ",paste(Plot_options, collapse = ", "),"." )
  }
  if(is.null(Plot_SettingsInfo)==FALSE){
    if(is.vector(Plot_SettingsInfo)==FALSE){
      stop("Plot_SettingsInfo must be a named vector or NULL.")
    }
  }
  if(is.null(Plot_SettingsFile)==FALSE & "Metabolite" %in% names(Plot_SettingsFile)==FALSE){
    stop("Check input. Plot_SettingsFile must contain a column named `Metabolite` including the metabolite names.")
  }
  if(is.vector(Plot_SettingsInfo)==TRUE & is.null(Plot_SettingsFile)==TRUE){
    if("color" %in% names(Plot_SettingsInfo)==TRUE){#If Plot_SettingsFile=NULL, check if Input_data contain required columns
      if(Plot_SettingsInfo[["color"]] %in% names(Input_data)==FALSE){
        stop("You have chosen Plot_SettingsInfo option that requires you to provide a DF Plot_SettingsFile or the required columns to be present in Input_data.")
      }
    }
    if("shape" %in% names(Plot_SettingsInfo)==TRUE){
      if(Plot_SettingsInfo[["shape"]] %in% names(Input_data)==FALSE){
        stop("You have chosen Plot_SettingsInfo option that requires you to provide a DF Plot_SettingsFile or the required columns to be present in Input_data.")
      }
    }
    if("individual" %in% names(Plot_SettingsInfo)==TRUE){
      if(Plot_SettingsInfo[["individual" ]] %in% names(Input_data)==FALSE){
        stop("You have chosen Plot_SettingsInfo option that requires you to provide a DF Plot_SettingsFile or the required columns to be present in Input_data.")
      }
    }
    Plot_SettingsFile <- Input_data
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
    } else if(is.vector(Plot_SettingsInfo)==FALSE & is.null(Plot_SettingsFile)==FALSE){
      stop("Plot_SettingsInfo must be named vector or NULL.")
    }

  if(Plot_Settings=="PEA" & is.vector(Plot_SettingsInfo)==FALSE){
    stop("You have chosen Plot_Settings=`PEA` that requires you to provide a vector for Plot_SettingsInfo.")
  }else if(Plot_Settings=="PEA" & is.null(Plot_SettingsFile)==TRUE){
    stop("You have chosen Plot_Settings=`PEA` that requires you to provide a DF Plot_SettingsFile including the pathways used for the enrichment analysis.")
  } else if(Plot_Settings=="PEA" & is.null(Plot_SettingsFile)==FALSE & is.null(Plot_SettingsFile)==FALSE){
  if("individual" %in% names(Plot_SettingsInfo)==FALSE | "PEA_score" %in% names(Plot_SettingsInfo)==FALSE | "PEA_stat" %in% names(Plot_SettingsInfo)==FALSE | "PEA_Pathway" %in% names(Plot_SettingsInfo)==FALSE){
      stop("You have chosen Plot_Settings=`PEA` that requires you to provide a vector for Plot_SettingsInfo including `individual`, `PEA_Pathway`, `PEA_stat` and `PEA_score`.")
  }
  }

  #3. Select Input_data columns and Plot_SettingsFile columns
  if(is.null(Plot_SettingsFile)==FALSE){
    common_columns <- intersect(colnames(Input_data), colnames(Plot_SettingsFile))#check for overlapping names
    common_columns <- setdiff(common_columns, "Metabolite")#remove metabolites
    Plot_SettingsFile <- Plot_SettingsFile%>%#rename those column since they otherwise will cause issues when we merge the DFs later
      dplyr::rename_at(vars(common_columns), ~ paste0(., "_PlotSettingsFile"))
  }

  #4. AdditionalInput_data
  if(Plot_Settings=="Compare" & is.data.frame(AdditionalInput_data)==TRUE){
    if(paste(x) %in% colnames(AdditionalInput_data)==TRUE & paste(y) %in% colnames(AdditionalInput_data)==TRUE){
      AdditionalInput_data <- AdditionalInput_data %>%
        dplyr::rename("Log2FC"=paste(x),
                      "p.adj"=paste(y))
      } else{
        stop("Check your AdditionalInput_data. The column name of x and/or y does not exist in AdditionalInput_data.")
      }
     if("Metabolite" %in% names(AdditionalInput_data)==FALSE){
       stop("Check input. AdditionalInput_data must contain a column named `Metabolite` including the metabolite names.")
     }
     if(length(AdditionalInput_data[duplicated(AdditionalInput_data$Metabolite), "Metabolite"]) > 0){
       doublons <- as.character(AdditionalInput_data[duplicated(AdditionalInput_data$Metabolite), "Metabolite"])#number of duplications
       AdditionalInput_data <-AdditionalInput_data[!duplicated(AdditionalInput_data$Metabolite),]#remove duplications
       warning("AdditionalInput_data contained duplicates based on Metabolite! Dropping duplicate IDs and kept only the first entry. You had ", length(doublons), " duplicates. Note that you should do this before running VizVolcano.")
     }
     #Combine DFs and add appropriate column names
     Input_data[,"comparison"]  <- as.character(paste(Comparison_name[["Input_data"]]))
     AdditionalInput_data[,"comparison"]  <- as.character(paste(Comparison_name[["AdditionalInput_data"]]))
     Input_Comparison <- rbind(Input_data,AdditionalInput_data)
   } else if(Plot_Settings=="Comparison" & is.data.frame(AdditionalInput_data)==FALSE){
    stop("If Plot_Settings=`Comparison` you have to provide a DF for AdditionalInput_data.")
   }

   if(Plot_Settings=="PEA" & is.data.frame(AdditionalInput_data)==FALSE){
      stop("If Plot_Settings=`PEA` you have to provide a DF for AdditionalInput_data including the results of an enrichment analysis.")
     } else if(Plot_Settings=="PEA" & is.data.frame(AdditionalInput_data)==TRUE){
       AdditionalInput_data <- AdditionalInput_data%>%
        dplyr::rename("PEA_score"=paste(Plot_SettingsInfo[["PEA_score"]]),
                      "PEA_stat"=paste(Plot_SettingsInfo[["PEA_stat"]]),
                      "PEA_Pathway"=paste(Plot_SettingsInfo[["PEA_Pathway"]]))
       }

  # 5. Check other plot-specific parameters:
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
  if(is.null(SelectLab)==FALSE & is.vector(SelectLab)==FALSE){
      stop("Check input. SelectedLab mus be either NULL or a vector.")
  }
    if(is.logical(Connectors) == FALSE){
      stop("Check input. The Connectors value should be either = TRUE if connectors from names to points are to be added to the plot or =FALSE if not.")
    }
    if (!is.null(Save_as_Plot)) {
      Save_as_Plot_options <- c("svg","pdf", "png")
      if(Save_as_Plot %in% Save_as_Plot_options == FALSE){
        stop("Check input. The selected Save_as_Plot option is not valid. Please select one of the following: ",paste(Save_as_Plot_options,collapse = ", "),"." )
      }
    }

  ##theme
    if(is.null(Theme)==FALSE){
      Theme_options <- c("theme_grey()", "theme_gray()", "theme_bw()", "theme_linedraw()", "theme_light()", "theme_dark()", "theme_minimal()", "theme_classic()", "theme_void()", "theme_test()")
      if (Theme %in% Theme_options == FALSE){
      stop("Theme option is incorrect. You can check for complete themes here: https://ggplot2.tidyverse.org/reference/ggtheme.html. Options are the following: ",paste(Theme_options, collapse = ", "),"." )
    }
    }

  ## Rename the x and y lab if the information has been passed:
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

  ## ------------ Create Output folders ----------- ##
  if (!is.null(Save_as_Plot)) {
    name <- paste0("MetaProViz_Results_",Sys.Date())
    WorkD <- getwd()
    Results_folder <- file.path(WorkD, name)
    if (!dir.exists(Results_folder)) {dir.create(Results_folder)} # Make Results folder
    Results_folder_plots_Volcano_folder = file.path(Results_folder, "Volcano")  # This searches for a folder called "Preprocessing" within the "Results" folder in the current working directory and if its not found it creates one
    if (!dir.exists(Results_folder_plots_Volcano_folder)) {dir.create(Results_folder_plots_Volcano_folder)}  # check and create folder
  }
  ############################################################################################################
  ## ----------- Make the  plot based on the chosen parameters ------------ ##
  #####--- 1. Standard
  if(Plot_Settings=="Standard"){############################################################################################################
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

          if("color" %in% names(Plot_SettingsInfo)==TRUE ){
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
          if("shape" %in% names(Plot_SettingsInfo)==TRUE){
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

          if("color" %in% names(Plot_SettingsInfo)==FALSE & "shape" %in% names(Plot_SettingsInfo)==FALSE){
            LegendPos<- "none"
          }

          #Prepare the Plot:
          Plot<- EnhancedVolcano::EnhancedVolcano(InputVolcano,
                                                  lab = InputVolcano$Metabolite,#Metabolite name
                                                  selectLab = SelectLab,
                                                  x = paste(x),
                                                  y = paste(y),
                                                  xlab  =xlab,
                                                  ylab =ylab,
                                                  pCutoff = pCutoff,
                                                  FCcutoff = FCcutoff,#Cut off Log2FC, automatically 2
                                                  pointSize = 3,
                                                  labSize = 3,
                                                  axisLabSize = 10,
                                                  titleLabSize = 12,
                                                  subtitleLabSize = 11,
                                                  captionLabSize = 10,
                                                  colCustom = keyvals,
                                                  shapeCustom = keyvalsshape,
                                                  colAlpha = 1,
                                                  title= paste(OutputPlotName, ": ", i, sep=""),
                                                  subtitle = Subtitle,
                                                  caption = paste0("Total = ", nrow(InputVolcano), " Metabolites"),
                                                  xlim =  c(min(InputVolcano[[x]][is.finite(InputVolcano[[x]] )])-0.2, max(InputVolcano[[x]][is.finite(InputVolcano[[x]])])+1.2),
                                                  ylim = c(0,(ceiling(-log10(Reduce(min,InputVolcano[[y]]))))),
                                                  cutoffLineType = "dashed",
                                                  cutoffLineCol = "black",
                                                  cutoffLineWidth = 0.5,
                                                  legendLabels=c(paste(x," < |", FCcutoff, "|"), paste(x," > |", FCcutoff, "|"), paste(y, ' < ', pCutoff) , paste(y, ' < ', pCutoff,' & ',x," < |", FCcutoff, "|")),
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

          #Set the total heights and widths
          PlotTitle <- paste(OutputPlotName, ": ", i, sep="")
          Plot_Sized <-  MetaProViz:::plotGrob_Volcano(Input=Plot, keyvals = keyvals, keyvalsshape = keyvalsshape, OutputPlotName = PlotTitle, Subtitle = Subtitle)
          Plot <-Plot_Sized[[3]]
          Plot <- ggplot2::ggplot() +
            annotation_custom(Plot)

          #save plot and get rid of extra signs before saving
          cleaned_i <- gsub("[[:space:],/\\\\]", "-", i)#removes empty spaces and replaces /,\ with -
          if (!is.null(Save_as_Plot)) {
            if(OutputPlotName ==""){
              ggsave(file=paste(Results_folder_plots_Volcano_folder,"/", "Volcano_",cleaned_i, ".",Save_as_Plot, sep=""), plot=Plot, width=Plot_Sized[[2]], height=Plot_Sized[[1]], unit="cm")
            }else{
              ggsave(file=paste(Results_folder_plots_Volcano_folder,"/", "Volcano_", OutputPlotName, "_",cleaned_i, ".",Save_as_Plot, sep=""), plot=Plot, width=Plot_Sized[[2]], height=Plot_Sized[[1]], unit="cm")
            }
          }
          ## Store the plot in the 'plots' list
          PlotList[[cleaned_i]] <- Plot
          plot(Plot)
        }
      }
      #invisible(PlotList)
      } else if("individual" %in% names(Plot_SettingsInfo)==FALSE){
        if(is.null(Plot_SettingsFile)==FALSE){
          InputVolcano  <- merge(x=Plot_SettingsFile,y=Input_data, by="Metabolite", all.x=TRUE)%>%
            na.omit()
          }else{
            InputVolcano  <- Input_data
            }

        if(nrow(InputVolcano)>=1){
          if("color" %in% names(Plot_SettingsInfo)==TRUE ){
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
          if("shape" %in% names(Plot_SettingsInfo)==TRUE){
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

          if("color" %in% names(Plot_SettingsInfo)==FALSE & "shape" %in% names(Plot_SettingsInfo)==FALSE){
            LegendPos<- "none"
          }

          #Prepare the Plot:
          Plot<- EnhancedVolcano::EnhancedVolcano(InputVolcano,
                                                  lab = InputVolcano$Metabolite,#Metabolite name
                                                  selectLab = SelectLab,
                                                  x = paste(x),
                                                  y = paste(y),
                                                  xlab  =xlab,
                                                  ylab =ylab,
                                                  pCutoff = pCutoff,
                                                  FCcutoff = FCcutoff,#Cut off Log2FC, automatically 2
                                                  pointSize = 3,
                                                  labSize = 3,
                                                  axisLabSize = 10,
                                                  titleLabSize = 12,
                                                  subtitleLabSize = 11,
                                                  captionLabSize = 10,
                                                  colCustom = keyvals,
                                                  shapeCustom = keyvalsshape,
                                                  colAlpha = 1,
                                                  title= paste(OutputPlotName),
                                                  subtitle = Subtitle,
                                                  caption = paste0("Total = ", nrow(InputVolcano), " Metabolites"),
                                                  xlim =  c(min(InputVolcano[[x]][is.finite(InputVolcano[[x]] )])-0.2, max(InputVolcano[[x]][is.finite(InputVolcano[[x]])])+1.2),
                                                  ylim = c(0,(ceiling(-log10(Reduce(min,InputVolcano[[y]]))))),
                                                  cutoffLineType = "dashed",
                                                  cutoffLineCol = "black",
                                                  cutoffLineWidth = 0.5,
                                                  legendLabels=c(paste(x," < |", FCcutoff, "|"), paste(x," > |", FCcutoff, "|"), paste(y, ' < ', pCutoff) , paste(y, ' < ', pCutoff,' & ',x," < |", FCcutoff, "|")),
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

          #Set the total heights and widths
          Plot_Sized <-  MetaProViz:::plotGrob_Volcano(Input=Plot, keyvals = keyvals, keyvalsshape = keyvalsshape, OutputPlotName = OutputPlotName, Subtitle = Subtitle)
          Plot <-Plot_Sized[[3]]
          Plot <- ggplot2::ggplot() +
            annotation_custom(Plot)

          #save plot and get rid of extra signs before saving
          if (!is.null(Save_as_Plot)) {
            if(OutputPlotName ==""){
              ggsave(file=paste(Results_folder_plots_Volcano_folder,"/", "Volcano." ,Save_as_Plot, sep=""), plot=Plot, width=Plot_Sized[[2]], height=Plot_Sized[[1]], unit="cm")
            }else{
              ggsave(file=paste(Results_folder_plots_Volcano_folder,"/", "Volcano_", OutputPlotName, ".",Save_as_Plot, sep=""), plot=Plot, width=Plot_Sized[[2]], height=Plot_Sized[[1]], unit="cm")
            }
          }
          #Plot
          plot(Plot)
          #invisible(Plot)
        }
      }
  #####--- 2. Condition
  } else if(Plot_Settings=="Compare"){############################################################################################################
    if("individual" %in% names(Plot_SettingsInfo)==TRUE){
      # Create the list of individual plots that should be made:
      IndividualPlots <- Plot_SettingsFile[!duplicated(Plot_SettingsFile$individual),]
      IndividualPlots <- IndividualPlots$individual

      PlotList <- list()#Empty list to store all the plots

      for (i in IndividualPlots){
        Plot_SettingsFile_Select <- subset(Plot_SettingsFile, individual == paste(i))
        InputVolcano  <- merge(x=Plot_SettingsFile_Select,y=Input_Comparison, by="Metabolite", all.x=TRUE)%>%
          na.omit()

        if(nrow(InputVolcano)>=1){
          #Prepare the colour scheme:
          if("color" %in% names(Plot_SettingsInfo)==TRUE){
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
          if("shape" %in% names(Plot_SettingsInfo)==TRUE & "color" %in% names(Plot_SettingsInfo)==FALSE){
            shape_select <- safe_shape_palette[1:length(unique(InputVolcano$shape))]

            keyvalsshape <- c()
            for(row in 1:nrow(InputVolcano)){
              sha <- shape_select[unique(InputVolcano$shape) %in% InputVolcano[row, "shape"]]
              names(sha) <- InputVolcano$shape[row]
              keyvalsshape <- c(keyvalsshape, sha)
            }
          } else if("shape" %in% names(Plot_SettingsInfo)==TRUE & "color" %in% names(Plot_SettingsInfo)==TRUE){
            #Here we have already used color from Plot_SettingsInfo and we need to use shape for the conditions
            message("For Plot_setting= `Consitions`we can only use colour or shape from Plot_SettingsFile. We ignore shape and use it to label the Comparison_name.")
            shape_select <- safe_shape_palette[1:length(unique(InputVolcano$comparison))]

            keyvalsshape <- c()
            for(row in 1:nrow(InputVolcano)){
              sha <- shape_select[unique(InputVolcano$comparison) %in% InputVolcano[row, "comparison"]]
              names(sha) <- InputVolcano$comparison[row]
              keyvalsshape <- c(keyvalsshape, sha)
            }
          } else if("shape" %in% names(Plot_SettingsInfo)==FALSE & "color" %in% names(Plot_SettingsInfo)==FALSE){
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
                                                  lab = InputVolcano$Metabolite,#Metabolite name
                                                  selectLab = SelectLab,
                                                  x = paste(x),
                                                  y = paste(y),
                                                  xlab  =xlab,
                                                  ylab =ylab,
                                                  pCutoff = pCutoff,
                                                  FCcutoff = FCcutoff,#Cut off Log2FC, automatically 2
                                                  pointSize = 3,
                                                  labSize = 3,
                                                  axisLabSize = 10,
                                                  titleLabSize = 12,
                                                  subtitleLabSize = 11,
                                                  captionLabSize = 10,
                                                  colCustom = keyvals,
                                                  shapeCustom = keyvalsshape,
                                                  colAlpha = 1,
                                                  title= paste(OutputPlotName, ": ", i, sep=""),
                                                  subtitle = Subtitle,
                                                  caption = paste0("Total = ", (nrow(InputVolcano)/2), " Metabolites"),
                                                  xlim =  c(min(InputVolcano[[x]][is.finite(InputVolcano[[x]] )])-0.2, max(InputVolcano[[x]][is.finite(InputVolcano[[x]])])+1.2),
                                                  ylim = c(0,(ceiling(-log10(Reduce(min,InputVolcano[[y]]))))),
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

          #Set the total heights and widths
          PlotTitle <- paste(OutputPlotName, ": ", i, sep="")
          Plot_Sized <-  MetaProViz:::plotGrob_Volcano(Input=Plot, keyvals = keyvals, keyvalsshape = keyvalsshape, OutputPlotName = PlotTitle, Subtitle = Subtitle)
          Plot <-Plot_Sized[[3]]
          Plot <- ggplot2::ggplot() +
            annotation_custom(Plot)

          #save plot and get rid of extra signs before saving
          cleaned_i <- gsub("[[:space:],/\\\\]", "-", i)#removes empty spaces and replaces /,\ with -
          if (!is.null(Save_as_Plot)) {
            if(OutputPlotName ==""){
              ggsave(file=paste(Results_folder_plots_Volcano_folder,"/", "Volcano_",cleaned_i, ".",Save_as_Plot, sep=""), plot=Plot,  width=Plot_Sized[[2]], height=Plot_Sized[[1]], unit="cm")
            }else{
              ggsave(file=paste(Results_folder_plots_Volcano_folder,"/", "Volcano_", OutputPlotName, "_",cleaned_i, ".",Save_as_Plot, sep=""), plot=Plot,  width=Plot_Sized[[2]], height=Plot_Sized[[1]], unit="cm")
            }
          }
          ## Store the plot in the 'plots' list
          PlotList[[cleaned_i]] <- Plot
          plot(Plot)
        }
      }
      #invisible(PlotList)
    } else if("individual" %in% names(Plot_SettingsInfo)==FALSE){
      if(is.null(Plot_SettingsFile)==FALSE){
        InputVolcano  <- merge(x=Plot_SettingsFile,y=Input_Comparison, by="Metabolite", all.x=TRUE)%>%
          na.omit()
      }else{
        InputVolcano  <- Input_Comparison
      }

      if(nrow(InputVolcano)>=1){
        #Prepare the colour scheme:
        if("color" %in% names(Plot_SettingsInfo)==TRUE){
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
        if("shape" %in% names(Plot_SettingsInfo)==TRUE & "color" %in% names(Plot_SettingsInfo)==FALSE){
          shape_select <- safe_shape_palette[1:length(unique(InputVolcano$shape))]

          keyvalsshape <- c()
          for(row in 1:nrow(InputVolcano)){
            sha <- shape_select[unique(InputVolcano$shape) %in% InputVolcano[row, "shape"]]
            names(sha) <- InputVolcano$shape[row]
            keyvalsshape <- c(keyvalsshape, sha)
          }
        } else if("shape" %in% names(Plot_SettingsInfo)==TRUE & "color" %in% names(Plot_SettingsInfo)==TRUE){
          #Here we have already used color from Plot_SettingsInfo and we need to use shape for the conditions
          message("For Plot_setting= `Consitions`we can only use colour or shape from Plot_SettingsFile. We ignore shape and use it to label the Comparison_name.")
          shape_select <- safe_shape_palette[1:length(unique(InputVolcano$comparison))]

          keyvalsshape <- c()
          for(row in 1:nrow(InputVolcano)){
            sha <- shape_select[unique(InputVolcano$comparison) %in% InputVolcano[row, "comparison"]]
            names(sha) <- InputVolcano$comparison[row]
            keyvalsshape <- c(keyvalsshape, sha)
          }
        } else if("shape" %in% names(Plot_SettingsInfo)==FALSE & "color" %in% names(Plot_SettingsInfo)==FALSE){
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
                                                lab = InputVolcano$Metabolite,#Metabolite name
                                                selectLab = SelectLab,
                                                x = paste(x),
                                                y = paste(y),
                                                xlab  =xlab,
                                                ylab =ylab,
                                                pCutoff = pCutoff,
                                                FCcutoff = FCcutoff,#Cut off Log2FC, automatically 2
                                                pointSize = 3,
                                                labSize = 3,
                                                axisLabSize = 10,
                                                titleLabSize = 12,
                                                subtitleLabSize = 11,
                                                captionLabSize = 10,
                                                colCustom = keyvals,
                                                shapeCustom = keyvalsshape,
                                                colAlpha = 1,
                                                title= paste(OutputPlotName),
                                                subtitle = Subtitle,
                                                caption = paste0("Total = ", (nrow(InputVolcano)/2), " Metabolites"),
                                                xlim =  c(min(InputVolcano[[x]][is.finite(InputVolcano[[x]] )])-0.2, max(InputVolcano[[x]][is.finite(InputVolcano[[x]])])+1.2),
                                                ylim = c(0,(ceiling(-log10(Reduce(min,InputVolcano[[y]]))))),
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

        #Set the total heights and widths
        Plot_Sized <-  MetaProViz:::plotGrob_Volcano(Input=Plot, keyvals = keyvals, keyvalsshape = keyvalsshape, OutputPlotName = OutputPlotName, Subtitle = Subtitle)
        Plot <-Plot_Sized[[3]]
        Plot <- ggplot2::ggplot() +
          annotation_custom(Plot)

        #save plot and get rid of extra signs before saving i
        if (!is.null(Save_as_Plot)) {
          if(OutputPlotName ==""){
            ggsave(file=paste(Results_folder_plots_Volcano_folder,"/", "Volcano." ,Save_as_Plot, sep=""), plot=Plot,  width=Plot_Sized[[2]], height=Plot_Sized[[1]], unit="cm")
          }else{
            ggsave(file=paste(Results_folder_plots_Volcano_folder,"/", "Volcano_", OutputPlotName, ".",Save_as_Plot, sep=""), plot=Plot,  width=Plot_Sized[[2]], height=Plot_Sized[[1]], unit="cm")
          }
        }
        #Plot
        plot(Plot)
        #invisible(Plot)
      }
    }
  #####--- 3. PEA
  } else if(Plot_Settings=="PEA"){############################################################################################################
    # Create the list of individual plots that should be made:
    IndividualPlots <- Plot_SettingsFile[!duplicated(Plot_SettingsFile$individual),]
    IndividualPlots <- IndividualPlots$individual

    PlotList <- list()#Empty list to store all the plots

    for (i in IndividualPlots){
      Plot_SettingsFile_Select <- subset(Plot_SettingsFile, individual == paste(i))
      InputVolcano  <- merge(x=Plot_SettingsFile_Select,y=Input_data, by="Metabolite", all.x=TRUE)%>%
        na.omit()

      AdditionalInput_data_Select<- subset(AdditionalInput_data, PEA_Pathway == paste(i)) #Select pathway we plot and use the score and stats

      if(nrow(InputVolcano)>=1){
        #Prepare the colour scheme:
        if("color" %in% names(Plot_SettingsInfo)==TRUE){
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
        if("shape" %in% names(Plot_SettingsInfo)==TRUE){
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

        if("color" %in% names(Plot_SettingsInfo)==FALSE & "shape" %in% names(Plot_SettingsInfo)==FALSE){
          LegendPos<- "none"
        }

        #Prepare the Plot:
        Plot<- EnhancedVolcano::EnhancedVolcano(InputVolcano,
                                                lab = InputVolcano$Metabolite,#Metabolite name
                                                selectLab = SelectLab,
                                                x = paste(x),
                                                y = paste(y),
                                                xlab  =xlab,
                                                ylab =ylab,
                                                pCutoff = pCutoff,
                                                FCcutoff = FCcutoff,#Cut off Log2FC, automatically 2
                                                pointSize = 3,
                                                labSize = 3,
                                                axisLabSize = 10,
                                                titleLabSize = 12,
                                                subtitleLabSize = 11,
                                                captionLabSize = 10,
                                                colCustom = keyvals,
                                                shapeCustom = keyvalsshape,
                                                colAlpha = 1,
                                                title= paste(OutputPlotName, ": ", i, sep=""),
                                                subtitle = paste(Plot_SettingsInfo[["PEA_score"]],"= ", AdditionalInput_data_Select$PEA_score, ", ",Plot_SettingsInfo[["PEA_stat"]] , "= ", AdditionalInput_data_Select$PEA_stat, sep=""),
                                                caption = paste0("Total = ", nrow(InputVolcano), " of ", nrow(Plot_SettingsFile_Select), " metabolites in pathway"),
                                                xlim =  c(min(InputVolcano[[x]][is.finite(InputVolcano[[x]] )])-0.2, max(InputVolcano[[x]][is.finite(InputVolcano[[x]])])+1.2),
                                                ylim = c(0,(ceiling(-log10(Reduce(min,InputVolcano[[y]]))))),
                                                cutoffLineType = "dashed",
                                                cutoffLineCol = "black",
                                                cutoffLineWidth = 0.5,
                                                legendLabels=c(paste(x," < |", FCcutoff, "|"), paste(x," > |", FCcutoff, "|"), paste(y, ' < ', pCutoff) , paste(y, ' < ', pCutoff,' & ',x," < |", FCcutoff, "|")),
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

        #Set the total heights and widths
        PlotTitle <- paste(OutputPlotName, ": ", i, sep="")
        Plot_Sized <-  MetaProViz:::plotGrob_Volcano(Input=Plot, keyvals = keyvals, keyvalsshape = keyvalsshape, OutputPlotName = PlotTitle, Subtitle = Subtitle)
        Plot <-Plot_Sized[[3]]

        # First we want to convert the plot back into a ggplot object:
        Plot <- ggplot2::ggplot() +
          annotation_custom(Plot)

        #save plot and get rid of extra signs before saving
        cleaned_i <- gsub("[[:space:],/\\\\]", "-", i)#removes empty spaces and replaces /,\ with -
        if (!is.null(Save_as_Plot)) {
          if(OutputPlotName ==""){
            ggsave(file=paste(Results_folder_plots_Volcano_folder,"/", "Volcano_",cleaned_i, ".",Save_as_Plot, sep=""), plot=Plot, width=Plot_Sized[[2]], height=Plot_Sized[[1]], unit="cm")
          }else{
            ggsave(file=paste(Results_folder_plots_Volcano_folder,"/", "Volcano_", OutputPlotName, "_",cleaned_i, ".",Save_as_Plot, sep=""), plot=Plot, width=Plot_Sized[[2]], height=Plot_Sized[[1]], unit="cm")
          }
        }
        ## Store the plot in the 'plots' list
        PlotList[[cleaned_i]] <- Plot
       # plot(Plot)
      }
    }
    #invisible(PlotList)
  }
  return(invisible(Plot))
}


##############################################################
### ### ### Volcano helper function: Internal Function ### ### ###
##############################################################

#' @param Input This is the ggplot object generated within the VizVolcano function.
#' @param keyvals Generated in VizVolcano
#' @param keyvalsshape Generated in VizVolcano
#' @param OutputPlotName Passed to VizVolcano
#' @param Subtitle
#'
#' @keywords Volcano helper function
#' @noRd
#'

plotGrob_Volcano <- function(Input, keyvals, keyvalsshape, OutputPlotName, Subtitle){
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

    if((OutputPlotName=="" | Subtitle=="")==FALSE){#Check how much width is needed for the figure title/subtitle
      Titles <- c(OutputPlotName, Subtitle)
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

    if(OutputPlotName=="" & Subtitle==""){
      plottable$heights[c(6)] <- unit(0.5,"cm")#controls margins --> Some space above the plot
      plottable$heights[c(1,2,3,4,5)] <- unit(0,"cm")#controls margins --> not needed
      plot_heights <- 11
    } else{
      plottable$heights[c(3)] <- unit(1,"cm")#controls margins --> OutputPlotName and subtitle
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

    if((OutputPlotName=="" | Subtitle=="")==FALSE){#Check how much width is needed for the figure title/subtitle
      Titles <- c(OutputPlotName, Subtitle)
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

    if(OutputPlotName=="" & Subtitle==""){
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
      plottable$heights[c(3)] <- unit(1,"cm")#controls margins --> OutputPlotName and subtitle
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

    if((OutputPlotName=="" | Subtitle=="")==FALSE){#Check how much width is needed for the figure title/subtitle
      Titles <- c(OutputPlotName, Subtitle)
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

    if(OutputPlotName=="" & Subtitle==""){
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
      plottable$heights[c(3)] <- unit(1,"cm")#controls margins --> OutputPlotName and subtitle
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
