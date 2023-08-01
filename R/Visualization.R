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

#' @param Plot_SettingsInfo \emph{Optional: } NULL or Named vector including at least one of those three information : c(color="ColumnName_Plot_SettingsFile", shape= "ColumnName_Plot_SettingsFile"). \strong{Default = NULL}
#' @param Plot_SettingsFile \emph{Optional: } DF which contains information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames. Column "Conditions" with information about the sample conditions (e.g. "N" and "T" or "Normal" and "Tumor"), can be used for feature filtering and colour coding in the PCA. Column "AnalyticalReplicate" including numerical values, defines technical repetitions of measurements, which will be summarised. Column "BiologicalReplicates" including numerical values. Please use the following names: "Conditions", "Biological_Replicates", "Analytical_Replicates".\strong{Default = NULL}
#' @param Input_data DF with a column "UniqueID" with unique sample identifiers  and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected. includes experimental design and outlier column.
#' @param color_palette \emph{Optional: } Provide customiced color-palette in vector format. \strong{Default = NULL}
#' @param shape_palette \emph{Optional: } Provide customiced shape-palette in vector format. \strong{Default = NULL}
#' @param Show_Loadings  \emph{Optional: } TRUE or FALSE for whether PCA loadings are also plotted on the PCA (biplot) \strong{Default = FALSE}
#' @param Scaling  \emph{Optional: } TRUE or FALSE for whether a data scaling is used \strong{Default = TRUE}
#' @param Theme \emph{Optional: } Selection of theme for plots from ggplot2. \strong{Default = NULL}
#' @param OutputPlotName \emph{Optional: } String which is added to the output files of the PCA \strong{Default = ""}
#' @param Save_as_Plot \emph{Optional: } Select the file type of output plots. Options are svg, png, pdf. \strong{Default = svg}
#' @param Save_as_Results \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt" \strong{default: "csv"}
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
                   OutputPlotName= '',
                   Save_as_Plot = "svg",
                   Save_as_Results = "csv"
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
  if("UniqueID" %in% names(Input_data)==FALSE){
    stop("Check input. Input_data must contain a column named `UniqueID` including unique sample names.")
  }else if(length(Input_data[duplicated(Input_data$UniqueID), "UniqueID"]) > 0){
    stop("Input_data contained duplicates based on UniqueID!")
  } else{
    Input_data_m<- Input_data%>%
      column_to_rownames("UniqueID")
    Test_num <- apply(Input_data_m, 2, function(x) is.numeric(x))
    if((any(Test_num) ==  FALSE) ==  TRUE){
      stop("Input_data needs to be of class numeric in all columns (except column `UniqueID`).")
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

  Save_as_Plot_options <- c("svg","pdf", "png")
  if(Save_as_Plot %in% Save_as_Plot_options == FALSE){
    stop("Check input. The selected Save_as_Plot option is not valid. Please select one of the folowing: ",paste(Save_as_Plot_options,collapse = ", "),"." )
  }

  ## ------------ Create Output folders ----------- ##
  name <- paste0("MetaProViz_Results_",Sys.Date())
  WorkD <- getwd()
  Results_folder <- file.path(WorkD, name)
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)} # Make Results folder
  Results_folder_plots_PCA_folder = file.path(Results_folder, "PCA")  # This searches for a folder called "Preprocessing" within the "Results" folder in the current working directory and if its not found it creates one
  if (!dir.exists(Results_folder_plots_PCA_folder)) {dir.create(Results_folder_plots_PCA_folder)}  # check and create folder

  ############################################################################################################
  ## ----------- Make the  plot based on the choosen parameters ------------ ##
  if(is.null(Plot_SettingsFile)==FALSE){
    InputPCA  <- merge(x=Plot_SettingsFile, y=Input_data, by="UniqueID", all.y=TRUE)%>%
      column_to_rownames("UniqueID")
  }else{
    InputPCA  <- Input_data%>%
      column_to_rownames("UniqueID")
  }

  #Prepare the color scheme:
  if("color" %in% names(Plot_SettingsInfo)==TRUE){
    #color that will be used
    color_select <- safe_colorblind_palette[1:length(unique(InputPCA$color))]

    # numeric scale or continuous
    if(is.numeric(InputPCA$color) == TRUE | is.integer(InputPCA$color) == TRUE){
      if(length(unique(InputPCA$color)) > 4){ # change this to change the number after which color becomes from distinct to continuous
        InputPCA$color <- as.numeric(InputPCA$color)
      }else{
        InputPCA$color <- as.factor(InputPCA$color)
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

  #Make the plot:
  PCA <- autoplot(prcomp(as.matrix(Input_data_m, scale. = as.logical(Scaling))),
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

  plot(PCA)

  #Prepare output DF
  loading_data <- prcomp(as.matrix(Input_data_m, scale. = as.logical(Scaling)))

  loading_data_table <-as.data.frame(loading_data$rotation)
  loading_data_table <- loading_data_table
  loading_data_table <- tibble::rownames_to_column(loading_data_table, "Metabolite")

  # Save output
  if(OutputPlotName ==""){
    write.csv(loading_data_table, paste(Results_folder_plots_PCA_folder,"/Loadings.csv", sep=""))
  }else{
    write.csv(loading_data_table, paste(Results_folder_plots_PCA_folder,"/", OutputPlotName, "_Loadings.csv", sep=""))
  }

  if (Save_as_Results == "xlsx"){
    if(OutputPlotName ==""){
      xlsPCA <-  paste(Results_folder_plots_PCA_folder,"/Loadings.xlsx", sep="")
      writexl::write_xlsx(loading_data_table,xlsPCA, col_names = TRUE)
    }else{
      xlsPCA <-  paste(Results_folder_plots_PCA_folder,"/", OutputPlotName, "_Loadings.xlsx", sep="")
      writexl::write_xlsx(loading_data_table,xlsPCA, col_names = TRUE)
    }
  }else if (Save_as_Results == "csv"){
    if(OutputPlotName ==""){
      csvPCA <-  paste(Results_folder_plots_PCA_folder,"/Loadings.csv", sep="")
      write.csv(loading_data_table,csvPCA)
    }else{
      csvPCA <-  paste(Results_folder_plots_PCA_folder,"/", OutputPlotName, "_Loadings.csv", sep="")
      write.csv(loading_data_table,csvPCA)
    }
  }else if (Save_as_Results == "txt"){
    if(OutputPlotName ==""){
      txtPCA <-  paste(Results_folder_plots_PCA_folder,"/Loadings.txt", sep="")
      write.csv(loading_data_table,txtPCA)
    }else{
      txtPCA <-  paste(Results_folder_plots_PCA_folder,"/", OutputPlotName, "_Loadings.txt", sep="")
      write.csv(loading_data_table,txtPCA)
    }
  }

  if(OutputPlotName ==""){
    ggsave(file=paste(Results_folder_plots_PCA_folder,"/", "PCA", OutputPlotName, ".",Save_as_Plot, sep=""), plot=PCA, width=10, height=8)
  }else{
    ggsave(file=paste(Results_folder_plots_PCA_folder,"/", "PCA_", OutputPlotName, ".",Save_as_Plot, sep=""), plot=PCA, width=10, height=8)
  }
}



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
#' @param Output_Name \emph{Optional: } String which is added to the output files of the plot. \strong{Default = ""}
#'
#' @param Comparison_name \emph{Optional: } Named vector including those information about the two datasets that are compared on the plots when choosing Plot_Settings= "Compare". \strong{Default = c(Input_data="Cond1", AdditionalInput_data= "Cond2")}
#' @param xlab \emph{Optional: } String to replace x-axis label in plot. \strong{Default = NULL}
#' @param ylab \emph{Optional: } String to replace y-axis label in plot. \strong{Default = NULL}
#' @param pCutoff \emph{Optional: } Number of the desired p value cutoff for assessing significance. \strong{Default = 0.05}
#' @param FCcutoff \emph{Optional: } Number of the desired log fold change cutoff for assessing significance. \strong{Default = 0.5}
#' @param Legend \emph{Optional: } Legend=="Pie" will plot a PieChart as the legend for color or Legend="Standard, plot the standard legend for color. \strong{Default = "Standard"}
#' @param color_palette \emph{Optional: } Provide customiced color-palette in vector format. \strong{Default = NULL}
#' @param shape_palette \emph{Optional: } Provide customiced shape-palette in vector format. \strong{Default = NULL}
#' @param Connectors \emph{Optional: } TRUE or FALSE for whether Connectors from names to points are to be added to the plot. \strong{Default =  FALSE}
#' @param Subtitle \emph{Optional: } \strong{Default = ""}
#' @param Theme \emph{Optional: } Selection of theme for plot. \strong{Default = NULL}
#'
#' @param Save_as_Plot \emph{Optional: } Select the file type of output plots. Options are svg, pdf or png. \strong{Default = "svg"}
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
                       Legend="Standard",
                       color_palette= NULL,
                       shape_palette=NULL,
                       Connectors=  FALSE,
                       Subtitle= "",
                       Theme= NULL,
                       Save_as_Plot= "svg"
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

   #3. AdditionalInput_data
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

  # 4. Check other plot-specific parameters:
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
    if(is.logical(Connectors) == FALSE){
      stop("Check input. The Connectors value should be either = TRUE if connectors from names to points are to be added to the plot or =FALSE if not.")
    }
    Save_as_Plot_options <- c("svg","pdf", "png")
    if(Save_as_Plot %in% Save_as_Plot_options == FALSE){
      stop("Check input. The selected Save_as_Plot option is not valid. Please select one of the following: ",paste(Save_as_Plot_options,collapse = ", "),"." )
    }

  #Legend=="Pie" or Legend="Standard --> If color is not provided and Legend=="Pie" this will be ignored! --> change parameter and give warning!
  #outputplotname
  #theme

  # Rename the x and y lab if the information has been passed:
  if(is.null(xlab)==TRUE){#use column name of x provided by user
    xlab <- bquote(.(as.symbol(x)))
    }else if(is.null(xlab)==FALSE){
    xlab <- bquote(.(as.symbol(xlab)))
    }

  if(is.null(xlab)==TRUE){#use column name of x provided by user
    ylab <- bquote(.(as.symbol(y)))
    }else if(is.null(ylab)==FALSE){
      ylab <- bquote(.(as.symbol(ylab)))
      }

  ## ------------ Create Output folders ----------- ##
  name <- paste0("MetaProViz_Results_",Sys.Date())
  WorkD <- getwd()
  Results_folder <- file.path(WorkD, name)
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)} # Make Results folder
  Results_folder_plots_Volcano_folder = file.path(Results_folder, "Volcano")  # This searches for a folder called "Preprocessing" within the "Results" folder in the current working directory and if its not found it creates one
  if (!dir.exists(Results_folder_plots_Volcano_folder)) {dir.create(Results_folder_plots_Volcano_folder)}  # check and create folder

  ############################################################################################################
  ## ----------- Make the  plot based on the choosen parameters ------------ ##
  #Check the plot type: Comparison/Standard/GSE

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

          #Prepare the colour scheme:
          if("color" %in% names(Plot_SettingsInfo)==TRUE & Legend=="Pie"){
            #Prepare Summary to make PieChart
            color_select <- safe_colorblind_palette[1:length(unique(InputVolcano$color))]

            Summary <- InputVolcano %>%
              group_by(InputVolcano$color) %>%
              summarise(percent = round(100 * n() / nrow(InputVolcano))) %>%
              mutate(csum = rev(cumsum(rev(percent))),
                     pos = percent/2 + lead(csum, 1),
                     pos = if_else(is.na(pos), percent/2, pos))%>%
              rename("Group"=1)
            Summary$Label <- paste(Summary$Group, " (", Summary$percent, "%)")
            Summary$Palette <- color_select

            #Make PieChart
            PieChart <- ggplot(Summary, aes(x="", y=percent, fill=Label))+
              geom_col(width = 1, color = 1, alpha=0.7) +
              coord_polar(theta = "y") +
              scale_fill_manual(values=Summary$Palette)+
              ggrepel::geom_label_repel(data = Summary,
                               aes(y = pos, label = paste0(Group, " (", percent, "%)")),
                               size = 3.5, nudge_x = 1, show.legend = FALSE, fill = "white") +
              guides(fill = guide_legend(title = "Group")) +
              theme_void()+
              theme(legend.position = "none")

            #Assign LegendParameter for EnahncedVolcano:
            LegendPos<- "none"
          } else if(Legend=="Standard"){
            LegendPos<- "right"
          }
          if("color" %in% names(Plot_SettingsInfo)==TRUE){
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
            shape_select <- safe_shape_palette[1:length(unique(InputVolcano$shape))]

            keyvalsshape <- c()
            for(row in 1:nrow(InputVolcano)){
              sha <- shape_select[unique(InputVolcano$shape) %in% InputVolcano[row, "shape"]]
              names(sha) <- InputVolcano$shape[row]
              keyvalsshape <- c(keyvalsshape, sha)
            }
            #Assign LegendParameter for EnahncedVolcano:
            LegendPos<- "right"
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
                                                  title= paste(OutputPlotName, ": ", i, sep=""),
                                                  subtitle = Subtitle,
                                                  caption = paste0("total = ", nrow(InputVolcano), " Metabolites"),
                                                  xlim =  c(min(InputVolcano$Log2FC[is.finite(InputVolcano$Log2FC )])-0.2, max(InputVolcano$Log2FC[is.finite(InputVolcano$Log2FC )])+1.2),
                                                  ylim = c(0,(ceiling(-log10(Reduce(min,InputVolcano$p.adj))))),
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

          #Add PieChart
          if("color" %in% names(Plot_SettingsInfo)==TRUE & Legend=="Pie"){
            Plot <- Plot+ annotation_custom(
              grob = ggplotGrob(PieChart),
              xmin = (max(InputVolcano$Log2FC[is.finite(InputVolcano$Log2FC )])+0.2), xmax =(max(InputVolcano$Log2FC[is.finite(InputVolcano$Log2FC )])+1.5),
              ymin = ((ceiling(-log10(Reduce(min,InputVolcano$p.adj))))+0.5), ymax =((ceiling(-log10(Reduce(min,InputVolcano$p.adj))))-1.5)
            )
          }

          #save plot and get rid of extra signs before saving
          cleaned_i <- gsub("[[:space:],/\\\\]", "-", i)#removes empty spaces and replaces /,\ with -
          if(OutputPlotName ==""){
            ggsave(file=paste(Results_folder_plots_Volcano_folder,"/", "Volcano_",cleaned_i, ".",Save_as_Plot, sep=""), plot=Plot, width=8, height=6)
          }else{
            ggsave(file=paste(Results_folder_plots_Volcano_folder,"/", "Volcano_", OutputPlotName, "_",cleaned_i, ".",Save_as_Plot, sep=""), plot=Plot, width=8, height=6)
          }
          ## Store the plot in the 'plots' list
          PlotList[[cleaned_i]] <- Plot
          plot(Plot)
        }
      }
      # Return PlotList into the environment to enable the user to view the plots directly
      #assign("VolcanoPlots", PlotList, envir=.GlobalEnv)
      # Combine plots into a single plot using facet_grid or patchwork::wrap_plots
      Return <- PlotList
      } else if("individual" %in% names(Plot_SettingsInfo)==FALSE){
        if(is.null(Plot_SettingsFile)==FALSE){
          InputVolcano  <- merge(x=Plot_SettingsFile,y=Input_data, by="Metabolite", all.x=TRUE)%>%
            na.omit()
          }else{
            InputVolcano  <- Input_data
            }

        if(nrow(InputVolcano)>=1){
          if("color" %in% names(Plot_SettingsInfo)==TRUE & Legend=="Pie"){
            #Prepare Summary to make PieChart
            color_select <- safe_colorblind_palette[1:length(unique(InputVolcano$color))]

            Summary <- InputVolcano %>%
              group_by(InputVolcano$color) %>%
              summarise(percent = round(100 * n() / nrow(InputVolcano))) %>%
              mutate(csum = rev(cumsum(rev(percent))),
                     pos = percent/2 + lead(csum, 1),
                     pos = if_else(is.na(pos), percent/2, pos))%>%
              rename("Group"=1)
            Summary$Label <- paste(Summary$Group, " (", Summary$percent, "%)")
            Summary$Palette <- color_select

            #Make PieChart
            PieChart <- ggplot(Summary, aes(x="", y=percent, fill=Label))+
              geom_col(width = 1, color = 1, alpha=0.7) +
              coord_polar(theta = "y") +
              scale_fill_manual(values=Summary$Palette)+
              ggrepel::geom_label_repel(data = Summary,
                                        aes(y = pos, label = paste0(Group, " (", percent, "%)")),
                                        size = 3.5, nudge_x = 1, show.legend = FALSE, fill = "white") +
              guides(fill = guide_legend(title = "Group")) +
              theme_void()+
              theme(legend.position = "none")

            #Assign LegendParameter for EnahncedVolcano:
            LegendPos<- "none"
          } else if(Legend=="Standard"){
            LegendPos<- "right"
          }
          if("color" %in% names(Plot_SettingsInfo)==TRUE ){
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
            shape_select <- safe_shape_palette[1:length(unique(InputVolcano$shape))]

            keyvalsshape <- c()
            for(row in 1:nrow(InputVolcano)){
              sha <- shape_select[unique(InputVolcano$shape) %in% InputVolcano[row, "shape"]]
              names(sha) <- InputVolcano$shape[row]
              keyvalsshape <- c(keyvalsshape, sha)
            }
            #Assign LegendParameter for EnahncedVolcano:
            LegendPos<- "top"
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
                                                  xlim =  c(min(InputVolcano$Log2FC[is.finite(InputVolcano$Log2FC )])-0.2, max(InputVolcano$Log2FC[is.finite(InputVolcano$Log2FC )])+1.2),
                                                  ylim = c(0,(ceiling(-log10(Reduce(min,InputVolcano$p.adj))))),
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

          #Add PieChart
          if("color" %in% names(Plot_SettingsInfo)==TRUE & Legend=="Pie"){
            Plot <- Plot+ annotation_custom(
              grob = ggplotGrob(PieChart),
              xmin = (max(InputVolcano$Log2FC[is.finite(InputVolcano$Log2FC)])+0.2), xmax =(max(InputVolcano$Log2FC[is.finite(InputVolcano$Log2FC)])+1.5),
              ymin = ((ceiling(-log10(Reduce(min,InputVolcano$p.adj))))+0.5), ymax =((ceiling(-log10(Reduce(min,InputVolcano$p.adj))))-1.5)
            )
          }

          #save plot and get rid of extra signs before saving
          if(OutputPlotName ==""){
            ggsave(file=paste(Results_folder_plots_Volcano_folder,"/", "Volcano." ,Save_as_Plot, sep=""), plot=Plot, width=8, height=6)
          }else{
            ggsave(file=paste(Results_folder_plots_Volcano_folder,"/", "Volcano_", OutputPlotName, ".",Save_as_Plot, sep=""), plot=Plot, width=8, height=6)
          }
          #Plot
          plot(Plot)
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
                                                  title= paste(OutputPlotName, ": ", i, sep=""),
                                                  subtitle = Subtitle,
                                                  caption = paste0("total = ", (nrow(InputVolcano)/2), " Metabolites"),
                                                  xlim =  c(min(InputVolcano$Log2FC[is.finite(InputVolcano$Log2FC )])-0.2, max(InputVolcano$Log2FC[is.finite(InputVolcano$Log2FC )])+1.2),
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
          #save plot and get rid of extra signs before saving
          cleaned_i <- gsub("[[:space:],/\\\\]", "-", i)#removes empty spaces and replaces /,\ with -
          if(OutputPlotName ==""){
            ggsave(file=paste(Results_folder_plots_Volcano_folder,"/", "Volcano_",cleaned_i, ".",Save_as_Plot, sep=""), plot=Plot, width=8, height=6)
          }else{
            ggsave(file=paste(Results_folder_plots_Volcano_folder,"/", "Volcano_", OutputPlotName, "_",cleaned_i, ".",Save_as_Plot, sep=""), plot=Plot, width=8, height=6)
          }
          ## Store the plot in the 'plots' list
          PlotList[[cleaned_i]] <- Plot
          plot(Plot)
        }
      }
      # Return PlotList into the environment to enable the user to view the plots directly
      #assign("VolcanoPlots", PlotList, envir=.GlobalEnv)
      # Combine plots into a single plot using facet_grid or patchwork::wrap_plots
      Return <- PlotList
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
                                                caption = paste0("total = ", (nrow(InputVolcano)/2), " Metabolites"),
                                                xlim =  c(min(InputVolcano$Log2FC[is.finite(InputVolcano$Log2FC )])-0.2, max(InputVolcano$Log2FC[is.finite(InputVolcano$Log2FC )])+1.2),
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
          ggsave(file=paste(Results_folder_plots_Volcano_folder,"/", "Volcano." ,Save_as_Plot, sep=""), plot=Plot, width=8, height=6)
        }else{
          ggsave(file=paste(Results_folder_plots_Volcano_folder,"/", "Volcano_", OutputPlotName, ".",Save_as_Plot, sep=""), plot=Plot, width=8, height=6)
        }
        #Plot
        plot(Plot)
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
                                                title= paste(OutputPlotName, ": ", i, sep=""),
                                                subtitle = paste(Plot_SettingsInfo[["PEA_score"]],"= ", AdditionalInput_data_Select$PEA_score, ", ",Plot_SettingsInfo[["PEA_stat"]] , "= ", AdditionalInput_data_Select$PEA_stat, sep=""),
                                                caption = paste0("total = ", nrow(InputVolcano), " metabolites of ", nrow(Plot_SettingsFile_Select), " metabolites in pathway"),
                                                xlim =  c(min(InputVolcano$Log2FC[is.finite(InputVolcano$Log2FC )])-0.2, max(InputVolcano$Log2FC[is.finite(InputVolcano$Log2FC )])+1.2),
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
        #save plot and get rid of extra signs before saving
        cleaned_i <- gsub("[[:space:],/\\\\]", "-", i)#removes empty spaces and replaces /,\ with -
        if(OutputPlotName ==""){
          ggsave(file=paste(Results_folder_plots_Volcano_folder,"/", "Volcano_",cleaned_i, ".",Save_as_Plot, sep=""), plot=Plot, width=8, height=6)
        }else{
          ggsave(file=paste(Results_folder_plots_Volcano_folder,"/", "Volcano_", OutputPlotName, "_",cleaned_i, ".",Save_as_Plot, sep=""), plot=Plot, width=8, height=6)
        }
        ## Store the plot in the 'plots' list
        PlotList[[cleaned_i]] <- Plot
        plot(Plot)
      }
    }
    # Return PlotList into the environment to enable the user to view the plots directly
    #assign("VolcanoPlots", PlotList, envir=.GlobalEnv)
    # Combine plots into a single plot using facet_grid or patchwork::wrap_plots
    Return <- PlotList

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
                     Save_as_Plot = pdf
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
  ### ### ### make output plot Save_as_Plot name ### ### ###
  Save_as_Plot_var <- Save_as_Plot
  Save_as_Plot= deparse(substitute(Save_as_Plot))


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


  Save_as_Plot_var(paste(Results_folder_plots_MetabolicCluster_folder,"/Metabolic_Clusters_",Condition1,"-versus-",Condition2,"_", Output_Name,  ".",Save_as_Plot, sep=""), width=12, height=9)
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
#' @param Plot_Settings \emph{Optional: } Choose between "Standard" (Input_data), "Compare" (plot two comparisons together Input_data and Input_data2) or "PEA" (Pathway Enrichment Analysis) \strong{Default = "Standard"}
#' @param Plot_SettingsInfo \emph{Optional: } NULL or Named vector including at least one of those three information for Plot_Settings="Standard" or "Compare": c(color="ColumnName_Plot_SettingsFile", shape= "ColumnName_Plot_SettingsFile", individual="ColumnName_Plot_SettingsFile"). For Plot_Settings="PEA" a named vector with c(PEA_Pathway="ColumnNameAdditionalInput_data", PEA_score="ColumnNameAdditionalInput_data", PEA_stat= "ColumnNameAdditionalInput_data", individual="Plot_SettingsFile), optionally you can additionally include c(color="ColumnName_Plot_SettingsFile", shape= "ColumnName_Plot_SettingsFile").\strong{Default = NULL}
#' @param Plot_SettingsFile \emph{Optional: } DF with column "Metabolite" including the Metabolite names (needs to match Metabolite names of Input_data) and other columns with required PlotSettingInfo. \strong{Default = NULL}
#' @param Input_data DF with column "Metabolite" including the Metabolite names, Log2FC, pvalue/padjusted values. Can also include additional columns with metadata usable for Plot_Setting_Info. Multiple Input data frames can be added using a list ie. list(df1, df2, df3)
#' @param OutputPlotName \emph{Optional: } String which is added to the output files of the plot. \strong{Default = ""}
#' @param test \emph{Optional: } String which selects pvalue or padj for significance. \strong{Default = padj}
#' @param Comparison_name \emph{Optional: } List including those information about the datasets that are compared on the plots when choosing Plot_Settings= "Compare". \strong{Default = list("Cond1", "Cond2")}
#' @param pCutoff \emph{Optional: } Number of the desired p value cutoff for assessing significance. \strong{Default = 0.05}
#' @param FCcutoff \emph{Optional: } Number of the desired log fold change cutoff for assessing significance. \strong{Default = 0.5}
#' @param Legend \emph{Optional: } Legend=="Pie" will plot a PieChart as the legend for color or Legend="Standard, plot the standard legend for color. \strong{Default = "Standard"}
#' @param Subtitle \emph{Optional: } \strong{Default = ""}
#' @param Theme \emph{Optional: } Selection of theme for plot. \strong{Default = NULL}
#'
#' @param Save_as_Plot \emph{Optional: } Select the file type of output plots. Options are svg, pdf or png. \strong{Default = "svg"}
#'
#' @keywords Lolipop plot, pathways
#' @export

# Helper function needed for adding column to pathway file defining if this metabolite is unique/multiple pathways



VizLolipop<- function(Plot_Settings="Standard",
                      Plot_SettingsInfo= NULL,
                      Plot_SettingsFile= NULL, # Input_pathways = NULL,
                      Input_data, # a dataframe of list of dataframes
                      AdditionalInput_data= NULL, # unsd only for PEA
                      x = "Log2FC",
                      y = "Metabolite",
                      OutputPlotName= "",
                      Comparison_name= c(Input_data="Cond1", AdditionalInput_data= "Cond2"),
                      Subtitle= "",
                      Theme= NULL,
                      Save_as_Plot = "svg",
                      parameter_size="Reverse" #or default "Standard"
){


  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse", "showtext", "cowplot")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  suppressMessages(library(tidyverse))

  ## ------------ Check Input files ----------- ##
  # 1. The input data:
  if(inherits(Input_data,"list") ==FALSE){
    temp <- list(Input_data)
    Input_data <- temp
  }
  if(inherits(Plot_SettingsFile,"list") ==FALSE){
    temp <- list(Plot_SettingsFile)
    Plot_SettingsFile <- temp
  }
  Plot_SettingsFile_List <- Plot_SettingsFile
  if(length(Plot_SettingsFile_List)==1){
    Plot_SettingsFile <- Plot_SettingsFile_List[[1]]
  }

  Flip <- c() # Vector to check if x or y in numeric
  for(i in 1:length(Input_data)){
    # i=2
    data <- Input_data[[i]]
    if(any(duplicated(row.names(data)))==TRUE){
      stop("Duplicated row.names of Input_data, whilst row.names must be unique")
    }
    if("Metabolite" %in% colnames(data) == FALSE){
      stop("Check input. Input_data must contain a column named `Metabolite` including the metabolite names.")
    }
    # Check if the next lines work correctly in case of duplicated metabolties
    if(length(data[duplicated(data$Metabolite), "Metabolite"]) > 0){
      doublons <- as.character(data[duplicated(data$Metabolite), "Metabolite"])#number of duplications
      data <-data[!duplicated(data$Metabolite),]#remove duplications
      warning("Input_data contained duplicates based on Metabolite! Dropping duplicate IDs and kept only the first entry. You had ", length(doublons), " duplicates. Note that you should do this before running VizLolipop.")
    }
    if(paste(x) %in% colnames(data)==FALSE | paste(y) %in% colnames(data)==FALSE){
      stop("Check your input. The column name of x and/ore y does not exist in Input_data.")
    }
    if (is.numeric(data[[x]]) && is.character(data[[y]])) {
      Flip[i] <-FALSE
    } else if (is.character(data[[x]]) && is.numeric(data[[y]])) {
      Flip[i] <- TRUE

    } else {
      stop("One of the x or y must by numeric and the other must be a character")
    }

  }

  if (sum(Flip) == length(Input_data)) {
    Flip=TRUE
    temp<- x
    x<-y
    y<- temp
  } else{
    Flip <- FALSE
  }


  # 2. The Plot_settings: Plot_Settings, Plot_SettingInfo and Plot_SettingFile
  Plot_options <- c("Standard", "Compare", "PEA")
  if (Plot_Settings %in% Plot_options == FALSE){
    stop("Plot_Settings option is incorrect. The allowed options are the following: ",paste(Plot_options, collapse = ", "),"." )
  }
  if(is.vector(Plot_SettingsInfo)==TRUE & is.null(Plot_SettingsFile)==TRUE){
    stop("You have chosen Plot_SettingsInfo option that requires you to provide a DF Plot_SettingsFile.")
  }
  if(Plot_Settings=="Compare" & "color" %in% names(Plot_SettingsInfo)==TRUE){
    warning("When Plot_Settings='Compare' color is used to colorcode the input datasets. Color was added in the Plot_SettingsInfo but it will be ignored. Please use only the size option.")
    Plot_SettingsInfo <- Plot_SettingsInfo[!names(Plot_SettingsInfo) == "color"]
  }

  if("color" %in% names(Plot_SettingsInfo)==TRUE){
    for(i in 1:length(Plot_SettingsFile_List)){
      data <- Plot_SettingsFile_List[[i]]
      if (Plot_SettingsInfo[["color"]] %in% colnames(data) == FALSE) {
        stop("You have chosen color = ",paste(Plot_SettingsInfo[["color"]]), ", ", paste(Plot_SettingsInfo[["color"]])," does not exist in the PlotSettingsFile."   )
      }
    }
  }
  if("size" %in% names(Plot_SettingsInfo)==TRUE){
    for(i in 1:length(Plot_SettingsFile_List)){
      data <- Plot_SettingsFile_List[[i]]
      if (Plot_SettingsInfo[["size"]] %in% colnames(data) == FALSE) {
        stop("You have chosen size = ",paste(Plot_SettingsInfo[["size"]]), ", ", paste(Plot_SettingsInfo[["size"]])," does not exist in the PlotSettingsFile."   )
      }
    }
  }
  if("label_dot" %in% names(Plot_SettingsInfo)==TRUE){
    for(i in 1:length(Plot_SettingsFile_List)){
      data <- Plot_SettingsFile_List[[i]]
      if (Plot_SettingsInfo[["label_dot"]] %in% colnames(data) == FALSE) {
        stop("You have chosen label_dot = ",paste(Plot_SettingsInfo[["label_dot"]]), ", ", paste(Plot_SettingsInfo[["label_dot"]])," does not exist in the PlotSettingsFile."   )
      }
    }
  }

  if(is.vector(Plot_SettingsInfo)==TRUE){
    if("color" %in% names(Plot_SettingsInfo)==TRUE & "size" %in% names(Plot_SettingsInfo)==TRUE){
      if((Plot_SettingsInfo[["size"]] == Plot_SettingsInfo[["color"]])==TRUE){
        Plot_SettingsFile$size <- Plot_SettingsFile[,paste(Plot_SettingsInfo[["color"]])]
        Plot_SettingsFile<- Plot_SettingsFile%>%
          dplyr::rename("color"=paste(Plot_SettingsInfo[["color"]]))
      }
      if((Plot_SettingsInfo[["size"]] == Plot_SettingsInfo[["color"]])==FALSE & "color" %in% names(Plot_SettingsInfo)==TRUE){
        Plot_SettingsFile <- Plot_SettingsFile%>%
          dplyr::rename("color"=paste(Plot_SettingsInfo[["color"]]))
      }
      if((Plot_SettingsInfo[["size"]] == Plot_SettingsInfo[["color"]])==FALSE & "size" %in% names(Plot_SettingsInfo)==TRUE){
        Plot_SettingsFile <- Plot_SettingsFile%>%
          dplyr::rename("size"=paste(Plot_SettingsInfo[["size"]]))
      }
      Plot_SettingsFile2 <- Plot_SettingsFile%>% select(Metabolite, color, size)########################################################
    } else if("color" %in% names(Plot_SettingsInfo)==TRUE & "size" %in% names(Plot_SettingsInfo)==FALSE){
      Plot_SettingsFile <- Plot_SettingsFile%>%
        dplyr::rename("color"=paste(Plot_SettingsInfo[["color"]]))
      Plot_SettingsFile2 <- Plot_SettingsFile%>% select(Metabolite, color)##########################################################
    } else if("color" %in% names(Plot_SettingsInfo)==FALSE & "size" %in% names(Plot_SettingsInfo)==TRUE){
      if(length(Plot_SettingsFile_List)==1){
        Plot_SettingsFile <- Plot_SettingsFile%>%
          dplyr::rename("size"=paste(Plot_SettingsInfo[["size"]]))
        Plot_SettingsFile2 <- Plot_SettingsFile%>% select(Metabolite,size) ##############################
      }else if(length(Plot_SettingsFile_List)>1){
        for(i in 1:length(Plot_SettingsFile)){
          file <- Plot_SettingsFile[[i]]
          file <- file%>%
            dplyr::rename("size"=paste(Plot_SettingsInfo[["size"]]))
          Plot_SettingsFile_List[[i]] <- file %>% select(Metabolite, size)
          Plot_SettingsFile2 <- Plot_SettingsFile_List[[i]] # Needs fix This is not needed it is here for the code not to break
        }
      }
    }


    if("individual" %in% names(Plot_SettingsInfo)==TRUE){
      if(length(Plot_SettingsFile_List)==1){
        Plot_SettingsFile <- Plot_SettingsFile%>%
          dplyr::rename("individual"=paste(Plot_SettingsInfo[["individual"]]))
        if("color" %in% names(Plot_SettingsInfo) | "size" %in% names(Plot_SettingsInfo)){
          Plot_SettingsFile3 <- Plot_SettingsFile %>% select(individual, Metabolite)
          Plot_SettingsFile <- merge(Plot_SettingsFile2, Plot_SettingsFile3, by="Metabolite")
        }
      }else if(length(Plot_SettingsFile_List)>1){
        for(i in 1:length(Plot_SettingsFile)){
          file <- Plot_SettingsFile[[i]]
          file <- file%>%
            dplyr::rename("individual"=paste(Plot_SettingsInfo[["individual"]]))

          if("size" %in% names(Plot_SettingsInfo)){
            Plot_SettingsFile_List[[i]] <- merge(Plot_SettingsFile_List[[i]] ,file %>% select(Metabolite, individual), by = "Metabolite")
            Plot_SettingsFile2 <- Plot_SettingsFile_List[[i]] # Needs fix This is not needed it is here for the code not to break

          }else{
            Plot_SettingsFile_List[[i]] <- file %>% select(Metabolite, individual)
            Plot_SettingsFile2 <- Plot_SettingsFile_List[[i]] # Needs fix This is not needed it is here for the code not to break
          }
        }
      }
    }else{
      Plot_SettingsFile <- Plot_SettingsFile2
    }
  }else if(is.vector(Plot_SettingsInfo)==FALSE & is.null(Plot_SettingsInfo)==FALSE){
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

  ## The next lines need checks/corrections
  # if(Plot_Settings=="PEA" & length(Input_data)==1){
  #   stop("If Plot_Settings=`PEA` you have to provide a DF for AdditionalInput_data including the results of an enrichment analysis.")
  # } else if(Plot_Settings=="PEA" & length(Input_data)>1){
  #   AdditionalInput_data <- AdditionalInput_data%>%
  #     dplyr::rename("PEA_score"=paste(Plot_SettingsInfo[["PEA_score"]]),
  #                   "PEA_stat"=paste(Plot_SettingsInfo[["PEA_stat"]]),
  #                   "PEA_Pathway"=paste(Plot_SettingsInfo[["PEA_Pathway"]]))
  # }

  # The next lines need checks/corrections
  if(Plot_Settings=="PEA" & is.null(AdditionalInput_data)==TRUE){
    stop("If Plot_Settings=`PEA` you have to provide a DF for AdditionalInput_data including the results of an enrichment analysis.")
  } else if(Plot_Settings=="PEA" & is.null(AdditionalInput_data)==FALSE){
    AdditionalInput_data <- AdditionalInput_data%>%
      dplyr::rename("PEA_score"=paste(Plot_SettingsInfo[["PEA_score"]]),
                    "PEA_stat"=paste(Plot_SettingsInfo[["PEA_stat"]]),
                    "PEA_Pathway"=paste(Plot_SettingsInfo[["PEA_Pathway"]]))
  }


  # 4. Check other plot-specific parameters:
  Save_as_Plot_options <- c("svg","pdf","png")
  if(Save_as_Plot %in% Save_as_Plot_options == FALSE){
    stop("Check input. The selected Save_as_Plot option is not valid. Please select one of the following: ",paste(Save_as_Plot_options,collapse = ", "),"." )
  }



  ## ------------ Create Output folders ----------- ##
  name <- paste0("MetaProViz_Results_",Sys.Date())
  WorkD <- getwd()
  Results_folder <- file.path(WorkD, name)
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)} # Make Results folder
  Results_folder_plots_Lolipop_folder = file.path(Results_folder, "Lolipop")  # This searches for a folder called "Preprocessing" within the "Results" folder in the current working directory and if its not found it creates one
  if (!dir.exists(Results_folder_plots_Lolipop_folder)) {dir.create(Results_folder_plots_Lolipop_folder)}  # check and create folder


  if(Plot_Settings=="Standard"){############################################################################################################
    Input_data <- Input_data[[1]]
    if("individual" %in% names(Plot_SettingsInfo)==TRUE){
      # Create the list of individual plots that should be made:
      IndividualPlots <- Plot_SettingsFile[!duplicated(Plot_SettingsFile$individual),]
      IndividualPlots <- IndividualPlots$individual

      PlotList <- list()#Empty list to store all the plots

      for (i in IndividualPlots){
        # i = IndividualPlots[1]
        Plot_SettingsInfo_indi <- Plot_SettingsInfo
        Plot_SettingsFile_Select <- subset(Plot_SettingsFile, individual == paste(i))
        InputLolipop  <- merge(x=Plot_SettingsFile_Select,y=Input_data, by="Metabolite", all.x=TRUE)%>%
          na.omit()

        #Select metabolites for the cut offs selected
        loli.data <- InputLolipop %>% mutate(names=Metabolite)

        if("size" %in% names(Plot_SettingsInfo_indi)==TRUE ){
          if(is.numeric(loli.data$size)==FALSE){ # run is color is discrete
            stop("Size can take only numeric values")
          }else{# color = continuous
            keyvalssize <- loli.data$size
          }
        } else{
          Plot_SettingsInfo_indi= c(Plot_SettingsInfo_indi,size="p.adj")
          keyvalssize <- loli.data$size
        }


        if("label_dot" %in% names(Plot_SettingsInfo)==TRUE ){
          label <-  Plot_SettingsInfo[["label_dot"]]
          loli.data[ Plot_SettingsInfo[["label_dot"]]] <- round(loli.data[[label]], digits = 3)
        }else{
          label <-  ""
        }

        p2 <- NULL
        if("color" %in% names(Plot_SettingsInfo_indi)==TRUE ){
          if(is.numeric(loli.data$color)==FALSE){ # run if color is discrete

            col_var_name <- Plot_SettingsInfo_indi[['color']]

            position <- which(names(loli.data)=="color" )
            names(loli.data)[position]<-"plot_color_variable"


            loli.data <- loli.data %>%
              arrange(plot_color_variable, get(x),Metabolite)

            loli.data_avg <- loli.data %>%
              arrange(plot_color_variable, get(x), Metabolite) %>%
              mutate(Metab_name = row_number()) %>%
              group_by(plot_color_variable) %>%
              mutate(
                avg = mean(get(x))
              ) %>%
              ungroup() %>%
              mutate(plot_color_variable = factor(plot_color_variable))


            loli_lines <-   loli.data_avg %>%
              arrange(plot_color_variable, Metabolite) %>%
              group_by(plot_color_variable) %>%
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
              ggplot(aes(Metab_name, get(x))) + # names in aes ro Metab_name
              geom_hline(
                data = tibble(y = -5:5),
                aes(yintercept = y),
                color = "grey82",
                size = .5 )

            p2 <- p2 + geom_segment(
              aes(
                xend = Metab_name,          # names
                yend = 0,#avg,
                color = plot_color_variable,
                #color = after_scale(colorspace::lighten(color, .2))
              ))

            p2 <- p2 + # geom_line( data = loli_lines, aes(x, y),  color = "grey40"  ) +
              geom_line(
                data = loli_lines,
                aes( x_group, y,
                     color = plot_color_variable,
                     #  color = after_scale(colorspace::darken(color, .2))
                ), size = 2.5) +  geom_point(aes(size = keyvalssize, color = plot_color_variable)
                )

            p2 <- p2 +theme(axis.text.y=element_blank(),
                            axis.ticks.y=element_blank()
            )

            if(Flip == TRUE){
              p2 <- p2 +theme(axis.text.x=element_blank(),
                              axis.ticks.x=element_blank())
              lab_pos_metab <- loli.data_avg[loli.data_avg[x]>0,]  %>% select(Metabolite, Metab_name, x)

              p2 <- p2+ annotate("text", x = lab_pos_metab$Metab_name, y = lab_pos_metab[[x]]+1.5, label = lab_pos_metab$Metabolite,angle = 90, size = 3)

              lab_neg_metab <-  loli.data_avg[loli.data_avg[x]<0,]  %>% select(Metabolite, Metab_name, x)
              p2<- p2+ annotate("text", x = lab_neg_metab$Metab_name, y = lab_neg_metab[[x]]-1.5, label = lab_neg_metab$Metabolite, angle = 90,size = 3)

              if("label_dot" %in% names(Plot_SettingsInfo)==TRUE ){
                dot_pos_metab <- loli.data_avg[loli.data_avg[x]>0,]  %>% select(Metabolite, Metab_name, x, label)
                p2<- p2+ annotate("text", x = dot_pos_metab$Metab_name, y = dot_pos_metab[[x]], label = dot_pos_metab[[label]], size = 3)

                dot_neg_metab <- loli.data_avg[loli.data_avg[x]<0,]  %>% select(Metabolite, Metab_name, x,label)
                p2<- p2+ annotate("text", x = dot_neg_metab$Metab_name, y = dot_neg_metab[[x]], label = dot_neg_metab[[label]], size = 3)

              }

              p2 <- p2+Theme
              p2 <- p2+ labs(color=col_var_name)+
                labs(size=Plot_SettingsInfo[['size']])

              p2 <- p2  + labs(title = paste(OutputPlotName),subtitle = Subtitle)+
                theme(plot.title = element_text(color = "black", size = 12, face = "bold"),
                      plot.subtitle = element_text(color = "black", size=10),
                      plot.caption = element_text(color = "black",size=9, face = "italic", hjust = 2.5))

            }else{
              p2<- p2 + coord_flip()
              p2 <- p2 +theme(axis.text.x=element_blank(),
                              axis.ticks.x=element_blank())
              lab_pos_metab <- loli.data_avg %>% filter(get(x)>0) %>% select(Metabolite, Metab_name, get(x))
              p2<- p2+ annotate("text", x = lab_pos_metab$Metab_name, y = lab_pos_metab[[x]]+1.5, label = lab_pos_metab$Metabolite, size = 3)

              lab_neg_metab <- loli.data_avg %>% filter(get(x)<0) %>% select(Metabolite, Metab_name, get(x))
              p2<- p2+ annotate("text", x = lab_neg_metab$Metab_name, y = lab_neg_metab[[x]]-1.5, label = lab_neg_metab$Metabolite, size = 3)

              if("label_dot" %in% names(Plot_SettingsInfo)==TRUE ){
                dot_pos_metab <- loli.data_avg[loli.data_avg[x]>0,]  %>% select(Metabolite, Metab_name, x, label)
                p2<- p2+ annotate("text", x = dot_pos_metab$Metab_name, y = dot_pos_metab[[x]], label = dot_pos_metab[[label]], size = 3)

                dot_neg_metab <- loli.data_avg[loli.data_avg[x]<0,]  %>% select(Metabolite, Metab_name, x,label)
                p2<- p2+ annotate("text", x = dot_neg_metab$Metab_name, y = dot_neg_metab[[x]], label = dot_neg_metab[[label]], size = 3)

              }

              p2 <- p2+Theme

              p2 <- p2  + labs(title = paste(OutputPlotName,": ", i ),subtitle = Subtitle) +
                theme(plot.title = element_text(color = "black", size = 12, face = "bold"),
                      plot.subtitle = element_text(color = "black", size=10),
                      plot.caption = element_text(color = "black",size=9, face = "italic", hjust = 2.5))

              p2 <- p2+ labs(color=col_var_name)+
                labs(size=Plot_SettingsInfo_indi[['size']])
            }

            lolipop_plot <-  p2
            if(parameter_size=="Reverse" & is.null(keyvalssize)==FALSE){
              lolipop_plot <-   lolipop_plot + scale_size(trans = 'reverse')
            }

            # Put back the correct name in the data df
            names(loli.data)[position]<- col_var_name

          }else{# color = continuous
            keyvals <- loli.data$color
          }
        }else{
          Plot_SettingsInfo_indi= c(Plot_SettingsInfo_indi, color="p.adj")
          keyvals <- loli.data$color
        }

        if(is.null(p2)==TRUE){

          loli.data$names <- as.factor(loli.data$names)
          loli.data[[x]]<- as.numeric(loli.data[[x]])
          loli.data$names <- reorder(loli.data$names, -loli.data[[x]])

          lolipop_plot <- ggplot(loli.data , aes(x = get(x), y = names,label=!!label)) +
            geom_segment(aes(x = 0, xend = get(x), y = names, yend = names)) +
            geom_point(aes(colour = keyvals, size = keyvalssize ))   +
            # scale_size_continuous(range = c(1,5))+# , trans = 'reverse') +

            scale_colour_gradient(low = "red", high = "blue")+
            geom_vline(xintercept = 0)+
            Theme+
            labs(color=Plot_SettingsInfo_indi[['color']])+
            labs(size=Plot_SettingsInfo_indi[['size']]) +
            ylab(y)+
            xlab(x)+
            labs(title = paste(OutputPlotName,": ", i ),subtitle = Subtitle)+
            theme(plot.title = element_text(color = "black", size = 12, face = "bold"),
                  plot.subtitle = element_text(color = "black", size=10),
                  plot.caption = element_text(color = "black",size=9, face = "italic", hjust = 2.5))

          if(parameter_size=="Reverse" & is.null(keyvalssize)==FALSE){
            lolipop_plot <-   lolipop_plot + scale_size(trans = 'reverse')
          }
          if(Flip==TRUE){
            lolipop_plot <- lolipop_plot + coord_flip() +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
          }
          if("label_dot" %in% names(Plot_SettingsInfo)==TRUE ){
            lolipop_plot <-   lolipop_plot + geom_text(color="black", size=2)
          }
        }

        #save plot and get rid of extra signs before saving
        cleaned_i <- gsub("[[:space:],/\\\\]", "-", i)#removes empty spaces and replaces /,\ with -
        if(OutputPlotName ==""){
          ggsave(file=paste(Results_folder_plots_Lolipop_folder,"/", "Lolipop_",cleaned_i, ".",Save_as_Plot, sep=""), plot=lolipop_plot, width=8, height=6)
        }else{
          ggsave(file=paste(Results_folder_plots_Lolipop_folder,"/", "Lolipop_", OutputPlotName, "_",cleaned_i, ".",Save_as_Plot, sep=""), plot=lolipop_plot, width=8, height=6)
        }
        ## Store the plot in the 'plots' list
        PlotList[[cleaned_i]] <- lolipop_plot
        plot(lolipop_plot)
      }
      # Return PlotList into the environment to enable the user to view the plots directly
      #assign("LolipopPlots", PlotList, envir=.GlobalEnv)
      # Combine plots into a single plot using facet_grid or patchwork::wrap_plots
      Return <- PlotList
    }
    else if("individual" %in% names(Plot_SettingsInfo)==FALSE){

      if(is.null(Plot_SettingsFile)==FALSE){
        InputLolipop  <- merge(x=Plot_SettingsFile,y=Input_data, by="Metabolite", all.x=TRUE)%>%
          na.omit()
      }else{
        InputLolipop  <- Input_data
      }

      InputLolipop<- InputLolipop %>% drop_na()

      #Select metabolites for the cut offs selected
      loli.data <- InputLolipop %>% mutate(names=Metabolite)

      if("size" %in% names(Plot_SettingsInfo)==TRUE ){
        if(is.numeric(loli.data$size)==FALSE){ # run is color is discrete
          stop("Size can take only numeric values")
        }else{# color = continuous
          keyvalssize <- loli.data$size
        }
      } else{
        Plot_SettingsInfo= c(Plot_SettingsInfo,size="p.adj")
        keyvalssize <- loli.data$size
      }


      if("label_dot" %in% names(Plot_SettingsInfo)==TRUE ){
        label <-  Plot_SettingsInfo[["label_dot"]]
        loli.data[ Plot_SettingsInfo[["label_dot"]]] <- round(loli.data[[label]], digits = 3)
      }else{
        label <-  ""
      }


      p2 <- NULL
      if("color" %in% names(Plot_SettingsInfo)==TRUE ){
        if(is.numeric(loli.data$color)==FALSE){ # run if color is discrete

          col_var_name <- Plot_SettingsInfo[['color']]

          position <- which(names(loli.data)=="color" )
          names(loli.data)[position]<-"plot_color_variable"


          loli.data <- loli.data %>%
            arrange(plot_color_variable,get(x), Metabolite)

          loli.data_avg <- loli.data %>%
            arrange(plot_color_variable,get(x), Metabolite) %>%
            mutate(Metab_name = row_number()) %>%
            group_by(plot_color_variable) %>%
            mutate(
              avg = mean(get(x))
            ) %>%
            ungroup() %>%
            mutate(plot_color_variable = factor(plot_color_variable))



          loli_lines <-   loli.data_avg %>%
            arrange(plot_color_variable, Metabolite) %>%
            group_by(plot_color_variable) %>%
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
            ggplot(aes(Metab_name, get(x))) + # names in aes ro Metab_name
            geom_hline(
              data = tibble(y = -5:5),
              aes(yintercept = y),
              color = "grey82",
              size = .5 )

          p2 <- p2 + geom_segment(
            aes(
              xend = Metab_name,          # names
              yend = 0,#avg,
              color = plot_color_variable,
              #color = after_scale(colorspace::lighten(color, .2))
            ))

          p2 <- p2 + # geom_line( data = loli_lines, aes(x, y),  color = "grey40"  ) +
            geom_line(
              data = loli_lines,
              aes( x_group, y,
                   color = plot_color_variable,
                   #  color = after_scale(colorspace::darken(color, .2))
              ), size = 2.5) +  geom_point(aes(size = keyvalssize, color = plot_color_variable)
              )

          if(Flip == TRUE){
            p2 <- p2 +theme(axis.text.x=element_blank(),
                            axis.ticks.x=element_blank())
            lab_pos_metab <- loli.data_avg[loli.data_avg[x]>0,]  %>% select(Metabolite, Metab_name, x)
            p2 <- p2+ annotate("text", x = lab_pos_metab$Metab_name, y = lab_pos_metab[[x]]+1.5, label = lab_pos_metab$Metabolite,angle = 90, size = 3)

            lab_neg_metab <- loli.data_avg[loli.data_avg[x]<0,]  %>% select(Metabolite, Metab_name, x)
            p2<- p2+ annotate("text", x = lab_neg_metab$Metab_name, y = lab_neg_metab[[x]]-1.5, label = lab_neg_metab$Metabolite, angle = 90,size = 3)

            if("label_dot" %in% names(Plot_SettingsInfo)==TRUE ){
              dot_pos_metab <- loli.data_avg[loli.data_avg[x]>0,]  %>% select(Metabolite, Metab_name, x, label)
              p2<- p2+ annotate("text", x = dot_pos_metab$Metab_name, y = dot_pos_metab[[x]], label = dot_pos_metab[[label]], size = 3)

              dot_neg_metab <- loli.data_avg[loli.data_avg[x]<0,]  %>% select(Metabolite, Metab_name, x,label)
              p2<- p2+ annotate("text", x = dot_neg_metab$Metab_name, y = dot_neg_metab[[x]], label = dot_neg_metab[[label]], size = 3)

            }

            p2 <- p2+Theme
            p2 <- p2+ labs(color=col_var_name)+
              labs(size=Plot_SettingsInfo[['size']])

            p2 <- p2  + labs(title = paste(OutputPlotName),subtitle = Subtitle)+
              theme(plot.title = element_text(color = "black", size = 12, face = "bold"),
                    plot.subtitle = element_text(color = "black", size=10),
                    plot.caption = element_text(color = "black",size=9, face = "italic", hjust = 2.5))

          }else{
            p2<- p2 + coord_flip()
            p2 <- p2 +theme(axis.text.y=element_blank(),
                            axis.ticks.y=element_blank()
            )
            lab_pos_metab <- loli.data_avg[loli.data_avg[x]>0,]  %>% select(Metabolite, Metab_name, x)
            p2<- p2+ annotate("text", x = lab_pos_metab$Metab_name, y = lab_pos_metab[[x]]+1.5, label = lab_pos_metab$Metabolite, size = 3)


            lab_neg_metab <- loli.data_avg[loli.data_avg[x]<0,]  %>% select(Metabolite, Metab_name, x)
            p2<- p2+ annotate("text", x = lab_neg_metab$Metab_name, y = lab_neg_metab[[x]]-1.5, label = lab_neg_metab$Metabolite, size = 3)

            # p2 <- p2+ annotate("text", x = max(lab_neg_metab$Metab_name)+ 7, y = 0, label = OutputPlotName, size = 5)

            if("label_dot" %in% names(Plot_SettingsInfo)==TRUE ){
              dot_pos_metab <- loli.data_avg[loli.data_avg[x]>0,]  %>% select(Metabolite, Metab_name, x, label)
              p2<- p2+ annotate("text", x = dot_pos_metab$Metab_name, y = dot_pos_metab[[x]], label = dot_pos_metab[[label]], size = 3)

              dot_neg_metab <- loli.data_avg[loli.data_avg[x]<0,]  %>% select(Metabolite, Metab_name, x,label)
              p2<- p2+ annotate("text", x = dot_neg_metab$Metab_name, y = dot_neg_metab[[x]], label = dot_neg_metab[[label]], size = 3)

            }


            p2 <- p2+Theme
            p2 <- p2+ labs(color=col_var_name)+
              labs(size=Plot_SettingsInfo[['size']])

            p2 <- p2  + labs(title = paste(OutputPlotName),subtitle = Subtitle)+
              theme(plot.title = element_text(color = "black", size = 12, face = "bold"),
                    plot.subtitle = element_text(color = "black", size=10),
                    plot.caption = element_text(color = "black",size=9, face = "italic", hjust = 2.5))

          }
          lolipop_plot <-  p2
          if(parameter_size=="Reverse" & is.null(keyvalssize)==FALSE){
            lolipop_plot <-   lolipop_plot + scale_size(trans = 'reverse')
          }

          # Put back the correct name in the data df
          names(loli.data)[position]<- col_var_name

        }else{# color = continuous
          keyvals <- loli.data$color
        }
      }else{
        Plot_SettingsInfo= c(Plot_SettingsInfo, color="p.adj")
        keyvals <- loli.data$color
      }

      if(is.null(p2)==TRUE){

        loli.data$names <- as.factor(loli.data$names)
        loli.data[[x]]<- as.numeric(loli.data[[x]])
        loli.data$names <- reorder(loli.data$names, -loli.data[[x]])

        lolipop_plot <- ggplot(loli.data , aes(x = get(x), y = names,label=!!label)) +
          geom_segment(aes(x = 0, xend = get(x), y = names, yend = names)) +
          geom_point(aes(colour = keyvals, size = keyvalssize ))   +
          # scale_size_continuous(range = c(1,5))+# , trans = 'reverse') +
          scale_colour_gradient(low = "red", high = "blue")+
          geom_vline(xintercept = 0)+
          Theme+
          labs(color=Plot_SettingsInfo[['color']])+
          labs(size=Plot_SettingsInfo[['size']]) +
          ylab(y)+
          xlab(x)+
          labs(title = OutputPlotName,subtitle = Subtitle)+
          theme(plot.title = element_text(color = "black", size = 12, face = "bold"),
                plot.subtitle = element_text(color = "black", size=10),
                plot.caption = element_text(color = "black",size=9, face = "italic", hjust = 2.5))

        if(parameter_size=="Reverse" & is.null(keyvalssize)==FALSE){
          lolipop_plot <- lolipop_plot + scale_size(trans = 'reverse')
        }
        if(Flip==TRUE){
          lolipop_plot <- lolipop_plot + coord_flip() +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
        }
        if("label_dot" %in% names(Plot_SettingsInfo)==TRUE ){
          lolipop_plot <-   lolipop_plot + geom_text(color="black", size=2)
        }
      }



      plot(lolipop_plot)


      if(OutputPlotName ==""){
        ggsave(file=paste(Results_folder_plots_Lolipop_folder,"/", "Lolipop", OutputPlotName, ".",Save_as_Plot, sep=""), plot=lolipop_plot, width=10, height=10)
      }else{
        ggsave(file=paste(Results_folder_plots_Lolipop_folder,"/", "Lolipop_", OutputPlotName, ".",Save_as_Plot, sep=""), plot=lolipop_plot, width=10, height=10)
      }
    }
  }else if(Plot_Settings=="Compare"){
    if("individual" %in% names(Plot_SettingsInfo)==TRUE){

      Combined_Input <- data.frame(matrix(ncol = 4, nrow = 0))
      comb.colnames <- c("Metabolite",paste(x), "Condition")
      colnames(Combined_Input) <- comb.colnames

      for (i in 1:length(Input_data)){
        Input_data[[i]]$Condition <- Comparison_name[[i]]
        data <-  Input_data[[i]]
        if(is.null(Plot_SettingsFile)==FALSE){
          data <- merge(Input_data[[i]], Plot_SettingsFile_List[[i]], by = "Metabolite")
        }
        Combined_Input <- rbind(Combined_Input, data)

      }

      if("label_dot" %in% names(Plot_SettingsInfo)==TRUE ){
        Plot_SettingsInfo[["label_dot"]]
        Combined_Input[ Plot_SettingsInfo[["label_dot"]]] <- round(Combined_Input[stat], digits = 4)
      }
      Combined_Input <- Combined_Input[is.finite(Combined_Input[[x]]),]


      # Create the list of individual plots that should be made:
      IndividualPlots <- Combined_Input[!duplicated(Combined_Input$individual),]
      IndividualPlots <- IndividualPlots$individual

      PlotList <- list()#Empty list to store all the plots

      for (i in IndividualPlots){
        # i = IndividualPlots[1]
        Plot_SettingsInfo_indi <- Plot_SettingsInfo
        # Plot_SettingsFile_Select <- subset(Combined_Input, individual == paste(i))
        # InputLolipop  <- merge(x=Plot_SettingsFile_Select,y=Input_data, by="Metabolite", all.x=TRUE)%>%
        #   na.omit()
        InputLolipop <- subset(Combined_Input, individual == paste(i))


        #Select metabolites for the cut offs selected
        loli.data <- InputLolipop %>% mutate(names=Metabolite)



        if("size" %in% names(Plot_SettingsInfo_indi)==TRUE ){
          if(is.numeric(loli.data$size)==FALSE){ # run is color is discrete
            stop("Size can take only numeric values")
          }else{# color = continuous
            keyvalssize <- loli.data$size
          }
        } else{
          Plot_SettingsInfo_indi= c(Plot_SettingsInfo_indi,size="p.adj")
          keyvalssize <- loli.data$size
        }

        if("label_dot" %in% names(Plot_SettingsInfo)==TRUE ){
          Combined_Input[ Plot_SettingsInfo[["label_dot"]]] <- round(Combined_Input[stat], digits = 4)
          label <-  Plot_SettingsInfo[["label_dot"]]
        }else{
          label <-  ""
        }

        if(Flip==TRUE){
          lolipop_plot <- ggplot(loli.data , aes(x = get(x), y = names , label=get(label))) +
            geom_segment(aes(x = 0, xend = get(x), y = names, yend = names)) +
            geom_point(aes(colour = Condition, size = keyvalssize ))   +
            #scale_size_continuous(range = c(1,5))+
            geom_vline(xintercept = 0) +
            coord_flip()+
            Theme+
            theme(plot.title = element_text(color = "black", size = 12, face = "bold"),
                  plot.subtitle = element_text(color = "black", size=10),
                  plot.caption = element_text(color = "black",size=9, face = "italic", hjust = 2.5),
                  axis.text.x = element_text(angle = 90, hjust = 1))+
            ylab(y)+
            xlab(x)+
            labs(size=Plot_SettingsInfo_indi[['size']])  +
            labs(title = paste(OutputPlotName,": ", i ),subtitle = Subtitle)

        }else{
          lolipop_plot <- ggplot(loli.data , aes(x = get(x), y = names , label=!!label)) +
            geom_segment(aes(x = 0, xend = get(x), y = names, yend = names)) +
            geom_point(aes(colour = Condition, size = keyvalssize ))   +
            # scale_size_continuous(range = c(1,5))+# , trans = 'reverse') +
            #   scale_colour_gradient(low = "red", high = "blue")+#, limits = c(0, max(loli.data[,stat]))) +
            geom_vline(xintercept = 0) +
            Theme+
            theme(plot.title = element_text(color = "black", size = 12, face = "bold"),
                  plot.subtitle = element_text(color = "black", size=10),
                  plot.caption = element_text(color = "black",size=9, face = "italic", hjust = 2.5))+
            ylab(x)+
            xlab(y)+
            # #  labs(color=Plot_SettingsInfo_indi[['color']]) +
            labs(size=Plot_SettingsInfo_indi[['size']])  +
            labs(title = paste(OutputPlotName,": ", i ),subtitle = Subtitle)
        }
        if(parameter_size=="Reverse" & is.null(keyvalssize)==FALSE){
          lolipop_plot <-   lolipop_plot + scale_size(trans = 'reverse')
        }
        if("label_dot" %in% names(Plot_SettingsInfo)==TRUE ){
          lolipop_plot <-   lolipop_plot + geom_text(color="black", size=2)
        }
        lolipop_plot


        #save plot and get rid of extra signs before saving
        cleaned_i <- gsub("[[:space:],/\\\\]", "-", i)#removes empty spaces and replaces /,\ with -
        if(OutputPlotName ==""){
          ggsave(file=paste(Results_folder_plots_Lolipop_folder,"/", "Lolipop_",cleaned_i, ".",Save_as_Plot, sep=""), plot=lolipop_plot, width=8, height=6)
        }else{
          ggsave(file=paste(Results_folder_plots_Lolipop_folder,"/", "Lolipop_", OutputPlotName, "_",cleaned_i, ".",Save_as_Plot, sep=""), plot=lolipop_plot, width=8, height=6)
        }
        ## Store the plot in the 'plots' list
        PlotList[[cleaned_i]] <- lolipop_plot
        plot(lolipop_plot)
      }
      # Return PlotList into the environment to enable the user to view the plots directly
      #assign("LolipopPlots", PlotList, envir=.GlobalEnv)
      # Combine plots into a single plot using facet_grid or patchwork::wrap_plots
      Return <- PlotList

    }
    else if("individual" %in% names(Plot_SettingsInfo)==FALSE){
      Combined_Input <- data.frame(matrix(ncol = 4, nrow = 0))
      comb.colnames <- c("Metabolite",paste(x), "Condition")
      colnames(Combined_Input) <- comb.colnames

      for (i in 1:length(Input_data)){
        Input_data[[i]]$Condition <- Comparison_name[[i]]
        data <-  Input_data[[i]]
        if(is.null(Plot_SettingsFile)==FALSE){
          data <- merge(Input_data[[i]], Plot_SettingsFile_List[[i]], by = "Metabolite")
        }
        Combined_Input <- rbind(Combined_Input, data)

      }

      if("size" %in% names(Plot_SettingsInfo)==TRUE ){
        if(is.numeric(Combined_Input$size)==FALSE){ # run is color is discrete
          stop("Size can take only numeric values")
        }else{# color = continuous
          keyvalssize <- Combined_Input$size
        }
      } else{
        Plot_SettingsInfo= c(Plot_SettingsInfo,size="p.adj")
        keyvalssize <- Combined_Input$size
      }

      if("label_dot" %in% names(Plot_SettingsInfo)==TRUE ){
        Combined_Input[ Plot_SettingsInfo[["label_dot"]]] <- round(Combined_Input[stat], digits = 4)
        label <-  Plot_SettingsInfo[["label_dot"]]
      }else{
        label <-  ""
      }

      # Remove the metabolite with inf in logFC because it messes the plot
      Combined_Input <- Combined_Input[is.finite(Combined_Input[[x]]),]


      if(Flip==TRUE){
        lolipop_plot <- ggplot(Combined_Input, aes(x=reorder(Metabolite, + get(x)), y=.data[[x]], label=!!label)) +
          geom_point(stat = 'identity', aes(size = keyvalssize, col = Condition))  +
          #   geom_text(color="black", size=2) +
          ylim(((Reduce(min,Combined_Input[[x]]))-0.5),((Reduce(max,Combined_Input[[x]]))+0.5)) +
          geom_hline(yintercept = 0) +
          Theme+
          theme(plot.title = element_text(color = "black", size = 12, face = "bold"),
                plot.subtitle = element_text(color = "black", size=10),
                plot.caption = element_text(color = "black",size=9, face = "italic", hjust = 2.5),
                axis.text.x = element_text(angle = 90, hjust = 1))+
          ylab(x)+
          xlab(y)+
          labs(size=Plot_SettingsInfo[['size']])  +
          labs(title = OutputPlotName,subtitle = Subtitle)

      }else{
        lolipop_plot <- ggplot(Combined_Input, aes(x=reorder(Metabolite,+ get(x)), y=.data[[x]], label=!!label)) +
          geom_point(stat = 'identity', aes(size = keyvalssize, col = Condition))  +
          #geom_text(color="black", size=2) +
          ylim(((Reduce(min,Combined_Input[[x]]))-0.5),((Reduce(max,Combined_Input[[x]]))+0.5)) +
          coord_flip()+
          geom_hline(yintercept = 0) +
          Theme+
          theme(plot.title = element_text(color = "black", size = 12, face = "bold"),
                plot.subtitle = element_text(color = "black", size=10),
                plot.caption = element_text(color = "black",size=9, face = "italic", hjust = 2.5))+
          ylab(y)+
          xlab(x)+
          labs(size=Plot_SettingsInfo[['size']])  +
          labs(title = OutputPlotName,subtitle = Subtitle)
      }
      if(parameter_size=="Reverse" & is.null(keyvalssize)==FALSE){
        lolipop_plot <-   lolipop_plot + scale_size(trans = 'reverse')
      }

      if("label_dot" %in% names(Plot_SettingsInfo)==TRUE ){
        lolipop_plot <-   lolipop_plot + geom_text(color="black", size=2)
      }




      if(OutputPlotName ==""){
        ggsave(file=paste(Results_folder_plots_Lolipop_folder,"/", "Lolipop", ".",Save_as_Plot, sep=""), plot=lolipop_plot, width=12, height=14)
      }else{
        ggsave(file=paste(Results_folder_plots_Lolipop_folder,"/", "Lolipop_", OutputPlotName, ".",Save_as_Plot, sep=""), plot=lolipop_plot, width=12, height=14)
      }
      plot(lolipop_plot)

    }


  }else if(Plot_Settings=="PEA"){# Code Missing
    Input_data <- Input_data[[1]]

    # Create the list of individual plots that should be made:
    IndividualPlots <- Plot_SettingsFile[!duplicated(Plot_SettingsFile$individual),]
    IndividualPlots <- IndividualPlots$individual

    PlotList <- list()#Empty list to store all the plots

    for (i in IndividualPlots){
      # i <- IndividualPlots[1]
      Plot_SettingsInfo_indi <- Plot_SettingsInfo
      Plot_SettingsFile_Select <- subset(Plot_SettingsFile, individual == paste(i))
      InputLolipop  <- merge(x=Plot_SettingsFile_Select,y=Input_data, by="Metabolite", all.x=TRUE)%>%
        na.omit()

      AdditionalInput_data_Select<- subset(AdditionalInput_data, PEA_Pathway == paste(i)) #Select pathway we plot and use the score and stats

      #Select metabolites for the cut offs selected
      loli.data <- InputLolipop %>% mutate(names=Metabolite)


      if("size" %in% names(Plot_SettingsInfo_indi)==TRUE ){
        if(is.numeric(loli.data$size)==FALSE){ # run is color is discrete
          stop("Size can take only numeric values")
        }else{# color = continuous
          keyvalssize <- loli.data$size
        }
      } else{
        Plot_SettingsInfo_indi= c(Plot_SettingsInfo_indi,size="p.adj")
        keyvalssize <- loli.data$size
      }


      if("label_dot" %in% names(Plot_SettingsInfo)==TRUE ){
        label <-  Plot_SettingsInfo[["label_dot"]]
        loli.data[ Plot_SettingsInfo[["label_dot"]]] <- round(loli.data[[label]], digits = 3)
      }else{
        label <-  ""
      }

      p2 <- NULL
      # add color here
      if("color" %in% names(Plot_SettingsInfo_indi)==TRUE ){
        if(is.numeric(loli.data$color)==FALSE){ # run if color is discrete

          col_var_name <- Plot_SettingsInfo_indi[['color']]

          position <- which(names(loli.data)=="color" )
          names(loli.data)[position]<-"plot_color_variable"


          loli.data <- loli.data %>%
            arrange(plot_color_variable, get(x),Metabolite)

          loli.data_avg <- loli.data %>%
            arrange(plot_color_variable, get(x), Metabolite) %>%
            mutate(Metab_name = row_number()) %>%
            group_by(plot_color_variable) %>%
            mutate(
              avg = mean(get(x))
            ) %>%
            ungroup() %>%
            mutate(plot_color_variable = factor(plot_color_variable))


          loli_lines <-   loli.data_avg %>%
            arrange(plot_color_variable, Metabolite) %>%
            group_by(plot_color_variable) %>%
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
            ggplot(aes(Metab_name, get(x))) + # names in aes ro Metab_name
            geom_hline(
              data = tibble(y = -5:5),
              aes(yintercept = y),
              color = "grey82",
              size = .5 )

          p2 <- p2 + geom_segment(
            aes(
              xend = Metab_name,          # names
              yend = 0,#avg,
              color = plot_color_variable,
              #color = after_scale(colorspace::lighten(color, .2))
            ))

          p2 <- p2 + # geom_line( data = loli_lines, aes(x, y),  color = "grey40"  ) +
            geom_line(
              data = loli_lines,
              aes( x_group, y,
                   color = plot_color_variable,
                   #  color = after_scale(colorspace::darken(color, .2))
              ), size = 2.5) +  geom_point(aes(size = keyvalssize, color = plot_color_variable)
              )

          p2 <- p2 +theme(axis.text.y=element_blank(),
                          axis.ticks.y=element_blank()
          )

          if(Flip == TRUE){
            p2 <- p2 +theme(axis.text.x=element_blank(),
                            axis.ticks.x=element_blank())
            lab_pos_metab <- loli.data_avg[loli.data_avg[x]>0,]  %>% select(Metabolite, Metab_name, x)

            p2 <- p2+ annotate("text", x = lab_pos_metab$Metab_name, y = lab_pos_metab[[x]]+1.5, label = lab_pos_metab$Metabolite,angle = 90, size = 3)

            lab_neg_metab <-  loli.data_avg[loli.data_avg[x]<0,]  %>% select(Metabolite, Metab_name, x)
            p2<- p2+ annotate("text", x = lab_neg_metab$Metab_name, y = lab_neg_metab[[x]]-1.5, label = lab_neg_metab$Metabolite, angle = 90,size = 3)

            if("label_dot" %in% names(Plot_SettingsInfo)==TRUE ){
              dot_pos_metab <- loli.data_avg[loli.data_avg[x]>0,]  %>% select(Metabolite, Metab_name, x, label)
              p2<- p2+ annotate("text", x = dot_pos_metab$Metab_name, y = dot_pos_metab[[x]], label = dot_pos_metab[[label]], size = 3)

              dot_neg_metab <- loli.data_avg[loli.data_avg[x]<0,]  %>% select(Metabolite, Metab_name, x,label)
              p2<- p2+ annotate("text", x = dot_neg_metab$Metab_name, y = dot_neg_metab[[x]], label = dot_neg_metab[[label]], size = 3)

            }

            p2 <- p2+Theme
            p2 <- p2+ labs(color=col_var_name)+
              labs(size=Plot_SettingsInfo[['size']])

            p2 <- p2  +  labs(title = paste(OutputPlotName, ": ", i, sep=""),
                              subtitle = paste(Plot_SettingsInfo[["PEA_score"]],"= ", AdditionalInput_data_Select$PEA_score, ", ",Plot_SettingsInfo[["PEA_stat"]] , "= ", AdditionalInput_data_Select$PEA_stat, sep=""),
                              caption = paste0("total = ", nrow(InputLolipop), " metabolites of ", nrow(Plot_SettingsFile_Select), " metabolites in pathway"))+
              theme(plot.title = element_text(color = "black", size = 12, face = "bold"),
                    plot.subtitle = element_text(color = "black", size=10),
                    plot.caption = element_text(color = "black",size=9, face = "italic", hjust = 2.5))

          }else{
            p2<- p2 + coord_flip()
            p2 <- p2 +theme(axis.text.x=element_blank(),
                            axis.ticks.x=element_blank())
            lab_pos_metab <- loli.data_avg %>% filter(get(x)>0) %>% select(Metabolite, Metab_name, get(x))
            p2<- p2+ annotate("text", x = lab_pos_metab$Metab_name, y = lab_pos_metab[[x]]+1.5, label = lab_pos_metab$Metabolite, size = 3)

            lab_neg_metab <- loli.data_avg %>% filter(get(x)<0) %>% select(Metabolite, Metab_name, get(x))
            p2<- p2+ annotate("text", x = lab_neg_metab$Metab_name, y = lab_neg_metab[[x]]-1.5, label = lab_neg_metab$Metabolite, size = 3)

            if("label_dot" %in% names(Plot_SettingsInfo)==TRUE ){
              dot_pos_metab <- loli.data_avg[loli.data_avg[x]>0,]  %>% select(Metabolite, Metab_name, x, label)
              p2<- p2+ annotate("text", x = dot_pos_metab$Metab_name, y = dot_pos_metab[[x]], label = dot_pos_metab[[label]], size = 3)

              dot_neg_metab <- loli.data_avg[loli.data_avg[x]<0,]  %>% select(Metabolite, Metab_name, x,label)
              p2<- p2+ annotate("text", x = dot_neg_metab$Metab_name, y = dot_neg_metab[[x]], label = dot_neg_metab[[label]], size = 3)

            }

            p2 <- p2+Theme

            p2 <- p2  +  labs(title = paste(OutputPlotName, ": ", i, sep=""),
                              subtitle = paste(Plot_SettingsInfo[["PEA_score"]],"= ", AdditionalInput_data_Select$PEA_score, ", ",Plot_SettingsInfo[["PEA_stat"]] , "= ", AdditionalInput_data_Select$PEA_stat, sep=""),
                              caption = paste0("total = ", nrow(InputLolipop), " metabolites of ", nrow(Plot_SettingsFile_Select), " metabolites in pathway"))+
              theme(plot.title = element_text(color = "black", size = 12, face = "bold"),
                    plot.subtitle = element_text(color = "black", size=10),
                    plot.caption = element_text(color = "black",size=9, face = "italic", hjust = 2.5))

            p2 <- p2+ labs(color=col_var_name)+
              labs(size=Plot_SettingsInfo_indi[['size']])
          }

          lolipop_plot <-  p2
          if(parameter_size=="Reverse" & is.null(keyvalssize)==FALSE){
            lolipop_plot <-   lolipop_plot + scale_size(trans = 'reverse')
          }

          # Put back the correct name in the data df
          names(loli.data)[position]<- col_var_name

        }else{# color = continuous
          keyvals <- loli.data$color
        }
      }else{
        Plot_SettingsInfo_indi= c(Plot_SettingsInfo_indi, color="p.adj")
        keyvals <- loli.data$color
      }

      if(is.null(p2)==TRUE){

        loli.data$names <- as.factor(loli.data$names)
        loli.data[[x]]<- as.numeric(loli.data[[x]])
        loli.data$names <- reorder(loli.data$names, -loli.data[[x]])

        lolipop_plot <- ggplot(loli.data , aes(x = get(x), y = names,label=!!label)) +
          geom_segment(aes(x = 0, xend = get(x), y = names, yend = names)) +
          geom_point(aes(colour = keyvals, size = keyvalssize ))   +
          # scale_size_continuous(range = c(1,5))+# , trans = 'reverse') +

          scale_colour_gradient(low = "red", high = "blue")+
          geom_vline(xintercept = 0)+
          Theme+
          labs(color=Plot_SettingsInfo_indi[['color']])+
          labs(size=Plot_SettingsInfo_indi[['size']]) +
          ylab(y)+
          xlab(x)+
          labs(title = paste(OutputPlotName, ": ", i, sep=""),
               subtitle = paste(Plot_SettingsInfo[["PEA_score"]],"= ", AdditionalInput_data_Select$PEA_score, ", ",Plot_SettingsInfo[["PEA_stat"]] , "= ", AdditionalInput_data_Select$PEA_stat, sep=""),
               caption = paste0("total = ", nrow(InputLolipop), " metabolites of ", nrow(Plot_SettingsFile_Select), " metabolites in pathway"))+
          theme(plot.title = element_text(color = "black", size = 12, face = "bold"),
                plot.subtitle = element_text(color = "black", size=10),
                plot.caption = element_text(color = "black",size=9, face = "italic", hjust = 2.5))

        if(parameter_size=="Reverse" & is.null(keyvalssize)==FALSE){
          lolipop_plot <-   lolipop_plot + scale_size(trans = 'reverse')
        }
        if(Flip==TRUE){
          lolipop_plot <- lolipop_plot + coord_flip() +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
        }
        if("label_dot" %in% names(Plot_SettingsInfo)==TRUE ){
          lolipop_plot <-   lolipop_plot + geom_text(color="black", size=2)
        }
        plot(lolipop_plot)
      }

      #save plot and get rid of extra signs before saving
      cleaned_i <- gsub("[[:space:],/\\\\]", "-", i)#removes empty spaces and replaces /,\ with -
      if(OutputPlotName ==""){
        ggsave(file=paste(Results_folder_plots_Lolipop_folder,"/", "Lolipop_",cleaned_i, ".",Save_as_Plot, sep=""), plot=lolipop_plot, width=8, height=6)
      }else{
        ggsave(file=paste(Results_folder_plots_Lolipop_folder,"/", "Lolipop_", OutputPlotName, "_",cleaned_i, ".",Save_as_Plot, sep=""), plot=lolipop_plot, width=8, height=6)
      }
      ## Store the plot in the 'plots' list
      PlotList[[cleaned_i]] <- lolipop_plot
      plot(lolipop_plot)
    }
    # Return PlotList into the environment to enable the user to view the plots directly
    #assign("LolipopPlots", PlotList, envir=.GlobalEnv)
    # Combine plots into a single plot using facet_grid or patchwork::wrap_plots
    Return <- PlotList

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
#' @param Plot_SettingsInfo  \emph{Optional: } NULL or Named vector  where you can include vectors or lists for annotation c(individual= "", color_Metab="ColumnName_Plot_SettingsFile_Metab", color_Sample= list("ColumnName_Plot_SettingsFile_Sample", "ColumnName_Plot_SettingsFile_Sample",...)).\strong{Default = NULL}
#' @param Plot_SettingsFile_Sample DF which contains information about the samples, which will be combined with your input data based on the unique sample identifiers. and other columns with required PlotSettingInfo.\strong{Default = NULL}
#' @param Plot_SettingsFile_Metab  \emph{Optional: } DF with column "Metabolite" including the Metabolite names (needs to match Metabolite names of Input_data) and other columns with required PlotSettingInfo. \strong{Default = NULL}
#' @param OutputPlotName \emph{Optional: } String which is added to the output files of the plot
#' @param SCALE \emph{Optional: } String with the information for scale row or column. \strong{Default = row}
#' @param Save_as_Plot \emph{Optional: } Select the file type of output plots. Options are svg, pdf or png. \strong{Default = "svg"}
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
                       Save_as_Plot = "svg"
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

    # Check if the next lines work correctly in case of duplicated metabolties (colnames)

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
          if( metab_color %in% colnames(Plot_SettingsFile_Metab)== FALSE){
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
          if( samp_color %in% colnames(Plot_SettingsFile_Sample)== FALSE){
            stop("You have chosen color_Sample = ",paste(samp_color), ", ", paste(samp_color)," does not exist in the Plot_SettingsFile_Sample"   )
          }
        }
      }
    }
  }

<<<<<<< HEAD
=======
  # 4. Check other plot-specific parameters:
  Save_as_Plot_options <- c("svg","pdf","png")
  if(Save_as_Plot %in% Save_as_Plot_options == FALSE){
    stop("Check input. The selected Save_as_Plot option is not valid. Please select one of the following: ",paste(Save_as_Plot_options,collapse = ", "),"." )
  }

>>>>>>> 8350d1eb04353fa3a28581f9478e5a968b19445d
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
      Features <- as.data.frame(t(data_path ))
      if(nrow(Features)>200){
        show_rownames <- FALSE
      }else{
        show_rownames <- TRUE
      }

      # Make the plot
      set.seed(1234)
      heatmap <- pheatmap::pheatmap(t(data_path),
                                    show_rownames = as.logical(show_rownames),
                                    clustering_method =  "complete",
                                    scale = SCALE,
                                    clustering_distance_rows = "correlation",
                                    annotation_col = col_annot,
                                    annotation_row = row_annot,
                                    main = paste(OutputPlotName, i, sep=" " ))

      #Width and height according to Sample and metabolite number
      # 1. Get the legend dimensions from the pheatmap object
      legend_width <- as.numeric(regmatches(heatmap$gtable$widths[4], gregexpr("[0-9.]+", heatmap$gtable$widths[4]))[[1]])
      width_cm <- sum((grid::convertX(unit(legend_width[1], "npc"), "cm", valueOnly = TRUE)),(grid::convertX(unit(legend_width[2], "bigpts"), "cm", valueOnly = TRUE)))
      legend_height <- as.numeric(regmatches(heatmap$gtable$heights[5], gregexpr("[0-9.]+", heatmap$gtable$heights[5]))[[1]])
      height_cm <- sum((grid::convertX(unit(legend_height[1], "npc"), "cm", valueOnly = TRUE)),(grid::convertX(unit(legend_height[2], "bigpts"), "cm", valueOnly = TRUE)))/2

      #2. Give value to width and heights
      plot_width <- (nrow(data_path ) * 0.2) + width_cm
      plot_height <- (nrow(Features) * 0.3) + height_cm

      #Save
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
    if(nrow(Features)>200){
      show_rownames <- FALSE
    }else{
      show_rownames <- TRUE
    }

    #Make the plot:
    set.seed(1234)
    heatmap <- pheatmap::pheatmap(t(data),
                                  show_rownames = as.logical(show_rownames),
                                  clustering_method =  "complete",
                                  scale = SCALE,
                                  clustering_distance_rows = "correlation",
                                  annotation_col = col_annot,
                                  annotation_row = row_annot,
                                  main = OutputPlotName)

    #Width and height according to Sample and metabolite number
    # 1. Get the legend dimensions from the pheatmap object
    legend_width <- as.numeric(regmatches(heatmap$gtable$widths[4], gregexpr("[0-9.]+", heatmap$gtable$widths[4]))[[1]])
    width_cm <- sum((grid::convertX(unit(legend_width[1], "npc"), "cm", valueOnly = TRUE)),(grid::convertX(unit(legend_width[2], "bigpts"), "cm", valueOnly = TRUE)))
    legend_height <- as.numeric(regmatches(heatmap$gtable$heights[5], gregexpr("[0-9.]+", heatmap$gtable$heights[5]))[[1]])
    height_cm <- sum((grid::convertX(unit(legend_height[1], "npc"), "cm", valueOnly = TRUE)),(grid::convertX(unit(legend_height[2], "bigpts"), "cm", valueOnly = TRUE)))/2

    #2. Give value to width and heights
    plot_width <- (nrow(data) * 0.2) + width_cm
    plot_height <- (nrow(Features) * 0.3) + height_cm

    #Save
    ggsave(file=paste(Results_folder_plots_Heatmaps_folder,"/", "Heatmap",OutputPlotName, ".", Save_as_Plot ,sep=""), plot=heatmap, width=plot_width, height= plot_height, units = "cm")
  }
}



##############################
### ### ### Plots  ### ### ###
##############################

#' @param Input_data DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Includes experimental design and outlier column.
#' @param Experimental_design DF which contains information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames. Column "Conditions" with information about the sample conditions (e.g. "N" and "T" or "Normal" and "Tumor"), can be used for feature filtering and colour coding in the PCA. Column "AnalyticalReplicate" including numerical values, defines technical repetitions of measurements, which will be summarised. Column "BiologicalReplicates" including numerical values. Please use the following names: "Conditions", "Biological_Replicates", "Analytical_Replicates".
#' @param Graprh_Style String with the information of the Graph style. Available options are Bar. Box and Violin  \strong{Default = Box}
#' @param Superplot \emph{Optional: } String with a Column name of the Experimental_design as string which is used to make the plots Superplots.
#' @param Output_Name \emph{Optional: } String which is added to the output files of the plot.
#' @param Output_plots String with plot save information. Available options are "Individual" for plots of each Individual metabolite and "Together" for a pdf containing all the plots. \strong{Default = Together}
#' @param Selected_Conditions Vector with names of selected Conditions for the plot. \strong{Default = NULL}
#' @param Selected_Comparisons Logical, TRUE to use t.test between the Selected_Conditions or FALSE. \strong{Default = NULL}
#' @param Theme \emph{Optional: } Selection of theme for plot. \strong{Default = theme_classic} ??
#' @param Save_as_Plot \emph{Optional: } Select the file type of output plots. Options are svg, pdf or png. \strong{Default = svg}
#'
#' @keywords Barplot, Boxplot, Violinplot, Superplot
#' @export

VizPlots <- function(Input_data,
                     Experimental_design,
                     Graprh_Style = "Box", # Bar, Box, Violin
                     Superplot = NULL,
                     OutputPlotName = "",
                     Output_plots = "Together",
                     Selected_Conditions = NULL,
                     Selected_Comparisons = NULL,
                     Theme = theme_classic(),
                     Save_as_Plot = "svg"
){

  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse", "ggplot2", "ggpubr", "ggbeeswarm")
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

  #2. Plot options
  if(Graprh_Style %in% c("Box", "Bar", "Violin") == FALSE){
    stop("Graprh_Style must be either Box, Bar ot Violin.")
  }
  if(is.null(Superplot)==FALSE){
    if(Superplot %in% colnames(Experimental_design) == FALSE){
      stop("Superplot is active. However, the column ", paste(Superplot), " is not the in Experimental_design." )
    }
  }

  Output_plots_options <- c("Individual", "Together")
  if (Output_plots %in% Output_plots_options == FALSE){
    stop("Check Input the Plot_pathways option is incorrect. The Allowed options are the following: ",paste(Output_plots_options,collapse = ", "),"." )
  }
  if("Conditions" %in% colnames(Experimental_design)==FALSE){
    stop("There is no column named `Conditions` in Input_data.")
  }

  #3. Comparison options
  if(is.null(Selected_Conditions)==FALSE){
    for (Conditions in Selected_Conditions){
      if(Conditions %in% Experimental_design$Conditions==FALSE){
        stop("Check Input. The Selected_Conditions were not found in the Conditions Column.")
      }
    }
  }

  if(is.null(Selected_Comparisons)==FALSE){
    for (Comp in Selected_Comparisons){
      if(Comp[1] %in% Experimental_design$Conditions ==FALSE){
        stop("Check Input. The Selected_Comparisons condition ",paste(Comp[1]), " is not found in the Conditions Column of the Experimental_design.")
      }
      if(Comp[2] %in% Experimental_design$Conditions ==FALSE){
        stop("Check Input. The Selected_Comparisons condition ",paste(Comp[2]), " is not found in the Conditions Column of the Experimental_design.")
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
  Results_folder_plots_Barplots_folder = file.path(Results_folder, paste(Graprh_Style, "plots", sep=""))
  if (!dir.exists(Results_folder_plots_Barplots_folder)) {dir.create(Results_folder_plots_Barplots_folder)}  # check and create folder



  Metabolite_Names <- colnames(data)
  data <- merge( Experimental_design[c("Conditions",Superplot)] ,data,, by=0)
  data <- column_to_rownames(data, "Row.names")

  # make a list for plotting all plots together
  outlier_plot_list <- list()
  k=1

  for (i in Metabolite_Names){
    # i = Metabolite_Names[1]

    dataMeans <- data %>%  select(i, Conditions) %>% group_by(Conditions) %>% summarise_at(vars(i), list(mean = mean, sd = sd)) %>% as.data.frame()
    names(dataMeans)[2] <- "Intensity"
    plotdata <- data %>%  select(i,Conditions, Superplot) %>%  group_by(Conditions)  %>% as.data.frame()
    names(plotdata)[1] <- c("Intensity")

    # Take only selected conditions
    if (is.null(Selected_Conditions) == "FALSE"){
      dataMeans <- dataMeans %>% filter(Conditions %in% Selected_Conditions)
      plotdata <- plotdata %>% filter(Conditions %in% Selected_Conditions)
    }


    # Make the Plot
    Plot <- ggplot(plotdata, aes(x = factor(Conditions), y = Intensity))#,fill = Superplot)

    # Add graph style
    if (Graprh_Style == "Bar"){
      Plot <- Plot+
        geom_bar(stat = "summary", fun = "mean", fill = "skyblue") +
        geom_errorbar(data = dataMeans, aes(x=Conditions, ymin=Intensity-sd, ymax=Intensity+sd), width=0.4, colour="black", alpha=0.9, size=0.5)
    } else if (Graprh_Style == "Violin"){
      Plot <- Plot+
        geom_violin(fill = "skyblue")  +
        geom_errorbar(data = dataMeans, aes(x=Conditions, ymin=Intensity-sd, ymax=Intensity+sd), width=0.4, colour="black", alpha=0.9, size=0.5)
    } else if (Graprh_Style == "Box"){
      Plot <- Plot +
        geom_boxplot(fill="skyblue")
    }

    # Add superplot
    if(is.null(Superplot)==FALSE){
      Plot <- Plot+ ggbeeswarm::geom_beeswarm(aes(x=Conditions,y=Intensity,color=as.factor(get(Superplot))),size=3)
    }


    if(is.null(Selected_Comparisons)== TRUE){
      Plot <- Plot + geom_errorbar(data = dataMeans, aes(x=Conditions, ymin=Intensity-sd, ymax=Intensity+sd), width=0.4, colour="black", alpha=0.9, size=0.5)
    }else{
      Plot <- Plot+ ggpubr::stat_compare_means(comparisons = Selected_Comparisons,
                                               label = "p.format", method = "t.test", hide.ns = TRUE, position = position_dodge(0.9), vjust = 0.25, show.legend = FALSE)
    }



    Plot <- Plot + theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "right")+xlab("Conditions")+ ylab("Mean Intensity")
    Plot <- Plot + Theme
    Plot <- Plot + ggtitle(paste(i))

    Plot

    if(Output_plots=="Individual"){

      i <- (gsub("/","_",i))#remove "/" cause this can not be safed in a PDF name
      i <- (gsub(":","_",i))
      i <-(gsub("\\*","",i))

      if(OutputPlotName ==""){
        ggsave(file=paste(Results_folder_plots_Barplots_folder, "/",i, ".",Save_as_Plot, sep=""), plot=Plot, width=10, height=8)
      }else{
        ggsave(file=paste(Results_folder_plots_Barplots_folder, "/",OutputPlotName,"_",i, ".",Save_as_Plot, sep=""), plot=Plot, width=10, height=8)
      }


    } else if(Output_plots=="Together"){

      plot(Plot)
      outlier_plot_list[[k]] <- recordPlot()
      dev.off()
      k=k+1
    }
  }

  if(Output_plots=="Together"){
    if(OutputPlotName ==""){
      pdf(file= paste(Results_folder_plots_Barplots_folder,"/",Graprh_Style, "plots", OutputPlotName,".pdf", sep = ""), onefile = TRUE )
    }else{
      pdf(file= paste(Results_folder_plots_Barplots_folder,"/",Graprh_Style, "plots", OutputPlotName,".pdf", sep = ""), onefile = TRUE )
      for (plot in outlier_plot_list){
        replayPlot(plot)
      }
      dev.off()
    }
  }
}



