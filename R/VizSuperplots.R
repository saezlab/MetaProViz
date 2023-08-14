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

##############################
### ### ### Superplots  ### ### ###
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

VizSuperplot <- function(Input_data,
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



