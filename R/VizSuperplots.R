## ---------------------------
##
## Script name: Visualization Superplots
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
#' This script allows you to perform different visualizations (bar, box, violin plots) using the results of the MetaProViz analysis


#' @param Input_data DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Includes experimental design and outlier column.
#' @param Input_SettingsFile DF which contains information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames. Column "Conditions" with information about the sample conditions (e.g. "N" and "T" or "Normal" and "Tumor"), can be used for feature filtering and colour coding in the PCA. Column "AnalyticalReplicate" including numerical values, defines technical repetitions of measurements, which will be summarised. Column "BiologicalReplicates" including numerical values. Please use the following names: "Conditions", "Biological_Replicates", "Analytical_Replicates".
#' @param Input_SettingsInfo Named vector including at least information on the conditions column: c(conditions="ColumnName_Plot_SettingsFile"). Additionally superplots can be made by adding superplot ="olumnName_Plot_SettingsFile", which are ususally biological replicates or patient IDs. \strong{Default = c(conditions="Conditions", superplot = NULL)}
#' @param Graph_Style String with the information of the Graph style. Available options are Bar. Box and Violin  \strong{Default = Box}
#' @param Output_Name \emph{Optional: } String which is added to the output files of the plot.
#' @param Individual_plots \emph{Optional: }  Logical, TRUE to save each plot individually. \strong{Default = FALSE}
#' @param Selected_Conditions Vector with names of selected Conditions for the plot. \strong{Default = NULL}
#' @param Selected_Comparisons Logical, TRUE to use t.test between the Selected_Conditions or FALSE. \strong{Default = NULL}
#' @param Theme \emph{Optional: } Selection of theme for plot. \strong{Default = theme_classic}
#' @param palette \emph{Optional: } Provide customiced color-palette in vector format. \strong{Default = NULL}
#' @param Save_as_Plot \emph{Optional: } Select the file type of output plots. Options are svg, pdf, png or NULL. \strong{Default = svg}
#'
#' @keywords Barplot, Boxplot, Violinplot, superplot
#' @export

Vizsuperplot <- function(Input_data,
                     Input_SettingsFile,
                     Input_SettingsInfo = c(conditions="Conditions", superplot = NULL),
                     Graph_Style = "Box", # Bar, Box, Violin
                     OutputPlotName = "",
                     Selected_Conditions = NULL,
                     Selected_Comparisons = NULL,
                     Individual_plots = FALSE,
                     Theme = theme_classic(),
                     palette = NULL,
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
  } else if("Conditions" %in% colnames(Input_SettingsFile)==FALSE){
    stop("There is no column named `Conditions` in Input_SettingsFile to obtain Condition1 and Condition2.")
  } else{
    Test_num <- apply(Input_data, 2, function(x) is.numeric(x))
    if((any(Test_num) ==  FALSE) ==  TRUE){
      stop("Input_data needs to be of class numeric")
    } else{
      Test_match <- merge(Input_SettingsFile, Input_data, by.x = "row.names",by.y = "row.names", all =  FALSE) # Do the unique IDs of the "Input_data" match the row names of the "Input_SettingsFile"?
      if(nrow(Test_match) ==  0){
        stop("row.names Input_data need to match row.names Input_SettingsFile.")
      } else(
        data <- Input_data
      )
    }
    Input_SettingsFile <- Input_SettingsFile
  }

  ## ------------ Check Input SettingsInfo ----------- ##
  #2. Input_SettingsInfo
  if(Input_SettingsInfo[["conditions"]] %in% Input_SettingsInfo==TRUE){
    if(Input_SettingsInfo[["conditions"]] %in% colnames(Input_SettingsFile)== TRUE){
      Input_SettingsFile<- Input_SettingsFile%>%
        dplyr::rename("Conditions"= paste(Input_SettingsInfo[["conditions"]]) )
    }else{# if true rename to Conditions
      stop("The ",Input_SettingsInfo[["conditions"]], " column selected as conditions in Input_SettingsInfo was not found in Input_SettingsFile. Please check your input.")
    }
  }else{
    stop("You have to provide a Input_SettingsInfo for conditions.")
  }

  if(Input_SettingsInfo[["superplot"]] %in% colnames(Input_SettingsFile)== TRUE){
    if(Input_SettingsInfo[["superplot"]] %in% Input_SettingsInfo==TRUE){
      Input_SettingsFile<- Input_SettingsFile%>%
        dplyr::rename("superplot"= paste(Input_SettingsInfo[["superplot"]]) )
    }else{
      stop("The ",Input_SettingsInfo[["superplot"]], " column selected for superplot in Input_SettingsInfo was not found in Input_SettingsFile. Please check your input.")
    }
  }


  #3. Plot options
  if(Graph_Style %in% c("Box", "Bar", "Violin") == FALSE){
    stop("Graph_Style must be either Box, Bar ot Violin.")
  }


  #4. Comparison options
  if(is.null(Selected_Conditions) == FALSE){
    for (Condition in Selected_Conditions){
      if(Condition %in% Input_SettingsFile$Conditions==FALSE){
        stop(paste0("Check Input. The Selected_Conditions ",Condition," were not found in the Conditions Column."))
      }
    }
  }

  if(is.null(Selected_Comparisons)==FALSE){
    for (Comp in Selected_Comparisons){
      if(is.null(Selected_Conditions)==FALSE){
        if(Selected_Conditions[Comp[1]] %in% Input_SettingsFile$Conditions ==FALSE){
        stop("Check Input. The Selected_Comparisons condition ",paste(Comp[1]), " is not found in the Conditions Column of the Input_SettingsFile.")
          }
        if(Selected_Conditions[Comp[2]] %in% Input_SettingsFile$Conditions ==FALSE){
        stop("Check Input. The Selected_Comparisons condition ",paste(Comp[2]), " is not found in the Conditions Column of the Input_SettingsFile.")
        }
      }
      }
  }


  #5. Check other plot-specific parameters:
  if (!is.null(Save_as_Plot)) {
    Save_as_Plot_options <- c("svg","pdf","png")
    if(Save_as_Plot %in% Save_as_Plot_options == FALSE){
      stop("Check input. The selected Save_as_Plot option is not valid. Please select one of the following: ",paste(Save_as_Plot_options,collapse = ", "),"." )
    }
  }

  #6. Check palette:
  if(is.null(palette)){
    palette <- "grey"
  }else{
    if(is.null(Selected_Conditions)==TRUE){
      if(length(palette) !=length(unique(Input_SettingsFile$Conditions)) ){
        stop(paste0("The palette colors used are not 1-to-1 with the groups in the ",Input_SettingsInfo[["conditions"]] ," column. Please check your Input."))
      }
    }else{
      if(length(palette) !=length(Selected_Conditions) ){
        stop("The Selected_Conditions and the palette used are not 1-to-1. Please check your Input.")
      }
    }


    }
  data <- Input_data


  ## ------------ Create Output folders ----------- ##
  if (!is.null(Save_as_Plot)) {
    name <- paste0("MetaProViz_Results_",Sys.Date())
    WorkD <- getwd()
    Results_folder <- file.path(WorkD, name)
    if (!dir.exists(Results_folder)) {dir.create(Results_folder)} # Make Results folder
    Results_folder_plots_Barplots_folder = file.path(Results_folder, paste(Graph_Style, "plots", sep=""))
    if (!dir.exists(Results_folder_plots_Barplots_folder)) {dir.create(Results_folder_plots_Barplots_folder)}  # check and create folder

    if(Individual_plots ==TRUE){
      Results_folder_plots_Barplots_folder_Individual = file.path(Results_folder_plots_Barplots_folder, "Individual")
      if (!dir.exists(Results_folder_plots_Barplots_folder_Individual)) {dir.create(Results_folder_plots_Barplots_folder_Individual)}  # check and create folder

    }
  }


  Metabolite_Names <- colnames(data)
  data <- merge( Input_SettingsFile[c("Conditions","superplot")] ,data, by=0)
  data <- column_to_rownames(data, "Row.names")

  # make a list for plotting all plots together
  plot_list <- list()
  k=1

  for (i in Metabolite_Names){
    suppressWarnings(dataMeans <- data %>%  select(i, Conditions) %>% group_by(Conditions) %>% summarise_at(vars(i), list(mean = mean, sd = sd)) %>% as.data.frame())
    names(dataMeans)[2] <- "Intensity"
    suppressWarnings(plotdata <- data %>%  select(i,Conditions, superplot) %>%  group_by(Conditions)  %>% as.data.frame() )
    names(plotdata)[1] <- c("Intensity")
    # Make conditions a factor
    plotdata$Conditions <- factor(plotdata$Conditions)

    # Take only selected conditions
    if (is.null(Selected_Conditions) == FALSE){
      dataMeans <- dataMeans %>% filter(Conditions %in% Selected_Conditions)
      plotdata <- plotdata %>% filter(Conditions %in% Selected_Conditions)
      plotdata$Conditions <- factor(plotdata$Conditions, levels = Selected_Conditions)
    }

    # Make the Plot
    Plot <- ggplot(plotdata, aes(x = Conditions, y = Intensity))

    # Add graph style and error bar
    data_summary <- function(x){#calculate error bar!
      m <- mean(x)
      ymin <- m-sd(x)
      ymax <- m+sd(x)
      return(c(y=m,ymin=ymin,ymax=ymax))
    }

    if (Graph_Style == "Bar"){
      Plot <- Plot+  geom_bar(stat = "summary", fun = "mean", fill = palette)+ stat_summary(fun.data=data_summary,
                                                                                            geom="errorbar", color="black", width=0.2)
    } else if (Graph_Style == "Violin"){
      Plot <- Plot+ geom_violin(fill = palette)+ stat_summary(fun.data=data_summary,
                                                              geom="errorbar", color="black", width=0.2)
    } else if (Graph_Style == "Box"){
      Plot <- Plot +  geom_boxplot(fill=palette,  width=0.5, position=position_dodge(width = 0.5))
    }

    # Add superplot
    if ("superplot" %in% names(Input_SettingsInfo)){
      Plot <- Plot+ ggbeeswarm::geom_beeswarm(aes(x=Conditions,y=Intensity,color=as.factor(superplot)),size=3)+
        labs(color=Input_SettingsInfo[["superplot"]], fill = Input_SettingsInfo[["superplot"]])
    }

    if(is.null(Selected_Comparisons)== TRUE){
      #Plot <- Plot + geom_errorbar(data = dataMeans, aes(x=Conditions, ymin=Intensity-sd, ymax=Intensity+sd), width=0.4, colour="black", alpha=0.9, size=0.5)
    }else{
      if(length(Selected_Comparisons)==1){# t-test
        Plot <- Plot+ ggpubr::stat_compare_means(comparisons = Selected_Comparisons,
                                                 label = "p.format", method = "t.test", hide.ns = TRUE,
                                               position = position_dodge(0.9), vjust = 0.25, show.legend = FALSE)
        Plot <- Plot +labs(caption = "p.val using pairwise t-test")
      }else{
      suppressMessages(Log2FCRes <- Log2FC(Input_data=data.frame("Intensity" = plotdata[,-c(2:3)]),
                            Input_SettingsFile=plotdata[,c(2:3)],
                            Input_SettingsInfo=c(conditions="Conditions"),
                            Save_as_Results = NULL,
                            Save_as_Plot=NULL,
                            Plot = FALSE))

        comparison_table <- data.frame(matrix(ncol = length(Selected_Comparisons), nrow = 2))
        for (k in seq_along(Selected_Comparisons)) {
          comparison_table[, k] <- sapply(Selected_Comparisons[[k]], function(x) Selected_Conditions[x])
          colnames(comparison_table)[k] <- paste0(Selected_Conditions[Selected_Comparisons[[k]][[1]]],"_vs_",  Selected_Conditions[Selected_Comparisons[[k]][[2]]])
        }

        STAT_C1vC2 <- AOV(Input_data=data.frame("Intensity" = plotdata[,-c(2:3)]),
                          conditions= plotdata[,c(2)],
                          Input_SettingsInfo=c(conditions="Conditions"),
                          STAT_padj="fdr",
                          Log2FC_table=Log2FCRes,
                          all_vs_all=TRUE,
                          comparisons=comparison_table)

        # make data for plot
        df <- data.frame(comparisons = names(STAT_C1vC2), stringsAsFactors = FALSE)
        split_names <- strsplit(df$comparisons, "_vs_")
        df$group1 <- sapply(split_names, function(x) x[1])
        df$group2 <- sapply(split_names, function(x) x[2])
        df$Log2FC <- sapply(STAT_C1vC2, function(x) x$Log2FC)
        df$p.adj <- round(sapply(STAT_C1vC2, function(x) x$p.adj),5)
        df$comparisons <- NULL
        df$y.position <-c(max(dataMeans$Intensity + 2*dataMeans$sd),
                          max(dataMeans$Intensity + 2*dataMeans$sd)+0.04* max(dataMeans$Intensity + 2*dataMeans$sd) ,
                          max(dataMeans$Intensity + 2*dataMeans$sd)+0.08* max(dataMeans$Intensity + 2*dataMeans$sd))
        # add stats to plot
        if(Graph_Style == "Bar"){
          Plot <- Plot +ggpubr::stat_pvalue_manual(df, hide.ns = FALSE, size = 3, tip.length = 0.01, step.increase=0.05)#http://rpkgs.datanovia.com/ggpubr/reference/stat_pvalue_manual.html
        }else{
          Plot <- Plot +ggpubr::stat_pvalue_manual(df, hide.ns = FALSE, size = 3, tip.length = 0.01, step.increase=0.01)#http://rpkgs.datanovia.com/ggpubr/reference/stat_pvalue_manual.html
        }
        Plot <- Plot +labs(caption = "p.adj using Anova and FDR")
     }
    }

    Plot <- Plot + Theme+ ggtitle(paste(i))
    Plot <- Plot + theme(legend.position = "right",plot.title = element_text(size=12, face = "bold"), axis.text.x = element_text(angle = 90, hjust = 1))+ xlab("Conditions")+ ylab("Normalized Intensity")

    # Make plot into nice format:
    Plot_Sized <-  MetaProViz:::plotGrob_Superplot(Input=Plot, Input_SettingsInfo=Input_SettingsInfo, Input_SettingsFile=Input_SettingsFile, MetaboliteName=i, Graph_Style=Graph_Style)
    Plot <- Plot_Sized[[3]]
    # First we want to convert the plot back into a ggplot object:
    Plot <- ggplot2::ggplot() +
      annotation_custom(Plot)

    if(Individual_plots==TRUE){

      i <- (gsub("/","_",i))#remove "/" cause this can not be safed in a PDF name
      i <- (gsub(":","_",i))
      i <-(gsub("\\*","",i))

      if (!is.null(Save_as_Plot)) {
        if(OutputPlotName ==""){
          ggsave(file=paste(Results_folder_plots_Barplots_folder_Individual, "/",i, ".",Save_as_Plot, sep=""), plot=Plot, width=10, height=8)
        }else{
          ggsave(file=paste(Results_folder_plots_Barplots_folder_Individual, "/",OutputPlotName,"_",i, ".",Save_as_Plot, sep=""), plot=Plot, width=10, height=8)
        }
      }
    }
    else{
      plot(Plot)
      plot_list[[k]] <- recordPlot()
      dev.off()
      k=k+1
    }
  }

  if(Individual_plots==FALSE){
    if (!is.null(Save_as_Plot)) {
      if(OutputPlotName ==""){
        pdf(file= paste(Results_folder_plots_Barplots_folder,"/",Graph_Style, "plots", OutputPlotName,".pdf", sep = ""), onefile = TRUE )
      }else{
        pdf(file= paste(Results_folder_plots_Barplots_folder,"/",Graph_Style, "plots", OutputPlotName,".pdf", sep = ""), onefile = TRUE )
        for (plot in plot_list){
          replayPlot(plot)
        }
        dev.off()
      }
    }
  }
}



#####################################################################
### ### ### Superplots helper function: Internal Function ### ### ###
#####################################################################

#' @param Input This is the ggplot object generated within the VizSuperplots function.
#' @param Input_SettingsInfo Passed to VizSuperplots
#' @param Input_SettingsFile Passed to VizSuperplots
#' @param MetaboliteName Passed to VizSuperplots
#' @param Graph_Style Passed to VizSuperplots
#'
#' @keywords PCA helper function
#' @noRd

plotGrob_Superplot <- function(Input, Input_SettingsInfo, Input_SettingsFile, MetaboliteName, Graph_Style){
  #------- Set the total heights and widths
  #we need ggplot_grob to edit the gtable of the ggplot object. Using this we can manipulate the gtable arguments directly.
  plottable <- ggplot2::ggplotGrob(Input) # Convert the plot to a gtable

  #-----widths (adapt for number of conditions)
  Number_Conditions <- Input_SettingsFile%>%
    dplyr::distinct(Conditions) %>%
    nrow()

  if(Graph_Style == "Bar"){
    plottable$widths[5] <- unit(Number_Conditions * 0.5, "cm")#controls x-axis
  }else{
   plottable$widths[5] <- unit(Number_Conditions * 1, "cm")#controls x-axis
  }

  plottable$widths[c(1)] <- unit(0.5,"cm")#controls margins --> y-axis label is there
  plottable$widths[c(4)] <- unit(2,"cm")#controls margins --> y-axis label is there
  plottable$widths[c(2,3)] <- unit(0,"cm")#controls margins --> not needed
  plottable$widths[c(6)] <- unit(1,"cm")#controls margins --> start Figure legend
  plottable$widths[c(7,8,10,11)] <- unit(0,"cm")#controls margins --> not needed
  plot_widths <- as.numeric(plottable$widths[5])+4

  plot(plottable)

  if(Input_SettingsInfo[["superplot"]] %in% Input_SettingsInfo==TRUE){#legend will be present!
    Value <- round(as.numeric(plottable$widths[9]),1) #plottable$widths[9] is a <unit/unit_v2> object and we can extract the extract the numeric part
    plot_widths <- plot_widths+Value
  }else{
    plottable$widths[c(9)] <- unit(0,"cm")
    plot_widths <- plot_widths
  }

  character_count <- nchar(MetaboliteName)#Check how much width is needed for the figure title/subtitle
  Titles_width <- (character_count*0.25)+0.8
  if(Titles_width>plot_widths){#If the title needs more space than the plot offers:
      plottable$widths[11] <- unit(Titles_width-plot_widths,"cm")#controls margins --> start Figure legend
      plot_widths <- Titles_width
      }

  #-----heigths
  plottable$heights[7] <- unit(8, "cm")#controls x-axis
  plottable$heights[c(8)] <- unit(1,"cm")#controls margins --> x-axis label
  plottable$heights[c(10)] <- unit(1,"cm")#controls margins --> Figure caption
  plottable$heights[c(6,9,11,12)] <- unit(0,"cm")#controls margins --> not needed
  plottable$heights[c(3)] <- unit(1,"cm")#controls margins --> Some space above the plot
  plottable$heights[c(1,2,4,5)] <- unit(0,"cm")#controls margins --> not needed

  if(Input_SettingsInfo[["superplot"]] %in% Input_SettingsInfo==TRUE){#legend will be present!
    #------- Legend heights
    Legend <- ggpubr::get_legend(Input) # Extract legend to adjust separately
    Legend_heights <- (round(as.numeric(Legend$heights[3]),1))+(round(as.numeric(Legend$heights[5]),1))
    if(Legend_heights>11){#If the legend requires more heights than the Plot
      Add <- (Legend_heights-11)/2
      plottable$heights[1] <- unit(Add,"cm")#controls margins --> Can be increased if Figure legend needs more space on the top
      plottable$heights[12] <- unit(Add,"cm")#controls margins --> Can be increased if Figure legend needs more space on the bottom
      plot_heights <- Legend_heights
    }else{
      plottable$heights[1] <- unit(0,"cm")#controls margins --> Can be increased if Figure legend needs more space on the top
      plottable$heights[12] <- unit(0,"cm")#controls margins --> Can be increased if Figure legend needs more space on the bottom
      plot_heights <- 10.5
    }
  }else{
    plot_heights <- 11
  }

  #plot_param <-c(plot_heights=plot_heights, plot_widths=plot_widths)
  Output<- list(plot_heights, plot_widths, plottable)
}


