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
#' @param Plot_SettingsFile DF which contains information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames. Column "Conditions" with information about the sample conditions (e.g. "N" and "T" or "Normal" and "Tumor"), can be used for feature filtering and colour coding in the PCA. Column "AnalyticalReplicate" including numerical values, defines technical repetitions of measurements, which will be summarised. Column "BiologicalReplicates" including numerical values. Please use the following names: "Conditions", "Biological_Replicates", "Analytical_Replicates".
#' @param Plot_SettingsInfo Named vector including at least information on the conditions column: c(conditions="ColumnName_Plot_SettingsFile"). Additionally superplots can be made by adding superplot ="olumnName_Plot_SettingsFile", which are ususally biological replicates or patient IDs. \strong{Default = c(conditions="Conditions", superplot = NULL)}
#' @param Graph_Style String with the information of the Graph style. Available options are Bar. Box and Violin  \strong{Default = Box}
#' @param STAT_pval \emph{Optional: } String which contains an abbreviation of the selected test to calculate p.value. For one-vs-one comparisons choose t.test or wilcox.test , for one-vs-all or all-vs-all comparison choose aov (=anova) or kruskal.test \strong{Default = NULL}
#' @param STAT_padj \emph{Optional: } String which contains an abbreviation of the selected p.adjusted test for p.value correction for multiple Hypothesis testing. Search: ?p.adjust for more methods:"BH", "fdr", "bonferroni", "holm", etc.\strong{Default = NULL}
#' @param Output_Name \emph{Optional: } String which is added to the output files of the plot.
#' @param Individual_plots \emph{Optional: }  Logical, TRUE to save each plot individually. \strong{Default = FALSE}
#' @param Selected_Conditions Vector with names of selected Conditions for the plot. \strong{Default = NULL}
#' @param Selected_Comparisons List of numeric vectors containing Condition pairs to compare based on the order of the Selected_Conditions vector. \strong{Default = NULL}
#' @param Theme \emph{Optional: } Selection of theme for plot. \strong{Default = theme_classic}
#' @param color_palette \emph{Optional: } Provide customized color_palette in vector format. \strong{Default = NULL}
#' @param color_palette_dot \emph{Optional: } Provide customized color_palette in vector format. \strong{Default = NULL}
#' @param Save_as_Plot \emph{Optional: } Select the file type of output plots. Options are svg, pdf, png or NULL. \strong{Default = svg}
#'
#' @keywords Barplot, Boxplot, Violinplot, superplot
#' @export

VizSuperplot <- function(Input_data,
                     Plot_SettingsFile,
                     Plot_SettingsInfo = c(conditions="Conditions", superplot = NULL),
                     Graph_Style = "Box", # Bar, Box, Violin
                     STAT_pval =NULL,
                     STAT_padj=NULL,
                     OutputPlotName = "",
                     Selected_Conditions = NULL,
                     Selected_Comparisons = NULL,
                     Individual_plots = FALSE,
                     Theme = theme_classic(),
                     color_palette = NULL,
                     color_palette_dot=NULL,
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
  } else if(Plot_SettingsInfo[["conditions"]] %in% colnames(Plot_SettingsFile)==FALSE){
    stop("There is no column named `Conditions` in Plot_SettingsFile to obtain Condition1 and Condition2.")
  } else{
    Test_num <- apply(Input_data, 2, function(x) is.numeric(x))
    if((any(Test_num) ==  FALSE) ==  TRUE){
      stop("Input_data needs to be of class numeric")
    } else{
      Test_match <- merge(Plot_SettingsFile, Input_data, by.x = "row.names",by.y = "row.names", all =  FALSE) # Do the unique IDs of the "Input_data" match the row names of the "Plot_SettingsFile"?
      if(nrow(Test_match) ==  0){
        stop("row.names Input_data need to match row.names Plot_SettingsFile.")
      } else(
        data <- Input_data
      )
    }
    Plot_SettingsFile <- Plot_SettingsFile
  }

  ## ------------ Check Input SettingsInfo ----------- ##
    #2. Plot_SettingsInfo
  if(Plot_SettingsInfo[["conditions"]] %in% Plot_SettingsInfo==TRUE){
    if(Plot_SettingsInfo[["conditions"]] %in% colnames(Plot_SettingsFile)== TRUE){
      Plot_SettingsFile<- Plot_SettingsFile%>%
        dplyr::rename("Conditions"= paste(Plot_SettingsInfo[["conditions"]]) )
    }else{# if true rename to Conditions
      stop("The ",Plot_SettingsInfo[["conditions"]], " column selected as conditions in Plot_SettingsInfo was not found in Plot_SettingsFile. Please check your input.")
    }
  }else{
    stop("You have to provide a Plot_SettingsInfo for conditions.")
  }

  if("superplot" %in% names(Plot_SettingsInfo)){
     if(Plot_SettingsInfo[["superplot"]] %in% colnames(Plot_SettingsFile)== TRUE){
       if(Plot_SettingsInfo[["superplot"]] %in% Plot_SettingsInfo==TRUE){
         Plot_SettingsFile<- Plot_SettingsFile%>%
           dplyr::rename("superplot"= paste(Plot_SettingsInfo[["superplot"]]) )
         }else{
           stop("The ",Plot_SettingsInfo[["superplot"]], " column selected for superplot in Plot_SettingsInfo was not found in Plot_SettingsFile. Please check your input.")
         }
     }
    }



  ## ------------ Check other plot parameters ----------- ##
  #3. Plot options
  if(Graph_Style %in% c("Box", "Bar", "Violin") == FALSE){
    stop("Graph_Style must be either Box, Bar ot Violin.")
  }

  #4. Comparison options
  if(is.null(Selected_Conditions) == FALSE){
    for (Condition in Selected_Conditions){
      if(Condition %in% Plot_SettingsFile$Conditions==FALSE){
        stop(paste0("Check Input. The Selected_Conditions ",Condition," were not found in the Conditions Column."))
      }
    }
  }

  if(is.null(Selected_Conditions) == TRUE & is.null(Selected_Comparisons)==FALSE){
    stop(paste0("The comparisons of the Selected_Comparisons parameter are selected from the Selected_Conditions vector. Comparisons were added while the Selected_Conditions vector is NULL. Please check your input."))
  }

  if(is.null(Selected_Comparisons)==FALSE){
    for (Comp in Selected_Comparisons){
      if(is.null(Selected_Conditions)==FALSE){
        if(Selected_Conditions[Comp[1]] %in% Plot_SettingsFile$Conditions ==FALSE){
        stop("Check Input. The Selected_Comparisons condition ",paste(Comp[1]), " is not found in the Conditions Column of the Plot_SettingsFile.")
          }
        if(Selected_Conditions[Comp[2]] %in% Plot_SettingsFile$Conditions ==FALSE){
        stop("Check Input. The Selected_Comparisons condition ",paste(Comp[2]), " is not found in the Conditions Column of the Plot_SettingsFile.")
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

  #6. Check color_palette:
  if(is.null(color_palette)){
    color_palette <- "grey"
  }else{
    if(is.null(Selected_Conditions)==TRUE){
      if(length(color_palette) !=length(unique(Plot_SettingsFile$Conditions)) ){
        stop(paste0("The color_palette colors used are not 1-to-1 with the groups in the ",Plot_SettingsInfo[["conditions"]] ," column. Please check your Input."))
      }
    }else{
      if(length(color_palette) !=length(Selected_Conditions) ){
        stop("The Selected_Conditions and the color_palette used are not 1-to-1. Please check your Input.")
      }
    }
  }

  #7. Check Selected_Comparisons & Selected_Conditions
  if(is.null(Selected_Conditions)){
    Number_Cond <- length(unique(tolower(Plot_SettingsFile[["Conditions"]])))
    if(Number_Cond<=2){
      MultipleComparison = FALSE
    }else{
      MultipleComparison = TRUE
    }
  }else if(length(Selected_Conditions)>2){
    MultipleComparison = TRUE
  }else if(length(Selected_Conditions)<=2){
    Number_Cond <- length(unique(tolower(Plot_SettingsFile[["Conditions"]])))
    if(Number_Cond<=2){
      MultipleComparison = FALSE
    }else{
      MultipleComparison = TRUE
    }
  }

  #8. Check Stat values:
  STAT_pval_options <- c("t.test", "wilcox.test", "aov", "kruskal.test")
  if(is.null(STAT_pval)==FALSE){
    if(STAT_pval %in% STAT_pval_options == FALSE & is.null(STAT_pval)==FALSE){
      stop("Check input. The selected STAT_pval option for Hypothesis testing is not valid. Please select NULL or one of the following: ",paste(STAT_pval_options,collapse = ", "),"." )
    }
    }

  if(is.null(STAT_pval)==FALSE){
    if(MultipleComparison == TRUE & (STAT_pval=="t.test" | STAT_pval=="wilcox.test")){
      warning("The selected STAT_pval option for Hypothesis testing,", STAT_pval, " is for one-versus-one comparison, but you have more than 2 conditions. Hence aov is performed.")
      STAT_pval <- "aov"
    }else if(MultipleComparison == FALSE & (STAT_pval=="aov" | STAT_pval=="kruskal.test")){
      warning("The selected STAT_pval option for Hypothesis testing,", STAT_pval, " is for multiple comparison, but you have only 2 conditions. Hence t.test is performed.")
      STAT_pval <- "t.test"
      }
    }




  if(is.null(STAT_pval)==TRUE & MultipleComparison == FALSE){
    STAT_pval <- "t.test"
  }

  if(is.null(STAT_pval)==TRUE & MultipleComparison == TRUE){
    STAT_pval <- "aov"
  }

  STAT_padj_options <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
  if(is.null(STAT_padj)==FALSE){
    if(STAT_padj %in% STAT_padj_options == FALSE){
      stop("Check input. The selected STAT_padj option for multiple Hypothesis testing correction is not valid. Please select NULL or one of the folowing: ",paste(STAT_padj_options,collapse = ", "),"." )
  }
  }

  if(is.null(STAT_padj)==TRUE){
    STAT_padj <- "fdr"
  }





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

  ## ------------ Create plots ----------- ##
  data <- Input_data
  Metabolite_Names <- colnames(data)
  if("superplot" %in% names(Plot_SettingsInfo)){
    data <- merge(Plot_SettingsFile[c("Conditions","superplot")] ,data, by=0)
    data <- column_to_rownames(data, "Row.names")
  }else{
    data <- merge(Plot_SettingsFile[c("Conditions")] ,data, by=0)
    data <- column_to_rownames(data, "Row.names")
  }


  # make a list for plotting all plots together
  PlotList <- list()#Empty list to store all the plots
  PlotList_adaptedGrid <- list()#Empty list to store all the plots

  for (i in Metabolite_Names){
    #Prepare the dfs:
    suppressWarnings(dataMeans <- data %>%
                       select(i, Conditions)
                     %>% group_by(Conditions)
                     %>% summarise_at(vars(i), list(mean = mean, sd = sd))
                     %>% as.data.frame())
    names(dataMeans)[2] <- "Intensity"

    if("superplot" %in% names(Plot_SettingsInfo)){
      suppressWarnings(plotdata <- data %>%
                       select(i,Conditions, superplot)
                     %>%  group_by(Conditions)
                     %>% as.data.frame() )
    }else{
      suppressWarnings(plotdata <- data %>%
                         select(i,Conditions)
                       %>%  group_by(Conditions)
                       %>% as.data.frame() )
    }
    names(plotdata)[1] <- c("Intensity")
    plotdata$Conditions <- factor(plotdata$Conditions)# Change conditions to factor

    # Take only selected conditions
    if(is.null(Selected_Conditions) == FALSE){
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
      Plot <- Plot+  geom_bar(stat = "summary", fun = "mean", fill = color_palette)+ stat_summary(fun.data=data_summary,
                                                                                            geom="errorbar", color="black", width=0.2)
    } else if (Graph_Style == "Violin"){
      Plot <- Plot+ geom_violin(fill = color_palette)+ stat_summary(fun.data=data_summary,
                                                              geom="errorbar", color="black", width=0.2)
    } else if (Graph_Style == "Box"){
      Plot <- Plot +  geom_boxplot(fill=color_palette,  width=0.5, position=position_dodge(width = 0.5))
    }

    # Add superplot
    if ("superplot" %in% names(Plot_SettingsInfo)){
      if(is.null(color_palette_dot)==FALSE){
        Plot <- Plot+ ggbeeswarm::geom_beeswarm(aes(x=Conditions,y=Intensity,color=as.factor(superplot)),size=3)+
          labs(color=Plot_SettingsInfo[["superplot"]], fill = Plot_SettingsInfo[["superplot"]])+
          scale_color_manual(values = color_palette_dot)
      }else{
        Plot <- Plot+ ggbeeswarm::geom_beeswarm(aes(x=Conditions,y=Intensity,color=as.factor(superplot)),size=3)+
          labs(color=Plot_SettingsInfo[["superplot"]], fill = Plot_SettingsInfo[["superplot"]])
      }
    }else{
      Plot <- Plot+ ggbeeswarm::geom_beeswarm(aes(x=Conditions,y=Intensity),size=2)
    }

    ####---- Add stats:
    if(STAT_pval=="t.test" | STAT_pval=="wilcox.test"){
      # One vs. One comparison: t-test
      if(is.null(Selected_Comparisons)==FALSE){
        Plot <- Plot+ ggpubr::stat_compare_means(comparisons = Selected_Comparisons,
                                                 label = "p.format", method = STAT_pval, hide.ns = TRUE,
                                               position = position_dodge(0.9), vjust = 0.25, show.legend = FALSE)
      }else{
        comparison <- unique(plotdata$Conditions)
        Plot <- Plot+ ggpubr::stat_compare_means(comparisons = comparison ,
                                                 label = "p.format", method = STAT_pval, hide.ns = TRUE,
                                                 position = position_dodge(0.9), vjust = 0.25, show.legend = FALSE)

      }
      Plot <- Plot +labs(caption = paste("p.val using pairwise ", STAT_pval))
      }else{
        #All-vs-All comparisons table:
        conditions <- Plot_SettingsFile$Conditions
        denominator <-unique(Plot_SettingsFile$Conditions)
        numerator <-unique(Plot_SettingsFile$Conditions)
        comparisons <- combn(unique(conditions), 2) %>% as.matrix()

        #Prepare Stat results using MetaProViz::DMA STAT helper functions
        if(STAT_pval=="aov"){
        STAT_C1vC2 <- MetaProViz:::AOV(Input_data=data.frame("Intensity" = plotdata[,-c(2:3)]),
                          conditions= plotdata[,c(2)],
                          Input_SettingsInfo=c(conditions="Conditions"),
                          STAT_padj=STAT_padj,
                          Log2FC_table=NULL,
                          all_vs_all=TRUE,
                          comparisons=comparisons)
        }else if(STAT_pval=="kruskal.test"){
          STAT_C1vC2 <-MetaProViz:::Kruskal(Input_data=data.frame("Intensity" = plotdata[,-c(2:3)]),
                                            conditions=plotdata[,c(2)],
                                            STAT_padj=STAT_padj,
                                            Log2FC_table=NULL,
                                            all_vs_all=all_vs_all,
                                            comparisons=comparisons)
        }

        #Prepare df to add stats to plot
        df <- data.frame(comparisons = names(STAT_C1vC2), stringsAsFactors = FALSE)%>%
          separate(comparisons, into=c("group1", "group2"), sep="_vs_", remove=FALSE)%>%
          unite(comparisons_rev, c("group2", "group1"), sep="_vs_", remove=FALSE)
        df$p.adj <- round(sapply(STAT_C1vC2, function(x) x$p.adj),5)
        df$y.position <-c(max(dataMeans$Intensity + 2*dataMeans$sd),
                          max(dataMeans$Intensity + 2*dataMeans$sd)+0.04* max(dataMeans$Intensity + 2*dataMeans$sd) ,
                          max(dataMeans$Intensity + 2*dataMeans$sd)+0.08* max(dataMeans$Intensity + 2*dataMeans$sd))

        # select stats based on comparison_table
        if(is.null(Selected_Comparisons)== FALSE){
          # Generate the comparisons
          df_select <- data.frame()
          for(comp in Selected_Comparisons){
            entry <- paste0(Selected_Conditions[comp[1]], "_vs_", Selected_Conditions[comp[2]])
            df_select <- rbind(df_select, data.frame(entry))
          }

          df_merge <- merge(df_select, df, by.x="entry", by.y="comparisons", all.x=TRUE)%>%
            column_to_rownames("entry")

          if(all(is.na(df_merge))==TRUE){#in case the reverse comparisons are needed
            df_merge <- merge(df_select, df, by.x="entry", by.y="comparisons_rev", all.x=TRUE)%>%
              column_to_rownames("entry")
          }
        }else{
          df_merge <- df[,-2]%>%
            column_to_rownames("comparisons")
          }


        # add stats to plot
        if(Graph_Style == "Bar"){
          Plot <- Plot +ggpubr::stat_pvalue_manual(df_merge, hide.ns = FALSE, size = 3, tip.length = 0.01, step.increase=0.05)#http://rpkgs.datanovia.com/ggpubr/reference/stat_pvalue_manual.html
        }else{
          Plot <- Plot +ggpubr::stat_pvalue_manual(df_merge, hide.ns = FALSE, size = 3, tip.length = 0.01, step.increase=0.01)#http://rpkgs.datanovia.com/ggpubr/reference/stat_pvalue_manual.html
        }
        Plot <- Plot +labs(caption = paste("p.adj using ", STAT_pval, "and", STAT_padj))
     }

    Plot <- Plot + Theme+ ggtitle(paste(i))
    Plot <- Plot + theme(legend.position = "right",plot.title = element_text(size=12, face = "bold"), axis.text.x = element_text(angle = 90, hjust = 1))+ xlab("")+ ylab("Normalized Intensity")

    ## Store the plot in the 'plots' list
    PlotList[[i]] <- Plot


    # Make plot into nice format:
    # MetaProViz:::
    Plot_Sized <-  MetaProViz:::plotGrob_Superplot(Input=Plot, Plot_SettingsInfo=Plot_SettingsInfo, Plot_SettingsFile=Plot_SettingsFile, MetaboliteName=i, Graph_Style=Graph_Style)
    Plot <- Plot_Sized[[3]]

    # First we want to convert the plot back into a ggplot object:
    Plot <- ggplot2::ggplot() +
      annotation_custom(Plot)
    Plot <-Plot + theme(panel.background = element_rect(fill = "transparent"))

    PlotList_adaptedGrid[[i]] <- Plot

    if(Individual_plots==TRUE){
      cleaned_i <- (gsub("/","_",i))#remove "/" cause this can not be safed in a PDF name
      cleaned_i <- (gsub(":","_",cleaned_i))
      cleaned_i <-(gsub("\\*","",cleaned_i))

      if (!is.null(Save_as_Plot)) {
        if(OutputPlotName ==""){
          ggsave(file=paste(Results_folder_plots_Barplots_folder_Individual, "/",cleaned_i, ".",Save_as_Plot, sep=""), plot=Plot, width=10, height=8)
        }else{
          ggsave(file=paste(Results_folder_plots_Barplots_folder_Individual, "/",OutputPlotName,"_",cleaned_i, ".",Save_as_Plot, sep=""), plot=Plot, width=10, height=8)
        }
      }
    }
  }

  if(Individual_plots==FALSE){
    if (!is.null(Save_as_Plot)){
      if(OutputPlotName ==""){
        pdf(file= paste(Results_folder_plots_Barplots_folder,"/",Graph_Style, "plots", OutputPlotName,".pdf", sep = ""), onefile = TRUE )
      }else{
        pdf(file= paste(Results_folder_plots_Barplots_folder,"/",Graph_Style, "plots","_", OutputPlotName,".pdf", sep = ""), onefile = TRUE )
      }
      dev.off()
    }
  }
  return(invisible(list("Plot"=PlotList,"Plot_Sized" = PlotList_adaptedGrid)))
}



#####################################################################
### ### ### Superplots helper function: Internal Function ### ### ###
#####################################################################

#' @param Input This is the ggplot object generated within the VizSuperplots function.
#' @param Plot_SettingsInfo Passed to VizSuperplots
#' @param Plot_SettingsFile Passed to VizSuperplots
#' @param MetaboliteName Passed to VizSuperplots
#' @param Graph_Style Passed to VizSuperplots
#'
#' @keywords PCA helper function
#' @noRd

plotGrob_Superplot <- function(Input, Plot_SettingsInfo, Plot_SettingsFile, MetaboliteName, Graph_Style){
  #------- Set the total heights and widths
  #we need ggplot_grob to edit the gtable of the ggplot object. Using this we can manipulate the gtable arguments directly.
  plottable <- ggplot2::ggplotGrob(Input) # Convert the plot to a gtable

  #-----widths (adapt for number of conditions)
  Number_Conditions <- Plot_SettingsFile%>%
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
  plottable$widths[c(6)] <- unit((Number_Conditions * 0.5)-1,"cm")#controls margins --> start Figure legend
  plottable$widths[c(7,8,10,11)] <- unit(0,"cm")#controls margins --> not needed
  plot_widths <- as.numeric(plottable$widths[5])+4

  plot(plottable)

  if("superplot" %in% Plot_SettingsInfo==TRUE){#legend will be present!
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

  if("superplot" %in% Plot_SettingsInfo==TRUE){#legend will be present!
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


