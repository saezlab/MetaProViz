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


#' @param InputData DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Includes experimental design and outlier column.
#' @param SettingsFile_Sample DF which contains information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames. Column "Conditions" with information about the sample conditions (e.g. "N" and "T" or "Normal" and "Tumor"), can be used for feature filtering and colour coding in the PCA. Column "AnalyticalReplicate" including numerical values, defines technical repetitions of measurements, which will be summarised. Column "BiologicalReplicates" including numerical values. Please use the following names: "Conditions", "Biological_Replicates", "Analytical_Replicates".
#' @param SettingsInfo Named vector including at least information on the conditions column: c(Conditions="ColumnName_SettingsFile_Sample"). Additionally Superplots can be made by adding Superplot ="ColumnName_SettingsFile_Sample", which are usually biological replicates or patient IDs. \strong{Default = c(Conditions="Conditions", Superplot = NULL)}
#' @param PlotType String with the information of the Graph style. Available options are Bar. Box and Violin  \strong{Default = Box}
#' @param PlotName \emph{Optional: } String which is added to the output files of the plot.
#' @param PlotConditions Vector with names of selected Conditions for the plot. Can also be used to order the Conditions in the way they should be displayed on the x-axis of the plot. \strong{Default = NULL}
#' @param StatComparisons List of numeric vectors containing Condition pairs to compare based on the order of the PlotConditions vector. \strong{Default = NULL}
#' @param StatPval \emph{Optional: } String which contains an abbreviation of the selected test to calculate p.value. For one-vs-one comparisons choose t.test or wilcox.test , for one-vs-all or all-vs-all comparison choose aov (=anova) or kruskal.test \strong{Default = NULL}
#' @param StatPadj \emph{Optional: } String which contains an abbreviation of the selected p.adjusted test for p.value correction for multiple Hypothesis testing. Search: ?p.adjust for more methods:"BH", "fdr", "bonferroni", "holm", etc.\strong{Default = NULL}
#' @param xlab \emph{Optional: } String to replace x-axis label in plot. \strong{Default = NULL}
#' @param ylab \emph{Optional: } String to replace y-axis label in plot. \strong{Default = NULL}
#' @param Theme \emph{Optional: } Selection of theme for plot, e.g. theme_grey(). You can check for complete themes here: https://ggplot2.tidyverse.org/reference/ggtheme.html. \strong{Default = NULL}
#' @param ColorPalette \emph{Optional: } Provide customized ColorPalette in vector format. \strong{Default = NULL}
#' @param ColorPalette_Dot \emph{Optional: } Provide customized ColorPalette in vector format. \strong{Default = NULL}
#' @param SaveAs_Plot \emph{Optional: } Select the file type of output plots. Options are svg, pdf, png or NULL. \strong{Default = svg}
#' @param PrintPlot \emph{Optional: } TRUE or FALSE, if TRUE plots are saved as an overview of the results. \strong{Default = TRUE}
#' @param FolderPath \emph{Optional:} Path to the folder the results should be saved at. \strong(Default = NULL)
#'
#' @keywords Barplot, Boxplot, Violinplot, Superplot
#' @export

VizSuperplot <- function(InputData,
                         SettingsFile_Sample,
                         SettingsInfo = c(Conditions="Conditions", Superplot = NULL),
                         PlotType = "Box", # Bar, Box, Violin
                         PlotName = "",
                         PlotConditions = NULL,
                         StatComparisons = NULL,
                         StatPval =NULL,
                         StatPadj=NULL,
                         xlab= NULL,
                         ylab= NULL,
                         Theme = NULL,
                         ColorPalette = NULL,
                         ColorPalette_Dot =NULL,
                         SaveAs_Plot = "svg",
                         PrintPlot=TRUE,
                         FolderPath = NULL
){

  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse", "ggplot2", "ggpubr", "ggbeeswarm")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  suppressMessages(library(tidyverse))

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
  if(is.null(SettingsInfo)==TRUE){
    stop("You must provide the column name for Conditions via SettingsInfo=c(Conditions=ColumnName) in order to plot the x-axis conditions.")
  }

  if(PlotType %in% c("Box", "Bar", "Violin") == FALSE){
    stop("PlotType must be either Box, Bar or Violin.")
  }

  if(is.null(PlotConditions) == FALSE){
    for (Condition in PlotConditions){
      if(Condition %in% SettingsFile_Sample[[SettingsInfo[["Conditions"]]]]==FALSE){
        stop(paste0("Check Input. The PlotConditions ",Condition," were not found in the Conditions Column."))
      }
    }
  }

  if(is.null(StatComparisons)==FALSE){
    for (Comp in StatComparisons){
      if(is.null(PlotConditions)==FALSE){
        if(PlotConditions[Comp[1]] %in% SettingsFile_Sample[[SettingsInfo[["Conditions"]]]] ==FALSE){
          stop("Check Input. The StatComparisons condition ",paste(Comp[1]), " is not found in the Conditions Column of the SettingsFile_Sample.")
        }
        if(PlotConditions[Comp[2]] %in% SettingsFile_Sample[[SettingsInfo[["Conditions"]]]] ==FALSE){
          stop("Check Input. The StatComparisons condition ",paste(Comp[2]), " is not found in the Conditions Column of the SettingsFile_Sample.")
        }
      }
    }
  }

  if(is.null(ColorPalette)){
    ColorPalette <- "grey"
  }


  ## ------------ Check Input SettingsInfo ----------- ##
  #7. Check StatComparisons & PlotConditions
  if(is.null(PlotConditions)){
    Number_Cond <- length(unique(tolower(SettingsFile_Sample[["Conditions"]])))
    if(Number_Cond<=2){
      MultipleComparison = FALSE
    }else{
      MultipleComparison = TRUE
    }
  }else if(length(PlotConditions)>2){
    MultipleComparison = TRUE
  }else if(length(PlotConditions)<=2){
    Number_Cond <- length(unique(tolower(SettingsFile_Sample[["Conditions"]])))
    if(Number_Cond<=2){
      MultipleComparison = FALSE
    }else{
      MultipleComparison = TRUE
    }
  }

  if(is.null(StatPval)==FALSE){
    if(MultipleComparison == TRUE & (StatPval=="t.test" | StatPval=="wilcox.test")){
      warning("The selected StatPval option for Hypothesis testing,", StatPval, " is for one-versus-one comparison, but you have more than 2 conditions. Hence aov is performed.")
      StatPval <- "aov"
    }else if(MultipleComparison == FALSE & (StatPval=="aov" | StatPval=="kruskal.test")){
      warning("The selected StatPval option for Hypothesis testing,", StatPval, " is for multiple comparison, but you have only 2 conditions. Hence t.test is performed.")
      StatPval <- "t.test"
      }
    }

  if(is.null(StatPval)==TRUE & MultipleComparison == FALSE){
    StatPval <- "t.test"
  }

  if(is.null(StatPval)==TRUE & MultipleComparison == TRUE){
    StatPval <- "aov"
  }

  STAT_padj_options <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
  if(is.null(StatPadj)==FALSE){
    if(StatPadj %in% STAT_padj_options == FALSE){
      stop("Check input. The selected StatPadj option for multiple Hypothesis testing correction is not valid. Please select NULL or one of the folowing: ",paste(STAT_padj_options,collapse = ", "),"." )
  }
  }

  if(is.null(StatPadj)==TRUE){
    StatPadj <- "fdr"
  }



  ## ------------ Create Results output folder ----------- ##
  if(is.null(SaveAs_Plot)==FALSE){
    Folder <- MetaProViz:::SavePath(FolderName=  paste(PlotType, "Plots", sep=""),
                                    FolderPath=FolderPath)
  }


  ###############################################################################################################################################################################################################
  ## ------------ Prepare Input ----------- ##
  SettingsFile_Sample<- SettingsFile_Sample%>%
    dplyr::rename("Conditions"= paste(SettingsInfo[["Conditions"]]) )

  if("Superplot" %in% names(SettingsInfo)){
    SettingsFile_Sample<- SettingsFile_Sample%>%
      dplyr::rename("Superplot"= paste(SettingsInfo[["Superplot"]]) )

    data <- merge(SettingsFile_Sample[c("Conditions","Superplot")] ,InputData, by=0)
    data <- column_to_rownames(data, "Row.names")
  }else{
    data <- merge(SettingsFile_Sample[c("Conditions")] ,InputData, by=0)
    data <- column_to_rownames(data, "Row.names")
  }

  # Rename the x and y lab if the information has been passed:
  if(is.null(xlab)==TRUE){#use column name of x provided by user
    xlab <- bquote(.(as.symbol(SettingsInfo[["Conditions"]])))
  }else if(is.null(xlab)==FALSE){
    xlab <- bquote(.(as.symbol(xlab)))
  }

  if(is.null(ylab)==TRUE){#use column name of x provided by user
    ylab <- bquote(.(as.symbol("Intensity")))
  }else if(is.null(ylab)==FALSE){
    ylab <- bquote(.(as.symbol(ylab)))
  }

  #Set the theme:
  if(is.null(Theme)==TRUE){
    Theme <- theme_classic()
  }

  ## ------------ Create plots ----------- ##
  # make a list for plotting all plots together
  PlotList <- list()#Empty list to store all the plots
  PlotList_adaptedGrid <- list()#Empty list to store all the plots

  for (i in colnames(InputData)){
    #Prepare the dfs:
    suppressWarnings(dataMeans <- data %>%
                       select(i, Conditions)
                     %>% group_by(Conditions)
                     %>% summarise_at(vars(i), list(mean = mean, sd = sd))
                     %>% as.data.frame())
    names(dataMeans)[2] <- "Intensity"

    if("Superplot" %in% names(SettingsInfo)){
      suppressWarnings(plotdata <- data %>%
                       select(i,Conditions, Superplot)
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
    if(is.null(PlotConditions) == FALSE){
      dataMeans <- dataMeans %>% filter(Conditions %in% PlotConditions)
      plotdata <- plotdata %>% filter(Conditions %in% PlotConditions)
      plotdata$Conditions <- factor(plotdata$Conditions, levels = PlotConditions)
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

    if (PlotType == "Bar"){
      Plot <- Plot+  geom_bar(stat = "summary", fun = "mean", fill = ColorPalette)+ stat_summary(fun.data=data_summary,
                                                                                            geom="errorbar", color="black", width=0.2)
    } else if (PlotType == "Violin"){
      Plot <- Plot+ geom_violin(fill = ColorPalette)+ stat_summary(fun.data=data_summary,
                                                              geom="errorbar", color="black", width=0.2)
    } else if (PlotType == "Box"){
      Plot <- Plot +  geom_boxplot(fill=ColorPalette,  width=0.5, position=position_dodge(width = 0.5))
    }

    # Add Superplot
    if ("Superplot" %in% names(SettingsInfo)){
      if(is.null(ColorPalette_Dot)==FALSE){
        Plot <- Plot+ ggbeeswarm::geom_beeswarm(aes(x=Conditions,y=Intensity,color=as.factor(Superplot)),size=3)+
          labs(color=SettingsInfo[["Superplot"]], fill = SettingsInfo[["Superplot"]])+
          scale_color_manual(values = ColorPalette_Dot)
      }else{
        Plot <- Plot+ ggbeeswarm::geom_beeswarm(aes(x=Conditions,y=Intensity,color=as.factor(Superplot)),size=3)+
          labs(color=SettingsInfo[["Superplot"]], fill = SettingsInfo[["Superplot"]])
      }
    }else{
      Plot <- Plot+ ggbeeswarm::geom_beeswarm(aes(x=Conditions,y=Intensity),size=2)
    }

    ####---- Add stats:
    if(StatPval=="t.test" | StatPval=="wilcox.test"){
      # One vs. One comparison: t-test
      if(is.null(StatComparisons)==FALSE){
        Plot <- Plot+ ggpubr::stat_compare_means(comparisons = StatComparisons,
                                                 label = "p.format", method = StatPval, hide.ns = TRUE,
                                               position = position_dodge(0.9), vjust = 0.25, show.legend = FALSE)
      }else{
        comparison <- unique(plotdata$Conditions)
        Plot <- Plot+ ggpubr::stat_compare_means(comparisons = comparison ,
                                                 label = "p.format", method = StatPval, hide.ns = TRUE,
                                                 position = position_dodge(0.9), vjust = 0.25, show.legend = FALSE)

      }
      Plot <- Plot +labs(caption = paste("p.val using pairwise ", StatPval))
      }else{
        #All-vs-All comparisons table:
        conditions <- SettingsFile_Sample$Conditions
        denominator <-unique(SettingsFile_Sample$Conditions)
        numerator <-unique(SettingsFile_Sample$Conditions)
        comparisons <- combn(unique(conditions), 2) %>% as.matrix()

        #Prepare Stat results using MetaProViz::DMA STAT helper functions
        if(StatPval=="aov"){
        STAT_C1vC2 <- MetaProViz:::AOV(InputData=data.frame("Intensity" = plotdata[,-c(2:3)]),
                          SettingsInfo=c(Conditions="Conditions", Numerator = unique(SettingsFile_Sample$Conditions), Denominator  = unique(SettingsFile_Sample$Conditions)),
                          SettingsFile_Sample= SettingsFile_Sample,
                          Log2FC_table=NULL)
        }else if(StatPval=="kruskal.test"){
          STAT_C1vC2 <-MetaProViz:::Kruskal(InputData=data.frame("Intensity" = plotdata[,-c(2:3)]),
                                            SettingsInfo=c(Conditions="Conditions", Numerator = unique(SettingsFile_Sample$Conditions), Denominator  = unique(SettingsFile_Sample$Conditions)),,
                                            SettingsFile_Sample= SettingsFile_Sample,
                                            Log2FC_table=NULL)
        }

        #Prepare df to add stats to plot
        df <- data.frame(comparisons = names(STAT_C1vC2), stringsAsFactors = FALSE)%>%
          separate(comparisons, into=c("group1", "group2"), sep="_vs_", remove=FALSE)%>%
          unite(comparisons_rev, c("group2", "group1"), sep="_vs_", remove=FALSE)
        df$p.adj <- round(sapply(STAT_C1vC2, function(x) x$p.adj),5)

        # Add the 'res' column by repeating 'position' to match the number of rows
        position <- c(max(dataMeans$Intensity + 2*dataMeans$sd),
                      max(dataMeans$Intensity + 2*dataMeans$sd)+0.04* max(dataMeans$Intensity + 2*dataMeans$sd) ,
                      max(dataMeans$Intensity + 2*dataMeans$sd)+0.08* max(dataMeans$Intensity + 2*dataMeans$sd))

        df <- df %>%
          mutate(y.position = rep(position, length.out = n()))

        # select stats based on comparison_table
        if(is.null(StatComparisons)== FALSE){
          # Generate the comparisons
          df_select <- data.frame()
          for(comp in StatComparisons){
            entry <- paste0(PlotConditions[comp[1]], "_vs_", PlotConditions[comp[2]])
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
        if(PlotType == "Bar"){
          Plot <- Plot +ggpubr::stat_pvalue_manual(df_merge, hide.ns = FALSE, size = 3, tip.length = 0.01, step.increase=0.05)#http://rpkgs.datanovia.com/ggpubr/reference/stat_pvalue_manual.html
        }else{
          Plot <- Plot +ggpubr::stat_pvalue_manual(df_merge, hide.ns = FALSE, size = 3, tip.length = 0.01, step.increase=0.01)#http://rpkgs.datanovia.com/ggpubr/reference/stat_pvalue_manual.html
        }
        Plot <- Plot +labs(caption = paste("p.adj using ", StatPval, "and", StatPadj))
     }

    Plot <- Plot + Theme+ labs(title = PlotName,
                                subtitle = i)# ggtitle(paste(i))
    Plot <- Plot + theme(legend.position = "right",plot.title = element_text(size=12, face = "bold"), axis.text.x = element_text(angle = 90, hjust = 1))+ xlab(xlab)+ ylab(ylab)

    ## Store the plot in the 'plots' list
    PlotList[[i]] <- Plot


    # Make plot into nice format:
    # MetaProViz:::
    Plot_Sized <-  MetaProViz:::plotGrob_Superplot(Input=Plot, SettingsInfo=SettingsInfo, SettingsFile_Sample=SettingsFile_Sample, MetaboliteName=i, PlotName=PlotName, PlotType=PlotType)
    Plot <- Plot_Sized[[3]]

    # First we want to convert the plot back into a ggplot object:
    Plot <- ggplot2::ggplot() +
      annotation_custom(Plot)
    Plot <-Plot + theme(panel.background = element_rect(fill = "transparent"))

    ####################################################################################################################################
    ## --------------- save -----------------##
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
                           FileName= paste(PlotType, "Plots_",PlotName, sep=""),
                           CoRe=FALSE,
                           PrintPlot=PrintPlot,
                           PlotHeight=Plot_Sized[[1]],
                           PlotWidth=Plot_Sized[[2]],
                           PlotUnit="cm")))
  }
  return(invisible(list("Plot"=PlotList,"Plot_Sized" = PlotList_adaptedGrid)))
}



#####################################################################
### ### ### Superplots helper function: Internal Function ### ### ###
#####################################################################

#' @param Input This is the ggplot object generated within the VizSuperplots function.
#' @param SettingsInfo Passed to VizSuperplots
#' @param SettingsFile_Sample Passed to VizSuperplots
#' @param MetaboliteName Passed to VizSuperplots
#' @param PlotName Passed to VizSuperplots
#' @param PlotType Passed to VizSuperplots
#'
#' @keywords PCA helper function
#' @noRd

plotGrob_Superplot <- function(Input,
                               SettingsInfo,
                               SettingsFile_Sample,
                               MetaboliteName,
                               PlotName,
                               PlotType){
  #------- Set the total heights and widths
  #we need ggplot_grob to edit the gtable of the ggplot object. Using this we can manipulate the gtable arguments directly.
  plottable <- ggplot2::ggplotGrob(Input) # Convert the plot to a gtable

  #-----widths (adapt for number of conditions)
  Number_Conditions <- SettingsFile_Sample%>%
    dplyr::distinct(Conditions) %>%
    nrow()

  if(PlotType == "Bar"){
    plottable$widths[5] <- unit(Number_Conditions * 0.5, "cm")#controls x-axis
  }else{
   plottable$widths[5] <- unit(Number_Conditions * 1, "cm")#controls x-axis
  }

  plottable$widths[c(1)] <- unit(0.5,"cm")#controls margins --> y-axis label is there
  plottable$widths[c(4)] <- unit(2,"cm")#controls margins --> y-axis label is there
  plottable$widths[c(2,3)] <- unit(0,"cm")#controls margins --> not needed
  plottable$widths[c(6)] <- unit((Number_Conditions * 0.5)-1,"cm")#controls margins --> start Figure legend
  plottable$widths[c(7,8,10)] <- unit(0,"cm")#controls margins --> not needed
  plottable$widths[c(11)] <- unit(0,"cm")#controls margins --> width
  plot_widths <- as.numeric(plottable$widths[5])+4
  #plot(plottable)

  if("Superplot" %in% names(SettingsInfo)==TRUE){#legend will be present!
    Value <- round(as.numeric(plottable$widths[9]),1) #plottable$widths[9] is a <unit/unit_v2> object and we can extract the extract the numeric part
    plot_widths <- plot_widths+Value
  }else{
    plottable$widths[c(9)] <- unit(0,"cm")
    plot_widths <- plot_widths
  }

  character_count_M <- nchar(MetaboliteName)#Check how much width is needed for the figure title/subtitle
  character_count_T <- nchar(PlotName)
  if (character_count_T >= character_count_M) {
    character_count <- character_count_T
  } else {
    character_count <- character_count_M
  }

  Titles_width <- (character_count*0.5)+0.8
  if(Titles_width>plot_widths){#If the title needs more space than the plot offers:
      plottable$widths[11] <- unit(character_count*0.125,"cm")#controls margins --> start Figure legend
      plot_widths <- plot_widths+character_count*0.125
      }

  #-----heigths
  plottable$heights[7] <- unit(8, "cm")#controls x-axis
  plottable$heights[c(8)] <- unit(1,"cm")#controls margins --> x-axis label
  plottable$heights[c(10)] <- unit(0.5,"cm")#controls margins --> Figure caption
  plottable$heights[c(6,9,11,12)] <- unit(0,"cm")#controls margins --> not needed (maybe 9 is needed in cases)
  plottable$heights[c(4)] <- unit(1,"cm")#controls margins --> Some space above the plot
  plottable$heights[c(1)] <- unit(1,"cm")#controls margins --> not needed
  plottable$heights[c(1,2,3,5)] <- unit(0,"cm")#controls margins --> not needed

  if("Superplot" %in% names(SettingsInfo)==TRUE){#legend will be present!
    #------- Legend heights
    Legend <- ggpubr::get_legend(Input) # Extract legend to adjust separately
    Legend_heights <- (round(as.numeric(Legend$heights[3]),1))+(round(as.numeric(Legend$heights[5]),1))
    if(Legend_heights>12){#If the legend requires more heights than the Plot
      Add <- (Legend_heights-12)/2
      plottable$heights[1] <- unit(Add,"cm")#controls margins --> Can be increased if Figure legend needs more space on the top
      plottable$heights[12] <- unit(Add,"cm")#controls margins --> Can be increased if Figure legend needs more space on the bottom
      plot_heights <- Legend_heights
    }else{
      plot_heights <- 12
    }
  }else{
    plot_heights <- 12
  }

  character_counts <- sapply(as.character(unique(SettingsFile_Sample$Conditions)), nchar)
  character_counts <- max(character_counts)# Get the longest character count
  if(character_counts>2){#If the title needs more space than the plot offers:
    plottable$heights[c(8)] <- unit(character_counts*0.21,"cm")
    plot_heights <- (plot_heights+character_counts*0.5)-1
  }

  #plot_param <-c(plot_heights=plot_heights, plot_widths=plot_widths)
  Output<- list(plot_heights, plot_widths, plottable)
}


