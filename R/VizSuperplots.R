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

#' Bar, Box or Violin plot in Superplot style visualization
#'
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
#' @param FolderPath \emph{Optional:} Path to the folder the results should be saved at. \strong{Default = NULL}
#'
#' @return List with two elements: Plot and Plot_Sized
#'
#' @examples
#' Intra <- ToyData("IntraCells_Raw")[,c(1:6)]
#' Res <- VizSuperplot(InputData=Intra[,-c(1:3)], SettingsFile_Sample=Intra[,c(1:3)], SettingsInfo = c(Conditions="Conditions", Superplot = NULL))
#'
#' @keywords Barplot, Boxplot, Violinplot, Superplot
#'
#' @importFrom ggplot2 ggplot theme geom_violin stat_summary geom_boxplot
#' @importFrom ggplot2 geom_bar labs scale_color_manual theme xlab ylab element_text
#' @importFrom ggpubr stat_pvalue_manual
#' @importFrom grid convertUnit
#' @importFrom dplyr rename select group_by summarise filter mutate n across
#' @importFrom tidyr separate unite
#' @importFrom magrittr %>% %<>%
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom ggbeeswarm geom_beeswarm
#' @importFrom logger log_trace log_info
#' @importFrom tidyselect all_of
#'
#' @export
#'
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
                         FolderPath = NULL){

  ## ------------ Create log file ----------- ##
  MetaProViz_Init()

  logger::log_info("VizSuperplot: Superplot visualization")

  ## ------------ Check Input files ----------- ##
  # HelperFunction `CheckInput`
  CheckInput(InputData=InputData,
                          SettingsFile_Sample=SettingsFile_Sample,
                          SettingsFile_Metab=NULL,
                          SettingsInfo=SettingsInfo,
                          SaveAs_Plot=SaveAs_Plot,
                          SaveAs_Table=NULL,
                          CoRe=FALSE,
                          PrintPlot= PrintPlot)

  # CheckInput` Specific
  if(is.null(SettingsInfo)==TRUE){
    message <- paste0("You must provide the column name for Conditions via SettingsInfo=c(Conditions=ColumnName) in order to plot the x-axis conditions.")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }

  if(PlotType %in% c("Box", "Bar", "Violin") == FALSE){
    message <- paste0("PlotType must be either Box, Bar or Violin.")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }

  if(is.null(PlotConditions) == FALSE){
    for (Condition in PlotConditions){
      if(Condition %in% SettingsFile_Sample[[SettingsInfo[["Conditions"]]]]==FALSE){
        message <- paste0("Check Input. The PlotConditions ",Condition," were not found in the Conditions Column.")
        logger::log_trace(paste("Error ", message, sep=""))
        stop(message)
      }
    }
  }

  if(is.null(StatComparisons)==FALSE){
    for (Comp in StatComparisons){
      if(is.null(PlotConditions)==FALSE){
        if(PlotConditions[Comp[1]] %in% SettingsFile_Sample[[SettingsInfo[["Conditions"]]]] ==FALSE){
          message <- paste0("Check Input. The StatComparisons condition ",Comp[1], " is not found in the Conditions Column of the SettingsFile_Sample.")
          logger::log_trace(paste("Error ", message, sep=""))
          stop(message)
        }
        if(PlotConditions[Comp[2]] %in% SettingsFile_Sample[[SettingsInfo[["Conditions"]]]] ==FALSE){
          message <- paste0("Check Input. The StatComparisons condition ",Comp[2], " is not found in the Conditions Column of the SettingsFile_Sample.")
          logger::log_trace(paste("Error ", message, sep=""))
          stop(message)
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
      message <- paste0("Check input. The selected StatPval option for Hypothesis testing,", StatPval, " is for multiple comparison, but you have only 2 conditions. Hence aov is performed.")
      logger::log_trace(paste("Warning ", message, sep=""))
      warning(message)
      StatPval <- "aov"
    }else if(MultipleComparison == FALSE & (StatPval=="aov" | StatPval=="kruskal.test")){
      message <- paste0("Check input. The selected StatPval option for Hypothesis testing,", StatPval, " is for multiple comparison, but you have only 2 conditions. Hence t.test is performed.")
      logger::log_trace(paste("Warning ", message, sep=""))
      warning(message)
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
      message <- paste0("Check input. The selected StatPadj option for multiple Hypothesis testing correction is not valid. Please select NULL or one of the folowing: ",paste(STAT_padj_options,collapse = ", "),"." )
      logger::log_trace(paste("Error ", message, sep=""))
      stop(message)
  }
  }

  if(is.null(StatPadj)==TRUE){
    StatPadj <- "fdr"
  }

  ## ------------ Create Results output folder ----------- ##
  if(is.null(SaveAs_Plot)==FALSE){
    Folder <- SavePath(FolderName=  paste(PlotType, "Plots", sep=""),
                                    FolderPath=FolderPath)
  }
  logger::log_info("VizSuperplot results saved at ", Folder)

  ###############################################################################################################################################################################################################
  ## ------------ Prepare Input ----------- ##
  SettingsFile_Sample<- SettingsFile_Sample%>%
    dplyr::rename("Conditions"= paste(SettingsInfo[["Conditions"]]) )

  if("Superplot" %in% names(SettingsInfo)){
    SettingsFile_Sample<- SettingsFile_Sample%>%
      dplyr::rename("Superplot"= paste(SettingsInfo[["Superplot"]]) )

    data <- merge(SettingsFile_Sample[c("Conditions","Superplot")] ,InputData, by=0)
    data <- tibble::column_to_rownames(data, "Row.names")
  }else{
    data <- merge(SettingsFile_Sample[c("Conditions")] ,InputData, by=0)
    data <- tibble::column_to_rownames(data, "Row.names")
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
    Theme <- ggplot2::theme_classic()
  }

  ## ------------ Create plots ----------- ##
  # make a list for plotting all plots together
  PlotList <- list()#Empty list to store all the plots
  PlotList_adaptedGrid <- list()#Empty list to store all the plots

  for (i in colnames(InputData)){
    #Prepare the dfs:
    suppressWarnings(
      dataMeans <-
        data %>%
        dplyr::select(i, Conditions) %>%
        dplyr::group_by(Conditions) %>%
        dplyr::summarise(
          dplyr::across(
            tidyselect::all_of(i),
            list(mean = mean, sd = sd)
          )
        ) %>%
        as.data.frame()
    )

    names(dataMeans)[2] <- "Intensity"

    if("Superplot" %in% names(SettingsInfo)){
      suppressWarnings(plotdata <- data %>%
                         dplyr::select(i,Conditions, Superplot)
                     %>%  dplyr::group_by(Conditions)
                     %>% as.data.frame() )
    }else{
      suppressWarnings(plotdata <- data %>%
                         dplyr::select(i,Conditions)
                       %>%  dplyr::group_by(Conditions)
                       %>% as.data.frame() )
    }
    names(plotdata)[1] <- c("Intensity")
    plotdata$Conditions <- factor(plotdata$Conditions)# Change conditions to factor

    # Take only selected conditions
    if(is.null(PlotConditions) == FALSE){
      dataMeans <- dataMeans %>% dplyr::filter(Conditions %in% PlotConditions)
      plotdata <- plotdata %>% dplyr::filter(Conditions %in% PlotConditions)
      plotdata$Conditions <- factor(plotdata$Conditions, levels = PlotConditions)
    }

    # Make the Plot
    Plot <- ggplot2::ggplot(plotdata, aes(x = Conditions, y = Intensity))

    # Add graph style and error bar
    data_summary <- function(x){#calculate error bar!
      m <- mean(x)
      ymin <- m-sd(x)
      ymax <- m+sd(x)
      return(c(y=m,ymin=ymin,ymax=ymax))
    }

    if (PlotType == "Bar"){
      Plot <- Plot+  ggplot2::geom_bar(stat = "summary", fun = "mean", fill = ColorPalette)+ ggplot2::stat_summary(fun.data=data_summary,
                                                                                            geom="errorbar", color="black", width=0.2)
    } else if (PlotType == "Violin"){
      Plot <- Plot+ ggplot2::geom_violin(fill = ColorPalette)+ ggplot2::stat_summary(fun.data=data_summary,
                                                              geom="errorbar", color="black", width=0.2)
    } else if (PlotType == "Box"){
      Plot <- Plot +  ggplot2::geom_boxplot(fill=ColorPalette,  width=0.5, position=position_dodge(width = 0.5))
    }

    # Add Superplot
    if ("Superplot" %in% names(SettingsInfo)){
      if(is.null(ColorPalette_Dot)==FALSE){
        Plot <- Plot+ ggbeeswarm::geom_beeswarm(aes(x=Conditions,y=Intensity,color=as.factor(Superplot)),size=3)+
          ggplot2::labs(color=SettingsInfo[["Superplot"]], fill = SettingsInfo[["Superplot"]])+
          ggplot2::scale_color_manual(values = ColorPalette_Dot)
      }else{
        Plot <- Plot+ ggbeeswarm::geom_beeswarm(aes(x=Conditions,y=Intensity,color=as.factor(Superplot)),size=3)+
          ggplot2::labs(color=SettingsInfo[["Superplot"]], fill = SettingsInfo[["Superplot"]])
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
      Plot <- Plot +ggplot2::labs(caption = paste("p.val using pairwise ", StatPval))
      }else{
        #All-vs-All comparisons table:
        conditions <- SettingsFile_Sample$Conditions
        denominator <-unique(SettingsFile_Sample$Conditions)
        numerator <-unique(SettingsFile_Sample$Conditions)
        comparisons <- combn(unique(conditions), 2) %>% as.matrix()

        #Prepare Stat results using DMA STAT helper functions
        if(StatPval=="aov"){
        STAT_C1vC2 <- AOV(InputData=data.frame("Intensity" = plotdata[,-c(2:3)]),
                          SettingsInfo=c(Conditions="Conditions", Numerator = unique(SettingsFile_Sample$Conditions), Denominator  = unique(SettingsFile_Sample$Conditions)),
                          SettingsFile_Sample= SettingsFile_Sample,
                          Log2FC_table=NULL)
        }else if(StatPval=="kruskal.test"){
          STAT_C1vC2 <-Kruskal(InputData=data.frame("Intensity" = plotdata[,-c(2:3)]),
                                            SettingsInfo=c(Conditions="Conditions", Numerator = unique(SettingsFile_Sample$Conditions), Denominator  = unique(SettingsFile_Sample$Conditions)),
                                            SettingsFile_Sample= SettingsFile_Sample,
                                            Log2FC_table=NULL)
        }

        #Prepare df to add stats to plot
        df <- data.frame(comparisons = names(STAT_C1vC2), stringsAsFactors = FALSE)%>%
          tidyr::separate(comparisons, into=c("group1", "group2"), sep="_vs_", remove=FALSE)%>%
          tidyr::unite(comparisons_rev, c("group2", "group1"), sep="_vs_", remove=FALSE)
        df$p.adj <- round(sapply(STAT_C1vC2, function(x) x$p.adj),5)

        # Add the 'res' column by repeating 'position' to match the number of rows
        position <- c(max(dataMeans$Intensity + 2*dataMeans$sd),
                      max(dataMeans$Intensity + 2*dataMeans$sd)+0.04* max(dataMeans$Intensity + 2*dataMeans$sd) ,
                      max(dataMeans$Intensity + 2*dataMeans$sd)+0.08* max(dataMeans$Intensity + 2*dataMeans$sd))

        df <- df %>%
          dplyr::mutate(y.position = rep(position, length.out = dplyr::n()))

        # select stats based on comparison_table
        if(is.null(StatComparisons)== FALSE){
          # Generate the comparisons
          df_select <- data.frame()
          for(comp in StatComparisons){
            entry <- paste0(PlotConditions[comp[1]], "_vs_", PlotConditions[comp[2]])
            df_select <- rbind(df_select, data.frame(entry))
          }

          df_merge <- merge(df_select, df, by.x="entry", by.y="comparisons", all.x=TRUE)%>%
            tibble::column_to_rownames("entry")

          if(all(is.na(df_merge))==TRUE){#in case the reverse comparisons are needed
            df_merge <- merge(df_select, df, by.x="entry", by.y="comparisons_rev", all.x=TRUE)%>%
              tibble::column_to_rownames("entry")
          }
        }else{
          df_merge <- df[,-2]%>%
            tibble::column_to_rownames("comparisons")
          }


        # add stats to plot
        if(PlotType == "Bar"){
          Plot <- Plot +ggpubr::stat_pvalue_manual(df_merge, hide.ns = FALSE, size = 3, tip.length = 0.01, step.increase=0.05)#http://rpkgs.datanovia.com/ggpubr/reference/stat_pvalue_manual.html
        }else{
          Plot <- Plot +ggpubr::stat_pvalue_manual(df_merge, hide.ns = FALSE, size = 3, tip.length = 0.01, step.increase=0.01)#http://rpkgs.datanovia.com/ggpubr/reference/stat_pvalue_manual.html
        }
        Plot <- Plot +ggplot2::labs(caption = paste("p.adj using ", StatPval, "and", StatPadj))
     }

    Plot <- Plot + Theme+ ggplot2::labs(title = PlotName,
                                subtitle = i)# ggtitle(paste(i))
    Plot <- Plot + ggplot2::theme(legend.position = "right",plot.title = ggplot2::element_text(size=12, face = "bold"), axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))+ ggplot2::xlab(xlab)+ ggplot2::ylab(ylab)

    ## Store the plot in the 'plots' list
    PlotList[[i]] <- Plot

    # Make plot into nice format:
    Plot_Sized <-  plotGrob_Superplot(InputPlot=Plot, SettingsInfo=SettingsInfo, SettingsFile_Sample=SettingsFile_Sample,  PlotName = PlotName, Subtitle = i, PlotType=PlotType)
    PlotHeight <- grid::convertUnit(Plot_Sized$height, 'cm', valueOnly = TRUE)
    PlotWidth <- grid::convertUnit(Plot_Sized$width, 'cm', valueOnly = TRUE)
    Plot_Sized %<>%
      {ggplot2::ggplot() + annotation_custom(.)} %>%
      add(theme(panel.background = element_rect(fill = "transparent")))

   ####################################################################################################################################
    ## --------------- save -----------------##
    cleaned_i <- gsub("[[:space:],/\\\\:*?\"<>|]", "-", i)#removes empty spaces and replaces /,\ with -
    PlotList_adaptedGrid[[cleaned_i]] <- Plot_Sized

    SaveList <- list()
    SaveList[[cleaned_i]] <- Plot_Sized
    #----- Save
    suppressMessages(suppressWarnings(
      SaveRes(InputList_DF=NULL,
                           InputList_Plot= SaveList,
                           SaveAs_Table=NULL,
                           SaveAs_Plot=SaveAs_Plot,
                           FolderPath= Folder,
                           FileName= paste(PlotType, "Plots_",PlotName, sep=""),
                           CoRe=FALSE,
                           PrintPlot=PrintPlot,
                           PlotHeight=PlotHeight,
                           PlotWidth=PlotWidth,
                           PlotUnit="cm")))
  }
  return(invisible(list("Plot"=PlotList,"Plot_Sized" = PlotList_adaptedGrid)))
}
