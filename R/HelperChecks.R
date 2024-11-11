## ---------------------------
##
## Script name: HelperFunctions
##
## Purpose of script: General helper functions to check function input and save results
##
## Author: Christina Schmidt
##
## Date Created: 2023-06-14
##
## Copyright (c) Christina Schmidt
## Email:
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------


################################################################################################
### ### ### Helper function: Internal Function to check function input ### ### ###
################################################################################################

#' Check input general parameters
#'
#' @param InputData Passed to MetaProViz functions. Usually DF which contains unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. But can also be differential expression results or other InputData.
#' @param InputData_Num  \emph{Optional: } If InputData must be numeric \strong{Default = TRUE}
#' @param SettingsFile_Sample \emph{Optional: } DF which contains information about the samples, which will be combined with the input data based on the unique sample identifiers used as rownames. If not avaliable can be set to NULL. \strong{Default = NULL}
#' @param SettingsFile_Metab \emph{Optional: } DF which contains information about the features. If not avaliable can be set to NULL. \strong{Default = NULL}
#' @param SettingsInfo \emph{Optional: } Passed to MetaProViz functions. Usually named vector containing the information about the names of the experimental parameters in SettingsFile_Sample, SettingsFile_Metab or InputData. \strong{Default = NULL}
#' @param SaveAs_Plot \emph{Optional: } Select the file type of output plots. Options are svg, png, pdf. If set to NULL, plots are not saved.\strong{Default = NULL}
#' @param SaveAs_Table \emph{Optional: } Select the file type of output table. Options are "csv", "xlsx", "txt". If set to NULL, plots are not saved. \strong{Default = NULL}
#' @param CoRe \emph{Optional: } If TRUE, a consumption-release experiment has been performed. If not avaliable can be set to NULL. \strong{Default = FALSE}
#' @param PrintPlot \emph{Optional: } If TRUE prints an overview of resulting plots. If not avaliable can be set to NULL. \strong{Default = FALSE}
#' @param Theme \emph{Optional: } Can be set for VizX functions. If not avaliable can be set to NULL. \strong{Default = NULL}
#' @param PlotSettings \emph{Optional: } Needs to be set for VizX functions. Options are "Sample", "Feature", Both". This refers to SettingsInfo color, shape, individual as for some plots we have both feature and sample settings. \strong{Default = NULL}
#'
#' @return returns warnings and errors if input is not correct
#'
#' @keywords Input check
#'
#' @importFrom logger log_trace
#'
#' @noRd
#'
CheckInput <- function(InputData,
                       InputData_Num=TRUE,
                       SettingsFile_Sample=NULL,
                       SettingsFile_Metab=NULL,
                       SettingsInfo=NULL,
                       SaveAs_Plot=NULL,
                       SaveAs_Table=NULL,
                       CoRe=FALSE,
                       PrintPlot=FALSE,
                       Theme=NULL,
                       PlotSettings=NULL){
  ############## Parameters valid for multiple MetaProViz functions

  #-------------InputData
  if(is.data.frame(InputData)==FALSE){
    message <- paste0("InputData should be a data.frame. It's currently a ", paste(class(InputData)), ".", sep = "")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }
  if(any(duplicated(row.names(InputData)))==TRUE){
    message <- paste0("Duplicated row.names of InputData, whilst row.names must be unique")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }

  if(InputData_Num==TRUE){
     Test_num <- apply(InputData, 2, function(x) is.numeric(x))
     if((any(Test_num) ==  FALSE) ==  TRUE){
       message <- paste0("InputData needs to be of class numeric")
       logger::log_trace(paste("Error ", message, sep=""))
       stop(message)
       }
  }

  if(sum(duplicated(colnames(InputData))) > 0){
    message <- paste0("InputData contained duplicates column names, whilst col.names must be unique.")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }

  #-------------SettingsFile
  if(is.null(SettingsFile_Sample)==FALSE){
    Test_match <- merge(SettingsFile_Sample, InputData, by = "row.names", all =  FALSE)
    if(nrow(Test_match) ==  0){
      message <- paste0("row.names InputData need to match row.names SettingsFile_Sample.")
      logger::log_trace(paste("Error ", message, sep=""))
      stop(message)
      }
  }

  if(is.null(SettingsFile_Metab)==FALSE){
    Test_match <- merge(SettingsFile_Metab, as.data.frame(t(InputData)), by = "row.names", all =  FALSE)
    if(nrow(Test_match) ==  0){
      stop("col.names InputData need to match row.names SettingsFile_Metab.")
    }
  }

  #-------------SettingsInfo
  if(is.vector(SettingsInfo)==FALSE & is.null(SettingsInfo)==FALSE){
    message <- paste0("SettingsInfo should be NULL or a vector. It's currently a ", paste(class(SettingsInfo), ".", sep = ""))
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }

  if(is.null(SettingsInfo)==FALSE){
    #Conditions
    if("Conditions" %in% names(SettingsInfo)){
      if(SettingsInfo[["Conditions"]] %in% colnames(SettingsFile_Sample)== FALSE){
        message <- paste0("The ", SettingsInfo[["Conditions"]], " column selected as Conditions in SettingsInfo was not found in SettingsFile. Please check your input.")
        logger::log_trace(paste("Error ", message, sep=""))
        stop(message)
      }
    }

    #Biological replicates
    if("Biological_Replicates" %in% names(SettingsInfo)){
      if(SettingsInfo[["Biological_Replicates"]] %in% colnames(SettingsFile_Sample)== FALSE){
        message <- paste0("The ", SettingsInfo[["Biological_Replicates"]], " column selected as Biological_Replicates in SettingsInfo was not found in SettingsFile_Sample. Please check your input.")
        logger::log_trace(paste("Error ", message, sep=""))
        stop(message)
      }
    }

    #Numerator
    if("Numerator" %in% names(SettingsInfo)==TRUE){
      if(SettingsInfo[["Numerator"]] %in% SettingsFile_Sample[[SettingsInfo[["Conditions"]]]]==FALSE){
        message <- paste0("The ",SettingsInfo[["Numerator"]], " column selected as numerator in SettingsInfo was not found in SettingsFile_Sample. Please check your input.")
        logger::log_trace(paste("Error ", message, sep=""))
        stop(message)
      }
    }

   #Denominator
    if("Denominator" %in% names(SettingsInfo)==TRUE){
      if(SettingsInfo[["Denominator"]] %in% SettingsFile_Sample[[SettingsInfo[["Conditions"]]]]==FALSE){
        message <- paste0("The ",SettingsInfo[["Denominator"]], " column selected as denominator in SettingsInfo was not found in SettingsFile_Sample. Please check your input.")
        logger::log_trace(paste("Error ", message, sep=""))
        stop(message)
      }
    }

    #Denominator & Numerator
    if("Denominator" %in% names(SettingsInfo)==FALSE  & "Numerator" %in% names(SettingsInfo) ==TRUE){
      message <- paste0("Check input. The selected denominator option is empty while ",paste(SettingsInfo[["Numerator"]])," has been selected as a numerator. Please add a denominator for 1-vs-1 comparison or remove the numerator for all-vs-all comparison.")
      logger::log_trace(paste("Error ", message, sep=""))
      stop(message)
    }

    #Superplot
    if("Superplot" %in% names(SettingsInfo)){
      if(SettingsInfo[["Superplot"]] %in% colnames(SettingsFile_Sample)== FALSE){
        message <- paste0("The ",SettingsInfo[["Superplot"]], " column selected as Superplot column in SettingsInfo was not found in SettingsFile_Sample. Please check your input.")
        logger::log_trace(paste("Error ", message, sep=""))
        stop(message)
      }
    }

    if(is.null(PlotSettings)==FALSE){
      if(PlotSettings== "Sample"){
        #Plot colour
        if("color" %in% names(SettingsInfo)){
          if(SettingsInfo[["color"]] %in% colnames(SettingsFile_Sample)== FALSE){
            message <- paste0("The ",SettingsInfo[["color"]], " column selected as color in SettingsInfo was not found in SettingsFile_Sample. Please check your input.")
            logger::log_trace(paste("Error ", message, sep=""))
            stop(message)
          }
        }

        #Plot shape
        if("shape" %in% names(SettingsInfo)){
          if(SettingsInfo[["shape"]] %in% colnames(SettingsFile_Sample)== FALSE){
            message <- paste0("The ",SettingsInfo[["shape"]], " column selected as shape in SettingsInfo was not found in SettingsFile_Sample. Please check your input.")
            logger::log_trace(paste("Error ", message, sep=""))
            stop(message)
          }
        }

        #Plot individual
        if("individual" %in% names(SettingsInfo)){
          if(SettingsInfo[["individual"]] %in% colnames(SettingsFile_Sample)== FALSE){
            message <- paste0("The ",SettingsInfo[["individual"]], " column selected as individual in SettingsInfo was not found in SettingsFile_Sample. Please check your input.")
            logger::log_trace(paste("Error ", message, sep=""))
            stop(message)
          }
        }
      }else if(PlotSettings== "Feature"){
        if("color" %in% names(SettingsInfo)){
          if(SettingsInfo[["color"]] %in% colnames(SettingsFile_Metab)== FALSE){
            message <- paste0("The ",SettingsInfo[["color"]], " column selected as color in SettingsInfo was not found in SettingsFile_Metab. Please check your input.")
            logger::log_trace(paste("Error ", message, sep=""))
            stop(message)
          }
        }

        #Plot shape
        if("shape" %in% names(SettingsInfo)){
          if(SettingsInfo[["shape"]] %in% colnames(SettingsFile_Metab)== FALSE){
            message <- paste0("The ",SettingsInfo[["shape"]], " column selected as shape in SettingsInfo was not found in SettingsFile_Metab. Please check your input.")
            logger::log_trace(paste("Error ", message, sep=""))
            stop(message)
          }
        }

        #Plot individual
        if("individual" %in% names(SettingsInfo)){
          if(SettingsInfo[["individual"]] %in% colnames(SettingsFile_Metab)== FALSE){
            message <- paste0("The ",SettingsInfo[["individual"]], " column selected as individual in SettingsInfo was not found in SettingsFile_Metab. Please check your input.")
            logger::log_trace(paste("Error ", message, sep=""))
            stop(message)
          }
        }
      }else if(PlotSettings== "Both"){
        #Plot colour sample
        if("color_Sample" %in% names(SettingsInfo)){
          if(SettingsInfo[["color_Sample"]] %in% colnames(SettingsFile_Sample)== FALSE){
            message <- paste0("The ",SettingsInfo[["color_Sample"]], " column selected as color_Sample in SettingsInfo was not found in SettingsFile_Sample. Please check your input.")
            logger::log_trace(paste("Error ", message, sep=""))
            stop(message)
          }
        }

        #Plot colour Metab
        if("color_Metab" %in% names(SettingsInfo)){
          if(SettingsInfo[["color_Metab"]] %in% colnames(SettingsFile_Metab)== FALSE){
            message <- paste0("The ",SettingsInfo[["color_Metab"]], " column selected as color_Metab in SettingsInfo was not found in SettingsFile_Metab. Please check your input.")
            logger::log_trace(paste("Error ", message, sep=""))
            stop(message)
          }
          if(sum(colnames(InputData) %in% SettingsFile_Metab$Metabolite) < length(InputData)  ){
            message <- paste0("The InputData contains metabolites not found in SettingsFile_Metab.")
            logger::log_trace(paste("Warning ", message, sep=""))
            warning(message)
          }
        }

       # Plot shape_metab
        if("shape_Metab" %in% names(SettingsInfo)){
          if(SettingsInfo[["shape_Metab"]] %in% colnames(SettingsFile_Metab)== FALSE){
            message <- paste0("The ",SettingsInfo[["shape_Metab"]], " column selected as shape_Metab in SettingsInfo was not found in SettingsFile_Metab. Please check your input.")
            logger::log_trace(paste("Error ", message, sep=""))
            stop(message)
          }
        }

        # Plot shape_metab
        if("shape_Sample" %in% names(SettingsInfo)){
          if(SettingsInfo[["shape_Sample"]] %in% colnames(SettingsFile_Metab)== FALSE){
            message <- paste0("The ",SettingsInfo[["shape_Sample"]], " column selected as shape_Metab in SettingsInfo was not found in SettingsFile_Sample. Please check your input.")
            logger::log_trace(paste("Error ", message, sep=""))
            stop(message)
          }
        }

        #Plot individual_Metab
        if("individual_Metab" %in% names(SettingsInfo)){
          if(SettingsInfo[["individual_Metab"]] %in% colnames(SettingsFile_Metab)== FALSE){
            message <- paste0("The ",SettingsInfo[["individual_Metab"]], " column selected as individual_Metab in SettingsInfo was not found in SettingsFile_Metab. Please check your input.")
            logger::log_trace(paste("Error ", message, sep=""))
            stop(message)
          }
        }

        #Plot individual_Sample
        if("individual_Sample" %in% names(SettingsInfo)){
          if(SettingsInfo[["individual_Sample"]] %in% colnames(SettingsFile_Sample)== FALSE){
            message <- paste0("The ",SettingsInfo[["individual_Sample"]], " column selected as individual_Sample in SettingsInfo was not found in SettingsFile_Sample. Please check your input.")
            logger::log_trace(paste("Error ", message, sep=""))
            stop(message)
          }
        }

      }

    }
   }

  #-------------SaveAs
  Save_as_Plot_options <- c("svg","pdf", "png")
  if(is.null(SaveAs_Plot)==FALSE){
    if(SaveAs_Plot %in% Save_as_Plot_options == FALSE){
      message <- paste0("Check input. The selected SaveAs_Plot option is not valid. Please select one of the folowwing: ",paste(Save_as_Plot_options,collapse = ", ")," or set to NULL if no plots should be saved.")
      logger::log_trace(paste("Error ", message, sep=""))
      stop(message)
  }
  }

  SaveAs_Table_options <- c("txt","csv", "xlsx", "RData")#RData = SummarizedExperiment (?)
  if(is.null(SaveAs_Table)==FALSE){
    if((SaveAs_Table %in% SaveAs_Table_options == FALSE)| (is.null(SaveAs_Table)==TRUE)){
      message <- paste0("Check input. The selected SaveAs_Table option is not valid. Please select one of the folowwing: ",paste(SaveAs_Table_options,collapse = ", ")," or set to NULL if no tables should be saved.")
      logger::log_trace(paste("Error ", message, sep=""))
      stop(message)
    }
  }

  #-------------CoRe
  if(is.logical(CoRe) == FALSE){
    message <- paste0("Check input. The CoRe value should be either =TRUE for preprocessing of Consuption/Release experiment or =FALSE if not.")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }

  #-------------Theme
  if(is.null(Theme)==FALSE){
    Theme_options <- c("theme_grey()", "theme_gray()", "theme_bw()", "theme_linedraw()", "theme_light()", "theme_dark()", "theme_minimal()", "theme_classic()", "theme_void()", "theme_test()")
    if (Theme %in% Theme_options == FALSE){
      message <- paste0("Check input. Theme option is incorrect. You can check for complete themes here: https://ggplot2.tidyverse.org/reference/ggtheme.html. Options are the following: ",paste(Theme_options, collapse = ", "),"." )
      logger::log_trace(paste("Error ", message, sep=""))
      stop(message)
    }
  }
  #------------- general
  if(is.logical(PrintPlot) == FALSE){
    message <- paste0("Check input. PrintPlot should be either =TRUE or =FALSE.")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }
}

################################################################################################
### ### ### PreProcessing helper function: Internal Function to check function input ### ### ###
################################################################################################

#' Check specific input parameters for PreProcessing()
#'
#' @param SettingsFile_Sample DF which contains information about the samples, which will be combined with the input data based on the unique sample identifiers used as rownames.
#' @param SettingsInfo Named vector containing the information about the names of the experimental parameters. c(Conditions="ColumnName_Plot_SettingsFile", Biological_Replicates="ColumnName_Plot_SettingsFile"). For CoRe = TRUE a CoRe_norm_factor = "Columnname_Input_SettingsFile" and CoRe_media = "Columnname_Input_SettingsFile", have to also be added.
#' @param  \emph{Optional: }CoRe If TRUE, a consumption-release experiment has been performed and the CoRe value will be calculated.\strong{Default = FALSE}
#' @param FeatureFilt \emph{Optional: }If NULL, no feature filtering is performed. If set to "Standard" then it applies the 80%-filtering rule (Bijlsma S. et al., 2006) on the metabolite features on the whole dataset. If is set to "Modified",filtering is done based on the different conditions, thus a column named "Conditions" must be provided in the Input_SettingsFile input file including the individual conditions you want to apply the filtering to (Yang, J et al., 2015). \strong{Default = "Standard"}
#' @param FeatureFilt_Value \emph{Optional: } Percentage of feature filtering. \strong{Default = 0.8}
#' @param TIC \emph{Optional: } If TRUE, Total Ion Count normalization is performed. \strong{Default = TRUE}
#' @param MVI \emph{Optional: } If TRUE, Missing Value Imputation (MVI) based on half minimum is performed \strong{Default = TRUE}
#' @param MVI_Percentage \emph{Optional: } Percentage 0-100 of imputed value based on the minimum value. \strong{Default = 50}
#' @param HotellinsConfidence \emph{Optional: } Defines the Confidence of Outlier identification in HotellingT2 test. Must be numeric.\strong{Default = 0.99}
#'
#' @return returns warnings and errors if input is not correct
#'
#' @keywords Input check for MetaProViz::PreProcessing
#'
#' @importFrom logger log_trace
#'
#' @noRd
#'
CheckInput_PreProcessing <- function(SettingsFile_Sample,
                                     SettingsInfo,
                                     CoRe=FALSE,
                                     FeatureFilt = "Modified",
                                     FeatureFilt_Value = 0.8,
                                     TIC = TRUE,
                                     MVI= TRUE,
                                     MVI_Percentage=50,
                                     HotellinsConfidence = 0.99){
  if(is.vector(SettingsInfo)==TRUE){
    #-------------SettingsInfo
    #CoRe
    if(CoRe == TRUE){   # parse CoRe normalisation factor
      message <- paste0("For Consumption Release experiment we are using the method from Jain M.  REF: Jain et. al, (2012), Science 336(6084):1040-4, doi: 10.1126/science.1218595.")
      logger::log_trace(paste("Message ", message, sep=""))
      message(message)
      if("CoRe_media" %in% names(SettingsInfo)){
        if(length(grep(SettingsInfo[["CoRe_media"]], SettingsFile_Sample[[SettingsInfo[["Conditions"]]]])) < 1){     # Check for CoRe_media samples
          message <- paste0("No CoRe_media samples were provided in the 'Conditions' in the SettingsFile_Sample. For a CoRe experiment control media samples without cells have to be measured and be added in the 'Conditions'
                            column labeled as 'CoRe_media' (see @param section). Please make sure that you used the correct labelling or whether you need CoRe = FALSE for your analysis")
          logger::log_trace(paste("Error ", message, sep=""))
          stop(message)
        }
      }

      if ("CoRe_norm_factor" %in% names(SettingsInfo)==FALSE){
        message <- paste0("No growth rate or growth factor provided for normalising the CoRe result, hence CoRe_norm_factor set to 1 for each sample")
        logger::log_trace(paste("Warning ", message, sep=""))
        warning(message)
      }
    }
  }

  #-------------General parameters
  Feature_Filtering_options <- c("Standard","Modified")
  if(FeatureFilt %in% Feature_Filtering_options == FALSE & is.null(FeatureFilt)==FALSE){
    message <- paste0("Check input. The selected FeatureFilt option is not valid. Please set to NULL or select one of the folowwing: ",paste(Feature_Filtering_options,collapse = ", "),"." )
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }
  if(is.numeric(FeatureFilt_Value) == FALSE |FeatureFilt_Value > 1 | FeatureFilt_Value < 0){
    message <- paste0("Check input. The selected FeatureFilt_Value should be numeric and between 0 and 1.")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }
  if(is.logical(TIC) == FALSE){
    message <- paste0("Check input. The TIC value should be either `TRUE` if TIC normalization is to be performed or `FALSE` if no data normalization is to be applied.")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }
  if(is.logical(MVI) == FALSE){
    message <- paste0("Check input. The MVI value should be either `TRUE` if mising value imputation should be performed or `FALSE` if not.")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }
  if(is.numeric(MVI_Percentage)== FALSE |HotellinsConfidence > 100 | HotellinsConfidence < 0){
    message <- paste0("Check input. The selected MVI_Percentage value should be numeric and between 0 and 100.")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }
  if( is.numeric(HotellinsConfidence)== FALSE |HotellinsConfidence > 1 | HotellinsConfidence < 0){
    message <- paste0("Check input. The selected HotellinsConfidence value should be numeric and between 0 and 1.")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }
}

################################################################################################
### ### ### DMA helper function: Internal Function to check function input ### ### ###
################################################################################################

#' Check input parameters of DMA
#'
#' @param InputData DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names.
#' @param SettingsFile_Sample DF which contains metadata information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames.
#' @param SettingsInfo Named vector including the information about the conditions column information on numerator or denominator c(Conditions="ColumnName_SettingsFile", Numerator = "ColumnName_SettingsFile", Denominator  = "ColumnName_SettingsFile"). \strong{Default = c(conditions="Conditions", numerator = NULL, denumerator = NULL)}
#' @param StatPval \emph{Optional: } String which contains an abbreviation of the selected test to calculate p.value. For one-vs-one comparisons choose t.test, wilcox.test, "chisq.test", "cor.test" or lmFit (=limma), for one-vs-all or all-vs-all comparison choose aov (=anova), welch(=welch anova), kruskal.test or lmFit (=limma) \strong{Default = "lmFit"}
#' @param StatPadj \emph{Optional: } String which contains an abbreviation of the selected p.adjusted test for p.value correction for multiple Hypothesis testing. Search: ?p.adjust for more methods:"BH", "fdr", "bonferroni", "holm", etc.\strong{Default = "fdr"}
#' @param VST TRUE or FALSE for whether to use variance stabilizing transformation on the data when linear modeling is used for hypothesis testing. \strong{Default = FALSE}
#' @param PerformShapiro TRUE or FALSE for whether to perform the shapiro.test and get informed about data distribution (normal versus not-normal distribution. \strong{Default = TRUE}
#' @param PerformBartlett TRUE or FALSE for whether to perform the bartlett.test. \strong{Default = TRUE}
#' @param Transform TRUE or FALSE. If TRUE we expect the data to be not log2 transformed and log2 transformation will be performed within the limma function and Log2FC calculation. If FALSE we expect the data to be log2 transformed as this impacts the Log2FC calculation and limma. \strong{Default= TRUE}
#'
#' @return returns warnings and errors if input is not correct
#'
#' @keywords Input check for MetaProViz::DMA
#'
#' @importFrom logger log_trace
#' @importFrom dplyr filter select_if
#' @importFrom magrittr %>%
#' @importFrom utils combn
#'
#' @noRd
#'
CheckInput_DMA <- function(InputData,
                           SettingsFile_Sample,
                           SettingsInfo= c(Conditions="Conditions", Numerator = NULL, Denominator  = NULL),
                           StatPval ="lmFit",
                           StatPadj="fdr",
                           VST=FALSE,
                           PerformShapiro =TRUE,
                           PerformBartlett =TRUE,
                           Transform=TRUE){

  #-------------SettingsInfo
  if(is.null(SettingsInfo)==TRUE){
    message <- paste0("You have to provide SettingsInfo's for Conditions.")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }

  ## ------------ Denominator/numerator ----------- ##
  # Denominator and numerator: Define if we compare one_vs_one, one_vs_all or all_vs_all.
  if("Denominator" %in% names(SettingsInfo)==FALSE  & "Numerator" %in% names(SettingsInfo) ==FALSE){
    # all-vs-all: Generate all pairwise combinations
    conditions = SettingsFile_Sample[[SettingsInfo[["Conditions"]]]]
    denominator <-unique(SettingsFile_Sample[[SettingsInfo[["Conditions"]]]])
    numerator <-unique(SettingsFile_Sample[[SettingsInfo[["Conditions"]]]])
    comparisons <- utils::combn(unique(conditions), 2) %>% as.matrix()
    #Settings:
    MultipleComparison = TRUE
    all_vs_all = TRUE
  }else if("Denominator" %in% names(SettingsInfo)==TRUE  & "Numerator" %in% names(SettingsInfo)==FALSE){
    #all-vs-one: Generate the pairwise combinations
    conditions = SettingsFile_Sample[[SettingsInfo[["Conditions"]]]]
    denominator <- SettingsInfo[["Denominator"]]
    numerator <-unique(SettingsFile_Sample[[SettingsInfo[["Conditions"]]]])
    # Remove denom from num
    numerator <- numerator[!numerator %in% denominator]
    comparisons  <- t(expand.grid(numerator, denominator)) %>% as.data.frame()
    #Settings:
    MultipleComparison = TRUE
    all_vs_all = FALSE
  }else if("Denominator" %in% names(SettingsInfo)==TRUE  & "Numerator" %in% names(SettingsInfo)==TRUE){
    # one-vs-one: Generate the comparisons
    denominator <- SettingsInfo[["Denominator"]]
    numerator <- SettingsInfo[["Numerator"]]
    comparisons <- matrix(c(SettingsInfo[["Denominator"]], SettingsInfo[["Numerator"]]))
    #Settings:
    MultipleComparison = FALSE
    all_vs_all = FALSE
  }

  ## ------------ Test statistics ----------- ##
  if(MultipleComparison==FALSE){
    STAT_pval_options <- c("t.test", "wilcox.test","chisq.test", "cor.test", "lmFit")
    if(StatPval %in% STAT_pval_options == FALSE){
      message <- paste0("Check input. The selected StatPval option for Hypothesis testing is not valid for multiple comparison (one-vs-all or all-vs-all). Please select one of the following: ",paste(STAT_pval_options,collapse = ", ")," or specify numerator and denumerator." )
      logger::log_trace(paste("Error ", message, sep=""))
      stop(message)
    }
  }else{
    STAT_pval_options <- c("aov", "kruskal.test", "welch" ,"lmFit")
    if(StatPval %in% STAT_pval_options == FALSE){
      message <- paste0("Check input. The selected StatPval option for Hypothesis testing is not valid for one-vs-one comparsion. Multiple comparison is selected. Please select one of the following: ",paste(STAT_pval_options,collapse = ", ")," or change numerator and denumerator." )
      logger::log_trace(paste("Error ", message, sep=""))
      stop(message)
    }
  }

  STAT_padj_options <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
  if(StatPadj %in% STAT_padj_options == FALSE){
    message <- paste0("Check input. The selected StatPadj option for multiple Hypothesis testing correction is not valid. Please select one of the folowing: ",paste(STAT_padj_options,collapse = ", "),"." )
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }

  ## ------------ Sample Numbers ----------- ##
  Num <- InputData %>%#Are sample numbers enough?
    dplyr::filter(SettingsFile_Sample[[SettingsInfo[["Conditions"]]]] %in% numerator) %>%
    dplyr::select_if(is.numeric)#only keep numeric columns with metabolite values
  Denom <- InputData %>%
    dplyr::filter(SettingsFile_Sample[[SettingsInfo[["Conditions"]]]] %in% denominator) %>%
    dplyr::select_if(is.numeric)

  if(nrow(Num)==1){
    message <- paste0("There is only one sample available for ", numerator, ", so no statistical test can be performed.")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  } else if(nrow(Denom)==1){
    message <- paste0("There is only one sample available for ", denominator, ", so no statistical test can be performed.")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }else if(nrow(Num)==0){
    message <- paste0("There is no sample available for ", numerator, ".")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }else if(nrow(Denom)==0){
    message <- paste0("There is no sample available for ", denominator, ".")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }

  ## ------------ Check Missingness ------------- ##
  Num_Miss <- replace(Num, Num==0, NA)
  Num_Miss <- Num_Miss[, (colSums(is.na(Num_Miss)) > 0), drop = FALSE]

  Denom_Miss <- replace(Denom, Denom==0, NA)
  Denom_Miss <- Denom_Miss[, (colSums(is.na(Denom_Miss)) > 0), drop = FALSE]

  if((ncol(Num_Miss)>0 & ncol(Denom_Miss)==0)){
    Metabolites_Miss <- colnames(Num_Miss)
    if(ncol(Num_Miss)<=10){
      message <- paste0("In `Numerator` ",paste0(toString(numerator)), ", NA/0 values exist in ", ncol(Num_Miss), " Metabolite(s): ", paste0(colnames(Num_Miss), collapse = ", "), ". Those metabolite(s) might return p.val= NA, p.adj.= NA, t.val= NA. The Log2FC = Inf, if all replicates are 0/NA.")
      logger::log_info(message)
      message(message)
    }else{
      message <- paste0("In `Numerator` ",paste0(toString(numerator)), ", NA/0 values exist in ", ncol(Num_Miss), " Metabolite(s).", " Those metabolite(s) might return p.val= NA, p.adj.= NA, t.val= NA. The Log2FC = Inf, if all replicates are 0/NA.")
      logger::log_info(message)
      message(message)
    }
  } else if(ncol(Num_Miss)==0 & ncol(Denom_Miss)>0){
    Metabolites_Miss <- colnames(Denom_Miss)
    if(ncol(Num_Miss)<=10){
      message <- paste0("In `Denominator` ",paste0(toString(denominator)), ", NA/0 values exist in ", ncol(Denom_Miss), " Metabolite(s): ", paste0(colnames(Denom_Miss), collapse = ", "), ". Those metabolite(s) might return p.val= NA, p.adj.= NA, t.val= NA. The Log2FC = Inf, if all replicates are 0/NA.")
      logger::log_info(message)
      message(message)
    }else{
      message <- paste0("In `Denominator` ",paste0(toString(denominator)), ", NA/0 values exist in ", ncol(Denom_Miss), " Metabolite(s).", " Those metabolite(s) might return p.val= NA, p.adj.= NA, t.val= NA. The Log2FC = Inf, if all replicates are 0/NA.")
      logger::log_info(message)
      message(message)#
    }
  } else if(ncol(Num_Miss)>0 & ncol(Denom_Miss)>0){
    Metabolites_Miss <- c(colnames(Num_Miss), colnames(Denom_Miss))
    Metabolites_Miss <- unique(Metabolites_Miss)

    message <- paste0("In `Numerator` ",paste0(toString(numerator)), ", NA/0 values exist in ", ncol(Num_Miss), " Metabolite(s).", " and in `denominator`",paste0(toString(denominator)), " ",ncol(Denom_Miss), " Metabolite(s).",
                      ". Those metabolite(s) might return p.val= NA, p.adj.= NA, t.val= NA. The Log2FC = Inf, if all replicates are 0/NA.")
    logger::log_info(message)
    message(message)
  } else{
    message <- paste0("There are no NA/0 values")
    logger::log_info(message)
    message(message)

    Metabolites_Miss <- c(colnames(Num_Miss), colnames(Denom_Miss))
    Metabolites_Miss <- unique(Metabolites_Miss)
  }

  #-------------General parameters
  if(is.logical(VST) == FALSE){
    message <- paste0("Check input. The VST value should be either =TRUE or =FALSE.")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }

  if(is.logical(PerformShapiro) == FALSE){
    message <- paste0("Check input. The Shapiro value should be either =TRUE or =FALSE.")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }
  if(is.logical(PerformBartlett) == FALSE){
    message <- paste0("Check input. The Bartlett value should be either =TRUE or =FALSE.")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }
  if(is.logical(Transform) == FALSE){
    message <- paste0("Check input. `Transform` should be either =TRUE or =FALSE.")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }

  #Settings <- list("comparisons"=comparisons, "MultipleComparison"=MultipleComparison, "all_vs_all"=all_vs_all, "Metabolites_Miss"=Metabolites_Miss, "denominator"=denominator, "numerator"=numerator)
  #return(invisible(Settings))
}

################################################################################################
### ### ### ORA helper function: Internal Function to check function input ### ### ###
################################################################################################

#' Check input parameters of ORA
#'
#' @param InputData Passed to main function PreProcessing()
#' @param SettingsInfo Passed to main function PreProcessing()
#'
#' @keywords Input check
#' @noRd
#'
#'
CheckInput_ORA <- function(InputData,
                           SettingsInfo,
                           RemoveBackground,
                           PathwayFile,
                           PathwayName,
                           minGSSize,
                           maxGSSize,
                           SaveAs_Table,
                           pCutoff,
                           PercentageCutoff
){
  # 1. The input data:
  if(class(InputData) != "data.frame"){
    stop("InputData should be a data.frame. It's currently a ", paste(class(InputData), ".",sep = ""))
  }
  if(any(duplicated(row.names(InputData)))==TRUE){
    stop("Duplicated row.names of InputData, whilst row.names must be unique")
  }

  # 2. Settings Columns:
  if(is.vector(SettingsInfo)==FALSE & is.null(SettingsInfo)==FALSE){
    stop("SettingsInfo should be NULL or a vector. It's currently a ", paste(class(SettingsInfo), ".", sep = ""))
  }

  if(is.null(SettingsInfo)==FALSE){
    #"ClusterColumn"
    if("ClusterColumn" %in% names(SettingsInfo)){
      if(SettingsInfo[["ClusterColumn"]] %in% colnames(InputData)== FALSE){
        stop("The ", SettingsInfo[["ClusterColumn"]], " column selected as ClusterColumn in SettingsInfo was not found in InputData. Please check your input.")
      }
    }

    #"BackgroundColumn"
    if("BackgroundColumn" %in% names(SettingsInfo)){
      if(SettingsInfo[["BackgroundColumn"]] %in% colnames(InputData)== FALSE){
        stop("The ", SettingsInfo[["BackgroundColumn"]], " column selected as BackgroundColumn in SettingsInfo was not found in InputData. Please check your input.")
      }
    }

    #"pvalColumn"
    if("pvalColumn" %in% names(SettingsInfo)){
      if(SettingsInfo[["pvalColumn"]] %in% colnames(InputData)== FALSE){
        stop("The ", SettingsInfo[["pvalColumn"]], " column selected as pvalColumn in SettingsInfo was not found in InputData. Please check your input.")
      }
    }

    #"PercentageColumn"
    if("PercentageColumn" %in% names(SettingsInfo)){
      if(SettingsInfo[["PercentageColumn"]] %in% colnames(InputData)== FALSE){
        stop("The ", SettingsInfo[["PercentageColumn"]], " column selected as PercentageColumn in SettingsInfo was not found in InputData. Please check your input.")
      }
    }


    #"PathwayTerm"
    if("PathwayTerm" %in% names(SettingsInfo)){
      if(SettingsInfo[["PathwayTerm"]] %in% colnames(PathwayFile)== FALSE){
        stop("The ", SettingsInfo[["PathwayTerm"]], " column selected as PathwayTerm in SettingsInfo was not found in PathwayFile. Please check your input.")
      }else{
        PathwayFile <- PathwayFile%>%
          dplyr::rename("term"=SettingsInfo[["PathwayTerm"]])
        PathwayFile$Description <- PathwayFile$term
      }
    }else{
      stop("SettingsInfo must provide the column name for PathwayTerm in PathwayFile")
    }

    # PathwayFeature
    if("PathwayFeature" %in% names(SettingsInfo)){
      if(SettingsInfo[["PathwayFeature"]] %in% colnames(PathwayFile)== FALSE){
        stop("The ", SettingsInfo[["PathwayFeature"]], " column selected as PathwayFeature in SettingsInfo was not found in PathwayFile. Please check your input.")
      }else{
        PathwayFile <- PathwayFile%>%
          dplyr::rename("gene"=SettingsInfo[["PathwayFeature"]])
      }
    }else{
      stop("SettingsInfo must provide the column name for PathwayFeature in PathwayFile")
    }

  }else{
    stop("you must provide SettingsInfo.")
  }

  # 3. General Settings
  if(is.character(PathwayName)==FALSE){
    stop("Check input. PathwayName must be a character of syntax 'example'.")
  }

  if(is.logical(RemoveBackground) == FALSE){
    stop("Check input. RemoveBackground value should be either =TRUE or = FALSE.")
  }

  if(is.numeric(minGSSize)== FALSE){
    stop("Check input. The selected minGSSize value should be numeric.")
  }

  if(is.numeric(maxGSSize)== FALSE){
    stop("Check input. The selected maxGSSize value should be numeric.")
  }

  SaveAs_Table_options <- c("txt","csv", "xlsx", "RData")#RData = SummarizedExperiment (?)
  if(is.null(SaveAs_Table)==FALSE){
    if((SaveAs_Table %in% SaveAs_Table_options == FALSE)| (is.null(SaveAs_Table)==TRUE)){
      stop("Check input. The selected SaveAs_Table option is not valid. Please select one of the folowwing: ",paste(SaveAs_Table_options,collapse = ", "),"." )
    }
  }

  if(is.null(pCutoff)== FALSE){
    if(is.numeric(pCutoff)== FALSE | pCutoff > 1 |  pCutoff < 0){
      stop("Check input. The selected Plot_pCutoff value should be numeric and between 0 and 1.")
    }
  }

  if(is.null(PercentageCutoff)== FALSE){
    if( is.numeric(PercentageCutoff)== FALSE  | PercentageCutoff > 100 | PercentageCutoff < 0){
      stop("Check input. The selected PercentageCutoff value should be numeric and between 0 and 100.")
    }
  }


  ## -------- Return Pathways ---------##
  return(invisible(PathwayFile))
}


################################################################################################
### ### ### MCA Helper function: Internal Function to check function input ### ### ###
################################################################################################

#' Check input parameters
#'
#' @param InputData_C1 Passed to main function MCA. If not avaliable can be set to NULL.
#' @param InputData_C2 Passed to main function MCA. If not avaliable can be set to NULL.
#' @param InputData_Intra Passed to main function MCA. If not avaliable can be set to NULL.
#' @param InputData_CoRe Passed to main function MCA. If not avaliable can be set to NULL.
#' @param SettingsInfo_C1 Passed to main function MCA. If not avaliable can be set to NULL.
#' @param SettingsInfo_C2 Passed to main function MCA. If not avaliable can be set to NULL.
#' @param SettingsInfo_Intra Passed to main function MCA. If not avaliable can be set to NULL.
#' @param SettingsInfo_CoRe Passed to main function MCA. If not avaliable can be set to NULL.
#' @param BackgroundMethod Passed to main function MCA.
#' @param FeatureID Passed to main function MCA.
#' @param SaveAs_Table Passed to main function PreProcessing(). If not avaliable can be set to NULL.
#'
#' @param Function Name of the MetaProViz Function that is checked.
#' @param InputList
#'
#'
#' @keywords Input check
#' @noRd
#'
#'

CheckInput_MCA <- function(InputData_C1,
                           InputData_C2,
                           InputData_CoRe,
                           InputData_Intra,
                           SettingsInfo_C1,
                           SettingsInfo_C2,
                           SettingsInfo_CoRe,
                           SettingsInfo_Intra,
                           BackgroundMethod,
                           FeatureID,
                           SaveAs_Table
){
  ## ------------ Create log file ----------- ##
  MetaProViz_Init()

  #------------- InputData
  if(is.null(InputData_C1)==FALSE){
    if(class(InputData_C1) != "data.frame"| class(InputData_C2) != "data.frame"){
      stop("InputData_C1 and InputData_C2 should be a data.frame. It's currently a ", paste(class(InputData_C1)), paste(class(InputData_C2)), ".",sep = "")
    }
    if(length(InputData_C1[duplicated(InputData_C1[[FeatureID]]), FeatureID]) > 0){
      stop("Duplicated FeatureIDs of InputData_C1, whilst features must be unique")
    }
    if(length(InputData_C2[duplicated(InputData_C2[[FeatureID]]), FeatureID]) > 0){
      stop("Duplicated FeatureIDs of InputData_C2, whilst features must be unique")
    }

  }else{
    if(class(InputData_Intra) != "data.frame"| class(InputData_CoRe) != "data.frame"){
      stop("InputData_Intra and InputData_CoRe should be a data.frame. It's currently a ", paste(class(InputData_Intra)), paste(class(InputData_CoRe)), ".",sep = "")
    }
    if(length(InputData_Intra[duplicated(InputData_Intra[[FeatureID]]), FeatureID]) > 0){
      stop("Duplicated FeatureIDs of InputData_Intra, whilst features must be unique")
    }
    if(length(InputData_CoRe[duplicated(InputData_CoRe[[FeatureID]]), FeatureID]) > 0){
      stop("Duplicated FeatureIDs of InputData_CoRe, whilst features must be unique")
    }
  }


  #------------- SettingsInfo
  if(is.null(SettingsInfo_C1)==FALSE){
    ## C1
    #ValueCol
    if("ValueCol" %in% names(SettingsInfo_C1)){
      if(SettingsInfo_C1[["ValueCol"]] %in% colnames(InputData_C1)== FALSE){
        stop("The ", SettingsInfo_C1[["ValueCol"]], " column selected as ValueCol in SettingsInfo_C1 was not found in InputData_C1. Please check your input.")
      }
    }
    #StatCol
    if("StatCol" %in% names(SettingsInfo_C1)){
      if(SettingsInfo_C1[["StatCol"]] %in% colnames(InputData_C1)== FALSE){
        stop("The ", SettingsInfo_C1[["StatCol"]], " column selected as StatCol in SettingsInfo_C1 was not found in InputData_C1. Please check your input.")
      }
    }

    ## C2
    #ValueCol
    if("ValueCol" %in% names(SettingsInfo_C2)){
      if(SettingsInfo_C2[["ValueCol"]] %in% colnames(InputData_C2)== FALSE){
        stop("The ", SettingsInfo_C2[["ValueCol"]], " column selected as ValueCol in SettingsInfo_C2 was not found in InputData_C2. Please check your input.")
      }
    }
    #StatCol
    if("StatCol" %in% names(SettingsInfo_C2)){
      if(SettingsInfo_C2[["StatCol"]] %in% colnames(InputData_C2)== FALSE){
        stop("The ", SettingsInfo_C2[["StatCol"]], " column selected as StatCol in SettingsInfo_C2 was not found in InputData_C2. Please check your input.")
      }
    }
  }else{
    ## Intra
    #ValueCol
    if("ValueCol" %in% names(SettingsInfo_Intra)){
      if(SettingsInfo_Intra[["ValueCol"]] %in% colnames(InputData_Intra)== FALSE){
        stop("The ", SettingsInfo_Intra[["ValueCol"]], " column selected as ValueCol in SettingsInfo_Intra was not found in InputData_Intra. Please check your input.")
      }
    }
    #StatCol
    if("StatCol" %in% names(SettingsInfo_Intra)){
      if(SettingsInfo_Intra[["StatCol"]] %in% colnames(InputData_Intra)== FALSE){
        stop("The ", SettingsInfo_Intra[["StatCol"]], " column selected as StatCol in SettingsInfo_Intra was not found in InputData_Intra. Please check your input.")
      }
    }

    ## CoRe
    #ValueCol
    if("ValueCol" %in% names(SettingsInfo_CoRe)){
      if(SettingsInfo_CoRe[["ValueCol"]] %in% colnames(InputData_CoRe)== FALSE){
        stop("The ", SettingsInfo_CoRe[["ValueCol"]], " column selected as ValueCol in SettingsInfo_CoRe was not found in InputData_CoRe. Please check your input.")
      }
    }
    #StatCol
    if("StatCol" %in% names(SettingsInfo_CoRe)){
      if(SettingsInfo_CoRe[["StatCol"]] %in% colnames(InputData_CoRe)== FALSE){
        stop("The ", SettingsInfo_CoRe[["StatCol"]], " column selected as StatCol in SettingsInfo_CoRe was not found in InputData_CoRe. Please check your input.")
      }
    }

    #StatCol
    if("DirectionCol" %in% names(SettingsInfo_CoRe)){
      if(SettingsInfo_CoRe[["DirectionCol"]] %in% colnames(InputData_CoRe)== FALSE){
        stop("The ", SettingsInfo_CoRe[["DirectionCol"]], " column selected as DirectionCol in SettingsInfo_CoRe was not found in InputData_CoRe. Please check your input.")
      }
    }

  }

  #------------- SettingsInfo Cutoffs:
  if(is.null(SettingsInfo_C1)==FALSE){
    if(is.na(as.numeric(SettingsInfo_C1[["StatCutoff"]])) == TRUE |as.numeric(SettingsInfo_C1[["StatCutoff"]]) > 1 | as.numeric(SettingsInfo_C1[["StatCutoff"]]) < 0){
      stop("Check input. The selected StatCutoff in SettingsInfo_C1 should be numeric and between 0 and 1.")
    }

    if(is.na(as.numeric(SettingsInfo_C2[["StatCutoff"]])) == TRUE |as.numeric(SettingsInfo_C2[["StatCutoff"]]) > 1 | as.numeric(SettingsInfo_C2[["StatCutoff"]]) < 0){
      stop("Check input. The selected StatCutoff in SettingsInfo_C2 should be numeric and between 0 and 1.")
    }

    if(is.na(as.numeric(SettingsInfo_C1[["ValueCutoff"]])) == TRUE){
      stop("Check input. The selected ValueCutoff in SettingsInfo_C1 should be numeric and between 0 and 1.")
    }

    if(is.na(as.numeric(SettingsInfo_C2[["ValueCutoff"]])) == TRUE){
      stop("Check input. The selected ValueCutoff in SettingsInfo_C2 should be numeric and between 0 and 1.")
    }

  }else{
    if(is.na(as.numeric(SettingsInfo_Intra[["StatCutoff"]])) == TRUE |as.numeric(SettingsInfo_Intra[["StatCutoff"]]) > 1 | as.numeric(SettingsInfo_Intra[["StatCutoff"]]) < 0){
      stop("Check input. The selected StatCutoff in SettingsInfo_Intra should be numeric and between 0 and 1.")
    }

    if(is.na(as.numeric(SettingsInfo_CoRe[["StatCutoff"]])) == TRUE |as.numeric(SettingsInfo_CoRe[["StatCutoff"]]) > 1 | as.numeric(SettingsInfo_CoRe[["StatCutoff"]]) < 0){
      stop("Check input. The selected StatCutoff in SettingsInfo_CoRe should be numeric and between 0 and 1.")
    }

    if(is.na(as.numeric(SettingsInfo_Intra[["ValueCutoff"]])) == TRUE){
      stop("Check input. The selected ValueCutoff in SettingsInfo_Intra should be numeric and between 0 and 1.")
    }

    if(is.na(as.numeric(SettingsInfo_CoRe[["ValueCutoff"]])) == TRUE){
      stop("Check input. The selected ValueCutoff in SettingsInfo_CoRe should be numeric and between 0 and 1.")
    }
  }

  #------------ NAs in data
  if(is.null(InputData_C1)==FALSE){
    if(nrow(InputData_C1[complete.cases(InputData_C1[[SettingsInfo_C1[["ValueCol"]]]], InputData_C1[[SettingsInfo_C1[["StatCol"]]]]), ]) < nrow(InputData_C1)){
      warning("InputData_C1 includes NAs in ", SettingsInfo_C1[["ValueCol"]], " and/or in ", SettingsInfo_C1[["StatCol"]], ". ", nrow(InputData_C1)- nrow(InputData_C1[complete.cases(InputData_C1[[SettingsInfo_C1[["ValueCol"]]]], InputData_C1[[SettingsInfo_C1[["StatCol"]]]]), ]) ," metabolites containing NAs are removed.")
    }

    if(nrow(InputData_C2[complete.cases(InputData_C2[[SettingsInfo_C2[["ValueCol"]]]], InputData_C2[[SettingsInfo_C2[["StatCol"]]]]), ]) < nrow(InputData_C2)){
      warning("InputData_C2 includes NAs in ", SettingsInfo_C2[["ValueCol"]], " and/or in", SettingsInfo_C2[["StatCol"]], ". ", nrow(InputData_C2)- nrow(InputData_C2[complete.cases(InputData_C2[[SettingsInfo_C2[["ValueCol"]]]], InputData_C2[[SettingsInfo_C2[["StatCol"]]]]), ]) ," metabolites containing NAs are removed.")
    }
  }else{
    if(nrow(InputData_Intra[complete.cases(InputData_Intra[[SettingsInfo_Intra[["ValueCol"]]]], InputData_Intra[[SettingsInfo_Intra[["StatCol"]]]]), ]) < nrow(InputData_Intra)){
      warning("InputData_Intra includes NAs in ", SettingsInfo_Intra[["ValueCol"]], " and/or in ", SettingsInfo_Intra[["StatCol"]], ". ", nrow(InputData_Intra)- nrow(InputData_Intra[complete.cases(InputData_Intra[[SettingsInfo_Intra[["ValueCol"]]]], InputData_Intra[[SettingsInfo_Intra[["StatCol"]]]]), ]) ," metabolites containing NAs are removed.")
    }

    if(nrow(InputData_CoRe[complete.cases(InputData_CoRe[[SettingsInfo_CoRe[["ValueCol"]]]], InputData_CoRe[[SettingsInfo_CoRe[["StatCol"]]]]), ]) < nrow(InputData_CoRe)){
      warning("InputData_CoRe includes NAs in ", SettingsInfo_CoRe[["ValueCol"]], " and/or in ", SettingsInfo_CoRe[["StatCol"]], ". ", nrow(InputData_CoRe)- nrow(InputData_CoRe[complete.cases(InputData_CoRe[[SettingsInfo_CoRe[["ValueCol"]]]], InputData_CoRe[[SettingsInfo_CoRe[["StatCol"]]]]), ]) ," metabolites containing NAs are removed.")
    }
  }

  #------------- BackgroundMethod
  if(is.null(SettingsInfo_C1)==FALSE){
    options <- c("C1|C2", "C1&C2", "C2", "C1" , "*")
    if(any(options %in% BackgroundMethod) == FALSE){
      stop("Check input. The selected BackgroundMethod option is not valid. Please select one of the folowwing: ",paste(options,collapse = ", "),"." )
    }
  }else{
    options <- c("Intra|CoRe", "Intra&CoRe", "CoRe", "Intra" , "*")
    if(any(options %in% BackgroundMethod) == FALSE){
      stop("Check input. The selected BackgroundMethod option is not valid. Please select one of the folowwing: ",paste(options,collapse = ", "),"." )
    }
  }

  #------------- SaveAs
  SaveAs_Table_options <- c("txt","csv", "xlsx", "RData")#RData = SummarizedExperiment (?)
  if(is.null(SaveAs_Table)==FALSE){
    if((SaveAs_Table %in% SaveAs_Table_options == FALSE)| (is.null(SaveAs_Table)==TRUE)){
      stop("Check input. The selected SaveAs_Table option is not valid. Please select one of the folowwing: ",paste(SaveAs_Table_options,collapse = ", "),"." )
    }
  }
}


