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
# REFACT: Description of each argument should start with its type; e.g.
# "@param x Character: name of the variable mapped to the x axis."
#' @param InputData Passed to main function Function().
#' @param InputData_Num  \emph{Optional: } If InputData must be numeric \strong{Default = TRUE}
#' @param SettingsFile_Sample \emph{Optional: } Passed to main function Function(). If not avaliable can be set to NULL. \strong{Default = NULL}
#' @param SettingsFile_Metab \emph{Optional: } Passed to main function Function(). If not avaliable can be set to NULL. \strong{Default = NULL}
#' @param SettingsInfo \emph{Optional: } Passed to main function Function() \strong{Default = NULL}
#' @param SaveAs_Plot \emph{Optional: } Passed to main function Function(). If not avaliable can be set to NULL.\strong{Default = NULL}
#' @param SaveAs_Table \emph{Optional: } Passed to main function Function(). If not avaliable can be set to NULL.\strong{Default = NULL}
#' @param CoRe \emph{Optional: } Passed to main function Function(). If not avaliable can be set to NULL. \strong{Default = FALSE}
#' @param PrintPlot \emph{Optional: } Passed to main function Function(). If not avaliable can be set to NULL. \strong{Default = FALSE}
#' @param Theme \emph{Optional: } Passed to main function Function(). If not avaliable can be set to NULL. \strong{Default = NULL}
#' @param PlotSettings \emph{Optional: } Needs to be set for VizX functions. Options are "Sample", "Feature", Both". This refers to SettingsInfo color, shape, individual as for some plots we have both feature and sample settings. \strong{Default = NULL}
#'
#'
#' @keywords Input check
#' @noRd
#'
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
  # REFACT: this will fail with "condition has length > 1" error, and also
  # gives wrong result. Use is.data.frame() instead.
  if(class(InputData) != "data.frame"){
    stop("InputData should be a data.frame. It's currently a ", paste(class(InputData), ".",sep = ""))
  }
  if(any(duplicated(row.names(InputData)))==TRUE){
    stop("Duplicated row.names of InputData, whilst row.names must be unique")
  }

  if(InputData_Num==TRUE){
     Test_num <- apply(InputData, 2, function(x) is.numeric(x))
     if((any(Test_num) ==  FALSE) ==  TRUE){
       stop("InputData needs to be of class numeric")
       }
  }


  if(sum(duplicated(colnames(InputData))) > 0){
    doublons <- as.character(colnames(InputData)[duplicated(colnames(InputData))])#number of duplications
    #data <-data[!duplicated(colnames(InputData)),]#remove duplications
    stop("InputData contained duplicates column names, whilst col.names must be unique.")
  }

  #-------------SettingsFile
  if(is.null(SettingsFile_Sample)==FALSE){
    Test_match <- merge(SettingsFile_Sample, InputData, by = "row.names", all =  FALSE)
    if(nrow(Test_match) ==  0){
        stop("row.names InputData need to match row.names SettingsFile_Sample.")
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
    stop("SettingsInfo should be NULL or a vector. It's currently a ", paste(class(SettingsInfo), ".", sep = ""))
  }

  if(is.null(SettingsInfo)==FALSE){
    #Conditions
    if("Conditions" %in% names(SettingsInfo)){
      if(SettingsInfo[["Conditions"]] %in% colnames(SettingsFile_Sample)== FALSE){
        stop("The ", SettingsInfo[["Conditions"]], " column selected as Conditions in SettingsInfo was not found in SettingsFile. Please check your input.")
      }
    }

    #Biological replicates
    if("Biological_Replicates" %in% names(SettingsInfo)){
      if(SettingsInfo[["Biological_Replicates"]] %in% colnames(SettingsFile_Sample)== FALSE){
        stop("The ",SettingsInfo[["Biological_Replicates"]], " column selected as Biological_Replicates in SettingsInfo was not found in SettingsFile_Sample. Please check your input.")
      }
    }

    #Numerator
    if("Numerator" %in% names(SettingsInfo)==TRUE){
      if(SettingsInfo[["Numerator"]] %in% SettingsFile_Sample[[SettingsInfo[["Conditions"]]]]==FALSE){
        stop("The ",SettingsInfo[["Numerator"]], " column selected as numerator in SettingsInfo was not found in SettingsFile_Sample. Please check your input.")
      }
    }

   #Denominator
    if("Denominator" %in% names(SettingsInfo)==TRUE){
      if(SettingsInfo[["Denominator"]] %in% SettingsFile_Sample[[SettingsInfo[["Conditions"]]]]==FALSE){
        stop("The ",SettingsInfo[["Denominator"]], " column selected as denominator in SettingsInfo was not found in SettingsFile_Sample. Please check your input.")
      }
    }

    #Denominator & Numerator
    if("Denominator" %in% names(SettingsInfo)==FALSE  & "Numerator" %in% names(SettingsInfo) ==TRUE){
      stop("Check input. The selected denominator option is empty while ",paste(SettingsInfo[["Numerator"]])," has been selected as a numerator. Please add a denominator for 1-vs-1 comparison or remove the numerator for all-vs-all comparison." )
    }

    #Superplot
    if("Superplot" %in% names(SettingsInfo)){
      if(SettingsInfo[["Superplot"]] %in% colnames(SettingsFile_Sample)== FALSE){
        stop("The ",SettingsInfo[["Superplot"]], " column selected as Superplot column in SettingsInfo was not found in SettingsFile_Sample. Please check your input.")
      }
    }

    if(is.null(PlotSettings)==FALSE){
      if(PlotSettings== "Sample"){
        #Plot colour
        if("color" %in% names(SettingsInfo)){
          if(SettingsInfo[["color"]] %in% colnames(SettingsFile_Sample)== FALSE){
            stop("The ",SettingsInfo[["color"]], " column selected as color in SettingsInfo was not found in SettingsFile_Sample. Please check your input.")
          }
        }

        #Plot shape
        if("shape" %in% names(SettingsInfo)){
          if(SettingsInfo[["shape"]] %in% colnames(SettingsFile_Sample)== FALSE){
            stop("The ",SettingsInfo[["shape"]], " column selected as shape in SettingsInfo was not found in SettingsFile_Sample. Please check your input.")
          }
        }

        #Plot individual
        if("individual" %in% names(SettingsInfo)){
          if(SettingsInfo[["individual"]] %in% colnames(SettingsFile_Sample)== FALSE){
            stop("The ",SettingsInfo[["individual"]], " column selected as individual in SettingsInfo was not found in SettingsFile_Sample. Please check your input.")
          }
        }
      }else if(PlotSettings== "Feature"){
        if("color" %in% names(SettingsInfo)){
          if(SettingsInfo[["color"]] %in% colnames(SettingsFile_Metab)== FALSE){
            stop("The ",SettingsInfo[["color"]], " column selected as color in SettingsInfo was not found in SettingsFile_Metab. Please check your input.")
          }
        }

        #Plot shape
        if("shape" %in% names(SettingsInfo)){
          if(SettingsInfo[["shape"]] %in% colnames(SettingsFile_Metab)== FALSE){
            stop("The ",SettingsInfo[["shape"]], " column selected as shape in SettingsInfo was not found in SettingsFile_Metab. Please check your input.")
          }
        }

        #Plot individual
        if("individual" %in% names(SettingsInfo)){
          if(SettingsInfo[["individual"]] %in% colnames(SettingsFile_Metab)== FALSE){
            stop("The ",SettingsInfo[["individual"]], " column selected as individual in SettingsInfo was not found in SettingsFile_Metab. Please check your input.")
          }
        }
      }else if(PlotSettings== "Both"){
        #Plot colour sample
        if("color_Sample" %in% names(SettingsInfo)){
          if(SettingsInfo[["color_Sample"]] %in% colnames(SettingsFile_Sample)== FALSE){
            stop("The ",SettingsInfo[["color_Sample"]], " column selected as color_Sample in SettingsInfo was not found in SettingsFile_Sample. Please check your input.")
          }
        }

        #Plot colour Metab
        if("color_Metab" %in% names(SettingsInfo)){
          if(SettingsInfo[["color_Metab"]] %in% colnames(SettingsFile_Metab)== FALSE){
            stop("The ",SettingsInfo[["color_Metab"]], " column selected as color_Metab in SettingsInfo was not found in SettingsFile_Metab. Please check your input.")
          }
          if(sum(colnames(InputData) %in% SettingsFile_Metab$Metabolite) < length(InputData)  ){
            warning("The InputData contains metabolites not found in SettingsFile_Metab.")
          }
        }

       # Plot shape_metab
        if("shape_Metab" %in% names(SettingsInfo)){
          if(SettingsInfo[["shape_Metab"]] %in% colnames(SettingsFile_Metab)== FALSE){
            stop("The ",SettingsInfo[["shape_Metab"]], " column selected as shape_Metab in SettingsInfo was not found in SettingsFile_Metab. Please check your input.")
          }
        }

        # Plot shape_metab
        if("shape_Sample" %in% names(SettingsInfo)){
          if(SettingsInfo[["shape_Sample"]] %in% colnames(SettingsFile_Metab)== FALSE){
            stop("The ",SettingsInfo[["shape_Sample"]], " column selected as shape_Metab in SettingsInfo was not found in SettingsFile_Sample. Please check your input.")
          }
        }

        #Plot individual_Metab
        if("individual_Metab" %in% names(SettingsInfo)){
          if(SettingsInfo[["individual_Metab"]] %in% colnames(SettingsFile_Metab)== FALSE){
            stop("The ",SettingsInfo[["individual_Metab"]], " column selected as individual_Metab in SettingsInfo was not found in SettingsFile_Metab. Please check your input.")
          }
        }

        #Plot individual_Sample
        if("individual_Sample" %in% names(SettingsInfo)){
          if(SettingsInfo[["individual_Sample"]] %in% colnames(SettingsFile_Sample)== FALSE){
            stop("The ",SettingsInfo[["individual_Sample"]], " column selected as individual_Sample in SettingsInfo was not found in SettingsFile_Sample. Please check your input.")
          }
        }

      }

    }
   }

  #-------------SaveAs
  Save_as_Plot_options <- c("svg","pdf", "png")
  if(is.null(SaveAs_Plot)==FALSE){
    if(SaveAs_Plot %in% Save_as_Plot_options == FALSE){
    stop("Check input. The selected SaveAs_Plot option is not valid. Please select one of the folowwing: ",paste(Save_as_Plot_options,collapse = ", ")," or set to NULL if no plots should be saved." )
  }
  }


  SaveAs_Table_options <- c("txt","csv", "xlsx", "RData")#RData = SummarizedExperiment (?)
  if(is.null(SaveAs_Table)==FALSE){
    if((SaveAs_Table %in% SaveAs_Table_options == FALSE)| (is.null(SaveAs_Table)==TRUE)){
      stop("Check input. The selected SaveAs_Table option is not valid. Please select one of the folowwing: ",paste(SaveAs_Table_options,collapse = ", "),"." )
    }
  }

  #-------------CoRe
  if(is.logical(CoRe) == FALSE){
    stop("Check input. The CoRe value should be either =TRUE for preprocessing of Consuption/Release experiment or =FALSE if not.")
  }

  #-------------Theme
  if(is.null(Theme)==FALSE){
    Theme_options <- c("theme_grey()", "theme_gray()", "theme_bw()", "theme_linedraw()", "theme_light()", "theme_dark()", "theme_minimal()", "theme_classic()", "theme_void()", "theme_test()")
    if (Theme %in% Theme_options == FALSE){
      stop("Theme option is incorrect. You can check for complete themes here: https://ggplot2.tidyverse.org/reference/ggtheme.html. Options are the following: ",paste(Theme_options, collapse = ", "),"." )
    }
  }
  #------------- general
  if(is.logical(PrintPlot) == FALSE){
    stop("Check input. PrintPlot should be either =TRUE or =FALSE.")
  }
}

################################################################################################
### ### ### PreProcessing helper function: Internal Function to check function input ### ### ###
################################################################################################

#' Check specific input parameters for PreProcessing()
#'
#' @param InputData Passed to main function PreProcessing()
#' @param SettingsFile_Sample Passed to main function PreProcessing()
#' @param SettingsInfo Passed to main function PreProcessing()
#' @param CoRe Passed to main function PreProcessing()
#' @param FeatureFilt Passed to main function PreProcessing()
#' @param FeatureFilt_Value Passed to main function PreProcessing()
#' @param TIC Passed to main function PreProcessing()
#' @param MVI Passed to main function PreProcessing()
#' @param MVI_Percentage Passed to main function PreProcessing()
#' @param HotellinsConfidence Passed to main function PreProcessing()
#'
#' @keywords Input check
#' @noRd
#'
#'

CheckInput_PreProcessing <- function(InputData,
                                     SettingsFile_Sample,
                                     SettingsInfo,
                                     CoRe,
                                     FeatureFilt,
                                     FeatureFilt_Value,
                                     TIC,
                                     MVI,
                                     MVI_Percentage,
                                     HotellinsConfidence){
  if(is.vector(SettingsInfo)==TRUE){
    #-------------SettingsInfo
    #CoRe
    if(CoRe == TRUE){   # parse CoRe normalisation factor
      message("For Consumption Release experiment we are using the method from Jain M.  REF: Jain et. al, (2012), Science 336(6084):1040-4, doi: 10.1126/science.1218595.")
      if("CoRe_media" %in% names(SettingsInfo)){
        if(length(grep(SettingsInfo[["CoRe_media"]], SettingsFile_Sample[[SettingsInfo[["Conditions"]]]])) < 1){     # Check for CoRe_media samples
          stop("No CoRe_media samples were provided in the 'Conditions' in the SettingsFile_Sample. For a CoRe experiment control media samples without cells have to be measured and be added in the 'Conditions'
           column labeled as 'CoRe_media' (see @param section). Please make sure that you used the correct labelling or whether you need CoRe = FALSE for your analysis")
        }
      }

      if ("CoRe_norm_factor" %in% names(SettingsInfo)==FALSE){
        warning("No growth rate or growth factor provided for normalising the CoRe result, hence CoRe_norm_factor set to 1 for each sample")
      }
    }
  }

  #-------------General parameters
  Feature_Filtering_options <- c("Standard","Modified")
  if(FeatureFilt %in% Feature_Filtering_options == FALSE & is.null(FeatureFilt)==FALSE){
    stop("Check input. The selected FeatureFilt option is not valid. Please set to NULL or select one of the folowwing: ",paste(Feature_Filtering_options,collapse = ", "),"." )
  }
  if(is.numeric(FeatureFilt_Value) == FALSE |FeatureFilt_Value > 1 | FeatureFilt_Value < 0){
    stop("Check input. The selected FeatureFilt_Value should be numeric and between 0 and 1.")
  }
  if(is.logical(TIC) == FALSE){
    stop("Check input. The TIC should be either `TRUE` if TIC normalization is to be performed or `FALSE` if no data normalization is to be applied.")
  }
  if(is.logical(MVI) == FALSE){
    stop("Check input. MVI value should be either `TRUE` if mising value imputation should be performed or `FALSE` if not.")
  }
  if(is.numeric(MVI_Percentage)== FALSE |HotellinsConfidence > 100 | HotellinsConfidence < 0){
    stop("Check input. The selected MVI_Percentage value should be numeric and between 0 and 100.")
  }
  if( is.numeric(HotellinsConfidence)== FALSE |HotellinsConfidence > 1 | HotellinsConfidence < 0){
    stop("Check input. The selected Filtering value should be numeric and between 0 and 1.")
  }

}


################################################################################################
### ### ### DMA helper function: Internal Function to check function input ### ### ###
################################################################################################

#' Check input parameters
#'
#' @param InputData Passed to main function PreProcessing()
#' @param SettingsFile_Sample Passed to main function PreProcessing()
#' @param SettingsInfo Passed to main function PreProcessing()
#' @param CoRe Passed to main function PreProcessing()
#' @param FeatureFilt Passed to main function PreProcessing()
#' @param FeatureFilt_Value Passed to main function PreProcessing()
#' @param TIC Passed to main function PreProcessing()
#' @param MVI Passed to main function PreProcessing()
#' @param MVI_Percentage Passed to main function PreProcessing()
#' @param HotellinsConfidence Passed to main function PreProcessing()
#'
#' @keywords Input check
#' @noRd
#'
#'

CheckInput_DMA <- function(InputData,
                           SettingsFile_Sample,
                           SettingsInfo,
                           StatPval,
                           StatPadj,
                           PerformShapiro,
                           PerformBartlett,
                           Transform){

  #-------------SettingsInfo
  if(is.null(SettingsInfo)==TRUE){
    stop("You have to provide a SettingsInfo for Conditions.") # If Numerator and/or Denominator = NULL, they are not in SettingsInfo!
  }

  ## ------------ Denominator/numerator ----------- ##
  # Denominator and numerator: Define if we compare one_vs_one, one_vs_all or all_vs_all.
  if("Denominator" %in% names(SettingsInfo)==FALSE  & "Numerator" %in% names(SettingsInfo) ==FALSE){
    # all-vs-all: Generate all pairwise combinations
    conditions = SettingsFile_Sample[[SettingsInfo[["Conditions"]]]]
    denominator <-unique(SettingsFile_Sample[[SettingsInfo[["Conditions"]]]])
    numerator <-unique(SettingsFile_Sample[[SettingsInfo[["Conditions"]]]])
    comparisons <- combn(unique(conditions), 2) %>% as.matrix()
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
      stop("Check input. The selected StatPval option for Hypothesis testing is not valid for multiple comparison (one-vs-all or all-vs-all). Please select one of the following: ",paste(STAT_pval_options,collapse = ", ")," or specify numerator and denumerator." )
    }
  }else{
    STAT_pval_options <- c("aov", "kruskal.test", "welch" ,"lmFit")
    if(StatPval %in% STAT_pval_options == FALSE){
      stop("Check input. The selected StatPval option for Hypothesis testing is not valid for one-vs-one comparsion. Multiple comparison is selected. Please select one of the following: ",paste(STAT_pval_options,collapse = ", ")," or change numerator and denumerator." )
    }
  }

  STAT_padj_options <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
  if(StatPadj %in% STAT_padj_options == FALSE){
    stop("Check input. The selected StatPadj option for multiple Hypothesis testing correction is not valid. Please select one of the folowing: ",paste(STAT_padj_options,collapse = ", "),"." )
  }

  ## ------------ Sample Numbers ----------- ##
  Num <- InputData %>%#Are sample numbers enough?
    filter(SettingsFile_Sample[[SettingsInfo[["Conditions"]]]] %in% numerator) %>%
    select_if(is.numeric)#only keep numeric columns with metabolite values
  Denom <- InputData %>%
    filter(SettingsFile_Sample[[SettingsInfo[["Conditions"]]]] %in% denominator) %>%
    select_if(is.numeric)

  if(nrow(Num)==1){
    stop("There is only one sample available for ", numerator, ", so no statistical test can be performed.")
  } else if(nrow(Denom)==1){
    stop("There is only one sample available for ", denominator, ", so no statistical test can be performed.")
  }else if(nrow(Num)==0){
    stop("There is no sample available for ", numerator, ".")
  }else if(nrow(Denom)==0){
    stop("There is no sample available for ", denominator, ".")
  }

  ## ------------ Check Missingness ------------- ##
  Num_Miss <- replace(Num, Num==0, NA)
  Num_Miss <- Num_Miss[, (colSums(is.na(Num_Miss)) > 0), drop = FALSE]

  Denom_Miss <- replace(Denom, Denom==0, NA)
  Denom_Miss <- Denom_Miss[, (colSums(is.na(Denom_Miss)) > 0), drop = FALSE]

  if((ncol(Num_Miss)>0 & ncol(Denom_Miss)==0)){
    Metabolites_Miss <- colnames(Num_Miss)
    if(ncol(Num_Miss)<=10){
      message("In `Numerator` ",paste0(toString(numerator)), ", NA/0 values exist in ", ncol(Num_Miss), " Metabolite(s): ", paste0(colnames(Num_Miss), collapse = ", "), ". Those metabolite(s) might return p.val= NA, p.adj.= NA, t.val= NA. The Log2FC = Inf, if all replicates are 0/NA.")
    }else{
      message("In `Numerator` ",paste0(toString(numerator)), ", NA/0 values exist in ", ncol(Num_Miss), " Metabolite(s).", " Those metabolite(s) mightl return p.val= NA, p.adj.= NA, t.val= NA. The Log2FC = Inf, if all replicates are 0/NA.")
    }
  } else if(ncol(Num_Miss)==0 & ncol(Denom_Miss)>0){
    Metabolites_Miss <- colnames(Denom_Miss)
    if(ncol(Num_Miss)<=10){
      message("In `Denominator` ",paste0(toString(denominator)), ", NA/0 values exist in ", ncol(Denom_Miss), " Metabolite(s): ", paste0(colnames(Denom_Miss), collapse = ", "), ". Those metabolite(s) might return p.val= NA, p.adj.= NA, t.val= NA. The Log2FC = Inf, if all replicates are 0/NA.")
    }else{
      message("In `Denominator` ",paste0(toString(denominator)), ", NA/0 values exist in ", ncol(Denom_Miss), " Metabolite(s).", " Those metabolite(s) might return p.val= NA, p.adj.= NA, t.val= NA. The Log2FC = Inf, if all replicates are 0/NA.")#
    }
  } else if(ncol(Num_Miss)>0 & ncol(Denom_Miss)>0){
    Metabolites_Miss <- c(colnames(Num_Miss), colnames(Denom_Miss))
    Metabolites_Miss <- unique(Metabolites_Miss)
    message("In `Numerator` ",paste0(toString(numerator)), ", NA/0 values exist in ", ncol(Num_Miss), " Metabolite(s).", " and in `denominator`",paste0(toString(denominator)), " ",ncol(Denom_Miss), " Metabolite(s).",". Those metabolite(s) might return p.val= NA, p.adj.= NA, t.val= NA. The Log2FC = Inf, if all replicates are 0/NA.")
  } else{
    message("There are no NA/0 values")
    Metabolites_Miss <- c(colnames(Num_Miss), colnames(Denom_Miss))
    Metabolites_Miss <- unique(Metabolites_Miss)
  }

  #-------------General parameters
  if(is.logical(PerformShapiro) == FALSE){
    stop("Check input. The Shapiro value should be either =TRUE or =FALSE.")
  }
  if(is.logical(PerformBartlett) == FALSE){
    stop("Check input. The Bartlett value should be either =TRUE or =FALSE.")
  }
  if(is.logical(Transform) == FALSE){
    stop("Check input. `Transform` should be either =TRUE or =FALSE.")
  }

  ## -------- Return settings ---------##
  Settings <- list("comparisons"=comparisons, "MultipleComparison"=MultipleComparison, "all_vs_all"=all_vs_all, "Metabolites_Miss"=Metabolites_Miss, "denominator"=denominator, "numerator"=numerator)
  return(invisible(Settings))
}


################################################################################################
### ### ### ORA helper function: Internal Function to check function input ### ### ###
################################################################################################

#' Check input parameters
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

