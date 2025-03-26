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
CheckInput <- function(
    se,
    ##InputData,
    InputData_Num = TRUE,
    ##SettingsFile_Sample = NULL,
    ##SettingsFile_Metab = NULL,
    SettingsInfo = NULL,
    SaveAs_Plot = NULL,
    SaveAs_Table = NULL,
    CoRe = FALSE,
    PrintPlot = FALSE,
    Theme = NULL,
    PlotSettings = NULL) {
  
    ############## Parameters valid for multiple MetaProViz functions
    
    ## obtain colnames and rownames from se
    cols_se <- colnames(se)
    rows_se <- rownames(se)
    cols_a <- colnames(assay(se))
    rows_a <- rownames(assay(se))
    cols_cD <- colnames(colData(se))
    rows_cD <- rownames(colData(se))
    cols_rD <- colnames(rowData(se))
    rows_rD <- rownames(rowData(se))

    #-------------InputData
    if (!is(se, "SummarizedExperiment")) {
        message <- paste0("InputData should be a SummarizedExperiment object. It is currently a ", 
                          paste(class(se)), ".")
        logger::log_trace(paste0("Error ", message))
        stop(message)
    }
    # if (!is.data.frame(InputData)) {
    #     message <- paste0("InputData should be a data.frame. It's currently a ", 
    #         paste(class(InputData)), ".")
    #     logger::log_trace(paste0("Error ", message))
    #     stop(message)
    # }
    if (any(duplicated(rows_se))) {
        message <- paste0("Duplicated rownames of se, whilst rownames must be unique")
        logger::log_trace(paste0("Error ", message))
        stop(message)
    }
    # if (any(duplicated(row.names(InputData)))) {
    #     message <- paste0("Duplicated row.names of InputData, whilst row.names must be unique")
    #     logger::log_trace(paste0("Error ", message))
    #     stop(message)
    # }

    if (InputData_Num) {
        Test_num <- apply(assay(se), 2, function(x) is.numeric(x))
        ##Test_num <- apply(InputData, 2, function(x) is.numeric(x))
        if (!any(Test_num)) {
            message <- paste0("InputData needs to be of class numeric")
            logger::log_trace(paste0("Error ", message))
            stop(message)
        }
    }

    if (any(duplicated(cols_se))) {
        message <- paste0("se contains duplicates column names, whilst colnames must be unique.")
        logger::log_trace(paste0("Error ", message))
        stop(message)
    }
    # if (sum(duplicated(colnames(InputData))) > 0) {
    #     message <- paste0("InputData contained duplicates column names, whilst col.names must be unique.")
    #     logger::log_trace(paste0("Error ", message))
    #     stop(message)
    # }

    #-------------SettingsFile
    if (!all(cols_se == rows_cD)) {
        message <- paste0("colnames assay(se) need to match rownames colData(se).")
        logger::log_trace(paste0("Error ", message))
        stop(message)
    }
    
    # if (!is.null(SettingsFile_Sample)) {
    #     Test_match <- merge(SettingsFile_Sample, InputData, 
    #         by = "row.names", all =  FALSE)
    #     if (nrow(Test_match) ==  0) {
    #         message <- paste0("row.names InputData need to match row.names SettingsFile_Sample.")
    #         logger::log_trace(paste0("Error ", message))
    #         stop(message)
    #     }
    # }

    if (!all(rows_a == rows_rD)) {
        message <- paste0("rownames assay(se) need to match rownames rowData(se).")
        logger::log_trace(paste0("Error ", message))
        stop(message)
    }
    # if (!is.null(SettingsFile_Metab)) {
    #     Test_match <- merge(SettingsFile_Metab, as.data.frame(t(InputData)), 
    #         by = "row.names", all =  FALSE)
    #     if (nrow(Test_match) ==  0) {
    #         stop("col.names InputData need to match row.names SettingsFile_Metab.")
    #     }
    # }

    #-------------SettingsInfo
    if (!is.vector(SettingsInfo) & !is.null(SettingsInfo)) {
        message <- paste0("SettingsInfo should be NULL or a vector. It's currently a ", 
            paste(class(SettingsInfo), "."))
        logger::log_trace(paste0("Error ", message))
        stop(message)
    }

    if (!is.null(SettingsInfo)) {
        
        ## Conditions
        if ("Conditions" %in% names(SettingsInfo)) {
            if (!SettingsInfo[["Conditions"]] %in% cols_cD) {
                message <- paste0("The ", SettingsInfo[["Conditions"]], 
                    " column selected as Conditions in SettingsInfo was not found in colData(se). Please check your input.")
                logger::log_trace(paste0("Error ", message))
                stop(message)
            }
        }

        ## Biological replicates
        if ("Biological_Replicates" %in% names(SettingsInfo)) {
            if (!SettingsInfo[["Biological_Replicates"]] %in% cols_cD) {
                message <- paste0("The ", SettingsInfo[["Biological_Replicates"]], 
                    " column selected as Biological_Replicates in SettingsInfo was not found in colData(se). Please check your input.")
                logger::log_trace(paste0("Error ", message))
                stop(message)
            }
        }

        ## Numerator
        if ("Numerator" %in% names(SettingsInfo)) {
            if (!SettingsInfo[["Numerator"]] %in% colData(se)[[SettingsInfo[["Conditions"]]]]) {
                message <- paste0("The ", SettingsInfo[["Numerator"]], 
                    " column selected as numerator in SettingsInfo was not found in colData(se). Please check your input.")
                logger::log_trace(paste0("Error ", message))
                stop(message)
            }
        }

        ## Denominator
        if ("Denominator" %in% names(SettingsInfo)) {
            if (!SettingsInfo[["Denominator"]] %in% colData(se)[[SettingsInfo[["Conditions"]]]]) {
                message <- paste0("The ", SettingsInfo[["Denominator"]], 
                    " column selected as denominator in SettingsInfo was not found in colData(se). Please check your input.")
                logger::log_trace(paste0("Error ", message))
                stop(message)
            }
        }

        ## Denominator & Numerator
        if (!("Denominator" %in% names(SettingsInfo))  & "Numerator" %in% names(SettingsInfo)) {
            message <- paste0("Check input. The selected denominator option is empty while ",
                paste(SettingsInfo[["Numerator"]]),
                    " has been selected as a numerator. Please add a denominator for 1-vs-1 comparison or remove the numerator for all-vs-all comparison.")
            logger::log_trace(paste0("Error ", message))
            stop(message)
        }

        ## Superplot
        if ("Superplot" %in% names(SettingsInfo)) {
            if (!SettingsInfo[["Superplot"]] %in% cols_cD) {
                message <- paste0("The ", SettingsInfo[["Superplot"]], 
                    " column selected as Superplot column in SettingsInfo was not found in colData(se). Please check your input.")
                logger::log_trace(paste0("Error ", message))
                stop(message)
            }
        }

        if (!is.null(PlotSettings)) {
            
            if (PlotSettings == "Sample") {
                
                ## plot colour
                if ("color" %in% names(SettingsInfo)) {
                    if (!SettingsInfo[["color"]] %in% cols_cD) {
                        message <- paste0("The ", SettingsInfo[["color"]], 
                            " column selected as color in SettingsInfo was not found in colData(se). Please check your input.")
                        logger::log_trace(paste0("Error ", message))
                        stop(message)
                    }
                }

                ## plot shape
                if ("shape" %in% names(SettingsInfo)) {
                    if (!SettingsInfo[["shape"]] %in% cols_cD) {
                        message <- paste0("The ", SettingsInfo[["shape"]], 
                            " column selected as shape in SettingsInfo was not found in colData(se). Please check your input.")
                        logger::log_trace(paste0("Error ", message))
                        stop(message)
                    }
                }

                ## plot individual
                if ("individual" %in% names(SettingsInfo)) {
                    if (!SettingsInfo[["individual"]] %in% cols_cD) {
                        message <- paste0("The ", SettingsInfo[["individual"]], " column selected as individual in SettingsInfo was not found in colData(se). Please check your input.")
                        logger::log_trace(paste0("Error ", message))
                        stop(message)
                    }
                } 
            } else if (PlotSettings == "Feature") {
                
                ## plot color
                if ("color" %in% names(SettingsInfo)) {
                    if (!SettingsInfo[["color"]] %in% cols_rD) {
                        message <- paste0("The ", SettingsInfo[["color"]], 
                            " column selected as color in SettingsInfo was not found in rowData(se). Please check your input.")
                        logger::log_trace(paste0("Error ", message))
                        stop(message)
                    }
                }

                ## plot shape
                if ("shape" %in% names(SettingsInfo)) {
                    if (!SettingsInfo[["shape"]] %in% cols_rD) {
                        message <- paste0("The ", SettingsInfo[["shape"]], 
                            " column selected as shape in SettingsInfo was not found in rowData(se). Please check your input.")
                        logger::log_trace(paste0("Error ", message))
                        stop(message)
                    }
                }

                ## Plot individual
                if ("individual" %in% names(SettingsInfo)) {
                    if (!SettingsInfo[["individual"]] %in% cols_rD) {
                        message <- paste0("The ", SettingsInfo[["individual"]], 
                            " column selected as individual in SettingsInfo was not found in rowData(se). Please check your input.")
                        logger::log_trace(paste0("Error ", message))
                        stop(message)
                    }
                }
            } else if (PlotSettings == "Both") {
        
                # plot colour sample
                if ("color_Sample" %in% names(SettingsInfo)) {
                    if (!SettingsInfo[["color_Sample"]] %in% cols_rD) {
                        message <- paste0("The ", SettingsInfo[["color_Sample"]], 
                            " column selected as color_Sample in SettingsInfo was not found in rowData(se). Please check your input.")
                        logger::log_trace(paste0("Error ", message))
                        stop(message)
                    }
                }

                ## plot colour Metab
                if ("color_Metab" %in% names(SettingsInfo)) {
                    if (!SettingsInfo[["color_Metab"]] %in% cols_rD) {
                        message <- paste0("The ", SettingsInfo[["color_Metab"]], 
                            " column selected as color_Metab in SettingsInfo was not found in rowData(se). Please check your input.")
                        logger::log_trace(paste0("Error ", message))
                        stop(message)
                    }
                    if (!all(rows_a == rows_rD)) {
                        message <- paste0("assay(se) has to contain the same metabolites as in rowData(se).")
                        logger::log_trace(paste0("Warning ", message))
                        warning(message)
                    }
                }

                ## plot shape_metab
                if ("shape_Metab" %in% names(SettingsInfo)) {
                    if (!SettingsInfo[["shape_Metab"]] %in% cols_rD) {
                        message <- paste0("The ", SettingsInfo[["shape_Metab"]], 
                            " column selected as shape_Metab in SettingsInfo was not found in rowData(se). Please check your input.")
                        logger::log_trace(paste0("Error ", message))
                        stop(message)
                    }
                }

                ## plot shape_metab
                if ("shape_Sample" %in% names(SettingsInfo)) {
                    if (!SettingsInfo[["shape_Sample"]] %in% cols_rD) {
                        message <- paste0("The ", SettingsInfo[["shape_Sample"]], 
                            " column selected as shape_Metab in SettingsInfo was not found in rowData(se). Please check your input.")
                        logger::log_trace(paste0("Error ", message))
                        stop(message)
                    }
                }

                ## plot individual_Metab
                if ("individual_Metab" %in% names(SettingsInfo)) {
                    if (!SettingsInfo[["individual_Metab"]] %in% cols_rD) {
                        message <- paste0("The ", SettingsInfo[["individual_Metab"]], 
                            " column selected as individual_Metab in SettingsInfo was not found in rowData(se). Please check your input.")
                        logger::log_trace(paste0("Error ", message))
                        stop(message)
                    }
                }

                ## plot individual_Sample
                if ("individual_Sample" %in% names(SettingsInfo)) {
                    if (!SettingsInfo[["individual_Sample"]] %in% cols_cD) {
                        message <- paste0("The ", SettingsInfo[["individual_Sample"]], 
                            " column selected as individual_Sample in SettingsInfo was not found in colData(se). Please check your input.")
                        logger::log_trace(paste0("Error ", message))
                        stop(message)
                    }
                }
            }
        }
    }

    #-------------SaveAs
    Save_as_Plot_options <- c("svg","pdf", "png") ## EDIT: this should be rewritten with match.arg
    if (!is.null(SaveAs_Plot)) {
        if (!SaveAs_Plot %in% Save_as_Plot_options) {
            message <- paste0("Check input. The selected SaveAs_Plot option is not valid. Please select one of the following: ", 
                paste(Save_as_Plot_options, collapse = ", "), 
                " or set to NULL if no plots should be saved.")
            logger::log_trace(paste0("Error ", message))
            stop(message)
        }
    }

    SaveAs_Table_options <- c("txt","csv", "xlsx", "RData")#RData = SummarizedExperiment (?) ## EDIT: this should be rewritten with match.arg
    if (!is.null(SaveAs_Table)) {
        if (!(SaveAs_Table %in% SaveAs_Table_options) | is.null(SaveAs_Table)) {
            message <- paste0("Check input. The selected SaveAs_Table option is not valid. Please select one of the following: ",
                paste(SaveAs_Table_options,collapse = ", "), 
                " or set to NULL if no tables should be saved.")
            logger::log_trace(paste0("Error ", message))
            stop(message)
        }
    }

    #-------------CoRe
    if (!is.logical(CoRe)) {
        message <- paste0("Check input. The CoRe value should be either TRUE for preprocessing of Consumption/Release experiment or FALSE if not.")
        logger::log_trace(paste0("Error ", message))
        stop(message)
    }

    #-------------Theme
    if (!is.null(Theme)) { ## EDIT: this should be rewritten with match.arg
        Theme_options <- c("theme_grey()", "theme_gray()", "theme_bw()", "theme_linedraw()", "theme_light()", "theme_dark()", "theme_minimal()", "theme_classic()", "theme_void()", "theme_test()")
        if (!Theme %in% Theme_options) {
            message <- paste0("Check input. Theme option is incorrect. You can check for complete themes here: https://ggplot2.tidyverse.org/reference/ggtheme.html. Options are the following: ",
                paste(Theme_options, collapse = ", "), "." )
            logger::log_trace(paste0("Error ", message))
            stop(message)
        }
    }
  
    #------------- general
    if (!is.logical(PrintPlot)) {
        message <- paste0("Check input. PrintPlot should be either TRUE or FALSE.")
        logger::log_trace(paste0("Error ", message))
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
CheckInput_PreProcessing <- function(se,
    SettingsInfo,
    CoRe = FALSE,
    FeatureFilt = "Modified",
    FeatureFilt_Value = 0.8,
    TIC = TRUE,
    MVI = TRUE,
    MVI_Percentage = 50,
    HotellinsConfidence = 0.99) {
  
    if (is.vector(SettingsInfo)) {
    
        #-------------SettingsInfo
        ## CoRe
        if (CoRe) {
            ## parse CoRe normalisation factor
            message <- paste0("For Consumption Release experiment we are using the method from Jain M.  REF: Jain et. al, (2012), Science 336(6084):1040-4, doi: 10.1126/science.1218595.")
            logger::log_trace(paste0("Message ", message))
            message(message)
            
            if ("CoRe_media" %in% names(SettingsInfo)) {
                if (length(grep(SettingsInfo[["CoRe_media"]], colData(se)[[SettingsInfo[["Conditions"]]]])) < 1) {     
                    ## check for CoRe_media samples
                    message <- paste0("No CoRe_media samples were provided in the 'Conditions' in colData(se). ",
                        "For a CoRe experiment control media samples without cells have to be measured and be added ",
                        "in the 'Conditions' column labeled as 'CoRe_media' (see @param section). Please make sure ",
                        "that you used the correct labelling or whether you need CoRe = FALSE for your analysis")
                    logger::log_trace(paste0("Error ", message))
                    stop(message)
                }
            }

            if (!"CoRe_norm_factor" %in% names(SettingsInfo)) {
                message <- paste0("No growth rate or growth factor provided for ",
                    "normalising the CoRe result, hence CoRe_norm_factor set to ",
                    "1 for each sample")
                logger::log_trace(paste0("Warning ", message))
                warning(message)
            }
        }
    }

    #-------------General parameters
    Feature_Filtering_options <- c("Standard", "Modified")
    if (!FeatureFilt %in% Feature_Filtering_options & !is.null(FeatureFilt)) {
        message <- paste0("Check input. The selected FeatureFilt option is not valid. ",
            "Please set to NULL or select one of the folowwing: ",
            paste(Feature_Filtering_options,collapse = ", "), "." )
        logger::log_trace(paste0("Error ", message))
        stop(message)
    }
    
    if (!is.numeric(FeatureFilt_Value) |FeatureFilt_Value > 1 | FeatureFilt_Value < 0) {
        message <- paste0("Check input. The selected FeatureFilt_Value should be numeric and between 0 and 1.")
        logger::log_trace(paste0("Error ", message))
        stop(message)
    }
    if (!is.logical(TIC)) {
        message <- paste0("Check input. The TIC value should be either `TRUE` if ",
            "TIC normalization is to be performed or `FALSE` if no data ",
            "normalization is to be applied.")
        logger::log_trace(paste0("Error ", message))
        stop(message)
    }
    if (!is.logical(MVI)) {
        message <- paste0("Check input. The MVI value should be either `TRUE` if ",
            "missing value imputation should be performed or `FALSE` if not.")
        logger::log_trace(paste0("Error ", message))
        stop(message)
    }
    if (!is.numeric(MVI_Percentage) | HotellinsConfidence > 100 | HotellinsConfidence < 0) {
        message <- paste0("Check input. The selected MVI_Percentage value should be numeric and between 0 and 100.")
        logger::log_trace(paste0("Error ", message))
        stop(message)
    }
    if (!is.numeric(HotellinsConfidence) | HotellinsConfidence > 1 | HotellinsConfidence < 0) {
        message <- paste0("Check input. The selected HotellinsConfidence value ",
            "should be numeric and between 0 and 1.")
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
#' @return Returns: 1. warnings and errors if input is not correct, 2. Settings
#'
#' @keywords Input check for MetaProViz::DMA
#'
#' @importFrom logger log_trace
#' @importFrom dplyr filter select_if
#' @importFrom magrittr %>%
#' @importFrom stats p.adjust.methods
#' @importFrom utils combn
#'
#' @noRd
#'
CheckInput_DMA <- function(
    se,
    #InputData,
    #SettingsFile_Sample,
    SettingsInfo = c(Conditions = "Conditions", Numerator = NULL, Denominator = NULL),
    StatPval = "lmFit", ## EDIT: show here the options and match with match.arg
    StatPadj = p.adjust.methods,
    VST = FALSE,
    PerformShapiro = TRUE,
    PerformBartlett = TRUE,
    Transform = TRUE) {

    ## match arguments
    StatPadj <- match.arg(StatPadj)
    
    ##-------------SettingsInfo
    if (is.null(SettingsInfo)) {
        message <- paste0("You have to provide SettingsInfo's for Conditions.")
        logger::log_trace(paste0("Error ", message))
        stop(message)
    }

    ## ------------ Denominator/numerator ----------- ##
    ## Denominator and numerator: Define if we compare one_vs_one, 
    ## one_vs_all or all_vs_all.
    if (!("Denominator" %in% names(SettingsInfo)) & !("Numerator" %in% names(SettingsInfo))) { ## EDIT: written several times in the funciton, write a function
        
        ## all-vs-all: Generate all pairwise combinations
        conditions <- colData(se)[[SettingsInfo[["Conditions"]]]]
        denominator <- unique(conditions)
        numerator <- unique(conditions)
        comparisons <- utils::combn(unique(conditions), 2) %>% 
            as.matrix()
        
        ## settings:
        MultipleComparison <- TRUE
        all_vs_all <- TRUE
        
    } else if ("Denominator" %in% names(SettingsInfo) & !("Numerator" %in% names(SettingsInfo))) {
        
        ## all-vs-one: Generate the pairwise combinations
        conditions = colData(se)[[SettingsInfo[["Conditions"]]]]
        denominator <- SettingsInfo[["Denominator"]]
        numerator <-unique(conditions)
        
        ## remove denom from num
        numerator <- numerator[!numerator %in% denominator]
        comparisons  <- t(expand.grid(numerator, denominator)) %>% 
            as.data.frame()
    
        ## settings:
        MultipleComparison <- TRUE
        all_vs_all <- FALSE
        
    } else if ("Denominator" %in% names(SettingsInfo)  & "Numerator" %in% names(SettingsInfo)) {
        
        ## one-vs-one: Generate the comparisons
        denominator <- SettingsInfo[["Denominator"]]
        numerator <- SettingsInfo[["Numerator"]]
        comparisons <- matrix(c(SettingsInfo[["Denominator"]], SettingsInfo[["Numerator"]]))
    
        ## settings:
        MultipleComparison <- FALSE
        all_vs_all <- FALSE
    }

    ## ------------ Test statistics ----------- ##
    if (!MultipleComparison) {
        STAT_pval_options <- c("t.test", "wilcox.test","chisq.test", "cor.test", "lmFit")
        
        if (!StatPval %in% STAT_pval_options) {
            message <- paste0("Check input. The selected StatPval option for ",
                "Hypothesis testing is not valid for multiple comparison ",
                "(one-vs-all or all-vs-all). Please select one of the following: ",
                paste(STAT_pval_options,collapse = ", "), 
                " or specify numerator and denumerator." )
            logger::log_trace(paste0("Error ", message))
            stop(message)
        }
    } else {
        STAT_pval_options <- c("aov", "kruskal.test", "welch" ,"lmFit")
        
        if (!StatPval %in% STAT_pval_options) {
            message <- paste0("Check input. The selected StatPval option for Hypothesis testing is not valid for one-vs-one comparsion. Multiple comparison is selected. Please select one of the following: ",paste(STAT_pval_options,collapse = ", ")," or change numerator and denumerator." )
            logger::log_trace(paste0("Error ", message))
            stop(message)
        }
    }

    STAT_padj_options <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
    if (!StatPadj %in% STAT_padj_options) {
        message <- paste0("Check input. The selected StatPadj option for ",
            "multiple Hypothesis testing correction is not valid. Please ",
            "select one of the folowing: ",
            paste(STAT_padj_options,collapse = ", "), "." )
        logger::log_trace(paste0("Error ", message))
        stop(message)
    }

    ## ------------ Sample Numbers ----------- ##
    Num <- assay(se)[, colData(se)[[SettingsInfo[["Conditions"]]]] %in% numerator]## %>%
        ## are sample numbers enough?
        ##dplyr::filter(colData(se)[[SettingsInfo[["Conditions"]]]] %in% numerator) %>%
        ## only keep numeric columns with metabolite values
        #dplyr::select_if (is.numeric)
    Denom <- assay(se)[, colData(se)[[SettingsInfo[["Conditions"]]]] %in% denominator]## %>%
        ##as.data.frame() |>
        #dplyr::select(colData(se)[[SettingsInfo[["Conditions"]]]] %in% denominator) %>%
        ##dplyr::select_if (is.numeric)

    if (ncol(Num) == 1) {
        message <- paste0("There is only one sample available for ", 
            numerator, ". No statistical test can be performed.")
        logger::log_trace(paste0("Error ", message))
        stop(message)
    } else if (ncol(Denom) == 1) {
        message <- paste0("There is only one sample available for ", 
            denominator, ". No statistical test can be performed.")
        logger::log_trace(paste0("Error ", message))
        stop(message)
    } else if (ncol(Num) == 0) {
        message <- paste0("There is no sample available for ", numerator, ".")
        logger::log_trace(paste0("Error ", message))
        stop(message)
    } else if (ncol(Denom) == 0) {
        message <- paste0("There is no sample available for ", denominator, ".")
        logger::log_trace(paste0("Error ", message))
        stop(message)
    }

    ## ------------ Check Missingness ------------- ##
    Num_Miss <- replace(Num, Num == 0, NA)
    Num_Miss <- Num_Miss[rowSums(is.na(Num_Miss)) > 0, , drop = FALSE]

    Denom_Miss <- replace(Denom, Denom == 0, NA)
    Denom_Miss <- Denom_Miss[rowSums(is.na(Denom_Miss)) > 0, , drop = FALSE]

    if (nrow(Num_Miss) > 0 & nrow(Denom_Miss) == 0) {
        Metabolites_Miss <- rownames(Num_Miss)
        if (nrow(Num_Miss) <= 10) {
            message <- paste0("In `Numerator` ", paste0(toString(numerator)), 
                ", NA/0 values exist in ", nrow(Num_Miss), " Metabolite(s): ", 
                paste0(rownames(Num_Miss), collapse = ", "), 
                ". Those metabolite(s) might return p.val= NA, p.adj.= NA, ",
                "t.val= NA. The Log2FC = Inf, if all replicates are 0/NA.")
            logger::log_info(message)
            message(message)
        } else {
            message <- paste0("In `Numerator` ", paste0(toString(numerator)), 
                ", NA/0 values exist in ", nrow(Num_Miss), " Metabolite(s).", 
                " Those metabolite(s) might return p.val= NA, p.adj.= NA, ",
                "t.val= NA. The Log2FC = Inf, if all replicates are 0/NA.")
            logger::log_info(message)
            message(message)
        }
    } else if (nrow(Num_Miss) == 0 & nrow(Denom_Miss) > 0) {
        Metabolites_Miss <- rownames(Denom_Miss)
        if (nrow(Num_Miss) <= 10) {
            message <- paste0("In `Denominator` ", paste0(toString(denominator)), 
                ", NA/0 values exist in ", nrow(Denom_Miss), " Metabolite(s): ", 
                paste0(rownames(Denom_Miss), collapse = ", "), 
                ". Those metabolite(s) might return p.val= NA, p.adj.= NA, ",
                "t.val= NA. The Log2FC = Inf, if all replicates are 0/NA.")
            logger::log_info(message)
            message(message)
        } else {
            message <- paste0("In `Denominator` ", 
                paste0(toString(denominator)), 
                ", NA/0 values exist in ", nrow(Denom_Miss), " Metabolite(s).", 
                " Those metabolite(s) might return p.val= NA, p.adj.= NA, ",
                "t.val= NA. The Log2FC = Inf, if all replicates are 0/NA.")
            logger::log_info(message)
            message(message)
        }
    } else if (nrow(Num_Miss) > 0 & nrow(Denom_Miss) > 0) {
        Metabolites_Miss <- c(rownames(Num_Miss), rownames(Denom_Miss))
        Metabolites_Miss <- unique(Metabolites_Miss)

        message <- paste0("In `Numerator` ", paste0(toString(numerator)), 
            ", NA/0 values exist in ", nrow(Num_Miss), " Metabolite(s).", 
            " and in `denominator`",paste0(toString(denominator)), " ",
            ncol(Denom_Miss), " Metabolite(s).", 
            ". Those metabolite(s) might return p.val= NA, p.adj.= NA, t.val= NA. ",
            "The Log2FC = Inf, if all replicates are 0/NA.")
        logger::log_info(message)
        message(message)
    } else {
        message <- paste0("There are no NA/0 values")
        logger::log_info(message)
        message(message)

        Metabolites_Miss <- c(rownames(Num_Miss), rownames(Denom_Miss))
        Metabolites_Miss <- unique(Metabolites_Miss)
    }

    #-------------General parameters
    if (!is.logical(VST)) {
        message <- paste0("Check input. The VST value should be either TRUE or FALSE.")
        logger::log_trace(paste0("Error ", message))
        stop(message)
    }

    if (!is.logical(PerformShapiro)) {
        message <- paste0("Check input. The Shapiro value should be either TRUE or FALSE.")
        logger::log_trace(paste0("Error ", message))
        stop(message)
    }
    if (!is.logical(PerformBartlett)) {
        message <- paste0("Check input. The Bartlett value should be either TRUE or FALSE.")
        logger::log_trace(paste0("Error ", message))
        stop(message)
    }
    if (!is.logical(Transform)) {
        message <- paste0("Check input. `Transform` should be either =TRUE or =FALSE.")
        logger::log_trace(paste0("Error ", message))
        stop(message)
    }

    Settings <- list(
        "comparisons" = comparisons, "MultipleComparison" = MultipleComparison, 
        "all_vs_all" = all_vs_all, "Metabolites_Miss" = Metabolites_Miss, 
        "denominator" = denominator, "numerator" = numerator)
    
    ## return
    invisible(Settings)
}

################################################################################################
### ### ### ORA helper function: Internal Function to check function input ### ### ###
################################################################################################

#' Check input parameters of ORA
#'
#' @param InputData DF with metabolite names/metabolite IDs as row names. Metabolite names/IDs need to match the identifier type (e.g. HMDB IDs) in the PathwayFile.
#' @param SettingsInfo \emph{Optional: } Pass ColumnName of the column including parameters to use for pCutoff and PercentageCutoff, ColumnName for PathwayFile. For MetaProViz::ClusterORA also BackgroundColumn. \strong{c(pvalColumn="p.adj", PercentageColumn="t.val", PathwayTerm= "term", PathwayFeature= "Metabolite")}
#' @param pCutoff \emph{Optional: } p-adjusted value cutoff from ORA results. Must be a numeric value. \strong{default: 0.05}
#' @param PercentageCutoff \emph{Optional: } Percentage cutoff of metabolites that should be considered for ORA. Selects Top/Bottom % of selected PercentageColumn, usually t.val or Log2FC \strong{default: 10}
#' @param PathwayFile DF that must include column "term" with the pathway name, column "Metabolite" with the Metabolite name or ID and column "Description" with pathway description that will be depicted on the plots.
#' @param PathwayName \emph{Optional: } Name of the PathwayFile used \strong{default: ""}
#' @param minGSSize \emph{Optional: } minimum group size in ORA \strong{default: 10}
#' @param maxGSSize \emph{Optional: } maximum group size in ORA \strong{default: 1000}
#' @param SaveAs_Table \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt" \strong{default: "csv"}
#' @param RemoveBackground For MetaProViz::ClusterORA the Background Settings are passed, for MetaProViz::StandardORA set to FALSE.
#'
#' @return Returns: 1. warnings and errors if input is not correct, 2. Pathway file
#'
#' @keywords Input check MetaProViz::StandardORA and MetaProViz::ClusterORA
#'
#' @importFrom logger log_trace
#' @importFrom dplyr rename
#'
#' @noRd
#'
CheckInput_ORA <- function(
    se,
    ##InputData,
    SettingsInfo = c(pvalColumn = "p.adj", PercentageColumn = "t.val", 
        PathwayTerm = "term", PathwayFeature = "Metabolite"),
    pCutoff = 0.05, 
    PercentageCutoff = 10,
    PathwayFile,
    PathwayName = "",
    minGSSize = 10,
    maxGSSize = 1000 ,
    SaveAs_Table = "csv", ## EDIT: name here the options and use match.arg
    RemoveBackground) {
    
    ## obtain colnames and rownames from assay(se)
    cols_a <- colnames(assay(se))
    rows_a <- rownames(assay(se))
    
    ## 1. The input data:
    if (is(assay(se), "matrix")) {
        message <- paste0(
            "assay(se) should be a matrix. It's currently a ", 
            paste(class(assay(se)), "."))
        logger::log_trace(paste0("Error ", message))
        stop(message)
    }
    if (any(duplicated(rows_a))) {
        message <- paste0("Duplicated rownames of assay(se), whilst rownames must be unique")
        logger::log_trace(paste0("Error ", message))
        stop(message)
    }

    # 2. Settings Columns:
    if (!is.vector(SettingsInfo) & !is.null(SettingsInfo)) {
        message <- paste0(
            "SettingsInfo should be NULL or a vector. It's currently a ", 
            paste0(class(SettingsInfo), "."))
        logger::log_trace(paste0("Error ", message))
        stop(message)
    }

    if (!is.null(SettingsInfo)) {
        ## "ClusterColumn"
        if ("ClusterColumn" %in% names(SettingsInfo)) {
            if (!SettingsInfo[["ClusterColumn"]] %in% cols_a) {
                message <- paste0("The ", SettingsInfo[["ClusterColumn"]], 
                    " column selected as ClusterColumn in SettingsInfo was not ",
                    "found in assay(se). Please check your input.")
                logger::log_trace(paste0("Error ", message))
                stop(message)
            }
        }

        ## "BackgroundColumn"
        if ("BackgroundColumn" %in% names(SettingsInfo)) {
            if (!SettingsInfo[["BackgroundColumn"]] %in% cols_a) {
                message <- paste0("The ", SettingsInfo[["BackgroundColumn"]], 
                    " column selected as BackgroundColumn in SettingsInfo was ",
                    "not found in assay(se). Please check your input.")
                logger::log_trace(paste0("Error ", message))
                stop(message)
            }
        }

        ## "pvalColumn"
        if ("pvalColumn" %in% names(SettingsInfo)) {
            if (!SettingsInfo[["pvalColumn"]] %in% cols_a) {
                message <- paste0("The ", SettingsInfo[["pvalColumn"]], 
                    " column selected as pvalColumn in SettingsInfo was not ",
                    "found in assay(se). Please check your input.")
                logger::log_trace(paste0("Error ", message))
                stop(message)
            }
        }

        ## "PercentageColumn"
        if ("PercentageColumn" %in% names(SettingsInfo)) {
            if (!SettingsInfo[["PercentageColumn"]] %in% cols_a) {
                message <- paste0("The ", SettingsInfo[["PercentageColumn"]], 
                    " column selected as PercentageColumn in SettingsInfo was ",
                    "not found in assay(se). Please check your input.")
                logger::log_trace(paste0("Error ", message))
                stop(message)
            }
        }

        ## "PathwayTerm"
        if ("PathwayTerm" %in% names(SettingsInfo)) {
            if (!SettingsInfo[["PathwayTerm"]] %in% colnames(PathwayFile)) {
                message <- paste0("The ", SettingsInfo[["PathwayTerm"]], 
                    " column selected as PathwayTerm in SettingsInfo was not ",
                    "found in PathwayFile. Please check your input.")
                logger::log_trace(paste0("Error ", message))
                stop(message)
            } else {
                PathwayFile <- PathwayFile%>%
                    dplyr::rename("term" = SettingsInfo[["PathwayTerm"]])
                PathwayFile$Description <- PathwayFile$term
            }
        } else {
            message <- paste0("SettingsInfo must provide the column name for PathwayTerm in PathwayFile")
            logger::log_trace(paste0("Error ", message))
            stop(message)
        }

        ## PathwayFeature
        if ("PathwayFeature" %in% names(SettingsInfo)) {
            if (!SettingsInfo[["PathwayFeature"]] %in% colnames(PathwayFile)) {
                message <- paste0("The ", SettingsInfo[["PathwayFeature"]], 
                    " column selected as PathwayFeature in SettingsInfo was not ",
                    "found in PathwayFile. Please check your input.")
                logger::log_trace(paste0("Error ", message))
                stop(message)
            } else {
                PathwayFile <- PathwayFile %>%
                    dplyr::rename("gene" = SettingsInfo[["PathwayFeature"]])
            }
        } else {
            message <- paste0("SettingsInfo must provide the column name for PathwayFeature in PathwayFile")
            logger::log_trace(paste0("Error ", message))
            stop(message)
        }

    } else {
        message <- paste0("You must provide SettingsInfo.")
        logger::log_trace(paste0("Error ", message))
        stop(message)
    }

    ## 3. General Settings
    if (!is.character(PathwayName)) {
        message <- paste0("Check input. PathwayName must be a character of syntax 'example'.")
        logger::log_trace(paste0("Error ", message))
        stop(message)
    }

    if (!is.logical(RemoveBackground)) {
        message <- paste0("Check input. RemoveBackground value should be either TRUE or FALSE.")
        logger::log_trace(paste0("Error ", message))
        stop(message)
    }

    if (!is.numeric(minGSSize)) {
        message <- paste0("Check input. The selected minGSSize value should be numeric.")
        logger::log_trace(paste0("Error ", message))
        stop(message)
    }

    if (!is.numeric(maxGSSize)) {
        message <- paste0("Check input. The selected maxGSSize value should be numeric.")
        logger::log_trace(paste0("Error ", message))
        stop(message)
    }

    SaveAs_Table_options <- c("txt","csv", "xlsx")
    if (!is.null(SaveAs_Table)) {
        if (!(SaveAs_Table %in% SaveAs_Table_options)| is.null(SaveAs_Table)) {
            message <- paste0("Check input. The selected SaveAs_Table option is not valid. Please select one of the folowing: ",paste(SaveAs_Table_options,collapse = ", "),"." )
            logger::log_trace(paste0("Error ", message))
            stop(message)
        }
    }

    if (!is.null(pCutoff)) {
        if (!is.numeric(pCutoff) | pCutoff > 1 |  pCutoff < 0) {
            message <- paste0("Check input. The selected pCutoff value should be numeric and between 0 and 1.")
            logger::log_trace(paste0("Error ", message))
            stop(message)
        }
    }

    if (!is.null(PercentageCutoff)) {
        if (!is.numeric(PercentageCutoff) | PercentageCutoff > 100 | PercentageCutoff < 0) {
            message <- paste0("Check input. The selected PercentageCutoff value should be numeric and between 0 and 100.")
            logger::log_trace(paste0("Error ", message))
            stop(message)
        }
    }

    ## -------- Return Pathways ---------##
    invisible(PathwayFile)
}


################################################################################################
### ### ### MCA Helper function: Internal Function to check function input ### ### ###
################################################################################################

#' Check input parameters of MCA
#'
#' @param InputData_Intra For MetaProViz::MCA_CoRe, otherwise NULL. DF for your data (results from e.g. DMA) containing metabolites in rows with corresponding Log2FC and stat (p-value, p.adjusted) value columns.
#' @param InputData_CoRe For MetaProViz::MCA_CoRe, otherwise NULL. DF for your data (results from e.g. DMA) containing metabolites in rows with corresponding Log2FC and stat (p-value, p.adjusted) value columns. Here we additionally require
#' @param SettingsInfo_Intra For MetaProViz::MCA_CoRe, otherwise NULL. Pass ColumnNames and Cutoffs for the intracellular metabolomics including the value column (e.g. Log2FC, Log2Diff, t.val, etc) and the stats column (e.g. p.adj, p.val). This must include: c(ValueCol=ColumnName_InputData_Intra,StatCol=ColumnName_InputData_Intra, StatCutoff= NumericValue, ValueCutoff=NumericValue)
#' @param SettingsInfo_CoRe  For MetaProViz::MCA_CoRe, otherwise NULL. Pass ColumnNames and Cutoffs for the consumption-release metabolomics including the direction column, the value column (e.g. Log2Diff, t.val, etc) and the stats column (e.g. p.adj, p.val). This must include: c(DirectionCol= ColumnName_InputData_CoRe,ValueCol=ColumnName_InputData_CoRe,StatCol=ColumnName_InputData_CoRe, StatCutoff= NumericValue, ValueCutoff=NumericValue)
#' @param InputData_C1 For MetaProViz::MCA_2Cond, otherwise NULL. DF for your data (results from e.g. DMA) containing metabolites in rows with corresponding Log2FC and stat (p-value, p.adjusted) value columns.
#' @param InputData_C2 For MetaProViz::MCA_2Cond, otherwise NULL. DF for your data (results from e.g. DMA) containing metabolites in rows with corresponding Log2FC and stat (p-value, p.adjusted) value columns.
#' @param SettingsInfo_C1  For MetaProViz::MCA_2Cond, otherwise NULL. Pass ColumnNames and Cutoffs for condition 1 including the value column (e.g. Log2FC, Log2Diff, t.val, etc) and the stats column (e.g. p.adj, p.val). This must include: c(ValueCol=ColumnName_InputData_C1,StatCol=ColumnName_InputData_C1, StatCutoff= NumericValue, ValueCutoff=NumericValue)
#' @param SettingsInfo_C2 For MetaProViz::MCA_2Cond, otherwise NULL. Pass ColumnNames and Cutoffs for condition 2 includingthe value column (e.g. Log2FC, Log2Diff, t.val, etc) and the stats column (e.g. p.adj, p.val). This must include: c(ValueCol=ColumnName_InputData_C2,StatCol=ColumnName_InputData_C2, StatCutoff= NumericValue, ValueCutoff=NumericValue)
#' @param FeatureID \emph{Optional: } Column name of Column including the Metabolite identifiers. This MUST BE THE SAME in each of your Input files. \strong{Default="Metabolite"}
#' @param BackgroundMethod \emph{Optional: } Background method `Intra|CoRe, Intra&CoRe, CoRe, Intra or * \strong{Default="Intra&CoRe"}
#' @param SaveAs_Table \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt" \strong{default: "csv"}
#'
#' @return Returns warnings and errors if input is not correct
#'
#' @keywords Input check MetaProViz::MCA_2Cond and MetaProViz::MCA_CoRe
#'
#' @importFrom logger log_trace
#'
#' @noRd
#'
CheckInput_MCA <- function(InputData_C1,
    InputData_C2,
    InputData_CoRe,
    InputData_Intra,
    SettingsInfo_C1,
    SettingsInfo_C2,
    SettingsInfo_CoRe,
    SettingsInfo_Intra,
    FeatureID = "Metabolite",
    BackgroundMethod = "Intra&CoRe",
    SaveAs_Table = "csv") { ## EDIT: provide options here

    ## ------------ Create log file ----------- ##
    MetaProViz_Init()

    #------------- InputData ## EDIT: this can be simplified, when you are checking for C1 and C2 for the same attributes, write a function and apply on C1 and C2 separately
    if (!is.null(InputData_C1)) {
        if (class(InputData_C1) != "data.frame"| class(InputData_C2) != "data.frame") {
            message <- paste0(
                "InputData_C1 and InputData_C2 should be a data.frame. It's currently a ", 
                paste(class(InputData_C1)), paste(class(InputData_C2)), ".")
            logger::log_trace(paste0("Error ", message))
            stop(message)
        }
        if (length(InputData_C1[duplicated(InputData_C1[[FeatureID]]), FeatureID]) > 0) {
            message <- paste0("Duplicated FeatureIDs of InputData_C1, whilst features must be unique")
            logger::log_trace(paste0("Error ", message))
            stop(message)
        }
        if (length(InputData_C2[duplicated(InputData_C2[[FeatureID]]), FeatureID]) > 0) {
            message <- paste0("Duplicated FeatureIDs of InputData_C2, whilst features must be unique")
            logger::log_trace(paste0("Error ", message))
            stop(message)
        }

    } else {
        if (class(InputData_Intra) != "data.frame"| class(InputData_CoRe) != "data.frame") {
            message <- paste0("InputData_Intra and InputData_CoRe should be a data.frame. It's currently a ", paste(class(InputData_Intra)), paste(class(InputData_CoRe)), ".",sep = "")
            logger::log_trace(paste0("Error ", message))
            stop(message)
        }
        if (length(InputData_Intra[duplicated(InputData_Intra[[FeatureID]]), FeatureID]) > 0) {
            message <- paste0("Duplicated FeatureIDs of InputData_Intra, whilst features must be unique")
            logger::log_trace(paste0("Error ", message))
            stop(message)
        }
        if (length(InputData_CoRe[duplicated(InputData_CoRe[[FeatureID]]), FeatureID]) > 0) {
            message <- paste0("Duplicated FeatureIDs of InputData_CoRe, whilst features must be unique")
            logger::log_trace(paste0("Error ", message))
            stop(message)
        }
    }

    #------------- SettingsInfo
    if (!is.null(SettingsInfo_C1)) {
        ## C1
        ## ValueCol
        if ("ValueCol" %in% names(SettingsInfo_C1)) {
            if (!SettingsInfo_C1[["ValueCol"]] %in% colnames(InputData_C1)) {
                message <- paste0("The ", SettingsInfo_C1[["ValueCol"]], " column selected as ValueCol in SettingsInfo_C1 was not found in InputData_C1. Please check your input.")
                logger::log_trace(paste0("Error ", message))
                stop(message)
            }
        }
        ## StatCol
        if ("StatCol" %in% names(SettingsInfo_C1)) {
            if (!SettingsInfo_C1[["StatCol"]] %in% colnames(InputData_C1)) {
                message <- paste0("The ", SettingsInfo_C1[["StatCol"]], " column selected as StatCol in SettingsInfo_C1 was not found in InputData_C1. Please check your input.")
                logger::log_trace(paste0("Error ", message))
                stop(message)
            }
        }

        ## C2
        ## ValueCol
        if ("ValueCol" %in% names(SettingsInfo_C2)) {
            if (!SettingsInfo_C2[["ValueCol"]] %in% colnames(InputData_C2)) {
                message <- paste0("The ", SettingsInfo_C2[["ValueCol"]], " column selected as ValueCol in SettingsInfo_C2 was not found in InputData_C2. Please check your input.")
                logger::log_trace(paste0("Error ", message))
                stop(message)
            }
        }
        ## StatCol
        if ("StatCol" %in% names(SettingsInfo_C2)) {
            if (!SettingsInfo_C2[["StatCol"]] %in% colnames(InputData_C2)) {
                message <- paste0("The ", SettingsInfo_C2[["StatCol"]], " column selected as StatCol in SettingsInfo_C2 was not found in InputData_C2. Please check your input.")
                logger::log_trace(paste0("Error ", message))
                stop(message)
            }
        }
    } else {
        ## Intra
        ## ValueCol
        if ("ValueCol" %in% names(SettingsInfo_Intra)) {
            if (!SettingsInfo_Intra[["ValueCol"]] %in% colnames(InputData_Intra)) {
                message <- paste0("The ", SettingsInfo_Intra[["ValueCol"]], " column selected as ValueCol in SettingsInfo_Intra was not found in InputData_Intra. Please check your input.")
                logger::log_trace(paste0("Error ", message))
                stop(message)
            }
        }
        ## StatCol
        if ("StatCol" %in% names(SettingsInfo_Intra)) {
            if (!SettingsInfo_Intra[["StatCol"]] %in% colnames(InputData_Intra)) {
                message <- paste0("The ", SettingsInfo_Intra[["StatCol"]], " column selected as StatCol in SettingsInfo_Intra was not found in InputData_Intra. Please check your input.")
                logger::log_trace(paste0("Error ", message))
                stop(message)
            }
        }

        ## CoRe
        ## ValueCol
        if ("ValueCol" %in% names(SettingsInfo_CoRe)) {
            if (!SettingsInfo_CoRe[["ValueCol"]] %in% colnames(InputData_CoRe)) {
                message <- paste0("The ", SettingsInfo_CoRe[["ValueCol"]], " column selected as ValueCol in SettingsInfo_CoRe was not found in InputData_CoRe. Please check your input.")
                logger::log_trace(paste0("Error ", message))
                stop(message)
            }
        }
        ## StatCol
        if ("StatCol" %in% names(SettingsInfo_CoRe)) {
            if (!SettingsInfo_CoRe[["StatCol"]] %in% colnames(InputData_CoRe)) {
                message <- paste0("The ", SettingsInfo_CoRe[["StatCol"]], " column selected as StatCol in SettingsInfo_CoRe was not found in InputData_CoRe. Please check your input.")
                logger::log_trace(paste0("Error ", message))
                stop(message)
            }
        }

        ## StatCol
        if ("DirectionCol" %in% names(SettingsInfo_CoRe)) {
            if (!SettingsInfo_CoRe[["DirectionCol"]] %in% colnames(InputData_CoRe)) {
                message <- paste0("The ", SettingsInfo_CoRe[["DirectionCol"]], " column selected as DirectionCol in SettingsInfo_CoRe was not found in InputData_CoRe. Please check your input.")
                logger::log_trace(paste0("Error ", message))
                stop(message)
            }
        }
    }

    #------------- SettingsInfo Cutoffs: ## EDIT: this can be simplified, when you are checking for C1 and C2 for the same attributes, write a function and apply on C1 and C2 separately
    if (!is.null(SettingsInfo_C1)) {
        if (is.na(as.numeric(SettingsInfo_C1[["StatCutoff"]])) |as.numeric(SettingsInfo_C1[["StatCutoff"]]) > 1 | as.numeric(SettingsInfo_C1[["StatCutoff"]]) < 0) {
            message <- paste0("Check input. The selected StatCutoff in SettingsInfo_C1 should be numeric and between 0 and 1.")
            logger::log_trace(paste0("Error ", message))
            stop(message)
        }

        if (is.na(as.numeric(SettingsInfo_C2[["StatCutoff"]])) |as.numeric(SettingsInfo_C2[["StatCutoff"]]) > 1 | as.numeric(SettingsInfo_C2[["StatCutoff"]]) < 0) {
            message <- paste0("Check input. The selected StatCutoff in SettingsInfo_C2 should be numeric and between 0 and 1.")
            logger::log_trace(paste0("Error ", message))
            stop(message)
        }

        if (is.na(as.numeric(SettingsInfo_C1[["ValueCutoff"]]))) {
            message <- paste0("Check input. The selected ValueCutoff in SettingsInfo_C1 should be numeric and between 0 and 1.")
            logger::log_trace(paste0("Error ", message))
            stop(message)
        }

        if (is.na(as.numeric(SettingsInfo_C2[["ValueCutoff"]]))) {
            message <- paste0("Check input. The selected ValueCutoff in SettingsInfo_C2 should be numeric and between 0 and 1.")
            logger::log_trace(paste0("Error ", message))
            stop(message)
        }

    } else {
        if (is.na(as.numeric(SettingsInfo_Intra[["StatCutoff"]])) | as.numeric(SettingsInfo_Intra[["StatCutoff"]]) > 1 | as.numeric(SettingsInfo_Intra[["StatCutoff"]]) < 0) {
            message <- paste0("Check input. The selected StatCutoff in SettingsInfo_Intra should be numeric and between 0 and 1.")
            logger::log_trace(paste0("Error ", message))
            stop(message)
        }

        if (is.na(as.numeric(SettingsInfo_CoRe[["StatCutoff"]])) |as.numeric(SettingsInfo_CoRe[["StatCutoff"]]) > 1 | as.numeric(SettingsInfo_CoRe[["StatCutoff"]]) < 0) {
            message <- paste0("Check input. The selected StatCutoff in SettingsInfo_CoRe should be numeric and between 0 and 1.")
            logger::log_trace(paste0("Error ", message))
            stop(message)
        }

        if (is.na(as.numeric(SettingsInfo_Intra[["ValueCutoff"]]))) {
            message <- paste0("Check input. The selected ValueCutoff in SettingsInfo_Intra should be numeric and between 0 and 1.")
            logger::log_trace(paste0("Error ", message))
            stop(message)
        }

        if (is.na(as.numeric(SettingsInfo_CoRe[["ValueCutoff"]]))) {
            message <- paste0("Check input. The selected ValueCutoff in SettingsInfo_CoRe should be numeric and between 0 and 1.")
            logger::log_trace(paste0("Error ", message))
            stop(message)
        }
    }

    #------------ NAs in data
    if (!is.null(InputData_C1)) {
        if (nrow(InputData_C1[complete.cases(InputData_C1[[SettingsInfo_C1[["ValueCol"]]]], InputData_C1[[SettingsInfo_C1[["StatCol"]]]]), ]) < nrow(InputData_C1)) {
            message <- paste0("InputData_C1 includes NAs in ", SettingsInfo_C1[["ValueCol"]], " and/or in ", SettingsInfo_C1[["StatCol"]], ". ", nrow(InputData_C1) - nrow(InputData_C1[complete.cases(InputData_C1[[SettingsInfo_C1[["ValueCol"]]]], InputData_C1[[SettingsInfo_C1[["StatCol"]]]]), ]), " metabolites containing NAs are removed.")
            logger::log_trace(paste0("Warning ", messag))
            warning(message)
        }

        if (nrow(InputData_C2[complete.cases(InputData_C2[[SettingsInfo_C2[["ValueCol"]]]], InputData_C2[[SettingsInfo_C2[["StatCol"]]]]), ]) < nrow(InputData_C2)) {
            message <- paste0("InputData_C2 includes NAs in ", SettingsInfo_C2[["ValueCol"]], " and/or in", SettingsInfo_C2[["StatCol"]], ". ", nrow(InputData_C2) - nrow(InputData_C2[complete.cases(InputData_C2[[SettingsInfo_C2[["ValueCol"]]]], InputData_C2[[SettingsInfo_C2[["StatCol"]]]]), ]), " metabolites containing NAs are removed.")
            logger::log_trace(paste0("Warning ", message))
            warning(message)
        }
    } else {
        if (nrow(InputData_Intra[complete.cases(InputData_Intra[[SettingsInfo_Intra[["ValueCol"]]]], InputData_Intra[[SettingsInfo_Intra[["StatCol"]]]]), ]) < nrow(InputData_Intra)) {
            message <- paste0("InputData_Intra includes NAs in ", SettingsInfo_Intra[["ValueCol"]], " and/or in ", SettingsInfo_Intra[["StatCol"]], ". ", nrow(InputData_Intra) - nrow(InputData_Intra[complete.cases(InputData_Intra[[SettingsInfo_Intra[["ValueCol"]]]], InputData_Intra[[SettingsInfo_Intra[["StatCol"]]]]), ]), " metabolites containing NAs are removed.")
            logger::log_trace(paste0("Warning ", message))
            warning(message)
        }

        if (nrow(InputData_CoRe[complete.cases(InputData_CoRe[[SettingsInfo_CoRe[["ValueCol"]]]], InputData_CoRe[[SettingsInfo_CoRe[["StatCol"]]]]), ]) < nrow(InputData_CoRe)) {
            message <- paste0("InputData_CoRe includes NAs in ", SettingsInfo_CoRe[["ValueCol"]], " and/or in ", SettingsInfo_CoRe[["StatCol"]], ". ", nrow(InputData_CoRe) - nrow(InputData_CoRe[complete.cases(InputData_CoRe[[SettingsInfo_CoRe[["ValueCol"]]]], InputData_CoRe[[SettingsInfo_CoRe[["StatCol"]]]]), ]), " metabolites containing NAs are removed.")
            logger::log_trace(paste0("Warning ", message))
            warning(message)
        }
    }

    #------------- BackgroundMethod
    if (!is.null(SettingsInfo_C1)) {
        options <- c("C1|C2", "C1&C2", "C2", "C1" , "*")
        if (!any(options %in% BackgroundMethod)) {
            message <- paste0("Check input. The selected BackgroundMethod option is not valid. Please select one of the folowwing: ", paste(options, collapse = ", "), "." )
            logger::log_trace(paste0("Error ", message))
            stop(message)
        }
    } else {
        options <- c("Intra|CoRe", "Intra&CoRe", "CoRe", "Intra" , "*")
        if (!any(options %in% BackgroundMethod)) {
            message <- paste0("Check input. The selected BackgroundMethod option is not valid. Please select one of the folowwing: ", paste(options, collapse = ", "), "." )
            logger::log_trace(paste0("Error ", message))
            stop(message)
        }
    }

    #------------- SaveAs
    SaveAs_Table_options <- c("txt", "csv", "xlsx")
    if (!is.null(SaveAs_Table)) {
        if (!(SaveAs_Table %in% SaveAs_Table_options)| is.null(SaveAs_Table)) {
            message <- paste0("Check input. The selected SaveAs_Table option is not valid. Please select one of the folowwing: ", paste(SaveAs_Table_options, collapse = ", "), "." )
            logger::log_trace(paste0("Error ", message))
            stop(message)
        }
    }
}


