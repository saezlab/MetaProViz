# !/usr/bin/env Rscript

#
#  This file is part of the `MetaProViz` R package
#
#  Copyright 2022-2025
#  Saez Lab, Heidelberg University
#
#  Authors: see the file `README.md`
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file `LICENSE` or copy at
#      https://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: https://saezlab.github.io/MetaProViz
#  Git repo: https://github.com/saezlab/MetaProViz
#


#
# Helper function: Internal Function to check function input
#

#' Check input general parameters
#'
#' @param data Passed to MetaProViz functions. Usually DF which contains unique sample
#'     identifiers as row names and metabolite numerical values in columns with
#'     metabolite identifiers as column names. But can also be differential
#'     expression results or other data.
#' @param data_num \emph{Optional: } If data must be numeric \strong{Default = TRUE}
#' @param metadata_sample \emph{Optional: } DF which contains information about the samples, which
#'     will be combined with the input data based on the unique sample
#'     identifiers used as rownames. If not avaliable can be set to NULL.
#'     \strong{Default = NULL}
#' @param metadata_feature \emph{Optional: } DF which contains information about the features. If
#'     not avaliable can be set to NULL. \strong{Default = NULL}
#' @param metadata_info \emph{Optional: } Passed to MetaProViz functions. Usually named vector
#'     containing the information about the names of the experimental
#'     parameters in metadata_sample, metadata_feature or data. \strong{Default
#'     = NULL}
#' @param save_plot \emph{Optional: } Select the file type of output plots. Options are svg,
#'     png, pdf. If set to NULL, plots are not saved.\strong{Default = NULL}
#' @param save_table \emph{Optional: } Select the file type of output table. Options are
#'     "csv", "xlsx", "txt". If set to NULL, plots are not saved.
#'     \strong{Default = NULL}
#' @param core \emph{Optional: } If TRUE, a consumption-release experiment has been
#'     performed. If not avaliable can be set to NULL. \strong{Default = FALSE}
#' @param print_plot \emph{Optional: } If TRUE prints an overview of resulting plots. If not
#'     avaliable can be set to NULL. \strong{Default = FALSE}
#' @param theme \emph{Optional: } Can be set for VizX functions. If not avaliable can be
#'     set to NULL. \strong{Default = NULL}
#' @param plot_types \emph{Optional: } Needs to be set for VizX functions. Options are
#'     "Sample", "Feature", Both". This refers to metadata_info color, shape,
#'     individual as for some plots we have both feature and sample settings.
#'     \strong{Default = NULL}
#'
#' @return returns warnings and errors if input is not correct
#'
#' @importFrom logger log_trace
#' @noRd
check_param <- function(
    data,
    data_num=TRUE,
    metadata_sample=NULL,
    metadata_feature=NULL,
    metadata_info=NULL,
    save_plot=NULL,
    save_table=NULL,
    core=FALSE,
    print_plot=FALSE,
    theme=NULL,
    plot_types=NULL
) {
    # # ## # ## # ## # ### Parameters valid for multiple MetaProViz functions

    # -------------data
    if (is.data.frame(data) == FALSE) {
        message <-
            paste0(
                "data should be a data.frame. It's currently a ",
                paste(class(data)),
                ".",
                sep = ""
            )
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }
    if (any(duplicated(row.names(data))) == TRUE) {
        message <- paste0("Duplicated row.names of data, whilst row.names must be unique")
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }

    if (data_num == TRUE) {
        Test_num <- apply(data, 2, function(x) is.numeric(x))
    if ((any(Test_num) ==  FALSE) ==  TRUE) {
        message <- paste0("data needs to be of class numeric")
        log_trace(paste("Error ", message, sep =""))
        stop(message)
        }
    }

    if (sum(duplicated(colnames(data))) > 0) {
        message <- paste0("data contained duplicates column names, whilst col.names must be unique.")
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }

    # -------------SettingsFile
    if (is.null(metadata_sample) == FALSE) {
        Test_match <-
            merge(
                metadata_sample,
                data,
                by = "row.names",
                all =  FALSE
            )
    if (nrow(Test_match) ==  0) {
        message <- paste0("row.names data need to match row.names metadata_sample.")
        log_trace(paste("Error ", message, sep =""))
        stop(message)
        }
    }

    if (is.null(metadata_feature) == FALSE) {
        Test_match <-
            merge(
                metadata_feature,
                as.data.frame(t(data)),
                by = "row.names",
                all =  FALSE
            )
    if (nrow(Test_match) ==  0) {
        stop("col.names data need to match row.names metadata_feature.")
    }
    }

    # -------------metadata_info
    if (is.vector(metadata_info) == FALSE & is.null(metadata_info) == FALSE) {
        message <-
            paste0(
                "metadata_info should be NULL or a vector. It's currently a ",
                paste(class(metadata_info), ".", sep = "")
            )
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }

    if (is.null(metadata_info) == FALSE) {
    # Conditions
    if ("Conditions" %in% names(metadata_info)) {
        if (metadata_info[["Conditions"]] %in% colnames(metadata_sample) == FALSE) {
            message <-
                paste0(
                    "The ",
                    metadata_info[["Conditions"]],
                    " column selected as Conditions in metadata_info was not found in SettingsFile. Please check your input."
                )
            log_trace(paste("Error ", message, sep =""))
            stop(message)
        }
    }

    # Biological replicates
    if ("Biological_Replicates" %in% names(metadata_info)) {
        if (metadata_info[["Biological_Replicates"]] %in% colnames(metadata_sample) == FALSE) {
            message <-
                paste0(
                    "The ",
                    metadata_info[["Biological_Replicates"]],
                    " column selected as Biological_Replicates in metadata_info was not found in metadata_sample. Please check your input."
                )
            log_trace(paste("Error ", message, sep =""))
            stop(message)
        }
    }

    # Numerator
    if ("Numerator" %in% names(metadata_info) == TRUE) {
        if (metadata_info[["Numerator"]] %in% metadata_sample[[metadata_info[["Conditions"]]]] == FALSE) {
            message <-
                paste0(
                    "The ",
                    metadata_info[["Numerator"]],
                    " column selected as numerator in metadata_info was not found in metadata_sample. Please check your input."
                )
            log_trace(paste("Error ", message, sep =""))
            stop(message)
        }
    }

    # Denominator
    if ("Denominator" %in% names(metadata_info) == TRUE) {
        if (metadata_info[["Denominator"]] %in% metadata_sample[[metadata_info[["Conditions"]]]] == FALSE) {
            message <-
                paste0(
                    "The ",
                    metadata_info[["Denominator"]],
                    " column selected as denominator in metadata_info was not found in metadata_sample. Please check your input."
                )
            log_trace(paste("Error ", message, sep =""))
            stop(message)
        }
    }

    # Denominator & Numerator
    if ("Denominator" %in% names(metadata_info) == FALSE  & "Numerator" %in% names(metadata_info) == TRUE) {
        message <-
            paste0(
                "Check input. The selected denominator option is empty while ",
                paste(metadata_info[["Numerator"]]),
                " has been selected as a numerator. Please add a denominator for 1-vs-1 comparison or remove the numerator for all-vs-all comparison."
            )
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }

    # Superplot
    if ("Superplot" %in% names(metadata_info)) {
        if (metadata_info[["Superplot"]] %in% colnames(metadata_sample) == FALSE) {
            message <-
                paste0(
                    "The ",
                    metadata_info[["Superplot"]],
                    " column selected as Superplot column in metadata_info was not found in metadata_sample. Please check your input."
                )
            log_trace(paste("Error ", message, sep =""))
            stop(message)
        }
    }

    if (is.null(plot_types) == FALSE) {
        if (plot_types == "Sample") {
        # Plot colour
        if ("color" %in% names(metadata_info)) {
            if (metadata_info[["color"]] %in% colnames(metadata_sample) == FALSE) {
                message <-
                    paste0(
                        "The ",
                        metadata_info[["color"]],
                        " column selected as color in metadata_info was not found in metadata_sample. Please check your input."
                    )
                log_trace(paste("Error ", message, sep =""))
                stop(message)
            }
        }

        # Plot shape
        if ("shape" %in% names(metadata_info)) {
            if (metadata_info[["shape"]] %in% colnames(metadata_sample) == FALSE) {
                message <-
                    paste0(
                        "The ",
                        metadata_info[["shape"]],
                        " column selected as shape in metadata_info was not found in metadata_sample. Please check your input."
                    )
                log_trace(paste("Error ", message, sep =""))
                stop(message)
            }
        }

        # Plot individual
        if ("individual" %in% names(metadata_info)) {
            if (metadata_info[["individual"]] %in% colnames(metadata_sample) == FALSE) {
                message <-
                    paste0(
                        "The ",
                        metadata_info[["individual"]],
                        " column selected as individual in metadata_info was not found in metadata_sample. Please check your input."
                    )
                log_trace(paste("Error ", message, sep =""))
                stop(message)
            }
        }
        }else if (plot_types == "Feature") {
        if ("color" %in% names(metadata_info)) {
            if (metadata_info[["color"]] %in% colnames(metadata_feature) == FALSE) {
                message <-
                    paste0(
                        "The ",
                        metadata_info[["color"]],
                        " column selected as color in metadata_info was not found in metadata_feature. Please check your input."
                    )
                log_trace(paste("Error ", message, sep =""))
                stop(message)
            }
        }

        # Plot shape
        if ("shape" %in% names(metadata_info)) {
            if (metadata_info[["shape"]] %in% colnames(metadata_feature) == FALSE) {
                message <-
                    paste0(
                        "The ",
                        metadata_info[["shape"]],
                        " column selected as shape in metadata_info was not found in metadata_feature. Please check your input."
                    )
                log_trace(paste("Error ", message, sep =""))
                stop(message)
            }
        }

        # Plot individual
        if ("individual" %in% names(metadata_info)) {
            if (metadata_info[["individual"]] %in% colnames(metadata_feature) == FALSE) {
                message <-
                    paste0(
                        "The ",
                        metadata_info[["individual"]],
                        " column selected as individual in metadata_info was not found in metadata_feature. Please check your input."
                    )
                log_trace(paste("Error ", message, sep =""))
                stop(message)
            }
        }
        }else if (plot_types == "Both") {
        # Plot colour sample
        if ("color_Sample" %in% names(metadata_info)) {
            if (metadata_info[["color_Sample"]] %in% colnames(metadata_sample) == FALSE) {
                message <-
                    paste0(
                        "The ",
                        metadata_info[["color_Sample"]],
                        " column selected as color_Sample in metadata_info was not found in metadata_sample. Please check your input."
                    )
                log_trace(paste("Error ", message, sep =""))
                stop(message)
            }
        }

        # Plot colour Metab
        if ("color_Metab" %in% names(metadata_info)) {
            if (metadata_info[["color_Metab"]] %in% colnames(metadata_feature) == FALSE) {
                message <-
                    paste0(
                        "The ",
                        metadata_info[["color_Metab"]],
                        " column selected as color_Metab in metadata_info was not found in metadata_feature. Please check your input."
                    )
                log_trace(paste("Error ", message, sep =""))
                stop(message)
            }
            if (sum(colnames(data) %in% metadata_feature$Metabolite) < length(data)  ) {
                message <- paste0("The data contains metabolites not found in metadata_feature.")
                log_trace(paste("Warning ", message, sep =""))
                warning(message)
            }
        }

        # Plot shape_metab
        if ("shape_Metab" %in% names(metadata_info)) {
            if (metadata_info[["shape_Metab"]] %in% colnames(metadata_feature) == FALSE) {
                message <-
                    paste0(
                        "The ",
                        metadata_info[["shape_Metab"]],
                        " column selected as shape_Metab in metadata_info was not found in metadata_feature. Please check your input."
                    )
                log_trace(paste("Error ", message, sep =""))
                stop(message)
            }
        }

        # Plot shape_metab
        if ("shape_Sample" %in% names(metadata_info)) {
            if (metadata_info[["shape_Sample"]] %in% colnames(metadata_feature) == FALSE) {
                message <-
                    paste0(
                        "The ",
                        metadata_info[["shape_Sample"]],
                        " column selected as shape_Metab in metadata_info was not found in metadata_sample. Please check your input."
                    )
                log_trace(paste("Error ", message, sep =""))
                stop(message)
            }
        }

        # Plot individual_Metab
        if ("individual_Metab" %in% names(metadata_info)) {
            if (metadata_info[["individual_Metab"]] %in% colnames(metadata_feature) == FALSE) {
                message <-
                    paste0(
                        "The ",
                        metadata_info[["individual_Metab"]],
                        " column selected as individual_Metab in metadata_info was not found in metadata_feature. Please check your input."
                    )
                log_trace(paste("Error ", message, sep =""))
                stop(message)
            }
        }

        # Plot individual_Sample
        if ("individual_Sample" %in% names(metadata_info)) {
            if (metadata_info[["individual_Sample"]] %in% colnames(metadata_sample) == FALSE) {
                message <-
                    paste0(
                        "The ",
                        metadata_info[["individual_Sample"]],
                        " column selected as individual_Sample in metadata_info was not found in metadata_sample. Please check your input."
                    )
                log_trace(paste("Error ", message, sep =""))
                stop(message)
            }
        }

        }

    }
    }

    # -------------SaveAs
    Save_as_Plot_options <- c("svg","pdf", "png")
    if (is.null(save_plot) == FALSE) {
        if (save_plot %in% Save_as_Plot_options == FALSE) {
        message <-
            paste0(
                "Check input. The selected save_plot option is not valid. Please select one of the folowwing: ",
                paste(Save_as_Plot_options,collapse = ", "),
                " or set to NULL if no plots should be saved."
            )
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }
    }

    save_table_options <-
        c(
            "txt",
            "csv",
            "xlsx",
            "Rdata") # Rdata = SummarizedExperiment (?
        )
    if (is.null(save_table) == FALSE) {
        if ((save_table %in% save_table_options == FALSE) | (is.null(save_table) == TRUE)) {
        message <-
            paste0(
                "Check input. The selected save_table option is not valid. Please select one of the folowwing: ",
                paste(save_table_options,collapse = ", "),
                " or set to NULL if no tables should be saved."
            )
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }
    }

    # -------------core
    if (is.logical(core) == FALSE) {
        message <- paste0("Check input. The core value should be either = TRUE for preprocessing of Consuption/Release experiment or = FALSE if not.")
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }

    # -------------theme
    if (is.null(theme) == FALSE) {
        theme_options <-
            c(
                "theme_grey()",
                "theme_gray()",
                "theme_bw()",
                "theme_linedraw()",
                "theme_light()",
                "theme_dark()",
                "theme_minimal()",
                "theme_classic()",
                "theme_void()",
                "theme_test()"
            )
    if (theme %in% theme_options == FALSE) {
        message <-
            paste0(
                "Check input. theme option is incorrect. You can check for complete themes here: https://ggplot2.tidyverse.org/reference/ggtheme.html. Options are the following: ",
                paste(theme_options, collapse = ", "),
                "."
            )
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }
    }
    # ------------- general
    if (is.logical(print_plot) == FALSE) {
        message <- paste0("Check input. print_plot should be either = TRUE or = FALSE.")
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }
}


#
# processing helper function: Internal Function to check function input
#

#' Check specific input parameters for processing()
#'
#' @param metadata_sample DF which contains information about the samples, which will be combined
#'     with the input data based on the unique sample identifiers used as
#'     rownames.
#' @param metadata_info Named vector containing the information about the names of the
#'     experimental parameters. c(Conditions="ColumnName_Plot_SettingsFile",
#'     Biological_Replicates="ColumnName_Plot_SettingsFile"). For core = TRUE a
#'     core_norm_factor = "Columnname_Input_SettingsFile" and core_media =
#'     "Columnname_Input_SettingsFile", have to also be added.
#' @param \emph{Optional: }core If TRUE, a consumption-release experiment has been performed and
#'     the core value will be calculated.\strong{Default = FALSE}
#' @param featurefilt \emph{Optional: }If NULL, no feature filtering is performed. If set to
#'     "Standard" then it applies the 80%-filtering rule (Bijlsma S. et al.,
#'     2006) on the metabolite features on the whole dataset. If is set to
#'     "Modified",filtering is done based on the different conditions, thus a
#'     column named "Conditions" must be provided in the Input_SettingsFile
#'     input file including the individual conditions you want to apply the
#'     filtering to (Yang, J et al., 2015). \strong{Default = "Standard"}
#' @param cutoff_featurefilt \emph{Optional: } percentage of feature filtering. \strong{Default =
#'     0.8}
#' @param tic \emph{Optional: } If TRUE, total Ion Count normalization is performed.
#'     \strong{Default = TRUE}
#' @param mvi \emph{Optional: } If TRUE, Missing Value Imputation (mvi) based on half
#'     minimum is performed \strong{Default = TRUE}
#' @param mvi_percentage \emph{Optional: } percentage 0-100 of imputed value based on the minimum
#'     value. \strong{Default = 50}
#' @param hotellins_confidence \emph{Optional: } Defines the Confidence of Outlier identification in
#'     HotellingT2 test. Must be numeric.\strong{Default = 0.99}
#'
#' @return returns warnings and errors if input is not correct
#'
#' @importFrom logger log_trace
#' @noRd
check_param_processing <- function(
    metadata_sample,
    metadata_info,
    core=FALSE,
    featurefilt = "Modified",
    cutoff_featurefilt = 0.8,
    tic = TRUE,
    mvi= TRUE,
    mvi_percentage=50,
    hotellins_confidence = 0.99
) {
    if (is.vector(metadata_info) == TRUE) {
    # -------------metadata_info
    # core
    if (core == TRUE) {   # parse core normalisation factor
        message <- paste0("For Consumption Release experiment we are using the method from Jain M.  REF: Jain et. al, (2012), Science 336(6084):1040-4, doi: 10.1126/science.1218595.")
        log_trace(paste("Message ", message, sep =""))
        message(message)
        if ("core_media" %in% names(metadata_info)) {
            if (length(grep(metadata_info[["core_media"]], metadata_sample[[metadata_info[["Conditions"]]]])) < 1) {     # Check for core_media samples
            message <- paste0("No core_media samples were provided in the 'Conditions' in the metadata_sample. For a core experiment control media samples without cells have to be measured and be added in the 'Conditions'
                            column labeled as 'core_media' (see @param section). Please make sure that you used the correct labelling or whether you need core = FALSE for your analysis")
            log_trace(paste("Error ", message, sep =""))
            stop(message)
        }
        }

        if ("core_norm_factor" %in% names(metadata_info) == FALSE) {
            message <- paste0("No growth rate or growth factor provided for normalising the core result, hence core_norm_factor set to 1 for each sample")
            log_trace(paste("Warning ", message, sep =""))
            warning(message)
        }
    }
    }

    # -------------General parameters
    Feature_Filtering_options <- c("Standard","Modified")
    if (featurefilt %in% Feature_Filtering_options == FALSE & is.null(featurefilt) == FALSE) {
        message <-
            paste0(
                "Check input. The selected featurefilt option is not valid. Please set to NULL or select one of the folowwing: ",
                paste(Feature_Filtering_options,collapse = ", "),
                "."
            )
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }
    if (is.numeric(cutoff_featurefilt) == FALSE | cutoff_featurefilt > 1 | cutoff_featurefilt < 0) {
        message <- paste0("Check input. The selected cutoff_featurefilt should be numeric and between 0 and 1.")
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }
    if (is.logical(tic) == FALSE) {
        message <- paste0("Check input. The tic value should be either `TRUE` if tic normalization is to be performed or `FALSE` if no data normalization is to be applied.")
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }
    if (is.logical(mvi) == FALSE) {
        message <- paste0("Check input. The mvi value should be either `TRUE` if mising value imputation should be performed or `FALSE` if not.")
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }
    if (is.numeric(mvi_percentage) == FALSE | hotellins_confidence > 100 | hotellins_confidence < 0) {
        message <- paste0("Check input. The selected mvi_percentage value should be numeric and between 0 and 100.")
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }
    if ( is.numeric(hotellins_confidence) == FALSE | hotellins_confidence > 1 | hotellins_confidence < 0) {
        message <- paste0("Check input. The selected hotellins_confidence value should be numeric and between 0 and 1.")
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }
}


#
# dma helper function: Internal Function to check function input
#

#' Check input parameters of dma
#'
#' @param data DF with unique sample identifiers as row names and metabolite numerical
#'     values in columns with metabolite identifiers as column names.
#' @param metadata_sample DF which contains metadata information about the samples, which will be
#'     combined with your input data based on the unique sample identifiers
#'     used as rownames.
#' @param metadata_info Named vector including the information about the conditions column
#'     information on numerator or denominator
#'     c(Conditions="ColumnName_SettingsFile", Numerator =
#'     "ColumnName_SettingsFile", Denominator  = "ColumnName_SettingsFile").
#'     \strong{Default = c(conditions="Conditions", numerator = NULL,
#'     denumerator = NULL)}
#' @param pval \emph{Optional: } String which contains an abbreviation of the selected
#'     test to calculate p.value. For one-vs-one comparisons choose t.test,
#'     wilcox.test, "chisq.test", "cor.test" or lmFit (=limma), for one-vs-all
#'     or all-vs-all comparison choose aov (=anova), welch(=welch anova),
#'     kruskal.test or lmFit (=limma) \strong{Default = "lmFit"}
#' @param padj \emph{Optional: } String which contains an abbreviation of the selected
#'     p.adjusted test for p.value correction for multiple Hypothesis testing.
#'     Search: ?p.adjust for more methods:"BH", "fdr", "bonferroni", "holm",
#'     etc.\strong{Default = "fdr"}
#' @param vst TRUE or FALSE for whether to use variance stabilizing transformation on
#'     the data when linear modeling is used for hypothesis testing.
#'     \strong{Default = FALSE}
#' @param shapiro TRUE or FALSE for whether to perform the shapiro.test and get informed
#'     about data distribution (normal versus not-normal distribution.
#'     \strong{Default = TRUE}
#' @param bartlett TRUE or FALSE for whether to perform the bartlett.test. \strong{Default
#'     = TRUE}
#' @param transform TRUE or FALSE. If TRUE we expect the data to be not log2 transformed and
#'     log2 transformation will be performed within the limma function and
#'     Log2FC calculation. If FALSE we expect the data to be log2 transformed
#'     as this impacts the Log2FC calculation and limma. \strong{Default= TRUE}
#'
#' @return Returns: 1. warnings and errors if input is not correct, 2. Settings
#'
#' @importFrom logger log_trace log_info
#' @importFrom dplyr filter select_if
#' @importFrom magrittr %>%
#' @importFrom utils combn
#' @noRd
check_param_dma <- function(
    data,
    metadata_sample,
    metadata_info= c(Conditions="Conditions", Numerator = NULL, Denominator  = NULL),
    pval ="lmFit",
    padj="fdr",
    vst=FALSE,
    shapiro =TRUE,
    bartlett =TRUE,
    transform=TRUE
) {

    # -------------metadata_info
    if (is.null(metadata_info) == TRUE) {
        message <- paste0("You have to provide metadata_info's for Conditions.")
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }

    ## ------------ Denominator/numerator ----------- ##
    # Denominator and numerator: Define if we compare one_vs_one, one_vs_all or all_vs_all.
    if ("Denominator" %in% names(metadata_info) == FALSE  & "Numerator" %in% names(metadata_info) == FALSE) {
    # all-vs-all: Generate all pairwise combinations
        conditions <- metadata_sample[[metadata_info[["Conditions"]]]]
        denominator <- unique(metadata_sample[[metadata_info[["Conditions"]]]])
        numerator <- unique(metadata_sample[[metadata_info[["Conditions"]]]])
        comparisons <- combn(unique(conditions), 2) %>% as.matrix()
    # Settings:
        MultipleComparison <- TRUE
        all_vs_all <- TRUE
    }else if ("Denominator" %in% names(metadata_info) == TRUE  & "Numerator" %in% names(metadata_info) == FALSE) {
    # all-vs-one: Generate the pairwise combinations
    conditions <- metadata_sample[[metadata_info[["Conditions"]]]]
    denominator <- metadata_info[["Denominator"]]
    numerator <- unique(metadata_sample[[metadata_info[["Conditions"]]]])
    # Remove denom from num
    numerator <- numerator[!numerator %in% denominator]
    comparisons  <- t(expand.grid(numerator, denominator)) %>% as.data.frame()
    # Settings:
    MultipleComparison <- TRUE
    all_vs_all <- FALSE
    }else if ("Denominator" %in% names(metadata_info) == TRUE  & "Numerator" %in% names(metadata_info) == TRUE) {
    # one-vs-one: Generate the comparisons
    denominator <- metadata_info[["Denominator"]]
    numerator <- metadata_info[["Numerator"]]
    comparisons <- matrix(c(metadata_info[["Denominator"]], metadata_info[["Numerator"]]))
    # Settings:
    MultipleComparison <- FALSE
    all_vs_all <- FALSE
    }

    ## ------------ Test statistics ----------- ##
    if (MultipleComparison == FALSE) {
        STAT_pval_options <-
            c(
                "t.test",
                "wilcox.test",
                "chisq.test",
                "cor.test",
                "lmFit"
            )
    if (pval %in% STAT_pval_options == FALSE) {
        message <-
            paste0(
                "Check input. The selected pval option for Hypothesis testing is not valid for multiple comparison (one-vs-all or all-vs-all). Please select one of the following: ",
                paste(STAT_pval_options,collapse = ", "),
                " or specify numerator and denumerator."
            )
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }
    } else {
    STAT_pval_options <- c("aov", "kruskal.test", "welch" ,"lmFit")
    if (pval %in% STAT_pval_options == FALSE) {
        message <-
            paste0(
                "Check input. The selected pval option for Hypothesis testing is not valid for one-vs-one comparsion. Multiple comparison is selected. Please select one of the following: ",
                paste(STAT_pval_options,collapse = ", "),
                " or change numerator and denumerator."
            )
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }
    }

    STAT_padj_options <-
        c(
            "holm",
            "hochberg",
            "hommel",
            "bonferroni",
            "BH",
            "BY",
            "fdr",
            "none"
        )
    if (padj %in% STAT_padj_options == FALSE) {
        message <-
            paste0(
                "Check input. The selected padj option for multiple Hypothesis testing correction is not valid. Please select one of the folowing: ",
                paste(STAT_padj_options,collapse = ", "),
                "."
            )
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }

    ## ------------ Sample Numbers ----------- ##
    Num <- data %>% #Are sample numbers enough?
    filter(metadata_sample[[metadata_info[["Conditions"]]]] %in% numerator) %>%
    select_if(is.numeric) # only keep numeric columns with metabolite values
    Denom <- data %>%
    filter(metadata_sample[[metadata_info[["Conditions"]]]] %in% denominator) %>%
    select_if(is.numeric)

    if (nrow(Num) == 1) {
        message <-
            paste0(
                "There is only one sample available for ",
                numerator,
                ", so no statistical test can be performed."
            )
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    } else if (nrow(Denom) == 1) {
    message <-
        paste0(
            "There is only one sample available for ",
            denominator,
            ", so no statistical test can be performed."
        )
    log_trace(paste("Error ", message, sep =""))
    stop(message)
    }else if (nrow(Num) == 0) {
    message <- paste0("There is no sample available for ", numerator, ".")
    log_trace(paste("Error ", message, sep =""))
    stop(message)
    }else if (nrow(Denom) == 0) {
    message <- paste0("There is no sample available for ", denominator, ".")
    log_trace(paste("Error ", message, sep =""))
    stop(message)
    }

    ## ------------ Check Missingness ------------- ##
    Num_Miss <- replace(Num, Num == 0, NA)
    Num_Miss <- Num_Miss[, (colSums(is.na(Num_Miss)) > 0), drop = FALSE]

    Denom_Miss <- replace(Denom, Denom == 0, NA)
    Denom_Miss <- Denom_Miss[, (colSums(is.na(Denom_Miss)) > 0), drop = FALSE]

    if ((ncol(Num_Miss)>0 & ncol(Denom_Miss) == 0)) {
        Metabolites_Miss <- colnames(Num_Miss)
    if (ncol(Num_Miss)<=10) {
        message <-
            paste0(
                "In `Numerator` ",
                paste0(toString(numerator)),
                ", NA/0 values exist in ",
                ncol(Num_Miss),
                " Metabolite(s): ",
                paste0(colnames(Num_Miss), collapse = ", "),
                ". Those metabolite(s) might return p.val = NA, p.adj.= NA, t.val = NA. The Log2FC = Inf, if all replicates are 0/NA."
            )
        log_info(message)
        message(message)
    } else {
        message <-
            paste0(
                "In `Numerator` ",
                paste0(toString(numerator)),
                ", NA/0 values exist in ",
                ncol(Num_Miss),
                " Metabolite(s).",
                " Those metabolite(s) might return p.val = NA, p.adj.= NA, t.val = NA. The Log2FC = Inf, if all replicates are 0/NA."
            )
        log_info(message)
        message(message)
    }
    } else if (ncol(Num_Miss) == 0 & ncol(Denom_Miss)>0) {
    Metabolites_Miss <- colnames(Denom_Miss)
    if (ncol(Num_Miss)<=10) {
        message <-
            paste0(
                "In `Denominator` ",
                paste0(toString(denominator)),
                ", NA/0 values exist in ",
                ncol(Denom_Miss),
                " Metabolite(s): ",
                paste0(colnames(Denom_Miss), collapse = ", "),
                ". Those metabolite(s) might return p.val = NA, p.adj.= NA, t.val = NA. The Log2FC = Inf, if all replicates are 0/NA."
            )
        log_info(message)
        message(message)
    } else {
        message <-
            paste0(
                "In `Denominator` ",
                paste0(toString(denominator)),
                ", NA/0 values exist in ",
                ncol(Denom_Miss),
                " Metabolite(s).",
                " Those metabolite(s) might return p.val = NA, p.adj.= NA, t.val = NA. The Log2FC = Inf, if all replicates are 0/NA."
            )
        log_info(message)
        message(message)#
    }
    } else if (ncol(Num_Miss)>0 & ncol(Denom_Miss)>0) {
    Metabolites_Miss <- c(colnames(Num_Miss), colnames(Denom_Miss))
    Metabolites_Miss <- unique(Metabolites_Miss)

    message <-
        paste0(
            "In `Numerator` ",
            paste0(toString(numerator)),
            ", NA/0 values exist in ",
            ncol(Num_Miss),
            " Metabolite(s).",
            " and in `denominator`",
            paste0(toString(denominator)),
            " ",
            ncol(Denom_Miss),
            " Metabolite(s
        ).",
                        ". Those metabolite(s) might return p.val = NA, p.adj.= NA, t.val = NA. The Log2FC = Inf, if all replicates are 0/NA.")
    log_info(message)
    message(message)
    } else {
    message <- paste0("There are no NA/0 values")
    log_info(message)
    message(message)

    Metabolites_Miss <- c(colnames(Num_Miss), colnames(Denom_Miss))
    Metabolites_Miss <- unique(Metabolites_Miss)
    }

    # -------------General parameters
    if (is.logical(vst) == FALSE) {
        message <- paste0("Check input. The vst value should be either = TRUE or = FALSE.")
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }

    if (is.logical(shapiro) == FALSE) {
        message <- paste0("Check input. The shapiro value should be either = TRUE or = FALSE.")
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }
    if (is.logical(bartlett) == FALSE) {
        message <- paste0("Check input. The bartlett value should be either = TRUE or = FALSE.")
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }
    if (is.logical(transform) == FALSE) {
        message <- paste0("Check input. `transform` should be either = TRUE or = FALSE.")
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }

    Settings <-
        list(
            "comparisons"= comparisons,
            "MultipleComparison"= MultipleComparison,
            "all_vs_all"= all_vs_all,
            "Metabolites_Miss"= Metabolites_Miss,
            "denominator"= denominator,
            "numerator"= numerator
        )
    return(invisible(Settings))
}


#
# ORA helper function: Internal Function to check function input
#

#' Check input parameters of ORA
#'
#' @param data DF with metabolite names/metabolite IDs as row names. Metabolite
#'     names/IDs need to match the identifier type (e.g. HMDB IDs) in the
#'     input_pathway.
#' @param metadata_info \emph{Optional: } Pass ColumnName of the column including parameters to
#'     use for cutoff_stat and cutoff_percentage, ColumnName for input_pathway.
#'     For cluster_ora also BackgroundColumn. \strong{c(pvalColumn="p.adj",
#'     percentageColumn="t.val", PathwayTerm= "term", PathwayFeature=
#'     "Metabolite")}
#' @param cutoff_stat \emph{Optional: } p-adjusted value cutoff from ORA results. Must be a
#'     numeric value. \strong{default: 0.05}
#' @param cutoff_percentage \emph{Optional: } percentage cutoff of metabolites that should be
#'     considered for ORA. Selects top/Bottom % of selected percentageColumn,
#'     usually t.val or Log2FC \strong{default: 10}
#' @param input_pathway DF that must include column "term" with the pathway name, column
#'     "Metabolite" with the Metabolite name or ID and column "Description"
#'     with pathway description that will be depicted on the plots.
#' @param pathway_name \emph{Optional: } Name of the input_pathway used \strong{default: ""}
#' @param min_gssize \emph{Optional: } minimum group size in ORA \strong{default: 10}
#' @param max_gssize \emph{Optional: } maximum group size in ORA \strong{default: 1000}
#' @param save_table \emph{Optional: } File types for the analysis results are: "csv",
#'     "xlsx", "txt" \strong{default: "csv"}
#' @param remove_background For cluster_ora the Background Settings are passed, for standard_ora set
#'     to FALSE.
#'
#' @return Returns: 1. warnings and errors if input is not correct, 2. Pathway file
#'
#' @importFrom logger log_trace
#' @importFrom dplyr rename
#' @noRd
check_param_ora <- function(
    data,
    metadata_info=c(pvalColumn="p.adj", percentageColumn="t.val", PathwayTerm= "term", PathwayFeature= "Metabolite"),
    cutoff_stat=0.05,
    cutoff_percentage=10,
    input_pathway,
    pathway_name="",
    min_gssize=10,
    max_gssize=1000,
    save_table="csv",
    remove_background
) {
    # 1. The input data:
    if (!inherits(data, "data.frame")) {
        message <-
            paste0(
                "data should be a data.frame. It's currently a ",
                paste(class(data), ".",sep = "")
            )
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }
    if (any(duplicated(row.names(data))) == TRUE) {
        message <- paste0("Duplicated row.names of data, whilst row.names must be unique")
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }

    # 2. Settings Columns:
    if (is.vector(metadata_info) == FALSE & is.null(metadata_info) == FALSE) {
        message <-
            paste0(
                "metadata_info should be NULL or a vector. It's currently a ",
                paste(class(metadata_info), ".", sep = "")
            )
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }

    if (is.null(metadata_info) == FALSE) {
    # "ClusterColumn"
    if ("ClusterColumn" %in% names(metadata_info)) {
        if (metadata_info[["ClusterColumn"]] %in% colnames(data) == FALSE) {
            message <-
                paste0(
                    "The ",
                    metadata_info[["ClusterColumn"]],
                    " column selected as ClusterColumn in metadata_info was not found in data. Please check your input."
                )
            log_trace(paste("Error ", message, sep =""))
            stop(message)
        }
    }

    # "BackgroundColumn"
    if ("BackgroundColumn" %in% names(metadata_info)) {
        if (metadata_info[["BackgroundColumn"]] %in% colnames(data) == FALSE) {
            message <-
                paste0(
                    "The ",
                    metadata_info[["BackgroundColumn"]],
                    " column selected as BackgroundColumn in metadata_info was not found in data. Please check your input."
                )
            log_trace(paste("Error ", message, sep =""))
            stop(message)
        }
    }

    # "pvalColumn"
    if ("pvalColumn" %in% names(metadata_info)) {
        if (metadata_info[["pvalColumn"]] %in% colnames(data) == FALSE) {
            message <-
                paste0(
                    "The ",
                    metadata_info[["pvalColumn"]],
                    " column selected as pvalColumn in metadata_info was not found in data. Please check your input."
                )
            log_trace(paste("Error ", message, sep =""))
            stop(message)
        }
    }

    # "percentageColumn"
    if ("percentageColumn" %in% names(metadata_info)) {
        if (metadata_info[["percentageColumn"]] %in% colnames(data) == FALSE) {
            message <-
                paste0(
                    "The ",
                    metadata_info[["percentageColumn"]],
                    " column selected as percentageColumn in metadata_info was not found in data. Please check your input."
                )
            log_trace(paste("Error ", message, sep =""))
            stop(message)
        }
    }

    # "PathwayTerm"
    if ("PathwayTerm" %in% names(metadata_info)) {
        if (metadata_info[["PathwayTerm"]] %in% colnames(input_pathway) == FALSE) {
            message <-
                paste0(
                    "The ",
                    metadata_info[["PathwayTerm"]],
                    " column selected as PathwayTerm in metadata_info was not found in input_pathway. Please check your input."
                )
            log_trace(paste("Error ", message, sep =""))
            stop(message)
        } else {
        input_pathway <- input_pathway %>%
            rename("term"= metadata_info[["PathwayTerm"]])
        input_pathway$Description <- input_pathway$term
        }
    } else {
        message <- paste0("metadata_info must provide the column name for PathwayTerm in input_pathway")
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }

    # PathwayFeature
    if ("PathwayFeature" %in% names(metadata_info)) {
        if (metadata_info[["PathwayFeature"]] %in% colnames(input_pathway) == FALSE) {
            message <-
                paste0(
                    "The ",
                    metadata_info[["PathwayFeature"]],
                    " column selected as PathwayFeature in metadata_info was not found in input_pathway. Please check your input."
                )
            log_trace(paste("Error ", message, sep =""))
            stop(message)
        } else {
        input_pathway <- input_pathway %>%
            rename("gene"= metadata_info[["PathwayFeature"]])
        }
    } else {
        message <- paste0("metadata_info must provide the column name for PathwayFeature in input_pathway")
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }

    } else {
    message <- paste0("You must provide metadata_info.")
    log_trace(paste("Error ", message, sep =""))
    stop(message)
    }

    # 3. General Settings
    if (is.character(pathway_name) == FALSE) {
        message <- paste0("Check input. pathway_name must be a character of syntax 'example'.")
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }

    if (is.logical(remove_background) == FALSE) {
        message <- paste0("Check input. remove_background value should be either = TRUE or = FALSE.")
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }

    if (is.numeric(min_gssize) == FALSE) {
        message <- paste0("Check input. The selected min_gssize value should be numeric.")
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }

    if (is.numeric(max_gssize) == FALSE) {
        message <- paste0("Check input. The selected max_gssize value should be numeric.")
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }

    save_table_options <- c("txt","csv", "xlsx")
    if (is.null(save_table) == FALSE) {
        if ((save_table %in% save_table_options == FALSE) | (is.null(save_table) == TRUE)) {
        message <-
            paste0(
                "Check input. The selected save_table option is not valid. Please select one of the folowing: ",
                paste(save_table_options,collapse = ", "),
                "."
            )
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }
    }

    if (is.null(cutoff_stat) == FALSE) {
        if (is.numeric(cutoff_stat) == FALSE | cutoff_stat > 1 |  cutoff_stat < 0) {
        message <- paste0("Check input. The selected cutoff_stat value should be numeric and between 0 and 1.")
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }
    }

    if (is.null(cutoff_percentage) == FALSE) {
        if ( is.numeric(cutoff_percentage) == FALSE  | cutoff_percentage > 100 | cutoff_percentage < 0) {
        message <- paste0("Check input. The selected cutoff_percentage value should be numeric and between 0 and 100.")
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }
    }

    ## -------- Return Pathways --------- # #
    return(invisible(input_pathway))
}


#
# MCA Helper function: Internal Function to check function input
#

#' Check input parameters of MCA
#'
#' @param data_intra For mca_core, otherwise NULL. DF for your data (results from e.g. dma)
#'     containing metabolites in rows with corresponding Log2FC and stat
#'     (p-value, p.adjusted) value columns.
#' @param data_core For mca_core, otherwise NULL. DF for your data (results from e.g. dma)
#'     containing metabolites in rows with corresponding Log2FC and stat
#'     (p-value, p.adjusted) value columns. Here we additionally require
#' @param metadata_info_intra For mca_core, otherwise NULL. Pass ColumnNames and Cutoffs for the
#'     intracellular metabolomics including the value column (e.g. Log2FC,
#'     Log2Diff, t.val, etc) and the stats column (e.g. p.adj, p.val). This
#'     must include:
#'     c(ValueCol=ColumnName_data_intra,StatCol=ColumnName_data_intra,
#'     cutoff_stat= NumericValue, ValueCutoff=NumericValue)
#' @param metadata_info_core For mca_core, otherwise NULL. Pass ColumnNames and Cutoffs for the
#'     consumption-release metabolomics including the direction column, the
#'     value column (e.g. Log2Diff, t.val, etc) and the stats column (e.g.
#'     p.adj, p.val). This must include: c(DirectionCol=
#'     ColumnName_data_core,ValueCol=ColumnName_data_core,StatCol=ColumnName_data_core,
#'     cutoff_stat= NumericValue, ValueCutoff=NumericValue)
#' @param data_c1 For mca_2cond, otherwise NULL. DF for your data (results from e.g. dma)
#'     containing metabolites in rows with corresponding Log2FC and stat
#'     (p-value, p.adjusted) value columns.
#' @param data_c2 For mca_2cond, otherwise NULL. DF for your data (results from e.g. dma)
#'     containing metabolites in rows with corresponding Log2FC and stat
#'     (p-value, p.adjusted) value columns.
#' @param metadata_info_c1 For mca_2cond, otherwise NULL. Pass ColumnNames and Cutoffs for
#'     condition 1 including the value column (e.g. Log2FC, Log2Diff, t.val,
#'     etc) and the stats column (e.g. p.adj, p.val). This must include:
#'     c(ValueCol=ColumnName_data_c1,StatCol=ColumnName_data_c1, cutoff_stat=
#'     NumericValue, ValueCutoff=NumericValue)
#' @param metadata_info_c2 For mca_2cond, otherwise NULL. Pass ColumnNames and Cutoffs for
#'     condition 2 includingthe value column (e.g. Log2FC, Log2Diff, t.val,
#'     etc) and the stats column (e.g. p.adj, p.val). This must include:
#'     c(ValueCol=ColumnName_data_c2,StatCol=ColumnName_data_c2, cutoff_stat=
#'     NumericValue, ValueCutoff=NumericValue)
#' @param feature \emph{Optional: } Column name of Column including the Metabolite
#'     identifiers. This MUST BE THE SAME in each of your Input files.
#'     \strong{Default="Metabolite"}
#' @param method_background \emph{Optional: } Background method `Intra |core, Intra&core, core,
#'     Intra or * \strong{Default="Intra&core"}
#' @param save_table \emph{Optional: } File types for the analysis results are: "csv",
#'     "xlsx", "txt" \strong{default: "csv"}
#'
#' @return Returns warnings and errors if input is not correct
#'
#' @importFrom logger log_trace
#' @noRd
check_param_mca <- function(
    data_c1,
    data_c2,
    data_core,
    data_intra,
    metadata_info_c1,
    metadata_info_c2,
    metadata_info_core,
    metadata_info_intra,
    feature= "Metabolite",
    method_background="Intra&core",
    save_table = "csv"
) {
    ## ------------ Create log file ----------- ##
    metaproviz_init()

    # ------------- data
    if (is.null(data_c1) == FALSE) {
        if (!inherits(data_c1, "data.frame") | !inherits(data_c2, "data.frame")) {
        message <-
            paste0(
                "data_c1 and data_c2 should be a data.frame. It's currently a ",
                paste(class(data_c1)),
                paste(class(data_c2)),
                ".",
                sep = ""
            )
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }
    if (length(data_c1[duplicated(data_c1[[feature]]), feature]) > 0) {
        message <- paste0("Duplicated features of data_c1, whilst features must be unique")
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }
    if (length(data_c2[duplicated(data_c2[[feature]]), feature]) > 0) {
        message <- paste0("Duplicated features of data_c2, whilst features must be unique")
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }

    } else {
    if (!inherits(data_intra, "data.frame") | !inherits(data_core, "data.frame")) {
        message <-
            paste0(
                "data_intra and data_core should be a data.frame. It's currently a ",
                paste(class(data_intra)),
                paste(class(data_core)),
                ".",
                sep = ""
            )
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }
    if (length(data_intra[duplicated(data_intra[[feature]]), feature]) > 0) {
        message <- paste0("Duplicated features of data_intra, whilst features must be unique")
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }
    if (length(data_core[duplicated(data_core[[feature]]), feature]) > 0) {
        message <- paste0("Duplicated features of data_core, whilst features must be unique")
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }
    }

    # ------------- metadata_info
    if (is.null(metadata_info_c1) == FALSE) {
    ## C1
    # ValueCol
    if ("ValueCol" %in% names(metadata_info_c1)) {
        if (metadata_info_c1[["ValueCol"]] %in% colnames(data_c1) == FALSE) {
            message <-
                paste0(
                    "The ",
                    metadata_info_c1[["ValueCol"]],
                    " column selected as ValueCol in metadata_info_c1 was not found in data_c1. Please check your input."
                )
            log_trace(paste("Error ", message, sep =""))
            stop(message)
        }
    }
    # StatCol
    if ("StatCol" %in% names(metadata_info_c1)) {
        if (metadata_info_c1[["StatCol"]] %in% colnames(data_c1) == FALSE) {
            message <-
                paste0(
                    "The ",
                    metadata_info_c1[["StatCol"]],
                    " column selected as StatCol in metadata_info_c1 was not found in data_c1. Please check your input."
                )
            log_trace(paste("Error ", message, sep =""))
            stop(message)
        }
    }

    ## C2
    # ValueCol
    if ("ValueCol" %in% names(metadata_info_c2)) {
        if (metadata_info_c2[["ValueCol"]] %in% colnames(data_c2) == FALSE) {
            message <-
                paste0(
                    "The ",
                    metadata_info_c2[["ValueCol"]],
                    " column selected as ValueCol in metadata_info_c2 was not found in data_c2. Please check your input."
                )
            log_trace(paste("Error ", message, sep =""))
            stop(message)
        }
    }
    # StatCol
    if ("StatCol" %in% names(metadata_info_c2)) {
        if (metadata_info_c2[["StatCol"]] %in% colnames(data_c2) == FALSE) {
            message <-
                paste0(
                    "The ",
                    metadata_info_c2[["StatCol"]],
                    " column selected as StatCol in metadata_info_c2 was not found in data_c2. Please check your input."
                )
            log_trace(paste("Error ", message, sep =""))
            stop(message)
        }
    }
    } else {
    ## Intra
    # ValueCol
    if ("ValueCol" %in% names(metadata_info_intra)) {
        if (metadata_info_intra[["ValueCol"]] %in% colnames(data_intra) == FALSE) {
            message <-
                paste0(
                    "The ",
                    metadata_info_intra[["ValueCol"]],
                    " column selected as ValueCol in metadata_info_intra was not found in data_intra. Please check your input."
                )
            log_trace(paste("Error ", message, sep =""))
            stop(message)
        }
    }
    # StatCol
    if ("StatCol" %in% names(metadata_info_intra)) {
        if (metadata_info_intra[["StatCol"]] %in% colnames(data_intra) == FALSE) {
            message <-
                paste0(
                    "The ",
                    metadata_info_intra[["StatCol"]],
                    " column selected as StatCol in metadata_info_intra was not found in data_intra. Please check your input."
                )
            log_trace(paste("Error ", message, sep =""))
            stop(message)
        }
    }

    ## core
    # ValueCol
    if ("ValueCol" %in% names(metadata_info_core)) {
        if (metadata_info_core[["ValueCol"]] %in% colnames(data_core) == FALSE) {
            message <-
                paste0(
                    "The ",
                    metadata_info_core[["ValueCol"]],
                    " column selected as ValueCol in metadata_info_core was not found in data_core. Please check your input."
                )
            log_trace(paste("Error ", message, sep =""))
            stop(message)
        }
    }
    # StatCol
    if ("StatCol" %in% names(metadata_info_core)) {
        if (metadata_info_core[["StatCol"]] %in% colnames(data_core) == FALSE) {
            message <-
                paste0(
                    "The ",
                    metadata_info_core[["StatCol"]],
                    " column selected as StatCol in metadata_info_core was not found in data_core. Please check your input."
                )
            log_trace(paste("Error ", message, sep =""))
            stop(message)
        }
    }

    # StatCol
    if ("DirectionCol" %in% names(metadata_info_core)) {
        if (metadata_info_core[["DirectionCol"]] %in% colnames(data_core) == FALSE) {
            message <-
                paste0(
                    "The ",
                    metadata_info_core[["DirectionCol"]],
                    " column selected as DirectionCol in metadata_info_core was not found in data_core. Please check your input."
                )
            log_trace(paste("Error ", message, sep =""))
            stop(message)
        }
    }

    }

    # ------------- metadata_info Cutoffs:
    if (is.null(metadata_info_c1) == FALSE) {
        if (is.na(as.numeric(metadata_info_c1[["cutoff_stat"]])) == TRUE | as.numeric(metadata_info_c1[["cutoff_stat"]]) > 1 | as.numeric(metadata_info_c1[["cutoff_stat"]]) < 0) {
        message <- paste0("Check input. The selected cutoff_stat in metadata_info_c1 should be numeric and between 0 and 1.")
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }

    if (is.na(as.numeric(metadata_info_c2[["cutoff_stat"]])) == TRUE | as.numeric(metadata_info_c2[["cutoff_stat"]]) > 1 | as.numeric(metadata_info_c2[["cutoff_stat"]]) < 0) {
        message <- paste0("Check input. The selected cutoff_stat in metadata_info_c2 should be numeric and between 0 and 1.")
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }

    if (is.na(as.numeric(metadata_info_c1[["ValueCutoff"]])) == TRUE) {
        message <- paste0("Check input. The selected ValueCutoff in metadata_info_c1 should be numeric and between 0 and 1.")
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }

    if (is.na(as.numeric(metadata_info_c2[["ValueCutoff"]])) == TRUE) {
        message <- paste0("Check input. The selected ValueCutoff in metadata_info_c2 should be numeric and between 0 and 1.")
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }

    } else {
    if (is.na(as.numeric(metadata_info_intra[["cutoff_stat"]])) == TRUE | as.numeric(metadata_info_intra[["cutoff_stat"]]) > 1 | as.numeric(metadata_info_intra[["cutoff_stat"]]) < 0) {
        message <- paste0("Check input. The selected cutoff_stat in metadata_info_intra should be numeric and between 0 and 1.")
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }

    if (is.na(as.numeric(metadata_info_core[["cutoff_stat"]])) == TRUE | as.numeric(metadata_info_core[["cutoff_stat"]]) > 1 | as.numeric(metadata_info_core[["cutoff_stat"]]) < 0) {
        message <- paste0("Check input. The selected cutoff_stat in metadata_info_core should be numeric and between 0 and 1.")
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }

    if (is.na(as.numeric(metadata_info_intra[["ValueCutoff"]])) == TRUE) {
        message <- paste0("Check input. The selected ValueCutoff in metadata_info_intra should be numeric and between 0 and 1.")
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }

    if (is.na(as.numeric(metadata_info_core[["ValueCutoff"]])) == TRUE) {
        message <- paste0("Check input. The selected ValueCutoff in metadata_info_core should be numeric and between 0 and 1.")
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }
    }

    # ------------ NAs in data
    if (is.null(data_c1) == FALSE) {
        if (nrow(data_c1[complete.cases(data_c1[[metadata_info_c1[["ValueCol"]]]], data_c1[[metadata_info_c1[["StatCol"]]]]), ]) < nrow(data_c1)) {
        message <-
            paste0(
                "data_c1 includes NAs in ",
                metadata_info_c1[["ValueCol"]],
                " and/or in ",
                metadata_info_c1[["StatCol"]],
                ". ",
                nrow(data_c1)- nrow(data_c1[complete.cases(data_c1[[metadata_info_c1[["ValueCol"]]]], data_c1[[metadata_info_c1[["StatCol"]]]]), ]),
                " metabolites containing NAs are removed."
            )
        log_trace(paste("Warning ", message, sep =""))
        warning(message)
    }

    if (nrow(data_c2[complete.cases(data_c2[[metadata_info_c2[["ValueCol"]]]], data_c2[[metadata_info_c2[["StatCol"]]]]), ]) < nrow(data_c2)) {
        message <-
            paste0(
                "data_c2 includes NAs in ",
                metadata_info_c2[["ValueCol"]],
                " and/or in",
                metadata_info_c2[["StatCol"]],
                ". ",
                nrow(data_c2)- nrow(data_c2[complete.cases(data_c2[[metadata_info_c2[["ValueCol"]]]], data_c2[[metadata_info_c2[["StatCol"]]]]), ]),
                " metabolites containing NAs are removed."
            )
        log_trace(paste("Warning ", message, sep =""))
        warning(message)
    }
    } else {
    if (nrow(data_intra[complete.cases(data_intra[[metadata_info_intra[["ValueCol"]]]], data_intra[[metadata_info_intra[["StatCol"]]]]), ]) < nrow(data_intra)) {
        message <-
            paste0(
                "data_intra includes NAs in ",
                metadata_info_intra[["ValueCol"]],
                " and/or in ",
                metadata_info_intra[["StatCol"]],
                ". ",
                nrow(data_intra)- nrow(data_intra[complete.cases(data_intra[[metadata_info_intra[["ValueCol"]]]], data_intra[[metadata_info_intra[["StatCol"]]]]), ]),
                " metabolites containing NAs are removed."
            )
        log_trace(paste("Warning ", message, sep =""))
        warning(message)
    }

    if (nrow(data_core[complete.cases(data_core[[metadata_info_core[["ValueCol"]]]], data_core[[metadata_info_core[["StatCol"]]]]), ]) < nrow(data_core)) {
        message <-
            paste0(
                "data_core includes NAs in ",
                metadata_info_core[["ValueCol"]],
                " and/or in ",
                metadata_info_core[["StatCol"]],
                ". ",
                nrow(data_core)- nrow(data_core[complete.cases(data_core[[metadata_info_core[["ValueCol"]]]], data_core[[metadata_info_core[["StatCol"]]]]), ]),
                " metabolites containing NAs are removed."
            )
        log_trace(paste("Warning ", message, sep =""))
        warning(message)
    }
    }

    # ------------- method_background
    if (is.null(metadata_info_c1) == FALSE) {
        options <- c("C1 | C2", "C1&C2", "C2", "C1" , "*")
    if (any(options %in% method_background) == FALSE) {
        message <-
            paste0(
                "Check input. The selected method_background option is not valid. Please select one of the folowwing: ",
                paste(options,collapse = ", "),
                "."
            )
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }
    } else {
    options <- c("Intra | core", "Intra&core", "core", "Intra" , "*")
    if (any(options %in% method_background) == FALSE) {
        message <-
            paste0(
                "Check input. The selected method_background option is not valid. Please select one of the folowwing: ",
                paste(options,collapse = ", "),
                "."
            )
        log_trace(paste("Error ", message, sep =""))
        stop(message)
    }
    }

    # ------------- SaveAs
    save_table_options <- c("txt","csv", "xlsx")
    if (is.null(save_table) == FALSE) {
        if ((save_table %in% save_table_options == FALSE) | (is.null(save_table) == TRUE)) {
        message <-
            paste0(
                "Check input. The selected save_table option is not valid. Please select one of the folowwing: ",
                paste(save_table_options,collapse = ", "),
                "."
            )
        log_trace(paste("Error ", message, sep =""))
        stop(  message)
    }
    }
}


