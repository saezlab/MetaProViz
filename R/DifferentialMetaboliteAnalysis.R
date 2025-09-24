#!/usr/bin/env Rscript

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




########################################################
### ### ### Differential Metabolite Analysis ### ### ###
########################################################

#' This function allows you to perform differential metabolite analysis to
#' obtain a Log2FC, pval, padj and tval comparing two or multiple conditions.
#'
#' @param data DF with unique sample identifiers as row names and metabolite
#'        numerical values in columns with metabolite identifiers as column
#'        names. Use NA for metabolites that were not detected.
#' @param metadata_sample DF which contains metadata information about the
#'        samples, which will be combined with your input data based on the
#'        unique sample identifiers used as rownames.
#' @param metadata_info  \emph{Optional: } Named vector including the
#'        information about the conditions column information on numerator or
#'        denominator c(Conditions="ColumnName_SettingsFile", Numerator =
#'        "ColumnName_SettingsFile", Denominator  = "ColumnName_SettingsFile").
#'        Denominator and Numerator will specify which comparison(s) will be
#'        done (one-vs-one, all-vs-one, all-vs-all), e.g. Denominator=NULL and
#'        Numerator =NULL selects all the condition and performs multiple
#'        comparison all-vs-all. Log2FC are obtained by dividing the numerator
#'        by the denominator, thus positive Log2FC values mean higher expression
#'        in the numerator. \strong{Default = c(conditions="Conditions",
#'        numerator = NULL, denumerator = NULL)}
#' @param pval \emph{Optional: } String which contains an abbreviation of the
#'        selected test to calculate p.value. For one-vs-one comparisons choose
#'        t.test, wilcox.test, "chisq.test", "cor.test" or lmFit (=limma), for
#'        one-vs-all or all-vs-all comparison choose aov (=anova), welch(=welch
#'        anova), kruskal.test or lmFit (=limma) \strong{Default = "lmFit"}
#' @param padj \emph{Optional: } String which contains an abbreviation of the
#'        selected p.adjusted test for p.value correction for multiple
#'        Hypothesis testing. Search: ?p.adjust for more methods:"BH", "fdr",
#'        "bonferroni", "holm", etc.\strong{Default = "fdr"}
#' @param metadata_feature \emph{Optional: } DF which contains the metadata
#'        information , i.e. pathway information, retention time,..., for each
#'        metabolite. The row names must match the metabolite names in the
#'        columns of the data. \strong{Default = NULL}
#' @param core \emph{Optional: } TRUE or FALSE for whether a Consumption/Release
#'        input is used. \strong{Default = FALSE}
#' @param vst TRUE or FALSE for whether to use variance stabilizing
#'        transformation on the data when linear modeling is used for hypothesis
#'        testing. \strong{Default = FALSE}
#' @param shapiro TRUE or FALSE for whether to perform the shapiro.test and get
#'        informed about data distribution (normal versus not-normal
#'        distribution. \strong{Default = TRUE}
#' @param bartlett TRUE or FALSE for whether to perform the bartlett.test.
#'        \strong{Default = TRUE}
#' @param transform TRUE or FALSE. If TRUE we expect the data to be not log2
#'        transformed and log2 transformation will be performed within the limma
#'        function and Log2FC calculation. If FALSE we expect the data to be
#'        log2 transformed as this impacts the Log2FC calculation and limma.
#'        \strong{Default= TRUE}
#' @param save_plot \emph{Optional: } Select the file type of output plots.
#'        Options are svg, png, pdf. \strong{Default = svg}
#' @param save_table \emph{Optional: } File types for the analysis results are:
#'        "csv", "xlsx", "txt". \strong{Default = "csv"}
#' @param print_plot \emph{Optional: } TRUE or FALSE, if TRUE Volcano plot is
#'        saved as an overview of the results. \strong{Default = TRUE}
#' @param path \emph{Optional:} Path to the folder the results should be saved
#'        at. \strong{Default = NULL}
#'
#' @return Dependent on parameter settings, list of lists will be returned for
#'         dma (DF of each comparison), shapiro (Includes DF and Plot), bartlett
#'         (Includes DF and Histogram), vst (Includes DF and Plot) and
#'         VolcanoPlot (Plots of each comparison).
#'
#' @examples
#' Intra <- intracell_raw[-c(49:58), ] %>% column_to_rownames("Code")
#' ResI <- dma(
#'     data = Intra[, -c(1:3)],
#'     metadata_sample = Intra[, c(1:3)],
#'     metadata_info = c(
#'         Conditions = "Conditions", Numerator = NULL, Denominator = "HK2"
#'     )
#' )
#'
#' @keywords Differential Metabolite Analysis, Multiple Hypothesis testing,
#'           Normality testing
#'
#' @importFrom magrittr %>%
#' @importFrom logger log_info
#' @importFrom dplyr full_join rename
#' @importFrom logger log_info
#' @importFrom purrr map reduce
#' @importFrom tibble column_to_rownames rownames_to_column
#'
#' @export
dma <- function(
    data,
    metadata_sample,
    metadata_info = c(
        Conditions = "Conditions",
        Numerator = NULL,
        Denominator = NULL
    ) ,
    pval = "lmFit",
    padj = "fdr",
    metadata_feature = NULL,
    core = FALSE,
    vst = FALSE,
    shapiro = TRUE,
    bartlett = TRUE,
    transform = TRUE,
    save_plot = "svg",
    save_table = "csv",
    print_plot = TRUE,
    path = NULL
    ) {

    ## ------------ Create log file ----------- ##
    metaproviz_init()

    log_info("dma: Differential metabolite analysis.")

    ## ------------ Check Input files ----------- ##
    # HelperFunction `check_param`
    check_param(
        data = data,
        metadata_sample = metadata_sample,
        metadata_feature = metadata_feature,
        metadata_info = metadata_info,
        save_plot = save_plot,
        save_table = save_table,
        core = core,
        print_plot = print_plot
    )

    # HelperFunction `check_param` Specific
    Settings <-
        check_param_dma(
            data = data,
            metadata_sample = metadata_sample,
            metadata_info = metadata_info,
            pval = pval,
            padj = padj,
            shapiro = shapiro,
            bartlett = bartlett,
            vst = vst,
            transform = transform
        )

    ## ------------ Create Results output folder ----------- ##
    if (is.null(save_plot) == FALSE | is.null(save_table) == FALSE) {
        folder <- save_path(
            folder_name = "dma",
            path = path
        )

        if (shapiro == TRUE) {
            Subfolder_S <- file.path(folder, "shapiro")
            if (!dir.exists(Subfolder_S)) {
                dir.create(Subfolder_S)
            }
        }

        if (bartlett == TRUE) {
          Subfolder_B <- file.path(folder, "bartlett")
            if (!dir.exists(Subfolder_B)) {
                dir.create(Subfolder_B)
            }
        }

        if (vst == TRUE) {
            Subfolder_V <- file.path(folder, "vst")
            if (!dir.exists(Subfolder_V)) {
                dir.create(Subfolder_V)
            }
        }
    }

    ############################################################################
    ## ------------ Check hypothesis test assumptions ----------- ##
    # 1. Normality
    if (shapiro == TRUE) {
        if (length(Settings[["Metabolites_Miss"]] >= 1)) {

          message(
              paste0(
                "There are NA's/0s in the data. This can impact the output of ",
                "the SHapiro-Wilk test for all metabolites that include NAs/0s."
              )
          )
        }
        tryCatch(
            {
            Shapiro_output <-
              suppressWarnings(
                  shapiro(
                      data = data,
                      metadata_sample = metadata_sample,
                      metadata_info = metadata_info,
                      pval = pval,
                      qqplots = FALSE
                  )
              )
            },
            error = function(e) {
                message(
                    paste0(
                      "Error occurred during shapiro that performs ",
                      "the shapiro-Wilk test. Message:"
                    ),
                    conditionMessage(e)
                )
            }
        )
    }

    # 2. Variance homogeneity
    if (bartlett == TRUE) {
        # if we only have two conditions, which can happen even tough multiple
        # comparison (C1 versus C2 and C2 versus C1 is done)
        if (Settings[["MultipleComparison"]] == TRUE) {
            UniqueConditions <-
                metadata_sample %>%
                    subset(
                        metadata_sample[[metadata_info[["Conditions"]]]] %in% Settings[["numerator"]] | metadata_sample[[metadata_info[["Conditions"]]]] %in% Settings[["denominator"]],
                        select = c(metadata_info[["Conditions"]])
                    )
            UniqueConditions <-
                unique(UniqueConditions[[metadata_info[["Conditions"]]]])

            if (length(UniqueConditions) > 2) {
                tryCatch(
                    {
                    Bartlett_output <-
                        suppressWarnings(
                              bartlett(
                                  data = data,
                                  metadata_sample = metadata_sample,
                                  metadata_info = metadata_info
                              )
                        )
                    },
                    error = function(e) {
                        message(
                            paste0(
                                "Error occurred during bartlett that performs ",
                                "the bartlett test. Message:"
                            ),
                            conditionMessage(e)
                        )
                    }
                )
            }
        }
    }


    ############################################################################
    #### Prepare the data ######
    # 1. Metabolite names:
    savedMetaboliteNames <- data.frame("InputName" = colnames(data))
    savedMetaboliteNames$Metabolite <- paste0("M", seq(1, length(colnames(data))))
    colnames(data) <- savedMetaboliteNames$Metabolite

    ############################################################################
    ############### Calculate Log2FC, pval, padj, tval and add additional info
    ############### ###############
    log2fc_table <- log2fc(
        data = data,
        metadata_sample = metadata_sample,
        metadata_info = metadata_info,
        core = core,
        transform = transform
    )

    ############################################################################
    ############### Perform Hypothesis testing ###############
    if (Settings[["MultipleComparison"]] == FALSE) {
        if (pval == "lmFit") {
            STAT_C1vC2 <-
                dma_stat_limma(
                    data = data,
                    metadata_sample = metadata_sample,
                    metadata_info = metadata_info,
                    padj = padj,
                    log2fc_table = log2fc_table,
                    core = core,
                    transform = transform
                )
        } else {
            STAT_C1vC2 <-
                dma_stat_single(
                    data = data,
                    metadata_sample = metadata_sample,
                    metadata_info = metadata_info,
                    log2fc_table = log2fc_table,
                    pval = pval,
                    padj = padj
                )
        }
    } else { # MultipleComparison = TRUE
        # Correct data heteroscedasticity
        if (pval != "lmFit" & vst == TRUE) {
            vst_res <- vst(data)
            data <- vst_res[["DFs"]][["Corrected_data"]]
        }

        if (Settings[["all_vs_all"]] == TRUE) {
            message(
                paste0(
                    "No conditions were specified as numerator or denumerator. ",
                    "Performing multiple testing `all-vs-all` using"
                ),
                paste(pval),
                "."
            )
        } else { # for 1 vs all
            message(
                "No condition was specified as numerator and ",
                Settings[["denominator"]],
                paste0(
                    "was selected as a denominator. Performing ",
                    "multiple testing `all-vs-one` using"
                ),
                paste(pval),
                "."
            )
        }

        if (pval == "aov") {
            STAT_C1vC2 <-
                aov(
                    data = data,
                    metadata_sample = metadata_sample,
                    metadata_info = metadata_info,
                    log2fc_table = log2fc_table
                )
        } else if (pval == "kruskal.test") {
            STAT_C1vC2 <-
                kruskal(
                    data = data,
                    metadata_sample = metadata_sample,
                    metadata_info = metadata_info,
                    log2fc_table = log2fc_table,
                    padj = padj
                )
        } else if (pval == "welch") {
            STAT_C1vC2 <-
                welch(
                    data = data,
                    metadata_sample = metadata_sample,
                    metadata_info = metadata_info,
                    log2fc_table = log2fc_table
                )
        } else if (pval == "lmFit") {
            STAT_C1vC2 <- dma_stat_limma(
            data = data,
            metadata_sample = metadata_sample,
            metadata_info = metadata_info,
            padj = padj,
            log2fc_table = log2fc_table,
            core = core,
            transform = transform
            )
        }
    }


    ############################################################################
    ###############  Add the previous metabolite names back ###############
    DMA_Output <-
        lapply(
            STAT_C1vC2, function(df) {
                merged_df <-
                    merge(
                        savedMetaboliteNames,
                        df,
                        by = "Metabolite",
                        all.y = TRUE
                    )
                # remove the names we used as part of the function and
                # add back the input names.
                merged_df <- merged_df[, -1] %>%
                    rename("Metabolite" = 1)
                return(merged_df)
            }
        )

    ############################################################################
    ###############  Add the metabolite Metadata if available ###############
    if (is.null(metadata_feature) == FALSE) {
        DMA_Output <-
            lapply(
                DMA_Output, function(df) {
                    merged_df <-
                        merge(
                            df,
                            metadata_feature %>% rownames_to_column("Metabolite"),
                            by = "Metabolite",
                            all.x = TRUE
                        )
                    return(merged_df)
                }
            )
    }

    ############################################################################
    ###############  For core=TRUE create summary of Feature_metadata  #########
    if (core == TRUE) {
        df_list_selected <-
            map(
                names(DMA_Output), function(df_name) {
                    df <- DMA_Output[[df_name]] # Extract the dataframe

                    # Extract the dynamic column name
                    # Find the column that starts with "core_"
                    core_col <- grep("^core_", names(df), value = TRUE)
                    # Filter only columns where the part after "core_"
                    # is in valid_conditions
                    core_col <-
                        core_col[str_remove(core_col, "^core_") %in%
                            unique(
                                metadata_sample[[metadata_info[["Conditions"]]]]
                                )
                        ]

                    # Select only the relevant columns
                    df_selected <- df %>%
                        select(Metabolite, all_of(core_col))

                    return(df_selected)
                }
            )

        # Merge all dataframes by "Metabolite"
        merged_df <-
            reduce(
                df_list_selected,
                full_join,
                by = "Metabolite"
            )
        # It is likely we have duplications that cause .x, .y, .x.x, .y.y, etc. to
        # be added to the column names. We only keep one column (.x)
        names(merged_df) <- gsub("\\.x$", "", names(merged_df))

        Feature_Metadata <- merged_df %>%
            # Now we remove all other columns with .x.x, .y.y, etc.
            select(-all_of(grep("\\.[xy]+$", names(merged_df), value = TRUE)))

        if (is.null(metadata_feature) == FALSE) { # Add to Metadata file:
            Feature_Metadata <-
                merge(
                    metadata_feature %>% rownames_to_column("Metabolite"),
                    Feature_Metadata,
                    by = "Metabolite",
                    all.x = TRUE
                )
        }
    }

    ############################################################################
    ###############  Plots ###############
    if (core == TRUE) {
        x <- "Log2(Distance)"
        VolPlot_metadata_info <- c(color = "core")
        VolPlot_SettingsFile <- DMA_Output
    } else {
        x <- "Log2FC"
        VolPlot_metadata_info <- NULL
        VolPlot_SettingsFile <- NULL
    }

    volplotList <- list()
    for (DF in names(DMA_Output)) { # DF = names(DMA_Output)[2]
        Volplotdata <- DMA_Output[[DF]]

        if (core == TRUE) {
            VolPlot_SettingsFile <-
                DMA_Output[[DF]] %>%
                    column_to_rownames("Metabolite")
        }

        dev.new()
        VolcanoPlot <- invisible(
            viz_volcano(
                plot_types = "Standard",
                data = Volplotdata %>% column_to_rownames("Metabolite"),
                metadata_info = VolPlot_metadata_info,
                metadata_feature = VolPlot_SettingsFile,
                y = "p.adj",
                x = x,
                plot_name = DF,
                subtitle = bquote(italic("Differential Metabolite Analysis")),
                save_plot = NULL
            )
        )

        ## Remove special characters and replace spaces with underscores
        DF_save <- gsub("[^A-Za-z0-9._-]", "_", DF)
        volplotList[[DF_save]] <- VolcanoPlot[["Plot_Sized"]][[1]]

        dev.off()
    }

    ############################################################################
    ## ----- Save and Return
    DMA_Output_List <- list()
    # Here we make a list in which we will save the outputs:
    if (shapiro == TRUE & exists("Shapiro_output") == TRUE) {
        suppressMessages(
            suppressWarnings(
                save_res(
                    inputlist_df = Shapiro_output[["DF"]],
                    inputlist_plot = Shapiro_output[["Plot"]][["Distributions"]],
                    save_table = save_table,
                    save_plot = save_plot,
                    path = Subfolder_S,
                    file_name = "ShapiroTest",
                    core = core,
                    print_plot = print_plot
                )
            )
    )

      DMA_Output_List <- list("ShapiroTest" = Shapiro_output)
    }

    if (bartlett == TRUE & exists("Bartlett_output") == TRUE) {
        suppressMessages(
            suppressWarnings(
                save_res(
                    inputlist_df = Bartlett_output[["DF"]],
                    inputlist_plot = Bartlett_output[["Plot"]],
                    save_table = save_table,
                    save_plot = save_plot,
                    path = Subfolder_B,
                    file_name = "BartlettTest",
                    core = core,
                    print_plot = print_plot
                )
            )
        )

      DMA_Output_List <-
        c(
            DMA_Output_List,
            list("BartlettTest" = Bartlett_output)
        )
    }

    if (vst == TRUE & exists("vst_res") == TRUE) {
        suppressMessages(
            suppressWarnings(
                save_res(
                    inputlist_df = vst_res[["DF"]],
                    inputlist_plot = vst_res[["Plot"]],
                    save_table = save_table,
                    save_plot = save_plot,
                    path = Subfolder_V,
                    file_name = "vst_res",
                    core = core,
                    print_plot = print_plot
                )
            )
        )

      DMA_Output_List <- c(DMA_Output_List, list("vstres" = Bartlett_output))
    }

    if (core == TRUE) {
        suppressMessages(
            suppressWarnings(
                save_res(
                    inputlist_df = list("Feature_Metadata" = Feature_Metadata),
                    inputlist_plot = NULL,
                    save_table = save_table,
                    save_plot = NULL,
                    path = folder,
                    file_name = "dma",
                    core = core,
                    print_plot = print_plot
                )
            )
        )

        DMA_Output_List <-
            c(
                DMA_Output_List,
                list("Feature_Metadata" = Feature_Metadata)
            )
    }

    suppressMessages(
        suppressWarnings(
            save_res(
                # This needs to be a list, also for single comparisons
                inputlist_df = DMA_Output,
                inputlist_plot = volplotList,
                save_table = save_table,
                save_plot = save_plot,
                path = folder,
                file_name = "dma",
                core = core,
                print_plot = print_plot
            )
        )
    )

    DMA_Output_List <-
        c(
            DMA_Output_List,
            list("dma" = DMA_Output, "VolcanoPlot" = volplotList)
        )

    return(invisible(DMA_Output_List))
}


###############################
### ### ### Log2FC  ### ### ###
###############################

#' This helper function calculates the Log2(FoldChange) or in case of
#' core Log2(Distance).
#'
#' @param data DF with unique sample identifiers as row names and metabolite
#'        numerical values in columns with metabolite identifiers as column
#'        names. Use NA for metabolites that were not detected.
#' @param metadata_sample DF which contains metadata information about the
#'        samples, which will be combined with your input data based on the
#'        unique sample identifiers used as rownames.
#' @param metadata_info \emph{Optional: } Named vector including the information
#'        about the conditions column information on numerator or denominator
#'        c(Conditions="ColumnName_SettingsFile", Numerator =
#'        "ColumnName_SettingsFile", Denominator  = "ColumnName_SettingsFile").
#'        Denominator and Numerator will specify which comparison(s) will be
#'        done (one-vs-one, all-vs-one, all-vs-all), e.g. Denominator=NULL and
#'        Numerator =NULL selects all the condition and performs multiple
#'        comparison all-vs-all. Log2FC are obtained by dividing the numerator
#'        by the denominator, thus positive Log2FC values mean higher expression
#'        in the numerator. \strong{Default = c(conditions="Conditions",
#'        numerator = NULL, denumerator = NULL)}
#' @param core \emph{Optional: } TRUE or FALSE for whether a Consumption/Release
#'        input is used \strong{default = FALSE}
#' @param transform \emph{Optional: } If TRUE we expect the data to be not log2
#'        transformed and log2 transformation will be performed within the limma
#'        function and Log2FC calculation. If FALSE we expect the data to be
#'        log2 transformed as this impacts the Log2FC calculation and
#'        limma.\strong{default = TRUE}
#'
#' @return List of DFs named after comparison (e.g. Tumour versus Normal) with
#'         Log2FC or Log2(Distance) column and column with feature names
#'
#' @keywords Log2FC, core, Distance
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter mutate rename select_if summarise_all
#' @importFrom gtools foldchange2logratio
#' @importFrom tibble rownames_to_column
#'
#' @noRd
#'
log2fc <- function(
    data,
    metadata_sample,
    metadata_info = c(
        Conditions = "Conditions",
        Numerator = NULL,
        Denominator = NULL
    ),
    core = FALSE,
    transform = TRUE
    ) {
    ## ------------ Create log file ----------- ##
    metaproviz_init()

    # ------------ Assignments ----------- ##
    if ("Denominator" %in% names(metadata_info) == FALSE & "Numerator" %in% names(metadata_info) == FALSE) {
        # all-vs-all: Generate all pairwise combinations
        conditions <- metadata_sample[[metadata_info[["Conditions"]]]]
        denominator <- unique(metadata_sample[[metadata_info[["Conditions"]]]])
        numerator <- unique(metadata_sample[[metadata_info[["Conditions"]]]])
        comparisons <- combn(unique(conditions), 2) %>% as.matrix()
        # Settings:
        MultipleComparison <- TRUE
        all_vs_all <- TRUE
    } else if ("Denominator" %in% names(metadata_info) == TRUE & "Numerator" %in% names(metadata_info) == FALSE) {
        # all-vs-one: Generate the pairwise combinations
        conditions <- metadata_sample[[metadata_info[["Conditions"]]]]
        denominator <- metadata_info[["Denominator"]]
        numerator <- unique(metadata_sample[[metadata_info[["Conditions"]]]])
        # Remove denom from num
        numerator <- numerator[!numerator %in% denominator]
        comparisons <- t(expand.grid(numerator, denominator)) %>% as.data.frame()
        # Settings:
        MultipleComparison <- TRUE
        all_vs_all <- FALSE
    } else if ("Denominator" %in% names(metadata_info) == TRUE & "Numerator" %in% names(metadata_info) == TRUE) {
        # one-vs-one: Generate the comparisons
        denominator <- metadata_info[["Denominator"]]
        numerator <- metadata_info[["Numerator"]]
        comparisons <- matrix(c(numerator, denominator))
        # Settings:
        MultipleComparison <- FALSE
        all_vs_all <- FALSE
    }

    ## ------------ Check Missingness ------------- ##
    Num <- data %>%
        filter(
            metadata_sample[[metadata_info[["Conditions"]]]] %in% numerator
            ) %>%
            select_if(is.numeric)
    Denom <- data %>%
        filter(
            metadata_sample[[metadata_info[["Conditions"]]]] %in% denominator
        ) %>%
        select_if(is.numeric)

    Num_Miss <- Num[, colSums(Num == 0) > 0, drop = FALSE]
    Denom_Miss <- Denom[, colSums(Denom == 0) > 0, drop = FALSE]
    Metabolites_Miss <- unique(c(colnames(Num_Miss), colnames(Denom_Miss)))

    ## ------------ Denominator/numerator ----------- ##
    # Denominator and numerator: Define if we compare one_vs_one, one_vs_all or
    # all_vs_all.
    if ("Denominator" %in% names(metadata_info) == FALSE & "Numerator" %in% names(metadata_info) == FALSE) {
        MultipleComparison <- TRUE
    } else if ("Denominator" %in% names(metadata_info) == TRUE & "Numerator" %in% names(metadata_info) == FALSE) {
        MultipleComparison <- TRUE
    } else if ("Denominator" %in% names(metadata_info) == TRUE & "Numerator" %in% names(metadata_info) == TRUE) {
        MultipleComparison <- FALSE
    }

    ############################################################################
    ## ----------------- Log2FC ----------------------------
    log2fc_table <- list() # Create an empty list to store results data frames
    for (column in 1:dim(comparisons)[2]) {
        # Numerator
        C1 <- data %>%
            filter(
                metadata_sample[[metadata_info[["Conditions"]]]] %in% comparisons[1, column]
            ) %>%
            # only keep numeric columns with metabolite values
            select_if(is.numeric)
        # Deniminator
        C2 <- data %>%
            filter(
                metadata_sample[[metadata_info[["Conditions"]]]] %in% comparisons[2, column]
            ) %>%
            select_if(is.numeric)

        ## ------------  Calculate Log2FC ----------- ##
        # For C1_Mean and C2_Mean use 0 to obtain values, leading to Log2FC=NA if
        # mean = 0 (If one value is NA, the mean will be NA even though all other
        # values are available.)
        C1_Zero <- C1
        C1_Zero[is.na(C1_Zero)] <- 0
        Mean_C1 <- C1_Zero %>%
        summarise_all("mean")

        C2_Zero <- C2
        C2_Zero[is.na(C2_Zero)] <- 0
        Mean_C2 <- C2_Zero %>%
        summarise_all("mean")

    # Calculate absolute distance between the means. log2 transform
        # and add sign (-/+):
        if (core == TRUE) {
            # core values can be negative and positive, which can does not allow us to
            # calculate a Log2FC.
            Mean_C1_t <- as.data.frame(t(Mean_C1)) %>%
                rownames_to_column("Metabolite")
            Mean_C2_t <- as.data.frame(t(Mean_C2)) %>%
                rownames_to_column("Metabolite")
            Mean_Merge <- merge(
                Mean_C1_t,
                Mean_C2_t,
                by = "Metabolite",
                all = TRUE
            ) %>%
                rename(
                "C1" = 2,
                "C2" = 3
                )

            # Deal with NA/0s
            # Column to enable the check if mean values of 0 are due to missing
            # values (NA/0) and not by coincidence
            Mean_Merge$`NA/0` <- Mean_Merge$Metabolite %in% Metabolites_Miss

            if (any(
                (Mean_Merge$`NA/0` == FALSE & Mean_Merge$C1 == 0)  |
                (Mean_Merge$`NA/0` == FALSE & Mean_Merge$C2 == 0)
            ) == TRUE) {
                Mean_Merge <- Mean_Merge %>%
                    mutate(C1 = case_when(
                        # Here we have a "true" 0 value due to 0/NAs
                        # in the input data
                        C2 == 0 & `NA/0` == TRUE ~ paste(C1),
                        # Here we have a "true" 0 value due to 0/NAs
                        # in the input data
                        C1 == 0 & `NA/0` == TRUE ~ paste(C1),
                        # Here we have a "false" 0 value that occured at random
                        # and not due to 0/NAs in the input data, hence we add
                        # the constant +1
                        C2 == 0 & `NA/0` == FALSE ~ paste(C1 + 1),
                        # Here we have a "false" 0 value that occured at random
                        # and not due to 0/NAs in the input data, hence we add
                        # the constant +1
                        C1 == 0 & `NA/0` == FALSE ~ paste(C1 + 1),
                        TRUE ~ paste(C1)
                    )) %>%
                    mutate(C2 = case_when(
                        # Here we have a "true" 0 value due to 0/NAs in the
                        # input data
                        C1 == 0 & `NA/0` == TRUE ~ paste(C2),
                        # Here we have a "true" 0 value due to 0/NAs in the
                        # input data
                        C2 == 0 & `NA/0` == TRUE ~ paste(C2),
                        # Here we have a "false" 0 value that occured at random
                        # and not due to 0/NAs in the input data, hence we
                        # add the constant +1
                        C1 == 0 & `NA/0` == FALSE ~ paste(C2 + 1),
                        # Here we have a "false" 0 value that occured at random
                        # and not due to 0/NAs in the input data, hence we
                        # add the constant +1
                        C2 == 0 & `NA/0` == FALSE ~ paste(C2 + 1),
                        TRUE ~ paste(C2)
                    )) %>%
                    mutate(C1 = as.numeric(C1)) %>%
                    mutate(C2 = as.numeric(C2))

                X <- Mean_Merge %>%
                    subset(
                        (Mean_Merge$`NA/0` == FALSE & Mean_Merge$C1 == 0)  |
                        (Mean_Merge$`NA/0` == FALSE & Mean_Merge$C2 == 0)
                    )
                    message(
                        "We added +1 to the mean value of metabolite(s) ",
                        paste0(X$Metabolite, collapse = ", "),
                        paste0(
                            ", since the mean of the replicate values where 0. ",
                            "This was not due to missing values (NA/0)."
                        )
                    )
            }

            # Add the distance column:
            Mean_Merge$`Log2(Distance)` <-
                log2(abs(Mean_Merge$C1 - Mean_Merge$C2))

            # Now we can adapt the values to take into account the distance
            Mean_Merge <- Mean_Merge %>%
                mutate(`Log2(Distance)` = case_when(
                    # If C1>C2 the distance stays positive to reflect
                    # that C1 > C2
                    C1 > C2 ~ paste(`Log2(Distance)` * +1),
                    # If C1<C2 the distance gets a negative sign to reflect
                    # that C1 < C2
                    C1 < C2 ~ paste(`Log2(Distance)` * -1),
                    TRUE ~ "NA"
                    )
                ) %>%
                mutate(`Log2(Distance)` = as.numeric(`Log2(Distance)`))

            # Add additional information:
            temp1 <- Mean_C1
            temp2 <- Mean_C2
            # Add Info of core:
            core_info <- rbind(temp1, temp2, rep(0, length(temp1)))
            for (i in 1:length(temp1)) {
                if (temp1[i] > 0 & temp2[i] > 0) {
                    core_info[3, i] <- "Released"
                } else if (temp1[i] < 0 & temp2[i] < 0) {
                    core_info[3, i] <- "Consumed"
                } else if (temp1[i] > 0 & temp2[i] < 0) {
                    core_info[3, i] <-
                        paste(
                            "Released in",
                            comparisons[1,
                            column],
                            "and Consumed",
                            comparisons[2,
                            column],
                            sep = " "
                        )
                } else if (temp1[i] < 0 & temp2[i] > 0) {
                    core_info[3, i] <-
                        paste(
                            "Consumed in",
                            comparisons[1,
                            column],
                            " and Released",
                            comparisons[2,
                            column],
                            sep = " "
                        )
                } else {
                    core_info[3, i] <- "No Change"
                }
            }

            core_info <- t(core_info) %>% as.data.frame()
            core_info <- rownames_to_column(core_info, "Metabolite")
            names(core_info)[2] <-
                paste(
                    "Mean", comparisons[1, column], sep = "_"
                    )
            names(core_info)[3] <-
                paste(
                    "Mean", comparisons[2, column], sep = "_"
                )
            names(core_info)[4] <- "core_specific"

            core_info <- core_info %>%
                mutate(
                    core = case_when(
                        core_specific == "Released" ~ "Released",
                        core_specific == "Consumed" ~ "Consumed",
                        TRUE ~ "Released/Consumed"
                    )
                ) %>%
                mutate(!!paste(
                    "core_", comparisons[1, column], sep = ""
                    ) := case_when(
                        core_specific == "Released" ~ "Released",
                        core_specific == "Consumed" ~ "Consumed",
                        core_specific == paste(
                            "Consumed in",
                            comparisons[1, column],
                            " and Released",
                            comparisons[2, column],
                            sep = " "
                        ) ~ "Consumed",
                        core_specific == paste(
                            "Released in",
                            comparisons[1, column],
                            "and Consumed",
                            comparisons[2, column],
                            sep = " "
                        ) ~ "Released",
                        TRUE ~ "NA"
                )
                ) %>%
                mutate(!!paste(
                    "core_",
                    comparisons[2, column],
                    sep = ""
                ) := case_when(
                    core_specific == "Released" ~ "Released",
                    core_specific == "Consumed" ~ "Consumed",
                    core_specific == paste(
                        "Consumed in",
                        comparisons[1, column],
                        " and Released",
                        comparisons[2, column],
                        sep = " "
                    ) ~ "Released",
                    core_specific == paste(
                        "Released in",
                        comparisons[1, column],
                        "and Consumed",
                        comparisons[2, column],
                        sep = " "
                    ) ~ "Consumed",
                    TRUE ~ "NA"
                    )
                )


            Log2FC_C1vC2 <- merge(
                Mean_Merge[, c(1, 5)],
                core_info[, c(1, 2, 6, 3, 7, 4:5)],
                by = "Metabolite",
                all.x = TRUE
            )

            # Add info on Input:
            temp3 <- as.data.frame(t(C1)) %>%
                rownames_to_column("Metabolite")
            temp4 <- as.data.frame(t(C2)) %>%
                rownames_to_column("Metabolite")
            temp_3a4 <- merge(temp3, temp4, by = "Metabolite", all = TRUE)
            Log2FC_C1vC2 <- merge(
                Log2FC_C1vC2,
                temp_3a4,
                by = "Metabolite",
                all.x = TRUE
            )

            # Return DFs
            ## Make reverse DF
            Log2FC_C2vC1 <- Log2FC_C1vC2
            Log2FC_C2vC1$`Log2(Distance)` <- Log2FC_C2vC1$`Log2(Distance)` * -1

            ## Name them
            if (MultipleComparison == TRUE) {
                logname <- paste(
                comparisons[1, column],
                comparisons[2, column],
                sep = "_vs_"
                )
                logname_reverse <- paste(
                comparisons[2, column],
                comparisons[1, column],
                sep = "_vs_"
                )

                # Store the data frame in the results list,
                # named after the contrast
                log2fc_table[[logname]] <- Log2FC_C1vC2
                log2fc_table[[logname_reverse]] <- Log2FC_C2vC1
            } else {
                log2fc_table <- Log2FC_C1vC2
            }
        } else if (core == FALSE) {
            # Mean values could be 0, which can not be used to calculate a
            # Log2FC and hence the Log2FC(A versus B)=(log2(A+x)-log2(B+x))
            # for A and/or B being 0, with x being set to 1
            Mean_C1_t <- as.data.frame(t(Mean_C1)) %>%
                rownames_to_column("Metabolite")
            Mean_C2_t <- as.data.frame(t(Mean_C2)) %>%
                rownames_to_column("Metabolite")
            Mean_Merge <-
                merge(
                    Mean_C1_t,
                    Mean_C2_t,
                    by = "Metabolite",
                    all = TRUE
                ) %>%
                rename(
                    "C1" = 2,
                    "C2" = 3
                )
            # Column to enable the check if mean values of 0 are due to missing
            # values (NA/0) and not by coincidence
            Mean_Merge$`NA/0` <- Mean_Merge$Metabolite %in% Metabolites_Miss

            Mean_Merge <- Mean_Merge %>%
                mutate(C1_Adapted = case_when(
                    # Here we have a "true" 0 value due to
                    # 0/NAs in the input data
                    C2 == 0 & `NA/0` == TRUE ~ paste(C1),
                    # Here we have a "true" 0 value due to
                    # 0/NAs in the input data
                    C1 == 0 & `NA/0` == TRUE ~ paste(C1),
                    # Here we have a "false" 0 value that occured at random
                    # and not due to 0/NAs in the input data,
                    # hence we add the constant +1
                    C2 == 0 & `NA/0` == FALSE ~ paste(C1 + 1),
                    # Here we have a "false" 0 value that occured at random
                    # and not due to 0/NAs in the input data,
                    # hence we add the constant +1
                    C1 == 0 & `NA/0` == FALSE ~ paste(C1 + 1),
                    TRUE ~ paste(C1)
                    )
                ) %>%
                mutate(C2_Adapted = case_when(
                    # Here we have a "true" 0 value due to
                    # 0/NAs in the input data
                    C1 == 0 & `NA/0` == TRUE ~ paste(C2),
                    # Here we have a "true" 0 value due to
                    # 0/NAs in the input data
                    C2 == 0 & `NA/0` == TRUE ~ paste(C2),
                    # Here we have a "false" 0 value that occured at random
                    # and not due to 0/NAs in the input data,
                    # hence we add the constant +1
                    C1 == 0 & `NA/0` == FALSE ~ paste(C2 + 1),
                    # Here we have a "false" 0 value that occured at random
                    # and not due to 0/NAs in the input data,
                    # hence we add the constant +1
                    C2 == 0 & `NA/0` == FALSE ~ paste(C2 + 1),
                    TRUE ~ paste(C2)
                    )
                ) %>%
                mutate(C1_Adapted = as.numeric(C1_Adapted)) %>%
                    mutate(C2_Adapted = as.numeric(C2_Adapted))

            if (any(
                (Mean_Merge$`NA/0` == FALSE & Mean_Merge$C1 == 0)  |
                (Mean_Merge$`NA/0` == FALSE & Mean_Merge$C2 == 0)
            ) == TRUE) {
                X <- Mean_Merge %>%
                subset(
                    (Mean_Merge$`NA/0` == FALSE & Mean_Merge$C1 == 0)  |
                    (Mean_Merge$`NA/0` == FALSE & Mean_Merge$C2 == 0)
                )
                message(
                    "We added +1 to the mean value of metabolite(s) ",
                    paste0(X$Metabolite, collapse = ", "),
                    paste0(
                        ", since the mean of the replicate values where 0. ",
                        "This was not due to missing values (NA/0)."
                    )
                )
            }

            # Calculate the Log2FC
            if (transform == TRUE) { # data are not log2 transformed
                # FoldChange
                Mean_Merge$FC_C1vC2 <-
                    Mean_Merge$C1_Adapted / Mean_Merge$C2_Adapted
                Mean_Merge$Log2FC <-
                    foldchange2logratio(
                        Mean_Merge$FC_C1vC2,
                        base = 2
                    )
            }

            # data has been log2 transformed and hence we need to take this into
            # account when calculating the log2FC
            if (transform == FALSE) {
                Mean_Merge$FC_C1vC2 <- "Empty"
                Mean_Merge$Log2FC <-
                    Mean_Merge$C1_Adapted - Mean_Merge$C2_Adapted
            }

            # Add info on Input:
            temp3 <- as.data.frame(t(C1)) %>%
                rownames_to_column("Metabolite")
            temp4 <- as.data.frame(t(C2)) %>%
                rownames_to_column("Metabolite")
            temp_3a4 <- merge(temp3, temp4, by = "Metabolite", all = TRUE)
            Log2FC_C1vC2 <- merge(
                Mean_Merge[,
                c(1, 8)],
                temp_3a4,
                by = "Metabolite",
                all.x = TRUE
            )

            # Return DFs
            ## Make reverse DF
            Log2FC_C2vC1 <- Log2FC_C1vC2
            Log2FC_C2vC1$Log2FC <- Log2FC_C2vC1$Log2FC * -1

            if (MultipleComparison == TRUE) {
                logname <- paste(
                    comparisons[1, column],
                    comparisons[2, column],
                    sep = "_vs_"
                )
                logname_reverse <- paste(
                    comparisons[2, column],
                    comparisons[1, column],
                    sep = "_vs_"
                )

                # Store the data frame in the results list, named after the contrast
                log2fc_table[[logname]] <- Log2FC_C1vC2
                log2fc_table[[logname_reverse]] <- Log2FC_C2vC1
            } else {
                log2fc_table <- Log2FC_C1vC2
            }
        }
    }
    return(invisible(log2fc_table))
}


################################################################################
### ### dma helper function: Internal Function to perform single comparison  ###
################################################################################

#' Calculate One-vs-One comparison statistics
#'
#' @param data DF with unique sample identifiers as row names and metabolite
#'        numerical values in columns with metabolite identifiers as column
#'        names. Use NA for metabolites that were not detected.
#' @param metadata_sample DF which contains metadata information about the
#'        samples, which will be combined with your input data based on the
#'        unique sample identifiers used as rownames.
#' @param metadata_info  Named vector including the information about the
#'        conditions column information on numerator or denominator
#'        c(Conditions="ColumnName_SettingsFile", Numerator =
#'        "ColumnName_SettingsFile", Denominator  = "ColumnName_SettingsFile").
#'        Denominator and Numerator will specify which comparison(s) will be
#'        done (here one-vs-one).
#' @param log2fc_table \emph{Optional: } This is a List of DFs including a
#'        column "MetaboliteID" and Log2FC or Log2(Distance). This is the output
#'        from MetaProViz:::log2fc. If NULL, the output statistics will not be
#'        added into the Log2FC/Log2(Distance) DFs. \strong{Default = NULL}
#' @param pval \emph{Optional: } String which contains an abbreviation of the
#'        selected test to calculate p.value. For one-vs-one comparisons choose
#'        t.test, wilcox.test, "chisq.test" or "cor.test", \strong{Default =
#'        "t.test"}
#' @param padj \emph{Optional: } String which contains an abbreviation of the
#'        selected p.adjusted test for p.value correction for multiple
#'        Hypothesis testing. Search: ?p.adjust for more methods:"BH", "fdr",
#'        "bonferroni", "holm", etc.\strong{Default = "fdr"}
#'
#' @return List of DFs named after comparison (e.g. tumour versus Normal) with
#'         p-value, t-value and adjusted p-value column and column with feature
#'         names
#'
#' @keywords Statistical testing, p-value, t-value
#'
#' @importFrom stats p.adjust
#' @importFrom dplyr summarise_all filter mutate rename select_if
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @noRd
dma_stat_single <- function(data,
                            metadata_sample,
                            metadata_info,
                            log2fc_table=NULL,
                            pval="t.test",
                            padj="fdr"){
  ## ------------ Create log file ----------- ##
  metaproviz_init()

  ## ------------ Check Missingness ------------- ##
  Num <- data %>%
    filter(metadata_sample[[metadata_info[["Conditions"]]]] %in% metadata_info[["Numerator"]]) %>%
    select_if(is.numeric)
  Denom <- data %>%
    filter(metadata_sample[[metadata_info[["Conditions"]]]] %in% metadata_info[["Denominator"]]) %>%
    select_if(is.numeric)

  Num_Miss <- Num[, colSums(Num == 0) > 0, drop = FALSE]
  Denom_Miss <- Denom[, colSums(Denom == 0) > 0, drop = FALSE]
  Metabolites_Miss <- unique(c(colnames(Num_Miss), colnames(Denom_Miss)))

  # Comparisons
  comparisons <- matrix(c(metadata_info[["Numerator"]], metadata_info[["Denominator"]]))

  ## ------------ Perform Hypothesis testing ----------- ##
  for(column in 1:dim(comparisons)[2]){
    C1 <- data %>% # Numerator
      filter(metadata_sample[[metadata_info[["Conditions"]]]] %in% comparisons[1,column]) %>%
      select_if(is.numeric)#only keep numeric columns with metabolite values
    C2 <- data %>% # Denominator
      filter(metadata_sample[[metadata_info[["Conditions"]]]] %in%  comparisons[2,column]) %>%
      select_if(is.numeric)
  }

  # For C1 and C2 we use 0, since otherwise we can not perform the statistical testing.
  C1[is.na(C1)] <- 0
  C2[is.na(C2)] <- 0

  #### 1. p.value and test statistics (=t.val)
  T_C1vC2 <-mapply(pval, x= as.data.frame(C2), y = as.data.frame(C1), SIMPLIFY = F)

  VecPVAL_C1vC2 <- c()
  VecTVAL_C1vC2 <- c()
  for(i in 1:length(T_C1vC2)){
    p_value <- unlist(T_C1vC2[[i]][3])
    t_value <- unlist(T_C1vC2[[i]])[1]   # Extract the t-value
    VecPVAL_C1vC2[i] <- p_value
    VecTVAL_C1vC2[i] <- t_value
  }
  Metabolite <- colnames(C2)
  PVal_C1vC2 <- data.frame(Metabolite, p.val = VecPVAL_C1vC2, t.val = VecTVAL_C1vC2)

  #we set p.val= NA, for metabolites that had 1 or more replicates with NA/0 values and remove them prior to p-value adjustment
  PVal_C1vC2$`NA/0` <- PVal_C1vC2$Metabolite %in% Metabolites_Miss
  PVal_C1vC2 <-PVal_C1vC2%>%
    mutate(p.val = case_when(`NA/0`== TRUE ~ NA,
                             TRUE ~ paste(VecPVAL_C1vC2)))
  PVal_C1vC2$p.val = as.numeric(as.character(PVal_C1vC2$p.val))

  #### 2. p.adjusted
  #Split data for p.value adjustment to exclude NA
  PVal_NA <- PVal_C1vC2[is.na(PVal_C1vC2$p.val), c(1:3)]
  PVal_C1vC2 <-PVal_C1vC2[!is.na(PVal_C1vC2$p.val), c(1:3)]

  #perform adjustment
  VecPADJ_C1vC2 <- stats::p.adjust((PVal_C1vC2[,2]),method = padj, n = length((PVal_C1vC2[,2]))) #p-adjusted
  Metabolite <- PVal_C1vC2[,1]
  PADJ_C1vC2 <- data.frame(Metabolite, p.adj = VecPADJ_C1vC2)
  STAT_C1vC2 <- merge(PVal_C1vC2,PADJ_C1vC2, by="Metabolite")

  #Add Metabolites that have p.val=NA back into the DF for completeness.
  if(nrow(PVal_NA)>0){
    PVal_NA$p.adj <- NA
    STAT_C1vC2 <- rbind(STAT_C1vC2, PVal_NA)
  }

  #Add Log2FC
  if(is.null(log2fc_table)==FALSE){
    STAT_C1vC2 <- merge(log2fc_table,STAT_C1vC2[,c(1:2,4,3)], by="Metabolite")
  }

  #order for t.value
  STAT_C1vC2 <- STAT_C1vC2[order(STAT_C1vC2$t.val,decreasing=TRUE),] # order the df based on the t-value

  #list
  results_list <- list()
  results_list[[paste(metadata_info[["Numerator"]], "_vs_", metadata_info[["Denominator"]])]] <- STAT_C1vC2

  return(invisible(results_list))
}


################################################################
### ### ### aov: Internal Function to perform Anova  ### ### ###
################################################################

#' This helper function to calculate One-vs-All or All-vs-All comparison statistics
#'
#' @param data DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected.
#' @param metadata_sample DF which contains metadata information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames.
#' @param metadata_info  \emph{Optional: } Named vector including the information about the conditions column information on numerator or denominator c(Conditions="ColumnName_SettingsFile", Numerator = "ColumnName_SettingsFile", Denominator  = "ColumnName_SettingsFile"). Denominator and Numerator will specify which comparison(s) will be done (Here all-vs-one, all-vs-all), e.g. Denominator=NULL and Numerator =NULL selects all the condition and performs multiple comparison all-vs-all. \strong{Default = c(conditions="Conditions", numerator = NULL, denumerator = NULL)}
#' @param log2fc_table \emph{Optional: } This is a List of DFs including a column "MetaboliteID" and Log2FC or Log2(Distance). This is the output from MetaProViz:::log2fc. If NULL, the output statistics will not be added into the Log2FC/Log2(Distance) DFs. \strong{Default = NULL}
#'
#' @return List of DFs named after comparison (e.g. tumour versus Normal) with p-value, t-value and adjusted p-value column and column with feature names
#'
#' @keywords Statistical testing, p-value, t-value
#'
#' @importFrom stats aov TukeyHSD
#' @importFrom dplyr rename
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @noRd
aov <- function(data,
               metadata_sample,
               metadata_info=c(Conditions="Conditions", Numerator = NULL, Denominator  = NULL),
               log2fc_table=NULL){

  ## ------------ Create log file ----------- ##
  metaproviz_init()

  ## ------------ Denominator/numerator ----------- ##
  # Denominator and numerator: Define if we compare one_vs_one, one_vs_all or all_vs_all.
  if("Denominator" %in% names(metadata_info)==FALSE  & "Numerator" %in% names(metadata_info) ==FALSE){
    # all-vs-all: Generate all pairwise combinations
    conditions = metadata_sample[[metadata_info[["Conditions"]]]]
    denominator <-unique(metadata_sample[[metadata_info[["Conditions"]]]])
    numerator <-unique(metadata_sample[[metadata_info[["Conditions"]]]])
    comparisons <- combn(unique(conditions), 2) %>% as.matrix()
    #Settings:
    MultipleComparison = TRUE
    all_vs_all = TRUE
  }else if("Denominator" %in% names(metadata_info)==TRUE  & "Numerator" %in% names(metadata_info)==FALSE){
    #all-vs-one: Generate the pairwise combinations
    conditions = metadata_sample[[metadata_info[["Conditions"]]]]
    denominator <- metadata_info[["Denominator"]]
    numerator <-unique(metadata_sample[[metadata_info[["Conditions"]]]])
    # Remove denom from num
    numerator <- numerator[!numerator %in% denominator]
    comparisons  <- t(expand.grid(numerator, denominator)) %>% as.data.frame()
    #Settings:
    MultipleComparison = TRUE
    all_vs_all = FALSE
  }

  #############################################################################################
  ## 1. Anova p.val
  aov.res= apply(data, 2, function(x) stats::aov(x ~ conditions))

  ## 2. Tukey test p.adj
  posthoc.res = lapply(aov.res, TukeyHSD, conf.level=0.95)
  Tukey_res <- do.call('rbind', lapply(posthoc.res, function(x) x[1][[1]][,'p adj'])) %>% as.data.frame()

  comps <-   paste(comparisons[1, ], comparisons[2, ], sep="-")# normal
  opp_comps <-  paste(comparisons[2, ], comparisons[1, ], sep="-")

  if(sum(opp_comps %in%  colnames(Tukey_res))>0){# if opposite comparisons is true
    for (comp in 1: length(opp_comps)){
      colnames(Tukey_res)[colnames(Tukey_res) %in% opp_comps[comp]] <-  comps[comp]
    }
  }

  ## 3. t.val
  Tukey_res_diff <- do.call('rbind', lapply(posthoc.res, function(x) x[1][[1]][,'diff'])) %>% as.data.frame()

  if (sum(opp_comps %in%  colnames(Tukey_res_diff))>0){# if oposite comparisons is true
    for (comp in 1: length(opp_comps)){
      colnames(Tukey_res_diff)[colnames(Tukey_res_diff) %in% opp_comps[comp]] <-  comps[comp]
    }
  }

  #Make output DFs:
  Pval_table <- Tukey_res
  Pval_table <- rownames_to_column(Pval_table,"Metabolite")

  Tval_table <- rownames_to_column(Tukey_res_diff,"Metabolite")

  common_col_names <- setdiff(names(Tukey_res_diff), "row.names")#Here we need to adapt for one_vs_all or all_vs_all

  results_list <- list()
  for(col_name in common_col_names){
    # Create a new data frame by merging the two data frames
    merged_df <- merge(Pval_table[,c("Metabolite",col_name)], Tval_table[,c("Metabolite",col_name)], by="Metabolite", all=TRUE)%>%
      rename("p.adj"=2,
                    "t.val"=3)

    #We need to add _vs_ into the comparison col_name
    pattern <- paste(conditions, collapse = "|")
    conditions_present <- unique(unlist(regmatches(col_name, gregexpr(pattern, col_name))))
    modified_col_name <- paste(conditions_present[1], "vs", conditions_present[2], sep = "_")

    # Add the new data frame to the list with the column name as the list element name
    results_list[[modified_col_name]] <- merged_df
  }

  # Merge the data frames in list1 and list2 based on the "Metabolite" column
  if(is.null(log2fc_table)==FALSE){
    list_names <-  names(results_list)

    merged_list <- list()
    for(name in list_names){
      # Check if the data frames exist in both lists
      if(name %in% names(results_list) && name %in% names(log2fc_table)){
        merged_df <- merge(results_list[[name]], log2fc_table[[name]], by = "Metabolite", all = TRUE)
        merged_df <- merged_df[,c(1,4,2:3,5:ncol(merged_df))]#reorder the columns
        merged_list[[name]] <- merged_df
      }
    }
  }else{
    merged_list <- results_list
  }

  # Make sure the right comparisons are returned:
  if(all_vs_all==TRUE){
    STAT_C1vC2 <- merged_list
  }else if(all_vs_all==FALSE){
    #remove the comparisons that are not needed:
    modified_df_list <- list()
    for(df_name in names(merged_list)){
      if(endsWith(df_name, metadata_info[["Denominator"]])){
        modified_df_list[[df_name]] <- merged_list[[df_name]]
      }
    }
    STAT_C1vC2 <- modified_df_list
  }

  return(invisible(STAT_C1vC2))
}

###########################################################################
### ### ### kruskal: Internal Function to perform kruskal test  ### ### ###
###########################################################################

#' This helper function to calculate One-vs-All or All-vs-All comparison statistics
#'
#' @param data DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected.
#' @param metadata_sample DF which contains metadata information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames.
#' @param metadata_info  \emph{Optional: } Named vector including the information about the conditions column information on numerator or denominator c(Conditions="ColumnName_SettingsFile", Numerator = "ColumnName_SettingsFile", Denominator  = "ColumnName_SettingsFile"). Denominator and Numerator will specify which comparison(s) will be done (Here all-vs-one, all-vs-all), e.g. Denominator=NULL and Numerator =NULL selects all the condition and performs multiple comparison all-vs-all. \strong{Default = c(conditions="Conditions", numerator = NULL, denumerator = NULL)}
#' @param log2fc_table \emph{Optional: } This is a List of DFs including a column "MetaboliteID" and Log2FC or Log2(Distance). This is the output from MetaProViz:::log2fc. If NULL, the output statistics will not be added into the Log2FC/Log2(Distance) DFs. \strong{Default = NULL}
#' @param padj \emph{Optional: } String which contains an abbreviation of the selected p.adjusted test for p.value correction for multiple Hypothesis testing. Search: ?p.adjust for more methods:"BH", "fdr", "bonferroni", "holm", etc.\strong{Default = "fdr"}
#'
#' @return List of DFs named after comparison (e.g. tumour versus Normal) with p-value, t-value and adjusted p-value column and column with feature names
#'
#' @keywords Statistical testing, p-value, t-value
#'
#' @importFrom stats kruskal.test
#' @importFrom rstatix dunn_test
#' @importFrom dplyr rename mutate_all mutate select
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column column_to_rownames
#' @noRd
kruskal <- function(data,
                   metadata_sample,
                   metadata_info=c(Conditions="Conditions", Numerator = NULL, Denominator  = NULL),
                   log2fc_table = NULL,
                   padj = "fdr"
){

  ## ------------ Create log file ----------- ##
  metaproviz_init()

  ## ------------ Denominator/numerator ----------- ##
  # Denominator and numerator: Define if we compare one_vs_one, one_vs_all or all_vs_all.
  if("Denominator" %in% names(metadata_info)==FALSE  & "Numerator" %in% names(metadata_info) ==FALSE){
    # all-vs-all: Generate all pairwise combinations
    conditions = metadata_sample[[metadata_info[["Conditions"]]]]
    denominator <-unique(metadata_sample[[metadata_info[["Conditions"]]]])
    numerator <-unique(metadata_sample[[metadata_info[["Conditions"]]]])
    comparisons <- combn(unique(conditions), 2) %>% as.matrix()
    #Settings:
    MultipleComparison = TRUE
    all_vs_all = TRUE
  }else if("Denominator" %in% names(metadata_info)==TRUE  & "Numerator" %in% names(metadata_info)==FALSE){
    #all-vs-one: Generate the pairwise combinations
    conditions = metadata_sample[[metadata_info[["Conditions"]]]]
    denominator <- metadata_info[["Denominator"]]
    numerator <-unique(metadata_sample[[metadata_info[["Conditions"]]]])
    # Remove denom from num
    numerator <- numerator[!numerator %in% denominator]
    comparisons  <- t(expand.grid(numerator, denominator)) %>% as.data.frame()
    #Settings:
    MultipleComparison = TRUE
    all_vs_all = FALSE
  }

  #############################################################################################
  # kruskal test (p.val)
  aov.res= apply(data,2,function(x) stats::kruskal.test(x~conditions))
  anova_res<-do.call('rbind', lapply(aov.res, function(x) {x["p.value"]}))
  anova_res <- as.matrix(mutate_all(as.data.frame(anova_res), function(x) as.numeric(as.character(x))))
  colnames(anova_res) = c("Kruskal_p.val")

  # Dunn test (p.adj)
  Dunndata <- data %>%
    mutate(conditions = conditions) %>%
    select(conditions, everything())%>%
    as.data.frame()

  # Applying a loop to obtain p.adj and t.val:
  Dunn_Pres<- data.frame(comparisons = paste(comparisons[1,],    comparisons[2,], sep = "_vs_" ))
  Dunn_Tres<- Dunn_Pres
  for(col in 2:dim(Dunndata)[2]){
    data = Dunndata[,c(1,col)]
    colnames(data)[2] <- gsub("^\\d+", "", colnames(data)[2])

    ## If a metabolite starts with number remove it
    formula <- as.formula(paste(colnames(data)[2], "~ conditions"))
    posthoc.res= dunn_test(data, formula, p.adjust.method = padj)

    pres <- data.frame(comparisons = c(paste(posthoc.res$group1, posthoc.res$group2, sep = "_vs_" ), paste(posthoc.res$group2, posthoc.res$group1, sep = "_vs_" )))
    pres[[colnames(Dunndata)[col] ]] <-  c(posthoc.res$p.adj, posthoc.res$p.adj )
    pres <- pres[pres$comparisons %in% Dunn_Pres$comparisons ,] # take only the comparisons selected
    Dunn_Pres <- merge(Dunn_Pres,pres,by="comparisons")

    tres <- data.frame(comparisons = c(paste(posthoc.res$group1, posthoc.res$group2, sep = "_vs_" ), paste(posthoc.res$group2, posthoc.res$group1, sep = "_vs_" )))
    tres[[colnames(Dunndata)[col] ]] <-  c(posthoc.res$statistic, -posthoc.res$statistic )
    tres <- tres[tres$comparisons %in% Dunn_Pres$comparisons ,]# take only the comparisons selected
    Dunn_Tres <- merge(Dunn_Tres,tres,by="comparisons")
  }

  #Make output DFs:
  Dunn_Pres <- column_to_rownames(Dunn_Pres, "comparisons")%>% t() %>% as.data.frame()
  Pval_table <- as.matrix(mutate_all(as.data.frame(Dunn_Pres), function(x) as.numeric(as.character(x)))) %>%
    as.data.frame()%>%
    rownames_to_column("Metabolite")

  Dunn_Tres <- column_to_rownames(Dunn_Tres, "comparisons")%>% t() %>% as.data.frame()
  Tval_table <- as.matrix(mutate_all(as.data.frame(Dunn_Tres), function(x) as.numeric(as.character(x))))%>%
    as.data.frame()%>%
    rownames_to_column("Metabolite")

  common_col_names <- setdiff(names(Dunn_Pres), "row.names")

  results_list <- list()
  for(col_name in common_col_names){
    # Create a new data frame by merging the two data frames
    merged_df <- merge(Pval_table[,c("Metabolite",col_name)], Tval_table[,c("Metabolite",col_name)], by="Metabolite", all=TRUE)%>%
      rename("p.adj"=2,
                    "t.val"=3)
    # Add the new data frame to the list with the column name as the list element name
    results_list[[col_name]] <- merged_df
  }

  # Merge the data frames in list1 and list2 based on the "Metabolite" column
  if(is.null(log2fc_table)==FALSE){
    merged_list <- list()
    for(name in common_col_names){
      # Check if the data frames exist in both lists
      if(name %in% names(results_list) && name %in% names(log2fc_table)){
        merged_df <- merge(results_list[[name]], log2fc_table[[name]], by = "Metabolite", all = TRUE)
        merged_df <- merged_df[,c(1,4,2:3,5:ncol(merged_df))]#reorder the columns
        merged_list[[name]] <- merged_df
      }
    }
    STAT_C1vC2 <- merged_list
  }else{
    STAT_C1vC2 <- results_list
    }

  return(invisible(STAT_C1vC2))
}



#############################################################################################
### ### ### welch: Internal Function to perform anova for unequal variance groups ### ### ###
#############################################################################################

#' Calculate One-vs-All or All-vs-All comparison statistics
#'
#' @param data DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected.
#' @param metadata_sample DF which contains metadata information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames.
#' @param metadata_info  \emph{Optional: } Named vector including the information about the conditions column information on numerator or denominator c(Conditions="ColumnName_SettingsFile", Numerator = "ColumnName_SettingsFile", Denominator  = "ColumnName_SettingsFile"). Denominator and Numerator will specify which comparison(s) will be done (Here all-vs-one, all-vs-all), e.g. Denominator=NULL and Numerator =NULL selects all the condition and performs multiple comparison all-vs-all. \strong{Default = c(conditions="Conditions", numerator = NULL, denumerator = NULL)}
#' @param log2fc_table \emph{Optional: } This is a List of DFs including a column "MetaboliteID" and Log2FC or Log2(Distance). This is the output from MetaProViz:::log2fc. If NULL, the output statistics will not be added into the Log2FC/Log2(Distance) DFs. \strong{Default = NULL}
#'
#' @return List of DFs named after comparison (e.g. tumour versus Normal) with p-value, t-value and adjusted p-value column and column with feature names
#'
#' @keywords Statistical testing, p-value, t-value
#'
#' @importFrom rstatix games_howell_test
#' @importFrom dplyr rename
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @noRd
welch <- function(data,
                 metadata_sample,
                 metadata_info=c(Conditions="Conditions", Numerator = NULL, Denominator  = NULL),
                 log2fc_table=NULL
){
  ## ------------ Create log file ----------- ##
  metaproviz_init()

  ## ------------ Denominator/numerator ----------- ##
  # Denominator and numerator: Define if we compare one_vs_one, one_vs_all or all_vs_all.
  if("Denominator" %in% names(metadata_info)==FALSE  & "Numerator" %in% names(metadata_info) ==FALSE){
    # all-vs-all: Generate all pairwise combinations
    conditions = metadata_sample[[metadata_info[["Conditions"]]]]
    denominator <-unique(metadata_sample[[metadata_info[["Conditions"]]]])
    numerator <-unique(metadata_sample[[metadata_info[["Conditions"]]]])
    comparisons <- combn(unique(conditions), 2) %>% as.matrix()
    #Settings:
    MultipleComparison = TRUE
    all_vs_all = TRUE
  }else if("Denominator" %in% names(metadata_info)==TRUE  & "Numerator" %in% names(metadata_info)==FALSE){
    #all-vs-one: Generate the pairwise combinations
    conditions = metadata_sample[[metadata_info[["Conditions"]]]]
    denominator <- metadata_info[["Denominator"]]
    numerator <-unique(metadata_sample[[metadata_info[["Conditions"]]]])
    # Remove denom from num
    numerator <- numerator[!numerator %in% denominator]
    comparisons  <- t(expand.grid(numerator, denominator)) %>% as.data.frame()
    #Settings:
    MultipleComparison = TRUE
    all_vs_all = FALSE
  }

  ###############################################################################################################
  ## 1. welch's ANOVA using oneway.test is not used by the Games post.hoc function
  #aov.res= apply(Input_data,2,function(x) oneway.test(x~conditions))
  games_data <- merge(metadata_sample, data, by=0)%>%
    rename("conditions"=metadata_info[["Conditions"]])
  games_data$conditions <- conditions
  posthoc.res.list <- list()

  ## 2. Games post hoc test
  for (col in names(data)){ # col = names(Input_data)[1]
    posthoc.res <- games_howell_test(data = games_data,detailed =TRUE, formula = as.formula(paste0(`col`, " ~ ", "conditions"))) %>% as.data.frame()

    result.df <- rbind(data.frame(p.adj = posthoc.res[,"p.adj"],
                                  t.val = posthoc.res[,"statistic"],
                                  row.names = paste(posthoc.res[["group1"]], posthoc.res[["group2"]], sep = "-")),
                       data.frame(p.adj = posthoc.res[,"p.adj"],
                                  t.val = -posthoc.res[,"statistic"],
                                  row.names = paste(posthoc.res[["group2"]], posthoc.res[["group1"]], sep = "-")))
    posthoc.res.list[[col]] <- result.df
  }
  Games_Pres <- do.call('rbind', lapply(posthoc.res.list, function(x) x[,'p.adj'])) %>% as.data.frame()
  colnames(Games_Pres) <- rownames(posthoc.res.list[[1]])
  comps <-   paste(comparisons[1, ], comparisons[2, ], sep="-")# normal
  Games_Pres <- Games_Pres[,colnames(Games_Pres) %in% comps] %>% rownames_to_column("Metabolite")
  # In case of p.adj =0 we change it to 10^-6
  Games_Pres[Games_Pres ==0] <- 0.000001

  ## 3. t.val
  Games_Tres <- do.call('rbind', lapply(posthoc.res.list, function(x) x[,'t.val'])) %>% as.data.frame()
  colnames(Games_Tres) <- rownames(posthoc.res.list[[1]])
  Games_Tres <- Games_Tres[,colnames(Games_Tres) %in% comps] %>% rownames_to_column("Metabolite")

  results_list <- list()
  for(col_name in colnames(Games_Pres)){
    # Create a new data frame by merging the two data frames
    merged_df <- merge(Games_Pres[,c("Metabolite",col_name)], Games_Tres[,c("Metabolite",col_name)], by="Metabolite", all=TRUE)%>%
      rename("p.adj"=2,
                    "t.val"=3)

    #We need to add _vs_ into the comparison col_name
    pattern <- paste(conditions, collapse = "|")
    conditions_present <- unique(unlist(regmatches(col_name, gregexpr(pattern, col_name))))
    modified_col_name <- paste(conditions_present[1], "vs", conditions_present[2], sep = "_")

    # Add the new data frame to the list with the column name as the list element name
    results_list[[modified_col_name]] <- merged_df
  }

  # Merge the data frames in list1 and list2 based on the "Metabolite" column
  if(is.null(log2fc_table)==FALSE){
    list_names <-  names(results_list)

    merged_list <- list()
    for(name in list_names){
      # Check if the data frames exist in both lists
      if(name %in% names(results_list) && name %in% names(log2fc_table)){
        merged_df <- merge(results_list[[name]], log2fc_table[[name]], by = "Metabolite", all = TRUE)
        merged_df <- merged_df[,c(1,4,2:3,5:ncol(merged_df))]#reorder the columns
        merged_list[[name]] <- merged_df
      }
      }
    STAT_C1vC2 <- merged_list
    }else{
      STAT_C1vC2 <-STAT_C1vC2 <- results_list
      }

  return(invisible(STAT_C1vC2))
}


##########################################################################################
### ### ### dma helper function: Internal Function to perform limma ### ### ###
##########################################################################################

#' This helper function to calculate One-vs-One, One-vs-All or All-vs-All comparison statistics
#'
#' @param data DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected.
#' @param metadata_sample DF which contains metadata information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames.
#' @param metadata_info  \emph{Optional: } Named vector including the information about the conditions column information on numerator or denominator c(Conditions="ColumnName_SettingsFile", Numerator = "ColumnName_SettingsFile", Denominator  = "ColumnName_SettingsFile"). Denominator and Numerator will specify which comparison(s) will be done (one-vs-all, all-vs-one, all-vs-all), e.g. Denominator=NULL and Numerator =NULL selects all the condition and performs multiple comparison all-vs-all. \strong{Default = c(conditions="Conditions", numerator = NULL, denumerator = NULL)}
#' @param log2fc_table \emph{Optional: } This is a List of DFs including a column "MetaboliteID" and Log2FC or Log2(Distance). This is the output from MetaProViz:::log2fc. If NULL, the output statistics will not be added into the Log2FC/Log2(Distance) DFs. \strong{Default = NULL}
#' @param padj \emph{Optional: } String which contains an abbreviation of the selected p.adjusted test for p.value correction for multiple Hypothesis testing. Search: ?p.adjust for more methods:"BH", "fdr", "bonferroni", "holm", etc.\strong{Default = "fdr"}
#' @param core \emph{Optional: } TRUE or FALSE for whether a Consumption/Release  input is used. \strong{Default = FALSE}
#' @param transform TRUE or FALSE. If TRUE we expect the data to be not log2 transformed and log2 transformation will be performed within the limma function and Log2FC calculation. If FALSE we expect the data to be log2 transformed as this impacts the Log2FC calculation and limma. \strong{Default= TRUE}
#'
#' @return List of DFs named after comparison (e.g. tumour versus Normal) with p-value, t-value and adjusted p-value column and column with feature names
#'
#' @keywords Statistical testing, p-value, t-value
#'
#' @importFrom limma lmFit makeContrasts contrasts.fit eBayes topTable
#' @importFrom dplyr rename arrange filter distinct
#' @importFrom tidyr separate unite
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column column_to_rownames
#' @noRd
dma_stat_limma <- function(data,
                           metadata_sample,
                           metadata_info=c(Conditions="Conditions", Numerator = NULL, Denominator  = NULL),
                           log2fc_table = NULL,
                           padj ="fdr",
                           core=FALSE,
                           transform= TRUE){

  ## ------------ Create log file ----------- ##
  metaproviz_init()

  ## ------------ Denominator/numerator ----------- ##
  # Denominator and numerator: Define if we compare one_vs_one, one_vs_all or all_vs_all.
  if("Denominator" %in% names(metadata_info)==FALSE  & "Numerator" %in% names(metadata_info) ==FALSE){
    MultipleComparison = TRUE
    all_vs_all = TRUE
  }else if("Denominator" %in% names(metadata_info)==TRUE  & "Numerator" %in% names(metadata_info)==FALSE){
    MultipleComparison = TRUE
    all_vs_all = FALSE
  }else if("Denominator" %in% names(metadata_info)==TRUE  & "Numerator" %in% names(metadata_info)==TRUE){
    MultipleComparison = FALSE
    all_vs_all = FALSE
  }

  ####------ Ensure that Input_data is ordered by conditions and sample names are the same as in Input_metadata_sample:
  targets <- metadata_sample%>%
    rownames_to_column("sample")
  targets<- targets[,c("sample", metadata_info[["Conditions"]])]%>%
    rename("condition"=2)%>%
    arrange(sample)#Order the column "sample" alphabetically
  targets$condition_limma_compatible <-make.names(targets$condition)#make appropriate condition names accepted by limma

  if(MultipleComparison==FALSE){
    #subset the data:
    targets<-targets%>%
      subset(condition==metadata_info[["Numerator"]] | condition==metadata_info[["Denominator"]])%>%
      arrange(sample)#Order the column "sample" alphabetically

    Limma_input <- data%>%tibble::rownames_to_column("sample")
    Limma_input <-merge(targets[,1:2],  Limma_input, by="sample", all.x=TRUE)
    Limma_input <- Limma_input[,-2]%>%
      arrange(sample)#Order the column "sample" alphabetically
  }else if(MultipleComparison==TRUE){
    Limma_input <- data%>%tibble::rownames_to_column("sample")%>%
      arrange(sample)#Order the column "sample" alphabetically
  }

  #Check if the order of the "sample" column is the same in both data frames
  if(identical(targets$sample, Limma_input$sample)==FALSE){
    stop("The order of the 'sample' column is different in both data frames. Please make sure that Input_metadata_sample and Input_data contain the same rownames and sample numbers.")
  }

  targets_limma <-targets[,-2]%>%
    rename("condition"="condition_limma_compatible")

  #We need to transpose the df to run limma. Also, if the data is not log2 transformed, we will not calculate the Log2FC as limma just substracts one condition from the other
  Limma_input <- as.data.frame(t(Limma_input%>%tibble::column_to_rownames("sample")))

  if(transform==TRUE){
    Limma_input <- log2(Limma_input) # communicate the log2 transformation --> how does limma deals with NA when calculating the change?
  }

  #### ------Run limma:
  ####  Make design matrix:
  fcond <- as.factor(targets_limma$condition)#all versus all

  design <- model.matrix(~0 + fcond)# Create the design matrix
  colnames(design) <- levels(fcond) # Give meaningful column names to the design matrix

  #### Fit the linear model
  fit <- limma::lmFit(Limma_input, design)

  ####  Make contrast matrix:
  if(all_vs_all ==TRUE & MultipleComparison==TRUE){
    unique_conditions <- levels(fcond)# Get unique conditions

    # Create an empty contrast matrix
    num_conditions <- length(unique_conditions)
    num_comparisons <- num_conditions * (num_conditions - 1) / 2
    cont.matrix <- matrix(0, nrow = num_comparisons, ncol = num_conditions)

    # Initialize an index for the column in the contrast matrix
    i <- 1

    # Initialize column and row names
    colnames(cont.matrix) <- unique_conditions
    rownames(cont.matrix) <- character(num_comparisons)

    # Loop through all pairwise combinations of unique conditions
    for (condition1 in 1:(num_conditions - 1)) {
      for (condition2 in (condition1 + 1):num_conditions) {
        # Create the pairwise comparison vector
        comparison <- rep(0, num_conditions)

        comparison[condition2] <- -1
        comparison[condition1] <- 1
        # Add the comparison vector to the contrast matrix
        cont.matrix[i, ] <- comparison
        # Set row name
        rownames(cont.matrix)[i] <- paste(unique_conditions[condition1], "_vs_", unique_conditions[condition2], sep="")
        i <- i + 1
      }
    }
    cont.matrix<- t(cont.matrix)
  }else if(all_vs_all ==FALSE & MultipleComparison==TRUE){
    unique_conditions <- levels(fcond)# Get unique conditions
    denominator  <- make.names(metadata_info[["Denominator"]])

    # Create an empty contrast matrix
    num_conditions <- length(unique_conditions)
    num_comparisons <- num_conditions - 1
    cont.matrix <- matrix(0, nrow = num_comparisons, ncol = num_conditions)


    # Initialize an index for the column in the contrast matrix
    i <- 1

    # Initialize column and row names
    colnames(cont.matrix) <- unique_conditions
    rownames(cont.matrix) <- character(num_comparisons)

    # Loop through all pairwise combinations of unique conditions
    for(condition in 2:num_conditions){
      # Create the pairwise comparison vector
      comparison <- rep(0, num_conditions)
      if(unique_conditions[1]== make.names(metadata_info[["Denominator"]])){
        comparison[1] <- -1
        comparison[condition] <- 1
        # Add the comparison vector to the contrast matrix
        cont.matrix[i, ] <- comparison
        # Set row name
        rownames(cont.matrix)[i] <- paste(unique_conditions[condition], "_vs_", unique_conditions[1], sep = "")
      }else{
        comparison[1] <- 1
        comparison[condition] <- -1
        # Add the comparison vector to the contrast matrix
        cont.matrix[i, ] <- comparison
        # Set row name
        rownames(cont.matrix)[i] <- paste(unique_conditions[1], "_vs_", unique_conditions[condition], sep = "")

      }
      i <- i + 1
    }
    cont.matrix<- t(cont.matrix)
  }else if(all_vs_all ==FALSE & MultipleComparison==FALSE){
    Name_Comp <- paste(make.names(metadata_info[["Numerator"]]), "-", make.names(metadata_info[["Denominator"]]), sep="")
    cont.matrix <- as.data.frame(makeContrasts(contrasts=Name_Comp, levels=colnames(design)))%>%
      rename(!!paste(make.names(metadata_info[["Numerator"]]), "_vs_", make.names(metadata_info[["Denominator"]]), sep="") := 1)
    cont.matrix <-as.matrix(cont.matrix)
  }

  # Fit the linear model with contrasts
  #fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)# Perform empirical Bayes moderation

  #### ------Extract results:
  contrast_names <- colnames(fit2$coefficients)  # Get all contrast names

  results_list <- list()# Create an empty list to store results data frames
  for (contrast_name in contrast_names) {
    # Extract results for the current contrast
    res.t <- topTable(fit2, coef=contrast_name, n=Inf, sort.by="n", adjust.method = padj)%>% # coef= the comparison the test is done for!
      rename("Log2FC"=1,
                    "t.val"=3,
                    "p.val"=4,
                    "p.adj"=5)

    res.t <- res.t%>%
      rownames_to_column("Metabolite")

    # Store the data frame in the results list, named after the contrast
    results_list[[contrast_name]] <- res.t
  }

  #Make the name_match_df
  name_match_df <- as.data.frame(names(results_list))%>%
    separate("names(results_list)", into=c("a", "b"), sep="_vs_", remove=FALSE)

  name_match_df <-merge(name_match_df, targets[,-c(1)] , by.x="a", by.y="condition_limma_compatible", all.x=TRUE)%>%
    rename("Condition1"=4)
  name_match_df <- merge(name_match_df, targets[,-c(1)] , by.x="b", by.y="condition_limma_compatible", all.x=TRUE)%>%
    rename("Condition2"=5)%>%
    unite("New", "Condition1", "Condition2", sep="_vs_", remove=FALSE)

  name_match_df<- name_match_df[,c(3,4)]%>%
    distinct(New, .keep_all = TRUE)

  results_list_new <- list()
  #Match the lists using name_match_df
  for(i in 1:nrow(name_match_df)){
    old_name <- name_match_df$`names(results_list)`[i]
    new_name <- name_match_df$New[i]
    results_list_new[[new_name]] <- results_list[[old_name]]
    }

  if(is.null(log2fc_table)==FALSE){
    if(core==TRUE){#If core=TRUE, we need to exchange the Log2FC with the Distance and we need to combine the lists
      # Merge the data frames in list1 and list2 based on the "Metabolite" column
      merged_list <- list()
      for(i in 1:nrow(name_match_df)){
        list_dfs <- name_match_df$New[i]

        # Check if the data frames exist in both lists
        if(list_dfs %in% names(results_list_new) && list_dfs %in% names(log2fc_table)){
          merged_df <- merge(results_list_new[[list_dfs]], log2fc_table[[list_dfs]], by = "Metabolite", all = TRUE)
          merged_list[[list_dfs]] <- merged_df
        }
      }
    STAT_C1vC2 <- merged_list
  }else{
    STAT_C1vC2 <- results_list_new
  }
  }

  #Add input data
  Cond <- metadata_sample%>%
    rownames_to_column("Code")

  InputReturn <- merge(Cond[,c("Code",metadata_info[["Conditions"]])], as.data.frame(t(Limma_input)),by.x="Code", by.y=0, all.y=TRUE)

  for(DFs in names(STAT_C1vC2)){
    parts <- unlist(strsplit(DFs, "_vs_"))
    C1 <- parts[1]
    C2 <- parts[2]
    InputReturn_Filt <- InputReturn%>%
      filter(get(metadata_info[["Conditions"]])==C1 | get(metadata_info[["Conditions"]])==C2)%>%
      column_to_rownames("Code")
    InputReturn_Filt <-as.data.frame(t(InputReturn_Filt[,-c(1)]))

    if(transform==TRUE){#Add prefix & suffix to each column since the data have been log2 transformed!
      colnames(InputReturn_Filt) <- paste0("log2(", colnames(InputReturn_Filt), ")")
      }

    InputReturn_Merge <- merge(STAT_C1vC2[[DFs]], InputReturn_Filt, by.x="Metabolite", by.y=0, all.x=TRUE)

    STAT_C1vC2[[DFs]] <- InputReturn_Merge
  }

  #return
  return(invisible(STAT_C1vC2))
}


#############################################################################################
### ### ### shapiro function: Internal Function to perform shapiro test and plots ### ### ###
#############################################################################################

#' Shapiro test and plots
#'
#' @param data DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected.
#' @param metadata_sample DF which contains metadata information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames.
#' @param metadata_info  \emph{Optional: } Named vector including the information about the conditions column information on numerator or denominator c(Conditions="ColumnName_SettingsFile", Numerator = "ColumnName_SettingsFile", Denominator  = "ColumnName_SettingsFile"). Denominator and Numerator will specify which comparison(s) will be done (one-vs-all, all-vs-one, all-vs-all), e.g. Denominator=NULL and Numerator =NULL selects all the condition and performs multiple comparison all-vs-all. \strong{Default = c(conditions="Conditions", numerator = NULL, denumerator = NULL)}
#' @param pval \emph{Optional: } String which contains an abbreviation of the selected test to calculate p.value. For one-vs-one comparisons choose t.test, wilcox.test, "chisq.test" or "cor.test", for one-vs-all or all-vs-all comparison choose aov (=annova), kruskal.test or lmFit (=limma) \strong{Default = "t-test"}
#' @param qqplots \emph {Optional: } TRUE or FALSE for whether QQ plots should be plotted  \strong{default = TRUE}
#'
#' @return List with tewo entries: DF (including the results DF) and Plots (including the Density and QQ plots)
#'
#' @keywords shapiro test,Normality testing, Density plot, QQplot
#'
#' @importFrom stats shapiro.test
#' @importFrom ggplot2 ggplot geom_histogram geom_density scale_x_continuous
#' @importFrom ggplot2 theme_minimal labs ggplot_build geom_qq geom_qq_line after_stat
#' @importFrom dplyr rename select_if filter
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column column_to_rownames
#' @noRd
shapiro <- function(data,
                   metadata_sample,
                   metadata_info=c(Conditions="Conditions", Numerator = NULL, Denominator  = NULL),
                   pval= "t-test",
                   qqplots=TRUE
){

  ## ------------ Create log file ----------- ##
  metaproviz_init()

  ## ------------- Checks --------------##
  if(grepl("[[:space:]()-./\\\\]", metadata_info[["Conditions"]])==TRUE){
    message("In metadata_info=c(Conditions= ColumnName): ColumnName contains special charaters, hence this is renamed.")
    ColumnNameCondition_clean <- gsub("[[:space:]()-./\\\\]", "_", metadata_info[["Conditions"]])
    metadata_sample <- metadata_sample%>%
      rename(!!paste(ColumnNameCondition_clean):= metadata_info[["Conditions"]])

    metadata_info[["Conditions"]] <- ColumnNameCondition_clean
  }

  ## ------------ Denominator/numerator ----------- ##
  # Denominator and numerator: Define if we compare one_vs_one, one_vs_all or all_vs_all.
  if("Denominator" %in% names(metadata_info)==FALSE  & "Numerator" %in% names(metadata_info) ==FALSE){
    # all-vs-all: Generate all pairwise combinations
    conditions = metadata_sample[[metadata_info[["Conditions"]]]]
    denominator <-unique(metadata_sample[[metadata_info[["Conditions"]]]])
    numerator <-unique(metadata_sample[[metadata_info[["Conditions"]]]])
    comparisons <- combn(unique(conditions), 2) %>% as.matrix()
    #Settings:
    MultipleComparison = TRUE
    all_vs_all = TRUE
   }else if("Denominator" %in% names(metadata_info)==TRUE  & "Numerator" %in% names(metadata_info)==FALSE){
    #all-vs-one: Generate the pairwise combinations
    conditions = metadata_sample[[metadata_info[["Conditions"]]]]
    denominator <- metadata_info[["Denominator"]]
    numerator <-unique(metadata_sample[[metadata_info[["Conditions"]]]])
    # Remove denom from num
    numerator <- numerator[!numerator %in% denominator]
    comparisons  <- t(expand.grid(numerator, denominator)) %>% as.data.frame()
    #Settings:
    MultipleComparison = TRUE
    all_vs_all = FALSE
  }else if("Denominator" %in% names(metadata_info)==TRUE  & "Numerator" %in% names(metadata_info)==TRUE){
    # one-vs-one: Generate the comparisons
    denominator <- metadata_info[["Denominator"]]
    numerator <- metadata_info[["Numerator"]]
    comparisons <- matrix(c(metadata_info[["Denominator"]], metadata_info[["Numerator"]]))
    #Settings:
    MultipleComparison = FALSE
    all_vs_all = FALSE
  }

  ################################################################################################################################################################################################
  ## ------------ Check data normality and statistical test chosen and generate Output DF----------- ##
  # Before Hypothesis testing, we have to decide whether to use a parametric or a non parametric test. We can test the data normality using the shapiro test.
  ##-------- First: Load the data and perform the shapiro.test on each metabolite across the samples of one condition. this needs to be repeated for each condition:
  #Prepare the input:
  Input_shaptest <- replace(data, is.na(data), 0)%>% #shapiro test can not handle NAs!
    filter(metadata_sample[[metadata_info[["Conditions"]]]] %in% numerator | metadata_sample[[metadata_info[["Conditions"]]]] %in% denominator)%>%
    select_if(is.numeric)
  temp<- sapply(Input_shaptest, function(x, na.rm = TRUE) var(x)) == 0#  we have to remove features with zero variance if there are any.
  temp <- temp[complete.cases(temp)]  # Remove NAs from temp
  columns_with_zero_variance <- names(temp[temp])# Extract column names where temp is TRUE

  if(length(Input_shaptest)==1){#handle a specific case where after filtering and selecting numeric variables, there's only one column left in Input_shaptest
    Input_shaptest <-data
  }else{
    if(length(columns_with_zero_variance)==0){
      Input_shaptest <-Input_shaptest
    }else{
      message("The following features have zero variance and are removed prior to performing the shaprio test: ",columns_with_zero_variance)
      Input_shaptest <- Input_shaptest[,!(names(Input_shaptest) %in% columns_with_zero_variance), drop = FALSE]#drop = FALSE argument is used to ensure that the subset operation doesn't simplify the result to a vector, preserving the data frame structure
    }
  }

  Input_shaptest_Cond <-merge(data.frame(Conditions = metadata_sample[, metadata_info[["Conditions"]], drop = FALSE]), Input_shaptest, by=0, all.y=TRUE)

  UniqueConditions <- metadata_sample%>%
    subset(metadata_sample[[metadata_info[["Conditions"]]]] %in% numerator | metadata_sample[[metadata_info[["Conditions"]]]] %in% denominator, select = c(metadata_info[["Conditions"]]))
  UniqueConditions <- unique(UniqueConditions[[metadata_info[["Conditions"]]]])

  #Generate the results
  shapiro_results <- list()
  for (i in UniqueConditions) {
    # Subset the data for the current condition
    subset_data <- Input_shaptest_Cond%>%
      column_to_rownames("Row.names")%>%
      subset(get(metadata_info[["Conditions"]]) == i, select = -c(1))

    #Check the sample size (shapiro.test(x) : sample size must be between 3 and 5000):
    if(nrow(subset_data)<3){
      warning("shapiro.test(x) : sample size must be between 3 and 5000. You have provided <3 Samples for condition ", i, ". Hence Shaprio test can not be performed for this condition.", sep="")
    }else if(nrow(subset_data)>5000){
      warning("shapiro.test(x) : sample size must be between 3 and 5000. You have provided >5000 Samples for condition ", i, ". Hence Shaprio test will not be performed for this condition.", sep="")
    }else{
      # Apply shapiro-Wilk test to each feature in the subset
     shapiro_results[[i]] <- as.data.frame(sapply(subset_data, function(x) shapiro.test(x)))
    }
  }

  if(nrow(subset_data)>=3 & nrow(subset_data)<=5000){
    #Make the output DF
    DF_shapiro_results <- as.data.frame(matrix(NA, nrow = length(UniqueConditions), ncol = ncol(Input_shaptest)))
    rownames(DF_shapiro_results) <- UniqueConditions
    colnames(DF_shapiro_results) <- colnames(Input_shaptest)
    for(k in 1:length(UniqueConditions)){
      for(l in 1:ncol(Input_shaptest)){
        DF_shapiro_results[k, l] <- shapiro_results[[UniqueConditions[k]]][[l]]$p.value
      }
    }
    colnames(DF_shapiro_results) <- paste("shapiro p.val(", colnames(DF_shapiro_results),")", sep = "")

    ##------ Second: Give feedback to the user if the chosen test fits the data distribution. The data are normal if the p-value of the shapiro.test > 0.05.
    Density_plots <- list()
    if(qqplots==TRUE){
      QQ_plots <- list()
    }
    for(x in 1:nrow(DF_shapiro_results)){
      #### Generate Results Table
      transpose <- as.data.frame(t(DF_shapiro_results[x,]))
      Norm <- format((round(sum(transpose[[1]] > 0.05)/nrow(transpose),4))*100, nsmall = 2) # percentage of normally distributed metabolites across samples
      NotNorm <- format((round(sum(transpose[[1]] < 0.05)/nrow(transpose),4))*100, nsmall = 2) # percentage of not-normally distributed metabolites across samples
      if(pval =="kruskal.test" | pval =="wilcox.test"){
        message("For the condition ", colnames(transpose) ," ", Norm, " % of the metabolites follow a normal distribution and ", NotNorm, " % of the metabolites are not-normally distributed according to the shapiro test. You have chosen ",paste(pval), ", which is for non parametric Hypothesis testing. `shapiro.test` ignores missing values in the calculation.")
      }else{
        message("For the condition ", colnames(transpose) ," ", Norm, " % of the metabolites follow a normal distribution and ", NotNorm, " % of the metabolites are not-normally distributed according to the shapiro test. You have chosen ",paste(pval), ", which is for parametric Hypothesis testing. `shapiro.test` ignores missing values in the calculation.")
      }

      # Assign the calculated values to the corresponding rows in result_df
      DF_shapiro_results$`Metabolites with normal distribution [%]`[x] <- Norm
      DF_shapiro_results$`Metabolites with not-normal distribution [%]`[x] <- NotNorm

      #Reorder the DF:
      all_columns <- colnames(DF_shapiro_results)
      include_columns <- c("Metabolites with normal distribution [%]", "Metabolites with not-normal distribution [%]")
      exclude_columns <- setdiff(all_columns, include_columns)
      DF_shapiro_results <- DF_shapiro_results[, c(include_columns, exclude_columns)]

      #### Make Group wise data distribution plot and QQ plots
      subset_data <- Input_shaptest_Cond%>%
        column_to_rownames("Row.names")%>%
        subset(get(metadata_info[["Conditions"]]) ==  colnames(transpose), select = -c(1))
      all_data <- unlist(subset_data)

      plot <- ggplot(data.frame(x = all_data), aes(x = x)) +
        geom_histogram(aes(y=after_stat(density)), binwidth=.5, colour="black", fill="white")  +
        geom_density(alpha = 0.2, fill = "grey45")

      density_values <- ggplot_build(plot)$data[[2]]

      plot <- ggplot(data.frame(x = all_data), aes(x = x)) +
        geom_histogram( aes(y=after_stat(density)), binwidth=.5, colour="black", fill="white") +
        geom_density(alpha=.2, fill="grey45") +
        scale_x_continuous(limits = c(0, density_values$x[max(which(density_values$scaled >= 0.1))]))

      density_values2 <- ggplot_build(plot)$data[[2]]

      suppressWarnings(sampleDist <- ggplot(data.frame(x = all_data), aes(x = x)) +
                         geom_histogram(aes(y=after_stat(density)), binwidth=.5, colour="black", fill="white") +
                         geom_density(alpha=.2, fill="grey45") +
                         scale_x_continuous(limits = c(0, density_values$x[max(which(density_values$scaled >= 0.1))])) +
                         theme_minimal()+
                         labs(title=paste("data distribution ",  colnames(transpose)), subtitle = paste(NotNorm, " of metabolites not normally distributed based on shapiro test"),x="Abundance", y = "Density")
      )

      Density_plots[[paste(colnames(transpose))]] <- sampleDist

      # QQ plots
      if(qqplots==TRUE){
        # Make folders
        conds <- unique(c(numerator, denominator))

        #QQ plots for each groups for each metabolite for normality visual check
        qq_plot_list <- list()
        for (col_name in colnames(subset_data)){
          qq_plot <- ggplot(data.frame(x = subset_data[[col_name]]), aes(sample = x)) +
            geom_qq() +
            geom_qq_line(color = "red") +
            labs(title = paste("QQPlot for", col_name),x = "Theoretical", y="Sample")+ theme_minimal()

          plot.new()
          plot(qq_plot)
          qq_plot_list[[col_name]] <-  recordPlot()

          col_name2 <- (gsub("/","_",col_name))#remove "/" cause this can not be safed in a PDF name
          col_name2 <- gsub("-", "", col_name2)
          col_name2 <- gsub("/", "", col_name2)
          col_name2 <- gsub(" ", "", col_name2)
          col_name2 <- gsub("\\*", "", col_name2)
          col_name2 <- gsub("\\+", "", col_name2)
          col_name2 <- gsub(",", "", col_name2)
          col_name2 <- gsub("\\(", "", col_name2)
          col_name2 <- gsub("\\)", "", col_name2)

          dev.off()
        }

        QQ_plots[[paste(colnames(transpose))]] <- qq_plot_list
      }
    }

    ######################################
    ##-------- Return
    #Here we make a list
    if(qqplots==TRUE){
      Shapiro_output_list <- list("DF" = list("Shapiro_result"=DF_shapiro_results%>%tibble::rownames_to_column("Code")),"Plot"=list( "Distributions"=Density_plots, "QQ_plots" = QQ_plots))
    }else{
      Shapiro_output_list <- list("DF" = list("Shapiro_result"=DF_shapiro_results%>%tibble::rownames_to_column("Code")),"Plot"=list( "Distributions"=Density_plots))
    }

    suppressWarnings(invisible(return(Shapiro_output_list)))
    }
}



###########################################################################################
### ### ### bartlett function: Internal Function to perform bartlett test and plots ### ###
###########################################################################################

#' Bartlett test for variance homogeneity check across groups
#'
#' @param data DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected.
#' @param metadata_sample DF which contains metadata information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames.
#' @param metadata_info  \emph{Optional: } Named vector including the information about the conditions column information on numerator or denominator c(Conditions="ColumnName_SettingsFile", Numerator = "ColumnName_SettingsFile", Denominator  = "ColumnName_SettingsFile"). Denominator and Numerator will specify which comparison(s) will be done (one-vs-all, all-vs-one, all-vs-all), e.g. Denominator=NULL and Numerator =NULL selects all the condition and performs multiple comparison all-vs-all. \strong{Default = c(conditions="Conditions", numerator = NULL, denumerator = NULL)}
#'
#' @return List with two entries: DF (including the results DF) and Plots (including the  histogramm plot)
#'
#' @keywords bartlett test,Normality testing, Density plot, QQplot
#'
#' @importFrom stats bartlett.test
#' @importFrom ggplot2 ggplot geom_histogram geom_density ggtitle xlab geom_vline
#' @importFrom dplyr mutate
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @noRd
bartlett <- function(data,
                    metadata_sample,
                    metadata_info){

  ## ------------ Create log file ----------- ##
  metaproviz_init()

  ################################################################################################################################################################################################

  conditions = metadata_sample[[metadata_info[["Conditions"]]]]

  # Use Bartletts test
  bartlett_res =  apply(data,2,function(x) bartlett.test(x~conditions))

  #Make the output DF
  DF_bartlett_results <- as.data.frame(matrix(NA, nrow = ncol(data)), ncol = 1)
  rownames(DF_bartlett_results) <- colnames(data)
  colnames(DF_bartlett_results) <- "bartlett p.val"

  for(l in 1:length(bartlett_res)){
    DF_bartlett_results[l, 1] <-bartlett_res[[l]]$p.value
  }
  DF_bartlett_results <- DF_bartlett_results %>% mutate(`Var homogeneity`= case_when(`bartlett p.val`< 0.05~ FALSE,
                                                                                     `bartlett p.val`>=0.05 ~ TRUE))
  # if p<0.05 then unequal variances
  message("For ",round(sum(DF_bartlett_results$`Var homogeneity`)/  nrow(DF_bartlett_results), digits = 4) * 100, "% of metabolites the group variances are equal.")

  DF_bartlett_results <- DF_bartlett_results %>% rownames_to_column("Metabolite") %>% relocate("Metabolite")
  DF_Bartlett_results_out <- DF_bartlett_results

  #### Plots:
  #Make density plots
  Bartlettplot <- ggplot(data.frame(x = DF_Bartlett_results_out), aes(x =DF_Bartlett_results_out$`bartlett p.val`)) +
    geom_histogram(aes(y=..density..), colour="black", fill="white")  +
    geom_density(alpha = 0.2, fill = "grey45")+
    ggtitle("bartlett's test p.value distribution") +
    xlab("p.value")+
    geom_vline(aes(xintercept = 0.05, color="darkred"))

  Bartlett_output_list<- list("DF"=list("Bartlett_result"=DF_Bartlett_results_out) , "Plot"=list("Histogram"=Bartlettplot))

  suppressWarnings(invisible(return(Bartlett_output_list)))

}


################################################################
### ### ### Variance stabilizing transformation function ### ###
################################################################

#' Variance stabilizing transformation (vst)
#'
#' @param data data frame with unique sample identifiers
# not true: no need to have row names, they are not used here, and in general,
# it is not a good practice to use row names
#'     as row names and metabolite numerical values in columns with metabolite
#'     identifiers as column names. Use NA for metabolites that were not detected.
#'
#' @return List with two entries: DF (including the results DF) and Plots (including the scedasticity_plot)
#'
#' @keywords Heteroscedasticity, variance stabilizing transformation
#'
#' @importFrom tidyr pivot_longer
#' @importFrom stats lm
#' @importFrom ggplot2 ggplot geom_point theme_bw scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous xlab ylab geom_abline ggtitle
#' @importFrom ggplot2 geom_smooth aes
#' @importFrom patchwork wrap_plots
#' @importFrom dplyr summarise group_by
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @noRd
vst <- function(data){

  ## ------------ Create log file ----------- ##
  metaproviz_init()

  # model the mean and variance relationship on the data
  het.data <-
    data %>%
    pivot_longer(-1L, names_to = 'variable') %>%
    group_by(variable) %>% # make a dataframe to save the values
    summarise(mean=mean(value), sd=sd(value))
  het.data$lm <- 1 # add a common group for the lm function to account for the whole data together

  invisible(het_plot <-
              ggplot(het.data, aes(x = mean, y = sd)) +
              geom_point() +
              theme_bw() +
              scale_x_continuous(trans='log2') +
              scale_y_continuous(trans='log2') +
              xlab("log(mean)") +
              ylab("log(sd)") +
              geom_abline(intercept = 0, slope = 1)  +
              ggtitle(" data heteroscedasticity")  +
              geom_smooth(aes(group=lm),method='lm', formula= y~x, color = "red"))

  # select data
  prevst.data <- het.data
  prevst.data$mean <- log(prevst.data$mean)
  prevst.data$sd <- log(prevst.data$sd)

  # calculate the slope of the log data
  data.fit <- lm(sd~mean, prevst.data)
  coef(data.fit)

  # Make the vst transformation
  data.vst <- as.data.frame(data^(1-coef(data.fit)['mean'][1]))

  # Heteroscedasticity visual check again
  het.vst.data <-
    data.vst %>%
    pivot_longer(-1L, names_to = 'variable') %>%
    group_by(variable) %>% # make a dataframe to save the values
    summarise(mean=mean(value), sd=sd(value))
  het.vst.data$lm <- 1 # add a common group for the lm function to account for the whole data together

  # plot variable stadard deviation as a function of the mean
  invisible(hom_plot <-
              ggplot(het.vst.data,  aes(x = mean, y = sd)) +
              geom_point() +
              theme_bw() +
              scale_x_continuous(trans='log2') +
              scale_y_continuous(trans='log2') +
              xlab("log(mean)") +
              ylab("log(sd)") +
              geom_abline(intercept = 0)  +
              ggtitle("vst transformed data")  +
              geom_smooth(aes(group=lm),method='lm', formula= y~x, color = "red"))

  invisible(scedasticity_plot <- wrap_plots(het_plot,hom_plot))

  return(invisible(list("DFs" = list("Corrected_data" = data.vst), "Plots" = list("scedasticity_plot" = scedasticity_plot))))
}

