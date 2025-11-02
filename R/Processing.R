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
# Metabolomics pre-processing
#

#' Data preprocessing and normalization
#'
#' Applies modular normalization including 80% filtering rule, total-ion count
#' normalization, missing value imputation, and outlier detection using
#' Hotelling's T2 test.
#'
#' @param data SummarizedExperiment or data frame. If SummarizedExperiment,
#'     metadata_sample is extracted from colData. If data frame, provide unique
#'     sample identifiers as row names and metabolite numerical values in
#'     columns with metabolite identifiers as column names. Use NA for
#'     undetected metabolites.
#' @param metadata_sample Data frame (optional). Only required if data is not a
#'     SummarizedExperiment. Contains information about samples, combined with
#'     input data based on unique sample identifiers used as row names. Must
#'     contain Conditions column. If experiment has no multiple conditions,
#'     assign all samples to same condition. Default: NULL.
#' @param metadata_info Named character vector (optional). Contains names of experimental
#'     parameters: c(Conditions="ColumnName",
#'     Biological_Replicates="ColumnName"). Column "Conditions" (mandatory)
#'     contains sample conditions (e.g., "N"/"T" or "Normal"/"Tumor"), used for
#'     feature filtering and PCA color coding. Column "BiologicalReplicates"
#'     (optional) contains numerical values. For core=TRUE, must also add
#'     core_norm_factor="ColumnName" and core_media="ColumnName". Column
#'     core_norm_factor is used for normalization; core_media specifies media
#'     controls in Conditions. Default: c(Conditions="Conditions").
#' @param featurefilt Character (optional). If NULL, no feature filtering is performed. If
#'     "Standard", applies 80% filtering rule (Bijlsma et al., 2006) on
#'     metabolite features across whole dataset. If "Modified", filtering is
#'     done per condition and Conditions column must be provided (Yang et al.,
#'     2015). Default: "Standard".
#' @param cutoff_featurefilt Numeric (optional). Percentage threshold for feature filtering. Default:
#'     0.8.
#' @param tic Logical (optional). Whether total ion count normalization is performed.
#'     Default: TRUE.
#' @param mvi Logical (optional). Whether missing value imputation (MVI) based on half
#'     minimum is performed. Default: TRUE.
#' @param mvi_percentage Numeric (optional). Percentage (0-100) of imputed value based on minimum
#'     value. Default: 50.
#' @param hotellins_confidence Numeric (optional). Confidence level for outlier identification in
#'     Hotelling's T2 test. Default: 0.99.
#' @param core Logical (optional). Whether consumption-release experiment was performed
#'     and core value should be calculated. If TRUE, provide normalization
#'     factor column "core_norm_factor" in metadata_sample where Conditions
#'     column matches. The normalization factor must be numerical value from
#'     growth rate (growth curve) or growth factor (ratio of cell count/protein
#'     quantification at start vs. end point). Additionally, control media
#'     samples must be available in data and defined as "core_media" in
#'     Conditions column of metadata_sample. Default: FALSE.
#' @param save_plot Character (optional). File type of output plots: "svg", "png", "pdf". If
#'     NULL, plots are not saved. Default: "svg".
#' @param save_table Character (optional). File type of output table: "csv", "xlsx", "txt".
#'     If NULL, tables are not saved. Default: "csv".
#' @param print_plot Logical (optional). Whether to print overview of resulting plots.
#'     Default: TRUE.
#' @param path Character (optional). Path to folder where results should be saved.
#'     Default: NULL.
#'
#' @return List with two elements: DF (all output tables generated) and Plot (all
#'     plots generated).
#'
#' @examples
#' data(intracell_raw_se)
#' ResI <- processing(
#'     data = intracell_raw_se,
#'     metadata_info = c(
#'         Conditions = "Conditions",
#'         Biological_Replicates = "Biological_Replicates"
#'     )
#' )
#'
#' data(intracell_raw)
#' Intra <- intracell_raw %>% tibble::column_to_rownames("Code")
#' ResI <- processing(
#'     data = Intra[-c(49:58), -c(1:3)],
#'     metadata_sample = Intra[-c(49:58), c(1:3)],
#'     metadata_info = c(
#'         Conditions = "Conditions",
#'         Biological_Replicates = "Biological_Replicates"
#'     )
#' )
#'
#' data(medium_raw)
#' Media <- medium_raw %>% tibble::column_to_rownames("Code")
#' ResM <- processing(
#'     data = Media[-c(40:45), -c(1:3)],
#'     metadata_sample = Media[-c(40:45), c(1:3)],
#'     metadata_info = c(
#'         Conditions = "Conditions",
#'         Biological_Replicates = "Biological_Replicates",
#'         core_norm_factor = "GrowthFactor",
#'         core_media = "blank"
#'     ),
#'     core = TRUE
#' )
#'
#' @importFrom dplyr mutate_all full_join
#' @importFrom magrittr %>% %<>%
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @export
processing <- function(
        data,
        metadata_sample = NULL,
        metadata_info = c(Conditions = "Conditions"),
        featurefilt = "Modified",
        cutoff_featurefilt = 0.8,
        tic = TRUE,
        mvi = TRUE,
        mvi_percentage = 50,
        hotellins_confidence = 0.99,
        core = FALSE,
        save_plot = "svg",
        save_table = "csv",
        print_plot = TRUE,
        path = NULL
) {
    ## ------------ Create log file ----------- ##
    metaproviz_init()

    ## ------------- Check SummarizedExperiment file ---------- ##
    input_data <- data
    if (inherits(data, "SummarizedExperiment"))  {
        log_info('Processing input SummarizedExperiment object.')
        se_list <- process_se(data)
        data <- se_list$data
        metadata_sample <- se_list$metadata_sample
    }

    ## ------------------ Check Input ------------------- ##
    # HelperFunction `check_param`
    check_param(
        data = data,
        metadata_sample = metadata_sample,
        metadata_feature = NULL,
        metadata_info = metadata_info,
        save_plot = save_plot,
        save_table = save_table,
        core = core,
        print_plot = print_plot
    )

    # HelperFunction `check_param` Specific
    check_param_processing(
        metadata_sample = metadata_sample,
        metadata_info = metadata_info,
        core = core,
        featurefilt = featurefilt,
        cutoff_featurefilt = cutoff_featurefilt,
        tic = tic,
        mvi = mvi,
        mvi_percentage = mvi_percentage,
        hotellins_confidence = hotellins_confidence
    )

    ## ------------------  Create output folders  and path ------------------- ##
    if (!is.null(save_plot) | !is.null(save_table) ) {
        folder <- save_path(
            folder_name = "Processing",
            path = path
        )

        Subfolder_P <- file.path(folder, "processing")
        if (!dir.exists(Subfolder_P)) {dir.create(Subfolder_P)}
    }

    ## ------------------ Prepare the data ------------------- ##
    # data files:
    data <- as.data.frame(data) %>%
    mutate_all(~ ifelse(grepl("^0*(\\.0*)?$", as.character(.)), NA, .)) # Make sure all 0 are changed to NAs

    data <- as.data.frame(mutate_all(as.data.frame(data), function(x) as.numeric(as.character(x))))

    ##############################################################################
    ## ------------------ 1. Feature filtering ------------------- ##
    if (!is.null(featurefilt)) {
        data_Filtered <- feature_filtering(
            data = data,
            featurefilt = featurefilt,
            cutoff_featurefilt = cutoff_featurefilt,
            metadata_sample = metadata_sample,
            metadata_info = metadata_info,
            core = core
        )

        data_Filt <- data_Filtered[["DF"]]
    } else {
        data_Filt <- data
    }

    ## ------------------ 2. Missing value Imputation ------------------- ##
    if (mvi) {
        mviRes <- mvi_imputation(
            data = data_Filt,
            metadata_sample = metadata_sample,
            metadata_info = metadata_info,
            core = core,
            mvi_percentage = mvi_percentage
        )
    } else {
        mviRes <- data_Filt
    }

    ## ------------------  3. total Ion Current Normalization ------------------- ##
    if (tic) {
        # Perform tic
        ticRes_List <- tic_norm(
            data = mviRes,
            metadata_sample = metadata_sample,
            metadata_info = metadata_info,
            tic = tic
        )
        ticRes <- ticRes_List[["DF"]][["data_tic"]]

        # Add plots to PlotList
        PlotList <- list()
        PlotList[["RLAPlot"]] <- ticRes_List[["Plot"]][["RLA_BeforeticNorm"]]
        PlotList[["RLAPlot_ticnorm"]] <- ticRes_List[["Plot"]][["RLA_AfterticNorm"]]
        PlotList[["RLAPlot_BeforeAfter_ticnorm"]] <- ticRes_List[["Plot"]][["norm_plots"]]
    } else {
        ticRes <- mviRes

        # Add plots to PlotList
        RLAPlot_List <- tic_norm(
            data = mviRes,
            metadata_sample = metadata_sample,
            metadata_info = metadata_info,
            tic = tic
        )
        PlotList <- list()
        PlotList[["RLAPlot"]] <- RLAPlot_List[["Plot"]][["RLA_BeforeticNorm"]]
    }

    ## ------------------ 4. core media QC (blank) and normalization ------------------- ##
    if (core) {
        data_coreNorm <- core_norm(
            data = ticRes,
            metadata_sample = metadata_sample,
            metadata_info = metadata_info
        )

        ticRes <- data_coreNorm[["DF"]][["core_Norm"]]
    }

    # ------------------ Final Output:
    data_norm <- ticRes %>% as.data.frame()

    ##########################################################################
    ## ------------------ Sample outlier identification ------------------- ##
    OutlierRes <-  outlier_detection(
        data = data_norm,
        metadata_sample = metadata_sample,
        metadata_info = metadata_info,
        core = core,
        hotellins_confidence = hotellins_confidence
    )

    ##########################################################################
    ## ------------------ Return ------------------- ##
    ## ---- DFs
    DFList <-
        data %>%
        rownames_to_column('Sample') %>%
        full_join(
            metadata_sample %>% rownames_to_column('Sample'),
            by = 'Sample'
        ) %>%
        relocate(!!sym(metadata_info[["Conditions"]]), .before = 1L) %>%
        list(data_Rawdata = .)

    if (!is.null(featurefilt)) {  # Add metabolites that where removed as part of the feature filtering

        if (length(data_Filtered[["RemovedMetabolites"]]) == 0) {

            DFList$Filtered_metabolites <-
                as.data.frame(
                    list(
                        feature_filtering = c(featurefilt),
                        cutoff_featurefilt = c(cutoff_featurefilt),
                        RemovedMetabolites = c("None")
                    )
                )

        } else {

            DFList$Filtered_metabolites <-
                as.data.frame(
                    list(
                        feature_filtering = rep(featurefilt, length(data_Filtered[["RemovedMetabolites"]])),
                        cutoff_featurefilt = rep(cutoff_featurefilt, length(data_Filtered[["RemovedMetabolites"]])),
                        RemovedMetabolites = data_Filtered[["RemovedMetabolites"]]
                    )
                )

        }

    }

    DFList$Preprocessing_output <- OutlierRes[["DF"]][["data_outliers"]]
    DFList$data_Rawdata %<>% column_to_rownames('Sample')

    if (core) {
        if (is.null(data_coreNorm[["DF"]][["Contigency_table_core_blank"]])) {
            DFList_core <- list( "CV_core_blank" = data_coreNorm[["DF"]][["CV_core_blank"]])
        } else {
            DFList_core <-
                list(
                    "CV_core_blank" = data_coreNorm[["DF"]][["CV_core_blank"]],
                    "Variation_ContigencyTable_core_blank" = data_coreNorm[["DF"]][["Contigency_table_core_blank"]]
                )
        }
        DFList <- c(DFList, DFList_core)
    }

    ## ---- Plots
    if (tic) {
        PlotList <- c(ticRes_List[["Plot"]], OutlierRes[["Plot"]])
    } else {
        PlotList <- c(RLAPlot_List[["Plot"]], OutlierRes[["Plot"]])
    }

    if (core) {
        PlotList <- c(PlotList , data_coreNorm[["Plot"]])
    }

    Res_List <- list("DF" = DFList ,"Plot" = PlotList)

    # Save Plots and DFs
    # As row names are not saved we need to make row.names to column for the DFs that needs this:
    DFList[["data_Rawdata"]] %<>% rownames_to_column("Code")
    DFList[["Preprocessing_output"]] %<>% rownames_to_column("Code")

    save_res(inputlist_df = DFList,
        inputlist_plot = PlotList,
        save_table = save_table,
        save_plot = save_plot,
        path = Subfolder_P,
        file_name = "processing",
        core = core,
        print_plot = print_plot
    )


    # Return
    # If the input was an se file, we also want to return a se file:
    if (inherits(input_data, "SummarizedExperiment"))  {
        PreRes <- Res_List[["DF"]][["Preprocessing_output"]]

        # find the "Outliers" column and split
        out_idx <- which(colnames(PreRes) == "Outliers")
        coldata_df <- PreRes[, seq_len(out_idx), drop = FALSE]
        assay_df <- PreRes[, (out_idx + 1):ncol(PreRes), drop = FALSE]

        # coerce assay to numeric matrix
        assay_mat <- as.matrix(t(assay_df))

        # build the SummarizedExperiment
        se_new <- SummarizedExperiment(
            assays = list(data = assay_mat),
            colData = DataFrame(coldata_df)
        )

        #Make list of se files:
        Se_list <-
            list(
                "data_Rawdata" = input_data,
                "Preprocessing_output" = se_new
            )
        Res_List <-
            list(
                "SE" = Se_list,
                "DF" = Res_List[["DF"]],
                "Plot" = Res_List[["Plot"]]
            )
        #Return
        invisible(return(Res_List))
    }else{
        invisible(return(Res_List))
    }
}

#
# Merge analytical replicates function   # #
#

#' Merges the analytical replicates of an experiment
#'
#' @param data SummarizedExperiment (se) file including assay and colData. If se file
#'     is provided, metadata_sample is extracted from the colData of the se
#'     object. Alternatively provide a DF which contains unique sample
#'     identifiers as row names and metabolite numerical values in columns with
#'     metabolite identifiers as column names. Use NA for metabolites that were
#'     not detected.
#' @param metadata_sample \emph{Optional: } Only required if you did not provide se file in
#'     parameter data. Provide DF which contains information about the samples,
#'     which will be combined with the input data based on the unique sample
#'     identifiers used as rownames. Must contain column with Conditions. If
#'     you do not have multiple conditions in your experiment assign all
#'     samples into the same condition. \strong{Default = NULL}
#' @param metadata_info \emph{Optional: } Named vector including the Conditions and Replicates
#'     information: c(Conditions="ColumnNameConditions",
#'     Biological_Replicates="ColumnName_metadata_sample",
#'     Analytical_Replicates="ColumnName_metadata_sample").\strong{Default =
#'     NULL}
#' @param save_table \emph{Optional: } File types for the analysis results are: "csv",
#'     "xlsx", "txt", ot NULL \strong{default: "csv"}
#' @param path \emph{Optional:} Path to the folder the results should be saved at.
#'     \strong{default: NULL}
#'
#' @return DF with the merged analytical replicates
#'
#' @examples
#' data(intracell_raw_se)
#' Res <- replicate_sum(data = intracell_raw_se)
#'
#' data(intracell_raw)
#' Intra <- intracell_raw %>% tibble::column_to_rownames("Code")
#' Res <- replicate_sum(
#'     data = Intra[-c(49:58), -c(1:3)],
#'     metadata_sample = Intra[-c(49:58), c(1:3)]
#' )
#'
#' @importFrom dplyr mutate_all summarise_all select rename ungroup group_by
#' @importFrom magrittr %>% %<>%
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom rlang !! :=
#' @importFrom tidyr unite
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @export
replicate_sum <- function(
        data,
        metadata_sample = NULL,
        metadata_info = c(Conditions = "Conditions", Biological_Replicates = "Biological_Replicates", Analytical_Replicates = "Analytical_Replicates"),
        save_table = "csv",
        path = NULL
) {

    # NSE vs. R CMD check workaround
    Biological_Replicates <- Conditions <- Analytical_Replicates <- UniqueID <- NULL
    ## ------------ Create log file ----------- ##
    metaproviz_init()


    ## ------------- Check SummarizedExperiment file ---------- ##
    input_data <- data
    if (inherits(data, "SummarizedExperiment"))  {
        log_info('Processing input SummarizedExperiment object.')
        se_list <- process_se(data)
        data <- se_list$data
        metadata_sample <- se_list$metadata_sample
    }

    ## ------------------ Check Input ------------------- ##
    # HelperFunction `check_param`
    check_param(
        data = data,
        metadata_sample = metadata_sample,
        metadata_feature = NULL,
        metadata_info = metadata_info,
        save_plot = NULL,
        save_table = save_table,
        core = FALSE,
        print_plot = FALSE
    )

    # `check_param` Specific
    if (metadata_info[["Conditions"]] %in% colnames(metadata_sample)) {
        # Conditions <- data[[metadata_info[["Conditions"]] ]]
    } else {
        stop("Column `Conditions` is required.")
    }
    if (metadata_info[["Biological_Replicates"]] %in% colnames(metadata_sample)) {
        # Biological_Replicates <- data[[metadata_info[["Biological_Replicates"]]]]
    } else {
        stop("Column `Biological_Replicates` is required.")
    }
    if (metadata_info[["Analytical_Replicates"]] %in% colnames(metadata_sample)) {
        # Analytical_Replicates <- data[[metadata_info[["Analytical_Replicates"]]]]
    } else {
        stop("Column `Analytical_Replicates` is required.")
    }

    ## ------------ Create Results output folder ----------- ##
    if (!is.null(save_table) ) {
        folder <- save_path(
            folder_name = "Processing",
            path = path
        )
        Subfolder <- file.path(folder, "replicate_sum")
        if (!dir.exists(Subfolder)) {dir.create(Subfolder)}
    }

    ## ------------  Load data and process  ----------- ##
    Input <- merge(x = metadata_sample %>% select(!!metadata_info[["Conditions"]], !!metadata_info[["Biological_Replicates"]], !!metadata_info[["Analytical_Replicates"]]),
    y = data,
    by = "row.names") %>%
    column_to_rownames("Row.names") %>%
    rename(
            "Conditions" = metadata_info[["Conditions"]],
            "Biological_Replicates" = metadata_info[["Biological_Replicates"]],
            "Analytical_Replicates" = metadata_info[["Analytical_Replicates"]]
        )

    # Make the replicate Sums
    Input_data_numeric_summed <- as.data.frame(Input %>%
    group_by(Biological_Replicates, Conditions) %>%
    summarise_all("mean") %>% select(-Analytical_Replicates))

    # Make a number of merged replicates column
    nReplicates <-  Input %>%
    group_by(Biological_Replicates, Conditions) %>%
    summarise_all("max") %>%
    ungroup() %>%
    select(Analytical_Replicates, Biological_Replicates, Conditions) %>%
    rename("n_AnalyticalReplicates_Summed" = "Analytical_Replicates")

    Input_data_numeric_summed <-
        merge(
            nReplicates,
            Input_data_numeric_summed,
            by = c("Conditions","Biological_Replicates")
        ) %>%
        unite(
            UniqueID,
            c("Conditions","Biological_Replicates"),
            sep = "_",
            remove = FALSE
        ) %>% # Create a uniqueID
        column_to_rownames("UniqueID")# set UniqueID to rownames

    # --------------- return ------------------ # #
    save_res(
        inputlist_df = list("Sum_AnalyticalReplicates" = Input_data_numeric_summed %>% tibble::rownames_to_column("Code")),
        inputlist_plot = NULL,
        save_table = save_table,
        save_plot = NULL,
        path = Subfolder,
        file_name = "Sum_AnalyticalReplicates",
        core = FALSE,
        print_plot = FALSE
    )

    # Return
    if (inherits(input_data, "SummarizedExperiment")) {
        PreRes <- Input_data_numeric_summed %>% rownames_to_column("Code")

        # find the "Outliers" column and split
        out_idx <- which(colnames(PreRes) == "n_AnalyticalReplicates_Summed")
        coldata_df <- PreRes[, seq_len(out_idx), drop = FALSE]
        assay_df <- PreRes[, (out_idx + 1):ncol(PreRes), drop = FALSE]
        rownames(assay_df) <- rownames(Input_data_numeric_summed)

        # coerce assay to numeric matrix
        assay_mat <- as.matrix(t(assay_df))

        # build the SummarizedExperiment
        se_new <- SummarizedExperiment(
            assays = list(data = assay_mat),
            colData = DataFrame(coldata_df)
        )

        invisible(return(se_new))
    } else {
        invisible(return(Input_data_numeric_summed))
    }

}


#
# Metabolite detection estimation using pool samples
#

#' Find metabolites with high variability across total pool samples
#'
#' @param data SummarizedExperiment (se) file including assay and colData. If se file
#'     is provided, metadata_sample is extracted from the colData of the se
#'     object. Alternatively provide a DF which contains unique sample
#'     identifiers as row names and metabolite numerical values in columns with
#'     metabolite identifiers as column names. Use NA for metabolites that were
#'     not detected.
#' @param metadata_sample \emph{Optional: } Only required if you did not provide se file in
#'     parameter data. Provide DF which contains information about the samples,
#'     which will be combined with the input data based on the unique sample
#'     identifiers used as rownames. Must contain column with Conditions. If
#'     you do not have multiple conditions in your experiment assign all
#'     samples into the same condition. \strong{Default = NULL}
#' @param metadata_info \emph{Optional: } NULL or Named vector including the Conditions and
#'     PoolSample information (Name of the Conditions column and Name of the
#'     pooled samples in the Conditions in the Input_SettingsFile)  :
#'     c(Conditions="ColumnNameConditions, PoolSamples=NamePoolCondition. If no
#'     Conditions is added in the Input_metadata_info, it is assumed that the
#'     conditions column is named 'Conditions' in the Input_SettingsFile. ).
#'     \strong{Default = NULL}
#' @param cutoff_cv \emph{Optional: } Filtering cutoff for high variance metabolites using
#'     the Coefficient of Variation. \strong{Default = 30}
#' @param save_plot \emph{Optional: } Select the file type of output plots. Options are svg,
#'     png, pdf or NULL. \strong{Default = svg}
#' @param save_table \emph{Optional: } File types for the analysis results are: "csv",
#'     "xlsx", "txt", ot NULL \strong{default: "csv"}
#' @param print_plot \emph{Optional: } If TRUE prints an overview of resulting plots.
#'     \strong{Default = TRUE}
#' @param path \emph{Optional:} Path to the folder the results should be saved at.
#'     \strong{default: NULL}
#'
#' @return List with two elements: DF (including input and output table) and Plot
#'     (including all plots generated)
#'
#' @examples
#' data(intracell_raw_se)
#' Res <- pool_estimation(
#'     data = intracell_raw_se,
#'     metadata_info = c(PoolSamples = "Pool", Conditions = "Conditions")
#' )
#'
#' data(intracell_raw)
#' Intra <- intracell_raw %>% tibble::column_to_rownames("Code")
#' Res <- pool_estimation(
#'     data = Intra[, -c(1:3)],
#'     metadata_sample = Intra[, c(1:3)],
#'     metadata_info = c(PoolSamples = "Pool", Conditions = "Conditions")
#' )
#'
#' @importFrom dplyr case_when select rowwise mutate
#' @importFrom magrittr %>% %<>%
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom logger log_info log_trace log_error
#' @importFrom ggplot2 after_stat
#' @importFrom ggrepel geom_text_repel
#' @importFrom grDevices dev.new dev.off
#' @importFrom stats sd
#' @importFrom SummarizedExperiment colData SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @export
pool_estimation <- function(
        data,
        metadata_sample = NULL,
        metadata_info = NULL,
        cutoff_cv = 30,
        save_plot = "svg",
        save_table = "csv",
        print_plot = TRUE,
        path = NULL
) {

    # NSE vs. R CMD check workaround
    CV <- Sample_type <- density <- HighVar <- NULL

    ## ------------ Create log file ----------- ##
    metaproviz_init()

    log_info('Starting pool estimation.')

    ## ------------- Check SummarizedExperiment file ---------- ##
    input_data <- data
    if (inherits(data, "SummarizedExperiment"))  {
        log_info('Processing input SummarizedExperiment object.')
        se_list <- process_se(data)
        data <- se_list$data
        metadata_sample <- se_list$metadata_sample
    }

    ## ------------------ Check Input ------------------- ##
    # HelperFunction `check_param`
    check_param(
        data = data,
        metadata_sample = metadata_sample,
        metadata_feature = NULL,
        metadata_info = metadata_info,
        save_plot = save_plot,
        save_table = save_table,
        core = FALSE,
        print_plot = print_plot
    )

    # `check_param` Specific
    if (!is.null(metadata_sample)) {
        if ("Conditions" %in% names(metadata_info)) {
            if (!(metadata_info[["Conditions"]] %in% colnames(metadata_sample)) ) {
                msg <-
                    sprintf(
                        "You have chosen Conditions = %s, %s was not found in metadata_sample as column. Please insert the name of the experimental conditions as stated in the metadata_sample.",
                        metadata_info[["Conditions"]],
                        metadata_info[["Conditions"]]
                    )
                log_error(msg)
                stop(msg, call. = FALSE)
            }
        }
        if ("PoolSamples" %in% names(metadata_info)) {
            if (!(metadata_info[["PoolSamples"]] %in% metadata_sample[[metadata_info[["Conditions"]]]]) ) {
                msg <-
                    sprintf(
                        "You have chosen PoolSamples = %s, %s was not found in metadata_sample as sample condition. Please insert the name of the pool samples as stated in the Conditions column of the metadata_sample.",
                        metadata_info[["PoolSamples"]],
                        metadata_info[["PoolSamples"]]
                    )
                log_error(msg)
                stop(msg, call. = FALSE)
            }
        }
    }

    if (!is.numeric(cutoff_cv) | cutoff_cv < 0) {
        stop("Check input. The selected cutoff_cv value should be a positive numeric value.")
    }

    ## ------------------  Create output folders  and path ------------------- ##
    if (!is.null(save_plot) | !is.null(save_table) ) {
        folder <- save_path(
            folder_name = "Processing",
            path = path
        )

        Subfolder <- file.path(folder, "pool_estimation")
        log_info('Selected output directory: `%s`.', Subfolder)
        if (!dir.exists(Subfolder)) {
            log_trace('Creating directory: `%s`.', Subfolder)
            dir.create(Subfolder)
        }
    }

    ## ------------------ Prepare the data ------------------- ##
    # data files:
    if (is.null(metadata_sample)) {
        Pooldata <- data %>%
        mutate_all(~ ifelse(grepl("^0*(\\.0*)?$", as.character(.)), NA, .)) # Make sure all 0 are changed to NAs
    } else {
        Pooldata <- data[metadata_sample[[metadata_info[["Conditions"]]]] == metadata_info[["PoolSamples"]],] %>%
        mutate_all(~ ifelse(grepl("^0*(\\.0*)?$", as.character(.)), NA, .)) # Make sure all 0 are changed to NAs
    }

    #####################################################################
    ## ------------------ Coefficient of Variation ------------------- ##

    log_trace('Calculating coefficient of variation.')
    result_df <-
        apply(
            Pooldata,
            2,
            function(x) { (sd(x, na.rm = TRUE)/  mean(x, na.rm = TRUE))*100 }  ) %>% t() %>% as.data.frame(
        )
    rownames(result_df)[1] <- "CV"

    NAvector <-
        apply(
            Pooldata,
            2,
            function(x) {(sum(is.na(x))/length(x))*100 }
        )# Calculate the NAs

    # Create Output DF
    result_df_final <- result_df %>%
    t() %>% as.data.frame() %>% rowwise() %>%
    mutate(HighVar = CV > cutoff_cv) %>% as.data.frame()

    result_df_final$MissingValuepercentage <- NAvector

    rownames(result_df_final)<- colnames(data)
    result_df_final_out <- rownames_to_column(result_df_final,"Metabolite" )

    # Remove Metabolites from data based on cutoff_cv
    log_trace('Applying CV cut-off.')
    if (!is.null(metadata_sample)) {
        unstable_metabs <- rownames(result_df_final)[result_df_final[["HighVar_Metabs"]]]
        if (length(unstable_metabs)>0) {
            filtered_Input_data <- data %>% select(!unstable_metabs)
        } else {
            filtered_Input_data <- NULL
        }
    } else {
        filtered_Input_data <- NULL
    }

    ## ------------------ QC plots ------------------- ##
    # Start QC plot list
    log_info('Plotting QC plots.')
    PlotList <- list()

    # 1. Pool Sample PCA
    log_trace('Pool sample PCA.')
    dev.new()
    if (is.null(metadata_sample)) {
        pca_data <- Pooldata
        pca_QC_pool <- invisible(
            viz_pca(data = pca_data,
                    plot_name = "QC Pool samples",
                    save_plot = NULL)
        )
    } else {
        pca_data <-
            merge(
                metadata_sample %>% select(metadata_info[["Conditions"]]),
                data,
                by = 0
            ) %>%
            column_to_rownames("Row.names") %>%
            mutate(
                Sample_type = case_when(.data[[metadata_info[["Conditions"]]]] == metadata_info[["PoolSamples"]] ~ "Pool",
                                        TRUE ~ "Sample")
            )

        pca_QC_pool <- invisible(viz_pca(data = pca_data %>% select(-all_of(metadata_info[["Conditions"]]), -Sample_type),
        metadata_info = c(color = "Sample_type"),
        metadata_sample = pca_data,
        plot_name = "QC Pool samples",
        save_plot = NULL))
    }
    dev.off()
    PlotList [["PCAPlot_PoolSamples"]] <- pca_QC_pool[["Plot_Sized"]][["Plot_Sized"]]


    # 2. Histogram of CVs
    log_trace('CV histogram.')
    HistCV <- invisible(ggplot(result_df_final_out, aes(CV)) +
                            geom_histogram(
                                aes(y = after_stat(density)),
                                color = "black",
                                fill = "white"
                            )+
                            geom_vline(aes(xintercept = cutoff_cv),
                            color = "darkred", linetype = "dashed", size = 1)+
                            geom_density(alpha = .2, fill = "#FF6666") +
                            labs(
                                title = "CV for metabolites of Pool samples",
                                x = "Coefficient of variation (CV%)",
                                y = "Frequency"
                            )+
                            theme_classic())

    HistCV_Sized <-
        plotGrob_Processing(
            input_plot = HistCV,
            plot_name = "CV for metabolites of Pool samples",
            plot_type = "Hist"
        )
    PlotList [["Histogram_CV-PoolSamples"]] <- HistCV_Sized

    # 2. ViolinPlot of CVs
    log_trace('CV violin plot.')
    # Make Violin of CVs
    Plot_cv_result_df <- result_df_final_out %>%
    mutate(HighVar = ifelse((CV > cutoff_cv), paste("> CV", cutoff_cv, sep = ""), paste("< CV", cutoff_cv, sep = "")))

    ViolinCV <- invisible(ggplot( Plot_cv_result_df, aes(y = CV, x = HighVar, label = Plot_cv_result_df$Metabolite))+
    geom_violin(alpha = 0.5 , fill = "#FF6666")+
    geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5) +
    geom_text_repel(
        aes(label = ifelse(Plot_cv_result_df$CV > cutoff_cv,
        as.character(Plot_cv_result_df$Metabolite), '')),
        hjust = 0,
        vjust = 0,
        box.padding = 0.5,
        # space between text and point
        point.padding = 0.5,
        # space around points
        max.overlaps = Inf
        ) + # allow for many labels
        labs(
            title = "CV for metabolites of Pool samples",
            x = "Metabolites",
            y = "Coefficient of variation (CV%)"
        )+
        theme_classic())

    ViolinCV_Sized <-
        plotGrob_Processing(
            input_plot = ViolinCV,
            plot_name = "CV for metabolites of Pool samples",
            plot_type = "Violin"
        )

    PlotList [["ViolinPlot_CV-PoolSamples"]] <- ViolinCV_Sized

    ############################################################
    ## ------------------ Return and Save ------------------- ##
    # Save
    log_info('Preparing saved and returned data.')
    if (!is.null(filtered_Input_data)) {
        DF_list <-
            list(
                "data" = data,
                "Filtered_data" = filtered_Input_data,
                "CV" = result_df_final_out
            )
    } else {
        DF_list <- list("data" = data, "CV" = result_df_final_out)
    }
    ResList <- list("DF" = DF_list,"Plot" = PlotList)

    # Save
    DF_list[["data"]]<-  DF_list[["data"]] %>% tibble::rownames_to_column("Code")

    log_info(
        'Saving results: [save_table=%s, save_plot=%s, path=%s].',
        save_table,
        save_plot,
        Subfolder
    )
    save_res(
        inputlist_df = DF_list,
        inputlist_plot = PlotList,
        save_table = save_table,
        save_plot = save_plot,
        path = Subfolder,
        file_name = "pool_estimation",
        core = FALSE,
        print_plot = print_plot
    )

    log_info('Finished pool estimation.')

    if (inherits(input_data, "SummarizedExperiment"))  {
        assay_mat <- as.matrix(t(as.data.frame(ResList[["DF"]][["data"]])))

        coldata_df <- as.data.frame(colData(input_data))
        rowdata_df <- as.data.frame(ResList[["DF"]][["CV"]])%>% tibble::column_to_rownames("Metabolite")

        # build the SummarizedExperiment
        se_new <- SummarizedExperiment(
            assays = list(data = assay_mat),
            colData = DataFrame(coldata_df),
            rowData = DataFrame(rowdata_df)
        )

        #Add to ResList:
        ResList <-
            list(
                "SE" = se_new,
                "DF" = ResList[["DF"]],
                "Plot" = ResList[["Plot"]]
            )

        #Return
        invisible(return(ResList))

    } else {
        invisible(return(ResList))
    }

}


#
# processing helper function: feature_filtering
#

#' Filter features
#'
#' @param data DF which contains unique sample identifiers as row names and metabolite
#'     numerical values in columns with metabolite identifiers as column names.
#'     Use NA for metabolites that were not detected and consider converting
#'     any zeros to NA unless they are true zeros.
#' @param metadata_sample DF which contains information about the samples, which will be combined
#'     with the input data based on the unique sample identifiers used as
#'     rownames.
#' @param metadata_info Named vector containing the information about the names of the
#'     experimental parameters. c(Conditions="ColumnName_Plot_SettingsFile",
#'     Biological_Replicates="ColumnName_Plot_SettingsFile"). Column
#'     "Conditions" with information about the sample conditions (e.g. "N" and
#'     "T" or "Normal" and "Tumor"), can be used for feature filtering and
#'     colour coding in the PCA. Column "BiologicalReplicates" including
#'     numerical values. For core = TRUE add core_media =
#'     "Columnname_Input_SettingsFile", which specifies the name of the media
#'     controls in the Conditions.
#' @param core \emph{Optional: } If TRUE, a consumption-release experiment has been
#'     performed.Should not be normalised to media blank. Provide information
#'     about control media sample names via metadata_info "core_media" samples.
#'     \strong{Default = FALSE}
#' @param featurefilt \emph{Optional: } If NULL, no feature filtering is performed. If set to
#'     "Standard" then it applies the 80%-filtering rule (Bijlsma S. et al.,
#'     2006) on the metabolite features on the whole dataset. If is set to
#'     "Modified",filtering is done based on the different conditions, thus a
#'     column named "Conditions" must be provided in the Input_SettingsFile
#'     input file including the individual conditions you want to apply the
#'     filtering to (Yang, J et al., 2015). \strong{Default = Modified}
#' @param cutoff_featurefilt \emph{Optional: } percentage of feature filtering. \strong{Default =
#'     0.8}
#'
#' @return List with two elements: filtered matrix  and features filtered
#'
#' @examples
#' data(intracell_raw)
#' Intra <- intracell_raw %>% tibble::column_to_rownames("Code")
#' Res <- feature_filtering(
#'     data =
#'         Intra[-c(49:58), -c(1:3)] %>%
#'         dplyr::mutate_all(
#'             ~ ifelse(grepl("^0*(\\.0*)?$", as.character(.)), NA, .)
#'         ),
#'     metadata_sample = Intra[-c(49:58), c(1:3)],
#'     metadata_info = c(
#'         Conditions = "Conditions",
#'         Biological_Replicates = "Biological_Replicates"
#'     )
#' )
#'
#' @importFrom dplyr filter mutate_all
#' @importFrom magrittr %>% %<>%
#' @importFrom logger log_info log_trace
#' @noRd
feature_filtering <- function(
        data,
        metadata_sample,
        metadata_info,
        core = FALSE,
        featurefilt = "Modified",
        cutoff_featurefilt = 0.8
) {
    ## ------------ Create log file ----------- ##
    metaproviz_init()

    ## ------------------ Prepare the data ------------------- ##
    feat_filt_data <- as.data.frame(data) %>%
    mutate_all(~ ifelse(grepl("^0*(\\.0*)?$", as.character(.)), NA, .)) # Make sure all 0 are changed to NAs

    if (core) { # remove core_media samples for feature filtering
        feat_filt_data <- feat_filt_data %>% filter(!metadata_sample[[metadata_info[["Conditions"]]]] == metadata_info[["core_media"]])
        Feature_Filtering <- paste0(featurefilt, "_core")
    }

    ## ------------------ Perform filtering ------------------ ##
    if (featurefilt ==  "Modified") {
        message <-
            paste0(
                "feature_filtering: Here we apply the modified 80%-filtering rule that takes the class information (Column `Conditions`) into account, which additionally reduces the effect of missing values (REF: Yang et. al., (2015), doi: 10.3389/fmolb.2015.00004). ",
                "Filtering value selected: ",
                cutoff_featurefilt,
                sep = ""
            )
        log_info(message)
        message(message)
        if (core) {
            feat_filt_Conditions <- metadata_sample[[metadata_info[["Conditions"]]]][!metadata_sample[[metadata_info[["Conditions"]]]] == metadata_info[["core_media"]]]
        } else {
            feat_filt_Conditions <- metadata_sample[[metadata_info[["Conditions"]]]]
        }

        if (is.null(unique(feat_filt_Conditions))) {
            message("Conditions information is missing.")
            log_trace(message)
            stop(message)
        }
        if (length(unique(feat_filt_Conditions)) ==  1) {
            message("to perform the Modified feature filtering there have to be at least 2 different Conditions in the `Condition` column in the Experimental design. Consider using the Standard feature filtering option.")
            log_trace(message)
            stop(message)
        }

        miss <- c()
        split_Input <-
            split(
                feat_filt_data,
                feat_filt_Conditions
            ) # split data frame into a list of dataframes by condition

        for (m in split_Input) { # Select metabolites to be filtered for different conditions
            for (i in seq_len(ncol(m))) {
                if (length(which(is.na(m[,i]))) > (1-cutoff_featurefilt)*nrow(m))
                    miss <- append(miss,i)
            }
        }

        if (length(miss) ==  0) { #remove metabolites if any are found
            message("There where no metabolites exluded")
            filtered_matrix <- data
            feat_file_res <- "There where no metabolites exluded"
        } else {
            names <- unique(colnames(data)[miss])
            msg <-
                sprintf(
                    "%d metabolites where removed: %s",
                    length(unique(miss)),
                    paste0(names, collapse = ", ")
                )
            log_info(msg)
            message(msg)
            filtered_matrix <- data[,-miss]
        }
    }else if (featurefilt ==  "Standard") {
        message <-
            paste0(
                "feature_filtering: Here we apply the so-called 80%-filtering rule, which removes metabolites with missing values in more than 80% of samples (REF: Smilde et. al. (2005), Anal. Chem. 77, 6729-6736., doi:10.1021/ac051080y). ",
                "Filtering value selected:",
                cutoff_featurefilt
            )
        log_info(message)
        message(message)

        split_Input <- feat_filt_data

        miss <- c()
        for (i in seq_len(ncol(split_Input))) { # Select metabolites to be filtered for one condition
            if (length(which(is.na(split_Input[,i]))) > (1-cutoff_featurefilt)*nrow(split_Input))
                miss <- append(miss,i)
        }

        if (length(miss) ==  0) { #remove metabolites if any are found
            message <- paste0("feature_filtering: There where no metabolites exluded")
            log_info(message)
            message(message)

            filtered_matrix <- data
            feat_file_res <- "There where no metabolites exluded"
        } else {
            names <- unique(colnames(data)[miss])
            message <-
                paste0(
                    length(unique(miss)),
                    " metabolites where removed: ",
                    paste0(names, collapse = ", ")
                )
            log_info(message)
            message(message)
            filtered_matrix <- data[,-miss]
        }
    }

    ## ------------------ Return ------------------ ##
    features_filtered <- unique(colnames(data)[miss]) %>% as.vector()
    filtered_matrix <- as.data.frame(mutate_all(as.data.frame(filtered_matrix), function(x) as.numeric(as.character(x))))

    Filtered_results <-
        list(
            "DF" = filtered_matrix,
            "RemovedMetabolites" = features_filtered
        )
    invisible(return(Filtered_results))
}


#
# processing helper function: Missing Value imputation
#

#' Missing Value Imputation using half minimum value
#'
#' @param data DF which contains unique sample identifiers as row names and metabolite
#'     numerical values in columns with metabolite identifiers as column names.
#'     Use NA for metabolites that were not detected and consider converting
#'     any zeros to NA unless they are true zeros.
#' @param metadata_sample DF which contains information about the samples, which will be combined
#'     with the input data based on the unique sample identifiers used as
#'     rownames.
#' @param metadata_info Named vector containing the information about the names of the
#'     experimental parameters. c(Conditions="ColumnName_Plot_SettingsFile",
#'     Biological_Replicates="ColumnName_Plot_SettingsFile", core_media =
#'     "Columnname_Input_SettingsFile"). Column "Conditions" with information
#'     about the sample conditions, Column "BiologicalReplicates" including
#'     numerical values and Column "Columnname_Input_SettingsFile" is used to
#'     specify the name of the media controls in the Conditions.
#' @param core \emph{Optional: } If TRUE, a consumption-release experiment has been
#'     performed. Should not be normalised to media blank. Provide information
#'     about control media sample names via metadata_info "core_media"
#'     samples.\strong{Default = FALSE}
#' @param mvi_percentage \emph{Optional: } percentage 0-100 of imputed value based on the minimum
#'     value. \strong{Default = 50}
#'
#' @return DF with imputed values
#'
#' @examples
#' data(intracell_raw)
#' Intra <- intracell_raw %>% tibble::column_to_rownames("Code")
#' Res <- mvi_imputation(
#'     data =
#'         Intra[-c(49:58), -c(1:3)] %>%
#'         dplyr::mutate_all(
#'             ~ ifelse(grepl("^0*(\\.0*)?$", as.character(.)), NA, .)
#'         ),
#'     metadata_sample = Intra[-c(49:58), c(1:3)],
#'     metadata_info = c(
#'         Conditions = "Conditions",
#'         Biological_Replicates = "Biological_Replicates"
#'     )
#' )
#'
#' @importFrom dplyr select mutate group_by filter mutate_all
#' @importFrom magrittr %>% %<>%
#' @importFrom tibble column_to_rownames
#' @importFrom logger log_info
#' @importFrom tidyselect everything where
#' @noRd
mvi_imputation <- function(
        data,
        metadata_sample,
        metadata_info,
        core = FALSE,
        mvi_percentage = 50
) {

    # NSE vs. R CMD check workaround
    Conditions <- NULL
    ## ------------------ Prepare the data ------------------- ##
    filtered_matrix <- data %>%
    mutate_all(~ ifelse(grepl("^0*(\\.0*)?$", as.character(.)), NA, .)) # Make sure all 0 are changed to NAs

    ## ------------------ Perform mvi ------------------ ##
    # Do mvi for the samples
    message <- paste0("Missing Value Imputation: Missing value imputation is performed, as a complementary approach to address the missing value problem, where the missing values are imputing using the `half minimum value`. REF: Wei et. al., (2018), Reports, 8, 663, doi:https://doi.org/10.1038/s41598-017-19120-0")
    log_info(message)
    message(message)

    if (core) { # remove blank samples
        NA_removed_matrix <- filtered_matrix %>% filter(!metadata_sample[[metadata_info[["Conditions"]]]] == metadata_info[["core_media"]])

    } else {
        NA_removed_matrix <- filtered_matrix %>% as.data.frame()
    }

    NA_removed_matrix %<>%
        mutate(Conditions = if (core) { # If we have a CoRe experiment we need to remove the sample metainformation of the media blank samples!
            metadata_sample %>%
            filter(!metadata_sample[[metadata_info[["Conditions"]]]] == metadata_info[["core_media"]]) %>%
            select(!!sym(metadata_info[["Conditions"]]))
        } else {
            metadata_sample[[metadata_info[["Conditions"]]]]# If we have a standard experiments no samples need to be removed
        }) %>%
        group_by(Conditions) %>%
        mutate(
            across(
                .cols = everything(),
                .fns = ~ {
                    if (all(is.na(.x))) {
                        .x  # leaves as-is (here keep NAs)
                    } else {
                        ifelse(
                            is.na(.x),
                            min(.x, na.rm = TRUE) * mvi_percentage / 100,
                            .x
                        )
                    }
                },
                .names = "{.col}"
            )
        ) %>%
        ungroup() %>%
        select(-Conditions) %>%
        as.data.frame() %>%
        `rownames<-`(rownames(NA_removed_matrix))

    # Check for groups with only NAs
    na_check <-
        NA_removed_matrix %>%
        summarise(
            across(
                .cols = where(is.numeric),
                .fns = ~ any(is.na(.x)),
                .names = "{.col}"
            )
        ) %>%
        select(where(~ any(.)))

    if (ncol(na_check) > 0L) {

        msg <- sprintf(
            paste0(
                "Some features had all NA values in all samples of certain ",
                "conditions - hence no imputation was done and NA remains for: %s"
            ),
            paste(names(na_check), collapse = ", ")
        )

        log_info(msg)
        message(msg)

    }

    if (core) {
        replaceNAdf <- filtered_matrix %>% filter(metadata_sample[[metadata_info[["Conditions"]]]] == metadata_info[["core_media"]])

        # find metabolites with NA
        na_percentage <- colMeans(is.na(replaceNAdf)) * 100
        highNA_metabs <- na_percentage[na_percentage>20 & na_percentage<100]
        OnlyNA_metabs <- na_percentage[na_percentage == 100]

        # report metabolites with NA
        if (sum(na_percentage)>0) {
            message <- paste0("NA values were found in Control_media samples for metabolites. For metabolites including NAs mvi is performed unless all samples of a metabolite are NA.")
            log_info(message)
            message(message)
            if (sum(na_percentage>20 & na_percentage<100)>0) {
                message <-
                    paste0(
                        "Metabolites with high NA load (>20%) in Control_media samples are: ",
                        paste(names(highNA_metabs), collapse = ", "),
                        "."
                    )
                log_info(message)
                message(message)
            }
            if (sum(na_percentage == 100)>0) {
                message <-
                    paste0(
                        "Metabolites with only NAs (= 100%) in Control_media samples are: ",
                        paste(names(OnlyNA_metabs), collapse = ", "),
                        ". Those NAs are set zero as we consider them true zeros"
                    )
                log_info(message)
                message(message)
            }
        }

        # if all values are NA set to 0
        replaceNAdf_zero <- as.data.frame(lapply(replaceNAdf, function(x) if (all(is.na(x))) replace(x, is.na(x), 0) else x))
        colnames(replaceNAdf_zero) <-  colnames(replaceNAdf)
        rownames(replaceNAdf_zero) <-  rownames(replaceNAdf)

        # If there is at least 1 value use the half minimum per feature
        replaceNAdf_Zero_mvi <-
            apply(
                replaceNAdf_zero,
                2,
                function(x) {x[is.na(x)] <-  min(x, na.rm = TRUE
            )/2
        return(x)
        }) %>% as.data.frame()
        rownames(replaceNAdf_Zero_mvi) <-  rownames(replaceNAdf)

        # add the samples in the original dataframe
        filtered_matrix_res <- rbind(NA_removed_matrix, replaceNAdf_Zero_mvi)
    } else {
        filtered_matrix_res <- NA_removed_matrix
    }

    ## ------------------ Return ------------------ ##
    invisible(return(filtered_matrix_res))
}


#
# processing helper function: total ion Count Normalization
#

#' Total ion count normalization
#'
#' @param data DF which contains unique sample identifiers as row names and metabolite
#'     numerical values in columns with metabolite identifiers as column names.
#'     Use NA for metabolites that were not detected and consider converting
#'     any zeros to NA unless they are true zeros.
#' @param metadata_sample DF which contains information about the samples, which will be combined
#'     with the input data based on the unique sample identifiers used as
#'     rownames.
#' @param metadata_info Named vector containing the information about the names of the
#'     experimental parameters. c(Conditions="ColumnName_Plot_SettingsFile").
#' @param tic \emph{Optional: } If TRUE, total Ion Count normalization is performed.
#'     If FALSE, only RLA QC plots are returned. \strong{Default = TRUE}
#'
#' @return List with two elements: DF (including output table) and Plot (including
#'     all plots generated)
#'
#' @examples
#' data(intracell_raw)
#' Intra <- intracell_raw %>% tibble::column_to_rownames("Code")
#' Res <- tic_norm(
#'     data =
#'         Intra[-c(49:58), -c(1:3)] %>%
#'         dplyr::mutate_all(
#'             ~ ifelse(grepl("^0*(\\.0*)?$", as.character(.)), NA, .)
#'         ),
#'     metadata_sample = Intra[-c(49:58), c(1:3)],
#'     metadata_info = c(Conditions = "Conditions")
#' )
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom tidyr pivot_longer
#' @importFrom gridExtra grid.arrange
#' @importFrom ggplot2 ggplot geom_boxplot geom_hline labs theme_classic
#' @importFrom ggplot2 theme_minimal theme annotation_custom aes_string element_text
#' @importFrom logger log_info log_trace
#' @importFrom grDevices dev.new dev.off
#' @importFrom stats median
#' @importFrom tidyselect everything
#' @noRd
tic_norm <- function(
        data,
        metadata_sample,
        metadata_info,
        tic = TRUE
) {

    # NSE vs. R CMD check workaround
    median <- NULL
    ## ------------------ Prepare the data ------------------- ##
    NA_removed_matrix <- data
    NA_removed_matrix[is.na(NA_removed_matrix)] <- 0 # replace NA with 0

    ## ------------------ QC plot ------------------- ##
    # # # Before tic Normalization
    # # ## Log() transformation:
    log_NA_removed_matrix <- log(NA_removed_matrix) %>% t() %>% as.data.frame() # log tranforms the data
    nan_count <- sum(is.nan(as.matrix(log_NA_removed_matrix)))# Count NaN values (produced by log(0))
    if (nan_count > 0) {# Issue a custom warning if NaNs are present
        message <-
            paste(
                "For the RLA plot before/after tic normalisation we have to perform log() transformation. This resulted in",
                nan_count,
                "NaN values due to 0s in the data."
            )
        log_trace("Warning: ", message, sep = "")
        warning(message)
    }

    medians <- apply(log_NA_removed_matrix, 2, median) # get median
    RLA_data_raw <- log_NA_removed_matrix - medians   # Subtract the medians from each column
    RLA_data_long <-
        pivot_longer(
            RLA_data_raw,
            cols = everything(),
            names_to = "Group"
        )
    names(RLA_data_long)<- c("Samples", "Intensity")
    RLA_data_long <- as.data.frame(RLA_data_long)
    RLA_data_long <- RLA_data_long
    metadata_sample <- metadata_sample
    for (row in seq_len(nrow(RLA_data_long))) { # add conditions
        RLA_data_long[row, metadata_info[["Conditions"]]] <- metadata_sample[rownames(metadata_sample) %in%RLA_data_long[row,1],metadata_info[["Conditions"]]]
    }

    # Create the ggplot boxplot
    RLA_data_raw <-
        ggplot(
            RLA_data_long,
            aes_string(x = "Samples", y = "Intensity", color = metadata_info[["Conditions"]])
        ) +
        geom_boxplot() +
        geom_hline(yintercept = 0, color = "red", linetype = "solid") +
        labs(title = "Before tic Normalization")+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 90, hjust = 1))+
        theme(legend.position = "none")

    # RLA_data_raw_Sized <- plotGrob_Processing(input_plot = RLA_data_raw,plot_name= "Before tic Normalization", plot_type= "RLA")

    if (tic) {
        ## ------------------ Perform tic ------------------- ##
        message <- paste0("total Ion Count (tic) normalization: total Ion Count (tic) normalization is used to reduce the variation from non-biological sources, while maintaining the biological variation. REF: Wulff et. al., (2018), Advances in Bioscience and Biotechnology, 9, 339-351, doi:https://doi.org/10.4236/abb.2018.98022")
        log_info(message)
        message(message)

        RowSums <- rowSums(NA_removed_matrix)
        Median_RowSums <- median(RowSums) #This will built the median
        data_tic_Pre <-
            apply(
                NA_removed_matrix,
                2,
                function(i) i/RowSums
            ) #This is dividing the ion intensity by the total ion count
        data_tic <- data_tic_Pre*Median_RowSums #Multiplies with the median metabolite intensity
        data_tic <- as.data.frame(data_tic)

        ## ------------------ QC plot ------------------- ##
        # # # After tic normalization
        log_data_tic  <- log(data_tic) %>% t() %>% as.data.frame() # log tranforms the data
        medians <- apply(log_data_tic, 2, median)
        RLA_data_norm <- log_data_tic - medians   # Subtract the medians from each column
        RLA_data_long <-
            pivot_longer(
                RLA_data_norm,
                cols = everything(),
                names_to = "Group"
            )
        names(RLA_data_long)<- c("Samples", "Intensity")
        for (row in seq_len(nrow(RLA_data_long))) { # add conditions
            RLA_data_long[row, metadata_info[["Conditions"]]] <- metadata_sample[rownames(metadata_sample) %in%RLA_data_long[row,1],metadata_info[["Conditions"]]]
        }

        # Create the ggplot boxplot
        RLA_data_norm <-
            ggplot(
                RLA_data_long,
                aes_string(x = "Samples", y = "Intensity", color = metadata_info[["Conditions"]])
            ) +
            geom_boxplot() +
            geom_hline(yintercept = 0, color = "red", linetype = "solid") +
            labs(title = "After tic Normalization")+
            theme_classic() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))+
            theme(legend.position = "none")

        # RLA_data_norm_Sized <- plotGrob_Processing(input_plot = RLA_data_norm,plot_name= "After tic Normalization", plot_type= "RLA")

        # Combine Plots
        dev.new()
        norm_plots <- grid.arrange(RLA_data_raw+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+ theme(legend.position = "none"),
        RLA_data_norm+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+ theme(legend.position = "none"),
        ncol = 2)
        dev.off()
        norm_plots <- ggplot() +theme_minimal()+ annotation_custom(norm_plots)


        ## ------------------ Return ------------------ ##
        Output_list <-
            list(
                "DF" = list("data_tic" = as.data.frame(data_tic)),
                "Plot" = list( "norm_plots" = norm_plots, "RLA_AfterticNorm" = RLA_data_norm,  "RLA_BeforeticNorm" = RLA_data_raw )
            )
        invisible(return(Output_list))
    } else {
        ## ------------------ Return ------------------ ##
        Output_list <- list("Plot" = list("RLA_BeforeticNorm" = RLA_data_raw))
        invisible(return(Output_list))
    }
}


#
# processing helper function: core nomalisation
#

#' Normalize consumption release data
#'
#' @param data DF which contains unique sample identifiers as row names and metabolite
#'     numerical values in columns with metabolite identifiers as column names.
#'     Use NA for metabolites that were not detected and consider converting
#'     any zeros to NA unless they are true zeros.
#' @param metadata_sample DF which contains information about the samples, which will be combined
#'     with the input data based on the unique sample identifiers used as
#'     rownames.
#' @param metadata_info Named vector containing the information about the names of the
#'     experimental parameters. c(Conditions="ColumnName_Plot_SettingsFile",
#'     core_norm_factor = "Columnname_Input_SettingsFile", core_media =
#'     "Columnname_Input_SettingsFile"). Column core_norm_factor is used for
#'     normalization and core_media is used to specify the name of the media
#'     controls in the Conditions.
#'
#' @return List with two elements: DF (including output table) and Plot (including
#'     all plots generated)
#'
#' @examples
#' Media <-
#'     medium_raw %>%
#'     tibble::column_to_rownames("Code") %>%
#'     subset(!Conditions == "Pool") %>%
#'     dplyr::mutate_all(~ ifelse(grepl("^0*(\\.0*)?$", as.character(.)), NA, .))
#' Res <- core_norm(
#'     data = Media[, -c(1:3)],
#'     metadata_sample = Media[, c(1:3)],
#'     metadata_info = c(
#'         Conditions = "Conditions",
#'         core_norm_factor = "GrowthFactor",
#'         core_media = "blank"
#'     )
#' )
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggplot2 ggplot geom_histogram geom_vline geom_density labs
#' @importFrom ggplot2 after_stat theme_classic geom_violin geom_dotplot
#' @importFrom logger log_info log_trace
#' @importFrom dplyr case_when summarise_all mutate rowwise mutate_all select filter pull
#' @importFrom grDevices dev.new dev.off
#' @importFrom stats fisher.test sd var
#' @noRd
core_norm <- function(
        data,
        metadata_sample,
        metadata_info
) {

    # NSE vs. R CMD check workaround
    Sample_type <- CV <- density <- HighVar <- NULL
    ## ------------------ Prepare the data ------------------- ##
    data_tic <- data
    data_tic[is.na(data_tic)] <- 0

    ## ------------------ Perform QC ------------------- ##
    Conditions <- metadata_sample[[metadata_info[["Conditions"]]]]
    core_medias <- data_tic[grep(metadata_info[["core_media"]], Conditions),]

    if (dim(core_medias)[1] == 1) {
        message <- paste0("Only 1 core_media sample was found. Thus, the consistency of the core_media samples cannot be checked. It is assumed that the core_media samples are already summed.")
        log_trace(paste("Warning: ", message, sep = ""))
        warning(message)

        core_media_df <- core_medias %>% t() %>% as.data.frame()
        colnames(core_medias) <- "core_mediaMeans"
    } else {
        ##########################################################################
        ## ------------------ QC Plots
        PlotList <- list()
        # # -- 1. PCA Media_control
        media_pca_data <-
            merge(
                x = metadata_sample %>% select(metadata_info[["Conditions"]]),
                y = data_tic,
                by = 0
            ) %>%
            column_to_rownames("Row.names") %>%
            mutate(
                Sample_type = case_when(Conditions == metadata_info[["core_media"]] ~ "core_media",
                                        TRUE ~ "Sample")
            )

        media_pca_data[is.na( media_pca_data)] <- 0

        dev.new()
        pca_QC_media <- invisible(viz_pca(data = media_pca_data %>% dplyr::select(-metadata_info[["Conditions"]], -Sample_type),
        metadata_info = c(color = "Sample_type"),
        metadata_sample = media_pca_data,
        plot_name = "QC Media_samples",
        save_plot = NULL))
        dev.off()

        PlotList[["PCA_coreMediaSamples"]] <- pca_QC_media[["Plot_Sized"]][["QC Media_samples"]]

        # # -- 2. Metabolite Variance Histogram
        # Coefficient of Variation
        result_df <-
            apply(
                core_medias,
                2,
                function(x) { (sd(x, na.rm = TRUE)/  mean(x, na.rm = TRUE))*100 } ) %>% t() %>% as.data.frame(
            )
        result_df[1, is.na(result_df[1,])]<- 0
        rownames(result_df)[1] <- "CV"

        cutoff_cv <- 30
        result_df <- result_df %>% t() %>% as.data.frame() %>% rowwise() %>%
        mutate(HighVar = CV > cutoff_cv) %>% as.data.frame()
        rownames(result_df)<- colnames(core_medias)

        # calculate the NAs
        NAvector <-
            apply(
                core_medias,
                2,
                function(x) { (sum(is.na(x))/length(x))*100 }
            )
        result_df$MissingValuepercentage <- NAvector

        cv_result_df <- result_df

        HighVar_metabs <- sum(result_df$HighVar)
        if (HighVar_metabs>0) {
            message <-
                paste0(
                    HighVar_metabs,
                    " of variables have high variability (CV > 30) in the core_media control samples. Consider checking the pooled samples to decide whether to remove these metabolites or not."
                )
            log_info(message)
            message(message)
        }

        # Make histogram of CVs
        HistCV <- invisible(ggplot(cv_result_df, aes(CV)) +
                                geom_histogram(
                                    aes(y = after_stat(density)),
                                    color = "black",
                                    fill = "white"
                                )+
                                geom_vline(aes(xintercept = cutoff_cv),
                                color = "darkred", linetype = "dashed", linewidth = 1)+
                                geom_density(alpha = .2, fill = "#FF6666") +
                                labs(
                                    title = "CV for metabolites of control media samples (no cells)",
                                    x = "Coefficient of variation (CV)",
                                    y = "Frequency"
                                )+
                                theme_classic())

        HistCV_Sized <-
            plotGrob_Processing(
                input_plot = HistCV,
                plot_name = "CV for metabolites of control media samples (no cells)",
                plot_type = "Hist"
            )

        PlotList[["Histogram_coreMediaCV"]] <- HistCV_Sized

        # Make Violin of CVs
        Plot_cv_result_df <- cv_result_df %>%
        mutate(HighVar = ifelse(HighVar, "> CV 30", "< CV 30"))

        ViolinCV <- invisible(ggplot(Plot_cv_result_df, aes(y = CV, x = HighVar, label = row.names(cv_result_df)))+
        geom_violin(alpha = 0.5 , fill = "#FF6666")+
        geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5) +
        geom_text_repel(
            aes(label = ifelse(Plot_cv_result_df$CV > cutoff_cv,
            as.character(row.names(Plot_cv_result_df)), '')),
            hjust = 0,
            vjust = 0,
            box.padding = 0.5,
            # space between text and point
            point.padding = 0.5,
            # space around points
            max.overlaps = Inf
            ) + # allow for many labels
            labs(
                title = "CV for metabolites of control media samples (no cells)",
                x = "Metabolites",
                y = "Coefficient of variation (CV)"
            )+
            theme_classic())

        ViolinCV_Sized <-
            plotGrob_Processing(
                input_plot = ViolinCV,
                plot_name = "CV for metabolites of control media samples (no cells)",
                plot_type = "Violin"
            )
        PlotList[["core_Media_CV_Violin"]] <- ViolinCV_Sized

        ##########################################################################
        ## ------------------ Outlier testing
        if (dim(core_medias)[1]>=3) {
            Outlier_data <- core_medias
            Outlier_data <- Outlier_data %>% mutate_all(.funs = ~ FALSE)

            while (HighVar_metabs>0) {
                # remove the furthest value from the mean
                if (HighVar_metabs>1) {
                    max_var_pos <-  core_medias[,result_df$HighVar]  %>%
                    as.data.frame() %>%
                    mutate_all(.funs = ~ . - mean(., na.rm = TRUE)) %>%
                    summarise_all(.funs = ~ which.max(abs(.)))
                } else {
                    max_var_pos <-  core_medias[,result_df$HighVar]  %>%
                    as.data.frame() %>%
                    mutate_all(.funs = ~ . - mean(., na.rm = TRUE)) %>%
                    summarise_all(.funs = ~ which.max(abs(.)))
                    colnames(max_var_pos)<- colnames(core_medias)[result_df$HighVar]
                }

                # Remove rows based on positions
                for (i in seq_along(max_var_pos)) {
                    core_medias[max_var_pos[[i]],names(max_var_pos)[i]] <- NA
                    Outlier_data[max_var_pos[[i]],names(max_var_pos)[i]] <- TRUE
                }

                # ReCalculate coefficient of variation for each column in the filtered data
                result_df <-
                    apply(
                        core_medias,
                        2,
                        function(x) { sd(x, na.rm = TRUE)/  mean(x, na.rm = TRUE) } ) %>% t() %>% as.data.frame(
                    )
                result_df[1, is.na(result_df[1,])]<- 0
                rownames(result_df)[1] <- "CV"

                result_df <- result_df %>% t() %>% as.data.frame() %>% rowwise() %>%
                mutate(HighVar = CV > cutoff_cv) %>% as.data.frame()
                rownames(result_df)<- colnames(core_medias)

                HighVar_metabs <- sum(result_df$HighVar)
            }

            data_cont <- Outlier_data %>% t() %>% as.data.frame()

            # List to store results
            fisher_test_results <- list()
            large_contingency_table <-
                matrix(
                    0,
                    nrow = 2,
                    ncol = ncol(data_cont)
                )

            for (i in seq_along(colnames(data_cont))) {
                sample <- colnames(data_cont)[i]
                current_sample <- data_cont[, sample]

                contingency_table <- matrix(0, nrow = 2, ncol = 2)
                contingency_table[1, 1] <- sum(current_sample)
                contingency_table[2, 1] <- sum(!current_sample)
                contingency_table[1, 2] <- sum(rowSums(data_cont) - current_sample)
                contingency_table[2, 2] <- dim(data_cont %>% select(!all_of(sample)))[1]*dim(data_cont %>% select(!all_of(sample)))[2] -sum( rowSums(data_cont) - current_sample)

                # Fisher's exact test
                fisher_test_result <- fisher.test(contingency_table)
                fisher_test_results[[sample]] <- fisher_test_result

                # Calculate the sum of "TRUE" and "FALSE" for the current sample
                large_contingency_table[1, i] <- sum(current_sample)  # Sum of "TRUE"
                large_contingency_table[2, i] <- sum(!current_sample) # Sum of "FALSE"
            }

            # Convert the matrix into a data_contframe for better readability
            contingency_data_contframe <- as.data.frame(large_contingency_table)
            colnames(contingency_data_contframe) <- colnames(data_cont)
            rownames(contingency_data_contframe) <- c("HighVar", "Low_var")

            contingency_data_contframe <- contingency_data_contframe %>% mutate(total = rowSums(contingency_data_contframe))
            contingency_data_contframe <-
                rbind(
                    contingency_data_contframe,
                    total = colSums(contingency_data_contframe)
                )

            different_samples <- c()
            for (sample in colnames(data_cont)) {
                p_value <- fisher_test_results[[sample]]$p.value
                if (p_value < 0.05) {  # Adjust the significance level as needed
                    different_samples <- c(different_samples, sample)
                }
            }

            if (!is.null(different_samples)) {
                message <-
                    paste(
                        "The core_media samples ",
                        paste(different_samples, collapse = ", "),
                        " were found to be different from the rest. They will not be included in the sum of the core_media samples."
                    )
                log_trace("Warning: " , message, sep = "")
                warning(message)
            }
            # Filter the core_media samples
            core_medias <- core_medias %>% filter(!rownames(core_medias) %in% different_samples)
        } else {
            message <- paste0("Only >=2 blank samples available. Thus,we can not perform outlier testing for the blank samples.")
            log_trace(message)
            message(message)

        }
        core_media_df <- as.data.frame(data.frame("core_mediaMeans" = colMeans(core_medias, na.rm = TRUE)))
    }

    cv_result_df <- rownames_to_column(cv_result_df, "Metabolite")

    ##########################################################################
    # # ------------------------ Substract mean (media control) from samples
    message <- paste("core data are normalised by substracting mean (blank) from each sample and multiplying with the core_norm_factor")
    log_info(message)
    message(message)

    # # -- Check core_norm_factor
    if (("core_norm_factor" %in% names(metadata_info))) {
        core_norm_factor <-   metadata_sample %>% filter(!!as.name(metadata_info[["Conditions"]]) !=metadata_info[["core_media"]]) %>% select(metadata_info[["core_norm_factor"]]) %>% dplyr::pull()
        if (var(core_norm_factor) ==  0) {
            message <- paste("The growth rate or growth factor for normalising the core result, is the same for all samples")
            log_trace("Warning: " , message, sep = "")
            warning(message)
        }
    } else {
        core_norm_factor <- as.numeric(rep(1,dim(metadata_sample %>% filter(!!as.name(metadata_info[["Conditions"]]) !=metadata_info[["core_media"]]))[1]))
    }

    # Remove core_media samples from the data
    data_tic <- merge(metadata_sample, data_tic, by = "row.names") %>%
    filter(!!as.name(metadata_info[["Conditions"]]) !=metadata_info[["core_media"]]) %>%
    column_to_rownames("Row.names") %>%
    select(-seq_len(ncol(metadata_sample)))

    data_tic_coreNorm_Media <- as.data.frame(t( apply(t(data_tic),2, function(i) i-core_media_df$core_mediaMeans)))  #Subtract from each sample the core_media mean
    data_tic_coreNorm <- as.data.frame(apply(data_tic_coreNorm_Media, 2, function(i) i*core_norm_factor))

    # Remove core_media samples from the data
    # Input_SettingsFile <- Input_SettingsFile[Input_SettingsFile$Conditions !="core_media",]
    # Conditions <- Conditions[!Conditions =="core_media"]

    ##########################################################################
    # # ------------------------ Return Plots and data
    if (dim(core_medias)[1]>=3) {
        DF_list <-
            list(
                "CV_core_blank" = cv_result_df,
                "Contigency_table_core_blank" = contingency_data_contframe,
                "core_Norm" = data_tic_coreNorm
            )
    } else {
        DF_list <-
            list(
                "CV_core_blank" = cv_result_df,
                "core_Norm" = data_tic_coreNorm
            )
    }

    # Return
    Output_list <- list("DF" = DF_list,"Plot" = PlotList)
    invisible(return(Output_list))
}


#
# processing helper function: Outlier detection
#

#' outlier_detection
#'
#' @param data DF which contains unique sample identifiers as row names and metabolite
#'     numerical values in columns with metabolite identifiers as column names.
#'     Use NA for metabolites that were not detected and consider converting
#'     any zeros to NA unless they are true zeros.
#' @param metadata_sample DF which contains information about the samples, which will be combined
#'     with the input data based on the unique sample identifiers used as
#'     rownames.
#' @param metadata_info Named vector containing the information about the names of the
#'     experimental parameters. c(Conditions="ColumnName_Plot_SettingsFile",
#'     Biological_Replicates="ColumnName_Plot_SettingsFile", core_media =
#'     "Columnname_Input_SettingsFile"). Column "Conditions" with information
#'     about the sample conditions, Column "BiologicalReplicates" including
#'     numerical values and Column "Columnname_Input_SettingsFile" is used to
#'     specify the name of the media controls in the Conditions.
#' @param core \emph{Optional: } If TRUE, a consumption-release experiment has been
#'     performed. If not normalised yet, provide information about control
#'     media sample names via metadata_info "core_media" samples.
#'     \strong{Default = FALSE}
#' @param hotellins_confidence \emph{Optional: } Confidence level for Hotellin's T2 test.
#'     \strong{Default = 0.99}
#'
#' @return List with two elements: : DF (including output tables) and Plot
#'     (including all plots generated)
#'
#' @examples
#' data(intracell_raw)
#' Intra <- intracell_raw %>% tibble::column_to_rownames("Code")
#' Res <- outlier_detection(
#'     data =
#'         Intra[-c(49:58), -c(1:3)] %>%
#'         dplyr::mutate_all(
#'             ~ ifelse(grepl("^0*(\\.0*)?$", as.character(.)), NA, .)
#'         ),
#'     metadata_sample = Intra[-c(49:58), c(1:3)],
#'     metadata_info = c(
#'         Conditions = "Conditions",
#'         Biological_Replicates = "Biological_Replicates"
#'     )
#' )
#'
#' @importFrom dplyr relocate case_when mutate mutate_all select
#' @importFrom inflection uik
#' @importFrom factoextra fviz_screeplot
#' @importFrom ggplot2 ggplot theme_classic theme geom_vline annotate geom_line element_text
#' @importFrom ggplot2 geom_point geom_hline scale_y_continuous ggtitle scale_linetype_discrete
#' @importFrom qcc mqcc
#' @importFrom magrittr %>% %<>%
#' @importFrom tibble rownames_to_column
#' @importFrom hash values keys hash
#' @importFrom logger log_info log_trace
#' @importFrom grDevices dev.new dev.off
#' @importFrom stats relevel var
#' @noRd
outlier_detection <- function(
        data,
        metadata_sample,
        metadata_info,
        core = FALSE,
        hotellins_confidence = 0.99
) {

    # NSE vs. R CMD check workaround
    Samples <- `Group summary statistics` <- Outliers <- NULL
    # Message:
    message <- paste("Outlier detection: Identification of outlier samples is performed using Hotellin's T2 test to define sample outliers in a mathematical way (Confidence = 0.99 ~ p.val < 0.01) (REF: Hotelling, H. (1931), Annals of Mathematical Statistics. 2 (3), 360-378, doi:https://doi.org/10.1214/aoms/1177732979). ",
    "hotellins_confidence value selected: ", hotellins_confidence, sep = "")
    log_info(message)
    message(message)

    # Load the data:
    data_norm <- data %>%
    mutate_all(~ replace(., is.nan(.), 0))
    data_norm[is.na(data_norm)] <- 0 #replace NA with 0

    if (core) {
        Conditions <- metadata_sample[[metadata_info[["Conditions"]]]][!metadata_sample[[metadata_info[["Conditions"]]]] == metadata_info[["core_media"]]]
    } else {
        Conditions <- metadata_sample[[metadata_info[["Conditions"]]]]
    }


    # Prepare the lists to store the results:
    Outlier_filtering_loop <- 10 #Here we do 10 rounds of hotelling filtering
    sample_outliers <- list()
    scree_plot_list <- list()
    outlier_plot_list <- list()
    metabolite_zero_var_total_list <- list()
    zero_var_metab_warning <- FALSE


    # # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## #
    # # --------- Perform Outlier testing:
    for (loop in seq_len(Outlier_filtering_loop)) {
        # # --- Zero variance metabolites
        metabolite_var <- as.data.frame(apply(data_norm, 2, var) %>% t()) # calculate each metabolites variance
        metabolite_zero_var_list <- list(colnames(metabolite_var)[which(metabolite_var[1,] == 0)]) # takes the names of metabolites with zero variance and puts them in list

        if (sum(metabolite_var[1,] == 0) == 0) {
            metabolite_zero_var_total_list[loop] <- 0
        }else if (sum(metabolite_var[1,] == 0)>0) {
            metabolite_zero_var_total_list[loop] <- metabolite_zero_var_list
            zero_var_metab_warning <- TRUE # This is used later to print and save the zero variance metabolites if any are found.
        }

        for (metab in metabolite_zero_var_list) {  # Remove the metabolites with zero variance from the data to do PCA
            data_norm <- data_norm %>% select(-all_of(metab))
        }

        # # ---  PCA
        PCA.res <- prcomp(data_norm, center = TRUE, scale. = TRUE)
        outlier_PCA_data <- data_norm
        outlier_PCA_data[[metadata_info[["Conditions"]]]] <- Conditions

        dev.new()
        pca_outlier <- invisible(
            viz_pca(data = data_norm,
                    metadata_info = c(color = metadata_info[["Conditions"]]),
                    metadata_sample = outlier_PCA_data,
                    plot_name = paste("PCA outlier test filtering round ",loop),
                    save_plot = NULL)
        )

        if (loop == 1) {
            pca_outlierloop1 <- pca_outlier[["Plot_Sized"]][[1]]
        }
        outlier_plot_list[[paste("PCA_round",loop,sep = "")]] <- pca_outlier[["Plot_Sized"]][[1]]
        dev.off()

        # # --- Scree plot
        inflect_df <- as.data.frame(c(seq_along(PCA.res$sdev))) # get Scree plot values for inflection point calculation
        colnames(inflect_df) <- "x"
        inflect_df$y <- summary(PCA.res)$importance[2,]
        inflect_df$Cumulative <- summary(PCA.res)$importance[3,]
        screeplot_cumul <-
            format(
                round(inflect_df$Cumulative[seq_len(20)]*100, 1),
                nsmall = 1
            ) #make cumulative variation labels for plot
        knee <-
            uik(
                inflect_df$x,
                inflect_df$y
            ) # Calculate the knee and select optimal number of components
        npcs <- knee -1 #Note: we subtract 1 components from the knee cause the root of the knee is the PC that does not add something. npcs = 30

        # Make a scree plot with the selected component cut-off for HotellingT2 test
        screeplot <-
            fviz_screeplot(
                PCA.res,
                main = paste("PCA Explained variance plot filtering round ",loop, sep = ""
            ),
                                    addlabels = TRUE,
                                    ncp = 20,
                                    geom = c("bar", "line"),
                                    barfill = "grey",
                                    barcolor = "grey",
                                    linecolor = "black",linetype = 1) +
            theme_classic()+
            geom_vline(xintercept = npcs+0.5, linetype = 2, color = "red") +
            annotate(
                "text",
                x = c(seq_len(20)),
                y = -0.8,
                label = screeplot_cumul,
                col = "black",
                size = 1.75
            )

        # screeplot_Sized <- plotGrob_Processing(input_plot = screeplot,plot_name= paste("PCA Explained variance plot filtering round ",loop, sep = ""), plot_type= "Scree")

        if (loop == 1) {
            scree_outlierloop1 <- screeplot
        }
        dev.new()

        outlier_plot_list[[paste("ScreePlot_round",loop,sep = "")]] <- screeplot # save plot
        dev.off()

        # # --- HotellingT2 test for outliers
        data_hot <- as.matrix(PCA.res$x[,seq_len(npcs)])
        hotelling_qcc <-
            mqcc(
                data_hot,
                type = "T2.single",
                labels = rownames(data_hot),
                confidence.level = hotellins_confidence,
                title = paste("Outlier filtering via HotellingT2 test filtering round ",loop,", with ",hotellins_confidence, "% Confidence",  sep = ""),
                plot = FALSE
            )
        HotellingT2plot_data <- as.data.frame(hotelling_qcc$statistics)
        HotellingT2plot_data <-
            rownames_to_column(
                HotellingT2plot_data,
                "Samples"
            )
        colnames(
            HotellingT2plot_data) <- c("Samples",
            "Group summary statistics"
        )
        outlier <- HotellingT2plot_data %>% filter(HotellingT2plot_data$`Group summary statistics`>hotelling_qcc$limits[2])
        limits <- as.data.frame(hotelling_qcc$limits)
        legend <- colnames(HotellingT2plot_data[2])
        LegendTitle <- "Limits"

        HotellingT2plot <-
            ggplot(
                HotellingT2plot_data,
                aes(x = Samples, y = `Group summary statistics`, group = 1, fill = )
            )
        HotellingT2plot <- HotellingT2plot +
            geom_point(
                aes(x = Samples,y = `Group summary statistics`),
                color = 'blue',
                size = 2
            ) +
            geom_point(
                data = outlier,
                aes(x = Samples,y = `Group summary statistics`),
                color = 'red',
                size = 3
            ) +
            geom_line(linetype = 2)

        # draw the horizontal lines corresponding to the LCL,UCL
        HotellingT2plot <- HotellingT2plot +
            geom_hline(
                aes(yintercept = limits[,1]),
                color = "black",
                data = limits,
                show.legend = FALSE
            ) +
            geom_hline(
                aes(yintercept = limits[,2], linetype = "UCL"),
                color = "red",
                data = limits,
                show.legend = TRUE
            ) +
            scale_y_continuous(breaks = sort(c(ggplot_build(HotellingT2plot)$layout$panel_ranges[[1]]$y.major_source, c(limits[,1],limits[,2]))))

        HotellingT2plot <- HotellingT2plot +
            theme_classic()+
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
            ggtitle(paste("Hotelling ", hotelling_qcc$type ," test filtering round ",loop,", with ", 100 * hotelling_qcc$confidence.level,"% Confidence"))+
            scale_linetype_discrete(name = LegendTitle,)+
            theme(plot.title = element_text(size = 13))+ # , face = "bold")) +
            theme(axis.text = element_text(size = 7))
        # HotellingT2plot_Sized <- plotGrob_Processing(input_plot = HotellingT2plot,plot_name= paste("Hotelling ", hotelling_qcc$type ," test filtering round ",loop,", with ", 100 * hotelling_qcc$confidence.level,"% Confidence"), plot_type= "Hotellings")

        if (loop == 1) {
            hotel_outlierloop1 <- HotellingT2plot
        }
        dev.new()
        # plot(HotellingT2plot)
        outlier_plot_list[[paste("HotellingsPlot_round",loop,sep = "")]] <- HotellingT2plot
        dev.off()

        a <- loop
        if (core) {
            a <- paste0(a,"_core")
        }

        if (length(hotelling_qcc[["violations"]][["beyond.limits"]]) == 0) { # loop for outliers until no outlier is detected
            data_norm <- data_norm
            break
        }else if (length(hotelling_qcc[["violations"]][["beyond.limits"]]) == 1) {
            data_norm <- data_norm[-hotelling_qcc[["violations"]][["beyond.limits"]],]# filter the selected outliers from the data
            Conditions <- Conditions[-hotelling_qcc[["violations"]][["beyond.limits"]]]

            # Change the names of outliers in mqcc . Instead of saving the order number it saves the name
            hotelling_qcc[["violations"]][["beyond.limits"]][1] <-   rownames(data_hot)[hotelling_qcc[["violations"]][["beyond.limits"]][1]]
            sample_outliers[loop] <- list(hotelling_qcc[["violations"]][["beyond.limits"]])
        } else {
            data_norm <- data_norm[-hotelling_qcc[["violations"]][["beyond.limits"]],]
            Conditions <- Conditions[-hotelling_qcc[["violations"]][["beyond.limits"]]]

            # Change the names of outliers in mqcc . Instead of saving the order number it saves the name
            sm_out <- c() # list of outliers samples
            for (i in seq_along(hotelling_qcc[["violations"]][["beyond.limits"]])) {
                sm_out <-
                    append(
                        sm_out,
                        rownames(data_hot)[hotelling_qcc[["violations"]][["beyond.limits"]][i]]
                    )
            }
            sample_outliers[loop] <- list(sm_out )
        }
    }

    # # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## #
    # # -- Print Outlier detection results about samples and metabolites
    if (length(sample_outliers) > 0) {   # Print outlier samples
        message <- paste("There are possible outlier samples in the data")
        log_info(message)
        message(message) #This was a warning
        for (i in seq_along(sample_outliers)  ) {
            message <-
                paste(
                    "Filtering round ",
                    i,
                    " Outlier Samples: ",
                    paste( head(sample_outliers[[i]]) ," ")
                )
            log_info(message)
            message(message)
        }
    } else {
        message <- paste("No sample outliers were found")
        log_info(message)
        message(message)
    }

    # # --  Print Zero variance metabolites
    zero_var_metab_export_df <- data.frame(1,2)
    names(zero_var_metab_export_df) <- c("Filtering round","Metabolite")

    if (zero_var_metab_warning) {
        message <- paste("Metabolites with zero variance have been identified in the data. As scaling in PCA cannot be applied when features have zero variace, these metabolites are not taken into account for the outlier detection and the PCA plots.")
        log_trace("Warning: " , message, sep = "")
        warning(message)
    }

    count <- 1
    for (i in seq_along(metabolite_zero_var_total_list)) {
        if (metabolite_zero_var_total_list[[i]] != 0) {
            message <-
                paste(
                    "Filtering round ",
                    i,
                    ". Zero variance metabolites identified: ",
                    paste( metabolite_zero_var_total_list[[i]] ," ")
                )
            log_info(message)
            message(message)

            zero_var_metab_export_df[count,"Filtering round"] <- paste(i)
            zero_var_metab_export_df[count,"Metabolite"] <- paste(metabolite_zero_var_total_list[[i]])
            count <- count +1
        }
    }

    # # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## #
    # # ---- 1. Make Output DF
    total_outliers <- hash() # make a dictionary
    if (length(sample_outliers) > 0) { # Create columns with outliers to merge to output dataframe
        for (i in seq_along(sample_outliers)  ) {
            total_outliers[[paste("Outlier_filtering_round_",i, sep = "")]] <- sample_outliers[i]
        }
    }

    data_norm_filtered_full <- as.data.frame(replace(data, data == 0, NA))

    if (length(total_outliers) > 0) {  # add outlier information to the full output dataframe
        data_norm_filtered_full$Outliers <- "no"
        for (i in seq_along(total_outliers)) {
            for (k in seq_along( values(total_outliers)[i] ) ) {
                data_norm_filtered_full[as.character(values(total_outliers)[[i]]) , "Outliers"] <- keys(total_outliers)[i]
            }
        }
    } else {
        data_norm_filtered_full$Outliers <- "no"
    }

    data_norm_filtered_full <-
        relocate(Outliers) %>%  # Put Outlier columns in the front
        merge(
            metadata_sample,
            data_norm_filtered_full,
            by = 0  # add the design in the output df (merge by rownames/sample names
        )
    rownames(data_norm_filtered_full) <- data_norm_filtered_full$Row.names
    data_norm_filtered_full$Row.names <- c()

    # # -- 2.  Quality Control (QC) PCA
    Metadata_Sample <- data_norm_filtered_full %>%
        mutate(
            Outliers = case_when(Outliers == "no" ~ 'no',
            Outliers == "Outlier_filtering_round_1" ~ ' Outlier_filtering_round = 1',
            Outliers == "Outlier_filtering_round_2" ~ ' Outlier_filtering_round = 2',
            Outliers == "Outlier_filtering_round_3" ~ ' Outlier_filtering_round = 3',
            Outliers == "Outlier_filtering_round_4" ~ ' Outlier_filtering_round = 4',
            TRUE ~ 'Outlier_filtering_round = or > 5')
        )
    Metadata_Sample$Outliers <- relevel(as.factor(Metadata_Sample$Outliers), ref = "no")

    # 1. Shape Outliers
    if (length(sample_outliers)>0) {
        dev.new()
        pca_QC <- invisible(viz_pca(data = as.data.frame(data) %>% dplyr::select(-zero_var_metab_export_df$Metabolite),
                                    metadata_info = c(color = metadata_info[["Conditions"]], shape = "Outliers"),
                                    metadata_sample = Metadata_Sample ,
                                    plot_name = "Quality Control PCA Condition clustering and outlier check",
                                    save_plot = NULL))
        dev.off()
        outlier_plot_list[["QC_PCA_and_Outliers"]] <- pca_QC[["Plot_Sized"]][[1]]
    }

    # 2. Shape Biological replicates
    if ("Biological_Replicates" %in% names(metadata_info)) {
        dev.new()
        pca_QC_repl <- invisible(viz_pca(data = as.data.frame(data) %>% dplyr::select(-zero_var_metab_export_df$Metabolite),
        metadata_info = c(color = metadata_info[["Conditions"]], shape = metadata_info[["Biological_Replicates"]]),
        metadata_sample = Metadata_Sample,
        plot_name = "Quality Control PCA replicate spread check",
        save_plot = NULL))
        dev.off()

        outlier_plot_list[["QC_PCA_Replicates"]] <- pca_QC_repl[["Plot_Sized"]][[1]]
    }


    # # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## #
    # # --- Save and Return plots and DFs
    DF_list <-
        list(
            "Zero_variance_metabolites_core" = zero_var_metab_export_df,
            "data_outliers" = data_norm_filtered_full
        )

    # Return
    Output_list <- list("DF" = DF_list,"Plot" = outlier_plot_list)
    invisible(return(Output_list))
}
