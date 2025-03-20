## ---------------------------
##
## Script name: PreProcessing
##
## Purpose of script: Metabolomics (raw ion counts) pre-processing, normalization, outlier detection and QC plots
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



###################################################
### ### ### Metabolomics pre-processing ### ### ###
###################################################

#' Modularised Normalization: 80%-filtering rule, total-ion count normalization, missing value imputation and Outlier Detection: HotellingT2.
#'
#' @param InputData DF which contains unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected.
#' @param SettingsFile_Sample DF which contains information about the samples, which will be combined with the input data based on the unique sample identifiers used as rownames.
#' @param SettingsInfo  Named vector containing the information about the names of the experimental parameters. c(Conditions="ColumnName_Plot_SettingsFile", Biological_Replicates="ColumnName_Plot_SettingsFile"). Column "Conditions" with information about the sample conditions (e.g. "N" and "T" or "Normal" and "Tumor"), can be used for feature filtering and colour coding in the PCA. Column "BiologicalReplicates" including numerical values. For CoRe = TRUE a CoRe_norm_factor = "Columnname_Input_SettingsFile" and CoRe_media = "Columnname_Input_SettingsFile", have to also be added. Column CoRe_norm_factor is used for normalization and CoRe_media is used to specify the name of the media controls in the Conditions.
#' @param FeatureFilt \emph{Optional: }If NULL, no feature filtering is performed. If set to "Standard" then it applies the 80%-filtering rule (Bijlsma S. et al., 2006) on the metabolite features on the whole dataset. If is set to "Modified",filtering is done based on the different conditions, thus a column named "Conditions" must be provided in the Input_SettingsFile input file including the individual conditions you want to apply the filtering to (Yang, J et al., 2015). \strong{Default = "Standard"}
#' @param FeatureFilt_Value \emph{Optional: } Percentage of feature filtering. \strong{Default = 0.8}
#' @param TIC \emph{Optional: } If TRUE, Total Ion Count normalization is performed. \strong{Default = TRUE}
#' @param MVI \emph{Optional: } If TRUE, Missing Value Imputation (MVI) based on half minimum is performed \strong{Default = TRUE}
#' @param MVI_Percentage \emph{Optional: } Percentage 0-100 of imputed value based on the minimum value. \strong{Default = 50}
#' @param HotellinsConfidence \emph{Optional: } Defines the Confidence of Outlier identification in HotellingT2 test. Must be numeric.\strong{Default = 0.99}
#' @param CoRe \emph{Optional: } If TRUE, a consumption-release experiment has been performed and the CoRe value will be calculated. Please consider providing a Normalisation factor column called "CoRe_norm_factor" in your "Input_SettingsFile" DF, where the column "Conditions" matches. The normalisation factor must be a numerical value obtained from growth rate that has been obtained from a growth curve or growth factor that was obtained by the ratio of cell count/protein quantification at the start point to cell count/protein quantification at the end point.. Additionally control media samples have to be available in the "Input" DF and defined as "CoRe_media" samples in the "Conditions" column in the "Input_SettingsFile" DF. \strong{Default = FALSE}
#' @param SaveAs_Plot \emph{Optional: } Select the file type of output plots. Options are svg, png, pdf. If set to NULL, plots are not saved. \strong{Default = svg}
#' @param SaveAs_Table \emph{Optional: } Select the file type of output table. Options are "csv", "xlsx", "txt". If set to NULL, plots are not saved. \strong{Default = "csv"}
#' @param PrintPlot  \emph{Optional: } If TRUE prints an overview of resulting plots. \strong{Default = TRUE}
#' @param FolderPath \emph{Optional:} Path to the folder the results should be saved at. \strong{default: NULL}
#'
#' @return List with two elements: DF (including all output tables generated) and Plot (including all plots generated)
#'
#' @examples
#' ## load the data and mapping Info
#' Intra <- ToyData("IntraCells_Raw")
#' MappingInfo <- ToyData(Data = "Cells_MetaData")
#' Media <- ToyData("CultureMedia_Raw")
#' 
#' ## create SummarizedExperiment objects
#' ## se_intra
#' rD <- MappingInfo
#' cD <- Intra[-c(49:58), c(1:3)]
#' a <- t(Intra[-c(49:58), -c(1:3)])
#' 
#' ## obtain overlapping metabolites
#' metabolites <- intersect(rownames(a), rownames(rD))
#' rD <- rD[metabolites, ]
#' a <- a[metabolites, ]
#' se_intra <- SummarizedExperiment::SummarizedExperiment(assays = a, rowData = rD, colData = cD)
#' 
#' ## se_media
#' rD <- MappingInfo
#' cD <- Media[, 1:3]
#' a <- t(Media[, 4:ncol(Media)])
#' 
#' ## obtain overlapping metabolites
#' metabolites <- intersect(rownames(a), rownames(rD))
#' rD <- rD[metabolites, ]
#' a <- a[metabolites, ]
#' se_media <- SummarizedExperiment::SummarizedExperiment(assays = a, rowData = rD, colData = cD)
#' 
#' ## apply the functions
#' ResI <- PreProcessing(
#'     se = se_intra,
#'     SettingsInfo = c(Conditions = "Conditions", 
#'         Biological_Replicates = "Biological_Replicates"))
#'  
#' ResM <- PreProcessing(
#'     se = se_media,
#'     SettingsInfo = c(Conditions = "Conditions", 
#'         Biological_Replicates = "Biological_Replicates", 
#'         CoRe_norm_factor = "GrowthFactor", CoRe_media = "blank"),
#'     CoRe = TRUE)
#'
#' @keywords 80  percent filtering rule, Missing Value Imputation, Total Ion Count normalization, PCA, HotellingT2, multivariate quality control charts
#'
#' @importFrom dplyr mutate_all
#' @importFrom magrittr %>% %<>%
#' @importFrom tibble rownames_to_column column_to_rownames
#'
#' @export
#'
PreProcessing <- function(
    se,
    #InputData,
    #SettingsFile_Sample,
    SettingsInfo,
    FeatureFilt = "Modified",
    FeatureFilt_Value = 0.8,
    TIC = TRUE,
    MVI = TRUE,
    MVI_Percentage = 50,
    HotellinsConfidence = 0.99,
    CoRe = FALSE,
    SaveAs_Plot = "svg",
    SaveAs_Table = "csv",
    PrintPlot = TRUE,
    FolderPath = NULL) {
  
    ## ------------ Create log file ----------- ##
    MetaProViz_Init()

    ## obtain InputData and SettingsFile_Sample
    InputData <- assay(se) |>
        t() |>
        as.data.frame() 
    SettingsFile_Sample <- colData(se) |>
        as.data.frame()
    
    ## ------------------ Check Input ------------------- ##
    ## HelperFunction `CheckInput`
    CheckInput(
        se = se,
        ##InputData = InputData, SettingsFile_Sample = SettingsFile_Sample,
        #SettingsFile_Metab = NULL, 
        SettingsInfo = SettingsInfo,
        SaveAs_Plot = SaveAs_Plot, SaveAs_Table = SaveAs_Table,
        CoRe = CoRe, PrintPlot = PrintPlot)

    ## HelperFunction `CheckInput` Specific
    CheckInput_PreProcessing(
        se = se,
        #SettingsFile_Sample = SettingsFile_Sample,
        SettingsInfo = SettingsInfo, CoRe = CoRe, FeatureFilt = FeatureFilt,
        FeatureFilt_Value = FeatureFilt_Value, TIC = TIC, MVI = MVI,
        MVI_Percentage = MVI_Percentage, 
        HotellinsConfidence = HotellinsConfidence)

    ## ------------------  Create output folders  and path ------------------- ##
    if (!is.null(SaveAs_Plot) |!is.null(SaveAs_Table)) {
        Folder <- SavePath(FolderName = "Processing", FolderPath = FolderPath)

        SubFolder_P <- file.path(Folder, "PreProcessing")
        if (!dir.exists(SubFolder_P)) {
            dir.create(SubFolder_P)
        }
    }

    ## ------------------ Prepare the data ------------------- ##
    ## InputData files:
    ## make sure all 0 are changed to NAs
    #InputData <- as.data.frame(InputData) %>%
    #    dplyr::mutate_all(~ ifelse(grepl("^0*(\\.0*)?$", as.character(.)), NA, .)) ## EDIT: why not:
    assay(se)[assay(se) == 0] <- NA

    #InputData <- as.data.frame(
    #    dplyr::mutate_all(as.data.frame(InputData), function(x) 
    #        as.numeric(as.character(x))))

    ############################################################################
    ## ------------------ 1. Feature filtering ------------------- ##
    if (!is.null(FeatureFilt)) {
        l_filtered <- FeatureFiltering(
            se = se,## InputData = InputData,
            FeatureFilt = FeatureFilt,
            FeatureFilt_Value = FeatureFilt_Value,
            #SettingsFile_Sample = SettingsFile_Sample,
            SettingsInfo = SettingsInfo,
            CoRe = CoRe)

        se_Filt <- l_filtered[["data"]][["se"]]
    } else {
        se_Filt <- se
    }

    ## ------------------ 2. Missing value Imputation ------------------- ##
    if (MVI) {
        l_imputed <- MVImputation(se = se_Filt, ##InputData = InputData_Filt,
            ##SettingsFile_Sample = SettingsFile_Sample,
            SettingsInfo = SettingsInfo,
            CoRe = CoRe,
            MVI_Percentage = MVI_Percentage)
        se_MVI <- l_imputed[["data"]][["se"]]
    } else {
        se_MVI <- se_Filt
    }

    ## ----------------  3. Total Ion Current Normalization ----------------- ##
    if (TIC) {
        
        ## perform TIC normalization
        l_tic <- TICNorm(se = se_MVI,
            SettingsInfo = SettingsInfo,
            TIC = TIC)
        se_tic <- l_tic[["data"]][["se"]]

        ## add plots to PlotList
        PlotList <- list()
        PlotList[["RLAPlot"]] <- l_tic[["plot"]][["beforeTicNormalization"]]
        PlotList[["RLAPlot_TICnorm"]] <- l_tic[["plot"]][["afterTicNormalization"]]
        PlotList[["RLAPlot_BeforeAfter_TICnorm"]] <- l_tic[["plot"]][["combined"]]
        
    } else {
        ##se_tic <- se_MVI ## EDIT: could also use the SE object returned from TICNorm?

        ## perform TIC normalization (TIC = FALSE)
        l_notic <- TICNorm(se = se_tic, ## EDIT: could it have the same name l_tic?
            SettingsInfo = SettingsInfo,
            TIC = TIC)
        se_tic <- l_tic[["data"]][["se"]]
        
        ## add plots to PlotList
        PlotList <- list()
        PlotList[["RLAPlot"]] <- l_notic[["plot"]][["beforeTicNormalization"]]
    }

    ## ------------- 4. CoRe media QC (blank) and normalization ------------- ##
    if (CoRe) {
        l_CoReNorm <- CoReNorm(se = se_tic, ##InputData = TICRes,
            ##SettingsFile_Sample = SettingsFile_Sample,
            SettingsInfo = SettingsInfo)
    
        se_tic <- l_CoReNorm[["data"]][["se"]]
    }

    ## ------------------ Final Output:
    
    ############################################################################
    ## ------------------ Sample outlier identification ------------------- ##
    l_outlier <-  OutlierDetection(se = se_tic, ##InputData = data_norm,
        ##SettingsFile_Sample = SettingsFile_Sample,
        SettingsInfo = SettingsInfo,
        CoRe = CoRe,
        HotellinsConfidence = HotellinsConfidence)

    ## continue from here ...
    ############################################################################
    ## ------------------ Return ------------------- ##
    ## ---- DFs
    if (!is.null(FeatureFilt)) {
        
        ## add metabolites that where removed as part of the feature filtering
        if (length(l_filtered[["RemovedMetabolites"]]) == 0) {
            
            l <- list(
                "se_raw"= se,
                "Filtered_metabolites"= as.data.frame(
                    list(FeatureFiltering = c(FeatureFilt), ## EDIT: are the c() needed?
                        FeatureFilt_Value = c(FeatureFilt_Value),
                        RemovedMetabolites = c("None"))), ## EDIT: for simplicity why not only return here l_filtered[["RemovedMetabolites"]]? / have only one return not dependong on length(l_filtered[["RemovedMetabolites"]])?
                "se_preprocessed" = l_outlier[["data"]][["se"]])
        } else {
            l <- list(
                "se_raw"= se,
                "Filtered_metabolites"= as.data.frame(
                    list(
                        FeatureFiltering = rep(FeatureFilt, length(l_filtered[["RemovedMetabolites"]])),
                        FeatureFilt_Value = rep(FeatureFilt_Value, length(l_filtered[["RemovedMetabolites"]])),
                        RemovedMetabolites = l_filtered[["RemovedMetabolites"]])),
                "se_preprocessed" = l_outlier[["data"]][["se"]]) 
        }
    } else {
        l <- list(
            "se_raw"= se, 
            "se_processed" = l_outlier[["data"]][["se"]]) ## EDIT: it seems to me that only Filtered_metabolites is different here, better to put this in the if/else and assemble everything outside of it
    }

    if (CoRe) {
        if (is.null(l_CoReNorm[["Contigency_table_CoRe_blank"]])) {
            l_CoRe <- list(
                "CV_CoRe_blank"= l_CoReNorm[["CV_CoRe_blank"]])
        } else {
            l_CoRe <- list(
                "CV_CoRe_blank" = l_CoReNorm[["CV_CoRe_blank"]],
                "Variation_ContigencyTable_CoRe_blank" = l_CoReNorm[["Contigency_table_CoRe_blank"]])
        }
        l <- c(l, l_CoRe)
    }

    ## ---- Plots
    if (TIC) {
        l_plot <- c(l_tic[["plot"]], l_outlier[["plot"]])
    } else {
        l_plot <- c(l_notic[["plot"]], l_outlier[["plot"]])
    }

    if (CoRe) {
        l_plot <- c(l_plot , l_CoReNorm[["plot"]])
    }

    ## save Plots and DFs
    ## as row names are not saved we need to make row.names to column for 
    ## the DFs that needs this:
    ##DFList[["InputData_RawData"]] <- DFList[["InputData_RawData"]] %>% 
    ##    tibble::rownames_to_column("Code")
    ##DFList[["Preprocessing_output"]] <- DFList[["Preprocessing_output"]] %>% 
    ##    tibble::rownames_to_column("Code") ## EDIT: not needed since we return SummarizedExperiment object

  suppressMessages(suppressWarnings(
        SaveRes(data = l,
            plot = l_plot,
            SaveAs_Table =SaveAs_Table,
            SaveAs_Plot = SaveAs_Plot,
            FolderPath = SubFolder_P,
            FileName = "PreProcessing",
            CoRe = CoRe,
            PrintPlot = PrintPlot)))

  ## return
  invisible(list("data" = l, "plot" = l_plot))
}


############################################################
### ### ### Merge analytical replicates function ### ### ###
############################################################

#' Merges the analytical replicates of an experiment
#'
#' @param InputData DF which contains unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected.
#'#@param SettingsFile_Sample DF which contains information about the samples Column "Conditions", "Biological_replicates" and "Analytical_Replicates has to exist.
#' @param SettingsInfo  \emph{Optional: } Named vector including the Conditions and Replicates information: c(Conditions="ColumnNameConditions", Biological_Replicates="ColumnName_SettingsFile_Sample", Analytical_Replicates="ColumnName_SettingsFile_Sample").\strong{Default = NULL}
#' @param SaveAs_Table \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt", ot NULL \strong{default: "csv"}
#' @param FolderPath \emph{Optional:} Path to the folder the results should be saved at. \strong{default: NULL}
#'
#' @return DF with the merged analytical replicates
#'
#' @examples
#' ## load the data
#' Intra <- ToyData("IntraCells_Raw")
#' 
#' ## create SummarizedExperiment
#' a <- t(Intra[-c(49:58), -c(1:3)])
#' rD <- DataFrame(feature = rownames(a))
#' cD <- Intra[-c(49:58), c(1:3)]
#' se <- SummarizedExperiment(assay = a, rowData = rD, colData = cD)
#' 
#' ## apply the function
#' ReplicateSum(se = se,
#'     SettingsInfo = c(Conditions = "Conditions", 
#'         Biological_Replicates = "Biological_Replicates", 
#'         Analytical_Replicates = "Analytical_Replicates"))
#'
#' @keywords Analytical Replicate Merge
#'
#' @importFrom dplyr mutate_all summarise_all select rename ungroup group_by
#' @importFrom magrittr %>% %<>%
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom rlang !! :=
#' @importFrom tidyr unite
#'
#' @export
#'
ReplicateSum <- function(se, ##InputData, ## EDIT: the name of this function is not informative? summarizeAnalyticalReplicates?
    ##SettingsFile_Sample,
    SettingsInfo = c(Conditions = "Conditions", 
        Biological_Replicates = "Biological_Replicates", 
        Analytical_Replicates = "Analytical_Replicates"),
    SaveAs_Table = "csv", ## EDIT: list here the options
    FolderPath = NULL) {
    
    ## ------------ Create log file ----------- ##
    MetaProViz_Init()

    ## ------------------ Check Input ------------------- ##
    ## HelperFunction `CheckInput`
    CheckInput(se, ##InputData = InputData,
        ##SettingsFile_Sample = SettingsFile_Sample,
        SettingsFile_Metab = NULL,
        SettingsInfo = SettingsInfo,
        SaveAs_Plot = NULL,
        SaveAs_Table = SaveAs_Table,
        CoRe = FALSE,
        PrintPlot = FALSE)
    
    ## create object that will simplify the calculations
    cD <- colData(se) |>
        as.data.frame()

    ## `CheckInput` Specific
    if (SettingsInfo[["Conditions"]] %in% colnames(cD)) {
        ## Conditions <- InputData[[SettingsInfo[["Conditions"]] ]] ## EDIT: simplify
    } else {
        stop("Column `Conditions` is required.")
    }
    if (SettingsInfo[["Biological_Replicates"]] %in% colnames(cD)) {
        ## Biological_Replicates <- InputData[[SettingsInfo[["Biological_Replicates"]]]]
    } else {
        stop("Column `Biological_Replicates` is required.")
    }
    if (SettingsInfo[["Analytical_Replicates"]] %in% colnames(cD)) {
        #Analytical_Replicates <- InputData[[SettingsInfo[["Analytical_Replicates"]]]]
    } else {
        stop("Column `Analytical_Replicates` is required.")
    } ## EDIT: why have the if/else here, if for "if" nothing is done? 

    ## ------------ Create Results output folder ----------- ##
    if (!is.null(SaveAs_Table)) {
        Folder <- SavePath(FolderName = "Processing", FolderPath = FolderPath)
        SubFolder <- file.path(Folder, "ReplicateSum")
        if (!dir.exists(SubFolder)) {
            dir.create(SubFolder)
        }
    }

    ## ------------  Load data and process  ----------- ##
    ##se_merged <- merge(
    ###    x =  dplyr::select(as.data.frame(colData(se)), !!SettingsInfo[["Conditions"]], 
    ##        !!SettingsInfo[["Biological_Replicates"]], 
    ##        !!SettingsInfo[["Analytical_Replicates"]]),
    ##    y = t(assay(se)),
    ##    by = "row.names") %>%
    ##    tibble::column_to_rownames("Row.names") %>%
    ##    dplyr::rename("Conditions" = SettingsInfo[["Conditions"]],
    ##        "Biological_Replicates" = SettingsInfo[["Biological_Replicates"]],
    ##        "Analytical_Replicates" = SettingsInfo[["Analytical_Replicates"]]) ## EDIT: not needed with SE

    ## Make the replicate Sums
    assay_summed <- t(assay(se)) |> 
        as.data.frame() |>
        dplyr::group_by(
            Biological_Replicates = cD[[SettingsInfo[["Biological_Replicates"]]]], 
            Conditions = cD[[SettingsInfo[["Conditions"]]]]) %>%
        dplyr::summarise_all("mean") %>% 
        #dplyr::select(-Analytical_Replicates) %>% 
        as.data.frame()

    ## make a number of merged replicates column
    n_replicates <- cD |>
        dplyr::group_by(Biological_Replicates, Conditions) %>%
        dplyr::summarise_all("max") %>%
        dplyr::ungroup() %>%
        dplyr::select(Analytical_Replicates, Biological_Replicates, Conditions) %>%
        dplyr::rename("n_AnalyticalReplicates_Summed "= "Analytical_Replicates") |>
        dplyr::mutate(sample_ID = paste(n_replicates[["Conditions"]], 
            n_replicates[["Biological_Replicates"]], sep = "_")) |>
        tibble::column_to_rownames(var = "sample_ID")

    ## create SummarizedExperiment
    a_colnames <- paste(assay_summed[["Conditions"]], assay_summed[["Biological_Replicates"]], sep = "_")
    a <- assay_summed[, rownames(se)] |>
        t()
    colnames(a) <- a_colnames
    
    ## make sure that a has same column order than row order of n_replicates and
    ## same row order than row order of rowData(se)
    a <- a[, rownames(n_replicates)]
    a <- a[rownames(rowData(se)), ]
    se <- SummarizedExperiment(assay = a, rowData = rowData(se), colData = n_replicates)
    
    ##assay_summed <- merge(n_replicates, assay_summed, 
    ##        by = c("Conditions", "Biological_Replicates")) %>%
    ##    ## create a uniqueID
    ##    tidyr::unite(UniqueID, c("Conditions", "Biological_Replicates"), 
    ##        sep = "_", remove = FALSE) %>% 
    ##    ## set UniqueID to rownames
    ##    tibble::column_to_rownames("UniqueID")

    ##--------------- return ------------------##
    l <- list("se" = se)
    SaveRes(data = list("se" = se),
        plot = NULL,
        SaveAs_Table = SaveAs_Table,
        SaveAs_Plot = NULL,
        FolderPath = SubFolder,
        FileName = "Sum_AnalyticalReplicates",
        CoRe = FALSE,
        PrintPlot = FALSE)

    ## return
    invisible(list("data" = l))
}


##########################################################################
### ### ### Metabolite detection estimation using pool samples ### ### ###
##########################################################################

#' Find metabolites with high variability across total pool samples
#'
#' @param InputData DF which contains unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected. Can be either a full dataset or a dataset with only the pool samples.
#' @param SettingsFile_Sample  \emph{Optional: } DF which contains information about the samples when a full dataset is inserted as Input_data. Column "Conditions" with information about the sample conditions (e.g. "N" and "T" or "Normal" and "Tumor"), has to exist.\strong{Default = NULL}
#' @param SettingsInfo  \emph{Optional: } NULL or Named vector including the Conditions and PoolSample information (Name of the Conditions column and Name of the pooled samples in the Conditions in the Input_SettingsFile)  : c(Conditions="ColumnNameConditions, PoolSamples=NamePoolCondition. If no Conditions is added in the Input_SettingsInfo, it is assumed that the conditions column is named 'Conditions' in the Input_SettingsFile.). \strong{Default = NULL}
#' @param CutoffCV \emph{Optional: } Filtering cutoff for high variance metabolites using the Coefficient of Variation. \strong{Default = 30}
#' @param SaveAs_Plot \emph{Optional: } Select the file type of output plots. Options are svg, png, pdf or NULL. \strong{Default = svg}
#' @param SaveAs_Table \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt", ot NULL \strong{default: "csv"}
#' @param PrintPlot \emph{Optional: } If TRUE prints an overview of resulting plots. \strong{Default = TRUE}
#' @param FolderPath \emph{Optional:} Path to the folder the results should be saved at. \strong{default: NULL}
#'
#' @return List with two elements: DF (including input and output table) and Plot (including all plots generated)
#'
#' @examples
#' ## load the data
#' Intra <- ToyData("IntraCells_Raw")
#' 
#' ## create SummarizedExperiment
#' a <- t(Intra[, -c(1:3)])
#' rD <- DataFrame(feature = rownames(a))
#' cD <- Intra[, c(1:3)]
#' se <- SummarizedExperiment(assay = a, rowData = rD, colData = cD)
#' 
#' ## apply the function
#' PoolEstimation(se = se,
#'     SettingsInfo = c(Conditions = "Conditions", PoolSamples = "Pool"))
#'
#' @keywords Coefficient of Variation, high variance metabolites
#'
#' @importFrom dplyr case_when select rowwise mutate ungroup
#' @importFrom magrittr %>% %<>%
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom logger log_info log_trace
#'
#' @export
#'
PoolEstimation <- function(se, InputData,
    ##SettingsFile_Sample = NULL,
    SettingsInfo = NULL,
    CutoffCV = 30,
    SaveAs_Plot = "svg", ## EDIT: name the option here and use match.arg
    SaveAs_Table = "csv", ## EDIT: name the option here and use match.arg
    PrintPlot = TRUE,
    FolderPath = NULL) {

    ## ------------ Create log file ----------- ##
    MetaProViz_Init()

    #InputData <- assay(se)
    #SettingsFile_Sample <- colData()
    
    logger::log_info('Starting pool estimation.')
    
    ## ------------------ Check Input ------------------- ##
    # HelperFunction `CheckInput`
    CheckInput(se = se, 
        SettingsInfo = SettingsInfo,
        SaveAs_Plot = SaveAs_Plot,
        SaveAs_Table = SaveAs_Table,
        CoRe = FALSE,
        PrintPlot = PrintPlot)

    ## `CheckInput` Specific
    if ("Conditions" %in% names(SettingsInfo)) {
        if (!SettingsInfo[["Conditions"]] %in% colnames(colData(se))) {
            stop("You have chosen Conditions = ", 
                SettingsInfo[["Conditions"]], ", ", SettingsInfo[["Conditions"]],
                " was not found in SettingsFile_Sample as column. Please insert the name of the experimental conditions as stated in the SettingsFile_Sample."  )
        }
    }
    if ("PoolSamples" %in% names(SettingsInfo)) {
        if (!SettingsInfo[["PoolSamples"]] %in% colData(se)[[SettingsInfo[["Conditions"]]]]) {
            stop("You have chosen PoolSamples = ", 
                 SettingsInfo[["PoolSamples"]], ", ", 
                 SettingsInfo[["PoolSamples"]],
                " was not found in SettingsFile_Sample as sample condition. Please insert the name of the pool samples as stated in the Conditions column of the SettingsFile_Sample."  )
        }
    }

    if (!is.numeric(CutoffCV) | CutoffCV < 0) {
        stop("Check input. The selected CutoffCV value should be a positive numeric value.")
    }

    ## -----------------  create output folders  and path ------------------- ##
    if (!is.null(SaveAs_Plot) | !is.null(SaveAs_Table)) {
        Folder <- SavePath(FolderName = "Processing", FolderPath = FolderPath)

        SubFolder <- file.path(Folder, "PoolEstimation")
        logger::log_info("Selected output directory: `%s`.", SubFolder)
        if (!dir.exists(SubFolder)) {
            logger::log_trace("Creating directory: `%s`.", SubFolder)
            dir.create(SubFolder)
        }
    }

    ## ------------------ Prepare the data ------------------- ##
    ## InputData files:
    ##if (is.null(SettingsFile_Sample)) {
    ##    PoolData <- assay(se) %>%
    ##        ## make sure all 0 are changed to NAs
    ##        dplyr::mutate_all(~ ifelse(grepl("^0*(\\.0*)?$", as.character(.)), NA, .)) ## EDIT: this mutate_all should be simplified
    ##} else {
    PoolData <- assay(se)[, colData(se)[[ SettingsInfo[["Conditions"]]]] == SettingsInfo[["PoolSamples"]]]## %>%
    PoolData[PoolData == 0] <- NA
            ## Make sure all 0 are changed to NAs
            ##dplyr::mutate_all(~ ifelse(grepl("^0*(\\.0*)?$", as.character(.)), NA, .)) 
   ## }

    ############################################################################
    ## ------------------ Coefficient of Variation ------------------- ##
    logger::log_trace("Calculating coefficient of variation.")
    result_df <- apply(PoolData, 1,  
        function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE) * 100) %>% ## EDIT: use an external function to calculate CVs
        t() %>% 
        as.data.frame()
    rownames(result_df)[1] <- "CV"

    ## calculate the NAs
    NAvector <- apply(PoolData, 1, function(x) sum(is.na(x)) / length(x) * 100)

    ## create Output DF
    result_df_final <- result_df %>%
        t() %>% 
        as.data.frame() %>% 
        dplyr::rowwise() %>%
        dplyr::mutate(HighVar = CV > CutoffCV) %>%
        as.data.frame()
    result_df_final$MissingValuePercentage <- NAvector

    rownames(result_df_final) <- rownames(se)
    result_df_final_out <- tibble::rownames_to_column(result_df_final, "Metabolite")

    ## remove Metabolites from se based on CutoffCV and assign to se_filtered
    logger::log_trace('Applying CV cut-off.')
    #if (!is.null(SettingsFile_Sample)) {
    unstable_metabs <- rownames(result_df_final)[result_df_final[["HighVar"]]]
    if (length(unstable_metabs) > 0) {
        se_filtered <- se[!rownames(se) %in% unstable_metabs, ]
    ##} else {
    ##    filtered_Input_data <- NULL
    ##}
    #} else {
    #    filtered_Input_data <- NULL ## EDIT: preset filtered_Input_data <- NULL and delete the else statements
    #}

    ## ------------------ QC plots ------------------- ##
    ## start QC plot list
    logger::log_info('Plotting QC plots.')
    l_plot <- list()

    ## 1. Pool Sample PCA
    logger::log_trace('Pool sample PCA.')
    #dev.new() ## EDIT: not sure if this is needed
    #if (is.null(SettingsFile_Sample)) {
    #    pca_data <- PoolData
    #    pca_QC_pool <-invisible(VizPCA(
    #        InputData = pca_data,
    #        PlotName = "QC Pool samples",
    #        SaveAs_Plot =  NULL))
    #} else {
        #pca_data <- merge(
        #    dplyr::select(SettingsFile_Sample, SettingsInfo[["Conditions"]]), 
        #    InputData, by = 0) %>%
        #    tibble::column_to_rownames("Row.names") %>%
        #    dplyr::mutate(Sample_type = dplyr::case_when(
        #        .data[[SettingsInfo[["Conditions"]]]] == SettingsInfo[["PoolSamples"]] ~ "Pool",
        #        TRUE ~ "Sample"))

        ## add column Sample_type to se object
        se[["Sample_type"]] <- ifelse(
            colData(se)[[SettingsInfo[["Conditions"]]]] == SettingsInfo[["PoolSamples"]],
            "Pool", "Sample")
        
        ## run PCA
        pca_QC_pool <- invisible(
            VizPCA(
                se = se,
                #InputData = dplyr::select(pca_data, -all_of(SettingsInfo[["Conditions"]]), -Sample_type),
                SettingsInfo = c(color = "Sample_type"),
                ##SettingsFile_Sample = pca_data,
                PlotName = "QC Pool samples",
                SaveAs_Plot = NULL))
    }
    dev.off() ## EDIT: not sure if this is needed
    l_plot [["PCAPlot_PoolSamples"]] <- pca_QC_pool[["Plot_Sized"]][["Plot_Sized"]]


    ## 2. Histogram of CVs
    logger::log_trace('CV histogram.')
    HistCV <- suppressWarnings(invisible(
        ggplot(result_df_final_out, aes(CV)) +
            geom_histogram(aes(y = after_stat(density)), color = "black", 
                fill = "white") +
            geom_vline(aes(xintercept = CutoffCV),
               color = "darkred", linetype = "dashed", size = 1) +
            geom_density(alpha = .2, fill = "#FF6666") +
            labs(title = "CV for metabolites of Pool samples", 
                x = "Coefficient of variation (CV%)", y = "Frequency") + 
            theme_classic()))

    HistCV_Sized <- plotGrob_Processing(InputPlot =  HistCV, 
        PlotName = "CV for metabolites of Pool samples", PlotType = "Hist")
    l_plot[["Histogram_CV-PoolSamples"]] <- HistCV_Sized

    ## 2. ViolinPlot of CVs
    logger::log_trace('CV violin plot.')
    ## Make Violin of CVs
    Plot_cv_result_df <- result_df_final_out %>%
        dplyr::mutate(
            HighVar = ifelse(CV > CutoffCV, 
                paste("> CV", CutoffCV, sep = ""), 
                paste("< CV", CutoffCV, sep="")))

    ViolinCV <- invisible(
        ggplot(Plot_cv_result_df, 
            aes(y = CV, x = HighVar, label = Plot_cv_result_df$Metabolite)) +
            geom_violin(alpha = 0.5 , fill = "#FF6666") +
            geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5) +
            ggrepel::geom_text_repel(
                aes(label = ifelse(Plot_cv_result_df$CV > CutoffCV,
                    as.character(Plot_cv_result_df$Metabolite), '')),
                hjust = 0, vjust = 0,
                ## space between text and point
                box.padding = 0.5, 
                ## space around points
                point.padding = 0.5, 
                ## allow for many labels
                max.overlaps = Inf) +
            labs(title = "CV for metabolites of Pool samples", 
                x = "Metabolites", y = "Coefficient of variation (CV%)") +
            theme_classic())

    ViolinCV_Sized <- plotGrob_Processing(InputPlot = ViolinCV, 
        PlotName = "CV for metabolites of Pool samples", PlotType = "Violin")
    l_plot[["ViolinPlot_CV-PoolSamples"]] <- ViolinCV_Sized

    ############################################################################
    ## ------------------ return and save ------------------- ##
    ## save
    logger::log_info('Preparing saved and returned data.')
    if (length(unstable_metabs) > 0) {
        l_pool <- list(
            "se" = se, 
            "se_filtered" = se_filtered, 
            "CV" = result_df_final_out)
    } else {
        l_pool <- list(
            "se" = se, 
            "CV" = result_df_final_out) ## EDIT: could be simplified, define DF_list and add Filtered_InputData IF
    }

    ## save
    ##DF_list[["InputData"]] <-  DF_list[["InputData"]] %>%
    ##    tibble::rownames_to_column("Code")
    logger::log_info(
        "Saving results: [SaveAs_Table=%s, SaveAs_Plot=%s, FolderPath=%s].",
        SaveAs_Table,
        SaveAs_Plot,
        SubFolder
    )
    SaveRes(
        data = l_pool,
        plot = l_plot,
        SaveAs_Table = SaveAs_Table,
        SaveAs_Plot = SaveAs_Plot,
        FolderPath = SubFolder,
        FileName = "PoolEstimation",
        CoRe = FALSE,
        PrintPlot = PrintPlot)

    ## return
    l_res <- list("data" = l_pool, "plot" = l_plot)
    logger::log_info('Finished pool estimation.')
    invisible(l_res)
}

################################################################################################
### ### ### PreProcessing helper function: FeatureFiltering ### ### ###
################################################################################################

#' FeatureFiltering
#'
#' @param InputData DF which contains unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected and consider converting any zeros to NA unless they are true zeros.
#' @param SettingsFile_Sample DF which contains information about the samples, which will be combined with the input data based on the unique sample identifiers used as rownames.
#' @param SettingsInfo  Named vector containing the information about the names of the experimental parameters. c(Conditions="ColumnName_Plot_SettingsFile", Biological_Replicates="ColumnName_Plot_SettingsFile"). Column "Conditions" with information about the sample conditions (e.g. "N" and "T" or "Normal" and "Tumor"), can be used for feature filtering and colour coding in the PCA. Column "BiologicalReplicates" including numerical values. For CoRe = TRUE add CoRe_media = "Columnname_Input_SettingsFile", which specifies the name of the media controls in the Conditions.
#' @param CoRe \emph{Optional: } If TRUE, a consumption-release experiment has been performed.Should not be normalised to media blank. Provide information about control media sample names via SettingsInfo "CoRe_media" samples. \strong{Default = FALSE}
#' @param FeatureFilt \emph{Optional: } If NULL, no feature filtering is performed. If set to "Standard" then it applies the 80%-filtering rule (Bijlsma S. et al., 2006) on the metabolite features on the whole dataset. If is set to "Modified",filtering is done based on the different conditions, thus a column named "Conditions" must be provided in the Input_SettingsFile input file including the individual conditions you want to apply the filtering to (Yang, J et al., 2015). \strong{Default = Modified}
#' @param FeatureFilt_Value \emph{Optional: } Percentage of feature filtering. \strong{Default = 0.8}
#'
#' @return List with two elements: filtered matrix  and features filtered
#'
#' @examples
#' ## load the data
#' Intra <- ToyData("IntraCells_Raw")
#' 
#' ## create SummarizedExperiment
#' a <- t(Intra[-c(49:58), -c(1:3)])
#' rD <- DataFrame(feature = rownames(a))
#' cD <- Intra[-c(49:58), c(1:3)]
#' se <- SummarizedExperiment(assay = a, rowData = rD, colData = cD)
#' 
#' FeatureFiltering(se = se,
#'     SettingsInfo = c(Conditions = "Conditions", 
#'     Biological_Replicates = "Biological_Replicates"))
#'
#' @keywords feature filtering or modified feature filtering
#'
#' @importFrom dplyr filter mutate_all
#' @importFrom magrittr %>% %<>%
#' @importFrom logger log_info log_trace
#'
#' @noRd
#'
FeatureFiltering <- function(se, ##InputData,
    #SettingsFile_Sample,
    SettingsInfo,
    CoRe = FALSE,
    FeatureFilt = "Modified", ## EDIT: name options here and use match.arg
    FeatureFilt_Value = 0.8) {
  
    ## ------------ Create log file ----------- ##
    MetaProViz_Init()

    ## ------------------ Prepare the data ------------------- ##
    
    #feat_filt_data <- as.data.frame(InputData) %>%
    #    ## make sure all 0 are changed to NAs
    #    dplyr::mutate_all(~ ifelse(grepl("^0*(\\.0*)?$", as.character(.)), NA, .)) ## EDIT: this term is applied several times, make a function?
    se_feat_filt_data <- se
    assay(se_feat_filt_data)[assay(se_feat_filt_data) == 0] <- NA
    #assay(se)[assay(se) == 0] <- NA
    
    if (CoRe) { 
        ## remove CoRe_media samples for feature filtering
        se_feat_filt_data <- se_feat_filt_data[, 
            !colData(se_feat_filt_data)[[SettingsInfo[["Conditions"]]]] == SettingsInfo[["CoRe_media"]]]
        Feature_Filtering <- paste0(FeatureFilt, "_CoRe")
    }

    ## ------------------ Perform filtering ------------------ ##
    if (FeatureFilt == "Modified") {
        message <- paste0("FeatureFiltering: Here we apply the modified 80%-filtering rule that takes the class information (Column `Conditions`) into account, which additionally reduces the effect of missing values (REF: Yang et. al., (2015), doi: 10.3389/fmolb.2015.00004). ", 
            "Filtering value selected: ", FeatureFilt_Value)
        logger::log_info(message)
        message(message)
        
        ## obtain the updated Conditions after filtering
        feat_filt_Conditions <- colData(se_feat_filt_data)[[SettingsInfo[["Conditions"]]]]
        # if (CoRe) { 
        #     feat_filt_Conditions <- colData(se_feat_filt_data)[[SettingsInfo[["Conditions"]]]][!colData(se)[[SettingsInfo[["Conditions"]]]] == SettingsInfo[["CoRe_media"]]]
        # } else {
        #     feat_filt_Conditions <- colData(se)[[SettingsInfo[["Conditions"]]]]
        # }

        if (is.null(unique(feat_filt_Conditions))) {
            message("Conditions information is missing.")
            logger::log_trace(message)
            stop(message)
        }
        if (length(unique(feat_filt_Conditions)) ==  1) {
            message("To perform the Modified feature filtering there have to be at least 2 different Conditions in the `Condition` column in the Experimental design. Consider using the Standard feature filtering option.")
            logger::log_trace(message)
            stop(message)
        }

        miss <- c()
        ## split data frame into a list of dataframes by condition
        
        split_Input <- assay(se_feat_filt_data) |>
            t() |>
            as.data.frame() |>
            split(f = feat_filt_Conditions, drop = FALSE) 
        
        for (m in split_Input) { 
            ## select metabolites to be filtered for different conditions
            for (i in seq_len(ncol(m))) {
                if (length(which(is.na(m[, i]))) > (1 - FeatureFilt_Value) * nrow(m)) 
                    miss <- append(miss, i)
            }
        }

        if (length(miss) ==  0) { 
            ## remove metabolites if any are found
            message("There where no metabolites exluded")
            filtered_matrix <- assay(se)
            feat_file_res <- "There where no metabolites exluded"
        } else {
            names_filt <- unique(rownames(se)[miss])
            message(
                length(unique(miss)), " metabolites where removed: ", 
                paste0(names_filt, collapse = ", "))
            filtered_matrix <- assay(se)[-miss, ]
        }
    } else if (FeatureFilt ==  "Standard") {
        message <- paste0("FeatureFiltering: Here we apply the so-called 80%-filtering rule, which removes metabolites with missing values in more than 80% of samples (REF: Smilde et. al. (2005), Anal. Chem. 77, 6729–6736., doi:10.1021/ac051080y). ", 
            "Filtering value selected:", FeatureFilt_Value)
        logger::log_info(message)
        message(message)

        split_Input <- assay(se_feat_filt_data) |>
            t()

        miss <- c()
        for (i in seq_len(ncol(split_Input))) { 
            ## select metabolites to be filtered for one condition
            if (length(which(is.na(split_Input[, i]))) > (1 - FeatureFilt_Value) * nrow(split_Input))
                miss <- append(miss, i)
        }

        if (length(miss) ==  0) {
            ## remove metabolites if any are found
            message <- paste0("FeatureFiltering: There where no metabolites exluded")
            logger::log_info(message)
            message(message)

            filtered_matrix <- assay(se)
            feat_file_res <- "There where no metabolites exluded"
        } else {
            names_filt <- unique(rownames(se)[miss])
            message <- paste0(length(unique(miss)), 
                " metabolites where removed: ", paste0(names_filt, collapse = ", "))
            logger::log_info(message)
            message(message)
            filtered_matrix <- assay(se)[, -miss]
        }
    }

    ## ------------------ Return ------------------ ##
    features_filtered <- unique(rownames(se)[miss]) %>% 
        as.vector()
    #filtered_matrix <- dplyr::mutate_all( ## EDIT: why is this needed?
    #    as.data.frame(filtered_matrix), function(x) as.numeric(as.character(x))) |>
    #    as.matrix()

    ## update the SummarizedObject   
    se <- se[rownames(filtered_matrix), ]
    assay(se) <- filtered_matrix
    
    ## assemble the object to return
    l <- list(
        "se" = se,
        "assay" = filtered_matrix, 
        "RemovedMetabolites" = features_filtered)
    
    ## return
    invisible(list("data" = l))
}

################################################################################################
### ### ### PreProcessing helper function: Missing Value imputation ### ### ###
################################################################################################

#' Missing Value Imputation using half minimum value
#'
#' @param InputData DF which contains unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected and consider converting any zeros to NA unless they are true zeros.
#' @param SettingsFile_Sample DF which contains information about the samples, which will be combined with the input data based on the unique sample identifiers used as rownames.
#' @param SettingsInfo  Named vector containing the information about the names of the experimental parameters. c(Conditions="ColumnName_Plot_SettingsFile", Biological_Replicates="ColumnName_Plot_SettingsFile", CoRe_media = "Columnname_Input_SettingsFile"). Column "Conditions" with information about the sample conditions, Column "BiologicalReplicates" including numerical values and Column "Columnname_Input_SettingsFile" is used to specify the name of the media controls in the Conditions.
#' @param CoRe \emph{Optional: } If TRUE, a consumption-release experiment has been performed. Should not be normalised to media blank. Provide information about control media sample names via SettingsInfo "CoRe_media" samples.\strong{Default = FALSE}
#' @param MVI_Percentage \emph{Optional: } Percentage 0-100 of imputed value based on the minimum value. \strong{Default = 50}
#'
#' @return DF with imputed values
#'
#' @examples
#' ## load the data
#' Intra <- ToyData("IntraCells_Raw")
#' 
#' ## create SummarizedExperiment
#' a <- t(Intra[-c(49:58), -c(1:3)])
#' rD <- DataFrame(feature = rownames(a))
#' cD <- Intra[-c(49:58), c(1:3)]
#' se <- SummarizedExperiment(assay = a, rowData = rD, colData = cD)
#' 
#' ## apply the function
#' MVImputation(se = se, 
#'     SettingsInfo = c(Conditions = "Conditions", 
#'         Biological_Replicates = "Biological_Replicates"))
#'
#' @keywords Half minimum missing value imputation
#'
#' @importFrom dplyr mutate group_by 
#' @importFrom MatrixGenerics rowMins
#' @importFrom logger log_info
#'
#' @noRd
#'
MVImputation <- function(se, ##InputData,
    ##SettingsFile_Sample,
    SettingsInfo,
    CoRe = FALSE,
    MVI_Percentage = 50) {
  
    ## ------------------ Prepare the data ------------------- ##
    if (CoRe) {
        ## remove blank samples
        se <- se[, 
            !colData(se)[[SettingsInfo[["Conditions"]]]] == SettingsInfo[["CoRe_media"]]]
    }
    
    #se_NA_removed <- se
    filtered_matrix <- assay(se)
    filtered_matrix[filtered_matrix == 0] <- NA
        
    ## ------------------ Perform MVI ------------------ ##
    ## do MVI for the samples
    message <- paste0("Missing Value Imputation: Missing value imputation is performed, as a complementary approach to address the missing value problem, where the missing values are imputing using the `half minimum value`. REF: Wei et. al., (2018), Reports, 8, 663, doi:https://doi.org/10.1038/s41598-017-19120-0")
    logger::log_info(message)
    message(message)

    ## impute features
    for (feature in rownames(se)) {
        ##feature_data <- merge(t(assay(se_NA_removed[feature, ])), as.data.frame(colData(se_NA_removed)) %>% dplyr::select(Conditions), by = 0)
        
        feature_data <- cbind(
            t(assay(se)[feature, , drop = FALSE]), 
            colData(se))
        feature_compatible <- make.names(feature)

        imputed_feature_data <- feature_data %>%
            as.data.frame() |>
            dplyr::group_by(Conditions) %>%
            dplyr::mutate(across(all_of(feature_compatible), 
            ~ {
                if (all(is.na(.))) {
                    message <- paste0(
                        "For some conditions all measured samples are NA for " , 
                        feature, 
                        ". Hence we can not perform half-minimum value imputation ",
                        "per condition for this metabolite and will assume it ",
                        "is a true biological 0 in those cases.")
                    logger::log_info(message)
                    message(message)
                    ## Return NA if all values are missing
                    return(0) ## EDIT: this returns 0 instead of NA?
                } else {
                    return(replace(., is.na(.), min(., na.rm = TRUE) * (MVI_Percentage / 100)))
                }
            }))

        assay(se)[feature, ] <- imputed_feature_data[[feature_compatible]]
    }

    ## for CoRe 
    if (CoRe) {
        replaceNA <- filtered_matrix[, colnames(se)]
        
        ## find metabolites with NA
        na_percentage <- rowMeans(is.na(replaceNA)) * 100
        highNA_metabs <- na_percentage[na_percentage > 20 & na_percentage < 100] ## EDIT: I would describe in @description or @details that the value is hardcoded
        OnlyNA_metabs <- na_percentage[na_percentage == 100]

        ## report metabolites with NA
        if (sum(na_percentage) > 0) {
            message <- paste0("NA values were found in Control_media samples for metabolites. For metabolites including NAs MVI is performed unless all samples of a metabolite are NA.")
            logger::log_info(message)
            message(message)
            if (sum(na_percentage > 20 & na_percentage < 100) > 0) {
                message <- paste0(
                    "Metabolites with high NA load (>20%) in Control_media samples are: ",
                    paste(names(highNA_metabs), collapse = ", "), ".")
                logger::log_info(message)
                message(message)
            }
            if (sum(na_percentage == 100) > 0) {
                message <- paste0("Metabolites with only NAs (=100%) in Control_media samples are: ",
                    paste(names(OnlyNA_metabs), collapse = ", "), 
                    ". Those NAs are set zero as we consider them true zeros")
                logger::log_info(message)
                message(message)
            }
        }

        ## if all values are NA set to 0
        replaceNA_zero <- replaceNA
        replaceNA_zero[apply(replaceNA_zero, 1, 
            function(row) all(is.na(row))), ] <- 0

        ## if there is at least 1 value use the half minimum per feature
        replaceNA_zero_MVI <- replaceNA_zero
        half_min <- rowMins(replaceNA_zero_MVI, na.rm = TRUE) / 2
        replaceNA_zero_MVI <- replace(replaceNA_zero_MVI, 
            is.na(replaceNA_zero_MVI), 
            half_min[row(replaceNA_zero_MVI)[is.na(replaceNA_zero_MVI)]])
    
        ## add the samples in the original dataframe
        filtered_matrix <- replaceNA_zero_MVI
    } else { ## EDIT: why is this only done for CoRe = TRUE?
        filtered_matrix <- assay(se)
    }
    
    ## update the SummarizedObject   
    se <- se[rownames(filtered_matrix), ]
    assay(se) <- filtered_matrix
    
    ## assemble the object to return
    l_imputed <- list(
        "se" = se,
        "assay" = filtered_matrix)
    
    ## return
    invisible(list("data" = l_imputed))
}


################################################################################################
### ### ### PreProcessing helper function: Total ion Count Normalization ### ### ###
################################################################################################

#' Total ion count normalisazion
#'
#' @param InputData DF which contains unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected and consider converting any zeros to NA unless they are true zeros.
#' @param SettingsFile_Sample DF which contains information about the samples, which will be combined with the input data based on the unique sample identifiers used as rownames.
#' @param SettingsInfo  Named vector containing the information about the names of the experimental parameters. c(Conditions="ColumnName_Plot_SettingsFile").
#' @param TIC \emph{Optional: } If TRUE, Total Ion Count normalization is performed. If FALSE, only RLA QC plots are returned. \strong{Default = TRUE}
#'
#' @return List with two elements: DF (including output table) and Plot (including all plots generated)
#'
#' @examples
#' ## load data
#' Intra <- ToyData("IntraCells_Raw")
#' 
#' ## create SummarizedExperiment
#' a <- t(Intra[-c(49:58), -c(1:3)])
#' rD <- DataFrame(feature = rownames(a))
#' cD <- Intra[-c(49:58), c(1:3)]
#' se <- SummarizedExperiment(assay = a, rowData = rD, colData = cD)
#' 
#' ## apply the function
#' TICNorm(se = se, SettingsInfo = c(Conditions = "Conditions"))
#'
#' @keywords total ion count normalisation
#'
#' @importFrom tidyr pivot_longer
#' @importFrom gridExtra grid.arrange
#' @importFrom ggplot2 ggplot sym geom_boxplot geom_hline labs theme_classic theme_minimal theme annotation_custom aes_string
#' @importFrom logger log_info log_trace
#'
#' @noRd
#'
TICNorm <- function(se, ##InputData,
    ##SettingsFile_Sample,
    SettingsInfo,
    TIC = TRUE) {
    
    ## ------------------ Prepare the data ------------------- ##
    NA_removed <- assay(se)
    
    ## replace NA with 0
    NA_removed[is.na(NA_removed)] <- 0 

    ## ------------------ QC plot ------------------- ##
    ## before TIC Normalization
    ## log() transformation:
    log_NA_removed <- suppressWarnings(
        ## log tranform the data
        log(NA_removed)) 
    
    ## count NaN values (produced by log(0))
    nan_count <- sum(is.nan(as.matrix(log_NA_removed)))
    if (nan_count > 0) {
        ## issue a custom warning if NaNs are present
        message <- paste("For the RLA plot before/after TIC normalisation we have to perform log() transformation. This resulted in", 
            nan_count, 
            "NaN values due to 0s in the data.")
        logger::log_trace("Warning: ", message, sep="")
        warning(message)
    }

    ## get median
    median_tic <- apply(log_NA_removed, 2, median)
    
    ## Subtract the medians from each column
    data_raw <- log_NA_removed - median_tic   
    data_long <- tidyr::pivot_longer(as.data.frame(data_raw), 
        cols = everything(), names_to = "Samples", values_to = "Intensity")
    
    ## add Conditions from colData
    data_long <- merge(data_long, colData(se), 
        by.x = "Samples", by.y = "row.names")

    ## create the ggplot boxplot
    gg_data_raw <- ggplot2::ggplot(data_long, 
            ggplot2::aes(x = !!sym("Samples"), y = !!sym("Intensity"), 
                color = !!sym(SettingsInfo[["Conditions"]]))) +
        ggplot2::geom_boxplot() +
        ggplot2::geom_hline(yintercept = 0, color = "red", linetype = "solid") +
        ggplot2::labs(title = "Before TIC Normalization") +
        ggplot2::theme_classic() +
        ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        ggplot2::theme(legend.position = "none")

    if (TIC) {
        
        ## ------------------ Perform TIC normalization ------------------- ##
        message <- paste0("Total Ion Count (TIC) normalization: Total Ion ",
            "Count (TIC) normalization is used to reduce the variation from ",
            "non-biological sources, while maintaining the biological ",
            "variation. REF: Wulff et. al., (2018), Advances in Bioscience ",
            "and Biotechnology, 9, 339-351, ",
            "doi:https://doi.org/10.4236/abb.2018.98022")
        logger::log_info(message)
        message(message)

        tic <- colSums(NA_removed)
        
        ## built the median
        median_tic <- median(tic)
        
        ## divide the ion intensity by the total ion count and multiply with
        ## the median intensity
        tic_norm <- sweep(NA_removed, 2, STATS = tic, FUN = "/") * median_tic ## EDIT: is this what should be done?
        
        ## ------------------ QC plot ------------------- ##
        ### After TIC normalization
        log_tic_norm  <- suppressWarnings(
            ## log tranforms the data
            log(tic_norm)) 
        median_tic_norm <- apply(log_tic_norm, 2, median)
        
        ## Subtract the medians from each column
        data_norm <- sweep(log_tic_norm, 2, STATS = median_tic_norm, FUN = "-")
        data_long <- tidyr::pivot_longer(as.data.frame(data_norm), 
            cols = everything(), names_to = "Samples", values_to = "Intensity")
        
        ## add Conditons from colData
        data_long <- merge(data_long, 
            colData(se), by.x = "Samples", by.y = "row.names")
        
        ## Create the ggplot boxplot
        gg_data_norm <- ggplot2::ggplot(data_long, 
                ggplot2::aes(x = !!sym("Samples"), y = !!sym("Intensity"), 
                    color = !!sym(SettingsInfo[["Conditions"]]))) +
            ggplot2::geom_boxplot() +
            ggplot2::geom_hline(yintercept = 0, color = "red", linetype = "solid") +
            ggplot2::labs(title = "After TIC Normalization") +
            ggplot2::theme_classic() +
            ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
            ggplot2::theme(legend.position = "none")

        ## RLA_data_norm_Sized <- plotGrob_Processing(InputPlot = RLA_data_norm, PlotName= "After TIC Normalization", PlotType= "RLA")

        ## combine Plots
        dev.new() ## EDIT: is this needed?
        plots_combined <- suppressWarnings(gridExtra::grid.arrange(
            gg_data_raw + 
                ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
                ggplot2::theme(legend.position = "none"), 
            gg_data_norm + 
                ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
                ggplot2::theme(legend.position = "none"), ncol = 2))
        dev.off()
        plots_combined <- ggplot2::ggplot() + 
            ggplot2::theme_minimal() + 
            ggplot2::annotation_custom(plots_combined)
        
        ## update the SummarizedExperiment object
        assay(se) <- tic_norm

        ## create the object to return
        l_normalized <- list(
            "data" = list(
                "se" = se,
                "assay" = assay(se)
            ),
            "plot" = list(
                "beforeTicNormalization" = gg_data_raw,
                "afterTicNormalization" = gg_data_norm, 
                "combined" = plots_combined
            ))
        
    } else {
        ## create the object to return
        l_normalized <- list(
            "data" = list(
                "se" = se,
                "assay" = assay(se) ## EDIT: correct? this is just a pass through    
            ),
            "plot" = list(
                "beforeTicNormalization" = gg_data_raw))
    }
    
    ## return
    invisible(l_normalized)
}

################################################################################################
### ### ### PreProcessing helper function: CoRe nomalisation ### ### ###
################################################################################################

#' Consumption Release Normalisation
#'
#' @param InputData DF which contains unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected and consider converting any zeros to NA unless they are true zeros.
#' @param SettingsFile_Sample DF which contains information about the samples, which will be combined with the input data based on the unique sample identifiers used as rownames.
#' @param SettingsInfo  Named vector containing the information about the names of the experimental parameters. c(Conditions="ColumnName_Plot_SettingsFile", CoRe_norm_factor = "Columnname_Input_SettingsFile", CoRe_media = "Columnname_Input_SettingsFile"). Column CoRe_norm_factor is used for normalization and CoRe_media is used to specify the name of the media controls in the Conditions.
#'
#' @return List with two elements: DF (including output table) and Plot (including all plots generated)
#'
#' @examples
#' ## load the data
#' Media <- ToyData("CultureMedia_Raw") %>% 
#'     subset(!Conditions=="Pool")
#'     
#' ## create the SummarizedExperiment
#' a <- t(Media[, -c(1:3)])
#' rD <- DataFrame(feature = rownames(a))
#' cD <- Media[, c(1:3)]
#' se <- SummarizedExperiment(assay = a, rowData = rD, colData = cD)
#' 
#' ## apply the function
#' Res <- CoReNorm(se = se, 
#'     SettingsInfo = c(Conditions = "Conditions", 
#'         CoRe_norm_factor = "GrowthFactor", CoRe_media = "blank"))
#'
#' @keywords Consumption Release Metaqbolomics, Normalisation, Exometabolomics
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggplot2 ggplot geom_histogram geom_vline geom_density labs theme_classic geom_violin geom_dotplot
#' @importFrom logger log_info log_trace
#' @importFrom dplyr case_when summarise_all mutate rowwise mutate_all select filter pull
#'
#'
#' @noRd
#'
#'
CoReNorm <-function(
    se,
    #InputData,
    #SettingsFile_Sample,
    SettingsInfo) {
  
    ## ------------------ Prepare the data ------------------- ##
    Data_TIC <- assay(se)
    Data_TIC[is.na(Data_TIC)] <- 0

    ## ------------------ Perform QC ------------------- ##
    Conditions <- colData(se)[[SettingsInfo[["Conditions"]]]]
    CoRe_medias <- Data_TIC[, grep(SettingsInfo[["CoRe_media"]], Conditions)]

    if (ncol(CoRe_medias) == 1) {
        message <- paste0("Only 1 CoRe_media sample was found. Thus, the ",
            "consistency of the CoRe_media samples cannot be checked. It is ",
            "assumed that the CoRe_media samples are already summed.")
        logger::log_trace(paste0("Warning: ", message))
        warning(message)

        CoRe_media_df <- CoRe_medias
        colnames(CoRe_medias) <- "CoRe_mediaMeans"
        
    } else {
        ########################################################################
        ## ------------------ QC Plots
        PlotList <- list()
        
        ##-- 1. PCA Media_control
        media_pca_data <- merge(
            x =  colData(se)[, SettingsInfo[["Conditions"]], drop = FALSE], 
            y = t(Data_TIC), by = "row.names") %>%
            tibble::column_to_rownames("Row.names") %>%
            dplyr::mutate(Sample_type = dplyr::case_when(
                Conditions == SettingsInfo[["CoRe_media"]] ~ "CoRe_media",
                TRUE ~ "Sample"))

        ##media_pca_data[is.na(media_pca_data)] <- 0 ## EDIT: this is not needed when we have done it above
        se_tmp <- se
        assay(se_tmp) <- media_pca_data[, rownames(se)] |>
            t()
        se_tmp@colData <- media_pca_data[, !colnames(media_pca_data) %in% rownames(se)] |>
            DataFrame()
        
        dev.new()
        pca_QC_media <- invisible(
            VizPCA(se = se_tmp, #InputData = dplyr::select(media_pca_data, -SettingsInfo[["Conditions"]], -Sample_type), ## EDIT: why not just use t(Data_TIC)?
                SettingsInfo = c(color = "Sample_type"),
                #SettingsFile_Sample = media_pca_data,
                PlotName = "QC Media_samples",
                SaveAs_Plot =  NULL))
        dev.off()

        PlotList[["PCA_CoReMediaSamples"]] <- pca_QC_media[["Plot_Sized"]][["QC Media_samples"]]

        ##-- 2. Metabolite Variance Histogram
        ## Coefficient of Variation
        result_df <- apply(CoRe_medias, MARGIN = 1, ## EDIT: should this be on samples (MARGIN = 2) or metabolites (MARGIN = 1)
            function(x) sd(x, na.rm = TRUE) /  mean(x, na.rm = TRUE) * 100) %>% ## EDIT: I would write here a function to calculate CVs + test it
            t() %>% 
            as.data.frame()
        result_df[1, is.na(result_df[1, ])] <- 0
        rownames(result_df)[1] <- "CV"

        CutoffCV <- 30
        result_df <- result_df %>% 
            t() %>% 
            as.data.frame() %>% 
            dplyr::rowwise() %>%
            dplyr::mutate(HighVar = CV > CutoffCV) %>%
            as.data.frame()
        rownames(result_df)<- rownames(CoRe_medias) ## adjust to colnames if for samples

        ## calculate the NAs
        NAvector <- apply(CoRe_medias, 1, 
            function(x) sum(is.na(x)) / length(x) * 100) ## EDIT: I would write here a function to calculate #NAs + test it
        result_df$MissingValuePercentage <- NAvector

        cv_result_df <- result_df

        HighVar_metabs <- sum(result_df$HighVar)
        if (HighVar_metabs > 0) {
            message <- paste0(HighVar_metabs, 
                " of variables have high variability (CV > 30) in the CoRe_media control samples. Consider checking the pooled samples to decide whether to remove these metabolites or not.")
            logger::log_info(message)
            message(message)
        }

        ## Make histogram of CVs
        HistCV <- invisible(
            ggplot2::ggplot(cv_result_df, aes(CV)) +
                ggplot2::geom_histogram(aes(y = after_stat(density)), 
                    color = "black", fill = "white") +
                ggplot2::geom_vline(aes(xintercept = CutoffCV),
                    color = "darkred", linetype = "dashed", linewidth = 1) +
                ggplot2::geom_density(alpha = .2, fill = "#FF6666") +
                ggplot2::labs(
                    title = "CV for metabolites of control media samples (no cells)",
                    x = "Coefficient of variation (CV)", y = "Frequency") +
                ggplot2::theme_classic())

        HistCV_Sized <- plotGrob_Processing(
            InputPlot = HistCV, 
            PlotName = "CV for metabolites of control media samples (no cells)", 
            PlotType = "Hist")

        PlotList[["Histogram_CoReMediaCV"]] <- HistCV_Sized

        ## Make Violin of CVs
        Plot_cv_result_df <- cv_result_df %>%
            dplyr::mutate(HighVar = ifelse(HighVar, "> CV 30", "< CV 30"))

        ViolinCV <- invisible(
            ggplot2::ggplot(Plot_cv_result_df, 
                    aes(y=CV, x = HighVar, label = row.names(cv_result_df))) +
                ggplot2::geom_violin(alpha = 0.5 , fill = "#FF6666")+
                ggplot2::geom_dotplot(binaxis = "y", stackdir = "center", 
                    dotsize = 0.5) +
                ggrepel::geom_text_repel(
                    aes(label = ifelse(Plot_cv_result_df$CV > CutoffCV,
                        as.character(row.names(Plot_cv_result_df)), "")),
                    hjust = 0, vjust = 0,
                    ## space between text and point
                    box.padding = 0.5, 
                    ## space around points
                    point.padding = 0.5, 
                    ## allow for many labels
                    max.overlaps = Inf) + 
                ggplot2::labs(
                    title = "CV for metabolites of control media samples (no cells)",
                    x = "Metabolites", y = "Coefficient of variation (CV)") +
                ggplot2::theme_classic())

        ViolinCV_Sized <- plotGrob_Processing(
            InputPlot = ViolinCV, 
            PlotName = "CV for metabolites of control media samples (no cells)", 
            PlotType = "Violin")
        PlotList[["CoRe_Media_CV_Violin"]] <- ViolinCV_Sized

        ########################################################################
        ## ------------------ Outlier testing
        if (ncol(CoRe_medias) >= 3) {
            Outlier_data <- CoRe_medias
            Outlier_data <- Outlier_data %>% 
                as.data.frame() |>
                dplyr::mutate_all(.funs = ~ FALSE)

            while(HighVar_metabs > 0) {
                ## remove the furthest value from the mean
                if (HighVar_metabs > 1) { 
                    max_var_pos <-  CoRe_medias[result_df$HighVar, ]  %>%
                        t() |>
                        as.data.frame() %>%
                        dplyr::mutate_all(.funs = ~ . - mean(., na.rm = TRUE)) %>%
                        dplyr::summarise_all(.funs = ~ which.max(abs(.)))
                } else {
                    max_var_pos <-  CoRe_medias[result_df$HighVar]  %>%
                        t() |>
                        as.data.frame() %>%
                        dplyr::mutate_all(.funs = ~ . - mean(., na.rm = TRUE)) %>%
                        dplyr::summarise_all(.funs = ~ which.max(abs(.)))
                    colnames(max_var_pos)<- rownames(CoRe_medias)[result_df$HighVar] ## EDIT: is this the only difference? Then you can wrap "colnames(max_var_pos) <- ..." in the if/else and everything before the if/else
                }

                ## Remove rows based on positions
                for (i in seq_along(max_var_pos)) {
                    CoRe_medias[names(max_var_pos)[i], max_var_pos[[i]]] <- NA
                    Outlier_data[names(max_var_pos)[i], max_var_pos[[i]]] <- TRUE
                }
    
                ## recalculate coefficient of variation for each column in the 
                ## filtered data
                result_df <- apply(CoRe_medias, MARGIN = 1, 
                        function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)) %>% ## EDIT: use here a function instead (see above)
                    t() %>% 
                    as.data.frame() 
                result_df[1, is.na(result_df[1, ])] <- 0
                rownames(result_df)[1] <- "CV"
    
                result_df <- result_df %>% 
                    t() %>%
                    as.data.frame() %>% 
                    dplyr::rowwise() %>%
                    dplyr::mutate(HighVar = CV > CutoffCV) %>% 
                    as.data.frame()
                rownames(result_df)<- rownames(CoRe_medias)
    
                HighVar_metabs <- sum(result_df$HighVar)
            }

            data_cont <- Outlier_data %>% 
                #t() %>% 
                as.data.frame()

            ## list to store results
            fisher_test_results <- list()
            large_contingency_table <- matrix(0, nrow = 2, 
                ncol = ncol(data_cont))

            for (i in seq_along(colnames(data_cont))) {
                sample <- colnames(data_cont)[i]
                current_sample <- data_cont[, sample]

                contingency_table <- matrix(0, nrow = 2, ncol = 2)
                contingency_table[1, 1] <- sum(current_sample)
                contingency_table[2, 1] <- sum(!current_sample)
                contingency_table[1, 2] <- sum(rowSums(data_cont) - current_sample)
                contingency_table[2, 2] <- nrow(dplyr::select(data_cont, !all_of(sample))) * 
                    ncol(dplyr::select(data_cont, !all_of(sample))) - 
                    sum(rowSums(data_cont) - current_sample)

                ## Fisher's exact test
                fisher_test_result <- fisher.test(contingency_table)
                fisher_test_results[[sample]] <- fisher_test_result

                ## calculate the sum of "TRUE" and "FALSE" for the current sample
                ## sum of "TRUE"
                large_contingency_table[1, i] <- sum(current_sample)
                ## Sum of "FALSE"
                large_contingency_table[2, i] <- sum(!current_sample)
            }

            ## convert the matrix into a data_contframe for better readability
            contingency_data_contframe <- as.data.frame(large_contingency_table)
            colnames(contingency_data_contframe) <- colnames(data_cont)
            rownames(contingency_data_contframe) <- c("HighVar", "Low_var")

            contingency_data_contframe <- contingency_data_contframe %>% 
                dplyr::mutate(Total = rowSums(contingency_data_contframe))
            contingency_data_contframe <- rbind(contingency_data_contframe, 
                Total = colSums(contingency_data_contframe))

            different_samples <- c()
            for (sample in colnames(data_cont)) {
                p_value <- fisher_test_results[[sample]]$p.value
                if (p_value < 0.05) {
                    ## adjust the significance level as needed 
                    different_samples <- c(different_samples, sample) ## EDIT: where is the adjustment taking place?
                }
            }

            if (!is.null(different_samples)) {
                message <- paste("The CoRe_media samples ", 
                    paste(different_samples, collapse = ", "), 
                    " were found to be different from the rest. They will not be included in the sum of the CoRe_media samples.")
                logger::log_trace("Warning: " , message, sep = "")
                warning(message)
            }
            
            ## filter the CoRe_media samples
            CoRe_medias <- CoRe_medias %>% 
                as.data.frame() |>
                dplyr::select(-all_of(different_samples))
        } else {
            message <- paste0(
                "Only >=2 blank samples available. Thus,we can not perform outlier testing for the blank samples.")
            logger::log_trace(message)
            message(message)
        }
        CoRe_media_df <- as.data.frame( ## EDIT: why convert a data.frame to a data.frame? Can the as.data.frame be removed?
            data.frame("CoRe_mediaMeans" = rowMeans(CoRe_medias, na.rm = TRUE)))
    }

    cv_result_df <- tibble::rownames_to_column(cv_result_df, "Metabolite")

    ############################################################################
    ##------------------------ Substract mean (media control) from samples
    message <- paste0("CoRe data are normalised by substracting mean ",
        "(blank) from each sample and multiplying with the CoRe_norm_factor")
    logger::log_info(message)
    message(message)

    ##-- Check CoRe_norm_factor
    if ("CoRe_norm_factor" %in% names(SettingsInfo)) {
        CoRe_norm_factor <- colData(se) %>% 
            as.data.frame() |>
            dplyr::filter(!!as.name(SettingsInfo[["Conditions"]]) != SettingsInfo[["CoRe_media"]]) %>% 
            dplyr::select(SettingsInfo[["CoRe_norm_factor"]]) %>%
            dplyr::pull()
        
            if (var(CoRe_norm_factor) ==  0) {
                message <- paste("The growth rate or growth factor for normalising the CoRe result, is the same for all samples")
                logger::log_trace("Warning: " , message)
                warning(message)
            }
    } else {
        CoRe_norm_factor <- rep(1, 
            sum(colData(se)[, SettingsInfo[["Conditions"]]] != SettingsInfo[["CoRe_media"]]))
    }

    ## remove CoRe_media samples from the data
    Data_TIC <- merge(colData(se), t(Data_TIC), by = "row.names") %>%
        dplyr::filter(!!as.name(SettingsInfo[["Conditions"]]) != SettingsInfo[["CoRe_media"]]) %>%
        tibble::column_to_rownames("Row.names") %>%
        dplyr::select(-c(seq_len(ncol(colData(se)))))

    ## subtract from each sample the CoRe_media mean
    Data_TIC_CoReNorm_Media <- t(
        apply(t(Data_TIC), 2, function(i) i - CoRe_media_df$CoRe_mediaMeans)) ## EDIT: is this correct?
    Data_TIC_CoReNorm <- t(
        apply(Data_TIC_CoReNorm_Media, 2, function(i) i * CoRe_norm_factor)) ## EDIT: is this correct? 

    ## remove CoRe_media samples from the data
    # Input_SettingsFile <- Input_SettingsFile[Input_SettingsFile$Conditions!="CoRe_media",]
    # Conditions <- Conditions[!Conditions=="CoRe_media"]

    ## update the SummarizedExperiment object
    se <- se[, colnames(Data_TIC_CoReNorm)]
    assay(se) <- Data_TIC_CoReNorm
    
    ############################################################################
    ##------------------------ Return Plots and Data
    if (nrow(CoRe_medias) >= 3) {
        l_corenorm <- list(
            "data" = list(
                "CV_CoRe_blank" = cv_result_df, 
                "Contigency_table_CoRe_blank" = contingency_data_contframe, 
                "se" = se   
            ),
            "plot" = PlotList)
    } else {
        l_corenorm <- list(
            "data" = list(
                "CV_CoRe_blank" = cv_result_df, 
                "se" = se    
            ),
            "plot" = PlotList) ## EDIT: define this list before and update with "Contigency_table_CoRe_blank" if TRUE
    }
    
    ## return
    invisible(l_corenorm)
}


################################################################################################
### ### ### PreProcessing helper function: Outlier detection ### ### ###
################################################################################################

#' @name OutlierDetection
#'
#' @param InputData DF which contains unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected and consider converting any zeros to NA unless they are true zeros.
#' @param SettingsFile_Sample DF which contains information about the samples, which will be combined with the input data based on the unique sample identifiers used as rownames.
#' @param SettingsInfo  Named vector containing the information about the names of the experimental parameters. c(Conditions="ColumnName_Plot_SettingsFile", Biological_Replicates="ColumnName_Plot_SettingsFile", CoRe_media = "Columnname_Input_SettingsFile"). Column "Conditions" with information about the sample conditions, Column "BiologicalReplicates" including numerical values and Column "Columnname_Input_SettingsFile" is used to specify the name of the media controls in the Conditions.
#' @param CoRe \emph{Optional: } If TRUE, a consumption-release experiment has been performed. If not normalised yet, provide information about control media sample names via SettingsInfo "CoRe_media" samples. \strong{Default = FALSE}
#' @param HotellinsConfidence \emph{Optional: } Confidence level for Hotellin's T2 test. \strong{Default = 0.99}
#'
#' @return List with two elements: : DF (including output tables) and Plot (including all plots generated)
#'
#' @examples
#' ## load the data
#' Intra <- ToyData("IntraCells_Raw")
#' 
#' ## create SummarizedExperiment object
#' a <- t(Intra[-c(49:58), -c(1:3)])
#' rD <- DataFrame(feature = rownames(a))
#' cD <- Intra[-c(49:58), c(1:3)]
#' se <- SummarizedExperiment(assay = a, rowData = rD, colData = cD)
#' 
#' ## apply the function
#' Res <- OutlierDetection(se = se,
#'     SettingsInfo = c(Conditions = "Conditions", 
#'         Biological_Replicates = "Biological_Replicates"))
#'
#' @keywords Hotellins T2 outlier detection
#'
#' @importFrom dplyr relocate case_when mutate mutate_all select
#' @importFrom inflection uik
#' @importFrom factoextra fviz_screeplot
#' @importFrom ggplot2 ggplot theme_classic theme geom_vline annotate geom_line geom_point geom_hline scale_y_continuous ggtitle scale_linetype_discrete
#' @importFrom qcc mqcc
#' @importFrom magrittr %>% %<>%
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom hash values keys
#' @importFrom logger log_info log_trace
#'
#' @noRd
#'
OutlierDetection <-function(se, ##InputData,
    ##SettingsFile_Sample,
    SettingsInfo,
    CoRe = FALSE,
    HotellinsConfidence = 0.99) {
  
    ## create and log message
    message <- paste(
        "Outlier detection: Identification of outlier samples is performed ",
        "using Hotellin's T2 test to define sample outliers in a mathematical ",
        "way (Confidence = 0.99 ~ p.val < 0.01) (REF: Hotelling, H. (1931), ",
        "Annals of Mathematical Statistics. 2 (3), 360–378, ",
        "doi:https://doi.org/10.1214/aoms/1177732979). ",
        "HotellinsConfidence value selected: ", HotellinsConfidence, sep = "")
    logger::log_info(message)
    message(message)

    ## load the data
    data_norm <- assay(se)# %>%
        #dplyr::mutate_all(~ replace(., is.nan(.), 0)) ## EDIT: this is repeated and should be done via a function
    
    ## replace NA with 0
    data_norm[is.na(data_norm)] <- 0 ## EDIT: what is the difference to the call before?

    if (CoRe) {
        Conditions <- colData(se)[[SettingsInfo[["Conditions"]]]][!colData(se)[[SettingsInfo[["Conditions"]]]] == SettingsInfo[["CoRe_media"]]]
    } else {
        Conditions <- colData(se)[[SettingsInfo[["Conditions"]]]] ## EDIT: define Conditions as here outside the if/else and truncate if TRUE as above
    }

    ## prepare the lists to store the results
    ## do 10 rounds of hotelling filtering
    Outlier_filtering_loop <- 10 
    sample_outliers <- list()
    scree_plot_list <- list()
    outlier_plot_list <- list()
    metabolite_zero_var_total_list <- list()
    zero_var_metab_warning = FALSE

    #################################################
    ##--------- Perform Outlier testing:
    for (loop in seq_len(Outlier_filtering_loop)) {
        ##--- Zero variance metabolites
        # calculate each metabolites variance
        metabolite_var <- apply(data_norm, 1, var) %>% 
            t() %>% 
            as.data.frame()
        # take the names of metabolites with zero variance and puts them in list
        metabolite_zero_var_list <- list(
            colnames(metabolite_var)[which(metabolite_var[1, ] == 0)]) 

        if (sum(metabolite_var[1, ] == 0) == 0) {
            metabolite_zero_var_total_list[loop] <- 0
        } else if (sum(metabolite_var[1, ] == 0) > 0) {
            metabolite_zero_var_total_list[loop] <- metabolite_zero_var_list
            ## this is used later to print and save the zero variance 
            ## metabolites if any are found
            zero_var_metab_warning <- TRUE 
        }

        for (metab in metabolite_zero_var_list) {  
            ## remove the metabolites with zero variance from the data to do PCA
            data_norm <- data_norm[!rownames(data_norm) %in% metab]
        }

        ##---  PCA
        PCA.res <- prcomp(t(data_norm), center =  TRUE, scale. =  TRUE)
        se_tmp <- se[rownames(data_norm), ]
        assay(se_tmp) <- data_norm
        #outlier_PCA_data <- data_norm
        #outlier_PCA_data$Conditions <- Conditions

        dev.new()
        pca_outlier <- invisible(
            VizPCA(se = se_tmp, 
                SettingsInfo = c(color = SettingsInfo[["Conditions"]]),
                ##SettingsFile_Sample = outlier_PCA_data,
                PlotName = paste("PCA outlier test filtering round ", loop),
                SaveAs_Plot = NULL))

        if (loop == 1) {
            pca_outlierloop1 <- pca_outlier[["Plot_Sized"]][[1]]
        }
        outlier_plot_list[[paste0("PCA_round", loop)]] <- pca_outlier[["Plot_Sized"]][[1]]
        dev.off()

        ##--- Scree plot
        ## get Scree plot values for inflection point calculation
        inflect_df <- as.data.frame(seq_along(PCA.res$sdev)) 
        colnames(inflect_df) <- "x"
        inflect_df$y <- summary(PCA.res)$importance[2, ]
        inflect_df$Cumulative <- summary(PCA.res)$importance[3, ]
        
        ## make cumulative variation labels for plot
        screeplot_cumul <- format(round(
            inflect_df$Cumulative[seq_len(20)] * 100, 1), nsmall = 1) 
        
        ## Calculate the knee and select optimal number of components
        knee <- inflection::uik(inflect_df$x, inflect_df$y)
        ## subtract 1 components from the knee cause the root of the knee is 
        ## the PC that does not add something. npcs = 30
        npcs <- knee -1 

        ## make a scree plot with the selected component cut-off for 
        ## HotellingT2 test
        screeplot <- factoextra::fviz_screeplot(PCA.res, 
            main = paste0("PCA Explained variance plot filtering round ", loop),
            addlabels = TRUE, ncp = 20, geom = c("bar", "line"),
            barfill = "grey", barcolor = "grey", linecolor = "black",
            linetype = 1) +
            ggplot2::theme_classic()+
            ggplot2::geom_vline(xintercept = npcs + 0.5, linetype = 2, 
                color = "red") +
            ggplot2::annotate("text", x = c(1:20), y = -0.8, 
                label = screeplot_cumul, col = "black", size = 1.75)

        #screeplot_Sized <- plotGrob_Processing(InputPlot = screeplot, PlotName= paste("PCA Explained variance plot filtering round ",loop, sep = ""), PlotType= "Scree")

        if (loop == 1) {
            scree_outlierloop1 <- screeplot
        }
        dev.new()
        ## save plot
        outlier_plot_list[[paste("ScreePlot_round",loop,sep="")]] <- screeplot 
        dev.off()

        ##--- HotellingT2 test for outliers
        data_hot <- as.matrix(PCA.res$x[, seq_len(npcs)])
        hotelling_qcc <- qcc::mqcc(data_hot, type = "T2.single",
            labels = rownames(data_hot), 
            confidence.level = HotellinsConfidence, 
            title = paste0(
                "Outlier filtering via HotellingT2 test filtering round ", 
                loop, ", with ", HotellinsConfidence, "% Confidence"), 
            plot = FALSE)
        HotellingT2plot_data <- as.data.frame(hotelling_qcc$statistics) |>
            tibble::rownames_to_column("Samples")
        colnames(HotellingT2plot_data) <- c("Samples", "Group summary statistics")
        outlier <- HotellingT2plot_data %>% 
            dplyr::filter(HotellingT2plot_data$`Group summary statistics` > hotelling_qcc$limits[2]) ## EDIT: choose a name that does not need ticks
        limits <- as.data.frame(hotelling_qcc$limits)
        legend <- colnames(HotellingT2plot_data[2])
        LegendTitle <- "Limits"

        HotellingT2plot <- ggplot2::ggplot(HotellingT2plot_data, 
                aes(x = Samples, y = `Group summary statistics`, 
                    group = 1, fill =)) + ## EDIT: fill is missing
            ggplot2::geom_point(aes(x = Samples, y = `Group summary statistics`), 
                color = 'blue', size = 2) +
            ggplot2::geom_point(data = outlier, 
                aes(x = Samples, y = `Group summary statistics`), 
                color = 'red', size = 3) +
            ggplot2::geom_line(linetype = 2)

        ## draw the horizontal lines corresponding to the LCL, UCL
        HotellingT2plot <- HotellingT2plot +
            ggplot2::geom_hline(aes(yintercept = limits[, 1]), 
                color = "black", data = limits,  show.legend = FALSE) +
            ggplot2::geom_hline(aes(yintercept = limits[, 2], linetype = "UCL"),
                color = "red", data = limits, show.legend = TRUE) +
            ggplot2::scale_y_continuous(breaks = sort(
                c(ggplot_build(HotellingT2plot)$layout$panel_ranges[[1]]$y.major_source, 
                    c(limits[, 1],limits[, 2]))))

        HotellingT2plot <- HotellingT2plot +
            ggplot2::theme_classic() +
            ggplot2::theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
            ggplot2::ggtitle(
                paste("Hotelling ", hotelling_qcc$type , 
                    " test filtering round ", loop, ", with ", 
                    100 * hotelling_qcc$confidence.level, "% Confidence")) +
            ggplot2::scale_linetype_discrete(name = LegendTitle,) + 
            ggplot2::theme(plot.title = element_text(size = 13)) + #, face = "bold")) +
            ggplot2::theme(axis.text = element_text(size = 7))
        #HotellingT2plot_Sized <- plotGrob_Processing(InputPlot = HotellingT2plot, PlotName= paste("Hotelling ", hotelling_qcc$type ," test filtering round ",loop,", with ", 100 * hotelling_qcc$confidence.level,"% Confidence"), PlotType= "Hotellings")

        if (loop == 1) {
            hotel_outlierloop1 <- HotellingT2plot
        }
        dev.new()
        #plot(HotellingT2plot)
        outlier_plot_list[[paste0("HotellingsPlot_round", loop)]] <- HotellingT2plot
        dev.off()

        a <- loop
        if (CoRe) {
            a <- paste0(a, "_CoRe")
        }

        ## loop for outliers until no outlier is detected
        if (length(hotelling_qcc[["violations"]][["beyond.limits"]]) == 0) {
            data_norm <- data_norm
            break ## EDIT: why is break needed here?
        } else if (length(hotelling_qcc[["violations"]][["beyond.limits"]]) == 1) {
            ## filter the selected outliers from the data
            data_norm <- data_norm[, -hotelling_qcc[["violations"]][["beyond.limits"]] ]
            Conditions <- Conditions[-hotelling_qcc[["violations"]][["beyond.limits"]]]

            ## Change the names of outliers in mqcc, instead of saving the 
            ## order number it saves the name
            hotelling_qcc[["violations"]][["beyond.limits"]][1] <- rownames(data_hot)[hotelling_qcc[["violations"]][["beyond.limits"]][1]]
            sample_outliers[loop] <- list(hotelling_qcc[["violations"]][["beyond.limits"]])
        } else {
            data_norm <- data_norm[, -hotelling_qcc[["violations"]][["beyond.limits"]]]
            Conditions <- Conditions[-hotelling_qcc[["violations"]][["beyond.limits"]]]

            ## change the names of outliers in mqcc, instead of saving the 
            ## order number it saves the name
            sm_out <- c() # list of outliers samples
            for (i in seq_along(hotelling_qcc[["violations"]][["beyond.limits"]])) {
                sm_out <-  append(sm_out, 
                    rownames(data_hot)[hotelling_qcc[["violations"]][["beyond.limits"]][i]])
            }
            sample_outliers[loop] <- list(sm_out)
        }
    }

    #################################################
    ##-- Print Outlier detection results about samples and metabolites
    if (length(sample_outliers) > 0) {  
        ## print outlier samples
        message <- "There are possible outlier samples in the data."
        logger::log_info(message)
        message(message) #This was a warning
        for (i in seq_along(sample_outliers)) {
            message <- paste("Filtering round ", i, " Outlier Samples: ", 
                paste(head(sample_outliers[[i]]), " "))
            logger::log_info(message)
            message(message)
        }
    } else {
        message <- "No sample outliers were found."
        logger::log_info(message)
        message(message)
    }

    ##--  Print Zero variance metabolites
    zero_var_metab_export_df <- data.frame(1, 2)
    names(zero_var_metab_export_df) <- c("Filtering round", "Metabolite")

    if (zero_var_metab_warning) {
        message <- paste("Metabolites with zero variance have been identified in the data. As scaling in PCA cannot be applied when features have zero variace, these metabolites are not taken into account for the outlier detection and the PCA plots.")
        logger::log_trace("Warning: " , message, sep = "")
        warning(message)
    }

    count <- 1
    for (i in seq_along(metabolite_zero_var_total_list)) {
        if (metabolite_zero_var_total_list[[i]] != 0) {
            message <- paste("Filtering round ", i, 
                ". Zero variance metabolites identified: ", 
                paste(metabolite_zero_var_total_list[[i]], " "))
            logger::log_info(message)
            message(message)

            zero_var_metab_export_df[count, "Filtering round"] <- paste(i)
            zero_var_metab_export_df[count, "Metabolite"] <- paste(metabolite_zero_var_total_list[[i]])
            count <- count +1
        }
    }

    #############################################
    ##---- 1. Make Output DF
    ## make a dictionary
    total_outliers <- hash::hash() 
    if (length(sample_outliers) > 0) { 
        ## create columns with outliers to merge to output dataframe
        for (i in seq_along(sample_outliers) ) { ## EDIT: with seq_along probably you do not need the outer if
            total_outliers[[paste0("Outlier_filtering_round_", i)]] <- sample_outliers[i]
        }
    }

    data_norm_filtered_full <- assay(se) |>
        t() |>
        as.data.frame()
    data_norm_filtered_full[data_norm_filtered_full == 0] <- NA
    data_norm_filtered_full$Outliers <- "no"
    
    ## add outlier information to the full output dataframe
    for (i in seq_along(total_outliers)) { ## EDIT: with seq_along probably you do not need the outer if
        for (k in seq_along(hash::values(total_outliers)[i])) {
            data_norm_filtered_full[as.character(hash::values(total_outliers)[[i]]), "Outliers"] <- hash::keys(total_outliers)[i]
        }
    }

    ## put Outlier columns in the front
    data_norm_filtered_full <- data_norm_filtered_full %>% 
        dplyr::relocate(Outliers) 
    
    ## add the design in the output df (merge by rownames/sample names)
    data_norm_filtered_full <- merge(as.data.frame(colData(se)), 
        data_norm_filtered_full,  by = "row.names")
    data_norm_filtered_full <- tibble::column_to_rownames(data_norm_filtered_full, "Row.names")

    ##-- 2.  Quality Control (QC) PCA
    MetaData_Sample <- data_norm_filtered_full %>%
        dplyr::mutate(Outliers = dplyr::case_when(
            Outliers == "no" ~ 'no',
            Outliers == "Outlier_filtering_round_1" ~ ' Outlier_filtering_round = 1',
            Outliers == "Outlier_filtering_round_2" ~ ' Outlier_filtering_round = 2',
            Outliers == "Outlier_filtering_round_3" ~ ' Outlier_filtering_round = 3',
            Outliers == "Outlier_filtering_round_4" ~ ' Outlier_filtering_round = 4',
            TRUE ~ 'Outlier_filtering_round = or > 5'))
    MetaData_Sample$Outliers <- relevel(
        as.factor(MetaData_Sample$Outliers), ref = "no")

    ## define updated SummarizedExperiment for plotting
    se_tmp <- se[-zero_var_metab_export_df$Metabolite, ]
    se_tmp@colData <- MetaData_Sample[, !colnames(MetaData_Sample) %in% rownames(se)] |>
        as.data.frame() |>
        tibble::column_to_rownames(var = "Row.names") |>
        DataFrame()
    
    ## 1. Shape Outliers
    if (length(sample_outliers) > 0) {
        dev.new()
        pca_QC <- invisible(
            VizPCA(
                se = se_tmp,
                ##InputData = dplyr::select(as.data.frame(InputData), -zero_var_metab_export_df$Metabolite),
                SettingsInfo = c(color = SettingsInfo[["Conditions"]], shape = "Outliers"),
                ##SettingsFile_Sample = MetaData_Sample ,
                PlotName = "Quality Control PCA Condition clustering and outlier check",
                SaveAs_Plot = NULL))
        dev.off()
        outlier_plot_list[["QC_PCA_and_Outliers"]] <- pca_QC[["Plot_Sized"]][[1]]
    }

    ## 2. Shape Biological replicates
    if ("Biological_Replicates" %in% names(SettingsInfo)) {
        dev.new()
        pca_QC_repl <- invisible(
            VizPCA(
                se = se_tmp,
                #InputData = dplyr::select(as.data.frame(InputData), -zero_var_metab_export_df$Metabolite),
                SettingsInfo = c(color = SettingsInfo[["Conditions"]], shape = SettingsInfo[["Biological_Replicates"]]),
                #SettingsFile_Sample = MetaData_Sample,
                PlotName = "Quality Control PCA replicate spread check",
                SaveAs_Plot =  NULL))
        dev.off()
        outlier_plot_list[["QC_PCA_Replicates"]] <- pca_QC_repl[["Plot_Sized"]][[1]]
    }

    ## add data_norm_filtered_full to assay
    tmp <- data_norm_filtered_full |>
        tibble::column_to_rownames(var = "Row.names")
    assay(se_tmp) <- tmp[, rownames(se_tmp)] |>
        t()
    
    #############################################
    ##--- Save and Return plots and DFs
    l_outlier <- list(
        "data" = list(
            "se" = se_tmp,
            "Zero_variance_metabolites_CoRe" = zero_var_metab_export_df  
        ),
        "plot" = outlier_plot_list)#, 
        #"data_outliers" = data_norm_filtered_full) ## EDIT: correct? For PCA InputData is used?

    ## return
    invisible(l_outlier)
}
