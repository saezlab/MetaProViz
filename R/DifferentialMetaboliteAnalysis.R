## ---------------------------
##
## Script name: DMA
##
## Purpose of script: Differential Metabolomics Analysis
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



########################################################
### ### ### Differential Metabolite Analysis ### ### ###
########################################################

#' This function allows you to perform differential metabolite analysis to obtain a Log2FC, pval, padj and tval comparing two or multiple conditions.
#'
#' @param InputData DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected.
#' @param SettingsFile_Sample DF which contains metadata information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames.
#' @param SettingsInfo  \emph{Optional: } Named vector including the information about the conditions column information on numerator or denominator c(Conditions="ColumnName_SettingsFile", Numerator = "ColumnName_SettingsFile", Denominator  = "ColumnName_SettingsFile"). Denominator and Numerator will specify which comparison(s) will be done (one-vs-one, all-vs-one, all-vs-all), e.g. Denominator=NULL and Numerator =NULL selects all the condition and performs multiple comparison all-vs-all. Log2FC are obtained by dividing the numerator by the denominator, thus positive Log2FC values mean higher expression in the numerator. \strong{Default = c(conditions="Conditions", numerator = NULL, denumerator = NULL)}
#' @param StatPval \emph{Optional: } String which contains an abbreviation of the selected test to calculate p.value. For one-vs-one comparisons choose t.test, wilcox.test, "chisq.test", "cor.test" or lmFit (=limma), for one-vs-all or all-vs-all comparison choose aov (=anova), welch(=welch anova), kruskal.test or lmFit (=limma) \strong{Default = "lmFit"}
#' @param StatPadj \emph{Optional: } String which contains an abbreviation of the selected p.adjusted test for p.value correction for multiple Hypothesis testing. Search: ?p.adjust for more methods:"BH", "fdr", "bonferroni", "holm", etc.\strong{Default = "fdr"}
#' @param SettingsFile_Metab \emph{Optional: } DF which contains the metadata information , i.e. pathway information, retention time,..., for each metabolite. The row names must match the metabolite names in the columns of the InputData. \strong{Default = NULL}
#' @param CoRe \emph{Optional: } TRUE or FALSE for whether a Consumption/Release  input is used. \strong{Default = FALSE}
#' @param VST TRUE or FALSE for whether to use variance stabilizing transformation on the data when linear modeling is used for hypothesis testing. \strong{Default = FALSE}
#' @param PerformShapiro TRUE or FALSE for whether to perform the shapiro.test and get informed about data distribution (normal versus not-normal distribution. \strong{Default = TRUE}
#' @param PerformBartlett TRUE or FALSE for whether to perform the bartlett.test. \strong{Default = TRUE}
#' @param Transform TRUE or FALSE. If TRUE we expect the data to be not log2 transformed and log2 transformation will be performed within the limma function and Log2FC calculation. If FALSE we expect the data to be log2 transformed as this impacts the Log2FC calculation and limma. \strong{Default= TRUE}
#' @param SaveAs_Plot \emph{Optional: } Select the file type of output plots. Options are svg, png, pdf. \strong{Default = svg}
#' @param SaveAs_Table \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt". \strong{Default = "csv"}
#' @param PrintPlot \emph{Optional: } TRUE or FALSE, if TRUE Volcano plot is saved as an overview of the results. \strong{Default = TRUE}
#' @param FolderPath \emph{Optional:} Path to the folder the results should be saved at. \strong{Default = NULL}
#'
#' @return Dependent on parameter settings, list of lists will be returned for DMA (DF of each comparison), Shapiro (Includes DF and Plot), Bartlett (Includes DF and Histogram), VST (Includes DF and Plot) and VolcanoPlot (Plots of each comparison).
#'
#' @examples
#' ## load data
#' Intra <- MetaProViz::ToyData("IntraCells_Raw")
#' 
#' ## create SummarizedExperiment
#' a <- t(Intra[-c(49:58), -c(1:3)])
#' rD <- DataFrame(feature = rownames(a))
#' cD <- Intra[-c(49:58) , c(1:3)]
#' se <- SummarizedExperiment(assay = a, rowData = rD, colData = cD)
#' 
#' ## apply the function
#' ResI <- MetaProViz::DMA(se, 
#'     SettingsInfo = c(Conditions = "Conditions", 
#'         Numerator = NULL, Denominator  = "HK2"))
#'
#' @keywords Differential Metabolite Analysis, Multiple Hypothesis testing, Normality testing
#'
#' @importFrom dplyr rename
#' @importFrom magrittr %>%
#' @importFrom stats p.adjust.methods
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom purrr map reduce
#' @importFrom logger log_info
#'
#' @export
#'
DMA <-function(se, ##InputData,
    ##SettingsFile_Sample,
    SettingsInfo = c(Conditions = "Conditions", Numerator = NULL, Denominator = NULL),
    StatPval = c("lmFit", "aov", "krustal.test", "welch"),
    StatPadj = p.adjust.methods,
    ##SettingsFile_Metab = NULL,
    CoRe = FALSE,
    VST = FALSE,
    PerformShapiro = TRUE,
    PerformBartlett = TRUE,
    Transform = TRUE,
    SaveAs_Plot = "svg",
    SaveAs_Table = "csv",
    PrintPlot = TRUE,
    FolderPath = NULL) {

    ## ------------ Create log file ----------- ##
    MetaProViz_Init()

    logger::log_info("DMA: Differential metabolite analysis.")

    ## ------------ Check Input files ----------- ##
    StatPval <- match.arg(StatPval)
    StatPadj <- match.arg(StatPadj)
    
    ## HelperFunction `CheckInput`
    CheckInput(se = se, SettingsInfo = SettingsInfo,
        SaveAs_Plot = SaveAs_Plot, SaveAs_Table = SaveAs_Table, CoRe = CoRe,
        PrintPlot = PrintPlot)

    # HelperFunction `CheckInput` Specific
    Settings <- CheckInput_DMA(se, ##InputData = InputData,
        ##SettingsFile_Sample = SettingsFile_Sample, 
        SettingsInfo = SettingsInfo,
        StatPval = StatPval, StatPadj = StatPadj,
        PerformShapiro = PerformShapiro, PerformBartlett = PerformBartlett,
        VST = VST, Transform = Transform)

    ## ------------ Create Results output folder ----------- ##
    if (!is.null(SaveAs_Plot) | !is.null(SaveAs_Table)) {
        Folder <- SavePath(FolderName = "DMA", FolderPath = FolderPath)
    
        if (PerformShapiro) {
            SubFolder_S <- file.path(Folder, "Shapiro")
            if (!dir.exists(SubFolder_S)) 
                dir.create(SubFolder_S)
        }

        if (PerformBartlett) {
            SubFolder_B <- file.path(Folder, "Bartlett")
            if (!dir.exists(SubFolder_B)) 
                dir.create(SubFolder_B)
        }

        if (VST) {
            SubFolder_V <- file.path(Folder, "VST")
            if (!dir.exists(SubFolder_V))
                dir.create(SubFolder_V)
        }
    }

    ############################################################################
    ## ------------ Check hypothesis test assumptions ----------- ##
    ## 1. Normality
    if (PerformShapiro) {
        
        if (length(Settings[["Metabolites_Miss"]] >= 1)) {
            message("There are NA's/0s in the data. This can impact the output of the Shapiro-Wilk test for all metabolites that include NAs/0s.")#
        }
        
        tryCatch({
            Shapiro_output <- suppressWarnings(
                Shapiro(se = se, SettingsInfo = SettingsInfo, 
                    StatPval = StatPval, QQplots = FALSE))
            }, error = function(e) {
                message(
                    "Error occurred during Shapiro that performs the Shapiro-Wilk test. Message: ", 
                    conditionMessage(e))
            })
    }

    # 2. Variance homogeneity
    if (PerformBartlett) {
        if (Settings[["MultipleComparison"]]) { ##if we only have two conditions, which can happen even tough multiple comparison (C1 versus C2 and C2 versus C1 is done)
            
            UniqueConditions <- SettingsFile_Sample %>%
                subset(
                    SettingsFile_Sample[[SettingsInfo[["Conditions"]]]] %in% Settings[["numerator"]] | 
                        SettingsFile_Sample[[SettingsInfo[["Conditions"]]]] %in% Settings[["denominator"]], 
                    select = c(SettingsInfo[["Conditions"]]))
            UniqueConditions <- unique(
                UniqueConditions[[SettingsInfo[["Conditions"]]]])

            if (length(UniqueConditions) > 2) {
            
                tryCatch({
                    Bartlett_output <- suppressWarnings(
                        Bartlett(se = se, SettingsInfo = SettingsInfo))
                    },
                    error = function(e) {
                        message(
                            "Error occurred during Bartlett that performs the Bartlett test. Message: ", 
                            conditionMessage(e))
                    })
            }
        }
    }

    ############################################################################
    #### Prepare the data ######
    ## 1. Metabolite names:
    savedMetaboliteNames <-  data.frame("InputName" = rownames(se))
    savedMetaboliteNames$Metabolite <- paste0("M", 
        seq(1, nrow(se)))
    rownames(se) <- savedMetaboliteNames$Metabolite

    ############################################################################
    ######## Calculate Log2FC, pval, padj, tval and add additional info ########
    Log2FC_table <- Log2FC_fun(se, SettingsInfo = SettingsInfo,
        CoRe = CoRe, Transform = Transform)

    ############################################################################
    ############### Perform Hypothesis testing ###############
    if (!Settings[["MultipleComparison"]]) {
        
        if (StatPval == "lmFit") {
            STAT_C1vC2 <- DMA_Stat_limma(se, 
                SettingsInfo = SettingsInfo, StatPadj = StatPadj,
                Log2FC_table = Log2FC_table, CoRe = CoRe, Transform = Transform)

        } else {
            STAT_C1vC2 <- DMA_Stat_single(se, 
                SettingsInfo = SettingsInfo, StatPadj = StatPadj, 
                Log2FC_table = Log2FC_table, StatPval = StatPval)
        }
    } else { ## MultipleComparison = TRUE
    
        ## Correct data heteroscedasticity
        if (StatPval != "lmFit" & VST) {
            VST_res <- vst(se)
            se <- VST_res[["data"]][["se"]]
        }

        if (Settings[["all_vs_all"]]) {
            message("No conditions were specified as numerator or denumerator. Performing multiple testing `all-vs-all` using ", 
                paste(StatPval), ".")
        } else { ## for 1 vs all
            message("No condition was specified as numerator and ", 
                Settings[["denominator"]], 
                " was selected as a denominator. Performing multiple testing `all-vs-one` using ", 
                paste(StatPval), ".")
        }

        if (StatPval == "aov") {
            STAT_C1vC2 <- AOV(se, 
                SettingsInfo = SettingsInfo, Log2FC_table = Log2FC_table)
        } else if (StatPval == "kruskal.test") {
            STAT_C1vC2 <- Kruskal(se,
                SettingsInfo = SettingsInfo, Log2FC_table = Log2FC_table,
                StatPadj = StatPadj)
        } else if (StatPval == "welch") {
            STAT_C1vC2 <- Welch(se,
                SettingsInfo = SettingsInfo, Log2FC_table = Log2FC_table)
        } else if (StatPval == "lmFit") {
            STAT_C1vC2 <- DMA_Stat_limma(se, 
                SettingsInfo = SettingsInfo, Log2FC_table = Log2FC_table,
                StatPadj = StatPadj, CoRe = CoRe, Transform = Transform)
        }
    }

    ############################################################################
    ###############  Add the previous metabolite names back ###############
    DMA_Output <- lapply(STAT_C1vC2, function(df) {
        merged_df <- merge(savedMetaboliteNames, df, by = "Metabolite", all.y = TRUE)
        #merged_df[, -1] %>% 
            ## remove the names we used as part of the function and add back the input names.
        #    dplyr::rename("Metabolite" = 1)
    })


    ############################################################################
    ################  Add the metabolite Metadata if available #################
    if (!is.null(rowData(se))) {
        DMA_Output <- lapply(DMA_Output, function(df){
            merge(df, 
                tibble::rownames_to_column(as.data.frame(rowData(se)), "Metabolite"), 
                by = "Metabolite", all.x = TRUE)
      })
    }

    ############################################################################
    #############  For CoRe=TRUE create summary of Feature_metadata ############
    if (CoRe) {
        df_list_selected <- purrr::map(names(DMA_Output), function(df_name) {
            df <- DMA_Output[[df_name]] ## Extract the dataframe

            ## Extract the dynamic column name
            ## Find the column that starts with "CoRe_"
            core_col <- grep("^CoRe_", names(df), value = TRUE)  
            
            ## Filter only columns where the part after "CoRe_" is in valid_conditions
            core_col <- core_col[str_remove(core_col, "^CoRe_") %in% 
                unique(colData(se)[[SettingsInfo[["Conditions"]]]])]

            ## Select only the relevant columns and return
            df %>%
                select(Metabolite, all_of(core_col))
        })

        ## Merge all dataframes by "Metabolite"
        merged_df <- purrr::reduce(df_list_selected, full_join, 
            by = "Metabolite")
        
        ## It is likely we have duplications that cause .x, .y, .x.x, .y.y, etc. 
        ## to be added to the column names. We only keep one column (.x)
        names(merged_df) <- gsub("\\.x$", "", names(merged_df)) 

        ## Now we remove all other columns with .x.x, .y.y, etc.
        Feature_Metadata <- merged_df %>%
            select(-all_of(grep("\\.[xy]+$", names(merged_df), value = TRUE)))

        ## Add to Metadata file:
        if (!is.null(rowData(se))){
            Feature_Metadata <- merge(
                tibble::rownames_to_column(rowData(se), "Metabolite"), 
                Feature_Metadata , by = "Metabolite", all.x = TRUE)
        }
    }

    ############################################################################
    ###############  Plots ###############
    if (CoRe) {
        x <- "Log2(Distance)"
        VolPlot_SettingsInfo <- c(color = "CoRe")
        VolPlot_SettingsFile <- DMA_Output
    } else {
        x <- "Log2FC"
        VolPlot_SettingsInfo <- NULL
        VolPlot_SettingsFile <- NULL
    }

    volplotList = list()
    se_l <- list()
    for (DF in names(DMA_Output)) { # DF = names(DMA_Output)[2]
        Volplotdata <- DMA_Output[[DF]] |>
            column_to_rownames("Metabolite")

        cD <- data.frame(name = colnames(Volplotdata |> select(-feature)))
        rownames(cD) <- cD$name
        se_volcano <- SummarizedExperiment(
            assays = column_to_rownames(DMA_Output[[DF]], "Metabolite") |> select(-feature), 
            colData = cD,
            rowData = column_to_rownames(DMA_Output[[DF]], "Metabolite"))
        
        se_l[[DF]] <- se_volcano
        #if (CoRe) { ## EDIT: needed
        #    VolPlot_SettingsFile <- DMA_Output[[DF]] %>%
        #        tibble::column_to_rownames("Metabolite")
        #}

        dev.new()
        VolcanoPlot <- invisible(VizVolcano(PlotSettings = "Standard", ### continue from here
            se = se_volcano, ##InputData = tibble::column_to_rownames(Volplotdata, "Metabolite"),
            SettingsInfo = VolPlot_SettingsInfo,
            ##SettingsFile_Metab = VolPlot_SettingsFile,
            y = "p.adj", x = x, PlotName = DF,
            Subtitle = bquote(italic("Differential Metabolite Analysis")),
            SaveAs_Plot = NULL))

        ## Remove special characters and replace spaces with underscores
        DF_save <- gsub("[^A-Za-z0-9._-]", "_", DF)
        volplotList[[DF_save]]<- VolcanoPlot[["Plot_Sized"]][[1]]

        dev.off()
    }
    
    ## assign names to se_l
    names(se_l) <- names(DMA_Output)

    ############################################################################
    ##----- Save and Return
    ## make a list in which we will save the outputs
    DMA_Output_List <- list()
    
    if (PerformShapiro & exists("Shapiro_output")) {
        suppressMessages(suppressWarnings(
            SaveRes(data = Shapiro_output[["DF"]],
                plot = Shapiro_output[["Plot"]][["Distributions"]],
                SaveAs_Table = SaveAs_Table, SaveAs_Plot = SaveAs_Plot,
                FolderPath = SubFolder_S, FileName = "ShapiroTest",
                CoRe = CoRe, PrintPlot = PrintPlot)))
        DMA_Output_List <- list("ShapiroTest" = Shapiro_output)
    }

    if (PerformBartlett & exists("Bartlett_output")) {
        suppressMessages(suppressWarnings(
            SaveRes(data = Bartlett_output[["DF"]],
                plot = Bartlett_output[["Plot"]],
                SaveAs_Table = SaveAs_Table, SaveAs_Plot = SaveAs_Plot,
                FolderPath = SubFolder_B, FileName = "BartlettTest",
                CoRe = CoRe, PrintPlot = PrintPlot)))
        DMA_Output_List <- c(DMA_Output_List, 
            list("BartlettTest" = Bartlett_output))
    }

    if (VST & exists("VST_res")) {
        suppressMessages(suppressWarnings(
            SaveRes(data = VST_res[["DF"]],
                plot = VST_res[["Plot"]], SaveAs_Table = SaveAs_Table,
                SaveAs_Plot = SaveAs_Plot, FolderPath = SubFolder_V,
                FileName = "VST_res", CoRe = CoRe, PrintPlot = PrintPlot)))
        DMA_Output_List <- c(DMA_Output_List, list("VSTres" = Bartlett_output))
    }

    if (CoRe) {
        suppressMessages(suppressWarnings(
            SaveRes(data = list("Feature_Metadata" = Feature_Metadata),
                plot = NULL, SaveAs_Table = SaveAs_Table,
                SaveAs_Plot = NULL, FolderPath = Folder, FileName = "DMA",
                CoRe = CoRe, PrintPlot = PrintPlot)))
        DMA_Output_List <- c(DMA_Output_List, 
            list("Feature_Metadata" = Feature_Metadata))
    }

    suppressMessages(suppressWarnings(
        SaveRes(data = se_l, ##This needs to be a list, also for single comparisons
            plot = volplotList, SaveAs_Table = SaveAs_Table,
            SaveAs_Plot = SaveAs_Plot, FolderPath = Folder, FileName = "DMA",
            CoRe = CoRe, PrintPlot = PrintPlot)))
    DMA_Output_List <- c(DMA_Output_List, 
        list("DMA" = se_l, "VolcanoPlot" = volplotList))

    ## return
    invisible(DMA_Output_List)
}


###############################
### ### ### Log2FC  ### ### ###
###############################

#' This helper function calculates the Log2(FoldChange) or in case of CoRe Log2(Distance).
#'
#' @param InputData DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected.
#' @param SettingsFile_Sample DF which contains metadata information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames.
#' @param SettingsInfo \emph{Optional: } Named vector including the information about the conditions column information on numerator or denominator c(Conditions="ColumnName_SettingsFile", Numerator = "ColumnName_SettingsFile", Denominator  = "ColumnName_SettingsFile"). Denominator and Numerator will specify which comparison(s) will be done (one-vs-one, all-vs-one, all-vs-all), e.g. Denominator=NULL and Numerator =NULL selects all the condition and performs multiple comparison all-vs-all. Log2FC are obtained by dividing the numerator by the denominator, thus positive Log2FC values mean higher expression in the numerator. \strong{Default = c(conditions="Conditions", numerator = NULL, denumerator = NULL)}
#' @param CoRe \emph{Optional: } TRUE or FALSE for whether a Consumption/Release  input is used \strong{default = FALSE}
#' @param Transform \emph{Optional: } If TRUE we expect the data to be not log2 transformed and log2 transformation will be performed within the limma function and Log2FC calculation. If FALSE we expect the data to be log2 transformed as this impacts the Log2FC calculation and limma.\strong{default = TRUE}
#'
#' @return List of DFs named after comparison (e.g. Tumour versus Normal) with Log2FC or Log2(Distance) column and column with feature names
#'
#' @keywords Log2FC, CoRe, Distance
#'
#' @importFrom dplyr select_if filter rename mutate summarise_all
#' @importFrom magrittr %>%
#' @importFrom gtools foldchange2logratio
#' @importFrom tibble rownames_to_column
#'
#' @noRd
#'
Log2FC_fun <-function(se, ##InputData,
    ##SettingsFile_Sample,
    SettingsInfo=c(Conditions = "Conditions", Numerator = NULL, Denominator = NULL),
    CoRe=FALSE,
    Transform=TRUE) {
    
    ## ------------ Create log file ----------- ##
    MetaProViz_Init()

    ## ------------ Assignments ----------- ##
    if (!"Denominator" %in% names(SettingsInfo) & !"Numerator" %in% names(SettingsInfo)) {
        ## all-vs-all: Generate all pairwise combinations
        conditions <- colData(se)[[SettingsInfo[["Conditions"]]]]
        denominator <- unique(conditions)
        numerator <- denominator
        comparisons <- combn(unique(conditions), 2) %>% 
            as.matrix()
        
        ## Settings:
        MultipleComparison <- TRUE
        all_vs_all <- TRUE
    } else if ("Denominator" %in% names(SettingsInfo) & !"Numerator" %in% names(SettingsInfo)) {
        ##all-vs-one: Generate the pairwise combinations
        conditions <- colData(se)[[SettingsInfo[["Conditions"]]]]
        denominator <- SettingsInfo[["Denominator"]]
        numerator <- unique(conditions)
        
        ## Remove denom from num
        numerator <- numerator[!numerator %in% denominator]
        comparisons  <- t(expand.grid(numerator, denominator)) %>% 
            as.data.frame()
        
        ## Settings:
        MultipleComparison <- TRUE
        all_vs_all <- FALSE
    } else if (("Denominator" %in% names(SettingsInfo)) & ("Numerator" %in% names(SettingsInfo))) {
        ## one-vs-one: Generate the comparisons
        denominator <- SettingsInfo[["Denominator"]]
        numerator <- SettingsInfo[["Numerator"]]
        comparisons <- matrix(c(numerator, denominator))
        
        ## Settings:
        MultipleComparison <- FALSE
        all_vs_all <- FALSE
    }

    ## ------------ Check Missingness ------------- ##
    cols_num <- colData(se)[[SettingsInfo[["Conditions"]]]] %in% numerator
    Num <- t(assay(se)[, cols_num])
    cols_denom <- colData(se)[[SettingsInfo[["Conditions"]]]] %in% denominator
    Denom <- t(assay(se)[, cols_denom])

    Num_Miss <- replace(Num, Num == 0, NA)
    Num_Miss <- Num_Miss[, colSums(is.na(Num_Miss)) > 0, drop = FALSE]

    Denom_Miss <- replace(Denom, Denom == 0, NA)
    Denom_Miss <- Denom_Miss[, colSums(is.na(Denom_Miss)) > 0, drop = FALSE]

    if (ncol(Num_Miss) > 0 & ncol(Denom_Miss) == 0){
        Metabolites_Miss <- colnames(Num_Miss)
    } else if (ncol(Num_Miss) == 0 & ncol(Denom_Miss) > 0) {
        Metabolites_Miss <- colnames(Denom_Miss)
    } else if (ncol(Num_Miss) > 0 & ncol(Denom_Miss) > 0) {
        Metabolites_Miss <- c(colnames(Num_Miss), colnames(Denom_Miss))
        Metabolites_Miss <- unique(Metabolites_Miss)
    } else {
        Metabolites_Miss <- c(colnames(Num_Miss), colnames(Denom_Miss))
        Metabolites_Miss <- unique(Metabolites_Miss)
    }

    ## ------------ Denominator/numerator ----------- ##
    ## Denominator and numerator: Define if we compare one_vs_one, one_vs_all or all_vs_all.
    if (!"Denominator" %in% names(SettingsInfo) & !"Numerator" %in% names(SettingsInfo)) {
        MultipleComparison <- TRUE
    } else if ("Denominator" %in% names(SettingsInfo) & !"Numerator" %in% names(SettingsInfo)) {
        MultipleComparison <- TRUE
    } else if ("Denominator" %in% names(SettingsInfo) & "Numerator" %in% names(SettingsInfo)) {
        MultipleComparison <- FALSE
    }

    ############################################################################
    ## ----------------- Log2FC ----------------------------
    Log2FC_table <- list() ## Create an empty list to store results data frames
    
    for (column in seq_len(ncol(comparisons))) {
        cols_num <- colData(se)[[SettingsInfo[["Conditions"]]]] %in% comparisons[1, column]
        C1 <- t(assay(se)[, cols_num]) ## Numerator
        cols_denom <- colData(se)[[SettingsInfo[["Conditions"]]]] %in% comparisons[2, column]
        C2 <- t(assay(se)[, cols_denom]) ## Denominator

        ## ------------  Calculate Log2FC ----------- ##
        ## For C1_Mean and C2_Mean use 0 to obtain values, leading to 
        ## Log2FC=NA if mean = 0 (If one value is NA, the mean will be NA even 
        ## though all other values are available.)
        C1_Zero <- C1
        C1_Zero[is.na(C1_Zero)] <- 0
        Mean_C1 <- C1_Zero %>%
            as.data.frame() |>
            dplyr::summarise_all("mean")

        C2_Zero <- C2
        C2_Zero[is.na(C2_Zero)] <- 0
        Mean_C2 <- C2_Zero %>%
            as.data.frame() |>
            dplyr::summarise_all("mean")

        ## calculate absolute distance between the means. log2 transform and 
        ## add sign (-/+):
        ## CoRe values can be negative and positive, which can does not allow 
        ## us to calculate a Log2FC.
        
        
        ## Mean values could be 0, which can not be used to calculate a 
        ## Log2FC and hence the Log2FC(A versus B)=(log2(A+x)-log2(B+x)) 
        ## for A and/or B being 0, with x being set to 1
        Mean_C1_t <- as.data.frame(t(Mean_C1)) %>%
            tibble::rownames_to_column("Metabolite")
        Mean_C2_t <- as.data.frame(t(Mean_C2)) %>%
            tibble::rownames_to_column("Metabolite")
        Mean_Merge <- merge(Mean_C1_t, Mean_C2_t, by = "Metabolite", 
                all = TRUE) %>%
            dplyr::rename("C1" = 2, "C2" = 3)
        
        ## Deal with NA/0s
        ## create column to enable the check if mean values of 0 are due 
        ## to missing values (NA/0) and not by coincidence
        Mean_Merge$`NA/0` <- Mean_Merge$Metabolite %in% Metabolites_Miss
        
        if (CoRe) { ## EDIT: can this be simplified? Identify steps that are identical between CoRe and !CoRe and try to remove duplicated code
            
            if (any(
                (!Mean_Merge$`NA/0` & Mean_Merge$C1 == 0) | (!Mean_Merge$`NA/0` & Mean_Merge$C2 == 0))) {
                
                Mean_Merge <- Mean_Merge %>%
                    dplyr::mutate(
                        C1 = case_when(
                            ## here we have a "true" 0 value due to 0/NAs in 
                            ## the input data    
                            C2 == 0 & `NA/0` ~ paste(C1), 
                            ## here we have a "true" 0 value due to 0/NAs in 
                            ## the input data
                            C1 == 0 & `NA/0` ~ paste(C1),
                            ## here we have a "false" 0 value that occured at 
                            ## random and not due to 0/NAs in the input data, 
                            ## hence we add the constant +1
                            C2 == 0 & !`NA/0` ~ paste(C1 + 1),
                            ## here we have a "false" 0 value that occured at 
                            ## random and not due to 0/NAs in the input data, 
                            ## hence we add the constant +1
                            C1 == 0 & !`NA/0` ~ paste(C1 + 1),
                            TRUE ~ paste(C1))) %>%
                    dplyr::mutate(
                        C2 = case_when(
                            ## Here we have a "true" 0 value due to 0/NAs in 
                            ## the input data
                            C1 == 0 & `NA/0` ~ paste(C2),
                            ## here we have a "true" 0 value due to 0/NAs in 
                            ## the input data    
                            C2 == 0 & `NA/0` ~ paste(C2),
                            ## here we have a "false" 0 value that occured at 
                            ## random and not due to 0/NAs in the input data, 
                            ## hence we add the constant +1
                            C1 == 0 & !`NA/0` ~ paste(C2 + 1), ## EDIT: is this correct
                            ## here we have a "false" 0 value that occured at 
                            ## random and not due to 0/NAs in the input data, 
                            ## hence we add the constant +1
                            C2 == 0 & !`NA/0` ~ paste(C2 + 1),
                        TRUE ~ paste(C2)))%>%
                    dplyr::mutate(C1 = as.numeric(C1), C2 = as.numeric(C2))

                X <- Mean_Merge %>%
                    subset((!Mean_Merge$`NA/0` & Mean_Merge$C1 ==0) | 
                        (!Mean_Merge$`NA/0` & Mean_Merge$C2==0))
                message("We added +1 to the mean value of metabolite(s) ", 
                    paste0(X$Metabolite, collapse = ", "), 
                    ", since the mean of the replicate values where 0. This was not due to missing values (NA/0).")
            }

            ## Add the distance column:
            Mean_Merge$`Log2(Distance)` <- log2(
                abs(Mean_Merge$C1 - Mean_Merge$C2))

            Mean_Merge <- Mean_Merge %>%
                ## adapt the values to take into account the distance
                dplyr::mutate(`Log2(Distance)` = case_when( ## FEAT: Why not use sign()?
                    ## If C1>C2 the distance stays positive to reflect that C1 > C2
                    C1 > C2 ~ paste(`Log2(Distance)` * + 1),
                    ## If C1<C2 the distance gets a negative sign to reflect that C1 < C2
                    C1 < C2 ~ paste(`Log2(Distance)` * - 1),
                    TRUE ~ 'NA')) %>%
                dplyr::mutate(`Log2(Distance)` = as.numeric(`Log2(Distance)`))

            ##Add additional information:
            temp1 <- Mean_C1
            temp2 <- Mean_C2
            
            ## Add Info of CoRe:
            CoRe_info <- rbind(temp1, temp2, rep(0, length(temp1)))
            for (i in seq_along(temp1)) {
                if (temp1[i] > 0 & temp2[i] > 0) {
                    CoRe_info[3, i] <- "Released"
                } else if (temp1[i] < 0 & temp2[i] < 0) {
                    CoRe_info[3, i] <- "Consumed"
                } else if (temp1[i] > 0 & temp2[i] < 0) {
                    CoRe_info[3, i] <- paste("Released in", 
                        comparisons[1, column] , "and Consumed",
                        comparisons[2, column] , sep = " ")
                } else if (temp1[i] < 0 & temp2[i] > 0) {
                    CoRe_info[3, i] <- paste("Consumed in", 
                        comparisons[1, column] , 
                        " and Released",comparisons[2, column], sep = " ")
                } else {
                    CoRe_info[3, i] <- "No Change"
                }
            }

            CoRe_info <- t(CoRe_info) %>% 
                as.data.frame()
            CoRe_info <- rownames_to_column(CoRe_info, "Metabolite")
            names(CoRe_info)[2] <- paste("Mean",  comparisons[1, column], 
                sep = "_")
            names(CoRe_info)[3] <- paste("Mean",  comparisons[2, column], 
                sep = "_")
            names(CoRe_info)[4] <- "CoRe_specific"

            CoRe_info <- CoRe_info %>%
                dplyr::mutate(CoRe = case_when(
                    CoRe_specific == "Released" ~ 'Released',
                    CoRe_specific == "Consumed" ~ 'Consumed',
                    TRUE ~ 'Released/Consumed')) %>%
                dplyr::mutate(
                        !!paste("CoRe_", comparisons[1, column], sep="") := case_when(
                    CoRe_specific == "Released" ~ 'Released',
                    CoRe_specific == "Consumed" ~ 'Consumed',
                    CoRe_specific == paste("Consumed in", 
                        comparisons[1, column], 
                        " and Released", comparisons[2, column] , sep = " ") ~ 'Consumed',
                    CoRe_specific == paste("Released in", 
                        comparisons[1, column], "and Consumed",
                        comparisons[2, column], sep =" ") ~ 'Released',
                        TRUE ~ 'NA')) %>%
                dplyr::mutate(
                        !!paste("CoRe_", comparisons[2, column], sep = "") := case_when(
                    CoRe_specific == "Released" ~ 'Released',
                    CoRe_specific == "Consumed" ~ 'Consumed',
                    CoRe_specific == paste("Consumed in", 
                        comparisons[1, column], " and Released", 
                        comparisons[2, column] , sep = " ") ~ 'Released',
                    CoRe_specific == paste("Released in", 
                        comparisons[1, column], "and Consumed",
                        comparisons[2, column] , sep = " ") ~ 'Consumed',
                        TRUE ~ 'NA'))

            Log2FC_C1vC2 <-merge(Mean_Merge[, c(1, 5)], 
                CoRe_info[, c(1, 2, 6, 3, 7, 4:5)], by = "Metabolite", 
                all.x = TRUE)

            ## add info on Input
            temp3 <- as.data.frame(t(C1)) %>%
                tibble::rownames_to_column("Metabolite")
            temp4 <- as.data.frame(t(C2)) %>%
                tibble::rownames_to_column("Metabolite")
            temp_3a4 <- merge(temp3, temp4, by = "Metabolite", 
                all = TRUE)
            Log2FC_C1vC2 <- merge(Log2FC_C1vC2, temp_3a4, by = "Metabolite", 
                all.x = TRUE)

            ## Return DFs
            ## make reverse DF
            Log2FC_C2vC1 <- Log2FC_C1vC2
            Log2FC_C2vC1$`Log2(Distance)` <- Log2FC_C2vC1$`Log2(Distance)` * -1

            ## name them
            if (MultipleComparison) {
                logname <- paste(comparisons[1, column], comparisons[2, column],
                    sep = "_vs_")
                logname_reverse <- paste(comparisons[2, column], 
                    comparisons[1, column], sep = "_vs_")

            ## store the data frame in the results list, named after the contrast
            Log2FC_table[[logname]] <- Log2FC_C1vC2
            Log2FC_table[[logname_reverse]] <- Log2FC_C2vC1
            } else {
                Log2FC_table <- Log2FC_C1vC2
            }
        } else { ## !Core
            
            Mean_Merge <- Mean_Merge %>%
                dplyr::mutate(C1_Adapted = case_when(
                    ## here we have a "true" 0 value due to 0/NAs in the input data
                    C2 == 0 & `NA/0` ~ paste(C1),
                    ## here we have a "true" 0 value due to 0/NAs in the input data
                    C1 == 0 & `NA/0` ~ paste(C1),
                    ## here we have a "false" 0 value that occured at random 
                    ## and not due to 0/NAs in the input data, hence we add 
                    ## the constant +1
                    C2 == 0 & !`NA/0` ~ paste(C1 + 1),
                    ## here we have a "false" 0 value that occured at random 
                    ## and not due to 0/NAs in the input data, hence we add 
                    ## the constant +1
                    C1 == 0 & !`NA/0` ~ paste(C1 + 1),
                    TRUE ~ paste(C1))) %>%
                dplyr::mutate(C2_Adapted = case_when(
                    ## here we have a "true" 0 value due to 0/NAs in the input data
                    C1 == 0 & `NA/0` ~ paste(C2),
                    ## here we have a "true" 0 value due to 0/NAs in the input data
                    C2 == 0 & `NA/0` ~ paste(C2),
                    ## here we have a "false" 0 value that occured at random 
                    ## and not due to 0/NAs in the input data, hence we add 
                    ## the constant +1   
                    C1 == 0 & !`NA/0` ~ paste(C2 + 1),
                    ## here we have a "false" 0 value that occured at random 
                    ## and not due to 0/NAs in the input data, hence we add 
                    ## the constant +1
                    C2 == 0 & !`NA/0` ~ paste(C2 + 1),
                    TRUE ~ paste(C2))) %>%
                dplyr::mutate(C1_Adapted = as.numeric(C1_Adapted), 
                    C2_Adapted = as.numeric(C2_Adapted))

            if (any((!Mean_Merge$`NA/0` & Mean_Merge$C1 == 0) | 
                    (!Mean_Merge$`NA/0` & Mean_Merge$C2 == 0))) {
                X <- Mean_Merge %>%
                    subset((!Mean_Merge$`NA/0` & Mean_Merge$C1 == 0) | 
                        (!Mean_Merge$`NA/0` & Mean_Merge$C2 == 0))
                message("We added +1 to the mean value of metabolite(s) ", 
                    paste0(X$Metabolite, collapse = ", "), 
                    ", since the mean of the replicate values where 0. This was not due to missing values (NA/0).")
            }

            ## Calculate the Log2FC
            if (Transform) {
                ## data are not log2 transformed
                ## calculate FoldChange
                Mean_Merge$FC_C1vC2 <- Mean_Merge$C1_Adapted / Mean_Merge$C2_Adapted 
                Mean_Merge$Log2FC <- gtools::foldchange2logratio(
                    Mean_Merge$FC_C1vC2, base = 2)
            }

            if (!Transform) {
                ## data has been log2 transformed and hence we need to 
                ## take this into account when calculating the log2FC
                Mean_Merge$FC_C1vC2 <- "Empty"
                Mean_Merge$Log2FC <- Mean_Merge$C1_Adapted - Mean_Merge$C2_Adapted
            }

            ## Add info on Input:
            temp3 <- as.data.frame(t(C1)) %>%
                tibble::rownames_to_column("Metabolite")
            temp4 <- as.data.frame(t(C2)) %>%
                tibble::rownames_to_column("Metabolite")
            temp_3a4 <- merge(temp3, temp4, by = "Metabolite", all = TRUE)
            Log2FC_C1vC2 <- merge(Mean_Merge[, c(1, 8)], temp_3a4, 
                by = "Metabolite", all.x = TRUE)

            ## Return DFs
            ## Make reverse DF
            Log2FC_C2vC1 <- Log2FC_C1vC2
            Log2FC_C2vC1$Log2FC <- Log2FC_C2vC1$Log2FC * -1

            if (MultipleComparison) {
                logname <- paste(comparisons[1, column], comparisons[2, column],
                    sep = "_vs_")
                logname_reverse <- paste(comparisons[2, column], 
                    comparisons[1, column], sep = "_vs_")

            ## store the data frame in the results list, named after the contrast
            Log2FC_table[[logname]] <- Log2FC_C1vC2
            Log2FC_table[[logname_reverse]] <- Log2FC_C2vC1
            } else {
                Log2FC_table <- Log2FC_C1vC2
            }
        }
    }
  
    ## return
    invisible(Log2FC_table)
}



##########################################################################################
### ### ### DMA helper function: Internal Function to perform single comparison ### ### ###
##########################################################################################

#' This helper function to calculate One-vs-One comparison statistics
#'
#' @param InputData DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected.
#' @param SettingsFile_Sample DF which contains metadata information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames.
#' @param SettingsInfo  Named vector including the information about the conditions column information on numerator or denominator c(Conditions="ColumnName_SettingsFile", Numerator = "ColumnName_SettingsFile", Denominator  = "ColumnName_SettingsFile"). Denominator and Numerator will specify which comparison(s) will be done (here one-vs-one).
#' @param Log2FC_table \emph{Optional: } This is a List of DFs including a column "MetaboliteID" and Log2FC or Log2(Distance). This is the output from MetaProViz:::Log2FC_fun. If NULL, the output statistics will not be added into the Log2FC/Log2(Distance) DFs. \strong{Default = NULL}
#' @param StatPval \emph{Optional: } String which contains an abbreviation of the selected test to calculate p.value. For one-vs-one comparisons choose t.test, wilcox.test, "chisq.test" or "cor.test", \strong{Default = "t.test"}
#' @param StatPadj \emph{Optional: } String which contains an abbreviation of the selected p.adjusted test for p.value correction for multiple Hypothesis testing. Search: ?p.adjust for more methods:"BH", "fdr", "bonferroni", "holm", etc.\strong{Default = "fdr"}
#'
#' @return List of DFs named after comparison (e.g. tumour versus Normal) with p-value, t-value and adjusted p-value column and column with feature names
#'
#' @keywords Statistical testing, p-value, t-value
#'
#' @importFrom stats p.adjust p.adjust.methods
#' @importFrom dplyr select_if filter rename mutate summarise_all
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#'
#' @noRd
#'
DMA_Stat_single <- function(se, ##InputData,
    ##SettingsFile_Sample,
    SettingsInfo,
    Log2FC_table = NULL,
    StatPval = "t.test", ## FEAT: add options
    StatPadj= p.adjust.methods) {
  
    ## ------------ Create log file ----------- ##
    MetaProViz_Init()

    
    ## StatPval <- match.arg(StatPval) ## FEAT: match arguments
    StatPadj <- match.arg(StatPadj)
    
    ## ------------ Check Missingness ------------- ##
    cols_num <- colData(se)[[SettingsInfo[["Conditions"]]]] %in% SettingsInfo[["Numerator"]]
    Num <- t(assay(se)[, cols_num])
    cols_denom <- colData(se)[[SettingsInfo[["Conditions"]]]] %in% SettingsInfo[["Denominator"]]
    Denom <- t(assay(se)[, cols_denom])
    
    Num_Miss <- replace(Num, Num == 0, NA)
    Num_Miss <- Num_Miss[, colSums(is.na(Num_Miss)) > 0, drop = FALSE]

    Denom_Miss <- replace(Denom, Denom == 0, NA)
    Denom_Miss <- Denom_Miss[, colSums(is.na(Denom_Miss) > 0), drop = FALSE]

    if (ncol(Num_Miss) > 0 & ncol(Denom_Miss) == 0) {
        Metabolites_Miss <- colnames(Num_Miss)
    } else if (ncol(Num_Miss) == 0 & ncol(Denom_Miss) > 0) {
        Metabolites_Miss <- colnames(Denom_Miss)
    } else if (ncol(Num_Miss) > 0 & ncol(Denom_Miss) > 0) {
        Metabolites_Miss <- c(colnames(Num_Miss), colnames(Denom_Miss)) |>
            unique()
    } else{
        Metabolites_Miss <- c(colnames(Num_Miss), colnames(Denom_Miss)) |>
            unique()
    }

    ## Comparisons
    comparisons <- matrix(c(SettingsInfo[["Numerator"]], SettingsInfo[["Denominator"]]))

    ## ------------ Perform Hypothesis testing ----------- ##
    for (column in seq_len(ncol(comparisons))) {
        ## Numerator
        cols_num <- colData(se)[[SettingsInfo[["Conditions"]]]] %in% comparisons[1, column]
        C1 <- t(assay(se)[, cols_num])
        ## Denominator
        cols_denom <- colData(se)[[SettingsInfo[["Conditions"]]]] %in% comparisons[2, column]
        C2 <- t(assay(se)[, cols_denom])
    } ## EDIT: What is going on here? Should the for loop be extended to all downstream operations?

    ## For C1 and C2 we use 0, since otherwise we can not perform the statistical testing.
    C1[is.na(C1)] <- 0
    C2[is.na(C2)] <- 0

    #### 1. p.value and test statistics (=t.val)
    T_C1vC2 <- mapply(FUN = StatPval, x = as.data.frame(C2), y = as.data.frame(C1), 
        SIMPLIFY = FALSE)

    VecPVAL_C1vC2 <- c()
    VecTVAL_C1vC2 <- c()
    for (i in 1:length(T_C1vC2)) {
        ## extract p-values and t-values
        p_value <- unlist(T_C1vC2[[i]][3])
        t_value <- unlist(T_C1vC2[[i]])[1]## EDIT: I would do this rather by calling explicitly the name instead of indexing
        VecPVAL_C1vC2[i] <- p_value ## EDITs: why not write directly to VecPVAL_C1vC2 and VecTVAL_C1vC2?
        VecTVAL_C1vC2[i] <- t_value
    }
    Metabolite <- colnames(C2)
    PVal_C1vC2 <- data.frame(Metabolite, p.val = VecPVAL_C1vC2, t.val = VecTVAL_C1vC2)

    ## we set p.val= NA, for metabolites that had 1 or more replicates with 
    ## NA/0 values and remove them prior to p-value adjustment
    PVal_C1vC2$`NA/0` <- PVal_C1vC2$Metabolite %in% Metabolites_Miss
    PVal_C1vC2 <- PVal_C1vC2 %>%
        dplyr::mutate(p.val = case_when(
            `NA/0` ~ NA,
            TRUE ~ paste(VecPVAL_C1vC2)))
    PVal_C1vC2$p.val <- as.numeric(as.character(PVal_C1vC2$p.val))

    #### 2. p.adjusted
    ## Split data for p.value adjustment to exclude NA
    PVal_NA <- PVal_C1vC2[is.na(PVal_C1vC2$p.val), c(1:3)]
    PVal_C1vC2 <- PVal_C1vC2[!is.na(PVal_C1vC2$p.val), c(1:3)]

    ## perform adjustment
    VecPADJ_C1vC2 <- stats::p.adjust(PVal_C1vC2[, 2], method = StatPadj, 
        n = length((PVal_C1vC2[, 2]))) ## p-adjusted
    Metabolite <- PVal_C1vC2[, 1]
    PADJ_C1vC2 <- data.frame(Metabolite, p.adj = VecPADJ_C1vC2)
    STAT_C1vC2 <- merge(PVal_C1vC2, PADJ_C1vC2, by = "Metabolite")

    ## add Metabolites that have p.val=NA back into the DF for completeness.
    if (nrow(PVal_NA) > 0) {
        PVal_NA$p.adj <- NA
        STAT_C1vC2 <- rbind(STAT_C1vC2, PVal_NA)
    }

    ## Add Log2FC
    if (!is.null(Log2FC_table)){
        STAT_C1vC2 <- merge(Log2FC_table, 
            STAT_C1vC2[, c(1:2, 4, 3)], by = "Metabolite")
    }

    ## order the df based on the t-value
    STAT_C1vC2 <- STAT_C1vC2[order(STAT_C1vC2$t.val, decreasing = TRUE), ] 

    ## create list and store results
    results_list <- list()
    results_list[[paste(SettingsInfo[["Numerator"]], "_vs_", SettingsInfo[["Denominator"]])]] <- STAT_C1vC2

    ## return
    invisible(results_list)
}


################################################################
### ### ### AOV: Internal Function to perform Anova  ### ### ###
################################################################

#' This helper function to calculate One-vs-All or All-vs-All comparison statistics
#'
#' @param InputData DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected.
#' @param SettingsFile_Sample DF which contains metadata information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames.
#' @param SettingsInfo  \emph{Optional: } Named vector including the information about the conditions column information on numerator or denominator c(Conditions="ColumnName_SettingsFile", Numerator = "ColumnName_SettingsFile", Denominator  = "ColumnName_SettingsFile"). Denominator and Numerator will specify which comparison(s) will be done (Here all-vs-one, all-vs-all), e.g. Denominator=NULL and Numerator =NULL selects all the condition and performs multiple comparison all-vs-all. \strong{Default = c(conditions="Conditions", numerator = NULL, denumerator = NULL)}
#' @param Log2FC_table \emph{Optional: } This is a List of DFs including a column "MetaboliteID" and Log2FC or Log2(Distance). This is the output from MetaProViz:::Log2FC_fun. If NULL, the output statistics will not be added into the Log2FC/Log2(Distance) DFs. \strong{Default = NULL}
#'
#' @return List of DFs named after comparison (e.g. tumour versus Normal) with p-value, t-value and adjusted p-value column and column with feature names
#'
#' @keywords Statistical testing, p-value, t-value
#'
#' @importFrom stats aov TukeyHSD
#' @importFrom dplyr rename
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#'
#' @noRd
#'
AOV <-function(se, ##InputData,
    ##SettingsFile_Sample,
    SettingsInfo = c(Conditions = "Conditions", Numerator = NULL, Denominator = NULL),
    Log2FC_table = NULL){

    ## ------------ Create log file ----------- ##
    MetaProViz_Init()

    ## ------------ Denominator/numerator ----------- ##
    ## Denominator and numerator: Define if we compare one_vs_one, 
    ## one_vs_all or all_vs_all.
    if (!"Denominator" %in% names(SettingsInfo)  & !"Numerator" %in% names(SettingsInfo)) { ## EDIT: line 1001-1028 is replicated across several functions and should be written as a function
        
        ## all-vs-all: Generate all pairwise combinations
        conditions <- colData(se)[[SettingsInfo[["Conditions"]]]]
        denominator <-unique(conditions)
        numerator <- unique(conditions)
        comparisons <- combn(unique(conditions), 2) %>% 
            as.matrix()
        
        ## Settings:
        MultipleComparison <- TRUE
        all_vs_all <- TRUE
    } else if ("Denominator" %in% names(SettingsInfo)  & !"Numerator" %in% names(SettingsInfo)) {
        
        ## all-vs-one: Generate the pairwise combinations
        conditions <- colData(se)[[SettingsInfo[["Conditions"]]]]
        denominator <- SettingsInfo[["Denominator"]]
        numerator <- unique(conditions)
        
        ## remove denom from num
        numerator <- numerator[!numerator %in% denominator]
        comparisons  <- t(expand.grid(numerator, denominator)) %>% 
            as.data.frame()
        
        ## Settings:
        MultipleComparison <- TRUE
        all_vs_all <- FALSE
    }

    #############################################################################################
    ## 1. Anova p.val
    aov.res <- apply(assay(se), 1, function(x) stats::aov(x ~ conditions))

    ## 2. Tukey test p.adj
    posthoc.res <- lapply(aov.res, stats::TukeyHSD, conf.level = 0.95)
    Tukey_res <- do.call("rbind", 
        lapply(posthoc.res, function(x) x[1][[1]][, "p adj"])) %>% 
        as.data.frame()

    comps <-paste(comparisons[1, ], comparisons[2, ], sep = "-") ## normal
    opp_comps <- paste(comparisons[2, ], comparisons[1, ], sep = "-")

    ## if opposite comparisons is true
    if (sum(opp_comps %in%  colnames(Tukey_res)) > 0) {
        for (comp in seq_along(opp_comps)) {
            colnames(Tukey_res)[colnames(Tukey_res) %in% opp_comps[comp]] <- comps[comp]
        }
    }

    ## 3. t.val
    Tukey_res_diff <- do.call("rbind", 
        lapply(posthoc.res, function(x) x[1][[1]][, "diff"])) %>% ## EDIT: could be "recycled" from l1043?
        as.data.frame()

    ## if oposite comparisons is true
    if (sum(opp_comps %in% colnames(Tukey_res_diff)) > 0) {
        for (comp in 1: length(opp_comps)){
            colnames(Tukey_res_diff)[colnames(Tukey_res_diff) %in% opp_comps[comp]] <-  comps[comp]
        }
    }

    ## Make output DFs:
    Pval_table <- Tukey_res |>
        tibble::rownames_to_column("Metabolite")
    Tval_table <- tibble::rownames_to_column(Tukey_res_diff, "Metabolite")
    
    ## here we need to adapt for one_vs_all or all_vs_all
    common_col_names <- setdiff(names(Tukey_res_diff), "row.names")

    results_list <- list()
    for(col_name in common_col_names){
        ## create a new data frame by merging the two data frames
        merged_df <- merge(Pval_table[, c("Metabolite", col_name)], 
                Tval_table[,c("Metabolite", col_name)], 
                by = "Metabolite", all = TRUE) %>%
            dplyr::rename("p.adj" = 2, "t.val" = 3)

        ## we need to add _vs_ into the comparison col_name
        pattern <- paste(conditions, collapse = "|")
        conditions_present <- unique(
            unlist(regmatches(col_name, gregexpr(pattern, col_name))))
        modified_col_name <- paste(conditions_present[1], "vs", 
            conditions_present[2], sep = "_")

        ## add the new data frame to the list with the column name as the 
        ## list element name
        results_list[[modified_col_name]] <- merged_df
    }

    ## merge the data frames in list1 and list2 based on the "Metabolite" column
    if (!is.null(Log2FC_table)) {
        list_names <- names(results_list)

        merged_list <- list()
        for (name in list_names) {
            ## check if the data frames exist in both lists
            if (name %in% names(results_list) && name %in% names(Log2FC_table)) {
                merged_df <- merge(results_list[[name]], 
                    Log2FC_table[[name]], by = "Metabolite", all = TRUE)
        
                ## reorder the columns
                merged_df <- merged_df[,c(1, 4, 2:3, 5:ncol(merged_df))]
                merged_list[[name]] <- merged_df
            }
        }
    } else {
        merged_list <- results_list
    }

    ## make sure the right comparisons are returned:
    if (all_vs_all) {
        STAT_C1vC2 <- merged_list
    } else if (!all_vs_all) { ## EDIT: the additional else is not needed, just use else
        ## remove the comparisons that are not needed:
        modified_df_list <- list()
        for (df_name in names(merged_list)) {
            if (endsWith(df_name, SettingsInfo[["Denominator"]])) {
                modified_df_list[[df_name]] <- merged_list[[df_name]]
            }
        }
        STAT_C1vC2 <- modified_df_list
    }

    ## return
    invisible(STAT_C1vC2)
}

###########################################################################
### ### ### Kruskal: Internal Function to perform Kruskal test  ### ### ###
###########################################################################

#' This helper function to calculate One-vs-All or All-vs-All comparison statistics
#'
#' @param InputData DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected.
#' @param SettingsFile_Sample DF which contains metadata information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames.
#' @param SettingsInfo  \emph{Optional: } Named vector including the information about the conditions column information on numerator or denominator c(Conditions="ColumnName_SettingsFile", Numerator = "ColumnName_SettingsFile", Denominator  = "ColumnName_SettingsFile"). Denominator and Numerator will specify which comparison(s) will be done (Here all-vs-one, all-vs-all), e.g. Denominator=NULL and Numerator =NULL selects all the condition and performs multiple comparison all-vs-all. \strong{Default = c(conditions="Conditions", numerator = NULL, denumerator = NULL)}
#' @param Log2FC_table \emph{Optional: } This is a List of DFs including a column "MetaboliteID" and Log2FC or Log2(Distance). This is the output from MetaProViz:::Log2FC_fun. If NULL, the output statistics will not be added into the Log2FC/Log2(Distance) DFs. \strong{Default = NULL}
#' @param StatPadj \emph{Optional: } String which contains an abbreviation of the selected p.adjusted test for p.value correction for multiple Hypothesis testing. Search: ?p.adjust for more methods:"BH", "fdr", "bonferroni", "holm", etc.\strong{Default = "fdr"}
#'
#' @return List of DFs named after comparison (e.g. tumour versus Normal) with p-value, t-value and adjusted p-value column and column with feature names
#'
#' @keywords Statistical testing, p-value, t-value
#'
#' @importFrom stats kruskal.test p.adjust.methods
#' @importFrom rstatix dunn_test
#' @importFrom dplyr rename mutate_all mutate select
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column column_to_rownames
#'
#' @noRd
#'
Kruskal <- function(se, ##InputData,
    ##SettingsFile_Sample,
    SettingsInfo = c(Conditions = "Conditions", Numerator = NULL, Denominator = NULL),
    Log2FC_table = NULL,
    StatPadj = p.adjust.methods) {

    ## ------------ Create log file ----------- ##
    MetaProViz_Init()

    ## match arguments
    StatPadj <- match.arg(StatPadj)
    
    ## ------------ Denominator/numerator ----------- ##
    ## Denominator and numerator: Define if we compare one_vs_one, 
    ## one_vs_all or all_vs_all.
    if (!"Denominator" %in% names(SettingsInfo) & !"Numerator" %in% names(SettingsInfo)) { ## EDIT: this is replicated and should be written as a function
        
        ## all-vs-all: Generate all pairwise combinations
        conditions <- colData(se)[[SettingsInfo[["Conditions"]]]]
        denominator <- unique(conditions)
        numerator <- unique(conditions)
        comparisons <- combn(unique(conditions), 2) %>%
            as.matrix()
        
        ## Settings:
        MultipleComparison <- TRUE
        all_vs_all <- TRUE
        
    } else if ("Denominator" %in% names(SettingsInfo) & !"Numerator" %in% names(SettingsInfo)) {
        ## all-vs-one: Generate the pairwise combinations
        conditions <- colData(se)[[SettingsInfo[["Conditions"]]]]
        denominator <- SettingsInfo[["Denominator"]]
        numerator <- unique(conditions)
        
        ## Remove denom from num
        numerator <- numerator[!numerator %in% denominator]
        comparisons  <- t(expand.grid(numerator, denominator)) %>% 
            as.data.frame()
    
        ## Settings:
        MultipleComparison <- TRUE
        all_vs_all <- FALSE
    }

    #############################################################################################
    ## Kruskal test (p.val)
    aov.res <- apply(assay(se), 2, function(x) stats::kruskal.test(x ~ conditions))
    anova_res <- do.call("rbind", lapply(aov.res, function(x) x["p.value"]))
    anova_res <- as.matrix(dplyr::mutate_all(as.data.frame(anova_res), 
        function(x) as.numeric(as.character(x))))
    colnames(anova_res) = c("Kruskal_p.val")

    ## Dunn test (p.adj)
    Dunndata <- assay(se) %>%
        dplyr::mutate(conditions = conditions) %>%
        dplyr::select(conditions, everything()) %>%
        as.data.frame()

    ## applying a loop to obtain p.adj and t.val:
    Dunn_Pres <- data.frame(
        comparisons = paste(comparisons[1, ],comparisons[2, ], sep = "_vs_" ))
    Dunn_Tres <- Dunn_Pres
    
    for (col in seq_len(ncol(Dunndata))[-1]) {
        data <- Dunndata[, c(1, col)]
        colnames(data)[2] <- gsub("^\\d+", "", colnames(data)[2])

        ## If a metabolite starts with number remove it
        formula <- as.formula(paste(colnames(data)[2], "~ conditions"))
        posthoc.res= rstatix::dunn_test(data, formula, p.adjust.method = StatPadj)

        pres <- data.frame(
            comparisons = c(paste(posthoc.res$group1, posthoc.res$group2, sep = "_vs_" ), 
                paste(posthoc.res$group2, posthoc.res$group1, sep = "_vs_" )))
        pres[[colnames(Dunndata)[col]]] <-  c(posthoc.res$p.adj, posthoc.res$p.adj )
        
        ## take only the comparisons selected
        pres <- pres[pres$comparisons %in% Dunn_Pres$comparisons, ] 
        Dunn_Pres <- merge(Dunn_Pres, pres, by="comparisons")

        tres <- data.frame(
            comparisons <- c(paste(posthoc.res$group1, posthoc.res$group2, sep = "_vs_" ), 
                paste(posthoc.res$group2, posthoc.res$group1, sep = "_vs_" )))
        tres[[colnames(Dunndata)[col]]] <- c(posthoc.res$statistic, -posthoc.res$statistic)
        
        ## take only the comparisons selected
        tres <- tres[tres$comparisons %in% Dunn_Pres$comparisons, ] 
        Dunn_Tres <- merge(Dunn_Tres, tres, by = "comparisons")
    }

    ## Make output DFs:
    Dunn_Pres <- tibble::column_to_rownames(Dunn_Pres, "comparisons") %>% 
        t() %>% 
        as.data.frame()
    Pval_table <- as.matrix(dplyr::mutate_all(as.data.frame(Dunn_Pres), 
            function(x) as.numeric(as.character(x)))) %>%
        as.data.frame()%>%
        tibble::rownames_to_column("Metabolite")

    Dunn_Tres <- tibble::column_to_rownames(Dunn_Tres, "comparisons") %>% 
        t() %>% 
        as.data.frame()
    Tval_table <- as.matrix(dplyr::mutate_all(as.data.frame(Dunn_Tres), 
            function(x) as.numeric(as.character(x)))) %>%
        as.data.frame() %>%
        tibble::rownames_to_column("Metabolite")

    common_col_names <- setdiff(names(Dunn_Pres), "row.names")

    results_list <- list()
    for (col_name in common_col_names) {
        ## create a new data frame by merging the two data frames
        merged_df <- merge(Pval_table[,c("Metabolite", col_name)], 
            Tval_table[,c("Metabolite",col_name)], 
            by = "Metabolite", all = TRUE) %>%
            dplyr::rename("p.adj" = 2, "t.val" = 3)
        
        ## add the new data frame to the list with the column name as the list 
        ## element name
        results_list[[col_name]] <- merged_df
    }

    ## merge the data frames in list1 and list2 based on the "Metabolite" column ## EDIT: this is replicated and could be written as a funciton
    if (!is.null(Log2FC_table)) {
        merged_list <- list()
        for (name in common_col_names) {
            ## check if the data frames exist in both lists
            if (name %in% names(results_list) && name %in% names(Log2FC_table)) {
                merged_df <- merge(results_list[[name]], Log2FC_table[[name]], 
                    by = "Metabolite", all = TRUE)
                ## reorder the columns
                merged_df <- merged_df[,c(1, 4, 2:3, 5:ncol(merged_df))] ## EDIT: could be written directly to merged_list[[name]]
                merged_list[[name]] <- merged_df
            }
        }
        STAT_C1vC2 <- merged_list
    } else {
        STAT_C1vC2 <- results_list
    }

    ## return
    invisible(STAT_C1vC2)
}



#############################################################################################
### ### ### Welch: Internal Function to perform anova for unequal variance groups ### ### ###
#############################################################################################

#' This helper function to calculate One-vs-All or All-vs-All comparison statistics
#'
#' @param InputData DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected.
#' @param SettingsFile_Sample DF which contains metadata information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames.
#' @param SettingsInfo  \emph{Optional: } Named vector including the information about the conditions column information on numerator or denominator c(Conditions="ColumnName_SettingsFile", Numerator = "ColumnName_SettingsFile", Denominator  = "ColumnName_SettingsFile"). Denominator and Numerator will specify which comparison(s) will be done (Here all-vs-one, all-vs-all), e.g. Denominator=NULL and Numerator =NULL selects all the condition and performs multiple comparison all-vs-all. \strong{Default = c(conditions="Conditions", numerator = NULL, denumerator = NULL)}
#' @param Log2FC_table \emph{Optional: } This is a List of DFs including a column "MetaboliteID" and Log2FC or Log2(Distance). This is the output from MetaProViz:::Log2FC_fun. If NULL, the output statistics will not be added into the Log2FC/Log2(Distance) DFs. \strong{Default = NULL}
#'
#' @return List of DFs named after comparison (e.g. tumour versus Normal) with p-value, t-value and adjusted p-value column and column with feature names
#'
#' @keywords Statistical testing, p-value, t-value
#'
#' @importFrom rstatix games_howell_test
#' @importFrom dplyr rename
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#'
#' @noRd
#'
Welch <-function(##InputData,
    se,
    ##SettingsFile_Sample,
    SettingsInfo = c(Conditions = "Conditions", Numerator = NULL, Denominator  = NULL),
    Log2FC_table = NULL) {
    
    ## ------------ Create log file ----------- ##
    MetaProViz_Init()

    ## ------------ Denominator/numerator ----------- ##
    ## Denominator and numerator: Define if we compare one_vs_one, 
    ## one_vs_all or all_vs_all.
    if (!"Denominator" %in% names(SettingsInfo) & !"Numerator" %in% names(SettingsInfo)) { ## EDIT: this is replicated and should be written as a function
        
        ## all-vs-all: Generate all pairwise combinations
        conditions <- colData(se)[[SettingsInfo[["Conditions"]]]]
        denominator <- unique(conditions)
        numerator <- unique(conditions)
        comparisons <- combn(unique(conditions), 2) %>%
            as.matrix()
        
        ## Settings:
        MultipleComparison <- TRUE
        all_vs_all <- TRUE
        
    } else if ("Denominator" %in% names(SettingsInfo) & !"Numerator" %in% names(SettingsInfo)) {
        ## all-vs-one: Generate the pairwise combinations
        conditions <- colData(se)[[SettingsInfo[["Conditions"]]]]
        denominator <- SettingsInfo[["Denominator"]]
        numerator <-unique(conditions)
        
        ## Remove denom from num
        numerator <- numerator[!numerator %in% denominator]
        comparisons  <- t(expand.grid(numerator, denominator)) %>% 
            as.data.frame()
        
        ## Settings:
        MultipleComparison <- TRUE
        all_vs_all <- FALSE
    }

    ############################################################################
    ## 1. Welch's ANOVA using oneway.test is not used by the Games post.hoc 
    ## function
    #aov.res= apply(Input_data,2,function(x) oneway.test(x~conditions))
    games_data <- merge(colData(se), t(assay(se)), by = 0) %>%
        dplyr::rename("conditions" = SettingsInfo[["Conditions"]])
    games_data$conditions <- conditions
    posthoc.res.list <- list()

    ## 2. Games post hoc test
    for (row in rownames(se)) { # col = names(Input_data)[1]
        posthoc.res <- rstatix::games_howell_test(data = games_data, ## EDIT: why use :: ?? if you define those in the NAMESPACE, also applies to other instances in the package
                detailed =TRUE, 
                formula = as.formula(paste0(`row`, " ~ ", "conditions"))) %>% 
            as.data.frame()

        result.df <- rbind(
            data.frame(p.adj = posthoc.res[, "p.adj"],
                t.val = posthoc.res[, "statistic"],
                row.names = paste(posthoc.res[["group1"]], posthoc.res[["group2"]], sep = "-")),
            data.frame(p.adj = posthoc.res[,"p.adj"],
                t.val = -posthoc.res[,"statistic"],
                row.names = paste(posthoc.res[["group2"]], posthoc.res[["group1"]], sep = "-")))
        posthoc.res.list[[row]] <- result.df
    }
    
    Games_Pres <- do.call("rbind", lapply(posthoc.res.list, 
            function(x) x[, "p.adj"])) %>% 
        as.data.frame()
    colnames(Games_Pres) <- rownames(posthoc.res.list[[1]])
    comps <-   paste(comparisons[1, ], comparisons[2, ], sep = "-")# normal
    Games_Pres <- Games_Pres[, colnames(Games_Pres) %in% comps] %>% 
        tibble::rownames_to_column("Metabolite")
  
    ## in case of p.adj = 0 we change it to 10^-6
    Games_Pres[Games_Pres == 0] <- 0.000001

    ## 3. t.val
    Games_Tres <- do.call("rbind", 
            lapply(posthoc.res.list, function(x) x[, "t.val"])) %>% 
        as.data.frame()
    colnames(Games_Tres) <- rownames(posthoc.res.list[[1]])
    Games_Tres <- Games_Tres[, colnames(Games_Tres) %in% comps] %>% 
        tibble::rownames_to_column("Metabolite")

    results_list <- list()
    for (col_name in colnames(Games_Pres)) {
        ## create a new data frame by merging the two data frames
        merged_df <- merge(Games_Pres[, c("Metabolite", col_name)], 
                Games_Tres[, c("Metabolite", col_name)], 
                by = "Metabolite", all = TRUE) %>%
            dplyr::rename("p.adj" = 2, "t.val" = 3)

        ## add _vs_ into the comparison col_name
        pattern <- paste(conditions, collapse = "|")
        conditions_present <- unique(
            unlist(regmatches(col_name, gregexpr(pattern, col_name))))
        modified_col_name <- paste(conditions_present[1], "vs", 
            conditions_present[2], sep = "_")

        ## Add the new data frame to the list with the column name as the list 
        ##element name
        results_list[[modified_col_name]] <- merged_df
    }

    ## merge the data frames in list1 and list2 based on the "Metabolite" column ## EDIT: this is replicated and could be written as a funciton
    if (!is.null(Log2FC_table)) {
        list_names <-  names(results_list)

        merged_list <- list()
        for (name in list_names) {
            ## check if the data frames exist in both lists
            if (name %in% names(results_list) && name %in% names(Log2FC_table)) {
                merged_df <- merge(results_list[[name]], Log2FC_table[[name]], 
                    by = "Metabolite", all = TRUE)
                
                ## reorder the columns
                merged_df <- merged_df[, c(1, 4, 2:3, 5:ncol(merged_df))]
                merged_list[[name]] <- merged_df
             }
        }
        STAT_C1vC2 <- merged_list
    } else {
        STAT_C1vC2 <- results_list
    }

    ## return
    invisible(STAT_C1vC2)
}


##########################################################################################
### ### ### DMA helper function: Internal Function to perform limma ### ### ###
##########################################################################################

#' This helper function to calculate One-vs-One, One-vs-All or All-vs-All comparison statistics
#'
#' @param InputData DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected.
#' @param SettingsFile_Sample DF which contains metadata information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames.
#' @param SettingsInfo  \emph{Optional: } Named vector including the information about the conditions column information on numerator or denominator c(Conditions="ColumnName_SettingsFile", Numerator = "ColumnName_SettingsFile", Denominator  = "ColumnName_SettingsFile"). Denominator and Numerator will specify which comparison(s) will be done (one-vs-all, all-vs-one, all-vs-all), e.g. Denominator=NULL and Numerator =NULL selects all the condition and performs multiple comparison all-vs-all. \strong{Default = c(conditions="Conditions", numerator = NULL, denumerator = NULL)}
#' @param Log2FC_table \emph{Optional: } This is a List of DFs including a column "MetaboliteID" and Log2FC or Log2(Distance). This is the output from MetaProViz:::Log2FC_fun. If NULL, the output statistics will not be added into the Log2FC/Log2(Distance) DFs. \strong{Default = NULL}
#' @param StatPadj \emph{Optional: } String which contains an abbreviation of the selected p.adjusted test for p.value correction for multiple Hypothesis testing. Search: ?p.adjust for more methods:"BH", "fdr", "bonferroni", "holm", etc.\strong{Default = "fdr"}
#' @param CoRe \emph{Optional: } TRUE or FALSE for whether a Consumption/Release  input is used. \strong{Default = FALSE}
#' @param Transform TRUE or FALSE. If TRUE we expect the data to be not log2 transformed and log2 transformation will be performed within the limma function and Log2FC calculation. If FALSE we expect the data to be log2 transformed as this impacts the Log2FC calculation and limma. \strong{Default= TRUE}
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
#'
#' @noRd
#'
DMA_Stat_limma <- function(se, ##InputData,
    ##SettingsFile_Sample,
    SettingsInfo = c(Conditions = "Conditions", Numerator = NULL, Denominator  = NULL),
    Log2FC_table = NULL,
    StatPadj = p.adjust.methods,
    CoRe = FALSE,
    Transform = TRUE){

    ## ------------ Create log file ----------- ##
    MetaProViz_Init()
    
    StatPadj <- match.arg(StatPadj)

    ## ------------ Denominator/numerator ----------- ##
    ## Denominator and numerator: Define if we compare one_vs_one, one_vs_all 
    ## or all_vs_all.
    if (!"Denominator" %in% names(SettingsInfo) & !"Numerator" %in% names(SettingsInfo)) { ## EDIT: can be simplified, define MultiComarison = TRUE and only when '"Denominator" %in% names(SettingsInfo)' adjust the values
        MultipleComparison = TRUE
        all_vs_all = TRUE
    } else if ("Denominator" %in% names(SettingsInfo) & !"Numerator" %in% names(SettingsInfo)) {
        MultipleComparison = TRUE
        all_vs_all = FALSE
    } else if ("Denominator" %in% names(SettingsInfo) & "Numerator" %in% names(SettingsInfo)) {
        MultipleComparison = FALSE
        all_vs_all = FALSE
    }

    ## ensure that Input_data is ordered by conditions and sample names are 
    ## the same as in Input_SettingsFile_Sample:
    targets <- colData(se) |>
        as.data.frame() |>
        tibble::rownames_to_column("sample")
    targets <- targets[, c("sample", SettingsInfo[["Conditions"]])] %>%
        dplyr::rename("condition" = 2) |># %>%
        ### order the column "sample" alphabetically
        #dplyr::arrange(sample)
        ## make appropriate condition names accepted by limma
        mutate(condition_limma_compatible = make.names(condition))

    ## create Limma_input
    Limma_input <- assay(se) %>%
        t() |>
        as.data.frame() |>
        tibble::rownames_to_column("sample")
    
    ## update Limma_input if MultipleComparison == FALSE
    if (!MultipleComparison) {
        ## subset the data:
        targets <- targets %>%
            subset(condition == SettingsInfo[["Numerator"]] | condition == SettingsInfo[["Denominator"]]) #%>%
            ## order the column "sample" alphabetically ## EDIT: is this actually needed after the arrange(smaple) step above?
            #dplyr::arrange(sample)

        #Limma_input <- assay(se) %>%
        #    tibble::rownames_to_column("sample") 
        Limma_input <- merge(targets[, 1:2],  Limma_input, 
            by = "sample", all.x = TRUE) ## EDIT: more robust to use column names instead of column indices
        Limma_input <- Limma_input[, -2] %>%
            ## Order the column "sample" alphabetically
            arrange(sample)
    }# else if (MultipleComparison) { ## EDIT: the other if statement is not needed, only else
    #    Limma_input <- assay(se) %>%
    #        t() |>
    #        as.data.frame() |>
    #        tibble::rownames_to_column("sample") %>%
    #        dplyr::arrange(sample)#Order the column "sample" alphabetically
    #}

    ## check if the order of the "sample" column is the same in both data frames
    #if (!identical(targets$sample, Limma_input$sample)) {
    #    stop("The order of the 'sample' column is different in both data frames. Please make sure that Input_SettingsFile_Sample and Input_data contain the same rownames and sample numbers.")
    #}

    targets_limma <- targets |>
        dplyr::select(-c("condition")) %>%
        dplyr::rename("condition" = "condition_limma_compatible")

    ## we need to transpose the df to run limma. Also, if the data is not 
    ## log2 transformed, we will not calculate the Log2FC as limma just 
    ## substracts one condition from the other
    Limma_input <- tibble::column_to_rownames(Limma_input, "sample") |>
        t() |>
        as.data.frame()

    if (Transform) {
        ## communicate the log2 transformation --> how does limma deals with NA when calculating the change?
        Limma_input <- log2(Limma_input) 
    }

    #### ------Run limma:
    ####  Make design matrix:
    ## all versus all
    fcond <- as.factor(targets_limma$condition)

    ## create the design matrix
    design <- model.matrix(~ 0 + fcond)
    ## give meaningful column names to the design matrix
    colnames(design) <- levels(fcond)

    ## fit the linear model
    fit <- limma::lmFit(Limma_input, design)

    ## make contrast matrix:
    if (all_vs_all & MultipleComparison) {
        ## get unique conditions
        unique_conditions <- levels(fcond)

        ## create an empty contrast matrix
        num_conditions <- length(unique_conditions)
        num_comparisons <- num_conditions * (num_conditions - 1) / 2
        cont.matrix <- matrix(0, nrow = num_comparisons, ncol = num_conditions)

        ## initialize an index for the column in the contrast matrix
        i <- 1

        ## initialize column and row names
        colnames(cont.matrix) <- unique_conditions
        rownames(cont.matrix) <- character(num_comparisons)

        ## Loop through all pairwise combinations of unique conditions
        for (condition1 in seq_len(num_conditions - 1)) {
            for (condition2 in (condition1 + 1):num_conditions) {
        
                ## create the pairwise comparison vector
                comparison <- rep(0, num_conditions)
    
                comparison[condition2] <- -1
                comparison[condition1] <- 1
                ## add the comparison vector to the contrast matrix
                cont.matrix[i, ] <- comparison
                ## set row name
                rownames(cont.matrix)[i] <- paste(unique_conditions[condition1], 
                    "_vs_", unique_conditions[condition2], sep = "")
                i <- i + 1
            }
        }
        cont.matrix <- t(cont.matrix)
    } else if (!all_vs_all & MultipleComparison) {
        ## get unique conditions
        unique_conditions <- levels(fcond)
        denominator  <- make.names(SettingsInfo[["Denominator"]])
    
        ## create an empty contrast matrix
        num_conditions <- length(unique_conditions)
        num_comparisons <- num_conditions - 1
        cont.matrix <- matrix(0, nrow = num_comparisons, ncol = num_conditions)
        
        ## initialize an index for the column in the contrast matrix
        i <- 1
    
        ## initialize column and row names
        colnames(cont.matrix) <- unique_conditions
        rownames(cont.matrix) <- character(num_comparisons)
    
        ## loop through all pairwise combinations of unique conditions
        for(condition in seq_len(num_conditions)[-1]) {
            
            ## create the pairwise comparison vector
            comparison <- rep(0, num_conditions)
            
            if (unique_conditions[1] == make.names(SettingsInfo[["Denominator"]])) {
                comparison[1] <- -1
                comparison[condition] <- 1
                
                ## add the comparison vector to the contrast matrix
                cont.matrix[i, ] <- comparison
                
                ## set row name
                rownames(cont.matrix)[i] <- paste(unique_conditions[condition], 
                    "_vs_", unique_conditions[1], sep = "")
            } else {
                
                comparison[1] <- 1
                comparison[condition] <- -1
                
                ## add the comparison vector to the contrast matrix
                cont.matrix[i, ] <- comparison
                
                ## set row name
                rownames(cont.matrix)[i] <- paste(unique_conditions[1], 
                    "_vs_", unique_conditions[condition], sep = "")

            }
            i <- i + 1
        }
    
        cont.matrix <- t(cont.matrix)
    } else if (!all_vs_all & !MultipleComparison) {
    
        Name_Comp <- paste(make.names(SettingsInfo[["Numerator"]]), "-", make.names(SettingsInfo[["Denominator"]]), sep = "")
        cont.matrix <- as.data.frame(
            limma::makeContrasts(contrasts = Name_Comp, levels = colnames(design))) %>%
            dplyr::rename(!!paste(make.names(SettingsInfo[["Numerator"]]), 
                "_vs_", make.names(SettingsInfo[["Denominator"]]), sep = "") := 1)
        cont.matrix <- as.matrix(cont.matrix)
    }

    ## fit the linear model with contrasts
    fit2 <- limma::contrasts.fit(fit, cont.matrix)
  
    ## Perform empirical Bayes moderation
    fit2 <- limma::eBayes(fit2)

    #### ------Extract results:
    ## get all contrast names
    contrast_names <- colnames(fit2$coefficients)

    ## create an empty list to store results data frames
    results_list <- list()
    for (contrast_name in contrast_names) {
        ## extract results for the current contrast
        res.t <- limma::topTable(fit2, coef = contrast_name, n = Inf, ## EDIT: . should be avoided in object names, better use _ instead
                sort.by = "n", adjust.method = StatPadj) %>% # coef= the comparison the test is done for!
            dplyr::rename("Log2FC" = 1, "t.val" = 3, "p.val" = 4, "p.adj" = 5)

        res.t <- res.t %>%
            tibble::rownames_to_column("Metabolite")

        ## store the data frame in the results list, named after the contrast
        results_list[[contrast_name]] <- res.t
    }

    ## make the name_match_df
    name_match_df <- as.data.frame(names(results_list)) %>%
        tidyr::separate("names(results_list)", into = c("a", "b"), 
            sep = "_vs_", remove = FALSE)

    name_match_df <- merge(name_match_df, targets[, -c(1)] , by.x = "a", 
            by.y = "condition_limma_compatible", all.x = TRUE) %>%
        dplyr::rename("Condition1" = 4)
    name_match_df <- merge(name_match_df, targets[, -c(1)] , by.x = "b", ## EDIT: why assign here again to name_match_df and not do everything in one go
            by.y = "condition_limma_compatible", all.x = TRUE) %>%
        dplyr::rename("Condition2" = 5) %>%
        tidyr::unite("New", "Condition1", "Condition2", sep = "_vs_", 
            remove = FALSE)

    name_match_df<- name_match_df[, c(3,4)] %>%
        dplyr::distinct(New, .keep_all = TRUE)

    results_list_new <- list()
    
    ## match the lists using name_match_df
    for (i in seq_len(nrow(name_match_df))) {
        old_name <- name_match_df$`names(results_list)`[i]
        new_name <- name_match_df$New[i]
        results_list_new[[new_name]] <- results_list[[old_name]]
    }

    if (!is.null(Log2FC_table)) {
        if (CoRe) {
            ##If CoRe=TRUE, we need to exchange the Log2FC with the Distance 
            ## and we need to combine the lists
            ## Merge the data frames in list1 and list2 based on the 
            ## "Metabolite" column
            merged_list <- list()
            for (i in seq_len(nrow(name_match_df))) {
                list_dfs <- name_match_df$New[i]

                ## Check if the data frames exist in both lists
                if (list_dfs %in% names(results_list_new) && list_dfs %in% names(Log2FC_table)) {
                    merged_df <- merge(results_list_new[[list_dfs]], 
                        Log2FC_table[[list_dfs]], by = "Metabolite", all = TRUE)
                    merged_list[[list_dfs]] <- merged_df
                }
            }
            STAT_C1vC2 <- merged_list
        } else {
            STAT_C1vC2 <- results_list_new
        }
    }

    ## add input data
    Cond <- colData(se) |>
        as.data.frame() |>
        tibble::rownames_to_column("Code")

    limma_return <- merge(Cond[, c("Code", SettingsInfo[["Conditions"]])], 
        as.data.frame(t(Limma_input)), by.x = "Code", by.y = 0, all.y = TRUE)

    for (DFs in names(STAT_C1vC2)) {
        
        parts <- unlist(strsplit(DFs, "_vs_"))
        C1 <- parts[1]
        C2 <- parts[2]
        limma_return_filt <- limma_return %>%
            dplyr::filter(get(SettingsInfo[["Conditions"]]) == C1 | get(SettingsInfo[["Conditions"]]) == C2) %>%
            tibble::column_to_rownames("Code")
        limma_return_filt <- as.data.frame(t(limma_return_filt[, -c(1)]))

        if (Transform) {
            ## add prefix & suffix to each column since the data have been 
            ## log2 transformed!
            colnames(limma_return_filt) <- paste0("log2(", 
                colnames(limma_return_filt), ")")
        }

        limma_return_merge <- merge(STAT_C1vC2[[DFs]], limma_return_filt, 
            by.x = "Metabolite", by.y = 0, all.x = TRUE)

        STAT_C1vC2[[DFs]] <- limma_return_merge
    }

    ## return
    invisible(STAT_C1vC2)
}


#############################################################################################
### ### ### Shapiro function: Internal Function to perform Shapiro test and plots ### ### ###
#############################################################################################

#' This helper function to perform Shapiro test and plots
#'
#' @param InputData DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected.
#' @param SettingsFile_Sample DF which contains metadata information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames.
#' @param SettingsInfo  \emph{Optional: } Named vector including the information about the conditions column information on numerator or denominator c(Conditions="ColumnName_SettingsFile", Numerator = "ColumnName_SettingsFile", Denominator  = "ColumnName_SettingsFile"). Denominator and Numerator will specify which comparison(s) will be done (one-vs-all, all-vs-one, all-vs-all), e.g. Denominator=NULL and Numerator =NULL selects all the condition and performs multiple comparison all-vs-all. \strong{Default = c(conditions="Conditions", numerator = NULL, denumerator = NULL)}
#' @param StatPval \emph{Optional: } String which contains an abbreviation of the selected test to calculate p.value. For one-vs-one comparisons choose t.test, wilcox.test, "chisq.test" or "cor.test", for one-vs-all or all-vs-all comparison choose aov (=annova), kruskal.test or lmFit (=limma) \strong{Default = "t-test"}
#' @param QQplots \emph {Optional: } TRUE or FALSE for whether QQ plots should be plotted  \strong{default = TRUE}
#'
#' @return List with tewo entries: DF (including the results DF) and Plots (including the Density and QQ plots)
#'
#' @keywords Shapiro test,Normality testing, Density plot, QQplot
#'
#' @importFrom stats shapiro.test
#' @importFrom ggplot2 ggplot geom_histogram geom_density scale_x_continuous theme_minimal labs ggplot_build geom_qq geom_qq_line
#' @importFrom dplyr rename select_if filter
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column column_to_rownames
#'
#' @noRd
#'
Shapiro <- function(se, ##InputData,
    ##SettingsFile_Sample,
    SettingsInfo = c(Conditions = "Conditions", Numerator = NULL, Denominator  = NULL),
    StatPval = "t-test", ## EDIT: the options should be outlined here and match.arg should be used to improve robustness
    QQplots = TRUE) {

    ## ------------ Create log file ----------- ##
    MetaProViz_Init()

    ## ------------- Checks --------------##
    if (grepl("[[:space:]()-./\\\\]", SettingsInfo[["Conditions"]])) { ## EDIT: for such complex regex a comment should be added
        message("In SettingsInfo=c(Conditions= ColumnName): ColumnName contains special charaters, hence this is renamed.")
        ColumnNameCondition_clean <- gsub("[[:space:]()-./\\\\]", "_", SettingsInfo[["Conditions"]])
        SettingsFile_Sample <- SettingsFile_Sample %>%
            dplyr::rename(!!paste(ColumnNameCondition_clean) := SettingsInfo[["Conditions"]])
        SettingsInfo[["Conditions"]] <- ColumnNameCondition_clean
    }

    ## ------------ Denominator/numerator ----------- ##
    ## Denominator and numerator: Define if we compare one_vs_one, one_vs_all or all_vs_all.
    if (!"Denominator" %in% names(SettingsInfo) & !"Numerator" %in% names(SettingsInfo)) { ## EDIT: this is replicated across several functions and should be written as fct
        
        ## all-vs-all: Generate all pairwise combinations
        conditions <- colData(se)[[SettingsInfo[["Conditions"]]]]
        denominator <- unique(conditions)
        numerator <- unique(conditions)
        comparisons <- combn(unique(conditions), 2) %>% 
            as.matrix()
        
        ## settings:
        MultipleComparison <- TRUE
        all_vs_all <- TRUE
    } else if ("Denominator" %in% names(SettingsInfo)  & !"Numerator" %in% names(SettingsInfo)) {
        
        ## all-vs-one: Generate the pairwise combinations
        conditions <- colData(se)[[SettingsInfo[["Conditions"]]]]
        denominator <- SettingsInfo[["Denominator"]]
        numerator <- unique(conditions)
        
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
    
        ## Settings:
        MultipleComparison <- FALSE
        all_vs_all <- FALSE
    }

    ############################################################################
    ## Check data normality and statistical test chosen and generate Output DF
    ## Before Hypothesis testing, we have to decide whether to use a 
    ## parametric or a non parametric test. We can test the data normality 
    ## using the Shapiro test.
    ## First: Load the data and perform the shapiro.test on each metabolite 
    ## across the samples of one condition. this needs to be repeated for 
    ## each condition:
    ## prepare the input (Shapiro test can not handle NAs):
    cols <- colData(se)[[SettingsInfo[["Conditions"]]]] %in% numerator | 
        colData(se)[[SettingsInfo[["Conditions"]]]] %in% denominator
    se_shaptest <- se[, cols]
    assay(se_shaptest)[is.na(assay(se_shaptest))] <- 0
    
    ## we have to remove features with zero variance if there are any.
    ##temp <- sapply(Input_shaptest, function(x, na.rm = TRUE) var(x)) == 0
    rows <- apply(assay(se_shaptest), MARGIN = 1, function(row) var(row, na.rm = TRUE)) != 0 ## EDIT: is this easier?
    se_shaptest <- se_shaptest[rows, ]
    
    ## remove NAs from temp
    ##temp <- temp[complete.cases(temp)]
    
    ## extract column names where temp is TRUE
    ##columns_with_zero_variance <- names(temp[temp])

    #handle a specific case where after filtering and selecting numeric 
    ## variables, there's only one column left in Input_shaptest
    if (nrow(se_shaptest) == 1) {
        se_shaptest <- se
    } else {
        if (!any(rows)) {
            message(
                "The following features have zero variance and are removed prior to performing the shapiro test: ", 
                names(rows[!rows]))
        }
    }

    #Input_shaptest_Cond <- merge(
    #    data.frame(Conditions = SettingsFile_Sample[, SettingsInfo[["Conditions"]], drop = FALSE]), 
     #   Input_shaptest, by = 0, all.y = TRUE)

    UniqueConditions <- colData(se_shaptest)[, SettingsInfo[["Conditions"]]] |>
        unique()
    #UniqueConditions <- SettingsFile_Sample %>%
    #    subset(SettingsFile_Sample[[SettingsInfo[["Conditions"]]]] %in% numerator | SettingsFile_Sample[[SettingsInfo[["Conditions"]]]] %in% denominator, select = c(SettingsInfo[["Conditions"]]))
    #UniqueConditions <- unique(UniqueConditions[[SettingsInfo[["Conditions"]]]])

    ## Generate the results
    shapiro_results <- list()
    for (i in UniqueConditions) {
        ## Subset the data for the current condition
        cols_i <- colData(se_shaptest)[, SettingsInfo[["Conditions"]]] == i
        subset_data <- assay(se_shaptest)[, cols_i] |>
            t() |>
            as.data.frame()
        
        #Input_shaptest_Cond %>%
         #   tibble::column_to_rownames("Row.names") %>%
         #   subset(get(SettingsInfo[["Conditions"]]) == i, select = -c(1))

        ## check the sample size (shapiro.test(x) : sample size must be 
        ## between 3 and 5000):
        if (nrow(subset_data) < 3) {
            warning("shapiro.test(x) : sample size must be between 3 and 5000. You have provided <3 Samples for condition ", 
                i, 
                ". Hence Shaprio test can not be performed for this condition.", 
                sep = "")
        } else if (nrow(subset_data) > 5000) {
            warning("shapiro.test(x) : sample size must be between 3 and 5000. You have provided >5000 Samples for condition ", 
                i, 
                ". Hence Shaprio test will not be performed for this condition.", 
                sep = "")
        } else {
            ## apply Shapiro-Wilk test to each feature in the subset
            shapiro_results[[i]] <- as.data.frame(
                sapply(subset_data, function(x) stats::shapiro.test(x)))
        }
    }

    if (nrow(subset_data) >= 3 & nrow(subset_data) <= 5000) {
        ## make the output DF
        DF_shapiro_results <- as.data.frame(
            matrix(NA, nrow = length(UniqueConditions), ncol = nrow(se_shaptest)))
        rownames(DF_shapiro_results) <- UniqueConditions
        colnames(DF_shapiro_results) <- rownames(se_shaptest)
        for (k in seq_along(UniqueConditions)) {
            for (l in seq_len(ncol(se_shaptest))) {
                DF_shapiro_results[k, l] <- shapiro_results[[UniqueConditions[k]]][[l]]$p.value
            }
        }
        colnames(DF_shapiro_results) <- paste("Shapiro p.val(", 
            colnames(DF_shapiro_results), ")", sep = "")

        ## Second: Give feedback to the user if the chosen test fits the 
        ## data distribution. The data are normal if the p-value of the 
        ## shapiro.test > 0.05.
        Density_plots <- list()
        if (QQplots) {
            QQ_plots <- list()
        }
        for (x in seq_len(nrow(DF_shapiro_results))) {
            
            ## Generate Results Table
            transpose <- as.data.frame(t(DF_shapiro_results[x, ]))
            ## calculate percentage of normally distributed metabolites 
            ## across samples
            Norm <- format(
                round(sum(transpose[[1]] > 0.05) / nrow(transpose), 4) * 100, 
                nsmall = 2) 
            ## calculate percentage of not-normally distributed metabolites across samples
            NotNorm <- format(
                round(sum(transpose[[1]] < 0.05) / nrow(transpose), 4) * 100, 
                nsmall = 2)
            
            if (StatPval == "kruskal.test" | StatPval == "wilcox.test") {
                message("For the condition ", colnames(transpose) ," ", Norm, 
                    " % of the metabolites follow a normal distribution and ", 
                    NotNorm, 
                    " % of the metabolites are not-normally distributed according to the shapiro test. You have chosen ",
                    paste(StatPval), 
                    ", which is for non parametric Hypothesis testing. `shapiro.test` ignores missing values in the calculation.")
            } else {
                message("For the condition ", colnames(transpose) ," ", Norm, 
                    " % of the metabolites follow a normal distribution and ", 
                    NotNorm, 
                    " % of the metabolites are not-normally distributed according to the shapiro test. You have chosen ",
                    paste(StatPval), 
                    ", which is for parametric Hypothesis testing. `shapiro.test` ignores missing values in the calculation.")
            }

            ## assign the calculated values to the corresponding rows in result_df
            DF_shapiro_results$`Metabolites with normal distribution [%]`[x] <- Norm
            DF_shapiro_results$`Metabolites with not-normal distribution [%]`[x] <- NotNorm

            ## reorder the DF:
            all_columns <- colnames(DF_shapiro_results)
            include_columns <- c("Metabolites with normal distribution [%]", "Metabolites with not-normal distribution [%]")
            exclude_columns <- setdiff(all_columns, include_columns)
            DF_shapiro_results <- DF_shapiro_results[, c(include_columns, exclude_columns)]

            ## make Group wise data distribution plot and QQ plots
            cols <- colData(se_shaptest)[, SettingsInfo[["Conditions"]]] == colnames(transpose)
            subset_data <- assay(se_shaptest)[, cols]
            #subset_data <-  %>%
            #    tibble::column_to_rownames("Row.names") %>%
            #    subset(get(SettingsInfo[["Conditions"]]) ==  colnames(transpose), 
            #        select = -c(1))
            all_data <- unlist(subset_data)

            plot <- ggplot2::ggplot(data.frame(x = all_data), aes(x = x)) +
                ggplot2::geom_histogram(ggplot2::aes(y = after_stat(density)), 
                    binwidth = .5, colour = "black", fill = "white")  +
                ggplot2::geom_density(alpha = 0.2, fill = "grey45")

            density_values <- ggplot2::ggplot_build(plot)$data[[2]]

            plot <- ggplot2::ggplot(data.frame(x = all_data), aes(x = x)) + ## EDIT: why replot, what about plot + scale_x_continous(...)
                ggplot2::geom_histogram(ggplot2::aes(y = after_stat(density)), 
                    binwidth = .5, colour = "black", fill = "white") +
                ggplot2::geom_density(alpha = .2, fill = "grey45") +
                ggplot2::scale_x_continuous(limits = c(0, density_values$x[max(which(density_values$scaled >= 0.1))]))

            density_values2 <- ggplot2::ggplot_build(plot)$data[[2]]

            suppressWarnings(
                sampleDist <- ggplot2::ggplot(data.frame(x = all_data), aes(x = x)) + ## EDIT: why replot, what about plot + theme_minimal(...) + labs(...)
                    ggplot2::geom_histogram(aes(y = after_stat(density)), binwidth = .5, colour = "black", fill = "white") +
                    ggplot2::geom_density(alpha = .2, fill = "grey45") +
                    ggplot2::scale_x_continuous(limits = c(0, density_values$x[max(which(density_values$scaled >= 0.1))])) +
                    ggplot2::theme_minimal()+
                    ggplot2::labs(title=paste("Data distribution ",  colnames(transpose)), subtitle = paste(NotNorm, " of metabolites not normally distributed based on Shapiro test"),x="Abundance", y = "Density")
            )

            Density_plots[[paste(colnames(transpose))]] <- sampleDist

            # QQ plots
            if (QQplots) {
                ## make folders
                conds <- unique(c(numerator, denominator))

                ## QQ plots for each groups for each metabolite for normality visual check
                qq_plot_list <- list()
                for (row_name in rownames(subset_data)) {
                    qq_plot <- ggplot2::ggplot(
                            data.frame(x = subset_data[row_name,]), 
                            aes(sample = x)) +
                        ggplot2::geom_qq() +
                        ggplot2::geom_qq_line(color = "red") +
                        ggplot2::labs(title = paste("QQPlot for", row_name),
                            x = "Theoretical", y = "Sample") + 
                        theme_minimal()
                    plot.new()
                    plot(qq_plot)
                    qq_plot_list[[row_name]] <- recordPlot() ## EDIT: not sure if it works but could it be just qq_plot_list[[col_name]] <- qq_plot (delete plot.new, recordPlot)

                    row_name2 <- (gsub("/", "_", row_name))#remove "/" cause this can not be safed in a PDF name
                    row_name2 <- gsub("-", "", row_name2)
                    row_name2 <- gsub("/", "", row_name2)
                    row_name2 <- gsub(" ", "", row_name2)
                    row_name2 <- gsub("\\*", "", row_name2)
                    row_name2 <- gsub("\\+", "", row_name2)
                    row_name2 <- gsub(",", "", row_name2)
                    row_name2 <- gsub("\\(", "", row_name2)
                    row_name2 <- gsub("\\)", "", row_name2) ## EDIT: row_name2 not used in the function?

                    dev.off() ## EDIT: is this needed?
                }
                QQ_plots[[paste(colnames(transpose))]] <- qq_plot_list
            }
        }

        ######################################
        ##-------- Return
        
        ## here we make a list
        l_shapiro <- list(
            "data" = list(
                "Shapiro_result" = tibble::rownames_to_column(
                    DF_shapiro_results, "Code")),
            "plot" = list(
                "Distributions" = Density_plots))
        
        if (QQplots) {
            l_shapiro[["plot"]][["QQ_plots"]] <- QQ_plots
        }

        ## return
        suppressWarnings(invisible(l_shapiro))
    }
}

################################################################################
### Bartlett function: Internal Function to perform Bartlett test and plots  ###
################################################################################

#' This helper function perform the Bartlett test to check the homogeneity of variances across groups
#'
#' @param InputData DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected.
#' @param SettingsFile_Sample DF which contains metadata information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames.
#' @param SettingsInfo  \emph{Optional: } Named vector including the information about the conditions column information on numerator or denominator c(Conditions="ColumnName_SettingsFile", Numerator = "ColumnName_SettingsFile", Denominator  = "ColumnName_SettingsFile"). Denominator and Numerator will specify which comparison(s) will be done (one-vs-all, all-vs-one, all-vs-all), e.g. Denominator=NULL and Numerator =NULL selects all the condition and performs multiple comparison all-vs-all. \strong{Default = c(conditions="Conditions", numerator = NULL, denumerator = NULL)}
#'
#' @return List with two entries: DF (including the results DF) and Plots (including the  histogramm plot)
#'
#' @keywords Bartlett test,Normality testing, Density plot, QQplot
#'
#' @importFrom stats bartlett.test
#' @importFrom ggplot2 ggplot geom_histogram geom_density ggtitle xlab geom_vline
#' @importFrom dplyr mutate
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#'
#' @noRd
#'
Bartlett <-function(se, ##InputData,
    ##SettingsFile_Sample,
    SettingsInfo) {

    ## ------------ Create log file ----------- ##
    MetaProViz_Init()

    ############################################################################

    conditions <- colData(se)[[SettingsInfo[["Conditions"]]]]

    # Use Bartletts test
    bartlett_res <-  apply(assay(se), 1, 
        function(x) stats::bartlett.test(x ~ conditions))

    ## make the output DF
    DF_bartlett_results <- data.frame("Bartlett_pvalue" = rep(NA, nrow(se)))
    #matrix(NA, nrow = nrow(InputData)), ncol = 1)
    rownames(DF_bartlett_results) <- rownames(se)
    ##colnames(DF_bartlett_results) <- "Bartlett p.val"

    for (i in seq_along(bartlett_res)) {
        DF_bartlett_results[i, "Bartlett_pvalue"] <- bartlett_res[[i]]$p.value
    }
    DF_bartlett_results <- DF_bartlett_results %>% 
        dplyr::mutate(
            Variance_homogeneity = case_when(
                Bartlett_pvalue <  0.05 ~ FALSE,
                Bartlett_pvalue >= 0.05 ~ TRUE))
    
    ## if p < 0.05 then unequal variances
    message("For ", 
        round(
            sum(DF_bartlett_results$Variance_homogeneity) / nrow(DF_bartlett_results), 
            digits = 4) * 100, 
        "% of metabolites the group variances are equal.")

    DF_bartlett_results <- DF_bartlett_results %>% 
        tibble::rownames_to_column("Metabolite") %>% 
        relocate("Metabolite")

    #### Plots:
    ## Make density plots
    Bartlettplot <- ggplot2::ggplot(
            data.frame(x = DF_bartlett_results), 
            aes(x = DF_bartlett_results$Bartlett_pvalue)) +
        ggplot2::geom_histogram(aes(y = ..density..), 
            colour = "black", fill = "white")  +
        ggplot2::geom_density(alpha = 0.2, fill = "grey45")+
        ggplot2::ggtitle("Bartlett's test p.value distribution") +
        ggplot2::xlab("p.value")+
        ggplot2::geom_vline(aes(xintercept = 0.05, color = "darkred"))

    Bartlett_output_list <- list(
        "data" = list("Bartlett_result" = DF_bartlett_results), 
        "plot" = list("Histogram" = Bartlettplot))

    ## return
    suppressWarnings(invisible(Bartlett_output_list))

}


################################################################
### ### ### Variance stabilizing transformation function ### ###
################################################################

#' @title Variance stabilizing transformation (VST) on assay of SummarizedExperiment
#' 
#' @description This function performs a variance stabilizing transformation 
#' (VST) on the assay of a SummarizedExperiment.
#'
#' @param se SummarizedExperiment with unique sample identifiers as colnames 
#'     and metabolite numerical values in columns with metabolite identifiers 
#'     as rownames. Use NA for metabolites that were not detected.
#'
#' @return list with two entries: "data" (including the vst-adjusted assay) and 
#' "plot" (including the scedasticity_plot)
#'
#' @keywords Heteroscedasticity, variance stabilizing transformation
#'
#' @importFrom reshape2 melt
#' @importFrom stats lm
#' @importFrom ggplot2 ggplot geom_point theme_bw scale_x_continuous scale_y_continuous xlab ylab geom_abline ggtitle geom_smooth aes
#' @importFrom patchwork wrap_plots
#' @importFrom dplyr summarise group_by
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#'
#' @noRd
#'
vst <- function(se) {

    ## ------------ Create log file ----------- ##
    MetaProViz_Init()

    ## model the mean and variance relationship on the data
    suppressMessages(
        melted <- assay(se) |>
            t() |>
            reshape2::melt())
    het_data <- melted %>%
        ## make a dataframe to save the values
        dplyr::group_by(Var2) %>% ## EDIT: should this be according to metabolites? ## was: variable
        dplyr::summarise(mean = mean(value), sd = sd(value))
    ## add a common group for the lm function to account for the whole data together
    het_data$lm <- 1

    invisible(het_plot <- ggplot2::ggplot(het_data, ggplot2::aes(x = mean, y = sd)) +
        ggplot2::geom_point() +
        ggplot2::theme_bw() +
        ggplot2::scale_x_continuous(trans='log2') +
        ggplot2::scale_y_continuous(trans='log2') +
        ggplot2::xlab("log(mean)") +
        ggplot2::ylab("log(sd)") +
        ggplot2::geom_abline(intercept = 0, slope = 1)  +
        ggplot2::ggtitle(" Data heteroscedasticity")  +
        ggplot2::geom_smooth(ggplot2::aes(group = lm), method = "lm", 
            formula = y~x, color = "red"))

    # select data
    prevst_data <- het_data
    prevst_data$mean <- log(prevst_data$mean)
    prevst_data$sd <- log(prevst_data$sd)
    
    ## calculate the slope of the log data ## EDIT: why are you not using some prebuild function from e.g. vsn2?
    data_fit <- stats::lm(sd ~ mean, prevst_data)
    coef(data_fit)

    ## make the vst transformation
    data_vst <- as.data.frame(t(assay(se)) ^ (1 - coef(data_fit)["mean"][1]))

    ## heteroscedasticity visual check again
    suppressMessages(melted_vst <- reshape::melt(data_vst))
    het_vst_data <- melted_vst %>%
        ## make a dataframe to save the values
        dplyr::group_by(variable) %>%
        dplyr::summarise(mean = mean(value), sd = sd(value))
    ## add a common group for the lm function to account for the whole data together
    het_vst_data$lm <- 1 

    ## plot variable stadard deviation as a function of the mean
    invisible(hom_plot <- ggplot2::ggplot(het_vst_data,  
            ggplot2::aes(x = mean, y = sd)) +
        ggplot2::geom_point() +
        ggplot2::theme_bw() +
        ggplot2::scale_x_continuous(trans = "log2") +
        ggplot2::scale_y_continuous(trans = "log2") +
        ggplot2::xlab("log(mean)") +
        ggplot2::ylab("log(sd)") +
        ggplot2::geom_abline(intercept = 0)  +
        ggplot2::ggtitle("Vst transformed data")  +
        ggplot2::geom_smooth(ggplot2::aes(group = lm), method='lm', 
            formula = y ~ x, color = "red"))

    scedasticity_plot <- patchwork::wrap_plots(het_plot, hom_plot)

    ## assemble the object to return
    assay(se) <- t(data_vst)
    l <- list(
        "data" = list(
            "se" = se,
            "assay" = data_vst), 
        "plot" = list(
            "scedasticity_plot" = scedasticity_plot))

    ## return
    invisible(l)
}

