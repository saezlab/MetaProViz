## ---------------------------
##
## Script name: MetaAnalysis
##
## Purpose of script:
##
## Author: Christina Schmidt
##
## Date Created: 2024-01-10
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
#'

###############################################
### ### ### MetaAnalysis ### ### ###
###############################################

#' This function performs a PCA analysis on the input data and combines it with the sample metadata to perform an ANOVA test to identify significant differences between the groups.
#'
#' @param InputData DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected. includes experimental design and outlier column.
#' @param SettingsFile_Sample \emph{Optional: } DF which contains information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames. Column "Conditions" with information about the sample conditions (e.g. "N" and "T" or "Normal" and "Tumor"), can be used for feature filtering and colour coding in the PCA. Column "AnalyticalReplicate" including numerical values, defines technical repetitions of measurements, which will be summarised. Column "BiologicalReplicates" including numerical values. Please use the following names: "Conditions", "Biological_Replicates", "Analytical_Replicates".\strong{Default = NULL}
#' @param Scaling \emph{Optional: } TRUE or FALSE for whether a data scaling is used \strong{Default = TRUE}
#' @param Percentage \emph{Optional: } Percentage of top and bottom features to be displayed in the results summary. \strong{Default = 0.1}
#' @param StatCutoff \emph{Optional: } Cutoff for the adjusted p-value of the ANOVA test for the results summary and on the heatmap. \strong{Default = 0.05}
#' @param VarianceCutoff \emph{Optional: } Cutoff for the PCs variance that should be displayed on the heatmap. \strong{Default = 1}
#' @param SaveAs_Plot \emph{Optional: } Select the file type of output plots. Options are svg, png, pdf. \strong{Default = svg}
#' @param SaveAs_Table \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt". \strong{Default = "csv"}
#' @param PrintPlot \emph{Optional: } TRUE or FALSE, if TRUE Volcano plot is saved as an overview of the results. \strong{Default = TRUE}
#' @param FolderPath \emph{Optional:} Path to the folder the results should be saved at. \strong{default: NULL}
#'
#' @return List of DFs: prcomp results, loadings, Top-Bottom features, annova results, results summary
#'
#' @examples
#' Tissue_Norm <- ToyData("Tissue_Norm")
#' Res <- MetaAnalysis(InputData=Tissue_Norm[,-c(1:13)],
#'                                 SettingsFile_Sample= Tissue_Norm[,c(2,4:5,12:13)])
#'
#' @keywords PCA, annova, metadata
#'
#' @importFrom dplyr filter bind_rows rename mutate ungroup group_by summarise select arrange rowwise mutate_all distinct
#' @importFrom magrittr %>%
#' @importFrom stats as.formula aov TukeyHSD
#' @importFrom broom tidy
#' @importFrom tidyr separate_rows
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom logger log_info
#'
#' @export
#'
MetaAnalysis <- function(InputData,
    SettingsFile_Sample,
    Scaling = TRUE,
    Percentage = 0.1,
    StatCutoff = 0.05,
    VarianceCutoff = 1,
    SaveAs_Table = "csv",
    SaveAs_Plot = "svg",
    PrintPlot = TRUE,
    FolderPath = NULL) {
    #SettingInfo= c(MainSeparator = "TISSUE_TYPE), # enable this parameter in the function --> main separator: Often a combination of demographics is is of paricular interest, e.g. comparing "Tumour versus Normal" for early stage patients and for late stage patients independently. If this is the case, we can use the parameter `SettingsInfo` and provide the column name of our main separator.

    ## ------------ Create log file ----------- ##
    MetaProViz_Init()

    ############################################################################
    ## ------------ Check Input files ----------- ##
    ## HelperFunction `CheckInput`
    CheckInput(InputData = InputData, SettingsFile_Sample = SettingsFile_Sample,
        SettingsFile_Metab = NULL, SettingsInfo = NULL, 
        SaveAs_Plot = SaveAs_Plot, SaveAs_Table = SaveAs_Table,
        CoRe = FALSE, PrintPlot = PrintPlot)

    ## Specific checks: Check the column names of the demographics --> 
    ## need to be R usable (no empty spaces, -, etc.)
    if (!is.null(SettingsFile_Sample)) { 
        if(any(grepl('[^[:alnum:]]', colnames(SettingsFile_Sample)))) {
            
            ## Remove special characters in colnames
            colnames(SettingsFile_Sample) <- make.names(colnames(SettingsFile_Sample))
            
            ## Message
            message <- paste("The column names of the 'SettingsFile_Sample' contain special character that where removed.")
            logger::log_info(message)
            message(message)
        }
    }

    ## ------------ Create Results output folder ----------- ##
    if(!is.null(SaveAs_Plot) | !is.null(SaveAs_Table)){
        Folder <- SavePath(FolderName = "MetaAnalysis", 
            FolderPath = FolderPath)
    }

    ############################################################################
    ## ---------- Run prcomp ------------##
    ##--- 1. prcomp
    ## Get PCs
    PCA.res <- prcomp(InputData, center = TRUE, scale = Scaling)
    PCA.res_Info <- as.data.frame(PCA.res$x)

    ## Extract loadings for each PC
    PCA.res_Loadings <- as.data.frame(PCA.res$rotation) %>%
        tibble::rownames_to_column("FeatureID")

    ##--- 2. Merge with demographics
    PCA.res_Info  <- merge(
        x = tibble::rownames_to_column(SettingsFile_Sample, "UniqueID"),
        y = tibble::rownames_to_column(PCA.res_Info, "UniqueID"),
        by = "UniqueID", all.y = TRUE) %>%
        tibble::column_to_rownames("UniqueID")

    ##--- 3. convert columns that are not numeric to factor:
    ## Demographics are often non-numerical, categorical explanatory 
    ## variables, which is often stored as characters, sometimes integers
    PCA.res_Info[sapply(PCA.res_Info, is.character)] <- lapply(
        PCA.res_Info[sapply(PCA.res_Info, is.character)], as.factor)
    PCA.res_Info[sapply(PCA.res_Info, is.integer)] <- lapply(
        PCA.res_Info[sapply(PCA.res_Info, is.integer)], as.factor)

    ############################################################################
    ## ---------- STATS ------------##
    ## 1. Anova p.val
    ## Iterate through each combination of meta and PC columns
    MetaData <- names(SettingsFile_Sample)

    Stat_results <- list()

    for (meta_col in MetaData) {
        for (pc_col in colnames(PCA.res_Info)[grepl("^PC", colnames(PCA.res_Info))]) {
            ## Create a formula for ANOVA --> When constructing the ANOVA 
            ## formula, it's important to ensure that the response variable 
            ## (dependent variable) is numeric.
            Formula <- stats::as.formula(paste(pc_col, "~", meta_col, sep = ""))
            #pairwiseFormula <- as.formula(paste("pairwise ~" , meta_col, sep=""))

            ## Perform ANOVA
            anova_result <- stats::aov(Formula, data = PCA.res_Info)
            ## broom::tidy --> convert statistics into table
            anova_result_tidy <- broom::tidy(anova_result)
            ## exclude Residuals row
            anova_row <- dplyr::filter(anova_result_tidy, term != "Residuals")  

            ## Perform Tukey test
            tukey_result <- stats::TukeyHSD(anova_result)
            ## #broom::tidy --> convert statistics into table
            tukey_result_tidy <- broom::tidy(tukey_result) ## EDIT: what is this fct doign, checking ?broom::tidy it seems it is just a call to generics::tidy, could dependencies be reduced?

            #lm_result_p.val <- lsmeans::lsmeans((lm(Formula, data = PCA.res_Info)),  pairwiseFormula , adjust=NULL)# adjust=NULL leads to the p-value!
            #lm_result_p.adj <- lsmeans::lsmeans((lm(Formula, data = PCA.res_Info)), pairwiseFormula, adjust="FDR")
            #lm_result <- as.data.frame(lm_result_p.val$contrasts)
            #lm_result$p.adj <- as.data.frame(lm_result_p.adj$contrasts)[,"p.value"]

            ## combine results
            combined_result <- data.frame(
                tukeyHSD_Contrast = tukey_result_tidy$contrast,
                PC = pc_col,
                term = anova_row$term,
                anova_sumsq = anova_row$sumsq,
                anova_meansq = anova_row$meansq,
                anova_statistic = anova_row$statistic,
                anova_p.value = anova_row$p.value,
                tukeyHSD_p.adjusted = tukey_result_tidy$adj.p.value
            )

            ## store ANOVA result
            Stat_results[[paste(pc_col, meta_col, sep = "_")]] <- combined_result
        }
    }

    ## add into one DF
    Stat_results <- dplyr::bind_rows(Stat_results)

    ## Add explained variance into the table:
    prop_var_ex <- as.data.frame(
        ## To compute the proportion of variance explained by each component in 
        ## percent, we divide the variance by sum of total variance and 
        ## multiply by 100(variance=standard deviation ^2)
            ((PCA.res[["sdev"]]) ^ 2 / sum((PCA.res[["sdev"]]) ^ 2)) * 100) %>%
        tibble::rownames_to_column("PC") %>%
        dplyr::mutate(PC = paste("PC", PC, sep = "")) %>%
        dplyr::rename("Explained_Variance" = 2)

    Stat_results <- merge(Stat_results, prop_var_ex,  by = "PC", all.x = TRUE)

    ############################################################################
    ## ---------- Top/Bottom ------------##
    ## add top/bottom related features to this
    ## create a data frame for top and bottom features for each PC
    TopBottom_Features <- lapply(seq_len(ncol(PCA.res_Loadings))[-1], function(i) {
        ## Make input
        pc_loadings <- PCA.res_Loadings[, c("FeatureID", 
            names(PCA.res_Loadings)[i])]
        ## get top and bottom features
        n_features <- nrow(pc_loadings)
        n_selected <- round(Percentage * n_features)

        top_features <- as.data.frame(head(
            arrange(pc_loadings, desc(!!sym(names(pc_loadings)[2]))), 
                n_selected)$FeatureID) %>%
            dplyr::rename(
                !!paste0("Features_", "(Top", Percentage, "%)") := 1)
        bottom_features <- as.data.frame(head(
            arrange(pc_loadings, !!sym(names(pc_loadings)[2])), 
            n_selected)$FeatureID) %>% ## EDIT: the calculation is done twice, should be simplified
            dplyr::rename(
                !!paste0("Features_", "(Bottom", Percentage, "%)") := 1)

        ## return
        cbind(data.frame(PC = paste("PC", i, sep = "")), 
            top_features, bottom_features)
        }) %>%
        dplyr::bind_rows(.id = "PC") %>%
        dplyr::mutate(PC = paste0("PC", PC))


    ## ---------- Final DF 1------------##
    ## Add to results DF
    Stat_results <- merge(Stat_results,
                         TopBottom_Features%>% ## EDIT: not sure if its really straitghtforward to understand what is done here, could it be simplified / rewritten?
                           dplyr::group_by(PC) %>%
                           dplyr::summarise(across(everything(), ~ paste(unique(gsub(", ", "_", .)), collapse = ", "))) %>%
                           dplyr::ungroup(),
                         by = "PC",
                         all.x = TRUE)

    ## ---------- DF 2: Metabolites as row names ------------##
    Res_Top <- Stat_results %>%
        dplyr::filter(tukeyHSD_p.adjusted < StatCutoff) %>%
        ## separate 'Features (Top 0.1%)'
        tidyr::separate_rows(paste0("Features_", "(Top", Percentage, "%)"), sep = ", ") %>% 
        dplyr::rename("FeatureID":= paste0("Features_", "(Top", Percentage, "%)")) %>%
        dplyr::select(- paste0("Features_", "(Bottom", Percentage, "%)"))

    Res_Bottom <- Stat_results%>%
        dplyr::filter(tukeyHSD_p.adjusted < StatCutoff) %>%
        ## separate 'Features (Bottom 0.1%)'
        tidyr::separate_rows(paste0("Features_", "(Bottom", Percentage, "%)"), sep = ", ") %>%    
        dplyr::rename("FeatureID":= paste0("Features_", "(Bottom", Percentage, "%)")) %>%
        dplyr::select(- paste0("Features_", "(Top", Percentage, "%)"))

    Res <- rbind(Res_Top, Res_Bottom) %>%
        ## group by FeatureID and term
        dplyr::group_by(FeatureID, term) %>%            
        dplyr::summarise(
            ## concatenate unique PC entries with commas
            PC = paste(unique(PC), collapse = ", "), 
            ## sum Explained_Variance
            `Sum(Explained_Variance)` = sum(Explained_Variance, na.rm = TRUE)) %>% ## EDIT: I think colnames should be after applying make.names
        ## remove previous grouping
        dplyr::ungroup() %>%  
        ## group by FeatureID for MainDriver calculation
        dplyr::group_by(FeatureID) %>%  
        ## mark TRUE for the highest value
        dplyr::mutate(MainDriver = (`Sum(Explained_Variance)` == max(`Sum(Explained_Variance)`))) %>% 
        ## remove grouping
        dplyr::ungroup()  

    Res_summary <- Res%>%
        dplyr::group_by(FeatureID) %>%
        dplyr::summarise(
            ## concatenate all terms separated by commas
            term = paste(term, collapse = ", "),
            ## concatenate all Sum(Explained_Variance) values
            `Sum(Explained_Variance)` = paste(`Sum(Explained_Variance)`, collapse = ", "),  
            ## extract the term where MainDriver is TRUE
            MainDriver = paste(MainDriver, collapse = ", ")) %>%  
        dplyr::ungroup()%>%
        dplyr::rowwise() %>%
        dplyr::mutate(  
            ## extract the term where MainDriver is TRUE
            MainDriver_Term = ifelse("TRUE" %in% strsplit(MainDriver, ", ")[[1]],
                strsplit(term, ", ")[[1]][which(strsplit(MainDriver, ", ")[[1]] == "TRUE")[1]],
                NA),
            ## extract the Sum(Explained_Variance) where MainDriver is TRUE
            `MainDriver_Sum(VarianceExplained)` = ifelse("TRUE" %in% strsplit(MainDriver, ", ")[[1]],
                as.numeric(strsplit(`Sum(Explained_Variance)`, ", ")[[1]][which(strsplit(MainDriver, ", ")[[1]] == "TRUE")[1]]),
                NA)) %>%
        dplyr::ungroup()%>%
        dplyr::arrange(desc(`MainDriver_Sum(VarianceExplained)`))

    #############################################################################
    ## ---------- Plot ------------##
    ## plot DF
    Data_Heat <- Stat_results %>%
        ## filter for significant results
        dplyr::filter(tukeyHSD_p.adjusted < StatCutoff) %>%
        ## exclude Residuals row
        dplyr::filter(Explained_Variance > VarianceCutoff) %>%
        ## only keep unique term~PC combinations AND STATS
        dplyr::distinct(term, PC, .keep_all = TRUE) %>% 
        dplyr::select(term, PC, Explained_Variance)

    Data_Heat <- reshape2::dcast(Data_Heat, term ~ PC, value.var = "Explained_Variance") %>%
        tibble::column_to_rownames("term") %>%
        dplyr::mutate_all(~replace(., is.na(.), 0))

    if (nrow(Data_Heat) > 2) {

        ## Plot
        invisible(
            VizHeatmap(InputData = Data_Heat,
                PlotName = paste0("ExplainedVariance-bigger-", 
                    VarianceCutoff , "Percent_AND_p.adj-smaller", 
                    StatCutoff),
                Scale = "none",
                SaveAs_Plot = SaveAs_Plot,
                PrintPlot = PrintPlot,
                FolderPath = Folder))
    } else {
        message <- paste0("StatCutoff of ", StatCutoff, 
            " and VarianceCutoff of ", VarianceCutoff, 
            " do only return <= 2 cases, hence no heatmap is plotted.")
        logger::log_info("warning: ", message)
        warning(message)
    }

    #############################################################################
    ## ---------- Return ------------##
    ## make list of Output DFs: 1. prcomp results, 2. Loadings result, 
    ## 3. TopBottom Features, 4. AOV
    ResList <- list(res_prcomp = PCA.res_Info, res_loadings = PCA.res_Loadings, 
        res_TopBottomFeatures = TopBottom_Features, res_aov = Stat_results, 
        res_summary = Res_summary)

    ## Save the results
    SaveRes(
        InputList_DF = ResList,
        InputList_Plot = NULL,
        SaveAs_Table = SaveAs_Table,
        SaveAs_Plot = FALSE,
        FolderPath = Folder,
        FileName = "MetaAnalysis",
        CoRe = FALSE,
        PrintPlot = FALSE)

    ## return
    invisible(ResList)
}


###############################################
### ### ### MetaPriorKnowlegde ### ### ###
###############################################

#' Meta prior-knowledge
#'
#' @param InputData DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected. includes experimental design and outlier column.
#' @param SettingsFile_Sample \emph{Optional: } DF which contains information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames. Column "Conditions" with information about the sample conditions (e.g. "N" and "T" or "Normal" and "Tumor"), can be used for feature filtering and colour coding in the PCA. Column "AnalyticalReplicate" including numerical values, defines technical repetitions of measurements, which will be summarised. Column "BiologicalReplicates" including numerical values. Please use the following names: "Conditions", "Biological_Replicates", "Analytical_Replicates".\strong{Default = NULL}
#' @param SettingsInfo \emph{Optional: } NULL or vector with column names that should be used, i.e. c("Age", "gender", "Tumour-stage"). \strong{default: NULL}
#' @param SaveAs_Table \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt". \strong{Default = "csv"}
#' @param FolderPath \emph{Optional:} Path to the folder the results should be saved at. \strong{default: NULL}
#'
#' @return DF with prior knowledge based on patient metadata
#'
#' @examples
#' Tissue_Norm <- ToyData("Tissue_Norm")
#' Res <- MetaPK(InputData=Tissue_Norm[,-c(1:13)],
#'                           SettingsFile_Sample= Tissue_Norm[,c(2,4:5,12:13)])
#'
#' @keywords prior knowledge, metadata
#'
#' @importFrom magrittr %>%
#' @importFrom tidyr pivot_longer unite
#' @importFrom tibble rownames_to_column
#'
#' @export
#'
MetaPK <- function(InputData,
    SettingsFile_Sample,
    SettingsInfo = NULL,
    SaveAs_Table = "csv",
    FolderPath = NULL) {

    ## ------------ Create log file ----------- ##
    MetaProViz_Init()

    ### Enrichment analysis-based
    ## *we can make pathway file from metadata and use this to run enrichment 
    ## analysis.
    ## *At the moment we can just make the pathways and we need to include a 
    ## standard fishers exact test.
    ##*Merge with function above?
    # *Advantage/disadvantage: Not specific - with annova PC we get granularity 
    ##e.g. smoker-ExSmoker-PC5, whilst here we only get smoking in generall 
    ## as parameter

    ############################################################################
    ## ------------ Check Input files ----------- ##
    ## HelperFunction `CheckInput`

    ## ------------ Create Results output folder ----------- ##
    if (!is.null(SaveAs_Table)) {
        Folder <- SavePath(FolderName = "MetaAnalysis", FolderPath = FolderPath)
    }

    ############################################################################
    ## create Prior Knowledge file format to perform enrichment analysis ##
    ## use the Sample metadata for this:
    if (is.null(SettingsInfo)) {
        MetaData <- names(SettingsFile_Sample)
        SettingsFile_Sample_subset <- SettingsFile_Sample %>%
            tibble::rownames_to_column("SampleID")
    } else {
        MetaData <- SettingsInfo
        SettingsFile_Sample_subset <- SettingsFile_Sample[, MetaData, 
                drop = FALSE] %>%
            tibble::rownames_to_column("SampleID")
    }

    ## convert into a pathway DF
    Metadata_df <- SettingsFile_Sample_subset %>%
        tidyr::pivot_longer(cols = -SampleID, 
            names_to = "ColumnName", values_to = "ColumnEntry") %>%
        tidyr::unite("term", c("ColumnName", "ColumnEntry"), sep = "_")
    Metadata_df$mor <- 1

    ## run ULM using decoupleR (input=PCs from prcomp and pkn=Metadata_df)
    #ULM_res <- decoupleR::run_ulm(mat=as.matrix(PCA.res$x),
    #                              net=Metadata_df,
    #                              .source='term',
    #                              .target='SampleID',
    #                              .mor='mor',
    #                              minsize = 5)

    ## add explained variance into the table:
    #prop_var_ex <- as.data.frame(((PCA.res[["sdev"]])^2/sum((PCA.res[["sdev"]])^2))*100)%>%#To compute the proportion of variance explained by each component in percent, we divide the variance by sum of total variance and multiply by 100(variance=standard deviation ^2)
    #    rownames_to_column("PC")%>%
    #    mutate(PC = paste("PC", PC, sep=""))%>%
    #   dplyr::rename("Explained_Variance"=2)

    #ULM_res <- merge(ULM_res, prop_var_ex,  by.x="condition",by.y="PC",all.x=TRUE)

    ############################################################################
    ## ---------- save ------------##
    ## add to results DF
    Res <- list(MetaData_PriorKnowledge = Metadata_df)

    ## return
    Res
}

