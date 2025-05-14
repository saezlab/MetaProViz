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
#' @param data DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected. includes experimental design and outlier column.
#' @param metadata_sample \emph{Optional: } DF which contains information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames. Column "Conditions" with information about the sample conditions (e.g. "N" and "T" or "Normal" and "Tumor"), can be used for feature filtering and colour coding in the PCA. Column "AnalyticalReplicate" including numerical values, defines technical repetitions of measurements, which will be summarised. Column "BiologicalReplicates" including numerical values. Please use the following names: "Conditions", "Biological_Replicates", "Analytical_Replicates".\strong{Default = NULL}
#' @param Scaling \emph{Optional: } TRUE or FALSE for whether a data scaling is used \strong{Default = TRUE}
#' @param Percentage \emph{Optional: } Percentage of top and bottom features to be displayed in the results summary. \strong{Default = 0.1}
#' @param StatCutoff \emph{Optional: } Cutoff for the adjusted p-value of the ANOVA test for the results summary and on the heatmap. \strong{Default = 0.05}
#' @param VarianceCutoff \emph{Optional: } Cutoff for the PCs variance that should be displayed on the heatmap. \strong{Default = 1}
#' @param save_plot \emph{Optional: } Select the file type of output plots. Options are svg, png, pdf. \strong{Default = svg}
#' @param save_table \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt". \strong{Default = "csv"}
#' @param print_plot \emph{Optional: } TRUE or FALSE, if TRUE Volcano plot is saved as an overview of the results. \strong{Default = TRUE}
#' @param path \emph{Optional:} Path to the folder the results should be saved at. \strong{default: NULL}
#'
#' @return List of DFs: prcomp results, loadings, Top-Bottom features, annova results, results summary
#'
#' @examples
#' Tissue_Norm <- MetaProViz::ToyData("Tissue_Norm")
#' Res <- MetaProViz::MetaAnalysis(data=Tissue_Norm[,-c(1:13)],
#'                                 metadata_sample= Tissue_Norm[,c(2,4:5,12:13)])
#'
#' @keywords PCA, annova, metadata
#'
#' @importFrom dplyr filter bind_rows rename mutate ungroup group_by summarise select arrange rowwise mutate_all distinct
#' @importFrom magrittr %>%
#' @importFrom stats as.formula aov TukeyHSD
#' @importFrom broom tidy
#' @importFrom tidyr separate_rows pivot_wider
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom logger log_info
#'
#' @export
#'
MetaAnalysis <- function(data,
                         metadata_sample,
                         Scaling = TRUE,
                         Percentage = 0.1,
                         StatCutoff= 0.05,
                         VarianceCutoff=1,
                         save_table = "csv",
                         save_plot = "svg",
                         print_plot= TRUE,
                         path = NULL
                         #SettingInfo= c(MainSeparator = "TISSUE_TYPE), # enable this parameter in the function --> main separator: Often a combination of demographics is is of paricular interest, e.g. comparing "Tumour versus Normal" for early stage patients and for late stage patients independently. If this is the case, we can use the parameter `metadata_info` and provide the column name of our main separator.

){
  ## ------------ Create log file ----------- ##
  MetaProViz_Init()

  ################################################################################################################################################################################################
  ## ------------ Check Input files ----------- ##
  # HelperFunction `CheckInput`
  CheckInput(data=data,
                          metadata_sample=metadata_sample,
                          metadata_feature=NULL,
                          metadata_info=NULL,
                          save_plot=save_plot,
                          save_table=save_table,
                          core=FALSE,
                          print_plot= print_plot)

  # Specific checks: Check the column names of the demographics --> need to be R usable (no empty spaces, -, etc.)
  if(is.null(metadata_sample)==FALSE){
    if(any(grepl('[^[:alnum:]]', colnames(metadata_sample)))==TRUE){
      #Remove special characters in colnames
      colnames(metadata_sample) <- make.names(colnames(metadata_sample))
      #Message:
      message <- paste("The column names of the 'metadata_sample' contain special character that where removed.")
      logger::log_info(message)
      message(message)
    }
  }

  ## ------------ Create Results output folder ----------- ##
  if(is.null(save_plot)==FALSE |is.null(save_table)==FALSE){
    Folder <- SavePath(FolderName= "MetaAnalysis",
                                    path=path)
  }

  ###############################################################################################################################################################################################################
  ## ---------- Run prcomp ------------##
  #--- 1. prcomp
  #Get PCs
  PCA.res <- prcomp(data, center = TRUE, scale=Scaling)
  PCA.res_Info <- as.data.frame(PCA.res$x)

  # Extract loadings for each PC
  PCA.res_Loadings <- as.data.frame(PCA.res$rotation)%>%
    tibble::rownames_to_column("FeatureID")

  #--- 2. Merge with demographics
  PCA.res_Info  <- merge(x=metadata_sample%>% tibble::rownames_to_column("UniqueID"),
                         y=PCA.res_Info%>% tibble::rownames_to_column("UniqueID"),
                         by="UniqueID",
                         all.y=TRUE)%>%
    tibble::column_to_rownames("UniqueID")

  #--- 3. convert columns that are not numeric to factor:
  ## Demographics are often non-numerical, categorical explanatory variables, which is often stored as characters, sometimes integers
  PCA.res_Info[sapply(PCA.res_Info, is.character)] <- lapply(PCA.res_Info[sapply(PCA.res_Info, is.character)], as.factor)
  PCA.res_Info[sapply(PCA.res_Info, is.integer)] <- lapply(PCA.res_Info[sapply(PCA.res_Info, is.integer)], as.factor)

  ###############################################################################################################################################################################################################
  ## ---------- STATS ------------##
  ## 1. Anova p.val
  # Iterate through each combination of meta and PC columns
  MetaData <- names(metadata_sample)

  Stat_results <- list()

  for (meta_col in MetaData) {
    for (pc_col in colnames(PCA.res_Info)[grepl("^PC", colnames(PCA.res_Info))]) {
      Formula <- stats::as.formula(paste(pc_col, "~", meta_col, sep=""))# Create a formula for ANOVA --> When constructing the ANOVA formula, it's important to ensure that the response variable (dependent variable) is numeric.
      #pairwiseFormula <- as.formula(paste("pairwise ~" , meta_col, sep=""))

      anova_result <- stats::aov(Formula, data = PCA.res_Info)# Perform ANOVA
      anova_result_tidy <- broom::tidy(anova_result)#broom::tidy --> convert statistics into table
      anova_row <- dplyr::filter(anova_result_tidy, term != "Residuals")  # Exclude Residuals row

      tukey_result <- stats::TukeyHSD(anova_result)# Perform Tukey test
      tukey_result_tidy <- broom::tidy( tukey_result)#broom::tidy --> convert statistics into table

      #lm_result_p.val <- lsmeans::lsmeans((lm(Formula, data = PCA.res_Info)),  pairwiseFormula , adjust=NULL)# adjust=NULL leads to the p-value!
      #lm_result_p.adj <- lsmeans::lsmeans((lm(Formula, data = PCA.res_Info)), pairwiseFormula, adjust="FDR")
      #lm_result <- as.data.frame(lm_result_p.val$contrasts)
      #lm_result$p.adj <- as.data.frame(lm_result_p.adj$contrasts)[,"p.value"]

      # Combine results
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

      # Store ANOVA result
      Stat_results[[paste(pc_col, meta_col, sep="_")]] <- combined_result
    }
  }

  # Add into one DF
  Stat_results <- dplyr::bind_rows(Stat_results)

  #Add explained variance into the table:
  prop_var_ex <- as.data.frame(((PCA.res[["sdev"]])^2/sum((PCA.res[["sdev"]])^2))*100)%>%#To compute the proportion of variance explained by each component in percent, we divide the variance by sum of total variance and multiply by 100(variance=standard deviation ^2)
    tibble::rownames_to_column("PC")%>%
    dplyr::mutate(PC = paste("PC", PC, sep=""))%>%
    dplyr::rename("Explained_Variance"=2)

  Stat_results <- merge(Stat_results, prop_var_ex,  by="PC",all.x=TRUE)

  ###############################################################################################################################################################################################################
  ## ---------- Top/Bottom ------------##
  #Add top/bottom related features to this
  ## Create a data frame for top and bottom features for each PC
  TopBottom_Features <- lapply(2:ncol(PCA.res_Loadings), function(i){
     #Make input
     pc_loadings <- PCA.res_Loadings[, c("FeatureID", names(PCA.res_Loadings)[i])]
     #get top and bottom features
     n_features <- nrow(pc_loadings)
     n_selected <- round(Percentage * n_features)

     top_features <- as.data.frame(head(arrange(pc_loadings, desc(!!sym(names(pc_loadings)[2]))), n_selected)$FeatureID)%>%
       dplyr::rename(!!paste("Features_", "(Top", Percentage, "%)", sep=""):=1)
     bottom_features <- as.data.frame(head(arrange(pc_loadings, !!sym(names(pc_loadings)[2])), n_selected)$FeatureID) %>%
       dplyr::rename(!!paste("Features_", "(Bottom", Percentage, "%)", sep=""):=1)

     #Return
     res <- cbind(data.frame(PC=paste("PC", i, sep = "")), top_features, bottom_features)
     }) %>%
     dplyr::bind_rows(.id = "PC")%>%
     dplyr::mutate(PC = paste("PC", PC, sep=""))


   ## ---------- Final DF 1------------##
   ## Add to results DF
   Stat_results <- merge(Stat_results,
                         TopBottom_Features%>%
                           dplyr::group_by(PC) %>%
                           dplyr::summarise(across(everything(), ~ paste(unique(gsub(", ", "_", .)), collapse = ", "))) %>%
                           dplyr::ungroup(),
                         by="PC",
                         all.x=TRUE)

   ## ---------- DF 2: Metabolites as row names ------------##
   Res_Top <- Stat_results%>%
     dplyr::filter(tukeyHSD_p.adjusted < StatCutoff)%>%
     tidyr::separate_rows(paste("Features_", "(Top", Percentage, "%)", sep=""), sep = ", ")%>% # Separate 'Features (Top 0.1%)'
     dplyr::rename("FeatureID":= paste("Features_", "(Top", Percentage, "%)", sep=""))%>%
     dplyr::select(- paste("Features_", "(Bottom", Percentage, "%)", sep=""))

   Res_Bottom <- Stat_results%>%
     dplyr::filter(tukeyHSD_p.adjusted< StatCutoff)%>%
     tidyr::separate_rows(paste("Features_", "(Bottom", Percentage, "%)", sep=""), sep = ", ")%>%    # Separate 'Features (Bottom 0.1%)'
     dplyr::rename("FeatureID":= paste("Features_", "(Bottom", Percentage, "%)", sep=""))%>%
     dplyr::select(- paste("Features_", "(Top", Percentage, "%)", sep=""))

   Res <- rbind(Res_Top, Res_Bottom)%>%
     dplyr::group_by(FeatureID, term) %>%            # Group by FeatureID and term
     dplyr::summarise(
       PC = paste(unique(PC), collapse = ", "), # Concatenate unique PC entries with commas
       `Sum(Explained_Variance)` = sum(Explained_Variance, na.rm = TRUE)) %>% # Sum Explained_Variance
     dplyr::ungroup()%>%  # Remove previous grouping
     dplyr::group_by(FeatureID) %>%  # Group by FeatureID for MainDriver calculation
     dplyr::mutate(MainDriver = (`Sum(Explained_Variance)` == max(`Sum(Explained_Variance)`))) %>% # Mark TRUE for the highest value
     dplyr::ungroup()  # Remove grouping

   Res_summary <- Res%>%
     dplyr::group_by(FeatureID) %>%
     dplyr::summarise(
       term = paste(term, collapse = ", "),  # Concatenate all terms separated by commas
       `Sum(Explained_Variance)` = paste(`Sum(Explained_Variance)`, collapse = ", "),  # Concatenate all Sum(Explained_Variance) values
       MainDriver = paste(MainDriver, collapse = ", ")) %>%  # Extract the term where MainDriver is TRUE
     dplyr::ungroup()%>%
     dplyr::rowwise() %>%
     dplyr::mutate(  # Extract the term where MainDriver is TRUE
       MainDriver_Term = ifelse("TRUE" %in% strsplit(MainDriver, ", ")[[1]],
                                strsplit(term, ", ")[[1]][which(strsplit(MainDriver, ", ")[[1]] == "TRUE")[1]],
                                NA),
       # Extract the Sum(Explained_Variance) where MainDriver is TRUE
       `MainDriver_Sum(VarianceExplained)` = ifelse("TRUE" %in% strsplit(MainDriver, ", ")[[1]],
                                                    as.numeric(strsplit(`Sum(Explained_Variance)`, ", ")[[1]][which(strsplit(MainDriver, ", ")[[1]] == "TRUE")[1]]),
                                                    NA)) %>%
     dplyr::ungroup()%>%
     dplyr::arrange(desc(`MainDriver_Sum(VarianceExplained)`))

   ###############################################################################################################################################################################################################
   ## ---------- Plot ------------##
   # Plot DF
   Data_Heat <- Stat_results %>%
     dplyr::filter(tukeyHSD_p.adjusted < StatCutoff) %>%  # Filter for significant results
     dplyr::filter(Explained_Variance > VarianceCutoff) %>%  # Exclude Residuals row
     dplyr::distinct(term, PC, .keep_all = TRUE) %>%  # only keep unique term~PC combinations AND STATS
     dplyr::select(term, PC, Explained_Variance) %>%
     tidyr::pivot_wider(names_from = PC, values_from = Explained_Variance) %>%
     tibble::column_to_rownames("term") %>%
     dplyr::mutate_all(~replace(., is.na(.), 0L))

   if(nrow(Data_Heat) > 2L){

     #Plot
     invisible(VizHeatmap(data = Data_Heat,
                         plot_name = paste0("ExplainedVariance-bigger-", VarianceCutoff , "Percent_AND_p.adj-smaller", StatCutoff, sep=""),
                          Scale = "none",
                          save_plot = save_plot,
                          print_plot = print_plot,
                          path = Folder))

   }else{
     message <- paste0("StatCutoff of ", StatCutoff, " and VarianceCutoff of ", VarianceCutoff, " do only return <= 2 cases, hence no heatmap is plotted.")
     logger::log_info("warning: ", message)
     warning(message)
   }



   ###############################################################################################################################################################################################################
   ## ---------- Return ------------##
    # Make list of Output DFs: 1. prcomp results, 2. Loadings result, 3. TopBottom Features, 4. AOV
    ResList <- list(res_prcomp = PCA.res_Info, res_loadings = PCA.res_Loadings, res_TopBottomFeatures = TopBottom_Features, res_aov=Stat_results, res_summary=Res_summary)

    # Save the results
    SaveRes(InputList_DF=ResList,
                         InputList_Plot = NULL,
                         save_table=save_table,
                         save_plot=FALSE,
                         path= Folder,
                         FileName= "MetaAnalysis",
                         core=FALSE,
                         print_plot=FALSE)

    #Return
    invisible(return(ResList))
}


###############################################
### ### ### MetaPriorKnowlegde ### ### ###
###############################################

#' Meta prior-knowledge
#'
#' @param data DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected. includes experimental design and outlier column.
#' @param metadata_sample \emph{Optional: } DF which contains information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames. Column "Conditions" with information about the sample conditions (e.g. "N" and "T" or "Normal" and "Tumor"), can be used for feature filtering and colour coding in the PCA. Column "AnalyticalReplicate" including numerical values, defines technical repetitions of measurements, which will be summarised. Column "BiologicalReplicates" including numerical values. Please use the following names: "Conditions", "Biological_Replicates", "Analytical_Replicates".\strong{Default = NULL}
#' @param metadata_info \emph{Optional: } NULL or vector with column names that should be used, i.e. c("Age", "gender", "Tumour-stage"). \strong{default: NULL}
#' @param save_table \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt". \strong{Default = "csv"}
#' @param path \emph{Optional:} Path to the folder the results should be saved at. \strong{default: NULL}
#'
#' @return DF with prior knowledge based on patient metadata
#'
#' @examples
#' Tissue_Norm <- MetaProViz::ToyData("Tissue_Norm")
#' Res <- MetaProViz::MetaPK(data=Tissue_Norm[,-c(1:13)],
#'                           metadata_sample= Tissue_Norm[,c(2,4:5,12:13)])
#'
#' @keywords prior knowledge, metadata
#'
#' @importFrom magrittr %>%
#' @importFrom tidyr pivot_longer unite
#' @importFrom tibble rownames_to_column
#'
#' @export
#'
MetaPK <- function(data,
                   metadata_sample,
                   metadata_info=NULL,
                   save_table = "csv",
                   path = NULL){

  ## ------------ Create log file ----------- ##
  MetaProViz_Init()

  ### Enrichment analysis-based
  #*we can make pathway file from metadata and use this to run enrichment analysis.
  #*At the moment we can just make the pathways and we need to include a standard fishers exact test.
  #*Merge with function above?
  # *Advantage/disadvantage: Not specific - with annova PC we get granularity e.g. smoker-ExSmoker-PC5, whilst here we only get smoking in generall as parameter

  ################################################################################################################################################################################################
  ## ------------ Check Input files ----------- ##
  # HelperFunction `CheckInput`





  ## ------------ Create Results output folder ----------- ##
  if(is.null(save_table)==FALSE){
    Folder <- SavePath(FolderName= "MetaAnalysis",
                                    path=path)
  }

  ###############################################################################################################################################################################################################
  ## ---------- Create Prior Knowledge file format to perform enrichment analysis ------------##
  # Use the Sample metadata for this:
  if(is.null(metadata_info)==TRUE){
    MetaData <- names(metadata_sample)
    metadata_sample_subset <- metadata_sample%>%
      tibble::rownames_to_column("SampleID")
  }else{
    MetaData <- metadata_info
    metadata_sample_subset <- metadata_sample[, MetaData, drop = FALSE]%>%
      tibble::rownames_to_column("SampleID")
  }

  # Convert into a pathway DF
  Metadata_df <- metadata_sample_subset %>%
    tidyr::pivot_longer(cols = -SampleID, names_to = "ColumnName", values_to = "ColumnEntry")%>%
    tidyr::unite("term", c("ColumnName", "ColumnEntry"), sep = "_")
  Metadata_df$mor <- 1

  #Run ULM using decoupleR (input=PCs from prcomp and pkn=Metadata_df)
  #ULM_res <- decoupleR::run_ulm(mat=as.matrix(PCA.res$x),
  #                              net=Metadata_df,
  #                              .source='term',
  #                              .target='SampleID',
  #                              .mor='mor',
  #                              minsize = 5)

  #Add explained variance into the table:
  #prop_var_ex <- as.data.frame(((PCA.res[["sdev"]])^2/sum((PCA.res[["sdev"]])^2))*100)%>%#To compute the proportion of variance explained by each component in percent, we divide the variance by sum of total variance and multiply by 100(variance=standard deviation ^2)
  #  rownames_to_column("PC")%>%
   # mutate(PC = paste("PC", PC, sep=""))%>%
   # dplyr::rename("Explained_Variance"=2)

  #ULM_res <- merge(ULM_res, prop_var_ex,  by.x="condition",by.y="PC",all.x=TRUE)

  ###############################################################################################################################################################################################################
  ## ---------- Save ------------##
  # Add to results DF
  Res <- list(MetaData_PriorKnowledge = Metadata_df)

}

