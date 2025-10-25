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

#'

###############################################
### ### ### metadata_analysis ### ### ###
###############################################

#' This function performs a PCA analysis on the input data and combines it with the sample metadata to perform an ANOVA test to identify significant differences between the groups.
#'
#' @param data DF where rows are unique samples and columns are features, with
#'     numerical values in columns, and metabolite identifiers as column names. Use
#'     NA for metabolites that were not detected. Includes experimental design and
#'     outlier column.
#' @param metadata_sample \emph{Optional: } DF which contains information about
#'     the samples, which will be combined with your input data based on the
#'     join specification in `by`. Column "Conditions" with information
#'     about the sample conditions (e.g. "N" and "T" or "Normal" and "Tumor"), can
#'     be used for feature filtering and colour coding in the PCA. Column
#'     "AnalyticalReplicate" including numerical values, defines technical
#'     repetitions of measurements, which will be summarised. Column
#'     "BiologicalReplicates" including numerical values. Please use the following
#'     names: "Conditions", "Biological_Replicates",
#'     "Analytical_Replicates".\strong{Default = NULL}
#' @param scaling \emph{Optional: } TRUE or FALSE for whether a data scaling is
#' used \strong{Default = TRUE}
#' @param percentage \emph{Optional: } percentage of top and bottom features to
#' be displayed in the results summary. \strong{Default = 0.1}
#' @param cutoff_stat \emph{Optional: } Cutoff for the adjusted p-value of the
#' ANOVA test for the results summary and on the heatmap. \strong{Default =
#' 0.05}
#' @param cutoff_variance \emph{Optional: } Cutoff for the PCs variance that
#' should be displayed on the heatmap. \strong{Default = 1}
#' @param save_plot \emph{Optional: } Select the file type of output plots.
#' Options are svg, png, pdf. \strong{Default = svg}
#' @param save_table \emph{Optional: } File types for the analysis results are:
#' "csv", "xlsx", "txt". \strong{Default = "csv"}
#' @param print_plot \emph{Optional: } TRUE or FALSE, if TRUE Volcano plot is
#' saved as an overview of the results. \strong{Default = TRUE}
#' @param path \emph{Optional:} Path to the folder the results should be saved
#' at. \strong{default: NULL}
#'
#' @return List of DFs: prcomp results, loadings, top-Bottom features, annova
#' results, results summary
#'
#' @examples
#' data(tissue_norm)
#' Res <- metadata_analysis(
#'     data = tissue_norm[,-c(2:14)]%>%tibble::column_to_rownames("Code"),
#'     metadata_sample = tissue_norm[,c(1,3,5:6,13:14)]%>%tibble::column_to_rownames("Code")
#' )
#'
#' @keywords PCA, annova, metadata
#'
#' @importFrom dplyr filter bind_rows rename mutate ungroup group_by summarise
#' @importFrom dplyr select arrange rowwise mutate_all distinct left_join
#' @importFrom magrittr %>% %<>%
#' @importFrom stats as.formula aov TukeyHSD
#' @importFrom broom tidy
#' @importFrom tidyr separate_rows pivot_wider
#' @importFrom tibble rownames_to_column
#' @importFrom logger log_info
#' @importFrom purrr map_lgl
#'
#' @export
metadata_analysis <- function(data,
                         metadata_sample,
                         #by = NULL,#Join specification between `data` and `metadata_sample`. See the docs of \code{left_join} for details. \strong{Default = NULL}
                         scaling = TRUE,
                         percentage = 0.1,
                         cutoff_stat= 0.05,
                         cutoff_variance=1,
                         save_table = "csv",
                         save_plot = "svg",
                         print_plot= TRUE,
                         path = NULL
                         #SettingInfo= c(MainSeparator = "TISSUE_TYPE), #
                         #enable this parameter in the function --> main
                         #separator: Often a combination of demographics is is
                         #of paricular interest, e.g. comparing "Tumour versus
                         #Normal" for early stage patients and for late stage
                         #patients independently. If this is the case, we can
                         #use the parameter `metadata_info` and provide the
                         #column name of our main separator.

){
  ## ------------ Create log file ----------- ##
  metaproviz_init()

  ################################################################################################################################################################################################
  ## ------------ Check Input files ----------- ##
  # HelperFunction `check_param`
  check_param(data=data,
                          metadata_sample=metadata_sample,
                          metadata_feature=NULL,
                          metadata_info=NULL,
                          save_plot=save_plot,
                          save_table=save_table,
                          core=FALSE,
                          print_plot= print_plot)

  # Specific checks: Check the column names of the demographics --> need to be R usable (no empty spaces, -, etc.)
  if(!is.null(metadata_sample)){
    if(any(grepl('[^[:alnum:]]', colnames(metadata_sample)))){
      #Remove special characters in colnames
      colnames(metadata_sample) <- make.names(colnames(metadata_sample))
      #Message:
      message <- paste("The column names of the 'metadata_sample' contain special character that where removed.")
      log_info(message)
      message(message)
    }
  }

  ## ------------ Create Results output folder ----------- ##
  if(is.null(save_plot)==FALSE |is.null(save_table)==FALSE){
    folder <- save_path(folder_name= "MetadataAnalysis",
                                    path=path)
  }

  ###############################################################################################################################################################################################################
  ## ---------- Run prcomp ------------##
  #--- 1. prcomp
  #Get PCs
  PCA.res <- prcomp(data, center = TRUE, scale=scaling)
  PCA.res_Info <- as.data.frame(PCA.res$x)

  # Extract loadings for each PC
  PCA.res_Loadings <- as.data.frame(PCA.res$rotation) %>%
    rownames_to_column("feature")

  #--- 2. Merge with demographics
  #PCA.res_Info %<>% left_join(metadata_sample, by = by)
  PCA.res_Info  <- merge(x=metadata_sample%>% rownames_to_column("UniqueID"),
                         y=PCA.res_Info%>% rownames_to_column("UniqueID"),
                         by="UniqueID",
                         all.y=TRUE)%>%
    column_to_rownames("UniqueID")

  #--- 3. convert columns that are not numeric to factor:
  ## Demographics are often non-numerical, categorical explanatory variables, which is often stored as characters, sometimes integers
  char_cols <- map_lgl(PCA.res_Info, is.character)
  PCA.res_Info[char_cols] <- lapply(PCA.res_Info[char_cols], as.factor)
  int_cols <- map_lgl(PCA.res_Info, is.integer)
  PCA.res_Info[int_cols] <- lapply(PCA.res_Info[int_cols], as.factor)

  ###############################################################################################################################################################################################################
  ## ---------- STATS ------------##
  ## 1. Anova p.val
  # Iterate through each combination of meta and PC columns
  Metadata <- names(metadata_sample)

  Stat_results <- list()

  for (meta_col in Metadata) {
    for (pc_col in colnames(PCA.res_Info)[grepl("^PC", colnames(PCA.res_Info))]) {
      Formula <- as.formula(paste(pc_col, "~", meta_col, sep=""))# Create a formula for ANOVA --> When constructing the ANOVA formula, it's important to ensure that the response variable (dependent variable) is numeric.
      #pairwiseFormula <- as.formula(paste("pairwise ~" , meta_col, sep=""))

      anova_result <- stats::aov(Formula, data = PCA.res_Info)# Perform ANOVA
      anova_result_tidy <- tidy(anova_result)#tidy --> convert statistics into table
      anova_row <- filter(anova_result_tidy, term != "Residuals")  # Exclude Residuals row

      tukey_result <- TukeyHSD(anova_result)# Perform Tukey test
      tukey_result_tidy <- tidy( tukey_result)#tidy --> convert statistics into table

      #lm_result_p.val <- lsmeans((lm(Formula, data = PCA.res_Info)),  pairwiseFormula , adjust=NULL)# adjust=NULL leads to the p-value!
      #lm_result_p.adj <- lsmeans((lm(Formula, data = PCA.res_Info)), pairwiseFormula, adjust="FDR")
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
  Stat_results <- bind_rows(Stat_results)

  #Add explained variance into the table:
  prop_var_ex <- as.data.frame(((PCA.res[["sdev"]])^2/sum((PCA.res[["sdev"]])^2))*100)%>%#to compute the proportion of variance explained by each component in percent, we divide the variance by sum of total variance and multiply by 100(variance=standard deviation ^2)
    rownames_to_column("PC")%>%
    mutate(PC = paste("PC", PC, sep=""))%>%
    rename("Explained_Variance"=2)

  Stat_results <- merge(Stat_results, prop_var_ex,  by="PC",all.x=TRUE)

  ###############################################################################################################################################################################################################
  ## ---------- top/Bottom ------------##
  #Add top/bottom related features to this
  ## Create a data frame for top and bottom features for each PC
  topBottom_Features <- lapply(2:ncol(PCA.res_Loadings), function(i){
     #Make input
     pc_loadings <- PCA.res_Loadings[, c("feature", names(PCA.res_Loadings)[i])]
     #get top and bottom features
     n_features <- nrow(pc_loadings)
     n_selected <- round(percentage * n_features)

     top_features <- as.data.frame(head(arrange(pc_loadings, desc(!!sym(names(pc_loadings)[2]))), n_selected)$feature)%>%
       rename(!!paste("Features_", "(top", percentage, "%)", sep=""):=1)
     bottom_features <- as.data.frame(head(arrange(pc_loadings, !!sym(names(pc_loadings)[2])), n_selected)$feature) %>%
       rename(!!paste("Features_", "(Bottom", percentage, "%)", sep=""):=1)

     #Return
     res <- cbind(data.frame(PC=paste("PC", i, sep = "")), top_features, bottom_features)
     }) %>%
     bind_rows(.id = "PC")%>%
     mutate(PC = paste("PC", PC, sep=""))


   ## ---------- Final DF 1------------##
   ## Add to results DF
   Stat_results <- merge(Stat_results,
                         topBottom_Features%>%
                           group_by(PC) %>%
                           summarise(across(everything(), ~ paste(unique(gsub(", ", "_", .)), collapse = ", "))) %>%
                           ungroup(),
                         by="PC",
                         all.x=TRUE)

   ## ---------- DF 2: Metabolites as row names ------------##
   Res_top <- Stat_results%>%
     filter(tukeyHSD_p.adjusted < cutoff_stat)%>%
     separate_rows(paste("Features_", "(top", percentage, "%)", sep=""), sep = ", ")%>% # Separate 'Features (top 0.1%)'
     rename("feature":= paste("Features_", "(top", percentage, "%)", sep=""))%>%
     select(- paste("Features_", "(Bottom", percentage, "%)", sep=""))

   Res_Bottom <- Stat_results%>%
     filter(tukeyHSD_p.adjusted< cutoff_stat)%>%
     separate_rows(paste("Features_", "(Bottom", percentage, "%)", sep=""), sep = ", ")%>%    # Separate 'Features (Bottom 0.1%)'
     rename("feature":= paste("Features_", "(Bottom", percentage, "%)", sep=""))%>%
     select(- paste("Features_", "(top", percentage, "%)", sep=""))

   Res <- rbind(Res_top, Res_Bottom)%>%
     group_by(feature, term) %>%            # Group by feature and term
     summarise(
       PC = paste(unique(PC), collapse = ", "), # Concatenate unique PC entries with commas
       `Sum(Explained_Variance)` = sum(Explained_Variance, na.rm = TRUE)) %>% # Sum Explained_Variance
     ungroup()%>%  # Remove previous grouping
     group_by(feature) %>%  # Group by feature for MainDriver calculation
     mutate(MainDriver = (`Sum(Explained_Variance)` == max(`Sum(Explained_Variance)`))) %>% # Mark TRUE for the highest value
     ungroup()  # Remove grouping

   Res_summary <- Res%>%
     group_by(feature) %>%
     summarise(
       term = paste(term, collapse = ", "),  # Concatenate all terms separated by commas
       `Sum(Explained_Variance)` = paste(`Sum(Explained_Variance)`, collapse = ", "),  # Concatenate all Sum(Explained_Variance) values
       MainDriver = paste(MainDriver, collapse = ", ")) %>%  # Extract the term where MainDriver is TRUE
     ungroup()%>%
     rowwise() %>%
     mutate(  # Extract the term where MainDriver is TRUE
       MainDriver_Term = ifelse("TRUE" %in% strsplit(MainDriver, ", ")[[1]],
                                strsplit(term, ", ")[[1]][which(strsplit(MainDriver, ", ")[[1]] == "TRUE")[1]],
                                NA),
       # Extract the Sum(Explained_Variance) where MainDriver is TRUE
       `MainDriver_Sum(VarianceExplained)` = ifelse("TRUE" %in% strsplit(MainDriver, ", ")[[1]],
                                                    as.numeric(strsplit(`Sum(Explained_Variance)`, ", ")[[1]][which(strsplit(MainDriver, ", ")[[1]] == "TRUE")[1]]),
                                                    NA)) %>%
     ungroup()%>%
     arrange(desc(`MainDriver_Sum(VarianceExplained)`))

   ###############################################################################################################################################################################################################
   ## ---------- Plot ------------##
   # Plot DF
   data_Heat <- Stat_results %>%
     filter(tukeyHSD_p.adjusted < cutoff_stat) %>%  # Filter for significant results
     filter(Explained_Variance > cutoff_variance) %>%  # Exclude Residuals row
     distinct(term, PC, .keep_all = TRUE) %>%  # only keep unique term~PC combinations AND STATS
     select(term, PC, Explained_Variance) %>%
     pivot_wider(names_from = PC, values_from = Explained_Variance) %>%
     column_to_rownames("term") %>%
     mutate_all(~replace(., is.na(.), 0L))

   if(nrow(data_Heat) > 2L){

     #Plot
     invisible(viz_heatmap(data = data_Heat,
                         plot_name = paste0("ExplainedVariance-bigger-", cutoff_variance , "Percent_AND_p.adj-smaller", cutoff_stat, sep=""),
                          scale = "none",
                          save_plot = save_plot,
                          print_plot = print_plot,
                          path = folder))

   }else{
     message <- paste0("cutoff_stat of ", cutoff_stat, " and cutoff_variance of ", cutoff_variance, " do only return <= 2 cases, hence no heatmap is plotted.")
     log_info("warning: ", message)
     warning(message)
   }



   ###############################################################################################################################################################################################################
   ## ---------- Return ------------##
    # Make list of Output DFs: 1. prcomp results, 2. Loadings result, 3. topBottom Features, 4. aov
    ResList <- list(res_prcomp = PCA.res_Info, res_loadings = PCA.res_Loadings, res_topBottomFeatures = topBottom_Features, res_aov=Stat_results, res_summary=Res_summary)

    # Save the results
    save_res(inputlist_df=ResList,
                         inputlist_plot = NULL,
                         save_table=save_table,
                         save_plot=FALSE,
                         path= folder,
                         file_name= "metadata_analysis",
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
#' data(tissue_norm)
#' Tissue_Norm <- tissue_norm %>%tibble::column_to_rownames("Code")
#' Res <- meta_pk(data=Tissue_Norm[,-c(1:13)],
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
meta_pk <- function(data,
                   metadata_sample,
                   metadata_info=NULL,
                   save_table = "csv",
                   path = NULL){

  ## ------------ Create log file ----------- ##
  metaproviz_init()

  ### Enrichment analysis-based
  #*we can make pathway file from metadata and use this to run enrichment analysis.
  #*At the moment we can just make the pathways and we need to include a standard fishers exact test.
  #*Merge with function above?
  # *Advantage/disadvantage: Not specific - with annova PC we get granularity e.g. smoker-ExSmoker-PC5, whilst here we only get smoking in generall as parameter

  ################################################################################################################################################################################################
  ## ------------ Check Input files ----------- ##
  # HelperFunction `check_param`





  ## ------------ Create Results output folder ----------- ##
  if(is.null(save_table)==FALSE){
    folder <- save_path(folder_name= "MetadataAnalysis",
                                    path=path)
  }

  ###############################################################################################################################################################################################################
  ## ---------- Create Prior Knowledge file format to perform enrichment analysis ------------##
  # Use the Sample metadata for this:
  if(is.null(metadata_info)==TRUE){
    Metadata <- names(metadata_sample)
    metadata_sample_subset <- metadata_sample%>%
      rownames_to_column("SampleID")
  }else{
    Metadata <- metadata_info
    metadata_sample_subset <- metadata_sample[, Metadata, drop = FALSE]%>%
      rownames_to_column("SampleID")
  }

  # Convert into a pathway DF
  Metadata_df <- metadata_sample_subset %>%
    pivot_longer(cols = -SampleID, names_to = "ColumnName", values_to = "ColumnEntry")%>%
    unite("term", c("ColumnName", "ColumnEntry"), sep = "_")
  Metadata_df$mor <- 1

  #Run ULM using decoupleR (input=PCs from prcomp and pkn=Metadata_df)
  #ULM_res <- run_ulm(mat=as.matrix(PCA.res$x),
  #                              net=Metadata_df,
  #                              .source='term',
  #                              .target='SampleID',
  #                              .mor='mor',
  #                              minsize = 5)

  #Add explained variance into the table:
  #prop_var_ex <- as.data.frame(((PCA.res[["sdev"]])^2/sum((PCA.res[["sdev"]])^2))*100)%>%#to compute the proportion of variance explained by each component in percent, we divide the variance by sum of total variance and multiply by 100(variance=standard deviation ^2)
  #  rownames_to_column("PC")%>%
   # mutate(PC = paste("PC", PC, sep=""))%>%
   # rename("Explained_Variance"=2)

  #ULM_res <- merge(ULM_res, prop_var_ex,  by.x="condition",by.y="PC",all.x=TRUE)

  ###############################################################################################################################################################################################################
  ## ---------- Save ------------##
  # Add to results DF
  Res <- list(Metadata_prior_knowledge = Metadata_df)

}


