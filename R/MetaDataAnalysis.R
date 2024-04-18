## ---------------------------
##
## Script name:
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
#' This script allows you to perform


###############################################
### ### ### MetaData-PC correlation ### ### ###
###############################################

#' @param Input_SettingsInfo \emph{Optional: } NULL or Named vector including at least one of those three information : c(color="ColumnName_Plot_SettingsFile", shape= "ColumnName_Plot_SettingsFile"). \strong{Default = NULL}
#' @param Input_SettingsFile_Feature \emph{Optional: } DF which contains information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames. Column "Conditions" with information about the sample conditions (e.g. "N" and "T" or "Normal" and "Tumor"), can be used for feature filtering and colour coding in the PCA. Column "AnalyticalReplicate" including numerical values, defines technical repetitions of measurements, which will be summarised. Column "BiologicalReplicates" including numerical values. Please use the following names: "Conditions", "Biological_Replicates", "Analytical_Replicates".\strong{Default = NULL}
#' @param Input_SettingsInfo_Feature \emph{Optional: } DF which contains information about the features (=metabolites)
#' @param Input_data DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected. includes experimental design and outlier column.
#' @param Scaling \emph{Optional: } TRUE or FALSE for whether a data scaling is used \strong{Default = TRUE}
#' @param OutputPlotName \emph{Optional: } String which is added to the output files of the PCA \strong{Default = ""}
#' @param Save_as_Results \emph{Optional: } Select the file type of output plots. Options are svg, png, pdf or NULL. \strong{Default = svg}
#' @param Folder_Name {Optional:} String which is added to the resulting folder name \strong(Default = NULL)
#'
#' @keywords PCA
#' @export

Meta_Corr <- function(Input_SettingsInfo= NULL,
                   Input_SettingsFile,
                   Input_data,
                   InputSettingsFile_Feature,
                   Input_SettingsInfo_Feature,
                   Scaling = TRUE,
                   percentage = 0.1,#Top/Bottom features
                   OutputName= '',
                   Save_as_Results = "svg",
                   Folder_Name = NULL,
                   path = NULL

){

  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  suppressMessages(library(tidyverse))

  ## ------------ Check Input files ----------- ##
  # 1. The input data:
  if(any(duplicated(row.names(Input_data))) ==  TRUE){# Is the "Input_data" has unique IDs as row names and numeric values in columns?
    stop("Duplicated row.names of Input_data, whilst row.names must be unique")
  }else{
    Input_data_m<- Input_data
    Test_num <- apply(Input_data_m, 2, function(x) is.numeric(x))
    if((any(Test_num) ==  FALSE) ==  TRUE){
      stop("Input_data needs to be of class numeric in all columns.")
    }
  }

  # 2. The Plot_settings: Plot_Settings, Plot_SettingInfo and Plot_SettingFile
  if(is.null(Input_SettingsInfo)==TRUE){
    message("You have not provided Input_SettingsInfo, hence each column in DF Input_SettingsFile will be included.")
  }
  if(is.vector(Input_SettingsInfo)==FALSE & is.null(Input_SettingsInfo)==FALSE){
    stop("Input_SettingsInfo must be named vector or NULL.")
  }
  if(is.null(Input_SettingsInfo)==FALSE & is.null(Input_SettingsFile)==FALSE){
    if(all(Input_SettingsInfo %in% colnames(Input_SettingsFile))==FALSE){
      warning("Some or all columns selected with Input_SettingsInfo do not exist in Input_SettingsFile.")
    }
  }
  #Add for feature info and rename above for Sample!


  #3. Check other plot-specific parameters:



  #save
  if (!is.null(Save_as_Results)) {
    Save_as_Results_options <- c("svg","pdf", "png")
    if(Save_as_Results %in% Save_as_Results_options == FALSE){
      stop("Check input. The selected Save_as_Plot option is not valid. Please select one of the folowing: ",paste(Save_as_Results_options,collapse = ", "),"." )
    }
  }

  ## ------------ Create Output folders ----------- ##
  if (!is.null(Save_as_Results)) {
    if(is.null(Folder_Name)){
      name <- paste("MetaProViz_Results",Sys.Date(),sep = "_" )
    }else{
      if(grepl('[^[:alnum:]]', Folder_Name)){
        stop("The 'Folder_Name' must not contain any special character.")
      }else{
        name <- paste("MetaProViz_Results",Sys.Date(),Folder_Name,sep = "_" )
      }
    }
    WorkD <- getwd()
    Results_folder <- file.path(WorkD, name)
    if (!dir.exists(Results_folder)) {dir.create(Results_folder)} # Make Results folder
    Results_folder_PCA_folder = file.path(Results_folder, "Meta_Corr")  # This searches for a folder called "" within the "Results" folder in the current working directory and if its not found it creates one
    if (!dir.exists(Results_folder_PCA_folder)) {dir.create(Results_folder_PCA_folder)}  # check and create folder
  }

  ############################################################################################################
  ## ---------- Run prcomp ------------##
  #--- 1. prcomp
  #Get PCs
  PCA.res <- prcomp(Input_data, center = TRUE, scale=Scaling)
  PCA.res_Info <- as.data.frame(PCA.res$x)

  # Extract loadings for each PC
  PCA.res_Loadings <- as.data.frame(PCA.res$rotation)%>%
    rownames_to_column("FeatureID")
  if(is.null(InputSettingsFile_Feature)==FALSE){
    PCA.res_Loadings_Info <- merge(x=InputSettingsFile_Feature, y=PCA.res_Loadings, by.x=Input_SettingsInfo_Feature[["FeatureID"]], by.y="FeatureID", all.y=TRUE)
  }

  #--- 2. Merge with demographics
  if(is.null(Input_SettingsInfo)==TRUE){
    MetaData <- names(Input_SettingsFile)
    PCA.res_Info  <- merge(x=Input_SettingsFile%>%rownames_to_column("UniqueID") , y=PCA.res_Info%>%rownames_to_column("UniqueID"), by="UniqueID", all.y=TRUE)%>%
      column_to_rownames("UniqueID")
  }else{
    MetaData <- Input_SettingsInfo
    Input_SettingsFile_subset <- Input_SettingsFile[, MetaData, drop = FALSE]
    PCA.res_Info  <- merge(x=Input_SettingsFile_subset%>%rownames_to_column("UniqueID") , y=PCA.res_Info%>%rownames_to_column("UniqueID"), by="UniqueID", all.y=TRUE)%>%
      column_to_rownames("UniqueID")
  }

  #--- 3. convert columns that are not numeric to factor:
  ## Demographics are often non-numerical, categorical explanatory variables, which is often stored as characters, sometimes integers
  PCA.res_Info[sapply(PCA.res_Info, is.character)] <- lapply(PCA.res_Info[sapply(PCA.res_Info, is.character)], as.factor)
  PCA.res_Info[sapply(PCA.res_Info, is.integer)] <- lapply(PCA.res_Info[sapply(PCA.res_Info, is.integer)], as.factor)

  ## ---------- STATS ------------##
  ## 1. Anova p.val
  # Iterate through each combination of meta and PC columns
  Stat_results <- list()

  for (meta_col in MetaData) {
    for (pc_col in colnames(PCA.res_Info)[grepl("^PC", colnames(PCA.res_Info))]) {
      Formula <- as.formula(paste(pc_col, "~", meta_col, sep=""))# Create a formula for ANOVA --> When constructing the ANOVA formula, it's important to ensure that the response variable (dependent variable) is numeric.
      #pairwiseFormula <- as.formula(paste("pairwise ~" , meta_col, sep=""))

      anova_result <- aov(Formula, data = PCA.res_Info)# Perform ANOVA
      anova_result_tidy <- broom::tidy(anova_result)#broom::tidy --> convert statistics into table
      anova_row <- filter(anova_result_tidy, term != "Residuals")  # Exclude Residuals row

      tukey_result <- TukeyHSD(anova_result)# Perform Tukey test
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
    rownames_to_column("PC")%>%
    mutate(PC = paste("PC", PC, sep=""))%>%
    dplyr::rename("Explained_Variance"=2)

  Stat_results <- merge(Stat_results, prop_var_ex,  by="PC",all.x=TRUE)

  #Add top/bottom related features to this









  ## ---------- ULM ------------##
  # Use the Sample metadata for this:
  if(is.null(Input_SettingsInfo)==TRUE){
    MetaData <- names(Input_SettingsFile)
    Input_SettingsFile <- Input_SettingsFile%>%
      rownames_to_column("SampleID")
  }else{
    MetaData <- Input_SettingsInfo
    Input_SettingsFile_subset <- Input_SettingsFile[, MetaData, drop = FALSE]%>%
      rownames_to_column("SampleID")
  }

  # Convert into a pathway DF
  Metadata_df <- Input_SettingsFile_subset %>%
    pivot_longer(cols = -SampleID, names_to = "ColumnName", values_to = "ColumnEntry")%>%
    unite("term", c("ColumnName", "ColumnEntry"), sep = "_")
  Metadata_df$mor <- 1

  #Run ULM using decoupleR (input=PCs from prcomp and pkn=Metadata_df)
  ULM_res <- decoupleR::run_ulm(mat=as.matrix(PCA.res$x),
                     net=Metadata_df,
                     .source='term',
                     .target='SampleID',
                     .mor='mor',
                     minsize = 5)

  #Add explained variance into the table:
  prop_var_ex <- as.data.frame(((PCA.res[["sdev"]])^2/sum((PCA.res[["sdev"]])^2))*100)%>%#To compute the proportion of variance explained by each component in percent, we divide the variance by sum of total variance and multiply by 100(variance=standard deviation ^2)
    rownames_to_column("PC")%>%
    mutate(PC = paste("PC", PC, sep=""))%>%
    dplyr::rename("Explained_Variance"=2)

  ULM_res <- merge(ULM_res, prop_var_ex,  by.x="condition",by.y="PC",all.x=TRUE)

  #Add top/bottom related features (= e.g. metabolites) to this
  ## Create a data frame for top and bottom features for each PC
  TopBottom_Features <- lapply(2:ncol(PCA.res_Loadings), function(i){
     #Make input
     pc_loadings <- PCA.res_Loadings[, c("FeatureID", names(PCA.res_Loadings)[i])]
     #get top and bottom features
     n_features <- nrow(pc_loadings)
     n_selected <- round(percentage * n_features)

     top_features <- as.data.frame(head(arrange(pc_loadings, desc(!!sym(names(pc_loadings)[2]))), n_selected)$FeatureID)%>%
       dplyr::rename("FeatureID"=1)
     bottom_features <- as.data.frame(head(arrange(pc_loadings, !!sym(names(pc_loadings)[2])), n_selected)$FeatureID) %>%
       dplyr::rename("FeatureID"=1)

     #Add other ID types if applicable
     if(is.null(InputSettingsFile_Feature)==FALSE){
       top_features <- merge(x=InputSettingsFile_Feature, y=top_features, by.x=Input_SettingsInfo_Feature[["FeatureID"]], by.y="FeatureID", all.y=TRUE)%>%
         summarise_all(~toString(unique(.)))%>%
         rename_all(~paste(paste("Top", percentage*100, "%_", sep=""), ., sep = ""))

       bottom_features <- merge(x=InputSettingsFile_Feature, y=bottom_features, by.x=Input_SettingsInfo_Feature[["FeatureID"]], by.y="FeatureID", all.y=TRUE)%>%
         summarise_all(~toString(unique(.)))%>%
         rename_all(~paste(paste("Bottom", percentage*100, "%_", sep=""), ., sep = ""))
     }

    #Return
     res <- cbind(data.frame(PC=paste("PC", i, sep = "")), top_features, bottom_features)
     return(res)
     }) %>%
     bind_rows(.id = "PC")%>%
     mutate(PC = paste("PC", PC, sep=""))

   ## Add to results DF
   ULM_res <- merge(ULM_res, TopBottom_Features,  by.x="condition",by.y="PC",all.x=TRUE)
   Stat_results <- merge(Stat_results, TopBottom_Features,  by="PC",all.x=TRUE)

    ## ---------- Heatmap ------------##
    #x= source (=Demographics = ) y= condition (=PCs=Factors)
    #Make outside of this?


    ## ---------- Return ------------##
    # Make list of Output DFs: 1. prcomp results, 2. Loadings result, 3. Aov, 4. ULM
    output_list <- list()
    output_dfs <- list(res_prcomp = PCA.res_Info, res_loadings = PCA.res_Loadings_Info, res_ULM = ULM_res, res_aov=Stat_results)

    if(is.null(Save_as_Results)==FALSE){
      if(Save_as_Results == "xlsx"){
        writexl::write_xlsx(output_dfs, paste0(Results_folder_Preprocessing_folder, "/MetaData_Analysis", OutputName, ".xlsx", sep = ""))   # Save the DMA results table
        }
      if(Save_as_Results == "csv"){
        #for(i in output_dfs){
        #  write.csv(DMA_Output,paste0(Results_folder_Conditions,paste0(Results_folder_Preprocessing_folder, "/MetaData_Analysis", OutputName, ".csv", sep = ""))
        #}
        }
      if(Save_as_Results == "txt"){
          txtDMA <- file.path(Results_folder_Conditions,paste0("DMA_Output_",toString(numerator),"_vs_",toString(denominator), OutputName, ".txt"))
          write.table(DMA_Output,txtDMA, col.names = TRUE, row.names = FALSE) # save the DMA result DF
    }
    }

}
