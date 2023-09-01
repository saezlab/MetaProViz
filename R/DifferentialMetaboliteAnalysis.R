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

#' This script allows you to perform differential metabolite analysis to obtain a Log2FC, pval, padj and tval comparing two conditions.
#'
#' @param Input_data DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected. includes experimental Input_SettingsFile and outlier column.
#' @param Input_SettingsFile DF which contains information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames. Column "Conditions" with information about the sample conditions (e.g. "N" and "T" or "Normal" and "Tumor"), can be used for feature filtering and colour coding in the PCA. Column "AnalyticalReplicate" including numerical values, defines technical repetitions of measurements, which will be summarised. Column "BiologicalReplicates" including numerical values. Please use the following names: "Conditions", "Biological_Replicates", "Analytical_Replicates".
#' @param Conditon1 Input needs to contain a column named "Condition" including the numerator that will be compared to denominator, e.g. "KO".
#' @param Conditon2 Input needs to contain a column named "Condition" including denominator that is compared to numerator, e.g. "WT".
#' @param STAT_pval \emph{Optional: } String which contains an abbreviation of the selected test to calculate p.value (t.test or wilcox.test) \strong{"t-test"}
#' @param STAT_padj \emph{Optional: } String which contains an abbreviation of the selected p.adjusted test for p.value correction for multiple Hypothesis testing. Search: ?p.adjust for more methods:"BH", "fdr", "bonferroni", "holm", etc.\strong{"fdr"}
#' @param OutputName String which is added to the output files of the DMA.
#' @param Input_Pathways \emph{Optional: } DF which contains the pathway information for each metabolite. \strong{NULL}
#' @param CoRe \emph{Optional: } TRUE or FALSE for whether a Consumption/Release  input is used \strong{FALSE}
#' @param plot \emph{Optional: } TRUE or FALSE, if TRUE Volcano plot is saved as an overview of the results. \strong{TRUE}
#' @param Save_as_Plot \emph{Optional: } Select the file type of output plots. Options are svg, png, pdf. \strong{Default = svg}
#' @param Save_as_Results \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt" \strong{default: "xlsx"}
#'
#' @keywords Differential Metabolite Analysis, Multiple Hypothesis testing, Normality testing
#' @export


########################################################
### ### ### Differential Metabolite Analysis ### ### ###
########################################################

DMA <-function(Input_data,
               Input_SettingsFile,
               Input_SettingsInfo = c(conditions="Conditions", numerator = NULL, denumerator = NULL),
               STAT_pval ="kruskal.test",
               STAT_padj="fdr",
               Input_Pathways = NULL,
               OutputName='',
               CoRe=FALSE,
               Plot = TRUE,
               Save_as_Plot = "svg",
               Save_as_Results = "xlsx" # txt or csv
){

  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse", "gtools", "EnhancedVolcano")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0){
    install.packages(new.packages)
  }
  suppressMessages(library(tidyverse))

  ################################################################################################################################################################################################
  ## ------------ Check Input files ----------- ##
  #1. Input_data and Conditions

  if(any(duplicated(row.names(Input_data)))==TRUE){
    stop("Duplicated row.names of Input_data, whilst row.names must be unique")
  } else{
    Test_num <- apply(Input_data, 2, function(x) is.numeric(x))
    if((any(Test_num) ==  FALSE) ==  TRUE){
      stop("Input_data needs to be of class numeric")
    } else{
      Test_match <- merge(Input_SettingsFile, Input_data, by.x = "row.names",by.y = "row.names", all =  FALSE) # Do the unique IDs of the "Input_data" match the row names of the "Input_SettingsFile"?
      if(nrow(Test_match) ==  0){
        stop("row.names Input_data need to match row.names Input_SettingsFile.")
      } else{
        Input_data <- Input_data
      }
    }
  }



  if(Input_SettingsInfo[["conditions"]] %in% colnames(Input_SettingsFile)== FALSE){
    stop("The ",Input_SettingsInfo[["conditions"]], " column selected as Conditions in Input_SettingsInfo was not found in Input_SettingsFile. Please check your input.")
  }else{# if true rename to Conditions
    Input_SettingsFile<- Input_SettingsFile%>%
      dplyr::rename("Conditions"= paste(Input_SettingsInfo[["conditions"]]) )
  }

  if ("denominator" %in% names(Input_SettingsInfo)==FALSE  & "numerator" %in% names(Input_SettingsInfo) ==FALSE ){
    # Anova all-vs-all
    message("No conditions were specified as numerator or denumerator. Performing Anova all against all.")
    conditions =Input_SettingsFile$Conditions
    denominator <-unique( Input_SettingsFile$Conditions)
    numerator <-unique( Input_SettingsFile$Conditions)
    ANOVA = TRUE
    All = TRUE
    # Generate all pairwise combinations
    comparisons <- combn( unique(conditions), 2) %>% as.matrix()
  }else if ("denominator" %in% names(Input_SettingsInfo)==FALSE  & "numerator" %in% names(Input_SettingsInfo) ==TRUE ){
    stop("Check input. The selected denominator option is empty while ",paste(Input_SettingsInfo[["numerator"]])," has been selected as a numerator. Please add a denminator for 1-vs-1 comparison or remove the numerator for all-vs-all comparison." )
  }else if ("denominator" %in% names(Input_SettingsInfo)==TRUE  & "numerator" %in% names(Input_SettingsInfo) ==FALSE ){
    # Anova all-vs-one
    ANOVA = TRUE
  }else if ("denominator" %in% names(Input_SettingsInfo)==TRUE  & "numerator" %in% names(Input_SettingsInfo) ==TRUE ){
    # one-vs-one
    if(Input_SettingsInfo[["denominator"]] %in% Input_SettingsFile$Conditions  == FALSE){
      stop("The ",Input_SettingsInfo[["denominator"]], " column selected as denominator in Input_SettingsInfo was not found in Input_SettingsFile. Please check your input.")
    }else{
      denominator <- Input_SettingsInfo[["denominator"]]
    }
    if(Input_SettingsInfo[["numerator"]] %in% Input_SettingsFile$Conditions  == FALSE){
      stop("The ",Input_SettingsInfo[["numerator"]], " column selected as numerator in Input_SettingsInfo was not found in Input_SettingsFile. Please check your input.")
    }else{
      numerator <- Input_SettingsInfo[["numerator"]]
    }
    comparisons <- matrix(c(denominator, numerator))
    ANOVA = FALSE
  }

  # If anova is used change the column names
  if(ANOVA == TRUE){
    for (i in 1:length(colnames(Input_data))){
      metabolite_name <- colnames(Input_data)[i]
      metabolite_name <- gsub("-", "", metabolite_name)
      metabolite_name <- gsub("/", "", metabolite_name)
      metabolite_name <- gsub(" ", "", metabolite_name)
      metabolite_name <- gsub("\\*", "", metabolite_name)
      metabolite_name <- gsub("\\+", "", metabolite_name)
      metabolite_name <- gsub(",", "", metabolite_name)
      metabolite_name <- gsub("\\(", "", metabolite_name)
      metabolite_name <- gsub("\\)", "", metabolite_name)
      colnames(Input_data)[i] <- metabolite_name
    }
  }

  Num <- Input_data %>%
    filter(Input_SettingsFile$Conditions %in% numerator) %>%
    select_if(is.numeric)#only keep numeric columns with metabolite values
  Denom <- Input_data %>%
    filter(Input_SettingsFile$Conditions %in% denominator) %>%
    select_if(is.numeric)

  if(nrow(Num)==1){
    stop("There is only one sample available for ", numerator, ", so no statistical test can be performed.")
  } else if(nrow(Denom)==1){
    stop("There is only one sample available for ", denominator, ", so no statistical test can be performed.")
  }else if(nrow(Num)==0){
    stop("There is no sample available for ", numerator, ".")
  }else if(nrow(Denom)==0){
    stop("There is no sample available for ", denominator, ".")
  }

  #2. General parameters
  if(ANOVA==FALSE){
    STAT_pval_options <- c("t.test", "wilcox.test","chisq.test", "cor.test")
    if(STAT_pval %in% STAT_pval_options == FALSE){
      stop("Check input. The selected STAT_pval option for Hypothesis testing is not valid. Please select one of the folowwing: ",paste(STAT_pval_options,collapse = ", "),"." )
    }
  }else{
    STAT_pval_options <- c("anova", "kruskal.test")
    if(STAT_pval %in% STAT_pval_options == FALSE){
      stop("Check input. The selected STAT_pval option for Hypothesis testing is not valid. Multiple comparison is selected. Please select one of the folowwing: ",paste(STAT_pval_options,collapse = ", "),"." )
    }
  }
  STAT_padj_options <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
  if(STAT_padj %in% STAT_padj_options == FALSE){
    stop("Check input. The selected STAT_padj option for multiple Hypothesis testing correction is not valid. Please select one of the folowwing: ",paste(STAT_padj_options,collapse = ", "),"." )
  }
  if(is.logical(CoRe) == FALSE){
    stop("Check input. The CoRe value should be either =TRUE for analysis of Consuption/Release experiment or =FALSE if not.")
  }
  if(is.logical(Plot) == FALSE){
    stop("Check input. The plot value should be either =TRUE if a Volcano plot presenting the DMA results is to be exported or =FALSE if not.")
  }
  Save_as_Plot_options <- c("svg","pdf","png")
  if(Save_as_Plot %in% Save_as_Plot_options == FALSE){
    stop("Check input. The selected Save_as_Plot option is not valid. Please select one of the folowwing: ",paste(Save_as_Plot_options,collapse = ", "),"." )
  }
  Save_as_Results_options <- c("txt","csv", "xlsx" )
  if(Save_as_Results %in% Save_as_Results_options == FALSE){
    stop("Check input. The selected Save_as_Results option is not valid. Please select one of the folowwing: ",paste(Save_as_Results_options,collapse = ", "),"." )
  }


  #3. Input_Pathways
  if(is.null(Input_Pathways) == FALSE){
    if('Metabolite' %in% colnames(Input_Pathways) == FALSE){
      warning("The provided file Input_Pathways must have 2 columns named: `Metabolite` and `Pathway`.")
    }else if('Pathway' %in% colnames(Input_Pathways) == FALSE){
      warning("The provided file Input_Pathways must have 2 columns named: `Metabolite` and `Pathway`.")
    }
  }




  ## ------------ Check Missingness ------------- ##
  # If missing value imputation has not been performed the input data will most likely contain NA or 0 values for some metabolites, which will lead to Log2FC = NA.
  # Here we will check how many metabolites this affects in Num and Denom, and weather all replicates of a metabolite are affected.
  Num_Miss <- replace(Num, Num==0, NA)
  Num_Miss <- Num_Miss[, (colSums(is.na(Num_Miss)) > 0), drop = FALSE]

  Denom_Miss <- replace(Denom, Denom==0, NA)
  Denom_Miss <- Denom_Miss[, (colSums(is.na(Denom_Miss)) > 0), drop = FALSE]

  if((ncol(Num_Miss)>0 & ncol(Denom_Miss)==0)){
    message("In `numerator` ",paste0(toString(numerator)), ", NA/0 values exist in ", ncol(Num_Miss), " Metabolite(s): ", paste0(colnames(Num_Miss), collapse = ", "), ". Those metabolite(s) will return p.val= NA, p.adj.= NA, t.val= NA. The Log2FC = Inf, if all replicates are 0/NA.")
    Metabolites_Miss <- colnames(Num_Miss)
  } else if(ncol(Num_Miss)==0 & ncol(Denom_Miss)>0){
    message("In `denominator` ",paste0(toString(denominator)), ", NA/0 values exist in ", ncol(Denom_Miss), " Metabolite(s): ", paste0(colnames(Denom_Miss), collapse = ", "), ". Those metabolite(s) will return p.val= NA, p.adj.= NA, t.val= NA. The Log2FC = Inf, if all replicates are 0/NA.")
    Metabolites_Miss <- colnames(Denom_Miss)
  } else if(ncol(Num_Miss)>0 & ncol(Denom_Miss)>0){
    message("In `numerator` ",paste0(toString(numerator)), ", NA/0 values exist in ", ncol(Num_Miss), " Metabolite(s): ", paste0(colnames(Num_Miss), collapse = ", "), " and in `denominator`",paste0(toString(denominator)), " ",ncol(Denom_Miss), " Metabolite(s): ", paste0(colnames(Denom_Miss), collapse = ", "),". Those metabolite(s) will return p.val= NA, p.adj.= NA, t.val= NA. The Log2FC = Inf, if all replicates are 0/NA.")
    Metabolites_Miss <- c(colnames(Num_Miss), colnames(Denom_Miss))
    Metabolites_Miss <- unique(Metabolites_Miss)
  } else{
    message("There are no NA/0 values")
    Metabolites_Miss <- c(colnames(Num_Miss), colnames(Denom_Miss))
    Metabolites_Miss <- unique(Metabolites_Miss)
  }

  ## ------------ Check data normality and statistical test chosen ----------- ##
  # Before Hypothesis testing, we have to decide whether to use a parametric or a non parametric test. we can test the data normality using the Shapiro test.

  # 1. Load the data and perform the shapiro.test on each metabolite:
  Input_data_NA <- replace(Input_data, Input_data==0, NA)#Shapiro test ignores NAs!

  shaptest <-   Input_data_NA %>% # Select data
    filter(Input_SettingsFile$Conditions %in% numerator | Input_SettingsFile$Conditions %in% denominator)%>%
    select_if(is.numeric)

  temp<- as.vector(sapply(shaptest, function(x) var(x)) == 0)#  we have to remove features with zero variance if there are any.
  shaptest <- shaptest[,!temp]

  shaptestres <- as.data.frame(sapply(shaptest, function(x) shapiro.test(x))) # do the test for each metabolite
  shaptestres <- as.data.frame(t(shaptestres))

  # 2. Give feedback to the user if the chosen test fits the data distribution. The data are normal if the p-value of the shapiro.test > 0.05.
  Norm <- format((round(sum(shaptestres$p.value > 0.05)/dim(shaptest)[2],4))*100, nsmall = 2) # Percentage of normally distributed metabolites across samples
  NotNorm <- format((round(sum(shaptestres$p.value < 0.05)/dim(shaptest)[2],4))*100, nsmall = 2) # Percentage of not-normally distributed metabolites across samples
  message(Norm, " % of the metabolites follow a normal distribution and ", NotNorm, " % of the metabolites are not-normally distributed according to the shapiro test. `shapiro.test` ignores missing values in the calculation.")

  if(ANOVA==TRUE){
    if (Norm > 50 & STAT_pval =="kruskal.test"){
      warning(Norm, " % of the metabolites follow a normal distribution but 'kruskal.test' for non parametric Hypothesis testing was chosen. Please consider selecting a parametric test (anova) instead.")
    }else if(NotNorm > 50 & STAT_pval =="anova"){
      message(NotNorm, " % of the metabolites follow a not-normal distribution but 'anova' for parametric Hypothesis testing was chosen. Please consider selecting a non-parametric test (kruskal.test) instead.")
    }
  }else{
    if (Norm > 50 & STAT_pval =="wilcox.test"){
      warning(Norm, " % of the metabolites follow a normal distribution but 'wilcox.test' for non parametric Hypothesis testing was chosen. Please consider selecting a parametric test (t.test) instead.")
    }else if(NotNorm > 50 & STAT_pval =="t.test"){
      message(NotNorm, " % of the metabolites follow a not-normal distribution but 't.test' for parametric Hypothesis testing was chosen. Please consider selecting a non-parametric test (wilcox.test) instead.")
    }else if((STAT_pval =="wilcox-test" & nrow(Num)<5)|(STAT_pval =="wilcox-test" & nrow(Denom)<5)){# check number of samples for wilcoxons test
      warning("Number of samples measured per condition is <5 in at least one of the two conditions, which is small for using wilcox.test. Consider using another test.")
    }
  }

  ## ------------ Create Results output folder ----------- ##
  name <- paste0("MetaProViz_Results_",Sys.Date())
  WorkD <- getwd()
  Results_folder <- file.path(WorkD, name) # Make Results folder
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)}
  Results_folder_DMA_folder <- file.path(Results_folder,"DMA") # Make DMA results folder
  if (!dir.exists(Results_folder_DMA_folder)) {dir.create(Results_folder_DMA_folder)}
  Results_folder_Conditions <- file.path(Results_folder_DMA_folder,paste0(toString(numerator),"_vs_",toString(denominator))) # Make comparison folder
  if (!dir.exists(Results_folder_Conditions)) {dir.create(Results_folder_Conditions)}

  ##################

  Log2FC_table <- data.frame(Metabolite = colnames(Input_data))
  for (column in 1:dim(comparisons)[2]){

    C1 <- Input_data %>% # Numerator
      filter(Input_SettingsFile$Conditions %in% comparisons[2,column]) %>%
      select_if(is.numeric)#only keep numeric columns with metabolite values
    C2 <- Input_data %>% # Deniminator
      filter(Input_SettingsFile$Conditions %in%  comparisons[1,column]) %>%
      select_if(is.numeric)


    ################################################################################################################################################################################################
    ############### Calculate Log2FC, pval, padj, tval ###############

    ## ------------  Calculate Log2FC ----------- ##
    # For C1_Mean and C2_Mean use 0 to obtain values, leading to Log2FC=NA if mean = 0 (If one value is NA, the mean will be NA even though all other values are available.)
    C1_Zero <- C1
    C1_Zero[is.na(C1_Zero)] <- 0
    Mean_C1 <- C1_Zero %>%
      summarise_all("mean")

    C2_Zero <- C2
    C2_Zero[is.na(C2_Zero)] <- 0
    Mean_C2 <- C2_Zero %>%
      summarise_all("mean")

    if(CoRe==TRUE){#Calculate Log2FC by taking into account the distance between the means:
      #CoRe values can be negative and positive, which can lead to problems when calculating the Log2FC and hence the Log2FC(A versus B)=(log2(A+x)-log2(B+x)) for A or B being negative, with x being a constant that is adapted to the size range of the respective metabolite.
      Mean_C1_t <- as.data.frame(t(Mean_C1))%>%
        rownames_to_column("Metabolite")
      Mean_C2_t <- as.data.frame(t(Mean_C2))%>%
        rownames_to_column("Metabolite")
      Mean_Merge <-merge(Mean_C1_t, Mean_C2_t, by="Metabolite", all=TRUE)%>%
        rename("C1"=2,
               "C2"=3)
      Mean_Merge$`NA/0` <- Mean_Merge$Metabolite %in% Metabolites_Miss#Column to enable the check if mean values of 0 are due to missing values (NA/0) and not by coincidence

      Mean_Merge <- Mean_Merge%>%#Now we can adapt the values to take into account the distance
        mutate(C1_Adapted = case_when(C1 < 0 & C2 > 0 ~ paste(C1*-1),#Here we have a negative value and transform it into a positive value
                                      C1 > 0 & C2 < 0 ~ paste(C1+(C2*-2)),#Here we have a positive value, but the other condition has a negative value that is transformed to a positive value, hence we add the distance
                                      C1 > 0 & C2 > 0 ~ paste(C1),#Both values are positive, so no action needed
                                      C1 < 0 & C2 < 0 ~ paste(C1),#Both values are negative, so no action needed
                                      C2 == 0 & `NA/0`== TRUE ~ paste(C1),#Here we have a "true" 0 value due to 0/NAs in the input data
                                      C1 == 0 & `NA/0`== TRUE ~ paste(C1),#Here we have a "true" 0 value due to 0/NAs in the input data
                                      C2 == 0 & `NA/0`== FALSE ~ paste(C1+1),#Here we have a "false" 0 value that occured at random and not due to 0/NAs in the input data, hence we add the constant +1
                                      C1 == 0 & `NA/0`== FALSE ~ paste(C1+1),#Here we have a "false" 0 value that occured at random and not due to 0/NAs in the input data, hence we add the constant +1
                                      TRUE ~ 'NA'))%>%
        mutate(C2_Adapted = case_when(C2 < 0 & C1 > 0 ~ paste(C2*-1),#Here we have a negative value and transform it into a positive value
                                      C2 > 0 & C1 < 0 ~ paste(C2+(C1*-2)),#Here we have a positive value, but the other condition has a negative value that is transformed to a positive value, hence we add the distance
                                      C2 > 0 & C1 > 0 ~ paste(C2),#Both values are positive, so no action needed
                                      C2 < 0 & C1 < 0 ~ paste(C2),#Both values are negative, so no action needed
                                      C1 == 0 & `NA/0`== TRUE ~ paste(C2),#Here we have a "true" 0 value due to 0/NAs in the input data
                                      C2 == 0 & `NA/0`== TRUE ~ paste(C2),#Here we have a "true" 0 value due to 0/NAs in the input data
                                      C1 == 0 & `NA/0`== FALSE ~ paste(C2+1),#Here we have a "false" 0 value that occured at random and not due to 0/NAs in the input data, hence we add the constant +1
                                      C2 == 0 & `NA/0`== FALSE ~ paste(C2+1),#Here we have a "false" 0 value that occured at random and not due to 0/NAs in the input data, hence we add the constant +1
                                      TRUE ~ 'NA'))%>%
        mutate(C1_Adapted = as.numeric(C1_Adapted))%>%
        mutate(C2_Adapted = as.numeric(C2_Adapted))

      if(any((Mean_Merge$`NA/0`==FALSE & Mean_Merge$C1 ==0) | (Mean_Merge$`NA/0`==FALSE & Mean_Merge$C2==0))==TRUE){
        X <- Mean_Merge%>%
          subset((Mean_Merge$`NA/0`==FALSE & Mean_Merge$C1 ==0) | (Mean_Merge$`NA/0`==FALSE & Mean_Merge$C2==0))
        message("We added +1 to the mean value of metabolite(s) ", paste0(X$Metabolite, collapse = ", "), ", since the mean of the replicate values where 0. This was not due to missing values (NA/0).")
      }

      #Calculate the Log2FC
      Mean_Merge$FC_C1vC2 <- Mean_Merge$C1_Adapted/Mean_Merge$C2_Adapted #FoldChange
      Mean_Merge$Log2FC <- gtools::foldchange2logratio(Mean_Merge$FC_C1vC2, base=2)
      Log2FC_C1vC2 <-Mean_Merge[,c(1,8)]

    }else if(CoRe==FALSE){
      #Mean values could be 0, which can not be used to calculate a Log2FC and hence the Log2FC(A versus B)=(log2(A+x)-log2(B+x)) for A and/or B being 0, with x being set to 1
      Mean_C1_t <- as.data.frame(t(Mean_C1))%>%
        rownames_to_column("Metabolite")
      Mean_C2_t <- as.data.frame(t(Mean_C2))%>%
        rownames_to_column("Metabolite")
      Mean_Merge <-merge(Mean_C1_t, Mean_C2_t, by="Metabolite", all=TRUE)%>%
        rename("C1"=2,
               "C2"=3)
      Mean_Merge$`NA/0` <- Mean_Merge$Metabolite %in% Metabolites_Miss#Column to enable the check if mean values of 0 are due to missing values (NA/0) and not by coincidence

      Mean_Merge <- Mean_Merge%>%#Now we can adapt the values to take into account the distance
        mutate(C1_Adapted = case_when(C2 == 0 & `NA/0`== TRUE ~ paste(C1),#Here we have a "true" 0 value due to 0/NAs in the input data
                                      C1 == 0 & `NA/0`== TRUE ~ paste(C1),#Here we have a "true" 0 value due to 0/NAs in the input data
                                      C2 == 0 & `NA/0`== FALSE ~ paste(C1+1),#Here we have a "false" 0 value that occured at random and not due to 0/NAs in the input data, hence we add the constant +1
                                      C1 == 0 & `NA/0`== FALSE ~ paste(C1+1),#Here we have a "false" 0 value that occured at random and not due to 0/NAs in the input data, hence we add the constant +1
                                      TRUE ~ paste(C1)))%>%
        mutate(C2_Adapted = case_when(C1 == 0 & `NA/0`== TRUE ~ paste(C2),#Here we have a "true" 0 value due to 0/NAs in the input data
                                      C2 == 0 & `NA/0`== TRUE ~ paste(C2),#Here we have a "true" 0 value due to 0/NAs in the input data
                                      C1 == 0 & `NA/0`== FALSE ~ paste(C2+1),#Here we have a "false" 0 value that occured at random and not due to 0/NAs in the input data, hence we add the constant +1
                                      C2 == 0 & `NA/0`== FALSE ~ paste(C2+1),#Here we have a "false" 0 value that occured at random and not due to 0/NAs in the input data, hence we add the constant +1
                                      TRUE ~ paste(C2)))%>%
        mutate(C1_Adapted = as.numeric(C1_Adapted))%>%
        mutate(C2_Adapted = as.numeric(C2_Adapted))

      if(any((Mean_Merge$`NA/0`==FALSE & Mean_Merge$C1 ==0) | (Mean_Merge$`NA/0`==FALSE & Mean_Merge$C2==0))==TRUE){
        X <- Mean_Merge%>%
          subset((Mean_Merge$`NA/0`==FALSE & Mean_Merge$C1 ==0) | (Mean_Merge$`NA/0`==FALSE & Mean_Merge$C2==0))
        message("We added +1 to the mean value of metabolite(s) ", paste0(X$Metabolite, collapse = ", "), ", since the mean of the replicate values where 0. This was not due to missing values (NA/0).")
      }

      #Calculate the Log2FC
      Mean_Merge$FC_C1vC2 <- Mean_Merge$C1_Adapted/Mean_Merge$C2_Adapted #FoldChange
      Mean_Merge$Log2FC <- gtools::foldchange2logratio(Mean_Merge$FC_C1vC2, base=2)
      Log2FC_C1vC2 <-Mean_Merge[,c(1,8)]

      if(ANOVA == TRUE){
        logname <- paste("Log2FC", paste( comparisons[2,column], comparisons[1,column],sep="-"), sep="_")
        names(Log2FC_C1vC2)[2] <- logname
      }
      Log2FC_table <- merge(Log2FC_table, Log2FC_C1vC2, by= "Metabolite")


    }else{
      stop("Please choose CoRe= TRUE or CoRe=FALSE.")
    }
  }


  ## ------------ Perform Hypothesis testing ----------- ##
  if(ANOVA == FALSE){
    # For C1 and C2 we use 0, since otherwise we can not perform the statistical testing.
    C1[is.na(C1)] <- 0
    C2[is.na(C2)] <- 0

    #### 1. p.value
    T_C1vC2 <-mapply(STAT_pval, x= as.data.frame(C2), y = as.data.frame(C1), SIMPLIFY = F)

    VecPVAL_C1vC2 <- c()
    for(i in 1:length(T_C1vC2)){
      p_value <- unlist(T_C1vC2[[i]][3])
      VecPVAL_C1vC2[i] <- p_value
    }
    Metabolite <- colnames(C2)
    PVal_C1vC2 <- data.frame(Metabolite, VecPVAL_C1vC2)

    #we set p.val= NA, for metabolites that had 1 or more replicates with NA/0 values and remove them prior to p-value adjustment
    PVal_C1vC2$`NA/0` <- PVal_C1vC2$Metabolite %in% Metabolites_Miss
    PVal_C1vC2 <-PVal_C1vC2%>%
      mutate(p.val = case_when(`NA/0`== TRUE ~ NA,
                               TRUE ~ paste(VecPVAL_C1vC2)))
    PVal_C1vC2$p.val = as.numeric(as.character(PVal_C1vC2$p.val))

    #### 2. p. adjusted
    #Split data for p.value adjustment to exclude NA
    PVal_NA <- PVal_C1vC2[is.na(PVal_C1vC2$p.val), c(1,4)]
    PVal_C1vC2 <-PVal_C1vC2[!is.na(PVal_C1vC2$p.val), c(1,4)]

    #perform adjustment
    VecPADJ_C1vC2 <- p.adjust((PVal_C1vC2[,2]),method = STAT_padj, n = length((PVal_C1vC2[,2]))) #p-adjusted
    Metabolite <- PVal_C1vC2[,1]
    PADJ_C1vC2 <- data.frame(Metabolite, VecPADJ_C1vC2)
    STAT_C1vC2 <- merge(PVal_C1vC2,PADJ_C1vC2, by="Metabolite")
    STAT_C1vC2 <- merge(Log2FC_table,STAT_C1vC2, by="Metabolite")
    names(STAT_C1vC2)[names(STAT_C1vC2) == "VecPADJ_C1vC2"] <- "p.adj"

    #### 3. t.value
    STAT_C1vC2$t.val <- qnorm((1 - STAT_C1vC2$p.val / 2)) * sign(STAT_C1vC2$Log2FC) # calculate and add t-value
    STAT_C1vC2 <- STAT_C1vC2[order(STAT_C1vC2$t.val,decreasing=TRUE),] # order the df based on the t-value

    #Add Metabolites that have p.val=NA back into the DF for completeness.
    if(nrow(PVal_NA)>0){
      PVal_NA <- merge(Log2FC_table,PVal_NA, by="Metabolite", all.y=TRUE)
      PVal_NA$p.adj <- NA
      PVal_NA$t.val <- NA
      STAT_C1vC2 <- rbind(STAT_C1vC2, PVal_NA)
    }

  }else{ # ANOVA = TRUE

    # for 1 vs all
    conditions =as.factor(conditions)
    #conditions=relevel(conditions, ref = "HK2")

    if(STAT_pval=="anova"){
      ## 1. Anova
      aov.res= apply(Input_data,2,function(x) aov(x~conditions))
      anova.res = lapply(aov.res  , anova)
      anova_res<-do.call('rbind', lapply(anova.res, function(x) {x["Pr(>F)"][1,]}))
      anova_res = cbind(anova_res,p.adjust(anova_res,STAT_padj )) %>% as.data.frame()
      colnames(anova_res) = c('Anova_p.val','Anova_p.adj')
      anova_res$Anova_p.val <- NULL # remove p.val

      ## 2. Tukey test
      posthoc.res = lapply(aov.res, TukeyHSD, conf.level=0.95)
      Tukey_res <- do.call('rbind', lapply(posthoc.res, function(x) x[1][[1]][,'p adj']))
      Tukey_res <- as.data.frame(Tukey_res)
      names(Tukey_res) <- paste("p.adj", names(Tukey_res), sep = "_")

      Pval_table <- merge(anova_res, Tukey_res, by=0)
      names(Pval_table)[1] <- "Metabolite"

      ### put Log2FC and pval dataframes
      merged_table <- merge(Pval_table, Log2FC_table, by = "Metabolite", all = TRUE)

    }else{# STAT_pval = kruskal.test
      # Kruskal test
      aov.res= apply(Input_data,2,function(x) kruskal.test(x~conditions))
      anova_res<-do.call('rbind', lapply(aov.res, function(x) {x["p.value"]}))
      anova_res <- as.matrix(mutate_all(as.data.frame(anova_res), function(x) as.numeric(as.character(x))))
      colnames(anova_res) = c("Kruskal_p.val")
      # Dunn test
      Dunndata <- Input_data %>%
        mutate(conditions = conditions) %>%
        select(conditions, everything())%>% as.data.frame()
      # apply doesnt for so going for a loop
      Dunn_res<- data.frame( comparisons = paste(comparisons[2,],    comparisons[1,], sep = "_" ))
      for(col in 2:dim(Dunndata)[2]){
        data = Dunndata[,c(1,col)]
        colnames(data)[2] <- gsub("^\\d+", "", colnames(data)[2])

        ## If a metabolite starts with number remove it
        formula <- as.formula(paste(colnames(data)[2], "~ conditions"))
        posthoc.res= rstatix::dunn_test(data, formula, p.adjust.method = STAT_padj)

        res <- data.frame(comparisons =  paste(posthoc.res$group2, posthoc.res$group1, sep = "_" ))
        res[[colnames(Dunndata)[col] ]] <-  posthoc.res$p.adj
        Dunn_res <- merge(Dunn_res,res,by="comparisons")
      }
      Dunn_res <- column_to_rownames(Dunn_res, "comparisons")%>% t() %>% as.data.frame()
      colnames(Dunn_res) <- paste("p.adj", colnames(Dunn_res), sep = "_")

      Dunn_res <- as.matrix(mutate_all(as.data.frame(Dunn_res), function(x) as.numeric(as.character(x))))

      Pval_table <- merge(anova_res, Dunn_res, by=0)
      names(Pval_table)[1] <- "Metabolite"

      merged_table <- merge(Pval_table, Log2FC_table, by = "Metabolite", all = TRUE)


    }
    #### 3. t.value
    for (i in 1:(dim(Pval_table)[2]-2)){
      t.val_name <- gsub("p.adj", "t.val", colnames(merged_table[i+2]) )
      merged_table[[t.val_name]] <- qnorm((1 - merged_table[[i+2]] / 2)) * sign(merged_table[[i+2+dim(Pval_table)[2]-2]]) # calculate and add t-value
    }

    # re-order columns
    order<- c(1,2)
    for(i in 1: (dim(Pval_table)[2]-2)){
      order <- append(order, dim(Pval_table)[2]+i)
      order <- append(order, 2+i)
      order <- append(order, (dim(Pval_table)[2] +dim(Log2FC_table)[2]-1) +i)
    }
    merged_table <- merged_table[, order]
    STAT_C1vC2 <- merged_table

  }


  ################################################################################################################################################################################################

  ## ------------ Add information on groups to DMA results----------- ##
  if(CoRe==TRUE){#add consumption release info to the result
    temp1 <- Mean_C1
    temp2 <- Mean_C2
    #Add Info:
    CoRe_info <- rbind(temp1, temp2,rep(0,length(temp1)))
    for (i in 1:length(temp1)){
      if (temp1[i]>0 & temp2[i]>0){
        CoRe_info[3,i] <- "Released"
      }else if (temp1[i]<0 & temp2[i]<0){
        CoRe_info[3,i] <- "Consumed"
      }else if(temp1[i]>0 & temp2[i]<0){
        CoRe_info[3,i] <- paste("Released in" ,numerator , "and Consumed",denominator , sep=" ")
      } else if(temp1[i]<0 & temp2[i]>0){
        CoRe_info[3,i] <- paste("Consumed in" ,numerator , " and Released",denominator , sep=" ")
      }else{
        CoRe_info[3,i] <- "No Change"
      }
    }

    CoRe_info <- t(CoRe_info) %>% as.data.frame()
    CoRe_info <- rownames_to_column(CoRe_info, "Metabolite")
    names(CoRe_info)[2] <- paste("Mean", "CoRe", numerator, sep="_")
    names(CoRe_info)[3] <- paste("Mean", "CoRe", denominator, sep="_")
    names(CoRe_info)[4] <- "CoRe"
    STAT_C1vC2 <- merge(STAT_C1vC2,CoRe_info,by= "Metabolite")
    STAT_C1vC2 <- STAT_C1vC2[order(STAT_C1vC2$t.val,decreasing=TRUE),] # order the df based on the t-value
  }

  ## ------------ Add pathway information to DMA results ----------- ##
  if(is.null(Input_Pathways)!=TRUE & 'Metabolite' %in% colnames(Input_Pathways) & 'Pathway' %in% colnames(Input_Pathways)){
    STAT_C1vC2$Metabolite %in% Input_Pathways$Metabolite
    STAT_C1vC2<- merge(STAT_C1vC2,Input_Pathways,by="Metabolite", all.x=T)
    STAT_C1vC2$Pathway[  is.na(STAT_C1vC2$Pathway)] <- "unknown"    # Merge the pathways to DMA result. All non matching metabolites get NA that are changed into "unknown".
  }else{
    warning("No column Pathway was added to the output. The pathway data must have 2 columns named: Metabolite and Pathway")
  }
  DMA_Output <- STAT_C1vC2

  if(ANOVA==TRUE){
    if(All==TRUE){
      numerator <- "Anova_All"
      denominator<- "All"
    }
  }
  if (Save_as_Results == "xlsx"){
    xlsDMA <- file.path(Results_folder_Conditions,paste0("DMA_Output_",toString(numerator),"_vs_",toString(denominator),"_", OutputName, ".xlsx"))   # Save the DMA results table
    writexl::write_xlsx(DMA_Output,xlsDMA, col_names = TRUE) # save the DMA result DF
  }else if (Save_as_Results == "csv"){
    csvDMA <- file.path(Results_folder_Conditions,paste0("DMA_Output_",toString(numerator),"_vs_",toString(denominator),"_", OutputName, ".csv"))
    write.csv(DMA_Output,csvDMA) # save the DMA result DF
  }else if (Save_as_Results == "txt"){
    txtDMA <- file.path(Results_folder_Conditions,paste0("DMA_Output_",toString(numerator),"_vs_",toString(denominator),"_", OutputName, ".txt"))
    write.table(DMA_Output,txtDMA, col.names = TRUE, row.names = FALSE) # save the DMA result DF
  }



  if(Plot == TRUE){ # Make a simple Volcano plot
    VolcanoPlot <- invisible(MetaProViz::VizVolcano(Plot_Settings="Standard",
                                                    Input_data=DMA_Output,
                                                    y= "p.adj",
                                                    x= "Log2FC",
                                                    AdditionalInput_data= NULL,
                                                    OutputPlotName= paste0(toString(numerator)," versus ",toString(denominator)),
                                                    Comparison_name= c(Input_data="Cond1", AdditionalInput_data= "Cond2"),
                                                    xlab= NULL,#"~Log[2]~FC"
                                                    ylab= NULL,#"~-Log[10]~p.adj"
                                                    pCutoff= 0.05,
                                                    FCcutoff= 0.5,
                                                    color_palette= NULL,
                                                    shape_palette=NULL,
                                                    SelectLab= DMA_Output$Metabolite,
                                                    Connectors=  FALSE,
                                                    Subtitle=  bquote(italic("Differential Metabolite Analysis")),
                                                    Theme= NULL,
                                                    Save_as_Plot= NULL))

    OutputPlotName = paste0(OutputName,"_padj_",0.05,"Log2FC_",0.5)

    volcanoDMA <- file.path(Results_folder_Conditions,paste0( "Volcano_Plot_",toString(numerator),"-versus-",toString(denominator),"_", OutputPlotName,".",Save_as_Plot))
    ggsave(volcanoDMA,plot=VolcanoPlot, width=10, height=8) # save the voplcano plot

    plot(VolcanoPlot)
  }
  assign(paste0("DMA_",toString(numerator),"_vs_",toString(denominator)), DMA_Output, envir=.GlobalEnv)
}




