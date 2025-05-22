#' ## ---------------------------
##
## Script name: MCA
##
## Purpose of script: Metabolite Clustering Analysis generates clusters of metabolites based on regulatory rules.
##
## Author: Christina Schmidt
##
## Date Created: 2022-10-28
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

####################################################################
### ### ### Metabolite Clustering Analysis: 2 Conditions ### ### ###
####################################################################

#' This script performs metabolite clustering analysis and computes clusters of metabolites based on regulatory rules between conditions.
#'
#' @param data_c1 DF for your data (results from e.g. dma) containing metabolites in rows with corresponding Log2FC and stat (p-value, p.adjusted) value columns.
#' @param data_c2 DF for your data (results from e.g. dma) containing metabolites in rows with corresponding Log2FC and stat (p-value, p.adjusted) value columns.
#' @param metadata_info_c1  \emph{Optional: } Pass ColumnNames and Cutoffs for condition 1 including the value column (e.g. Log2FC, Log2Diff, t.val, etc) and the stats column (e.g. p.adj, p.val). This must include: c(ValueCol=ColumnName_data_c1,StatCol=ColumnName_data_c1, cutoff_stat= NumericValue, ValueCutoff=NumericValue) \strong{Default=c(ValueCol="Log2FC",StatCol="p.adj", cutoff_stat= 0.05, ValueCutoff=1)}
#' @param metadata_info_c2  \emph{Optional: } Pass ColumnNames and Cutoffs for condition 2 includingthe value column (e.g. Log2FC, Log2Diff, t.val, etc) and the stats column (e.g. p.adj, p.val). This must include: c(ValueCol=ColumnName_data_c2,StatCol=ColumnName_data_c2, cutoff_stat= NumericValue, ValueCutoff=NumericValue)\strong{Default=c(ValueCol="Log2FC",StatCol="p.adj", cutoff_stat= 0.05, ValueCutoff=1)}
#' @param feature \emph{Optional: } Column name of Column including the Metabolite identifiers. This MUST BE THE SAME in each of your Input files. \strong{Default="Metabolite"}
#' @param method_background \emph{Optional: } Background method C1|C2, C1&C2, C2, C1 or * \strong{Default="C1&C2"}
#' @param save_table \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt" \strong{Default = "csv"}
#' @param path \emph{Optional:} Path to the folder the results should be saved at. \strong{Default = NULL}
#'
#' @return List of two DFs: 1. summary of the cluster count and 2. the detailed information of each metabolites in the clusters.
#'
#' @examples
#' Intra <- intracell_raw%>%tibble::column_to_rownames("Code")
#' Input <- MetaProViz::dma(data=Intra[-c(49:58) ,-c(1:3)], metadata_sample=Intra[-c(49:58) , c(1:3)], metadata_info = c(Conditions = "Conditions", Numerator = NULL, Denominator  = "HK2"))
#'
#' Res <- MetaProViz::mca_2cond(data_c1 = Input[["dma"]][["786-O_vs_HK2"]],
#'                              data_c2 = Input[["dma"]][["786-M1A_vs_HK2"]])
#'
#' @keywords biological clustering
#'
#' @importFrom dplyr rename mutate case_when mutate_at count
#' @importFrom tidyr replace_na
#' @importFrom magrittr %>% %<>%
#' @importFrom tibble rownames_to_column column_to_rownames

#' @export
#'
mca_2cond <- function(data_c1,
                      data_c2,
                      metadata_info_c1=c(ValueCol="Log2FC",StatCol="p.adj", cutoff_stat= 0.05, ValueCutoff=1),
                      metadata_info_c2=c(ValueCol="Log2FC",StatCol="p.adj", cutoff_stat= 0.05, ValueCutoff=1),
                      feature = "Metabolite",
                      save_table = "csv",
                      method_background="C1&C2",
                      path=NULL
                      ){

  ################################################################################################################################################################################################
  ## ------------ Check Input files ----------- ##
  check_param_mca(data_c1=data_c1,
                              data_c2=data_c2,
                              data_core=NULL,
                              data_intra=NULL,
                              metadata_info_c1=metadata_info_c1,
                              metadata_info_c2=metadata_info_c2,
                              metadata_info_core=NULL,
                              metadata_info_intra=NULL,
                              method_background=method_background,
                              feature=feature,
                              save_table=save_table)

  ## ------------ Create Results output folder ----------- ##
  if(is.null(save_table)==FALSE){
    folder <- save_path(folder_name= "MCA2Cond",
                                    path=path)
  }

  ################################################################################################################################################################################################
  ## ------------ Prepare the Input -------- ##
  #Import the data and check columns (here the user will get an error if the column can not be renamed as it does not exists.)
  Cond1_DF <- as.data.frame(data_c1)%>%
    dplyr::rename("MetaboliteID"=paste(feature),
                  "ValueCol"=metadata_info_c1[["ValueCol"]],
                  "PadjCol"=metadata_info_c1[["StatCol"]])
  Cond1_DF <- Cond1_DF[complete.cases(Cond1_DF$ValueCol, Cond1_DF$PadjCol), ]

  Cond2_DF<- as.data.frame(data_c2)%>%
    dplyr::rename("MetaboliteID"=paste(feature),
                  "ValueCol"=metadata_info_c2[["ValueCol"]],
                  "PadjCol"=metadata_info_c2[["StatCol"]])
  Cond2_DF <- Cond2_DF[complete.cases(Cond2_DF$ValueCol, Cond2_DF$PadjCol), ]

  #Tag genes that are detected in each data layer
  Cond1_DF$Detected <- "TRUE"
  Cond2_DF$Detected <- "TRUE"

  ## ------------ Assign Groups -------- ##
  #Assign to Group based on individual Cutoff ("UP", "DOWN", "No Change")
  Cond1_DF <-Cond1_DF%>%
    dplyr::mutate(Cutoff = dplyr::case_when(Cond1_DF$PadjCol <= as.numeric(metadata_info_c1[["cutoff_stat"]]) & Cond1_DF$ValueCol > as.numeric(metadata_info_c1[["ValueCutoff"]]) ~ 'UP',
                              Cond1_DF$PadjCol <= as.numeric(metadata_info_c1[["cutoff_stat"]]) & Cond1_DF$ValueCol < - as.numeric(metadata_info_c1[["ValueCutoff"]]) ~ 'DOWN',
                              TRUE ~ 'No Change'))%>%
    dplyr::mutate(Cutoff_Specific = dplyr::case_when(Cutoff == "UP" ~ 'UP',
                                       Cutoff == "DOWN" ~ 'DOWN',
                                       Cutoff == "No Change" & Cond1_DF$PadjCol <= as.numeric(metadata_info_c1[["cutoff_stat"]]) & Cond1_DF$ValueCol > 0 ~ 'Significant Positive',
                                       Cutoff == "No Change" & Cond1_DF$PadjCol <= as.numeric(metadata_info_c1[["cutoff_stat"]]) & Cond1_DF$ValueCol < 0 ~ 'Significant Negative',
                                       Cutoff == "No Change" & Cond1_DF$PadjCol > as.numeric(metadata_info_c1[["cutoff_stat"]]) ~ 'Not Significant',
                                       TRUE ~ 'NA'))

  Cond2_DF <- Cond2_DF%>%
    dplyr::mutate(Cutoff = dplyr::case_when(Cond2_DF$PadjCol <= as.numeric(metadata_info_c2[["cutoff_stat"]]) & Cond2_DF$ValueCol > as.numeric(metadata_info_c2[["ValueCutoff"]]) ~ 'UP',
                              Cond2_DF$PadjCol <= as.numeric(metadata_info_c2[["cutoff_stat"]]) & Cond2_DF$ValueCol < - as.numeric(metadata_info_c2[["ValueCutoff"]]) ~ 'DOWN',
                              TRUE ~ 'No Change')) %>%
    dplyr::mutate(Cutoff_Specific = dplyr::case_when(Cutoff == "UP" ~ 'UP',
                                       Cutoff == "DOWN" ~ 'DOWN',
                                       Cutoff == "No Change" & Cond2_DF$PadjCol <= as.numeric(metadata_info_c2[["cutoff_stat"]]) & Cond2_DF$ValueCol > 0 ~ 'Significant Positive',
                                       Cutoff == "No Change" & Cond2_DF$PadjCol <= as.numeric(metadata_info_c2[["cutoff_stat"]]) & Cond2_DF$ValueCol < 0 ~ 'Significant Negative',
                                       Cutoff == "No Change" & Cond2_DF$PadjCol > as.numeric(metadata_info_c2[["cutoff_stat"]]) ~ 'Not Significant',
                                       TRUE ~ 'NA'))



  #Merge the dataframes together: Merge the supplied Cond1 and Cond2 dataframes together.
  ##Add prefix to column names to distinguish the different data types after merge
  colnames(Cond2_DF) <- paste0("Cond2_DF_", colnames(Cond2_DF))
  Cond2_DF <- Cond2_DF%>%
    dplyr::rename("MetaboliteID" = "Cond2_DF_MetaboliteID")

  colnames(Cond1_DF) <- paste0("Cond1_DF_", colnames(Cond1_DF))
  Cond1_DF <-Cond1_DF%>%
    dplyr::rename("MetaboliteID"="Cond1_DF_MetaboliteID")

  ##Merge
  MergeDF <- merge(Cond1_DF, Cond2_DF, by="MetaboliteID", all=TRUE)

  ##Mark the undetected genes in each data layer
  MergeDF<-MergeDF %>%
    dplyr::mutate_at(c("Cond2_DF_Detected","Cond1_DF_Detected"), ~tidyr::replace_na(.,"FALSE"))%>%
    dplyr::mutate_at(c("Cond2_DF_Cutoff","Cond1_DF_Cutoff"), ~tidyr::replace_na(.,"No Change"))%>%
    dplyr::mutate_at(c("Cond2_DF_Cutoff_Specific", "Cond1_DF_Cutoff_Specific"), ~tidyr::replace_na(.,"Not Detected"))

  #Apply Background filter (label genes that will be removed based on choosen background)
  if(method_background == "C1|C2"){# C1|C2 = Cond2 OR Cond1
    MergeDF <- MergeDF%>%
      dplyr::mutate(BG_method = dplyr::case_when(Cond1_DF_Detected=="TRUE" & Cond2_DF_Detected=="TRUE" ~ 'TRUE', #Cond1 & Cond2
                                   Cond1_DF_Detected=="TRUE" & Cond2_DF_Detected=="FALSE" ~ 'TRUE', # JustCond1
                                   Cond1_DF_Detected=="FALSE" & Cond2_DF_Detected=="TRUE" ~ 'TRUE', # Just Cond2
                                   TRUE ~ 'FALSE'))
  }else if(method_background == "C1&C2"){ # Cond2 AND Cond1
    MergeDF <- MergeDF%>%
      dplyr::mutate(BG_method = dplyr::case_when(Cond1_DF_Detected=="TRUE" & Cond2_DF_Detected=="TRUE" ~ 'TRUE', #Cond1 & Cond2
                                   TRUE ~ 'FALSE'))
  }else if(method_background == "C2"){ # Cond2 has to be there
    MergeDF <- MergeDF%>%
      dplyr::mutate(BG_method = dplyr::case_when(Cond1_DF_Detected=="TRUE" & Cond2_DF_Detected=="TRUE" ~ 'TRUE', #Cond1 & Cond2
                                   Cond1_DF_Detected=="FALSE" & Cond2_DF_Detected=="TRUE" ~ 'TRUE', # Just Cond2
                                   TRUE ~ 'FALSE'))
  }else if(method_background == "C1"){ #Cond1 has to be there
    MergeDF <- MergeDF%>%
      dplyr::mutate(BG_method = dplyr::case_when(Cond1_DF_Detected=="TRUE" & Cond2_DF_Detected=="TRUE" ~ 'TRUE', #Cond1 & Cond2
                                   Cond1_DF_Detected=="TRUE" & Cond2_DF_Detected=="FALSE" ~ 'TRUE', # JustCond1
                                   TRUE ~ 'FALSE'))
  }else if(method_background == "*"){ # Use all genes as the background
    MergeDF$BG_method <- "TRUE"
  }else{
    stop("Please use one of the following method_backgrounds: C1|C2, C1&C2, C2, C1, *")#error message
  }

  #Assign SiRCle cluster names to the genes
  MergeDF <- MergeDF%>%
    dplyr::mutate(RG1_Specific_Cond2 = dplyr::case_when(BG_method =="FALSE"~ 'Background = FALSE',
                                          Cond1_DF_Cutoff=="DOWN" & Cond2_DF_Cutoff_Specific=="DOWN" ~ 'Cond1 DOWN + Cond2 DOWN',#State 1
                                          Cond1_DF_Cutoff=="DOWN" & Cond2_DF_Cutoff_Specific=="Not Detected" ~ 'Cond1 DOWN + Cond2 Not Detected',#State 2
                                          Cond1_DF_Cutoff=="DOWN" & Cond2_DF_Cutoff_Specific=="Not Significant" ~ 'Cond1 DOWN + Cond2 Not Significant',#State 3
                                          Cond1_DF_Cutoff=="DOWN" & Cond2_DF_Cutoff_Specific=="Significant Negative" ~ 'Cond1 DOWN + Cond2 Significant Negative',#State 4
                                          Cond1_DF_Cutoff=="DOWN" & Cond2_DF_Cutoff_Specific=="Significant Positive" ~ 'Cond1 DOWN + Cond2 Significant Positive',#State 5
                                          Cond1_DF_Cutoff=="DOWN" & Cond2_DF_Cutoff_Specific=="UP" ~ 'Cond1 DOWN + Cond2 UP',#State 6

                                          Cond1_DF_Cutoff=="No Change" & Cond2_DF_Cutoff_Specific=="DOWN" ~ 'Cond1 No Change + Cond2 DOWN',#State 7
                                          Cond1_DF_Cutoff=="No Change" & Cond2_DF_Cutoff_Specific=="Not Detected" ~ 'Cond1 No Change + Cond2 Not Detected',#State 8
                                          Cond1_DF_Cutoff=="No Change" & Cond2_DF_Cutoff_Specific=="Not Significant" ~ 'Cond1 No Change + Cond2 Not Significant',#State 9
                                          Cond1_DF_Cutoff=="No Change" & Cond2_DF_Cutoff_Specific=="Significant Negative" ~ 'Cond1 No Change + Cond2 Significant Negative',#State 10
                                          Cond1_DF_Cutoff=="No Change" & Cond2_DF_Cutoff_Specific=="Significant Positive" ~ 'Cond1 No Change + Cond2 Significant Positive',#State 11
                                          Cond1_DF_Cutoff=="No Change" & Cond2_DF_Cutoff_Specific=="UP" ~ 'Cond1 No Change + Cond2 UP',#State 6

                                          Cond1_DF_Cutoff=="UP" & Cond2_DF_Cutoff_Specific=="DOWN" ~ 'Cond1 UP + Cond2 DOWN',#State 12
                                          Cond1_DF_Cutoff=="UP" & Cond2_DF_Cutoff_Specific=="Not Detected" ~ 'Cond1 UP + Cond2 Not Detected',#State 13
                                          Cond1_DF_Cutoff=="UP" & Cond2_DF_Cutoff_Specific=="Not Significant" ~ 'Cond1 UP + Cond2 Not Significant',#State 14
                                          Cond1_DF_Cutoff=="UP" & Cond2_DF_Cutoff_Specific=="Significant Negative" ~ 'Cond1 UP + Cond2 Significant Negative',#State 15
                                          Cond1_DF_Cutoff=="UP" & Cond2_DF_Cutoff_Specific=="Significant Positive" ~ 'Cond1 UP + Cond2 Significant Positive',#State 16
                                          Cond1_DF_Cutoff=="UP" & Cond2_DF_Cutoff_Specific=="UP" ~ 'Cond1 UP + Cond2 UP',#State 17
                                          TRUE ~ 'NA'))%>%
    dplyr::mutate(RG1_Specific_Cond1 = dplyr::case_when(BG_method =="FALSE"~ 'Background = FALSE',
                                          Cond2_DF_Cutoff=="DOWN" & Cond1_DF_Cutoff_Specific=="DOWN" ~ 'Cond2 DOWN + Cond1 DOWN',#State 1
                                          Cond2_DF_Cutoff=="DOWN" & Cond1_DF_Cutoff_Specific=="Not Detected" ~ 'Cond2 DOWN + Cond1 Not Detected',#State 2
                                          Cond2_DF_Cutoff=="DOWN" & Cond1_DF_Cutoff_Specific=="Not Significant" ~ 'Cond2 DOWN + Cond1 Not Significant',#State 3
                                          Cond2_DF_Cutoff=="DOWN" & Cond1_DF_Cutoff_Specific=="Significant Negative" ~ 'Cond2 DOWN + Cond1 Significant Negative',#State 4
                                          Cond2_DF_Cutoff=="DOWN" & Cond1_DF_Cutoff_Specific=="Significant Positive" ~ 'Cond2 DOWN + Cond1 Significant Positive',#State 5
                                          Cond2_DF_Cutoff=="DOWN" & Cond1_DF_Cutoff_Specific=="UP" ~ 'Cond2 DOWN + Cond1 UP',#State 6

                                          Cond2_DF_Cutoff=="No Change" & Cond1_DF_Cutoff_Specific=="DOWN" ~ 'Cond2 No Change + Cond1 DOWN',#State 7
                                          Cond2_DF_Cutoff=="No Change" & Cond1_DF_Cutoff_Specific=="Not Detected" ~ 'Cond2 No Change + Cond1 Not Detected',#State 8
                                          Cond2_DF_Cutoff=="No Change" & Cond1_DF_Cutoff_Specific=="Not Significant" ~ 'Cond2 No Change + Cond1 Not Significant',#State 9
                                          Cond2_DF_Cutoff=="No Change" & Cond1_DF_Cutoff_Specific=="Significant Negative" ~ 'Cond2 No Change + Cond1 Significant Negative',#State 10
                                          Cond2_DF_Cutoff=="No Change" & Cond1_DF_Cutoff_Specific=="Significant Positive" ~ 'Cond2 No Change + Cond1 Significant Positive',#State 11
                                          Cond2_DF_Cutoff=="No Change" & Cond1_DF_Cutoff_Specific=="UP" ~ 'Cond2 No Change + Cond1 UP',#State 6

                                          Cond2_DF_Cutoff=="UP" & Cond1_DF_Cutoff_Specific=="DOWN" ~ 'Cond2 UP + Cond1 DOWN',#State 12
                                          Cond2_DF_Cutoff=="UP" & Cond1_DF_Cutoff_Specific=="Not Detected" ~ 'Cond2 UP + Cond1 Not Detected',#State 13
                                          Cond2_DF_Cutoff=="UP" & Cond1_DF_Cutoff_Specific=="Not Significant" ~ 'Cond2 UP + Cond1 Not Significant',#State 14
                                          Cond2_DF_Cutoff=="UP" & Cond1_DF_Cutoff_Specific=="Significant Negative" ~ 'Cond2 UP + Cond1 Significant Negative',#State 15
                                          Cond2_DF_Cutoff=="UP" & Cond1_DF_Cutoff_Specific=="Significant Positive" ~ 'Cond2 UP + Cond1 Significant Positive',#State 16
                                          Cond2_DF_Cutoff=="UP" & Cond1_DF_Cutoff_Specific=="UP" ~ 'Cond2 UP + Cond1 UP',#State 17
                                          TRUE ~ 'NA'))%>%
    dplyr::mutate(RG1_All = dplyr::case_when(BG_method =="FALSE"~ 'Background = FALSE',
                               Cond1_DF_Cutoff_Specific=="DOWN" & Cond2_DF_Cutoff_Specific=="DOWN" ~ 'Cond1 DOWN + Cond2 DOWN',#State 1
                               Cond1_DF_Cutoff_Specific=="DOWN" & Cond2_DF_Cutoff_Specific=="Not Detected" ~ 'Cond1 DOWN + Cond2 Not Detected',#State 2
                               Cond1_DF_Cutoff_Specific=="DOWN" & Cond2_DF_Cutoff_Specific=="Not Significant" ~ 'Cond1 DOWN + Cond2 Not Significant',#State 3
                               Cond1_DF_Cutoff_Specific=="DOWN" & Cond2_DF_Cutoff_Specific=="Significant Negative" ~ 'Cond1 DOWN + Cond2 Significant Negative',#State 4
                               Cond1_DF_Cutoff_Specific=="DOWN" & Cond2_DF_Cutoff_Specific=="Significant Positive" ~ 'Cond1 DOWN + Cond2 Significant Positive',#State 5
                               Cond1_DF_Cutoff_Specific=="DOWN" & Cond2_DF_Cutoff_Specific=="UP" ~ 'Cond1 DOWN + Cond2 UP',#State 6

                               Cond1_DF_Cutoff_Specific=="UP" & Cond2_DF_Cutoff_Specific=="DOWN" ~ 'Cond1 UP + Cond2 DOWN',#State 12
                               Cond1_DF_Cutoff_Specific=="UP" & Cond2_DF_Cutoff_Specific=="Not Detected" ~ 'Cond1 UP + Cond2 Not Detected',#State 13
                               Cond1_DF_Cutoff_Specific=="UP" & Cond2_DF_Cutoff_Specific=="Not Significant" ~ 'Cond1 UP + Cond2 Not Significant',#State 14
                               Cond1_DF_Cutoff_Specific=="UP" & Cond2_DF_Cutoff_Specific=="Significant Negative" ~ 'Cond1 UP + Cond2 Significant Negative',#State 15
                               Cond1_DF_Cutoff_Specific=="UP" & Cond2_DF_Cutoff_Specific=="Significant Positive" ~ 'Cond1 UP + Cond2 Significant Positive',#State 16
                               Cond1_DF_Cutoff_Specific=="UP" & Cond2_DF_Cutoff_Specific=="UP" ~ 'Cond1 UP + Cond2 UP',#State 17

                               Cond1_DF_Cutoff_Specific=="Not Detected" & Cond2_DF_Cutoff_Specific=="DOWN" ~ 'Cond1 Not Detected + Cond2 DOWN',#State 12
                               Cond1_DF_Cutoff_Specific=="Not Detected" & Cond2_DF_Cutoff_Specific=="Not Detected" ~ 'Cond1 Not Detected + Cond2 Not Detected',#State 13
                               Cond1_DF_Cutoff_Specific=="Not Detected" & Cond2_DF_Cutoff_Specific=="Not Significant" ~ 'Cond1 Not Detected + Cond2 Not Significant',#State 14
                               Cond1_DF_Cutoff_Specific=="Not Detected" & Cond2_DF_Cutoff_Specific=="Significant Negative" ~ 'Cond1 Not Detected + Cond2 Significant Negative',#State 15
                               Cond1_DF_Cutoff_Specific=="Not Detected" & Cond2_DF_Cutoff_Specific=="Significant Positive" ~ 'Cond1 Not Detected + Cond2 Significant Positive',#State 16
                               Cond1_DF_Cutoff_Specific=="Not Detected" & Cond2_DF_Cutoff_Specific=="UP" ~ 'Cond1 Not Detected + Cond2 UP',#State 17

                               Cond1_DF_Cutoff_Specific=="Significant Negative" & Cond2_DF_Cutoff_Specific=="DOWN" ~ 'Cond1 Significant Negative + Cond2 DOWN',#State 12
                               Cond1_DF_Cutoff_Specific=="Significant Negative" & Cond2_DF_Cutoff_Specific=="Not Detected" ~ 'Cond1 Significant Negative + Cond2 Not Detected',#State 13
                               Cond1_DF_Cutoff_Specific=="Significant Negative" & Cond2_DF_Cutoff_Specific=="Not Significant" ~ 'Cond1 Significant Negative + Cond2 Not Significant',#State 14
                               Cond1_DF_Cutoff_Specific=="Significant Negative" & Cond2_DF_Cutoff_Specific=="Significant Negative" ~ 'Cond1 Significant Negative + Cond2 Significant Negative',#State 15
                               Cond1_DF_Cutoff_Specific=="Significant Negative" & Cond2_DF_Cutoff_Specific=="Significant Positive" ~ 'Cond1 Significant Negative + Cond2 Significant Positive',#State 16
                               Cond1_DF_Cutoff_Specific=="Significant Negative" & Cond2_DF_Cutoff_Specific=="UP" ~ 'Cond1 Significant Negative + Cond2 UP',#State 17

                               Cond1_DF_Cutoff_Specific=="Significant Positive" & Cond2_DF_Cutoff_Specific=="DOWN" ~ 'Cond1 Significant Positive + Cond2 DOWN',#State 12
                               Cond1_DF_Cutoff_Specific=="Significant Positive" & Cond2_DF_Cutoff_Specific=="Not Detected" ~ 'Cond1 Significant Positive + Cond2 Not Detected',#State 13
                               Cond1_DF_Cutoff_Specific=="Significant Positive" & Cond2_DF_Cutoff_Specific=="Not Significant" ~ 'Cond1 Significant Positive + Cond2 Not Significant',#State 14
                               Cond1_DF_Cutoff_Specific=="Significant Positive" & Cond2_DF_Cutoff_Specific=="Significant Negative" ~ 'Cond1 Significant Positive + Cond2 Significant Negative',#State 15
                               Cond1_DF_Cutoff_Specific=="Significant Positive" & Cond2_DF_Cutoff_Specific=="Significant Positive" ~ 'Cond1 Significant Positive + Cond2 Significant Positive',#State 16
                               Cond1_DF_Cutoff_Specific=="Significant Positive" & Cond2_DF_Cutoff_Specific=="UP" ~ 'Cond1 Significant Positive + Cond2 UP',#State 17

                               Cond1_DF_Cutoff_Specific=="Not Significant" & Cond2_DF_Cutoff_Specific=="DOWN" ~ 'Cond1 Not Significant + Cond2 DOWN',#State 12
                               Cond1_DF_Cutoff_Specific=="Not Significant" & Cond2_DF_Cutoff_Specific=="Not Detected" ~ 'Cond1 Not Significant + Cond2 Not Detected',#State 13
                               Cond1_DF_Cutoff_Specific=="Not Significant" & Cond2_DF_Cutoff_Specific=="Not Significant" ~ 'Cond1 Not Significant + Cond2 Not Significant',#State 14
                               Cond1_DF_Cutoff_Specific=="Not Significant" & Cond2_DF_Cutoff_Specific=="Significant Negative" ~ 'Cond1 Not Significant + Cond2 Significant Negative',#State 15
                               Cond1_DF_Cutoff_Specific=="Not Significant" & Cond2_DF_Cutoff_Specific=="Significant Positive" ~ 'Cond1 Not Significant + Cond2 Significant Positive',#State 16
                               Cond1_DF_Cutoff_Specific=="Not Significant" & Cond2_DF_Cutoff_Specific=="UP" ~ 'Cond1 Not Significant + Cond2 UP',#State 1
                               TRUE ~ 'NA'))%>%
    dplyr::mutate(RG2_Significant = dplyr::case_when(BG_method =="FALSE"~ 'Background = FALSE',
                                       Cond1_DF_Cutoff_Specific=="DOWN" & Cond2_DF_Cutoff_Specific=="DOWN" ~ 'core_DOWN',#State 1
                                       Cond1_DF_Cutoff_Specific=="DOWN" & Cond2_DF_Cutoff_Specific=="Not Detected" ~ 'Cond1_DOWN',#State 2
                                       Cond1_DF_Cutoff_Specific=="DOWN" & Cond2_DF_Cutoff_Specific=="Not Significant" ~ 'Cond1_DOWN',#State 3
                                       Cond1_DF_Cutoff_Specific=="DOWN" & Cond2_DF_Cutoff_Specific=="Significant Negative" ~ 'core_DOWN',#State 4
                                       Cond1_DF_Cutoff_Specific=="DOWN" & Cond2_DF_Cutoff_Specific=="Significant Positive" ~ 'Opposite',#State 5
                                       Cond1_DF_Cutoff_Specific=="DOWN" & Cond2_DF_Cutoff_Specific=="UP" ~ 'Opposite',#State 6

                                       Cond1_DF_Cutoff_Specific=="UP" & Cond2_DF_Cutoff_Specific=="DOWN" ~ 'Opposite',#State 12
                                       Cond1_DF_Cutoff_Specific=="UP" & Cond2_DF_Cutoff_Specific=="Not Detected" ~ 'Cond1_UP',#State 13
                                       Cond1_DF_Cutoff_Specific=="UP" & Cond2_DF_Cutoff_Specific=="Not Significant" ~ 'Cond1_UP',#State 14
                                       Cond1_DF_Cutoff_Specific=="UP" & Cond2_DF_Cutoff_Specific=="Significant Negative" ~ 'Opposite',#State 15
                                       Cond1_DF_Cutoff_Specific=="UP" & Cond2_DF_Cutoff_Specific=="Significant Positive" ~ 'core_UP',#State 16
                                       Cond1_DF_Cutoff_Specific=="UP" & Cond2_DF_Cutoff_Specific=="UP" ~ 'core_UP',#State 17

                                       Cond1_DF_Cutoff_Specific=="Not Detected" & Cond2_DF_Cutoff_Specific=="DOWN" ~ 'Cond2_DOWN',#State 12
                                       Cond1_DF_Cutoff_Specific=="Not Detected" & Cond2_DF_Cutoff_Specific=="Not Detected" ~ 'None',#State 13
                                       Cond1_DF_Cutoff_Specific=="Not Detected" & Cond2_DF_Cutoff_Specific=="Not Significant" ~ 'None',#State 14
                                       Cond1_DF_Cutoff_Specific=="Not Detected" & Cond2_DF_Cutoff_Specific=="Significant Negative" ~ 'None',#State 15
                                       Cond1_DF_Cutoff_Specific=="Not Detected" & Cond2_DF_Cutoff_Specific=="Significant Positive" ~ 'None',#State 16
                                       Cond1_DF_Cutoff_Specific=="Not Detected" & Cond2_DF_Cutoff_Specific=="UP" ~ 'Cond2_UP',#State 17

                                       Cond1_DF_Cutoff_Specific=="Significant Negative" & Cond2_DF_Cutoff_Specific=="DOWN" ~ 'core_DOWN',#State 12
                                       Cond1_DF_Cutoff_Specific=="Significant Negative" & Cond2_DF_Cutoff_Specific=="Not Detected" ~ 'None',#State 13
                                       Cond1_DF_Cutoff_Specific=="Significant Negative" & Cond2_DF_Cutoff_Specific=="Not Significant" ~ 'None',#State 14
                                       Cond1_DF_Cutoff_Specific=="Significant Negative" & Cond2_DF_Cutoff_Specific=="Significant Negative" ~ 'None',#State 15
                                       Cond1_DF_Cutoff_Specific=="Significant Negative" & Cond2_DF_Cutoff_Specific=="Significant Positive" ~ 'None',#State 16
                                       Cond1_DF_Cutoff_Specific=="Significant Negative" & Cond2_DF_Cutoff_Specific=="UP" ~ 'Opposite',#State 17

                                       Cond1_DF_Cutoff_Specific=="Significant Positive" & Cond2_DF_Cutoff_Specific=="DOWN" ~ 'Opposite',#State 12
                                       Cond1_DF_Cutoff_Specific=="Significant Positive" & Cond2_DF_Cutoff_Specific=="Not Detected" ~ 'None',#State 13
                                       Cond1_DF_Cutoff_Specific=="Significant Positive" & Cond2_DF_Cutoff_Specific=="Not Significant" ~ 'None',#State 14
                                       Cond1_DF_Cutoff_Specific=="Significant Positive" & Cond2_DF_Cutoff_Specific=="Significant Negative" ~ 'None',#State 15
                                       Cond1_DF_Cutoff_Specific=="Significant Positive" & Cond2_DF_Cutoff_Specific=="Significant Positive" ~ 'None',#State 16
                                       Cond1_DF_Cutoff_Specific=="Significant Positive" & Cond2_DF_Cutoff_Specific=="UP" ~ 'core_UP',#State 17

                                       Cond1_DF_Cutoff_Specific=="Not Significant" & Cond2_DF_Cutoff_Specific=="DOWN" ~ 'Cond2_DOWN',#State 12
                                       Cond1_DF_Cutoff_Specific=="Not Significant" & Cond2_DF_Cutoff_Specific=="Not Detected" ~ 'None',#State 13
                                       Cond1_DF_Cutoff_Specific=="Not Significant" & Cond2_DF_Cutoff_Specific=="Not Significant" ~ 'None',#State 14
                                       Cond1_DF_Cutoff_Specific=="Not Significant" & Cond2_DF_Cutoff_Specific=="Significant Negative" ~ 'None',#State 15
                                       Cond1_DF_Cutoff_Specific=="Not Significant" & Cond2_DF_Cutoff_Specific=="Significant Positive" ~ 'None',#State 16
                                       Cond1_DF_Cutoff_Specific=="Not Significant" & Cond2_DF_Cutoff_Specific=="UP" ~ 'Cond2_UP',#State 1
                                       TRUE ~ 'NA'))%>%
    dplyr::mutate(RG3_SignificantChange = dplyr::case_when(BG_method =="FALSE"~ 'Background = FALSE',
                                             Cond1_DF_Cutoff_Specific=="DOWN" & Cond2_DF_Cutoff_Specific=="DOWN" ~ 'core_DOWN',#State 1
                                             Cond1_DF_Cutoff_Specific=="DOWN" & Cond2_DF_Cutoff_Specific=="Not Detected" ~ 'Cond1_DOWN',#State 2
                                             Cond1_DF_Cutoff_Specific=="DOWN" & Cond2_DF_Cutoff_Specific=="Not Significant" ~ 'Cond1_DOWN',#State 3
                                             Cond1_DF_Cutoff_Specific=="DOWN" & Cond2_DF_Cutoff_Specific=="Significant Negative" ~ 'Cond1_DOWN',#State 4
                                             Cond1_DF_Cutoff_Specific=="DOWN" & Cond2_DF_Cutoff_Specific=="Significant Positive" ~ 'Cond1_DOWN',#State 5
                                             Cond1_DF_Cutoff_Specific=="DOWN" & Cond2_DF_Cutoff_Specific=="UP" ~ 'Opposite',#State 6

                                             Cond1_DF_Cutoff_Specific=="UP" & Cond2_DF_Cutoff_Specific=="DOWN" ~ 'Opposite',#State 12
                                             Cond1_DF_Cutoff_Specific=="UP" & Cond2_DF_Cutoff_Specific=="Not Detected" ~ 'Cond1_UP',#State 13
                                             Cond1_DF_Cutoff_Specific=="UP" & Cond2_DF_Cutoff_Specific=="Not Significant" ~ 'Cond1_UP',#State 14
                                             Cond1_DF_Cutoff_Specific=="UP" & Cond2_DF_Cutoff_Specific=="Significant Negative" ~ 'Cond1_UP',#State 15
                                             Cond1_DF_Cutoff_Specific=="UP" & Cond2_DF_Cutoff_Specific=="Significant Positive" ~ 'Cond1_UP',#State 16
                                             Cond1_DF_Cutoff_Specific=="UP" & Cond2_DF_Cutoff_Specific=="UP" ~ 'core_UP',#State 17

                                             Cond1_DF_Cutoff_Specific=="Not Detected" & Cond2_DF_Cutoff_Specific=="DOWN" ~ 'Cond2_DOWN',#State 12
                                             Cond1_DF_Cutoff_Specific=="Not Detected" & Cond2_DF_Cutoff_Specific=="Not Detected" ~ 'None',#State 13
                                             Cond1_DF_Cutoff_Specific=="Not Detected" & Cond2_DF_Cutoff_Specific=="Not Significant" ~ 'None',#State 14
                                             Cond1_DF_Cutoff_Specific=="Not Detected" & Cond2_DF_Cutoff_Specific=="Significant Negative" ~ 'None',#State 15
                                             Cond1_DF_Cutoff_Specific=="Not Detected" & Cond2_DF_Cutoff_Specific=="Significant Positive" ~ 'None',#State 16
                                             Cond1_DF_Cutoff_Specific=="Not Detected" & Cond2_DF_Cutoff_Specific=="UP" ~ 'Cond2_UP',#State 17

                                             Cond1_DF_Cutoff_Specific=="Significant Negative" & Cond2_DF_Cutoff_Specific=="DOWN" ~ 'Cond2_DOWN',#State 12
                                             Cond1_DF_Cutoff_Specific=="Significant Negative" & Cond2_DF_Cutoff_Specific=="Not Detected" ~ 'None',#State 13
                                             Cond1_DF_Cutoff_Specific=="Significant Negative" & Cond2_DF_Cutoff_Specific=="Not Significant" ~ 'None',#State 14
                                             Cond1_DF_Cutoff_Specific=="Significant Negative" & Cond2_DF_Cutoff_Specific=="Significant Negative" ~ 'None',#State 15
                                             Cond1_DF_Cutoff_Specific=="Significant Negative" & Cond2_DF_Cutoff_Specific=="Significant Positive" ~ 'None',#State 16
                                             Cond1_DF_Cutoff_Specific=="Significant Negative" & Cond2_DF_Cutoff_Specific=="UP" ~ 'Cond2_UP',#State 17

                                             Cond1_DF_Cutoff_Specific=="Significant Positive" & Cond2_DF_Cutoff_Specific=="DOWN" ~ 'Cond2_DOWN',#State 12
                                             Cond1_DF_Cutoff_Specific=="Significant Positive" & Cond2_DF_Cutoff_Specific=="Not Detected" ~ 'None',#State 13
                                             Cond1_DF_Cutoff_Specific=="Significant Positive" & Cond2_DF_Cutoff_Specific=="Not Significant" ~ 'None',#State 14
                                             Cond1_DF_Cutoff_Specific=="Significant Positive" & Cond2_DF_Cutoff_Specific=="Significant Negative" ~ 'None',#State 15
                                             Cond1_DF_Cutoff_Specific=="Significant Positive" & Cond2_DF_Cutoff_Specific=="Significant Positive" ~ 'None',#State 16
                                             Cond1_DF_Cutoff_Specific=="Significant Positive" & Cond2_DF_Cutoff_Specific=="UP" ~ 'Cond2_UP',#State 17

                                             Cond1_DF_Cutoff_Specific=="Not Significant" & Cond2_DF_Cutoff_Specific=="DOWN" ~ 'Cond2_DOWN',#State 12
                                             Cond1_DF_Cutoff_Specific=="Not Significant" & Cond2_DF_Cutoff_Specific=="Not Detected" ~ 'None',#State 13
                                             Cond1_DF_Cutoff_Specific=="Not Significant" & Cond2_DF_Cutoff_Specific=="Not Significant" ~ 'None',#State 14
                                             Cond1_DF_Cutoff_Specific=="Not Significant" & Cond2_DF_Cutoff_Specific=="Significant Negative" ~ 'None',#State 15
                                             Cond1_DF_Cutoff_Specific=="Not Significant" & Cond2_DF_Cutoff_Specific=="Significant Positive" ~ 'None',#State 16
                                             Cond1_DF_Cutoff_Specific=="Not Significant" & Cond2_DF_Cutoff_Specific=="UP" ~ 'Cond2_UP',#State 1
                                             TRUE ~ 'NA'))

  #Safe the DF and return the groupings
  ##MCA DF (Merged InputDF filtered for background with assigned MCA cluster names)
  MergeDF_Select1 <- MergeDF[, c("MetaboliteID", "Cond1_DF_Detected","Cond1_DF_ValueCol","Cond1_DF_PadjCol","Cond1_DF_Cutoff", "Cond1_DF_Cutoff_Specific", "Cond2_DF_Detected", "Cond2_DF_ValueCol","Cond2_DF_PadjCol","Cond2_DF_Cutoff", "Cond2_DF_Cutoff_Specific", "BG_method", "RG1_All", "RG2_Significant", "RG3_SignificantChange")]

  Cond2ValueCol_Unique<-paste("Cond2_DF_",metadata_info_c2[["ValueCol"]])
  Cond2PadjCol_Unique <-paste("Cond2_DF_",metadata_info_c2[["StatCol"]])
  Cond1ValueCol_Unique<-paste("Cond1_DF_",metadata_info_c1[["ValueCol"]])
  Cond1PadjCol_Unique <-paste("Cond1_DF_",metadata_info_c1[["StatCol"]])

  MergeDF_Select2<- subset(MergeDF, select=-c(Cond1_DF_Detected,Cond1_DF_Cutoff, Cond2_DF_Detected,Cond2_DF_Cutoff, Cond2_DF_Cutoff_Specific, BG_method, RG1_All, RG2_Significant, RG3_SignificantChange))%>%
    dplyr::rename(!!Cond2ValueCol_Unique :="Cond2_DF_ValueCol",#This syntax is needed since paste(MetaboliteID)="MetaboliteID" is not working in dyplr
                  !!Cond2PadjCol_Unique :="Cond2_DF_PadjCol",
                  !!Cond1ValueCol_Unique :="Cond1_DF_ValueCol",
                  !!Cond1PadjCol_Unique :="Cond1_DF_PadjCol")

  MergeDF_Rearrange <- merge(MergeDF_Select1, MergeDF_Select2, by="MetaboliteID")
  MergeDF_Rearrange <-MergeDF_Rearrange%>%
    dplyr::rename("Metabolite"="MetaboliteID")

  ##summary SiRCle clusters (number of genes assigned to each SiRCle cluster in each grouping)
  Clustersummary_RG1 <- MergeDF_Rearrange[,c("Metabolite", "RG1_All")]%>%
    dplyr::count(RG1_All, name="Number of Features")%>%
    dplyr::rename("SiRCle cluster Name"= "RG1_All")
  Clustersummary_RG1$`Regulation Grouping` <- "RG1_All"

  Clustersummary_RG2 <- MergeDF_Rearrange[,c("Metabolite", "RG2_Significant")]%>%
    dplyr::count(RG2_Significant, name="Number of Features")%>%
    dplyr::rename("SiRCle cluster Name"= "RG2_Significant")
  Clustersummary_RG2$`Regulation Grouping` <- "RG2_Significant"

  Clustersummary_RG3 <- MergeDF_Rearrange[,c("Metabolite", "RG3_SignificantChange")]%>%
    dplyr::count(RG3_SignificantChange, name="Number of Features")%>%
    dplyr::rename("SiRCle cluster Name"= "RG3_SignificantChange")
  Clustersummary_RG3$`Regulation Grouping` <- "RG3_SignificantChange"

  Clustersummary <- rbind(Clustersummary_RG1, Clustersummary_RG2,Clustersummary_RG3)
  Clustersummary <- Clustersummary[,c(3,1,2)]

  ## Rename feature
  MergeDF_Rearrange <-MergeDF_Rearrange%>%
    dplyr::rename(!!feature := "Metabolite")

  ######################################################################################################################################################################
  ##----- Save and Return
  #Here we make a list in which we will save the outputs:
  DF_List <- list("MCA_2Cond_summary"=Clustersummary, "MCA_2Cond_Results"=MergeDF_Rearrange)

  suppressMessages(suppressWarnings(
    save_res(inputlist_df=DF_List,
                         inputlist_plot= NULL,
                         save_table=save_table,
                         save_plot=NULL,
                         path= folder,
                         file_name= "mca_2cond",
                         core=FALSE,
                         print_plot=FALSE)))

  #Return:
  invisible(return(DF_List))
}


######################################################
### ### ### Metabolite Clustering Analysis ### ### ###
######################################################

#' This script performs metabolite clustering analysis and computes clusters of metabolites based on regulatory rules between Intracellular and culture media metabolomics (core experiment).
#'
#' @param data_intra DF for your data (results from e.g. dma) containing metabolites in rows with corresponding Log2FC and stat (p-value, p.adjusted) value columns.
#' @param data_core DF for your data (results from e.g. dma) containing metabolites in rows with corresponding Log2FC and stat (p-value, p.adjusted) value columns. Here we additionally require
#' @param metadata_info_intra  \emph{Optional: } Pass ColumnNames and Cutoffs for the intracellular metabolomics including the value column (e.g. Log2FC, Log2Diff, t.val, etc) and the stats column (e.g. p.adj, p.val). This must include: c(ValueCol=ColumnName_data_intra,StatCol=ColumnName_data_intra, cutoff_stat= NumericValue, ValueCutoff=NumericValue) \strong{Default=c(ValueCol="Log2FC",StatCol="p.adj", cutoff_stat= 0.05, ValueCutoff=1)}
#' @param metadata_info_core  \emph{Optional: } Pass ColumnNames and Cutoffs for the consumption-release metabolomics including the direction column, the value column (e.g. Log2Diff, t.val, etc) and the stats column (e.g. p.adj, p.val). This must include: c(DirectionCol= ColumnName_data_core,ValueCol=ColumnName_data_core,StatCol=ColumnName_data_core, cutoff_stat= NumericValue, ValueCutoff=NumericValue)\strong{Default=c(DirectionCol="core", ValueCol="Log2(Distance)",StatCol="p.adj", cutoff_stat= 0.05, ValueCutoff=1)}
#' @param feature \emph{Optional: } Column name of Column including the Metabolite identifiers. This MUST BE THE SAME in each of your Input files. \strong{Default="Metabolite"}
#' @param method_background \emph{Optional: } Background method `Intra|core, Intra&core, core, Intra or * \strong{Default="Intra&core"}
#' @param save_table \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt" \strong{default: "csv"}
#' @param path \emph{Optional:} Path to the folder the results should be saved at. \strong{default: NULL}
#'
#' @return List of two DFs: 1. summary of the cluster count and 2. the detailed information of each metabolites in the clusters.
#'
#' @examples
#'
#' Media <- medium_raw %>%tibble::column_to_rownames("Code")
#' ResM <- MetaProViz::processing(data = Media[-c(40:45) ,-c(1:3)],
#'                                   metadata_sample = Media[-c(40:45) ,c(1:3)] ,
#'                                   metadata_info = c(Conditions = "Conditions", Biological_Replicates = "Biological_Replicates", core_norm_factor = "GrowthFactor", core_media = "blank"),
#'                                   core=TRUE)
#'
#' MediaDMA <- MetaProViz::dma(data=ResM[["DF"]][["Preprocessing_output"]][ ,-c(1:4)],
#'                             metadata_sample=ResM[["DF"]][["Preprocessing_output"]][ , c(1:4)],
#'                             metadata_info = c(Conditions = "Conditions", Numerator = NULL, Denominator  = "HK2"),
#'                             pval ="aov",
#'                             core=TRUE)
#'
#' IntraDMA <- intracell_dma
#'
#' Res <- MetaProViz::mca_core(
#'     data_intra = as.data.frame(IntraDMA),
#'     data_core = as.data.frame(MediaDMA[["dma"]][["786-M1A_vs_HK2"]])
#'   )
#'
#' @keywords biological clustering
#'
#' @importFrom dplyr rename mutate case_when mutate_at count
#' @importFrom tidyr replace_na
#' @importFrom magrittr %>% %<>%
#' @importFrom tibble rownames_to_column column_to_rownames
#'
#' @export
#'
mca_core <- function(data_intra,
                     data_core,
                     metadata_info_intra=c(ValueCol="Log2FC",StatCol="p.adj", cutoff_stat= 0.05, ValueCutoff=1),
                     metadata_info_core=c(DirectionCol="core", ValueCol="Log2(Distance)",StatCol="p.adj", cutoff_stat= 0.05, ValueCutoff=1),
                     feature= "Metabolite",
                     save_table = "csv",
                     method_background="Intra&core",
                     path=NULL
                     ){

  ## ------------ Create log file ----------- ##
  metaproviz_init()

  ################################################################################################################################################################################################
  ## ------------ Check Input files ----------- ##
  check_param_mca(
      data_c1=NULL,
      data_c2=NULL,
      data_core= data_core,
      data_intra= data_intra,
      metadata_info_c1=NULL,
      metadata_info_c2=NULL,
      metadata_info_core=metadata_info_core,
      metadata_info_intra=metadata_info_intra,
      method_background=method_background,
      feature=feature,
      save_table=save_table
  )


  ## ------------ Create Results output folder ----------- ##
  if(is.null(save_table)==FALSE){
    folder <- save_path(folder_name= "MCAcore",
                                    path=path)
  }


  ################################################################################################################################################################################################
  ## ------------ Prepare the Input -------- ##
  #Import the data and check columns (here the user will get an error if the column can not be renamed as it does not exists.)
  core_DF <-  as.data.frame(data_core)%>%
    dplyr::rename("MetaboliteID"=paste(feature),
                  "ValueCol"=metadata_info_core[["ValueCol"]],
                  "PadjCol"=metadata_info_core[["StatCol"]],
                  "core_Direction"=metadata_info_core[["DirectionCol"]])
  core_DF <- core_DF[complete.cases(core_DF$ValueCol, core_DF$PadjCol), ]

  Intra_DF<- as.data.frame(data_intra)%>%
    dplyr::rename("MetaboliteID"=paste(feature),
                  "ValueCol"=metadata_info_intra[["ValueCol"]],
                  "PadjCol"=metadata_info_intra[["StatCol"]])
  Intra_DF <- Intra_DF[complete.cases(Intra_DF$ValueCol, Intra_DF$PadjCol), ]

  #Tag genes that are detected in each data layer
  core_DF$Detected <- "TRUE"
  Intra_DF$Detected <- "TRUE"

  #Assign to Group based on individual Cutoff ("UP", "DOWN", "No Change")
  core_DF <- core_DF%>%
    dplyr::mutate(Cutoff = dplyr::case_when(core_DF$PadjCol <= as.numeric(metadata_info_core[["cutoff_stat"]]) & core_DF$ValueCol > as.numeric(metadata_info_core[["ValueCutoff"]]) ~ 'UP',
                              core_DF$PadjCol <= as.numeric(metadata_info_core[["cutoff_stat"]]) & core_DF$ValueCol < - as.numeric(metadata_info_core[["ValueCutoff"]]) ~ 'DOWN',
                              TRUE ~ 'No Change')) %>%
    dplyr::mutate(Cutoff_Specific = dplyr::case_when(Cutoff == "UP" ~ 'UP',
                                       Cutoff == "DOWN" ~ 'DOWN',
                                       Cutoff == "No Change" & core_DF$PadjCol <= as.numeric(metadata_info_core[["cutoff_stat"]]) & core_DF$ValueCol > 0 ~ 'Significant Positive',
                                       Cutoff == "No Change" & core_DF$PadjCol <= as.numeric(metadata_info_core[["cutoff_stat"]]) & core_DF$ValueCol < 0 ~ 'Significant Negative',
                                       Cutoff == "No Change" & core_DF$PadjCol > as.numeric(metadata_info_core[["cutoff_stat"]]) ~ 'Not Significant',
                                       TRUE ~ 'NA'))

  Intra_DF <-Intra_DF%>%
    dplyr::mutate(Cutoff = dplyr::case_when(Intra_DF$PadjCol <= as.numeric(metadata_info_intra[["cutoff_stat"]]) & Intra_DF$ValueCol > as.numeric(metadata_info_intra[["ValueCutoff"]]) ~ 'UP',
                              Intra_DF$PadjCol <= as.numeric(metadata_info_intra[["cutoff_stat"]]) & Intra_DF$ValueCol < - as.numeric(metadata_info_intra[["ValueCutoff"]]) ~ 'DOWN',
                              TRUE ~ 'No Change'))%>%
    dplyr::mutate(Cutoff_Specific = dplyr::case_when(Cutoff == "UP" ~ 'UP',
                                       Cutoff == "DOWN" ~ 'DOWN',
                                       Cutoff == "No Change" & Intra_DF$PadjCol <= as.numeric(metadata_info_intra[["cutoff_stat"]]) & Intra_DF$ValueCol > 0 ~ 'Significant Positive',
                                       Cutoff == "No Change" & Intra_DF$PadjCol <= as.numeric(metadata_info_intra[["cutoff_stat"]]) & Intra_DF$ValueCol < 0 ~ 'Significant Negative',
                                       Cutoff == "No Change" & Intra_DF$PadjCol > as.numeric(metadata_info_intra[["cutoff_stat"]]) ~ 'Not Significant',
                                       TRUE ~ 'NA'))

  #Merge the dataframes together: Merge the supplied Intra and core dataframes together.
  ##Add prefix to column names to distinguish the different data types after merge
  colnames(core_DF) <- paste0("core_DF_", colnames(core_DF))
  core_DF <- core_DF%>%
    dplyr::rename("MetaboliteID" = "core_DF_MetaboliteID")

  colnames(Intra_DF) <- paste0("Intra_DF_", colnames(Intra_DF))
  Intra_DF <-Intra_DF%>%
    dplyr::rename("MetaboliteID"="Intra_DF_MetaboliteID")

  ##Merge
  MergeDF <- merge(Intra_DF, core_DF, by="MetaboliteID", all=TRUE)

  ##Mark the undetected genes in each data layer
  MergeDF<-MergeDF %>%
    dplyr::mutate_at(c("core_DF_Detected","Intra_DF_Detected"), ~tidyr::replace_na(.,"FALSE"))%>%
    dplyr::mutate_at(c("core_DF_Cutoff","Intra_DF_Cutoff"), ~tidyr::replace_na(.,"No Change"))%>%
    dplyr::mutate_at(c("core_DF_Cutoff_Specific", "Intra_DF_Cutoff_Specific"), ~tidyr::replace_na(.,"Not Detected"))%>%
    dplyr::mutate_at(c("core_DF_core_Direction"), ~tidyr::replace_na(.,"Not Detected"))

  #Apply Background filter (label metabolites that will be removed based on chosen background)
  if(method_background == "Intra|core"){# C1|C2 = core OR Intra
    MergeDF <- MergeDF%>%
      dplyr::mutate(BG_method = dplyr::case_when(Intra_DF_Detected=="TRUE" & core_DF_Detected=="TRUE" ~ 'TRUE', #Intra & core
                                   Intra_DF_Detected=="TRUE" & core_DF_Detected=="FALSE" ~ 'TRUE', # JustIntra
                                   Intra_DF_Detected=="FALSE" & core_DF_Detected=="TRUE" ~ 'TRUE', # Just core
                                   TRUE ~ 'FALSE'))
  }else if(method_background == "Intra&core"){ # core AND Intra
    MergeDF <- MergeDF%>%
      dplyr::mutate(BG_method = dplyr::case_when(Intra_DF_Detected=="TRUE" & core_DF_Detected=="TRUE" ~ 'TRUE', #Intra & core
                                   TRUE ~ 'FALSE'))
  }else if(method_background == "core"){ # core has to be there
    MergeDF <- MergeDF%>%
      dplyr::mutate(BG_method = dplyr::case_when(Intra_DF_Detected=="TRUE" & core_DF_Detected=="TRUE" ~ 'TRUE', #Intra & core
                                   Intra_DF_Detected=="FALSE" & core_DF_Detected=="TRUE" ~ 'TRUE', # Just core
                                   TRUE ~ 'FALSE'))
  }else if(method_background == "Intra"){ #Intra has to be there
    MergeDF <- MergeDF%>%
      dplyr::mutate(BG_method = dplyr::case_when(Intra_DF_Detected=="TRUE" & core_DF_Detected=="TRUE" ~ 'TRUE', #Intra & core
                                   Intra_DF_Detected=="TRUE" & core_DF_Detected=="FALSE" ~ 'TRUE', # JustIntra
                                   TRUE ~ 'FALSE'))
  }else if(method_background == "*"){ # Use all metabolites as the background
    MergeDF$BG_method <- "TRUE"
  }else{
    stop("Please use one of the following method_backgrounds: Intra|core, Intra&core, core, Intra, *")#error message
  }

  #Assign Metabolite cluster names to the metabolites
  MergeDF <- MergeDF%>%
        dplyr::mutate(RG1_All = dplyr::case_when(BG_method =="FALSE"~ 'Background = FALSE',
                               Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Released" ~ 'Intra DOWN + core DOWN_Released',
                               Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected"~ 'Intra DOWN + core Not Detected',
                               Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Released"~ 'Intra DOWN + core Not Significant_Released',
                               Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Released"~ 'Intra DOWN + core Significant Negative_Released',
                               Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Released"~ 'Intra DOWN + core Significant Positive_Released',
                               Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Released" ~ 'Intra DOWN + core UP_Released',

                               Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Released" ~ 'Intra UP + core DOWN_Released',
                               Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected" ~ 'Intra UP + core Not Detected',
                               Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Released" ~ 'Intra UP + core Not Significant_Released',
                               Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Released" ~ 'Intra UP + core Significant Negative_Released',
                               Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Released" ~ 'Intra UP + core Significant Positive_Released',
                               Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Released" ~ 'Intra UP + core UP_Released',

                               Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Released" ~ 'Intra Not Detected + core DOWN_Released',
                               Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected" ~ 'Intra Not Detected + core Not Detected',
                               Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Released" ~ 'Intra Not Detected + core Not Significant_Released',
                               Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Released" ~ 'Intra Not Detected + core Significant Negative_Released',
                               Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Released" ~ 'Intra Not Detected + core Significant Positive_Released',
                               Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Released" ~ 'Intra Not Detected + core UP_Released',

                               Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Released"~ 'Intra Significant Negative + core DOWN_Released',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected"~ 'Intra Significant Negative + core Not Detected',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Released"~ 'Intra Significant Negative + core Not Significant_Released',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Released"~ 'Intra Significant Negative + core Significant Negative_Released',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Released"~ 'Intra Significant Negative + core Significant Positive_Released',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Released"~ 'Intra Significant Negative + core UP_Released',

                               Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Released"~ 'Intra Significant Positive + core DOWN_Released',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected"~ 'Intra Significant Positive + core Not Detected',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Released"~ 'Intra Significant Positive + core Not Significant_Released',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Released"~ 'Intra Significant Positive + core Significant Negative_Released',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Released"~ 'Intra Significant Positive + core Significant Positive_Released',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Released"~ 'Intra Significant Positive + core UP_Released',

                               Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Released"~ 'Intra Not Significant + core DOWN_Released',
                               Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected"~ 'Intra Not Significant + core Not Detected',
                               Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Released"~ 'Intra Not Significant + core Not Significant_Released',
                               Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Released"~ 'Intra Not Significant + core Significant Negative_Released',
                               Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Released"~ 'Intra Not Significant + core Significant Positive_Released',
                               Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Released"~ 'Intra Not Significant + core UP_Released',

                               #Consumed:
                               Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Consumed" ~ 'Intra DOWN + core DOWN_Consumed',
                               Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected"~ 'Intra DOWN + core Not Detected',
                               Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Consumed"~ 'Intra DOWN + core Not Significant_Consumed',
                               Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Consumed"~ 'Intra DOWN + core Significant Negative_Consumed',
                               Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Consumed"~ 'Intra DOWN + core Significant Positive_Consumed',
                               Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Consumed" ~ 'Intra DOWN + core UP_Consumed',

                               Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Consumed" ~ 'Intra UP + core DOWN_Consumed',
                               Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected" ~ 'Intra UP + core Not Detected',
                               Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Consumed" ~ 'Intra UP + core Not Significant_Consumed',
                               Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Consumed" ~ 'Intra UP + core Significant Negative_Consumed',
                               Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Consumed" ~ 'Intra UP + core Significant Positive_Consumed',
                               Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Consumed" ~ 'Intra UP + core UP_Consumed',

                               Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Consumed" ~ 'Intra Not Detected + core DOWN_Consumed',
                               Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected" ~ 'Intra Not Detected + core Not Detected',
                               Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Consumed" ~ 'Intra Not Detected + core Not Significant_Consumed',
                               Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Consumed" ~ 'Intra Not Detected + core Significant Negative_Consumed',
                               Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Consumed" ~ 'Intra Not Detected + core Significant Positive_Consumed',
                               Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Consumed" ~ 'Intra Not Detected + core UP_Consumed',

                               Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Consumed"~ 'Intra Significant Negative + core DOWN_Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected"~ 'Intra Significant Negative + core Not Detected',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Consumed"~ 'Intra Significant Negative + core Not Significant_Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Consumed"~ 'Intra Significant Negative + core Significant Negative_Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Consumed"~ 'Intra Significant Negative + core Significant Positive_Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Consumed"~ 'Intra Significant Negative + core UP_Consumed',

                               Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Consumed"~ 'Intra Significant Positive + core DOWN_Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected"~ 'Intra Significant Positive + core Not Detected',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Consumed"~ 'Intra Significant Positive + core Not Significant_Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Consumed"~ 'Intra Significant Positive + core Significant Negative_Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Consumed"~ 'Intra Significant Positive + core Significant Positive_Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Consumed"~ 'Intra Significant Positive + core UP_Consumed',

                               Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Consumed"~ 'Intra Not Significant + core DOWN_Consumed',
                               Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected"~ 'Intra Not Significant + core Not Detected',
                               Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Consumed"~ 'Intra Not Significant + core Not Significant_Consumed',
                               Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Consumed"~ 'Intra Not Significant + core Significant Negative_Consumed',
                               Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Consumed"~ 'Intra Not Significant + core Significant Positive_Consumed',
                               Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Consumed"~ 'Intra Not Significant + core UP_Consumed',

                               #Released/Consumed (Consumed in one, released in the other)
                               Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Released/Consumed"~ 'Intra DOWN + core DOWN_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected"~ 'Intra DOWN + core Not Detected',
                               Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Released/Consumed"~ 'Intra DOWN + core Not Significant_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Released/Consumed"~ 'Intra DOWN + core Significant Negative_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Released/Consumed"~ 'Intra DOWN + core Significant Positive_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Released/Consumed"~ 'Intra DOWN + core UP_Released/Consumed',

                               Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Released/Consumed"~ 'Intra UP + core DOWN_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected" ~ 'Intra UP + core Not Detected',
                               Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Released/Consumed"~ 'Intra UP + core Not Significant_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Released/Consumed"~ 'Intra UP + core Significant Negative_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Released/Consumed"~ 'Intra UP + core Significant Positive_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Released/Consumed"~ 'Intra UP + core UP_Released/Consumed',

                               Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Released/Consumed"~ 'Intra Not Detected + core DOWN_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected" ~ 'Intra Not Detected + core Not Detected',
                               Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Released/Consumed"~ 'Intra Not Detected + core Not Significant_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Released/Consumed"~ 'Intra Not Detected + core Significant Negative_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Released/Consumed"~ 'Intra Not Detected + core Significant Positive_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Released/Consumed"~ 'Intra Not Detected + core UP_Released/Consumed',

                               Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Released/Consumed"~ 'Intra Significant Negative + core DOWN_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected"~ 'Intra Significant Negative + core Not Detected',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Released/Consumed"~ 'Intra Significant Negative + core Not Significant_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Released/Consumed"~ 'Intra Significant Negative + core Significant Negative_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Released/Consumed"~ 'Intra Significant Negative + core Significant Positive_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Released/Consumed"~ 'Intra Significant Negative + core UP_Released/Consumed',

                               Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Released/Consumed"~ 'Intra Significant Positive + core DOWN_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected"~ 'Intra Significant Positive + core Not Detected',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Released/Consumed"~ 'Intra Significant Positive + core Not Significant_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Released/Consumed"~ 'Intra Significant Positive + core Significant Negative_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Released/Consumed"~ 'Intra Significant Positive + core Significant Positive_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Released/Consumed"~ 'Intra Significant Positive + core UP_Released/Consumed',

                               Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Released/Consumed"~ 'Intra Not Significant + core DOWN_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected"~ 'Intra Not Significant + core Not Detected',
                               Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Released/Consumed"~ 'Intra Not Significant + core Not Significant_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Released/Consumed"~ 'Intra Not Significant + core Significant Negative_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Released/Consumed"~ 'Intra Not Significant + core Significant Positive_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Released/Consumed"~ 'Intra Not Significant + core UP_Released/Consumed',
                               TRUE ~ 'NA'))%>%
    dplyr::mutate(RG2_Significant = dplyr::case_when(BG_method =="FALSE"~ 'Background = FALSE',
                                       Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Released" ~ 'Both_DOWN (Released)',
                                       Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected"~ 'None',
                                       Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Released"~ 'None',
                                       Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Released"~ 'Both_DOWN (Released)',
                                       Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Released"~ 'Opposite (Released UP)',
                                       Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Released" ~ 'Opposite (Released UP)',

                                       Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Released" ~ 'Opposite (Released DOWN)',
                                       Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Released" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Released" ~ 'Opposite (Released UP)',
                                       Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Released" ~ 'Both_UP (Released)',
                                       Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Released" ~ 'Both_UP (Released)',

                                       Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Released" ~ 'core_DOWN (Released)',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Released" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Released" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Released" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Released" ~ 'core_UP (Released)',

                                       Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Released"~ 'Both_DOWN (Released)',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Released"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Released"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Released"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Released"~ 'Opposite (Released UP)',

                                       Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Released"~ 'Opposite (Released DOWN)',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Released"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Released"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Released"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Released"~ 'Both_UP (Released)',

                                       Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Released"~ 'core_DOWN (Released)',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Released"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Released"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Released"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Released"~ 'core_UP (Released)',

                                       #Consumed:
                                       Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Consumed" ~ 'Both_DOWN (Consumed)',
                                       Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected"~ 'None',
                                       Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Consumed"~ 'Both_DOWN (Consumed)',
                                       Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Consumed"~ 'Opposite (Consumed UP)',
                                       Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Consumed" ~ 'Opposite (Consumed UP)',

                                       Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Consumed" ~ 'Opposite (Consumed DOWN)',
                                       Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Consumed" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Consumed" ~ 'Opposite (Consumed UP)',
                                       Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Consumed" ~ 'Both_UP (Consumed)',
                                       Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Consumed" ~ 'Both_UP (Consumed)',

                                       Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Consumed" ~ 'core_DOWN (Consumed)',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Consumed" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Consumed" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Consumed" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Consumed" ~ 'core_UP (Consumed)',

                                       Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Consumed"~ 'Both_DOWN (Consumed)',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Consumed"~ 'Opposite (Consumed UP)',

                                       Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Consumed"~ 'Opposite (Consumed DOWN)',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Consumed"~ 'Both_UP (Consumed)',

                                       Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Consumed"~ 'core_DOWN (Consumed)',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Consumed"~ 'core_UP (Consumed)',

                                       #Consumed/Released (Consumed in one, released in the other)
                                       Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Released/Consumed" ~ 'Both_DOWN (Released/Consumed)',
                                       Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected"~ 'None',
                                       Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Released/Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Released/Consumed"~ 'Both_DOWN (Released/Consumed)',
                                       Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Released/Consumed"~ 'Opposite (Released/Consumed UP)',
                                       Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Released/Consumed" ~ 'Opposite (Released/Consumed UP)',

                                       Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Released/Consumed" ~ 'Opposite (Released/Consumed DOWN)',
                                       Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Released/Consumed" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Released/Consumed" ~ 'Opposite (Released/Consumed UP)',
                                       Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Released/Consumed" ~ 'Both_UP (Released/Consumed)',
                                       Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Released/Consumed" ~ 'Both_UP (Released/Consumed)',

                                       Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Released/Consumed" ~ 'core_DOWN (Released/Consumed)',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Released/Consumed" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Released/Consumed" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Released/Consumed" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Released/Consumed" ~ 'core_UP (Released/Consumed)',

                                       Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Released/Consumed"~ 'Both_DOWN (Released/Consumed)',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Released/Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Released/Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Released/Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Released/Consumed"~ 'Opposite (Released/Consumed UP)',

                                       Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Released/Consumed"~ 'Opposite (Released/Consumed DOWN)',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Released/Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Released/Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Released/Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Released/Consumed"~ 'Both_UP (Released/Consumed)',

                                       Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Released/Consumed"~ 'core_DOWN (Released/Consumed)',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Released/Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Released/Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Released/Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Released/Consumed"~ 'core_UP (Released/Consumed)',
                                       TRUE ~ 'NA'))%>%
    dplyr::mutate(RG3_Change = dplyr::case_when(BG_method =="FALSE"~ 'Background = FALSE',
                                  Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Released" ~ 'Both_DOWN (Released)',
                                  Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected"~ 'None',
                                  Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Released"~ 'None',
                                  Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Released"~ 'None',
                                  Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Released"~ 'None',
                                  Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Released" ~ 'Opposite (Released UP)',

                                  Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Released" ~ 'Opposite (Released DOWN)',
                                  Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Released" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Released" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Released" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Released" ~ 'Both_UP (Released)',

                                  Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Released" ~ 'core_DOWN (Released)',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Released" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Released" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Released" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Released" ~ 'core_UP (Released)',

                                  Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Released"~ 'core_DOWN (Released)',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Released"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Released"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Released"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Released"~ 'core_UP (Released)',

                                  Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Released"~ 'core_DOWN (Released)',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Released"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Released"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Released"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Released"~ 'core_UP (Released)',

                                  Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Released"~ 'core_DOWN (Released)',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Released"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Released"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Released"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Released"~ 'core_UP (Released)',

                                  #Consumed:
                                  Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Consumed" ~ 'Both_DOWN (Consumed)',
                                  Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected"~ 'None',
                                  Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Consumed" ~ 'Opposite (Consumed UP)',

                                  Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Consumed" ~ 'Opposite (Consumed DOWN)',
                                  Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Consumed" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Consumed" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Consumed" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Consumed" ~ 'Both_UP (Consumed)',

                                  Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Consumed" ~ 'core_DOWN (Consumed)',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Consumed" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Consumed" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Consumed" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Consumed" ~ 'core_UP (Consumed)',

                                  Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Consumed"~ 'core_DOWN (Consumed)',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Consumed"~ 'core_UP (Consumed)',

                                  Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Consumed"~ 'core_DOWN (Consumed)',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Consumed"~ 'core_UP (Consumed)',

                                  Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Consumed"~ 'core_DOWN (Consumed)',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Consumed"~ 'core_UP (Consumed)',

                                  #Consumed/Released (Consumed in one, released in the other)
                                  Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Released/Consumed" ~ 'Both_DOWN (Released/Consumed)',
                                  Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected"~ 'None',
                                  Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Released/Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Released/Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Released/Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="DOWN" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Released/Consumed" ~ 'Opposite (Released/Consumed UP)',

                                  Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Released/Consumed" ~ 'Opposite (Released/Consumed DOWN)',
                                  Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Released/Consumed" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Released/Consumed" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Released/Consumed" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="UP" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Released/Consumed" ~ 'Both_UP (Released/Consumed)',

                                  Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Released/Consumed" ~ 'core_DOWN (Released/Consumed)',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Released/Consumed" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Released/Consumed" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Released/Consumed" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Released/Consumed" ~ 'core_UP (Released/Consumed)',

                                  Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Released/Consumed"~ 'core_DOWN (Released/Consumed)',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Released/Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Released/Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Released/Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Released/Consumed"~ 'core_UP (Released/Consumed)',

                                  Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Released/Consumed"~ 'core_DOWN (Released/Consumed)',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Released/Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Released/Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Released/Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Released/Consumed"~ 'core_UP (Released/Consumed)',

                                  Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="DOWN" & core_DF_core_Direction=="Released/Consumed"~ 'core_DOWN (Released/Consumed)',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="Not Detected" & core_DF_core_Direction=="Not Detected"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="Not Significant" & core_DF_core_Direction=="Released/Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="Significant Negative" & core_DF_core_Direction=="Released/Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="Significant Positive" & core_DF_core_Direction=="Released/Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & core_DF_Cutoff_Specific=="UP" & core_DF_core_Direction=="Released/Consumed"~ 'core_UP (Released/Consumed)',
                                  TRUE ~ 'NA'))

  #Safe the DF and return the groupings
  ##MCA DF (Merged InputDF filtered for background with assigned MCA cluster names)
  MergeDF_Select1 <- MergeDF[, c("MetaboliteID", "Intra_DF_Detected","Intra_DF_ValueCol","Intra_DF_PadjCol","Intra_DF_Cutoff", "Intra_DF_Cutoff_Specific", "core_DF_Detected", "core_DF_ValueCol","core_DF_PadjCol","core_DF_Cutoff", "core_DF_Cutoff_Specific", "BG_method", "RG1_All", "RG2_Significant", "RG3_Change")]

  coreValueCol_Unique<-paste("Cond2_DF_",metadata_info_core[["ValueCol"]])
  corePadjCol_Unique <-paste("Cond2_DF_",metadata_info_core[["StatCol"]])
  IntraValueCol_Unique<-paste("Cond1_DF_",metadata_info_intra[["ValueCol"]])
  IntraPadjCol_Unique <-paste("Cond1_DF_",metadata_info_intra[["StatCol"]])

  MergeDF_Select2<- subset(MergeDF, select=-c(Intra_DF_Detected,Intra_DF_Cutoff, core_DF_Detected,core_DF_Cutoff, core_DF_Cutoff_Specific, BG_method, RG1_All, RG2_Significant, RG3_Change))%>%
    dplyr::rename(!!coreValueCol_Unique :="core_DF_ValueCol",#This syntax is needed since paste(MetaboliteID)="MetaboliteID" is not working in dyplr
                  !!corePadjCol_Unique :="core_DF_PadjCol",
                  !!IntraValueCol_Unique :="Intra_DF_ValueCol",
                  !!IntraPadjCol_Unique :="Intra_DF_PadjCol")

  MergeDF_Rearrange <- merge(MergeDF_Select1, MergeDF_Select2, by="MetaboliteID")
  MergeDF_Rearrange <-MergeDF_Rearrange%>%
    dplyr::rename("Metabolite"="MetaboliteID")


  ##summary SiRCle clusters (number of genes assigned to each SiRCle cluster in each grouping)
  Clustersummary_RG1 <- MergeDF_Rearrange[,c("Metabolite", "RG1_All")]%>%
    dplyr::group_by(RG1_All) %>%
    dplyr::mutate("Number of Features" = n()) %>%
    dplyr::distinct(RG1_All, .keep_all = TRUE) %>%
    dplyr::rename("SiRCle cluster Name" = "RG1_All")
  Clustersummary_RG1$`Regulation Grouping` <- "RG1_All"
  Clustersummary_RG1 <- Clustersummary_RG1[-c(1)]

  Clustersummary_RG2 <- MergeDF_Rearrange[,c("Metabolite", "RG2_Significant")]%>%
    dplyr::group_by(RG2_Significant) %>%
    dplyr::mutate("Number of Features" = n()) %>%
    dplyr::distinct(RG2_Significant, .keep_all = TRUE) %>%
    dplyr::rename("SiRCle cluster Name"= "RG2_Significant")
  Clustersummary_RG2$`Regulation Grouping` <- "RG2_Significant"
  Clustersummary_RG2 <- Clustersummary_RG2[-c(1)]

  Clustersummary_RG3 <- MergeDF_Rearrange[,c("Metabolite", "RG3_Change")]%>%
    dplyr::group_by(RG3_Change) %>%
    dplyr::mutate("Number of Features" = n()) %>%
    dplyr::distinct(RG3_Change, .keep_all = TRUE) %>%
    dplyr::rename("SiRCle cluster Name"= "RG3_Change")
  Clustersummary_RG3$`Regulation Grouping` <- "RG3_Change"
  Clustersummary_RG3 <- Clustersummary_RG3[-c(1)]

  Clustersummary <- rbind(Clustersummary_RG1, Clustersummary_RG2,Clustersummary_RG3)
  Clustersummary <- Clustersummary[,c(3,1,2)]

  ## Rename feature
  MergeDF_Rearrange <-MergeDF_Rearrange%>%
    dplyr::rename(!!feature := "Metabolite")

  ######################################################################################################################################################################
  ##----- Save and Return
  #Here we make a list in which we will save the outputs:
  DF_List <- list("MCA_core_summary"=Clustersummary, "MCA_core_Results"=MergeDF_Rearrange)

  suppressMessages(suppressWarnings(
    save_res(inputlist_df=DF_List,
                         inputlist_plot= NULL,
                         save_table=save_table,
                         save_plot=NULL,
                         path= folder,
                         file_name= "mca_2cond",
                         core=FALSE,
                         print_plot=FALSE)))

  invisible(return(DF_List))
}
