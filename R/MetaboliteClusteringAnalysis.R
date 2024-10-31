#' ## ---------------------------
##
## Script name: MCA
##
## Purpose of script: Metabolite Clustering Analysis generates clusters of metabolites based on regulatory rules or kmeans clustering.
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


#' MCA_2Cond
#'
#' This script allows you to perform metabolite clustering analysis and computes clusters of metabolites based on regulatory rules between two conditions (e.g. KO versus WT in Hypoxia = Cond1 and KO versus WT in Normoxia = Cond2).
#'
#' @param InputData_C1 DF for your data (results from e.g. DMA) containing metabolites in rows with corresponding Log2FC and stat (p-value, p.adjusted) value columns.
#' @param InputData_C2 DF for your data (results from e.g. DMA) containing metabolites in rows with corresponding Log2FC and stat (p-value, p.adjusted) value columns.
#' @param SettingsInfo_C1  \emph{Optional: } Pass ColumnNames and Cutoffs for condition 1 including the value column (e.g. Log2FC, Log2Diff, t.val, etc) and the stats column (e.g. p.adj, p.val). This must include: c(ValueCol=ColumnName_InputData_C1,StatCol=ColumnName_InputData_C1, StatCutoff= NumericValue, ValueCutoff=NumericValue) \strong{Default=c(ValueCol="Log2FC",StatCol="p.adj", StatCutoff= 0.05, ValueCutoff=1)}
#' @param SettingsInfo_C2  \emph{Optional: } Pass ColumnNames and Cutoffs for condition 2 includingthe value column (e.g. Log2FC, Log2Diff, t.val, etc) and the stats column (e.g. p.adj, p.val). This must include: c(ValueCol=ColumnName_InputData_C2,StatCol=ColumnName_InputData_C2, StatCutoff= NumericValue, ValueCutoff=NumericValue)\strong{Default=c(ValueCol="Log2FC",StatCol="p.adj", StatCutoff= 0.05, ValueCutoff=1)}
#' @param FeatureID \emph{Optional: } Column name of Column including the Metabolite identifiers. This MUST BE THE SAME in each of your Input files. \strong{Default="Metabolite"}
#' @param BackgroundMethod \emph{Optional: } Background method C1|C2, C1&C2, C2, C1 or * \strong{Default="C1&C2"}
#' @param SaveAs_Table \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt" \strong{default: "csv"}
#' @param FolderPath \emph{Optional:} Path to the folder the results should be saved at. \strong{default: NULL}
#' @return MCA an instance of the MetaProViz package
#' @export
#'

##################################################
### ### ### Metabolite Clustering Analysis: 2 Conditions ### ### ###
##################################################

MCA_2Cond <- function(InputData_C1,
                      InputData_C2,
                      SettingsInfo_C1=c(ValueCol="Log2FC",StatCol="p.adj", StatCutoff= 0.05, ValueCutoff=1),
                      SettingsInfo_C2=c(ValueCol="Log2FC",StatCol="p.adj", StatCutoff= 0.05, ValueCutoff=1),
                      FeatureID = "Metabolite",
                      SaveAs_Table = "csv",
                      BackgroundMethod="C1&C2",
                      FolderPath=NULL
                      ){
  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0){
    install.packages(new.packages)
  }
  suppressMessages(library(tidyverse))

  ################################################################################################################################################################################################
  ## ------------ Check Input files ----------- ##
  CheckInput_MCA(InputData_C1=InputData_C1,
                              InputData_C2=InputData_C2,
                              InputData_CoRe=NULL,
                              InputData_Intra=NULL,
                              SettingsInfo_C1=SettingsInfo_C1,
                              SettingsInfo_C2=SettingsInfo_C2,
                              SettingsInfo_CoRe=NULL,
                              SettingsInfo_Intra=NULL,
                              BackgroundMethod=BackgroundMethod,
                              FeatureID=FeatureID,
                              SaveAs_Table=SaveAs_Table)

  ## ------------ Create Results output folder ----------- ##
  if(is.null(SaveAs_Table)==FALSE){
    Folder <- SavePath(FolderName= "MCA2Cond",
                                    FolderPath=FolderPath)
  }

  ################################################################################################################################################################################################
  ## ------------ Prepare the Input -------- ##
  #Import the data and check columns (here the user will get an error if the column can not be renamed as it does not exists.)
  Cond1_DF <- as.data.frame(InputData_C1)%>%
    dplyr::rename("MetaboliteID"=paste(FeatureID),
                  "ValueCol"=SettingsInfo_C1[["ValueCol"]],
                  "PadjCol"=SettingsInfo_C1[["StatCol"]])
  Cond1_DF <- Cond1_DF[complete.cases(Cond1_DF$ValueCol, Cond1_DF$PadjCol), ]

  Cond2_DF<- as.data.frame(InputData_C2)%>%
    dplyr::rename("MetaboliteID"=paste(FeatureID),
                  "ValueCol"=SettingsInfo_C2[["ValueCol"]],
                  "PadjCol"=SettingsInfo_C2[["StatCol"]])
  Cond2_DF <- Cond2_DF[complete.cases(Cond2_DF$ValueCol, Cond2_DF$PadjCol), ]

  #Tag genes that are detected in each data layer
  Cond1_DF$Detected <- "TRUE"
  Cond2_DF$Detected <- "TRUE"


  ## ------------ Assign Groups -------- ##
  #Assign to Group based on individual Cutoff ("UP", "DOWN", "No Change")
  Cond1_DF <-Cond1_DF%>%
    mutate(Cutoff = case_when(Cond1_DF$PadjCol <= as.numeric(SettingsInfo_C1[["StatCutoff"]]) & Cond1_DF$ValueCol > as.numeric(SettingsInfo_C1[["ValueCutoff"]]) ~ 'UP',
                              Cond1_DF$PadjCol <= as.numeric(SettingsInfo_C1[["StatCutoff"]]) & Cond1_DF$ValueCol < - as.numeric(SettingsInfo_C1[["ValueCutoff"]]) ~ 'DOWN',
                              TRUE ~ 'No Change'))%>%
    mutate(Cutoff_Specific = case_when(Cutoff == "UP" ~ 'UP',
                                       Cutoff == "DOWN" ~ 'DOWN',
                                       Cutoff == "No Change" & Cond1_DF$PadjCol <= as.numeric(SettingsInfo_C1[["StatCutoff"]]) & Cond1_DF$ValueCol > 0 ~ 'Significant Positive',
                                       Cutoff == "No Change" & Cond1_DF$PadjCol <= as.numeric(SettingsInfo_C1[["StatCutoff"]]) & Cond1_DF$ValueCol < 0 ~ 'Significant Negative',
                                       Cutoff == "No Change" & Cond1_DF$PadjCol > as.numeric(SettingsInfo_C1[["StatCutoff"]]) ~ 'Not Significant',
                                       TRUE ~ 'NA'))

  Cond2_DF <- Cond2_DF%>%
    mutate(Cutoff = case_when(Cond2_DF$PadjCol <= as.numeric(SettingsInfo_C2[["StatCutoff"]]) & Cond2_DF$ValueCol > as.numeric(SettingsInfo_C2[["ValueCutoff"]]) ~ 'UP',
                              Cond2_DF$PadjCol <= as.numeric(SettingsInfo_C2[["StatCutoff"]]) & Cond2_DF$ValueCol < - as.numeric(SettingsInfo_C2[["ValueCutoff"]]) ~ 'DOWN',
                              TRUE ~ 'No Change')) %>%
    mutate(Cutoff_Specific = case_when(Cutoff == "UP" ~ 'UP',
                                       Cutoff == "DOWN" ~ 'DOWN',
                                       Cutoff == "No Change" & Cond2_DF$PadjCol <= as.numeric(SettingsInfo_C2[["StatCutoff"]]) & Cond2_DF$ValueCol > 0 ~ 'Significant Positive',
                                       Cutoff == "No Change" & Cond2_DF$PadjCol <= as.numeric(SettingsInfo_C2[["StatCutoff"]]) & Cond2_DF$ValueCol < 0 ~ 'Significant Negative',
                                       Cutoff == "No Change" & Cond2_DF$PadjCol > as.numeric(SettingsInfo_C2[["StatCutoff"]]) ~ 'Not Significant',
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
    mutate_at(c("Cond2_DF_Detected","Cond1_DF_Detected"), ~replace_na(.,"FALSE"))%>%
    mutate_at(c("Cond2_DF_Cutoff","Cond1_DF_Cutoff"), ~replace_na(.,"No Change"))%>%
    mutate_at(c("Cond2_DF_Cutoff_Specific", "Cond1_DF_Cutoff_Specific"), ~replace_na(.,"Not Detected"))

  #Apply Background filter (label genes that will be removed based on choosen background)
  if(BackgroundMethod == "C1|C2"){# C1|C2 = Cond2 OR Cond1
    MergeDF <- MergeDF%>%
      mutate(BG_Method = case_when(Cond1_DF_Detected=="TRUE" & Cond2_DF_Detected=="TRUE" ~ 'TRUE', #Cond1 & Cond2
                                   Cond1_DF_Detected=="TRUE" & Cond2_DF_Detected=="FALSE" ~ 'TRUE', # JustCond1
                                   Cond1_DF_Detected=="FALSE" & Cond2_DF_Detected=="TRUE" ~ 'TRUE', # Just Cond2
                                   TRUE ~ 'FALSE'))
  }else if(BackgroundMethod == "C1&C2"){ # Cond2 AND Cond1
    MergeDF <- MergeDF%>%
      mutate(BG_Method = case_when(Cond1_DF_Detected=="TRUE" & Cond2_DF_Detected=="TRUE" ~ 'TRUE', #Cond1 & Cond2
                                   TRUE ~ 'FALSE'))
  }else if(BackgroundMethod == "C2"){ # Cond2 has to be there
    MergeDF <- MergeDF%>%
      mutate(BG_Method = case_when(Cond1_DF_Detected=="TRUE" & Cond2_DF_Detected=="TRUE" ~ 'TRUE', #Cond1 & Cond2
                                   Cond1_DF_Detected=="FALSE" & Cond2_DF_Detected=="TRUE" ~ 'TRUE', # Just Cond2
                                   TRUE ~ 'FALSE'))
  }else if(BackgroundMethod == "C1"){ #Cond1 has to be there
    MergeDF <- MergeDF%>%
      mutate(BG_Method = case_when(Cond1_DF_Detected=="TRUE" & Cond2_DF_Detected=="TRUE" ~ 'TRUE', #Cond1 & Cond2
                                   Cond1_DF_Detected=="TRUE" & Cond2_DF_Detected=="FALSE" ~ 'TRUE', # JustCond1
                                   TRUE ~ 'FALSE'))
  }else if(BackgroundMethod == "*"){ # Use all genes as the background
    MergeDF$BG_Method <- "TRUE"
  }else{
    stop("Please use one of the following BackgroundMethods: C1|C2, C1&C2, C2, C1, *")#error message
  }

  #Assign SiRCle cluster names to the genes
  MergeDF <- MergeDF%>%
    mutate(RG1_Specific_Cond2 = case_when(BG_Method =="FALSE"~ 'Background = FALSE',
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
    mutate(RG1_Specific_Cond1 = case_when(BG_Method =="FALSE"~ 'Background = FALSE',
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
    mutate(RG1_All = case_when(BG_Method =="FALSE"~ 'Background = FALSE',
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
    mutate(RG2_Significant = case_when(BG_Method =="FALSE"~ 'Background = FALSE',
                                       Cond1_DF_Cutoff_Specific=="DOWN" & Cond2_DF_Cutoff_Specific=="DOWN" ~ 'Core_DOWN',#State 1
                                       Cond1_DF_Cutoff_Specific=="DOWN" & Cond2_DF_Cutoff_Specific=="Not Detected" ~ 'Cond1_DOWN',#State 2
                                       Cond1_DF_Cutoff_Specific=="DOWN" & Cond2_DF_Cutoff_Specific=="Not Significant" ~ 'Cond1_DOWN',#State 3
                                       Cond1_DF_Cutoff_Specific=="DOWN" & Cond2_DF_Cutoff_Specific=="Significant Negative" ~ 'Core_DOWN',#State 4
                                       Cond1_DF_Cutoff_Specific=="DOWN" & Cond2_DF_Cutoff_Specific=="Significant Positive" ~ 'Opposite',#State 5
                                       Cond1_DF_Cutoff_Specific=="DOWN" & Cond2_DF_Cutoff_Specific=="UP" ~ 'Opposite',#State 6

                                       Cond1_DF_Cutoff_Specific=="UP" & Cond2_DF_Cutoff_Specific=="DOWN" ~ 'Opposite',#State 12
                                       Cond1_DF_Cutoff_Specific=="UP" & Cond2_DF_Cutoff_Specific=="Not Detected" ~ 'Cond1_UP',#State 13
                                       Cond1_DF_Cutoff_Specific=="UP" & Cond2_DF_Cutoff_Specific=="Not Significant" ~ 'Cond1_UP',#State 14
                                       Cond1_DF_Cutoff_Specific=="UP" & Cond2_DF_Cutoff_Specific=="Significant Negative" ~ 'Opposite',#State 15
                                       Cond1_DF_Cutoff_Specific=="UP" & Cond2_DF_Cutoff_Specific=="Significant Positive" ~ 'Core_UP',#State 16
                                       Cond1_DF_Cutoff_Specific=="UP" & Cond2_DF_Cutoff_Specific=="UP" ~ 'Core_UP',#State 17

                                       Cond1_DF_Cutoff_Specific=="Not Detected" & Cond2_DF_Cutoff_Specific=="DOWN" ~ 'Cond2_DOWN',#State 12
                                       Cond1_DF_Cutoff_Specific=="Not Detected" & Cond2_DF_Cutoff_Specific=="Not Detected" ~ 'None',#State 13
                                       Cond1_DF_Cutoff_Specific=="Not Detected" & Cond2_DF_Cutoff_Specific=="Not Significant" ~ 'None',#State 14
                                       Cond1_DF_Cutoff_Specific=="Not Detected" & Cond2_DF_Cutoff_Specific=="Significant Negative" ~ 'None',#State 15
                                       Cond1_DF_Cutoff_Specific=="Not Detected" & Cond2_DF_Cutoff_Specific=="Significant Positive" ~ 'None',#State 16
                                       Cond1_DF_Cutoff_Specific=="Not Detected" & Cond2_DF_Cutoff_Specific=="UP" ~ 'Cond2_UP',#State 17

                                       Cond1_DF_Cutoff_Specific=="Significant Negative" & Cond2_DF_Cutoff_Specific=="DOWN" ~ 'Core_DOWN',#State 12
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
                                       Cond1_DF_Cutoff_Specific=="Significant Positive" & Cond2_DF_Cutoff_Specific=="UP" ~ 'Core_UP',#State 17

                                       Cond1_DF_Cutoff_Specific=="Not Significant" & Cond2_DF_Cutoff_Specific=="DOWN" ~ 'Cond2_DOWN',#State 12
                                       Cond1_DF_Cutoff_Specific=="Not Significant" & Cond2_DF_Cutoff_Specific=="Not Detected" ~ 'None',#State 13
                                       Cond1_DF_Cutoff_Specific=="Not Significant" & Cond2_DF_Cutoff_Specific=="Not Significant" ~ 'None',#State 14
                                       Cond1_DF_Cutoff_Specific=="Not Significant" & Cond2_DF_Cutoff_Specific=="Significant Negative" ~ 'None',#State 15
                                       Cond1_DF_Cutoff_Specific=="Not Significant" & Cond2_DF_Cutoff_Specific=="Significant Positive" ~ 'None',#State 16
                                       Cond1_DF_Cutoff_Specific=="Not Significant" & Cond2_DF_Cutoff_Specific=="UP" ~ 'Cond2_UP',#State 1
                                       TRUE ~ 'NA'))%>%
    mutate(RG3_SignificantChange = case_when(BG_Method =="FALSE"~ 'Background = FALSE',
                                             Cond1_DF_Cutoff_Specific=="DOWN" & Cond2_DF_Cutoff_Specific=="DOWN" ~ 'Core_DOWN',#State 1
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
                                             Cond1_DF_Cutoff_Specific=="UP" & Cond2_DF_Cutoff_Specific=="UP" ~ 'Core_UP',#State 17

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
  MergeDF_Select1 <- MergeDF[, c("MetaboliteID", "Cond1_DF_Detected","Cond1_DF_ValueCol","Cond1_DF_PadjCol","Cond1_DF_Cutoff", "Cond1_DF_Cutoff_Specific", "Cond2_DF_Detected", "Cond2_DF_ValueCol","Cond2_DF_PadjCol","Cond2_DF_Cutoff", "Cond2_DF_Cutoff_Specific", "BG_Method", "RG1_All", "RG2_Significant", "RG3_SignificantChange")]

  Cond2ValueCol_Unique<-paste("Cond2_DF_",SettingsInfo_C2[["ValueCol"]])
  Cond2PadjCol_Unique <-paste("Cond2_DF_",SettingsInfo_C2[["StatCol"]])
  Cond1ValueCol_Unique<-paste("Cond1_DF_",SettingsInfo_C1[["ValueCol"]])
  Cond1PadjCol_Unique <-paste("Cond1_DF_",SettingsInfo_C1[["StatCol"]])

  MergeDF_Select2<- subset(MergeDF, select=-c(Cond1_DF_Detected,Cond1_DF_Cutoff, Cond2_DF_Detected,Cond2_DF_Cutoff, Cond2_DF_Cutoff_Specific, BG_Method, RG1_All, RG2_Significant, RG3_SignificantChange))%>%
    dplyr::rename(!!Cond2ValueCol_Unique :="Cond2_DF_ValueCol",#This syntax is needed since paste(MetaboliteID)="MetaboliteID" is not working in dyplr
                  !!Cond2PadjCol_Unique :="Cond2_DF_PadjCol",
                  !!Cond1ValueCol_Unique :="Cond1_DF_ValueCol",
                  !!Cond1PadjCol_Unique :="Cond1_DF_PadjCol")

  MergeDF_Rearrange <- merge(MergeDF_Select1, MergeDF_Select2, by="MetaboliteID")
  MergeDF_Rearrange <-MergeDF_Rearrange%>%
    dplyr::rename("Metabolite"="MetaboliteID")

  ##Summary SiRCle clusters (number of genes assigned to each SiRCle cluster in each grouping)
  ClusterSummary_RG1 <- MergeDF_Rearrange[,c("Metabolite", "RG1_All")]%>%
    count(RG1_All, name="Number of Features")%>%
    dplyr::rename("SiRCle cluster Name"= "RG1_All")
  ClusterSummary_RG1$`Regulation Grouping` <- "RG1_All"

  ClusterSummary_RG2 <- MergeDF_Rearrange[,c("Metabolite", "RG2_Significant")]%>%
    count(RG2_Significant, name="Number of Features")%>%
    dplyr::rename("SiRCle cluster Name"= "RG2_Significant")
  ClusterSummary_RG2$`Regulation Grouping` <- "RG2_Significant"

  ClusterSummary_RG3 <- MergeDF_Rearrange[,c("Metabolite", "RG3_SignificantChange")]%>%
    count(RG3_SignificantChange, name="Number of Features")%>%
    dplyr::rename("SiRCle cluster Name"= "RG3_SignificantChange")
  ClusterSummary_RG3$`Regulation Grouping` <- "RG3_SignificantChange"

  ClusterSummary <- rbind(ClusterSummary_RG1, ClusterSummary_RG2,ClusterSummary_RG3)
  ClusterSummary <- ClusterSummary[,c(3,1,2)]

  ## Rename FeatureID
  MergeDF_Rearrange <-MergeDF_Rearrange%>%
    dplyr::rename(!!FeatureID := "Metabolite")

  ######################################################################################################################################################################
  ##----- Save and Return
  #Here we make a list in which we will save the outputs:
  DF_List <- list("MCA_2Cond_Summary"=ClusterSummary, "MCA_2Cond_Results"=MergeDF_Rearrange)

  suppressMessages(suppressWarnings(
    SaveRes(InputList_DF=DF_List,
                         InputList_Plot= NULL,
                         SaveAs_Table=SaveAs_Table,
                         SaveAs_Plot=NULL,
                         FolderPath= Folder,
                         FileName= "MCA_2Cond",
                         CoRe=FALSE,
                         PrintPlot=FALSE)))

  #Return:
  invisible(return(DF_List))
}




#' MCA_CoRe
#'
#' This script allows you to perform metabolite clustering analysis and computes clusters of metabolites based on regulatory rules between Intracellular and culture media metabolomics (CoRe experiment).
#'
#' @param InputData_Intra DF for your data (results from e.g. DMA) containing metabolites in rows with corresponding Log2FC and stat (p-value, p.adjusted) value columns.
#' @param InputData_CoRe DF for your data (results from e.g. DMA) containing metabolites in rows with corresponding Log2FC and stat (p-value, p.adjusted) value columns. Here we additionally require
#' @param SettingsInfo_Intra  \emph{Optional: } Pass ColumnNames and Cutoffs for the intracellular metabolomics including the value column (e.g. Log2FC, Log2Diff, t.val, etc) and the stats column (e.g. p.adj, p.val). This must include: c(ValueCol=ColumnName_InputData_Intra,StatCol=ColumnName_InputData_Intra, StatCutoff= NumericValue, ValueCutoff=NumericValue) \strong{Default=c(ValueCol="Log2FC",StatCol="p.adj", StatCutoff= 0.05, ValueCutoff=1)}
#' @param SettingsInfo_CoRe  \emph{Optional: } Pass ColumnNames and Cutoffs for the consumption-release metabolomics including the direction column, the value column (e.g. Log2Diff, t.val, etc) and the stats column (e.g. p.adj, p.val). This must include: c(DirectionCol= ColumnName_InputData_CoRe,ValueCol=ColumnName_InputData_CoRe,StatCol=ColumnName_InputData_CoRe, StatCutoff= NumericValue, ValueCutoff=NumericValue)\strong{Default=c(DirectionCol="CoRe", ValueCol="Log2(Distance)",StatCol="p.adj", StatCutoff= 0.05, ValueCutoff=1)}
#' @param FeatureID \emph{Optional: } Column name of Column including the Metabolite identifiers. This MUST BE THE SAME in each of your Input files. \strong{Default="Metabolite"}
#' @param BackgroundMethod \emph{Optional: } Background method `Intra|CoRe, Intra&CoRe, CoRe, Intra or * \strong{Default="Intra&CoRe"}
#' @param SaveAs_Table \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt" \strong{default: "csv"}
#' @param FolderPath \emph{Optional:} Path to the folder the results should be saved at. \strong{default: NULL}
#' @return MCA an instance of the MetaProViz package
#' @export
#'

##################################################
### ### ### Metabolite Clustering Analysis ### ### ###
##################################################

MCA_CoRe <- function(InputData_Intra,
                     InputData_CoRe,
                     SettingsInfo_Intra=c(ValueCol="Log2FC",StatCol="p.adj", StatCutoff= 0.05, ValueCutoff=1),
                     SettingsInfo_CoRe=c(DirectionCol="CoRe", ValueCol="Log2(Distance)",StatCol="p.adj", StatCutoff= 0.05, ValueCutoff=1),
                     FeatureID= "Metabolite",
                     SaveAs_Table = "csv",
                     BackgroundMethod="Intra&CoRe",
                     FolderPath=NULL
                     ){

  ## ------------ Create log file ----------- ##
  MetaProViz_Init()

  ################################################################################################################################################################################################
  ## ------------ Check Input files ----------- ##
  CheckInput_MCA(InputData_C1=NULL,
                              InputData_C2=NULL,
                              InputData_CoRe= InputData_CoRe,
                              InputData_Intra=InputData_Intra,
                              SettingsInfo_C1=NULL,
                              SettingsInfo_C2=NULL,
                              SettingsInfo_CoRe=SettingsInfo_CoRe,
                              SettingsInfo_Intra=SettingsInfo_Intra,
                              BackgroundMethod=BackgroundMethod,
                              FeatureID=FeatureID,
                              SaveAs_Table=SaveAs_Table)


  ## ------------ Create Results output folder ----------- ##
  if(is.null(SaveAs_Table)==FALSE){
    Folder <- SavePath(FolderName= "MCACoRe",
                                    FolderPath=FolderPath)
  }


  ################################################################################################################################################################################################
  ## ------------ Prepare the Input -------- ##
  #Import the data and check columns (here the user will get an error if the column can not be renamed as it does not exists.)
  CoRe_DF <-  as.data.frame(InputData_CoRe)%>%
    dplyr::rename("MetaboliteID"=paste(FeatureID),
                  "ValueCol"=SettingsInfo_CoRe[["ValueCol"]],
                  "PadjCol"=SettingsInfo_CoRe[["StatCol"]],
                  "CoRe_Direction"=SettingsInfo_CoRe[["DirectionCol"]])
  CoRe_DF <- CoRe_DF[complete.cases(CoRe_DF$ValueCol, CoRe_DF$PadjCol), ]

  Intra_DF<- as.data.frame(InputData_Intra)%>%
    dplyr::rename("MetaboliteID"=paste(FeatureID),
                  "ValueCol"=SettingsInfo_Intra[["ValueCol"]],
                  "PadjCol"=SettingsInfo_Intra[["StatCol"]])
  Intra_DF <- Intra_DF[complete.cases(Intra_DF$ValueCol, Intra_DF$PadjCol), ]

  #Tag genes that are detected in each data layer
  CoRe_DF$Detected <- "TRUE"
  Intra_DF$Detected <- "TRUE"

  #Assign to Group based on individual Cutoff ("UP", "DOWN", "No Change")
  CoRe_DF <- CoRe_DF%>%
    mutate(Cutoff = case_when(CoRe_DF$PadjCol <= as.numeric(SettingsInfo_CoRe[["StatCutoff"]]) & CoRe_DF$ValueCol > as.numeric(SettingsInfo_CoRe[["ValueCutoff"]]) ~ 'UP',
                              CoRe_DF$PadjCol <= as.numeric(SettingsInfo_CoRe[["StatCutoff"]]) & CoRe_DF$ValueCol < - as.numeric(SettingsInfo_CoRe[["ValueCutoff"]]) ~ 'DOWN',
                              TRUE ~ 'No Change')) %>%
    mutate(Cutoff_Specific = case_when(Cutoff == "UP" ~ 'UP',
                                       Cutoff == "DOWN" ~ 'DOWN',
                                       Cutoff == "No Change" & CoRe_DF$PadjCol <= as.numeric(SettingsInfo_CoRe[["StatCutoff"]]) & CoRe_DF$ValueCol > 0 ~ 'Significant Positive',
                                       Cutoff == "No Change" & CoRe_DF$PadjCol <= as.numeric(SettingsInfo_CoRe[["StatCutoff"]]) & CoRe_DF$ValueCol < 0 ~ 'Significant Negative',
                                       Cutoff == "No Change" & CoRe_DF$PadjCol > as.numeric(SettingsInfo_CoRe[["StatCutoff"]]) ~ 'Not Significant',
                                       TRUE ~ 'NA'))

  Intra_DF <-Intra_DF%>%
    mutate(Cutoff = case_when(Intra_DF$PadjCol <= as.numeric(SettingsInfo_Intra[["StatCutoff"]]) & Intra_DF$ValueCol > as.numeric(SettingsInfo_Intra[["ValueCutoff"]]) ~ 'UP',
                              Intra_DF$PadjCol <= as.numeric(SettingsInfo_Intra[["StatCutoff"]]) & Intra_DF$ValueCol < - as.numeric(SettingsInfo_Intra[["ValueCutoff"]]) ~ 'DOWN',
                              TRUE ~ 'No Change'))%>%
    mutate(Cutoff_Specific = case_when(Cutoff == "UP" ~ 'UP',
                                       Cutoff == "DOWN" ~ 'DOWN',
                                       Cutoff == "No Change" & Intra_DF$PadjCol <= as.numeric(SettingsInfo_Intra[["StatCutoff"]]) & Intra_DF$ValueCol > 0 ~ 'Significant Positive',
                                       Cutoff == "No Change" & Intra_DF$PadjCol <= as.numeric(SettingsInfo_Intra[["StatCutoff"]]) & Intra_DF$ValueCol < 0 ~ 'Significant Negative',
                                       Cutoff == "No Change" & Intra_DF$PadjCol > as.numeric(SettingsInfo_Intra[["StatCutoff"]]) ~ 'Not Significant',
                                       TRUE ~ 'NA'))

  #Merge the dataframes together: Merge the supplied Intra and CoRe dataframes together.
  ##Add prefix to column names to distinguish the different data types after merge
  colnames(CoRe_DF) <- paste0("CoRe_DF_", colnames(CoRe_DF))
  CoRe_DF <- CoRe_DF%>%
    dplyr::rename("MetaboliteID" = "CoRe_DF_MetaboliteID")

  colnames(Intra_DF) <- paste0("Intra_DF_", colnames(Intra_DF))
  Intra_DF <-Intra_DF%>%
    dplyr::rename("MetaboliteID"="Intra_DF_MetaboliteID")

  ##Merge
  MergeDF <- merge(Intra_DF, CoRe_DF, by="MetaboliteID", all=TRUE)

  ##Mark the undetected genes in each data layer
  MergeDF<-MergeDF %>%
    mutate_at(c("CoRe_DF_Detected","Intra_DF_Detected"), ~replace_na(.,"FALSE"))%>%
    mutate_at(c("CoRe_DF_Cutoff","Intra_DF_Cutoff"), ~replace_na(.,"No Change"))%>%
    mutate_at(c("CoRe_DF_Cutoff_Specific", "Intra_DF_Cutoff_Specific"), ~replace_na(.,"Not Detected"))%>%
    mutate_at(c("CoRe_DF_CoRe_Direction"), ~replace_na(.,"Not Detected"))

  #Apply Background filter (label metabolites that will be removed based on chosen background)
  if(BackgroundMethod == "Intra|CoRe"){# C1|C2 = CoRe OR Intra
    MergeDF <- MergeDF%>%
      mutate(BG_Method = case_when(Intra_DF_Detected=="TRUE" & CoRe_DF_Detected=="TRUE" ~ 'TRUE', #Intra & CoRe
                                   Intra_DF_Detected=="TRUE" & CoRe_DF_Detected=="FALSE" ~ 'TRUE', # JustIntra
                                   Intra_DF_Detected=="FALSE" & CoRe_DF_Detected=="TRUE" ~ 'TRUE', # Just CoRe
                                   TRUE ~ 'FALSE'))
  }else if(BackgroundMethod == "Intra&CoRe"){ # CoRe AND Intra
    MergeDF <- MergeDF%>%
      mutate(BG_Method = case_when(Intra_DF_Detected=="TRUE" & CoRe_DF_Detected=="TRUE" ~ 'TRUE', #Intra & CoRe
                                   TRUE ~ 'FALSE'))
  }else if(BackgroundMethod == "CoRe"){ # CoRe has to be there
    MergeDF <- MergeDF%>%
      mutate(BG_Method = case_when(Intra_DF_Detected=="TRUE" & CoRe_DF_Detected=="TRUE" ~ 'TRUE', #Intra & CoRe
                                   Intra_DF_Detected=="FALSE" & CoRe_DF_Detected=="TRUE" ~ 'TRUE', # Just CoRe
                                   TRUE ~ 'FALSE'))
  }else if(BackgroundMethod == "Intra"){ #Intra has to be there
    MergeDF <- MergeDF%>%
      mutate(BG_Method = case_when(Intra_DF_Detected=="TRUE" & CoRe_DF_Detected=="TRUE" ~ 'TRUE', #Intra & CoRe
                                   Intra_DF_Detected=="TRUE" & CoRe_DF_Detected=="FALSE" ~ 'TRUE', # JustIntra
                                   TRUE ~ 'FALSE'))
  }else if(BackgroundMethod == "*"){ # Use all metabolites as the background
    MergeDF$BG_Method <- "TRUE"
  }else{
    stop("Please use one of the following BackgroundMethods: Intra|CoRe, Intra&CoRe, CoRe, Intra, *")#error message
  }

  #Assign Metabolite cluster names to the metabolites
  MergeDF <- MergeDF%>%
        mutate(RG1_All = case_when(BG_Method =="FALSE"~ 'Background = FALSE',
                               Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Released" ~ 'Intra DOWN + CoRe DOWN_Released',
                               Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected"~ 'Intra DOWN + CoRe Not Detected',
                               Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Released"~ 'Intra DOWN + CoRe Not Significant_Released',
                               Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Released"~ 'Intra DOWN + CoRe Significant Negative_Released',
                               Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Released"~ 'Intra DOWN + CoRe Significant Positive_Released',
                               Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Released" ~ 'Intra DOWN + CoRe UP_Released',

                               Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Released" ~ 'Intra UP + CoRe DOWN_Released',
                               Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected" ~ 'Intra UP + CoRe Not Detected',
                               Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Released" ~ 'Intra UP + CoRe Not Significant_Released',
                               Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Released" ~ 'Intra UP + CoRe Significant Negative_Released',
                               Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Released" ~ 'Intra UP + CoRe Significant Positive_Released',
                               Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Released" ~ 'Intra UP + CoRe UP_Released',

                               Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Released" ~ 'Intra Not Detected + CoRe DOWN_Released',
                               Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected" ~ 'Intra Not Detected + CoRe Not Detected',
                               Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Released" ~ 'Intra Not Detected + CoRe Not Significant_Released',
                               Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Released" ~ 'Intra Not Detected + CoRe Significant Negative_Released',
                               Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Released" ~ 'Intra Not Detected + CoRe Significant Positive_Released',
                               Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Released" ~ 'Intra Not Detected + CoRe UP_Released',

                               Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Released"~ 'Intra Significant Negative + CoRe DOWN_Released',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected"~ 'Intra Significant Negative + CoRe Not Detected',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Released"~ 'Intra Significant Negative + CoRe Not Significant_Released',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Released"~ 'Intra Significant Negative + CoRe Significant Negative_Released',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Released"~ 'Intra Significant Negative + CoRe Significant Positive_Released',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Released"~ 'Intra Significant Negative + CoRe UP_Released',

                               Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Released"~ 'Intra Significant Positive + CoRe DOWN_Released',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected"~ 'Intra Significant Positive + CoRe Not Detected',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Released"~ 'Intra Significant Positive + CoRe Not Significant_Released',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Released"~ 'Intra Significant Positive + CoRe Significant Negative_Released',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Released"~ 'Intra Significant Positive + CoRe Significant Positive_Released',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Released"~ 'Intra Significant Positive + CoRe UP_Released',

                               Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Released"~ 'Intra Not Significant + CoRe DOWN_Released',
                               Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected"~ 'Intra Not Significant + CoRe Not Detected',
                               Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Released"~ 'Intra Not Significant + CoRe Not Significant_Released',
                               Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Released"~ 'Intra Not Significant + CoRe Significant Negative_Released',
                               Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Released"~ 'Intra Not Significant + CoRe Significant Positive_Released',
                               Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Released"~ 'Intra Not Significant + CoRe UP_Released',

                               #Consumed:
                               Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Consumed" ~ 'Intra DOWN + CoRe DOWN_Consumed',
                               Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected"~ 'Intra DOWN + CoRe Not Detected',
                               Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Consumed"~ 'Intra DOWN + CoRe Not Significant_Consumed',
                               Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Consumed"~ 'Intra DOWN + CoRe Significant Negative_Consumed',
                               Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Consumed"~ 'Intra DOWN + CoRe Significant Positive_Consumed',
                               Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Consumed" ~ 'Intra DOWN + CoRe UP_Consumed',

                               Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Consumed" ~ 'Intra UP + CoRe DOWN_Consumed',
                               Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected" ~ 'Intra UP + CoRe Not Detected',
                               Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Consumed" ~ 'Intra UP + CoRe Not Significant_Consumed',
                               Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Consumed" ~ 'Intra UP + CoRe Significant Negative_Consumed',
                               Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Consumed" ~ 'Intra UP + CoRe Significant Positive_Consumed',
                               Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Consumed" ~ 'Intra UP + CoRe UP_Consumed',

                               Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Consumed" ~ 'Intra Not Detected + CoRe DOWN_Consumed',
                               Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected" ~ 'Intra Not Detected + CoRe Not Detected',
                               Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Consumed" ~ 'Intra Not Detected + CoRe Not Significant_Consumed',
                               Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Consumed" ~ 'Intra Not Detected + CoRe Significant Negative_Consumed',
                               Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Consumed" ~ 'Intra Not Detected + CoRe Significant Positive_Consumed',
                               Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Consumed" ~ 'Intra Not Detected + CoRe UP_Consumed',

                               Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Consumed"~ 'Intra Significant Negative + CoRe DOWN_Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected"~ 'Intra Significant Negative + CoRe Not Detected',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Consumed"~ 'Intra Significant Negative + CoRe Not Significant_Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Consumed"~ 'Intra Significant Negative + CoRe Significant Negative_Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Consumed"~ 'Intra Significant Negative + CoRe Significant Positive_Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Consumed"~ 'Intra Significant Negative + CoRe UP_Consumed',

                               Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Consumed"~ 'Intra Significant Positive + CoRe DOWN_Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected"~ 'Intra Significant Positive + CoRe Not Detected',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Consumed"~ 'Intra Significant Positive + CoRe Not Significant_Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Consumed"~ 'Intra Significant Positive + CoRe Significant Negative_Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Consumed"~ 'Intra Significant Positive + CoRe Significant Positive_Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Consumed"~ 'Intra Significant Positive + CoRe UP_Consumed',

                               Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Consumed"~ 'Intra Not Significant + CoRe DOWN_Consumed',
                               Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected"~ 'Intra Not Significant + CoRe Not Detected',
                               Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Consumed"~ 'Intra Not Significant + CoRe Not Significant_Consumed',
                               Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Consumed"~ 'Intra Not Significant + CoRe Significant Negative_Consumed',
                               Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Consumed"~ 'Intra Not Significant + CoRe Significant Positive_Consumed',
                               Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Consumed"~ 'Intra Not Significant + CoRe UP_Consumed',

                               #Released/Consumed (Consumed in one, released in the other)
                               Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'Intra DOWN + CoRe DOWN_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected"~ 'Intra DOWN + CoRe Not Detected',
                               Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'Intra DOWN + CoRe Not Significant_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'Intra DOWN + CoRe Significant Negative_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'Intra DOWN + CoRe Significant Positive_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'Intra DOWN + CoRe UP_Released/Consumed',

                               Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'Intra UP + CoRe DOWN_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected" ~ 'Intra UP + CoRe Not Detected',
                               Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'Intra UP + CoRe Not Significant_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'Intra UP + CoRe Significant Negative_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'Intra UP + CoRe Significant Positive_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'Intra UP + CoRe UP_Released/Consumed',

                               Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'Intra Not Detected + CoRe DOWN_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected" ~ 'Intra Not Detected + CoRe Not Detected',
                               Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'Intra Not Detected + CoRe Not Significant_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'Intra Not Detected + CoRe Significant Negative_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'Intra Not Detected + CoRe Significant Positive_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'Intra Not Detected + CoRe UP_Released/Consumed',

                               Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'Intra Significant Negative + CoRe DOWN_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected"~ 'Intra Significant Negative + CoRe Not Detected',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'Intra Significant Negative + CoRe Not Significant_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'Intra Significant Negative + CoRe Significant Negative_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'Intra Significant Negative + CoRe Significant Positive_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'Intra Significant Negative + CoRe UP_Released/Consumed',

                               Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'Intra Significant Positive + CoRe DOWN_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected"~ 'Intra Significant Positive + CoRe Not Detected',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'Intra Significant Positive + CoRe Not Significant_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'Intra Significant Positive + CoRe Significant Negative_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'Intra Significant Positive + CoRe Significant Positive_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'Intra Significant Positive + CoRe UP_Released/Consumed',

                               Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'Intra Not Significant + CoRe DOWN_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected"~ 'Intra Not Significant + CoRe Not Detected',
                               Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'Intra Not Significant + CoRe Not Significant_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'Intra Not Significant + CoRe Significant Negative_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'Intra Not Significant + CoRe Significant Positive_Released/Consumed',
                               Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'Intra Not Significant + CoRe UP_Released/Consumed',
                               TRUE ~ 'NA'))%>%
    mutate(RG2_Significant = case_when(BG_Method =="FALSE"~ 'Background = FALSE',
                                       Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Released" ~ 'Both_DOWN (Released)',
                                       Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected"~ 'None',
                                       Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Released"~ 'None',
                                       Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Released"~ 'Both_DOWN (Released)',
                                       Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Released"~ 'Opposite (Released UP)',
                                       Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Released" ~ 'Opposite (Released UP)',

                                       Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Released" ~ 'Opposite (Released DOWN)',
                                       Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Released" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Released" ~ 'Opposite (Released UP)',
                                       Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Released" ~ 'Both_UP (Released)',
                                       Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Released" ~ 'Both_UP (Released)',

                                       Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Released" ~ 'CoRe_DOWN (Released)',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Released" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Released" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Released" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Released" ~ 'CoRe_UP (Released)',

                                       Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Released"~ 'Both_DOWN (Released)',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Released"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Released"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Released"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Released"~ 'Opposite (Released UP)',

                                       Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Released"~ 'Opposite (Released DOWN)',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Released"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Released"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Released"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Released"~ 'Both_UP (Released)',

                                       Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Released"~ 'CoRe_DOWN (Released)',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Released"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Released"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Released"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Released"~ 'CoRe_UP (Released)',

                                       #Consumed:
                                       Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Consumed" ~ 'Both_DOWN (Consumed)',
                                       Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected"~ 'None',
                                       Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Consumed"~ 'Both_DOWN (Consumed)',
                                       Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Consumed"~ 'Opposite (Consumed UP)',
                                       Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Consumed" ~ 'Opposite (Consumed UP)',

                                       Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Consumed" ~ 'Opposite (Consumed DOWN)',
                                       Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Consumed" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Consumed" ~ 'Opposite (Consumed UP)',
                                       Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Consumed" ~ 'Both_UP (Consumed)',
                                       Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Consumed" ~ 'Both_UP (Consumed)',

                                       Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Consumed" ~ 'CoRe_DOWN (Consumed)',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Consumed" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Consumed" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Consumed" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Consumed" ~ 'CoRe_UP (Consumed)',

                                       Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Consumed"~ 'Both_DOWN (Consumed)',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Consumed"~ 'Opposite (Consumed UP)',

                                       Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Consumed"~ 'Opposite (Consumed DOWN)',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Consumed"~ 'Both_UP (Consumed)',

                                       Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Consumed"~ 'CoRe_DOWN (Consumed)',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Consumed"~ 'CoRe_UP (Consumed)',

                                       #Consumed/Released (Consumed in one, released in the other)
                                       Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Released/Consumed" ~ 'Both_DOWN (Released/Consumed)',
                                       Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected"~ 'None',
                                       Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'Both_DOWN (Released/Consumed)',
                                       Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'Opposite (Released/Consumed UP)',
                                       Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Released/Consumed" ~ 'Opposite (Released/Consumed UP)',

                                       Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Released/Consumed" ~ 'Opposite (Released/Consumed DOWN)',
                                       Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Released/Consumed" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Released/Consumed" ~ 'Opposite (Released/Consumed UP)',
                                       Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Released/Consumed" ~ 'Both_UP (Released/Consumed)',
                                       Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Released/Consumed" ~ 'Both_UP (Released/Consumed)',

                                       Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Released/Consumed" ~ 'CoRe_DOWN (Released/Consumed)',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Released/Consumed" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Released/Consumed" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Released/Consumed" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Released/Consumed" ~ 'CoRe_UP (Released/Consumed)',

                                       Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'Both_DOWN (Released/Consumed)',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'Opposite (Released/Consumed UP)',

                                       Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'Opposite (Released/Consumed DOWN)',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'Both_UP (Released/Consumed)',

                                       Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'CoRe_DOWN (Released/Consumed)',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'CoRe_UP (Released/Consumed)',
                                       TRUE ~ 'NA'))%>%
    mutate(RG3_Change = case_when(BG_Method =="FALSE"~ 'Background = FALSE',
                                  Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Released" ~ 'Both_DOWN (Released)',
                                  Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected"~ 'None',
                                  Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Released"~ 'None',
                                  Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Released"~ 'None',
                                  Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Released"~ 'None',
                                  Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Released" ~ 'Opposite (Released UP)',

                                  Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Released" ~ 'Opposite (Released DOWN)',
                                  Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Released" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Released" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Released" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Released" ~ 'Both_UP (Released)',

                                  Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Released" ~ 'CoRe_DOWN (Released)',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Released" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Released" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Released" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Released" ~ 'CoRe_UP (Released)',

                                  Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Released"~ 'CoRe_DOWN (Released)',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Released"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Released"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Released"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Released"~ 'CoRe_UP (Released)',

                                  Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Released"~ 'CoRe_DOWN (Released)',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Released"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Released"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Released"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Released"~ 'CoRe_UP (Released)',

                                  Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Released"~ 'CoRe_DOWN (Released)',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Released"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Released"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Released"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Released"~ 'CoRe_UP (Released)',

                                  #Consumed:
                                  Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Consumed" ~ 'Both_DOWN (Consumed)',
                                  Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected"~ 'None',
                                  Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Consumed" ~ 'Opposite (Consumed UP)',

                                  Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Consumed" ~ 'Opposite (Consumed DOWN)',
                                  Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Consumed" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Consumed" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Consumed" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Consumed" ~ 'Both_UP (Consumed)',

                                  Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Consumed" ~ 'CoRe_DOWN (Consumed)',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Consumed" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Consumed" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Consumed" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Consumed" ~ 'CoRe_UP (Consumed)',

                                  Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Consumed"~ 'CoRe_DOWN (Consumed)',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Consumed"~ 'CoRe_UP (Consumed)',

                                  Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Consumed"~ 'CoRe_DOWN (Consumed)',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Consumed"~ 'CoRe_UP (Consumed)',

                                  Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Consumed"~ 'CoRe_DOWN (Consumed)',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Consumed"~ 'CoRe_UP (Consumed)',

                                  #Consumed/Released (Consumed in one, released in the other)
                                  Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Released/Consumed" ~ 'Both_DOWN (Released/Consumed)',
                                  Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected"~ 'None',
                                  Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Released/Consumed" ~ 'Opposite (Released/Consumed UP)',

                                  Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Released/Consumed" ~ 'Opposite (Released/Consumed DOWN)',
                                  Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Released/Consumed" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Released/Consumed" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Released/Consumed" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Released/Consumed" ~ 'Both_UP (Released/Consumed)',

                                  Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Released/Consumed" ~ 'CoRe_DOWN (Released/Consumed)',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Released/Consumed" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Released/Consumed" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Released/Consumed" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Released/Consumed" ~ 'CoRe_UP (Released/Consumed)',

                                  Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'CoRe_DOWN (Released/Consumed)',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'CoRe_UP (Released/Consumed)',

                                  Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'CoRe_DOWN (Released/Consumed)',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'CoRe_UP (Released/Consumed)',

                                  Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'CoRe_DOWN (Released/Consumed)',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_CoRe_Direction=="Not Detected"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_DF_CoRe_Direction=="Released/Consumed"~ 'CoRe_UP (Released/Consumed)',
                                  TRUE ~ 'NA'))

  #Safe the DF and return the groupings
  ##MCA DF (Merged InputDF filtered for background with assigned MCA cluster names)
  MergeDF_Select1 <- MergeDF[, c("MetaboliteID", "Intra_DF_Detected","Intra_DF_ValueCol","Intra_DF_PadjCol","Intra_DF_Cutoff", "Intra_DF_Cutoff_Specific", "CoRe_DF_Detected", "CoRe_DF_ValueCol","CoRe_DF_PadjCol","CoRe_DF_Cutoff", "CoRe_DF_Cutoff_Specific", "BG_Method", "RG1_All", "RG2_Significant", "RG3_Change")]

  CoReValueCol_Unique<-paste("Cond2_DF_",SettingsInfo_CoRe[["ValueCol"]])
  CoRePadjCol_Unique <-paste("Cond2_DF_",SettingsInfo_CoRe[["StatCol"]])
  IntraValueCol_Unique<-paste("Cond1_DF_",SettingsInfo_Intra[["ValueCol"]])
  IntraPadjCol_Unique <-paste("Cond1_DF_",SettingsInfo_Intra[["StatCol"]])

  MergeDF_Select2<- subset(MergeDF, select=-c(Intra_DF_Detected,Intra_DF_Cutoff, CoRe_DF_Detected,CoRe_DF_Cutoff, CoRe_DF_Cutoff_Specific, BG_Method, RG1_All, RG2_Significant, RG3_Change))%>%
    dplyr::rename(!!CoReValueCol_Unique :="CoRe_DF_ValueCol",#This syntax is needed since paste(MetaboliteID)="MetaboliteID" is not working in dyplr
                  !!CoRePadjCol_Unique :="CoRe_DF_PadjCol",
                  !!IntraValueCol_Unique :="Intra_DF_ValueCol",
                  !!IntraPadjCol_Unique :="Intra_DF_PadjCol")

  MergeDF_Rearrange <- merge(MergeDF_Select1, MergeDF_Select2, by="MetaboliteID")
  MergeDF_Rearrange <-MergeDF_Rearrange%>%
    dplyr::rename("Metabolite"="MetaboliteID")


  ##Summary SiRCle clusters (number of genes assigned to each SiRCle cluster in each grouping)
  ClusterSummary_RG1 <- MergeDF_Rearrange[,c("Metabolite", "RG1_All")]%>%
    group_by(RG1_All) %>%
    mutate("Number of Features" = n()) %>%
    distinct(RG1_All, .keep_all = TRUE) %>%
    dplyr::rename("SiRCle cluster Name" = "RG1_All")
  ClusterSummary_RG1$`Regulation Grouping` <- "RG1_All"
  ClusterSummary_RG1 <- ClusterSummary_RG1[-c(1)]

  ClusterSummary_RG2 <- MergeDF_Rearrange[,c("Metabolite", "RG2_Significant")]%>%
    group_by(RG2_Significant) %>%
    mutate("Number of Features" = n()) %>%
    distinct(RG2_Significant, .keep_all = TRUE) %>%
    dplyr::rename("SiRCle cluster Name"= "RG2_Significant")
  ClusterSummary_RG2$`Regulation Grouping` <- "RG2_Significant"
  ClusterSummary_RG2 <- ClusterSummary_RG2[-c(1)]

  ClusterSummary_RG3 <- MergeDF_Rearrange[,c("Metabolite", "RG3_Change")]%>%
    group_by(RG3_Change) %>%
    mutate("Number of Features" = n()) %>%
    distinct(RG3_Change, .keep_all = TRUE) %>%
    dplyr::rename("SiRCle cluster Name"= "RG3_Change")
  ClusterSummary_RG3$`Regulation Grouping` <- "RG3_Change"
  ClusterSummary_RG3 <- ClusterSummary_RG3[-c(1)]

  ClusterSummary <- rbind(ClusterSummary_RG1, ClusterSummary_RG2,ClusterSummary_RG3)
  ClusterSummary <- ClusterSummary[,c(3,1,2)]

  ## Rename FeatureID
  MergeDF_Rearrange <-MergeDF_Rearrange%>%
    dplyr::rename(!!FeatureID := "Metabolite")

  ######################################################################################################################################################################
  ##----- Save and Return
  #Here we make a list in which we will save the outputs:
  DF_List <- list("MCA_CoRe_Summary"=ClusterSummary, "MCA_CoRe_Results"=MergeDF_Rearrange)

  suppressMessages(suppressWarnings(
    SaveRes(InputList_DF=DF_List,
                         InputList_Plot= NULL,
                         SaveAs_Table=SaveAs_Table,
                         SaveAs_Plot=NULL,
                         FolderPath= Folder,
                         FileName= "MCA_2Cond",
                         CoRe=FALSE,
                         PrintPlot=FALSE)))

  invisible(return(DF_List))
}

########################################
### ### ### Import MCA rules ### ### ###
########################################

#' Imports MCA regulatory rules into environment
#'
#' @param Method Either "2Cond" or "CoRe" depending which regulatory rules you would like to load
#' @title MCA regulatory rules Import
#' @description Import and process .csv file to create toy data.
#' @importFrom utils read.csv
#' @return A data frame containing the toy data.
#' @export
#'
MCA_rules<- function(Method){
  ## ------------ Create log file ----------- ##
  MetaProViz_Init()

  # Read the .csv files
  Cond <- system.file("data", "MCA_2Cond.csv", package = "MetaProViz")
  Cond<- read.csv( Cond, check.names=FALSE)

  CoRe <- system.file("data", "MCA_CoRe.csv", package = "MetaProViz")
  CoRe<- read.csv(CoRe, check.names=FALSE)

  # Return the toy data into environment
  if(Method=="2Cond"){
    assign("MCA_2Cond", Cond, envir=.GlobalEnv)
  } else if(Method=="CoRe"){
    assign("MCA_CoRe", CoRe, envir=.GlobalEnv)
  } else{
    warning("Please choose the MCA regulatory rules you would like to load: 2Cond, CoRe")
  }
}


################################################################################################
### ### ### MCA Helper function: Internal Function to check function input ### ### ###
################################################################################################

#' Check input parameters
#'
#' @param InputData_C1 Passed to main function MCA. If not avaliable can be set to NULL.
#' @param InputData_C2 Passed to main function MCA. If not avaliable can be set to NULL.
#' @param InputData_Intra Passed to main function MCA. If not avaliable can be set to NULL.
#' @param InputData_CoRe Passed to main function MCA. If not avaliable can be set to NULL.
#' @param SettingsInfo_C1 Passed to main function MCA. If not avaliable can be set to NULL.
#' @param SettingsInfo_C2 Passed to main function MCA. If not avaliable can be set to NULL.
#' @param SettingsInfo_Intra Passed to main function MCA. If not avaliable can be set to NULL.
#' @param SettingsInfo_CoRe Passed to main function MCA. If not avaliable can be set to NULL.
#' @param BackgroundMethod Passed to main function MCA.
#' @param FeatureID Passed to main function MCA.
#' @param SaveAs_Table Passed to main function PreProcessing(). If not avaliable can be set to NULL.
#'
#' @param Function Name of the MetaProViz Function that is checked.
#' @param InputList
#'
#'
#' @keywords Input check
#' @noRd
#'
#'

CheckInput_MCA <- function(InputData_C1,
                           InputData_C2,
                           InputData_CoRe,
                           InputData_Intra,
                           SettingsInfo_C1,
                           SettingsInfo_C2,
                           SettingsInfo_CoRe,
                           SettingsInfo_Intra,
                           BackgroundMethod,
                           FeatureID,
                           SaveAs_Table
                           ){
  ## ------------ Create log file ----------- ##
  MetaProViz_Init()

  #------------- InputData
  if(is.null(InputData_C1)==FALSE){
    if(class(InputData_C1) != "data.frame"| class(InputData_C2) != "data.frame"){
    stop("InputData_C1 and InputData_C2 should be a data.frame. It's currently a ", paste(class(InputData_C1)), paste(class(InputData_C2)), ".",sep = "")
    }
    if(length(InputData_C1[duplicated(InputData_C1[[FeatureID]]), FeatureID]) > 0){
      stop("Duplicated FeatureIDs of InputData_C1, whilst features must be unique")
    }
    if(length(InputData_C2[duplicated(InputData_C2[[FeatureID]]), FeatureID]) > 0){
      stop("Duplicated FeatureIDs of InputData_C2, whilst features must be unique")
    }

  }else{
    if(class(InputData_Intra) != "data.frame"| class(InputData_CoRe) != "data.frame"){
      stop("InputData_Intra and InputData_CoRe should be a data.frame. It's currently a ", paste(class(InputData_Intra)), paste(class(InputData_CoRe)), ".",sep = "")
    }
    if(length(InputData_Intra[duplicated(InputData_Intra[[FeatureID]]), FeatureID]) > 0){
      stop("Duplicated FeatureIDs of InputData_Intra, whilst features must be unique")
    }
    if(length(InputData_CoRe[duplicated(InputData_CoRe[[FeatureID]]), FeatureID]) > 0){
      stop("Duplicated FeatureIDs of InputData_CoRe, whilst features must be unique")
    }
  }


  #------------- SettingsInfo
  if(is.null(SettingsInfo_C1)==FALSE){
    ## C1
    #ValueCol
    if("ValueCol" %in% names(SettingsInfo_C1)){
      if(SettingsInfo_C1[["ValueCol"]] %in% colnames(InputData_C1)== FALSE){
        stop("The ", SettingsInfo_C1[["ValueCol"]], " column selected as ValueCol in SettingsInfo_C1 was not found in InputData_C1. Please check your input.")
      }
    }
    #StatCol
    if("StatCol" %in% names(SettingsInfo_C1)){
      if(SettingsInfo_C1[["StatCol"]] %in% colnames(InputData_C1)== FALSE){
        stop("The ", SettingsInfo_C1[["StatCol"]], " column selected as StatCol in SettingsInfo_C1 was not found in InputData_C1. Please check your input.")
      }
    }

    ## C2
    #ValueCol
    if("ValueCol" %in% names(SettingsInfo_C2)){
      if(SettingsInfo_C2[["ValueCol"]] %in% colnames(InputData_C2)== FALSE){
        stop("The ", SettingsInfo_C2[["ValueCol"]], " column selected as ValueCol in SettingsInfo_C2 was not found in InputData_C2. Please check your input.")
      }
    }
    #StatCol
    if("StatCol" %in% names(SettingsInfo_C2)){
      if(SettingsInfo_C2[["StatCol"]] %in% colnames(InputData_C2)== FALSE){
        stop("The ", SettingsInfo_C2[["StatCol"]], " column selected as StatCol in SettingsInfo_C2 was not found in InputData_C2. Please check your input.")
      }
    }
  }else{
    ## Intra
    #ValueCol
    if("ValueCol" %in% names(SettingsInfo_Intra)){
      if(SettingsInfo_Intra[["ValueCol"]] %in% colnames(InputData_Intra)== FALSE){
        stop("The ", SettingsInfo_Intra[["ValueCol"]], " column selected as ValueCol in SettingsInfo_Intra was not found in InputData_Intra. Please check your input.")
      }
    }
    #StatCol
    if("StatCol" %in% names(SettingsInfo_Intra)){
      if(SettingsInfo_Intra[["StatCol"]] %in% colnames(InputData_Intra)== FALSE){
        stop("The ", SettingsInfo_Intra[["StatCol"]], " column selected as StatCol in SettingsInfo_Intra was not found in InputData_Intra. Please check your input.")
      }
    }

    ## CoRe
    #ValueCol
    if("ValueCol" %in% names(SettingsInfo_CoRe)){
      if(SettingsInfo_CoRe[["ValueCol"]] %in% colnames(InputData_CoRe)== FALSE){
        stop("The ", SettingsInfo_CoRe[["ValueCol"]], " column selected as ValueCol in SettingsInfo_CoRe was not found in InputData_CoRe. Please check your input.")
      }
    }
    #StatCol
    if("StatCol" %in% names(SettingsInfo_CoRe)){
      if(SettingsInfo_CoRe[["StatCol"]] %in% colnames(InputData_CoRe)== FALSE){
        stop("The ", SettingsInfo_CoRe[["StatCol"]], " column selected as StatCol in SettingsInfo_CoRe was not found in InputData_CoRe. Please check your input.")
      }
    }

    #StatCol
    if("DirectionCol" %in% names(SettingsInfo_CoRe)){
      if(SettingsInfo_CoRe[["DirectionCol"]] %in% colnames(InputData_CoRe)== FALSE){
        stop("The ", SettingsInfo_CoRe[["DirectionCol"]], " column selected as DirectionCol in SettingsInfo_CoRe was not found in InputData_CoRe. Please check your input.")
      }
    }

  }

  #------------- SettingsInfo Cutoffs:
  if(is.null(SettingsInfo_C1)==FALSE){
    if(is.na(as.numeric(SettingsInfo_C1[["StatCutoff"]])) == TRUE |as.numeric(SettingsInfo_C1[["StatCutoff"]]) > 1 | as.numeric(SettingsInfo_C1[["StatCutoff"]]) < 0){
      stop("Check input. The selected StatCutoff in SettingsInfo_C1 should be numeric and between 0 and 1.")
    }

    if(is.na(as.numeric(SettingsInfo_C2[["StatCutoff"]])) == TRUE |as.numeric(SettingsInfo_C2[["StatCutoff"]]) > 1 | as.numeric(SettingsInfo_C2[["StatCutoff"]]) < 0){
      stop("Check input. The selected StatCutoff in SettingsInfo_C2 should be numeric and between 0 and 1.")
    }

    if(is.na(as.numeric(SettingsInfo_C1[["ValueCutoff"]])) == TRUE){
      stop("Check input. The selected ValueCutoff in SettingsInfo_C1 should be numeric and between 0 and 1.")
    }

    if(is.na(as.numeric(SettingsInfo_C2[["ValueCutoff"]])) == TRUE){
      stop("Check input. The selected ValueCutoff in SettingsInfo_C2 should be numeric and between 0 and 1.")
    }

  }else{
    if(is.na(as.numeric(SettingsInfo_Intra[["StatCutoff"]])) == TRUE |as.numeric(SettingsInfo_Intra[["StatCutoff"]]) > 1 | as.numeric(SettingsInfo_Intra[["StatCutoff"]]) < 0){
      stop("Check input. The selected StatCutoff in SettingsInfo_Intra should be numeric and between 0 and 1.")
    }

    if(is.na(as.numeric(SettingsInfo_CoRe[["StatCutoff"]])) == TRUE |as.numeric(SettingsInfo_CoRe[["StatCutoff"]]) > 1 | as.numeric(SettingsInfo_CoRe[["StatCutoff"]]) < 0){
      stop("Check input. The selected StatCutoff in SettingsInfo_CoRe should be numeric and between 0 and 1.")
    }

    if(is.na(as.numeric(SettingsInfo_Intra[["ValueCutoff"]])) == TRUE){
      stop("Check input. The selected ValueCutoff in SettingsInfo_Intra should be numeric and between 0 and 1.")
    }

    if(is.na(as.numeric(SettingsInfo_CoRe[["ValueCutoff"]])) == TRUE){
      stop("Check input. The selected ValueCutoff in SettingsInfo_CoRe should be numeric and between 0 and 1.")
    }
  }

  #------------ NAs in data
  if(is.null(InputData_C1)==FALSE){
    if(nrow(InputData_C1[complete.cases(InputData_C1[[SettingsInfo_C1[["ValueCol"]]]], InputData_C1[[SettingsInfo_C1[["StatCol"]]]]), ]) < nrow(InputData_C1)){
      warning("InputData_C1 includes NAs in ", SettingsInfo_C1[["ValueCol"]], " and/or in ", SettingsInfo_C1[["StatCol"]], ". ", nrow(InputData_C1)- nrow(InputData_C1[complete.cases(InputData_C1[[SettingsInfo_C1[["ValueCol"]]]], InputData_C1[[SettingsInfo_C1[["StatCol"]]]]), ]) ," metabolites containing NAs are removed.")
    }

    if(nrow(InputData_C2[complete.cases(InputData_C2[[SettingsInfo_C2[["ValueCol"]]]], InputData_C2[[SettingsInfo_C2[["StatCol"]]]]), ]) < nrow(InputData_C2)){
      warning("InputData_C2 includes NAs in ", SettingsInfo_C2[["ValueCol"]], " and/or in", SettingsInfo_C2[["StatCol"]], ". ", nrow(InputData_C2)- nrow(InputData_C2[complete.cases(InputData_C2[[SettingsInfo_C2[["ValueCol"]]]], InputData_C2[[SettingsInfo_C2[["StatCol"]]]]), ]) ," metabolites containing NAs are removed.")
    }
  }else{
    if(nrow(InputData_Intra[complete.cases(InputData_Intra[[SettingsInfo_Intra[["ValueCol"]]]], InputData_Intra[[SettingsInfo_Intra[["StatCol"]]]]), ]) < nrow(InputData_Intra)){
      warning("InputData_Intra includes NAs in ", SettingsInfo_Intra[["ValueCol"]], " and/or in ", SettingsInfo_Intra[["StatCol"]], ". ", nrow(InputData_Intra)- nrow(InputData_Intra[complete.cases(InputData_Intra[[SettingsInfo_Intra[["ValueCol"]]]], InputData_Intra[[SettingsInfo_Intra[["StatCol"]]]]), ]) ," metabolites containing NAs are removed.")
    }

    if(nrow(InputData_CoRe[complete.cases(InputData_CoRe[[SettingsInfo_CoRe[["ValueCol"]]]], InputData_CoRe[[SettingsInfo_CoRe[["StatCol"]]]]), ]) < nrow(InputData_CoRe)){
      warning("InputData_CoRe includes NAs in ", SettingsInfo_CoRe[["ValueCol"]], " and/or in ", SettingsInfo_CoRe[["StatCol"]], ". ", nrow(InputData_CoRe)- nrow(InputData_CoRe[complete.cases(InputData_CoRe[[SettingsInfo_CoRe[["ValueCol"]]]], InputData_CoRe[[SettingsInfo_CoRe[["StatCol"]]]]), ]) ," metabolites containing NAs are removed.")
    }
  }

 #------------- BackgroundMethod
  if(is.null(SettingsInfo_C1)==FALSE){
    options <- c("C1|C2", "C1&C2", "C2", "C1" , "*")
    if(any(options %in% BackgroundMethod) == FALSE){
        stop("Check input. The selected BackgroundMethod option is not valid. Please select one of the folowwing: ",paste(options,collapse = ", "),"." )
    }
  }else{
    options <- c("Intra|CoRe", "Intra&CoRe", "CoRe", "Intra" , "*")
    if(any(options %in% BackgroundMethod) == FALSE){
      stop("Check input. The selected BackgroundMethod option is not valid. Please select one of the folowwing: ",paste(options,collapse = ", "),"." )
    }
  }

  #------------- SaveAs
  SaveAs_Table_options <- c("txt","csv", "xlsx", "RData")#RData = SummarizedExperiment (?)
  if(is.null(SaveAs_Table)==FALSE){
    if((SaveAs_Table %in% SaveAs_Table_options == FALSE)| (is.null(SaveAs_Table)==TRUE)){
      stop("Check input. The selected SaveAs_Table option is not valid. Please select one of the folowwing: ",paste(SaveAs_Table_options,collapse = ", "),"." )
    }
  }
}
