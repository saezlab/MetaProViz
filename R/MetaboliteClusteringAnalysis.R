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
#' @param Cond1_File DF for your data (results from e.g. DMA) containing metabolites in rows with corresponding Log2FC and stat (p-value, p.adjusted) value columns.
#' @param Cond2_File DF for your data (results from e.g. DMA) containing metabolites in rows with corresponding Log2FC and stat (p-value, p.adjusted) value columns.
#' @param MetaboliteID Column name of Column including the Metabolite identifiers. This MUST BE THE SAME in each of your Input files.
#' @param Cond1ValueCol Column name of Log2FC in Cond1File
#' @param Cond1PadjCol Column name of adjusted p-value in Cond1File. Can also be you p-value column if you want to use this instead.
#' @param Cond2ValueCol  Column name of Log2FC in Cond2File
#' @param Cond2PadjCol Column name of adjusted p-value in Cond2File. Can also be you p-value column if you want to use this instead.
#' @param Cond1_padj_cutoff  \emph{Optional: } adjusted p-value cutoff for Cond1File. \strong{Default=0.05}
#' @param Cond1_FC_cutoff \emph{Optional: } Log2FC cutoff for Cond1File. \strong{Default=0.5}
#' @param Cond2_padj_cutoff \emph{Optional: } adjusted p-value cutoff for Cond2File. \strong{Default=0.05}
#' @param Cond2_FC_cutoff \emph{Optional: } Log2FC cutoff for Cond2File. \strong{Default=0.5}
#' @param backgroundMethod \emph{Optional: } Background method C1|C2, C1&C2, C2, C1 or * \strong{Default="C1&C2"}
#' @param outputFileName \emph{Optional: } Output filename \strong{Default=SiRCle_RCM.csv}
#' @param Save_as_Results \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt" \strong{default: "xlsx"}
#' @param Folder_Name {Optional:} String which is added to the resulting folder name \strong(Default = NULL)
#' @return MCA an instance of the MetaProViz package
#' @export
#'

##################################################
### ### ### Metabolite Clustering Analysis: 2 Conditions ### ### ###
##################################################

MCA_2Cond <- function(Cond1_File,
                      Cond2_File,
                      MetaboliteID= "Metabolite",
                      Cond1ValueCol="Log2FC",
                      Cond1PadjCol="p.adj",
                      Cond2ValueCol="Log2FC",
                      Cond2PadjCol="p.adj",
                      Cond1_padj_cutoff= 0.05,
                      Cond2_padj_cutoff = 0.05,
                      Cond1_FC_cutoff= 1,
                      Cond2_FC_cutoff = 1,
                      Save_as_Results = "xlsx", # txt or csv
                      backgroundMethod="C1&C2",
                      OutputFileName = "MCA_2Cond_",
                      Folder_Name = NULL
                      ){
  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse", "alluvial")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0){
    install.packages(new.packages)
  }
  suppressMessages(library(tidyverse))


  ## ------------ Check Input files ----------- ##

  if( is.numeric(Cond1_padj_cutoff)== FALSE |Cond1_padj_cutoff > 1 | Cond1_padj_cutoff < 0 | is.numeric(Cond2_padj_cutoff)== FALSE |Cond2_padj_cutoff > 1 | Cond2_padj_cutoff < 0){
    stop("Check input. The selected Cond_padj_cutoff value should be numeric and between 0 and 1.")
  }
  if( is.numeric(Cond1_FC_cutoff)== FALSE  | Cond1_FC_cutoff < 0 | is.numeric(Cond2_FC_cutoff)== FALSE  | Cond2_FC_cutoff < 0){
    stop("Check input. The selected Cond_FC_cutoff value should be numeric and positive (absolute value).")
  }

  #Import the data and check columns (here the user will get an error if the column can not be renamed as it does not exists.)
  Cond2_DF <- as.data.frame(Cond2_File)%>%
    dplyr::rename("MetaboliteID"=paste(MetaboliteID),
                  "ValueCol"=paste(Cond2ValueCol),
                  "PadjCol"=paste(Cond2PadjCol))
  Cond1_DF<- as.data.frame(Cond1_File)%>%
    dplyr::rename("MetaboliteID"=paste(MetaboliteID),
                  "ValueCol"=paste(Cond1ValueCol),
                  "PadjCol"=paste(Cond1PadjCol))

  #First check for duplicates in "MetaboliteID" and drop if there are any
  if(length(Cond2_DF[duplicated(Cond2_DF$MetaboliteID), "MetaboliteID"]) > 0){
    doublons <- as.character(Cond2_DF[duplicated(Cond2_DF$MetaboliteID), "MetaboliteID"])#number of duplications
    Cond2_DF <-Cond2_DF[!duplicated(Cond2_DF$MetaboliteID),]#remove duplications
    warning("Cond2 dataset contained duplicates based on MetaboliteID! Dropping duplicate IDs and kept only the first entry. You had ", length(doublons), " duplicates. Note that you should do this before running MCA.")
  }
  if(length(Cond1_DF[duplicated(Cond1_DF$MetaboliteID), "MetaboliteID"]) > 0){
    doublons <- as.character(Cond1_DF[duplicated(Cond1_DF$MetaboliteID), "MetaboliteID"])#number of duplications
    Cond1_DF <-Cond1_DF[!duplicated(Cond1_DF$MetaboliteID),]#remove duplications
    warning("Cond1 dataset contained duplicates based on MetaboliteID! Dropping duplicate IDs and kept only the first entry. You had ", length(doublons), " duplicates. Note that you should do this before running MCA.")
  }

  ## -------- Create Results output folder ---------- ##
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
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)}
  Results_folder_MCA_folder = file.path(Results_folder, "MCA") # select name of result directory
  if (!dir.exists(Results_folder_MCA_folder)) {dir.create(Results_folder_MCA_folder)}  # check and create folder



  ## ------------ Assign Groups -------- ##
  #Tag genes that are detected in each data layer
  Cond2_DF$Detected <- "TRUE"
  Cond1_DF$Detected <- "TRUE"

  #Assign to Group based on individual Cutoff ("UP", "DOWN", "No Change")
  Cond2_DF <- Cond2_DF%>%
    mutate(Cutoff = case_when(Cond2_DF$PadjCol < Cond2_padj_cutoff & Cond2_DF$ValueCol > Cond2_FC_cutoff ~ 'UP',
                              Cond2_DF$PadjCol < Cond2_padj_cutoff & Cond2_DF$ValueCol < -Cond2_FC_cutoff ~ 'DOWN',
                              TRUE ~ 'No Change')) %>%
    mutate(Cutoff_Specific = case_when(Cutoff == "UP" ~ 'UP',
                                       Cutoff == "DOWN" ~ 'DOWN',
                                       Cutoff == "No Change" & Cond2_DF$PadjCol < Cond2_padj_cutoff & Cond2_DF$ValueCol > 0 ~ 'Significant Positive',
                                       Cutoff == "No Change" & Cond2_DF$PadjCol < Cond2_padj_cutoff & Cond2_DF$ValueCol < 0 ~ 'Significant Negative',
                                       Cutoff == "No Change" & Cond2_DF$PadjCol > Cond2_padj_cutoff ~ 'Not Significant',
                                       TRUE ~ 'FALSE'))

  Cond1_DF <-Cond1_DF%>%
    mutate(Cutoff = case_when(Cond1_DF$PadjCol <Cond1_padj_cutoff &Cond1_DF$ValueCol >Cond1_FC_cutoff ~ 'UP',
                              Cond1_DF$PadjCol <Cond1_padj_cutoff &Cond1_DF$ValueCol < -Cond1_FC_cutoff ~ 'DOWN',
                              TRUE ~ 'No Change'))%>%
    mutate(Cutoff_Specific = case_when(Cutoff == "UP" ~ 'UP',
                                       Cutoff == "DOWN" ~ 'DOWN',
                                       Cutoff == "No Change" & Cond1_DF$PadjCol < Cond1_padj_cutoff & Cond1_DF$ValueCol > 0 ~ 'Significant Positive',
                                       Cutoff == "No Change" & Cond1_DF$PadjCol < Cond1_padj_cutoff & Cond1_DF$ValueCol < 0 ~ 'Significant Negative',
                                       Cutoff == "No Change" & Cond1_DF$PadjCol > Cond1_padj_cutoff ~ 'Not Significant',
                                       TRUE ~ 'FALSE'))

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
  if(backgroundMethod == "C1|C2"){# C1|C2 = Cond2 OR Cond1
    MergeDF <- MergeDF%>%
      mutate(BG_Method = case_when(Cond1_DF_Detected=="TRUE" & Cond2_DF_Detected=="TRUE" ~ 'TRUE', #Cond1 & Cond2
                                   Cond1_DF_Detected=="TRUE" & Cond2_DF_Detected=="FALSE" ~ 'TRUE', # JustCond1
                                   Cond1_DF_Detected=="FALSE" & Cond2_DF_Detected=="TRUE" ~ 'TRUE', # Just Cond2
                                   TRUE ~ 'FALSE'))
  }else if(backgroundMethod == "C1&C2"){ # Cond2 AND Cond1
    MergeDF <- MergeDF%>%
      mutate(BG_Method = case_when(Cond1_DF_Detected=="TRUE" & Cond2_DF_Detected=="TRUE" ~ 'TRUE', #Cond1 & Cond2
                                   TRUE ~ 'FALSE'))
  }else if(backgroundMethod == "C2"){ # Cond2 has to be there
    MergeDF <- MergeDF%>%
      mutate(BG_Method = case_when(Cond1_DF_Detected=="TRUE" & Cond2_DF_Detected=="TRUE" ~ 'TRUE', #Cond1 & Cond2
                                   Cond1_DF_Detected=="FALSE" & Cond2_DF_Detected=="TRUE" ~ 'TRUE', # Just Cond2
                                   TRUE ~ 'FALSE'))
  }else if(backgroundMethod == "C1"){ #Cond1 has to be there
    MergeDF <- MergeDF%>%
      mutate(BG_Method = case_when(Cond1_DF_Detected=="TRUE" & Cond2_DF_Detected=="TRUE" ~ 'TRUE', #Cond1 & Cond2
                                   Cond1_DF_Detected=="TRUE" & Cond2_DF_Detected=="FALSE" ~ 'TRUE', # JustCond1
                                   TRUE ~ 'FALSE'))
  }else if(backgroundMethod == "*"){ # Use all genes as the background
    MergeDF$BG_Method <- "TRUE"
  }else{
    stop("Please use one of the following backgroundMethods: C1|C2, C1&C2, C2, C1, *")#error message
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

  Cond2ValueCol_Unique<-paste("Cond2_DF_",Cond2ValueCol)
  Cond2PadjCol_Unique <-paste("Cond2_DF_",Cond2PadjCol)
  Cond1ValueCol_Unique<-paste("Cond1_DF_",Cond1ValueCol)
  Cond1PadjCol_Unique <-paste("Cond1_DF_",Cond1PadjCol)

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
    count(RG1_All, name="Number of Genes")%>%
    dplyr::rename("SiRCle cluster Name"= "RG1_All")
  ClusterSummary_RG1$`Regulation Grouping` <- "RG1_All"

  ClusterSummary_RG2 <- MergeDF_Rearrange[,c("Metabolite", "RG2_Significant")]%>%
    count(RG2_Significant, name="Number of Genes")%>%
    dplyr::rename("SiRCle cluster Name"= "RG2_Significant")
  ClusterSummary_RG2$`Regulation Grouping` <- "RG2_Significant"

  ClusterSummary_RG3 <- MergeDF_Rearrange[,c("Metabolite", "RG3_SignificantChange")]%>%
    count(RG3_SignificantChange, name="Number of Genes")%>%
    dplyr::rename("SiRCle cluster Name"= "RG3_SignificantChange")
  ClusterSummary_RG3$`Regulation Grouping` <- "RG3_SignificantChange"

  ClusterSummary <- rbind(ClusterSummary_RG1, ClusterSummary_RG2,ClusterSummary_RG3)
  ClusterSummary <- ClusterSummary[,c(3,1,2)]


  if (Save_as_Results == "xlsx"){
    xlsMCA <-  paste0(Results_folder_MCA_folder,"/MCA_Output_",OutputFileName,"_Summary.xlsx", sep = "") # Save the DMA results table
    writexl::write_xlsx(ClusterSummary,xlsMCA, col_names = TRUE) # save the DMA result DF

    xlsMCA <-  paste0(Results_folder_MCA_folder,"/MCA_Output_",OutputFileName,".xlsx", sep = "") # Save the DMA results table
    writexl::write_xlsx(MergeDF_Rearrange,xlsMCA, col_names = TRUE) # save the DMA result DF
  }else if (Save_as_Results == "csv"){
    csvMCA <-  paste0(Results_folder_MCA_folder,"/MCA_Output_",OutputFileName,"_Summary.csv", sep = "")
    write.csv(ClusterSummary,csvMCA) # save the DMA result DF

    csvMCA <-  paste0(Results_folder_MCA_folder,"/MCA_Output_",OutputFileName,".csv", sep = "")
    write.csv(MergeDF_Rearrange,csvMCA) # save the DMA result DF
  }else if (Save_as_Results == "txt"){
    txtMCA <-  paste0(Results_folder_MCA_folder,"/MCA_Output_",OutputFileName,"_Summary.txt", sep = "")
    write.table(ClusterSummary,txtMCA, col.names = TRUE, row.names = FALSE) # save the DMA result DF

    txtMCA <-  paste0(Results_folder_MCA_folder,"/MCA_Output_",OutputFileName,".txt", sep = "")
    write.table(MergeDF_Rearrange,txtMCA, col.names = TRUE, row.names = FALSE) # save the DMA result DF
  }

  #Return:
  invisible(return(list("DF"=list("MCA_result"=MergeDF_Rearrange, "ClusterSummary"= ClusterSummary))))
}




#' MCA_CoRe
#'
#' This script allows you to perform metabolite clustering analysis and computes clusters of metabolites based on regulatory rules between Intracellular and culture media metabolomics (CoRe experiment).
#'
#' @param Intra_File DF for your data (results from e.g. DMA) containing metabolites in rows with corresponding Log2FC and stat (p-value, p.adjusted) value columns.
#' @param CoRe_File DF for your data (results from e.g. DMA) containing metabolites in rows with corresponding Log2FC and stat (p-value, p.adjusted) value columns. Here we additionally require
#' @param MetaboliteID Column name of Column including the Metabolite identifiers. This MUST BE THE SAME in each of your Input files.
#' @param IntraValueCol Column name of Log2FC in IntraFile
#' @param IntraPadjCol Column name of adjusted p-value in IntraFile. Can also be you p-value column if you want to use this instead.
#' @param CoReValueCol  Column name of Log2FC in CoReFile
#' @param CoRePadjCol Column name of adjusted p-value in CoReFile. Can also be you p-value column if you want to use this instead.
#' @param Intra_padj_cutoff  \emph{Optional: } adjusted p-value cutoff for IntraFile. \strong{Default=0.05}
#' @param Intra_FC_cutoff \emph{Optional: } Log2FC cutoff for IntraFile. \strong{Default=0.5}
#' @param CoRe_padj_cutoff \emph{Optional: } adjusted p-value cutoff for CoReFile. \strong{Default=0.05}
#' @param CoRe_FC_cutoff \emph{Optional: } Log2FC cutoff for CoReFile. \strong{Default=0.5}
#' @param backgroundMethod \emph{Optional: } Background method `Intra|CoRe, Intra&CoRe, CoRe, Intra or * \strong{Default="C1&C2"}
#' @param outputFileName \emph{Optional: } Output filename \strong{Default=SiRCle_RCM.csv}
#' @param Folder_Name {Optional:} String which is added to the resulting folder name \strong(Default = NULL)
#' @return MCA an instance of the MetaProViz package
#' @export
#'

##################################################
### ### ### Metabolite Clustering Analysis ### ### ###
##################################################

MCA_CoRe <- function(Intra_File,
                     CoRe_File,
                     MetaboliteID= "Metabolite",
                     IntraValueCol="Log2FC",
                     IntraPadjCol="p.adj",
                     CoReValueCol="Log2FC",
                     CoReDirectionCol="CoRe_info",
                     CoRePadjCol="p.adj",
                     Intra_padj_cutoff= 0.05,
                     CoRe_padj_cutoff = 0.05,
                     Intra_FC_cutoff= 1,
                     Save_as_Results = "xlsx",
                     CoRe_FC_cutoff = 1,
                     backgroundMethod="Intra&CoRe",
                     OutputFileName = "MCA_",
                     Folder_Name = NULL){
  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse", "alluvial")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0){
    install.packages(new.packages)
  }
  suppressMessages(library(tidyverse))


  ## ------------ Check Input files ----------- ##

  if( is.numeric(Intra_padj_cutoff)== FALSE |Intra_padj_cutoff > 1 | Intra_padj_cutoff < 0 | is.numeric(CoRe_padj_cutoff)== FALSE |CoRe_padj_cutoff > 1 | CoRe_padj_cutoff < 0){
    stop("Check input. The selected Cond_padj_cutoff value should be numeric and between 0 and 1.")
  }
  if( is.numeric(Intra_FC_cutoff)== FALSE  | Intra_FC_cutoff < 0 | is.numeric(CoRe_FC_cutoff)== FALSE  | CoRe_FC_cutoff < 0){
    stop("Check input. The selected Cond_FC_cutoff value should be numeric and positive (absolute value).")
  }

  #Import the data and check columns (here the user will get an error if the column can not be renamed as it does not exists.)
  CoRe_DF <- as.data.frame(CoRe_File)%>%
    dplyr::rename("MetaboliteID"=paste(MetaboliteID),
                  "ValueCol"=paste(CoReValueCol),
                  "PadjCol"=paste(CoRePadjCol),
                  "CoRe_Direction"=paste(CoReDirectionCol))
  Intra_DF<- as.data.frame(Intra_File)%>%
    dplyr::rename("MetaboliteID"=paste(MetaboliteID),
                  "ValueCol"=paste(IntraValueCol),
                  "PadjCol"=paste(IntraPadjCol))

   #First check for duplicates in "MetaboliteID" and drop if there are any
  if(length(CoRe_DF[duplicated(CoRe_DF$MetaboliteID), "MetaboliteID"]) > 0){
    doublons <- as.character(CoRe_DF[duplicated(CoRe_DF$MetaboliteID), "MetaboliteID"])#number of duplications
    CoRe_DF <-CoRe_DF[!duplicated(CoRe_DF$MetaboliteID),]#remove duplications
    warning("CoRe dataset contained duplicates based on MetaboliteID! Dropping duplicate IDs and kept only the first entry. You had ", length(doublons), " duplicates.")
    warning("Note that you should do this before running MCA.")
  }
  if(length(Intra_DF[duplicated(Intra_DF$MetaboliteID), "MetaboliteID"]) > 0){
    doublons <- as.character(Intra_DF[duplicated(Intra_DF$MetaboliteID), "MetaboliteID"])#number of duplications
    Intra_DF <-Intra_DF[!duplicated(Intra_DF$MetaboliteID),]#remove duplications
    warning("Intra dataset contained duplicates based on MetaboliteID! Dropping duplicate IDs and kept only the first entry. You had ", length(doublons), " duplicates.")
    warning("Note that you should do this before running MCA.")
  }


  ## -------- Create Results output folder ---------- ##
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
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)}
  Results_folder_MCA_folder = file.path(Results_folder, "MCA") # select name of result directory
  if (!dir.exists(Results_folder_MCA_folder)) {dir.create(Results_folder_MCA_folder)}  # check and create folder



  ## ------------ Assign Groups -------- ##
  #Tag genes that are detected in each data layer
  CoRe_DF$Detected <- "TRUE"
  Intra_DF$Detected <- "TRUE"

  #Assign to Group based on individual Cutoff ("UP", "DOWN", "No Change")
  CoRe_DF <- CoRe_DF%>%
    mutate(Cutoff = case_when(CoRe_DF$PadjCol < CoRe_padj_cutoff & CoRe_DF$ValueCol > CoRe_FC_cutoff ~ 'UP',
                              CoRe_DF$PadjCol < CoRe_padj_cutoff & CoRe_DF$ValueCol < -CoRe_FC_cutoff ~ 'DOWN',
                              TRUE ~ 'No Change')) %>%
    mutate(Cutoff_Specific = case_when(Cutoff == "UP" ~ 'UP',
                                       Cutoff == "DOWN" ~ 'DOWN',
                                       Cutoff == "No Change" & CoRe_DF$PadjCol < CoRe_padj_cutoff & CoRe_DF$ValueCol > 0 ~ 'Significant Positive',
                                       Cutoff == "No Change" & CoRe_DF$PadjCol < CoRe_padj_cutoff & CoRe_DF$ValueCol < 0 ~ 'Significant Negative',
                                       Cutoff == "No Change" & CoRe_DF$PadjCol > CoRe_padj_cutoff ~ 'Not Significant',
                                       TRUE ~ 'FALSE'))

  Intra_DF <-Intra_DF%>%
    mutate(Cutoff = case_when(Intra_DF$PadjCol <Intra_padj_cutoff &Intra_DF$ValueCol >Intra_FC_cutoff ~ 'UP',
                              Intra_DF$PadjCol <Intra_padj_cutoff &Intra_DF$ValueCol < -Intra_FC_cutoff ~ 'DOWN',
                              TRUE ~ 'No Change'))%>%
    mutate(Cutoff_Specific = case_when(Cutoff == "UP" ~ 'UP',
                                       Cutoff == "DOWN" ~ 'DOWN',
                                       Cutoff == "No Change" & Intra_DF$PadjCol < Intra_padj_cutoff & Intra_DF$ValueCol > 0 ~ 'Significant Positive',
                                       Cutoff == "No Change" & Intra_DF$PadjCol < Intra_padj_cutoff & Intra_DF$ValueCol < 0 ~ 'Significant Negative',
                                       Cutoff == "No Change" & Intra_DF$PadjCol > Intra_padj_cutoff ~ 'Not Significant',
                                       TRUE ~ 'FALSE'))

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
    mutate_at(c("CoRe_Direction"), ~replace_na(.,"Not Detected"))

  #Apply Background filter (label metabolites that will be removed based on chosen background)
  if(backgroundMethod == "Intra|CoRe"){# C1|C2 = CoRe OR Intra
    MergeDF <- MergeDF%>%
      mutate(BG_Method = case_when(Intra_DF_Detected=="TRUE" & CoRe_DF_Detected=="TRUE" ~ 'TRUE', #Intra & CoRe
                                   Intra_DF_Detected=="TRUE" & CoRe_DF_Detected=="FALSE" ~ 'TRUE', # JustIntra
                                   Intra_DF_Detected=="FALSE" & CoRe_DF_Detected=="TRUE" ~ 'TRUE', # Just CoRe
                                   TRUE ~ 'FALSE'))
  }else if(backgroundMethod == "Intra&CoRe"){ # CoRe AND Intra
    MergeDF <- MergeDF%>%
      mutate(BG_Method = case_when(Intra_DF_Detected=="TRUE" & CoRe_DF_Detected=="TRUE" ~ 'TRUE', #Intra & CoRe
                                   TRUE ~ 'FALSE'))
  }else if(backgroundMethod == "CoRe"){ # CoRe has to be there
    MergeDF <- MergeDF%>%
      mutate(BG_Method = case_when(Intra_DF_Detected=="TRUE" & CoRe_DF_Detected=="TRUE" ~ 'TRUE', #Intra & CoRe
                                   Intra_DF_Detected=="FALSE" & CoRe_DF_Detected=="TRUE" ~ 'TRUE', # Just CoRe
                                   TRUE ~ 'FALSE'))
  }else if(backgroundMethod == "Intra"){ #Intra has to be there
    MergeDF <- MergeDF%>%
      mutate(BG_Method = case_when(Intra_DF_Detected=="TRUE" & CoRe_DF_Detected=="TRUE" ~ 'TRUE', #Intra & CoRe
                                   Intra_DF_Detected=="TRUE" & CoRe_DF_Detected=="FALSE" ~ 'TRUE', # JustIntra
                                   TRUE ~ 'FALSE'))
  }else if(backgroundMethod == "*"){ # Use all metabolites as the background
    MergeDF$BG_Method <- "TRUE"
  }else{
    stop("Please use one of the following backgroundMethods: Intra|CoRe, Intra&CoRe, CoRe, Intra, *")#error message
  }

  #Assign Metabolite cluster names to the metabolites
  MergeDF <- MergeDF%>%
        mutate(RG1_All = case_when(BG_Method =="FALSE"~ 'Background = FALSE',
                               Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Released" ~ 'Intra DOWN + CoRe DOWN_Released',
                               Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected"~ 'Intra DOWN + CoRe Not Detected',
                               Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Released"~ 'Intra DOWN + CoRe Not Significant_Released',
                               Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Released"~ 'Intra DOWN + CoRe Significant Negative_Released',
                               Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Released"~ 'Intra DOWN + CoRe Significant Positive_Released',
                               Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Released" ~ 'Intra DOWN + CoRe UP_Released',

                               Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Released" ~ 'Intra UP + CoRe DOWN_Released',
                               Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected" ~ 'Intra UP + CoRe Not Detected',
                               Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Released" ~ 'Intra UP + CoRe Not Significant_Released',
                               Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Released" ~ 'Intra UP + CoRe Significant Negative_Released',
                               Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Released" ~ 'Intra UP + CoRe Significant Positive_Released',
                               Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Released" ~ 'Intra UP + CoRe UP_Released',

                               Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Released" ~ 'Intra Not Detected + CoRe DOWN_Released',
                               Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected" ~ 'Intra Not Detected + CoRe Not Detected',
                               Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Released" ~ 'Intra Not Detected + CoRe Not Significant_Released',
                               Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Released" ~ 'Intra Not Detected + CoRe Significant Negative_Released',
                               Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Released" ~ 'Intra Not Detected + CoRe Significant Positive_Released',
                               Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Released" ~ 'Intra Not Detected + CoRe UP_Released',

                               Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Released"~ 'Intra Significant Negative + CoRe DOWN_Released',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected"~ 'Intra Significant Negative + CoRe Not Detected',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Released"~ 'Intra Significant Negative + CoRe Not Significant_Released',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Released"~ 'Intra Significant Negative + CoRe Significant Negative_Released',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Released"~ 'Intra Significant Negative + CoRe Significant Positive_Released',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Released"~ 'Intra Significant Negative + CoRe UP_Released',

                               Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Released"~ 'Intra Significant Positive + CoRe DOWN_Released',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected"~ 'Intra Significant Positive + CoRe Not Detected',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Released"~ 'Intra Significant Positive + CoRe Not Significant_Released',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Released"~ 'Intra Significant Positive + CoRe Significant Negative_Released',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Released"~ 'Intra Significant Positive + CoRe Significant Positive_Released',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Released"~ 'Intra Significant Positive + CoRe UP_Released',

                               Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Released"~ 'Intra Not Significant + CoRe DOWN_Released',
                               Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected"~ 'Intra Not Significant + CoRe Not Detected',
                               Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Released"~ 'Intra Not Significant + CoRe Not Significant_Released',
                               Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Released"~ 'Intra Not Significant + CoRe Significant Negative_Released',
                               Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Released"~ 'Intra Not Significant + CoRe Significant Positive_Released',
                               Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Released"~ 'Intra Not Significant + CoRe UP_Released',

                               #Consumed:
                               Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Consumed" ~ 'Intra DOWN + CoRe DOWN_Consumed',
                               Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected"~ 'Intra DOWN + CoRe Not Detected',
                               Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Consumed"~ 'Intra DOWN + CoRe Not Significant_Consumed',
                               Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Consumed"~ 'Intra DOWN + CoRe Significant Negative_Consumed',
                               Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Consumed"~ 'Intra DOWN + CoRe Significant Positive_Consumed',
                               Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Consumed" ~ 'Intra DOWN + CoRe UP_Consumed',

                               Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Consumed" ~ 'Intra UP + CoRe DOWN_Consumed',
                               Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected" ~ 'Intra UP + CoRe Not Detected',
                               Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Consumed" ~ 'Intra UP + CoRe Not Significant_Consumed',
                               Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Consumed" ~ 'Intra UP + CoRe Significant Negative_Consumed',
                               Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Consumed" ~ 'Intra UP + CoRe Significant Positive_Consumed',
                               Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Consumed" ~ 'Intra UP + CoRe UP_Consumed',

                               Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Consumed" ~ 'Intra Not Detected + CoRe DOWN_Consumed',
                               Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected" ~ 'Intra Not Detected + CoRe Not Detected',
                               Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Consumed" ~ 'Intra Not Detected + CoRe Not Significant_Consumed',
                               Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Consumed" ~ 'Intra Not Detected + CoRe Significant Negative_Consumed',
                               Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Consumed" ~ 'Intra Not Detected + CoRe Significant Positive_Consumed',
                               Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Consumed" ~ 'Intra Not Detected + CoRe UP_Consumed',

                               Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Consumed"~ 'Intra Significant Negative + CoRe DOWN_Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected"~ 'Intra Significant Negative + CoRe Not Detected',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Consumed"~ 'Intra Significant Negative + CoRe Not Significant_Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Consumed"~ 'Intra Significant Negative + CoRe Significant Negative_Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Consumed"~ 'Intra Significant Negative + CoRe Significant Positive_Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Consumed"~ 'Intra Significant Negative + CoRe UP_Consumed',

                               Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Consumed"~ 'Intra Significant Positive + CoRe DOWN_Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected"~ 'Intra Significant Positive + CoRe Not Detected',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Consumed"~ 'Intra Significant Positive + CoRe Not Significant_Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Consumed"~ 'Intra Significant Positive + CoRe Significant Negative_Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Consumed"~ 'Intra Significant Positive + CoRe Significant Positive_Consumed',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Consumed"~ 'Intra Significant Positive + CoRe UP_Consumed',

                               Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Consumed"~ 'Intra Not Significant + CoRe DOWN_Consumed',
                               Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected"~ 'Intra Not Significant + CoRe Not Detected',
                               Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Consumed"~ 'Intra Not Significant + CoRe Not Significant_Consumed',
                               Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Consumed"~ 'Intra Not Significant + CoRe Significant Negative_Consumed',
                               Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Consumed"~ 'Intra Not Significant + CoRe Significant Positive_Consumed',
                               Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Consumed"~ 'Intra Not Significant + CoRe UP_Consumed',

                               #Consumed/Released (Consumed in one, released in the other)
                               Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Consumed/Released" ~ 'Intra DOWN + CoRe DOWN_Consumed/Released',
                               Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected"~ 'Intra DOWN + CoRe Not Detected',
                               Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Consumed/Released"~ 'Intra DOWN + CoRe Not Significant_Consumed/Released',
                               Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Consumed/Released"~ 'Intra DOWN + CoRe Significant Negative_Consumed/Released',
                               Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Consumed/Released"~ 'Intra DOWN + CoRe Significant Positive_Consumed/Released',
                               Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Consumed/Released" ~ 'Intra DOWN + CoRe UP_Consumed/Released',

                               Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Consumed/Released" ~ 'Intra UP + CoRe DOWN_Consumed/Released',
                               Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected" ~ 'Intra UP + CoRe Not Detected',
                               Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Consumed/Released" ~ 'Intra UP + CoRe Not Significant_Consumed/Released',
                               Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Consumed/Released" ~ 'Intra UP + CoRe Significant Negative_Consumed/Released',
                               Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Consumed/Released" ~ 'Intra UP + CoRe Significant Positive_Consumed/Released',
                               Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Consumed/Released" ~ 'Intra UP + CoRe UP_Consumed/Released',

                               Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Consumed/Released" ~ 'Intra Not Detected + CoRe DOWN_Consumed/Released',
                               Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected" ~ 'Intra Not Detected + CoRe Not Detected',
                               Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Consumed/Released" ~ 'Intra Not Detected + CoRe Not Significant_Consumed/Released',
                               Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Consumed/Released" ~ 'Intra Not Detected + CoRe Significant Negative_Consumed/Released',
                               Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Consumed/Released" ~ 'Intra Not Detected + CoRe Significant Positive_Consumed/Released',
                               Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Consumed/Released" ~ 'Intra Not Detected + CoRe UP_Consumed/Released',

                               Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Consumed/Released"~ 'Intra Significant Negative + CoRe DOWN_Consumed/Released',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected"~ 'Intra Significant Negative + CoRe Not Detected',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Consumed/Released"~ 'Intra Significant Negative + CoRe Not Significant_Consumed/Released',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Consumed/Released"~ 'Intra Significant Negative + CoRe Significant Negative_Consumed/Released',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Consumed/Released"~ 'Intra Significant Negative + CoRe Significant Positive_Consumed/Released',
                               Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Consumed/Released"~ 'Intra Significant Negative + CoRe UP_Consumed/Released',

                               Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Consumed/Released"~ 'Intra Significant Positive + CoRe DOWN_Consumed/Released',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected"~ 'Intra Significant Positive + CoRe Not Detected',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Consumed/Released"~ 'Intra Significant Positive + CoRe Not Significant_Consumed/Released',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Consumed/Released"~ 'Intra Significant Positive + CoRe Significant Negative_Consumed/Released',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Consumed/Released"~ 'Intra Significant Positive + CoRe Significant Positive_Consumed/Released',
                               Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Consumed/Released"~ 'Intra Significant Positive + CoRe UP_Consumed/Released',

                               Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Consumed/Released"~ 'Intra Not Significant + CoRe DOWN_Consumed/Released',
                               Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected"~ 'Intra Not Significant + CoRe Not Detected',
                               Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Consumed/Released"~ 'Intra Not Significant + CoRe Not Significant_Consumed/Released',
                               Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Consumed/Released"~ 'Intra Not Significant + CoRe Significant Negative_Consumed/Released',
                               Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Consumed/Released"~ 'Intra Not Significant + CoRe Significant Positive_Consumed/Released',
                               Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Consumed/Released"~ 'Intra Not Significant + CoRe UP_Consumed/Released',
                               TRUE ~ 'NA'))%>%
    mutate(RG2_Significant = case_when(BG_Method =="FALSE"~ 'Background = FALSE',
                                       Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Released" ~ 'Both_DOWN (Released)',
                                       Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected"~ 'None',
                                       Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Released"~ 'None',
                                       Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Released"~ 'Both_DOWN (Released)',
                                       Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Released"~ 'Opposite (Released UP)',
                                       Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Released" ~ 'Opposite (Released UP)',

                                       Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Released" ~ 'Opposite (Released DOWN)',
                                       Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Released" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Released" ~ 'Opposite (Released UP)',
                                       Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Released" ~ 'Both_UP (Released)',
                                       Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Released" ~ 'Both_UP (Released)',

                                       Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Released" ~ 'CoRe_DOWN (Released)',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Released" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Released" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Released" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Released" ~ 'CoRe_UP (Released)',

                                       Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Released"~ 'Both_DOWN (Released)',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Released"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Released"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Released"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Released"~ 'Opposite (Released UP)',

                                       Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Released"~ 'Opposite (Released DOWN)',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Released"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Released"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Released"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Released"~ 'Both_UP (Released)',

                                       Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Released"~ 'CoRe_DOWN (Released)',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Released"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Released"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Released"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Released"~ 'CoRe_UP (Released)',

                                       #Consumed:
                                       Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Consumed" ~ 'Both_DOWN (Consumed)',
                                       Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected"~ 'None',
                                       Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Consumed"~ 'Both_DOWN (Consumed)',
                                       Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Consumed"~ 'Opposite (Consumed UP)',
                                       Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Consumed" ~ 'Opposite (Consumed UP)',

                                       Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Consumed" ~ 'Opposite (Consumed DOWN)',
                                       Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Consumed" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Consumed" ~ 'Opposite (Consumed UP)',
                                       Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Consumed" ~ 'Both_UP (Consumed)',
                                       Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Consumed" ~ 'Both_UP (Consumed)',

                                       Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Consumed" ~ 'CoRe_DOWN (Consumed)',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Consumed" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Consumed" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Consumed" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Consumed" ~ 'CoRe_UP (Consumed)',

                                       Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Consumed"~ 'Both_DOWN (Consumed)',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Consumed"~ 'Opposite (Consumed UP)',

                                       Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Consumed"~ 'Opposite (Consumed DOWN)',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Consumed"~ 'Both_UP (Consumed)',

                                       Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Consumed"~ 'CoRe_DOWN (Consumed)',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Consumed"~ 'CoRe_UP (Consumed)',

                                       #Consumed/Released (Consumed in one, released in the other)
                                       Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Released/Consumed" ~ 'Both_DOWN (Released/Consumed)',
                                       Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected"~ 'None',
                                       Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Released/Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Released/Consumed"~ 'Both_DOWN (Released/Consumed)',
                                       Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Released/Consumed"~ 'Opposite (Released/Consumed UP)',
                                       Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Released/Consumed" ~ 'Opposite (Released/Consumed UP)',

                                       Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Released/Consumed" ~ 'Opposite (Released/Consumed DOWN)',
                                       Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Released/Consumed" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Released/Consumed" ~ 'Opposite (Released/Consumed UP)',
                                       Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Released/Consumed" ~ 'Both_UP (Released/Consumed)',
                                       Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Released/Consumed" ~ 'Both_UP (Released/Consumed)',

                                       Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Released/Consumed" ~ 'CoRe_DOWN (Released/Consumed)',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Released/Consumed" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Released/Consumed" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Released/Consumed" ~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Released/Consumed" ~ 'CoRe_UP (Released/Consumed)',

                                       Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Released/Consumed"~ 'Both_DOWN (Released/Consumed)',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Released/Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Released/Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Released/Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Released/Consumed"~ 'Opposite (Released/Consumed UP)',

                                       Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Released/Consumed"~ 'Opposite (Released/Consumed DOWN)',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Released/Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Released/Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Released/Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Released/Consumed"~ 'Both_UP (Released/Consumed)',

                                       Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Released/Consumed"~ 'CoRe_DOWN (Released/Consumed)',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Released/Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Released/Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Released/Consumed"~ 'None',
                                       Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Released/Consumed"~ 'CoRe_UP (Released/Consumed)',
                                       TRUE ~ 'NA'))%>%
    mutate(RG3_Change = case_when(BG_Method =="FALSE"~ 'Background = FALSE',
                                  Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Released" ~ 'Both_DOWN (Released)',
                                  Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected"~ 'None',
                                  Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Released"~ 'None',
                                  Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Released"~ 'None',
                                  Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Released"~ 'None',
                                  Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Released" ~ 'Opposite (Released UP)',

                                  Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Released" ~ 'Opposite (Released DOWN)',
                                  Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Released" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Released" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Released" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Released" ~ 'Both_UP (Released)',

                                  Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Released" ~ 'CoRe_DOWN (Released)',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Released" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Released" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Released" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Released" ~ 'CoRe_UP (Released)',

                                  Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Released"~ 'CoRe_DOWN (Released)',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Released"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Released"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Released"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Released"~ 'CoRe_UP (Released)',

                                  Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Released"~ 'CoRe_DOWN (Released)',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Released"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Released"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Released"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Released"~ 'CoRe_UP (Released)',

                                  Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Released"~ 'CoRe_DOWN (Released)',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Released"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Released"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Released"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Released"~ 'CoRe_UP (Released)',

                                  #Consumed:
                                  Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Consumed" ~ 'Both_DOWN (Consumed)',
                                  Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected"~ 'None',
                                  Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Consumed" ~ 'Opposite (Consumed UP)',

                                  Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Consumed" ~ 'Opposite (Consumed DOWN)',
                                  Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Consumed" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Consumed" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Consumed" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Consumed" ~ 'Both_UP (Consumed)',

                                  Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Consumed" ~ 'CoRe_DOWN (Consumed)',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Consumed" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Consumed" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Consumed" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Consumed" ~ 'CoRe_UP (Consumed)',

                                  Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Consumed"~ 'CoRe_DOWN (Consumed)',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Consumed"~ 'CoRe_UP (Consumed)',

                                  Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Consumed"~ 'CoRe_DOWN (Consumed)',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Consumed"~ 'CoRe_UP (Consumed)',

                                  Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Consumed"~ 'CoRe_DOWN (Consumed)',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Consumed"~ 'CoRe_UP (Consumed)',

                                  #Consumed/Released (Consumed in one, released in the other)
                                  Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Released/Consumed" ~ 'Both_DOWN (Released/Consumed)',
                                  Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected"~ 'None',
                                  Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Released/Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Released/Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Released/Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="DOWN" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Released/Consumed" ~ 'Opposite (Released/Consumed UP)',

                                  Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Released/Consumed" ~ 'Opposite (Released/Consumed DOWN)',
                                  Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Released/Consumed" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Released/Consumed" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Released/Consumed" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="UP" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Released/Consumed" ~ 'Both_UP (Released/Consumed)',

                                  Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Released/Consumed" ~ 'CoRe_DOWN (Released/Consumed)',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Released/Consumed" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Released/Consumed" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Released/Consumed" ~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Detected" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Released/Consumed" ~ 'CoRe_UP (Released/Consumed)',

                                  Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Released/Consumed"~ 'CoRe_DOWN (Released/Consumed)',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Released/Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Released/Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Released/Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Negative" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Released/Consumed"~ 'CoRe_UP (Released/Consumed)',

                                  Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Released/Consumed"~ 'CoRe_DOWN (Released/Consumed)',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Released/Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Released/Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Released/Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Significant Positive" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Released/Consumed"~ 'CoRe_UP (Released/Consumed)',

                                  Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="DOWN" & CoRe_Direction=="Released/Consumed"~ 'CoRe_DOWN (Released/Consumed)',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Not Detected" & CoRe_Direction=="Not Detected"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Not Significant" & CoRe_Direction=="Released/Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Significant Negative" & CoRe_Direction=="Released/Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="Significant Positive" & CoRe_Direction=="Released/Consumed"~ 'None',
                                  Intra_DF_Cutoff_Specific=="Not Significant" & CoRe_DF_Cutoff_Specific=="UP" & CoRe_Direction=="Released/Consumed"~ 'CoRe_UP (Released/Consumed)',
                                  TRUE ~ 'NA'))

  #Safe the DF and return the groupings
  ##MCA DF (Merged InputDF filtered for background with assigned MCA cluster names)
  MergeDF_Select1 <- MergeDF[, c("MetaboliteID", "Intra_DF_Detected","Intra_DF_ValueCol","Intra_DF_PadjCol","Intra_DF_Cutoff", "Intra_DF_Cutoff_Specific", "CoRe_DF_Detected", "CoRe_DF_ValueCol","CoRe_DF_PadjCol","CoRe_DF_Cutoff", "CoRe_DF_Cutoff_Specific", "BG_Method", "RG1_All", "RG2_Significant", "RG3_SignificantChange")]

  CoReValueCol_Unique<-paste("CoRe_DF_",CoReValueCol)
  CoRePadjCol_Unique <-paste("CoRe_DF_",CoRePadjCol)
  IntraValueCol_Unique<-paste("Intra_DF_",IntraValueCol)
  IntraPadjCol_Unique <-paste("Intra_DF_",IntraPadjCol)

  MergeDF_Select2<- subset(MergeDF, select=-c(Intra_DF_Detected,Intra_DF_Cutoff, CoRe_DF_Detected,CoRe_DF_Cutoff, CoRe_DF_Cutoff_Specific, BG_Method, RG1_All, RG2_Significant, RG3_SignificantChange))%>%
    dplyr::rename(!!CoReValueCol_Unique :="CoRe_DF_ValueCol",#This syntax is needed since paste(MetaboliteID)="MetaboliteID" is not working in dyplr
                  !!CoRePadjCol_Unique :="CoRe_DF_PadjCol",
                  !!IntraValueCol_Unique :="Intra_DF_ValueCol",
                  !!IntraPadjCol_Unique :="Intra_DF_PadjCol")

  MergeDF_Rearrange <- merge(MergeDF_Select1, MergeDF_Select2, by="MetaboliteID")
  MergeDF_Rearrange <-MergeDF_Rearrange%>%
    dplyr::rename("Metabolite"="MetaboliteID")


  ##Summary SiRCle clusters (number of genes assigned to each SiRCle cluster in each grouping)
  ClusterSummary_RG1 <- MergeDF_Rearrange[,c("Metabolite", "RG1_All")]%>%
    count(RG1_All, name="Number of Genes")%>%
    dplyr::rename("SiRCle cluster Name"= "RG1_All")
  ClusterSummary_RG1$`Regulation Grouping` <- "RG1_All"

  ClusterSummary_RG2 <- MergeDF_Rearrange[,c("Metabolite", "RG2_Significant")]%>%
    count(RG2_Significant, name="Number of Genes")%>%
    dplyr::rename("SiRCle cluster Name"= "RG2_Significant")
  ClusterSummary_RG2$`Regulation Grouping` <- "RG2_Significant"

  ClusterSummary_RG3 <- MergeDF_Rearrange[,c("Metabolite", "RG3_SignificantChange")]%>%
    count(RG3_SignificantChange, name="Number of Genes")%>%
    dplyr::rename("SiRCle cluster Name"= "RG3_SignificantChange")
  ClusterSummary_RG3$`Regulation Grouping` <- "RG3_SignificantChange"

  ClusterSummary <- rbind(ClusterSummary_RG1, ClusterSummary_RG2,ClusterSummary_RG3)
  ClusterSummary <- ClusterSummary[,c(3,1,2)]


  if (Save_as_Results == "xlsx"){
    xlsMCA <-  paste0(Results_folder_MCA_folder,"/MCA_Output_",OutputFileName,"_Summary_CoRe.xlsx", sep = "") # Save the DMA results table
    writexl::write_xlsx(ClusterSummary,xlsMCA, col_names = TRUE) # save the DMA result DF

    xlsMCA <-  paste0(Results_folder_MCA_folder,"/MCA_Output_",OutputFileName,"_CoRe.xlsx", sep = "") # Save the DMA results table
    writexl::write_xlsx(MergeDF_Rearrange,xlsMCA, col_names = TRUE) # save the DMA result DF
  }else if (Save_as_Results == "csv"){
    csvMCA <-  paste0(Results_folder_MCA_folder,"/MCA_Output_",OutputFileName,"_Summary_CoRe.csv", sep = "")
    write.csv(ClusterSummary,csvMCA) # save the DMA result DF

    csvMCA <-  paste0(Results_folder_MCA_folder,"/MCA_Output_",OutputFileName,"_CoRe.csv", sep = "")
    write.csv(MergeDF_Rearrange,csvMCA) # save the DMA result DF
  }else if (Save_as_Results == "txt"){
    txtMCA <-  paste0(Results_folder_MCA_folder,"/MCA_Output_",OutputFileName,"_Summary_CoRe.txt", sep = "")
    write.table(ClusterSummary,txtMCA, col.names = TRUE, row.names = FALSE) # save the DMA result DF

    txtMCA <-  paste0(Results_folder_MCA_folder,"/MCA_Output_",OutputFileName,"_CoRe.txt", sep = "")
    write.table(MergeDF_Rearrange,txtMCA, col.names = TRUE, row.names = FALSE) # save the DMA result DF
  }


  invisible(return(list("DF"=list("MCA_result"=MergeDF_Rearrange, "ClusterSummary"= ClusterSummary))))
}



