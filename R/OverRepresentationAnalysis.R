#'## ---------------------------
##
## Script name: Over representation Analysis (ORA)
##
## Purpose of script: Run ORA on MetaProViz metabolite clusters from MCA or diffeential results from DMA
##
## Author: Christina Schmidt
##
## Date Created: 2023-07-03
##
## Copyright (c)
## Email:
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------
#'
#' This script allows you


#################################
### ### ### MC_ORA ### ### ###
#################################
#' MC_ORA
#'
#' Uses enricher to run ORA on each of the metabolite cluster from any of the MCA functions using a pathway list
#'
#' @param Input_data Input DF with column "Metabolite", which needs to match the identifier type in column ""Metabolite" in the PathwayFile.
#' @param MetabCluster_lab \emph{Optional: } Provide column name for the metabolite cluster labels in the Input_data.\strong{default: "RG2_Significant"}
#' @param RemoveBackground\emph{Optional: } If TRUE, column "BG_Method" needs to be in Input_data, which includes TRUE/FALSE for each metabolite to fall into background based on the chosen Background method for e.g. MCA_2Cond are removed from the universe. \strong{default: TRUE}
#' @param PathwayFile DF that must include column "term" with the pathway name, column "Metabolite" with the Metabolite name or ID and column "Description" with pathway description that will be depicted on the plots.
#' @param PathwayName \emph{Optional: } Name of the pathway list used \strong{default: ""}
#' @param minGSSize \emph{Optional: } minimum group size in ORA \strong{default: 10}
#' @param maxGSSize \emph{Optional: } maximum group size in ORA \strong{default: 1000}
#' @param Save_as_Plot \emph{Optional: } If Save_as_Plot=NULL no plots will be saved. Otherwise, file types for the figures are: "svg", "pdf", "png" \strong{default: "pdf"}
#' @param Save_as_Results \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt" \strong{default: "csv"}
#' @param Folder_Name {Optional:} String which is added to the resulting folder name \strong(Default = NULL)
#'
#' @return Saves results as individual .csv files.
#' @export

MC_ORA <- function(Input_data,
                   MetabCluster_lab="RG2_Significant",
                   RemoveBackground=TRUE,
                   PathwayFile,
                   PathwayName="",
                   minGSSize=10,
                   maxGSSize=1000 ,
                   Save_as_Results= "csv",
                   Folder_Name = NULL
                   ){
  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse","clusterProfiler", "ggupset")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  suppressMessages(library(tidyverse))
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0){
    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(new.packages)
  }
  ## ------------ Check Input files ----------- ##
  # 1. The input data:
  if("Metabolite" %in% names(Input_data)==FALSE){
    stop("Check input. Input_data must contain a column named `Metabolite` including the metabolite names.")
  }
  if(length(Input_data[duplicated(Input_data$Metabolite), "Metabolite"]) > 0){
    doublons <- as.character(Input_data[duplicated(Input_data$Metabolite), "Metabolite"])#number of duplications
    Input_data <-Input_data[!duplicated(Input_data$Metabolite),]#remove duplications
    warning("Input_data contained duplicates based on Metabolite! Dropping duplicate IDs and kept only the first entry. You had ", length(doublons), " duplicates. Note that you should do this before running VizVolcano.")
  }
  if(MetabCluster_lab %in% names(Input_data)==FALSE){
    stop("MetabCluster_lab must be column in Input_data.")
    }
 if(is.logical(RemoveBackground) == FALSE){
    stop("Check input. RemoveBackground value should be either =TRUE or = FALSE.")
  } else if(RemoveBackground==TRUE & "BG_Method" %in% names(Input_data)==FALSE){
    stop("Check input. Input_data must contain a column named `BG_Method` with entries TRUE or FALSE.")
  }


  # 2. PathwayFile and name
  if("Metabolite" %in% names(PathwayFile)==FALSE){
    stop("Check input. PathwayFile must contain a column named `Metabolite` including the metabolite names or IDs.")
  } else if("Metabolite" %in% names(PathwayFile)==TRUE){
    PathwayFile<-PathwayFile%>%
      dplyr::rename("gene"="Metabolite")
  }
  if("term" %in% names(PathwayFile)==FALSE){
    stop("Check input. PathwayFile must contain a column named `term` with pathway names.")
  }
  if("Description" %in% names(PathwayFile)==FALSE){
    stop("Check input. PathwayFile must contain a column named `Description` with description of the pathway.")
  }
  if(is.character(PathwayName)==FALSE){
    stop("Check input. PathwayName must be a character of syntax 'example'.")
  }

  # 3. General Settings
  if(is.numeric(minGSSize)== FALSE){
    stop("Check input. The selected minGSSize value should be numeric.")
  }
  if(is.numeric(maxGSSize)== FALSE){
    stop("Check input. The selected maxGSSize value should be numeric.")
  }
  Save_as_Results_options <- c("txt","csv", "xlsx" )
  if(Save_as_Results %in% Save_as_Results_options == FALSE){
    stop("Check input. The selected Save_as_Results option is not valid. Please select one of the folowwing: ",paste(Save_as_Results_options,collapse = ", "),"." )
  }

  ## ------------ Create Output folders ----------- ##
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
  Results_folder_MC_ORA = file.path(Results_folder, "MC_ORA")  # This searches for a folder within the "Results" folder in the current working directory and if its not found it creates one
  if (!dir.exists( Results_folder_MC_ORA)) {dir.create( Results_folder_MC_ORA)}  # check and create folder

  ############################################################################################################
  ## ------------ Run ----------- ##
  # open the data
  if(RemoveBackground==TRUE){
    df <- subset(Input_data, !Input_data$BG_Method == "FALSE")
  } else{
    df <- Input_data
  }

  #Select universe
  allMetabolites <- as.character(df$Metabolite)

  #Select MCA clusters
  grps_labels <- unlist(unique(df[[MetabCluster_lab]]))

  #Load Pathways
  Pathway <- PathwayFile
  Term2gene <- Pathway[,c("term", "gene")]# term and MetaboliteID (MetaboliteID= gene as syntax required for enricher)
  term2name <- Pathway[,c("term", "Description")]# term and description

  #Add the number of genes present in each pathway
  Pathway$Count <- 1
  Pathway_Mean <- aggregate(Pathway$Count, by=list(term=Pathway$term), FUN=sum)
  names(Pathway_Mean)[names(Pathway_Mean) == "x"] <- "Metabolites_in_Pathway"
  Pathway <- merge(x= Pathway[,-4], y=Pathway_Mean,by="term", all.x=TRUE)

  df_list = list()# Make an empty list to store the created DFs
  #Run ORA
  for(g in grps_labels){
    grpMetabolites <- subset(df, df[[MetabCluster_lab]] == g)
    message("Number of metabolite in cluster `",g, "`: ", nrow(grpMetabolites), sep="")

    clusterGo <- clusterProfiler::enricher(gene=as.character(grpMetabolites$Metabolite),
                                           pvalueCutoff = 1,
                                           pAdjustMethod = "BH",
                                           universe = allMetabolites,
                                           minGSSize=minGSSize,
                                           maxGSSize=maxGSSize,
                                           qvalueCutoff = 1,
                                           TERM2GENE=Term2gene ,
                                           TERM2NAME = term2name)
    clusterGoSummary <- data.frame(clusterGo)
    if(!(dim(clusterGoSummary)[1] == 0)){
      #Add pathway information (% of genes in pathway detected)
      clusterGoSummary <- merge(x= clusterGoSummary[,-2], y=Pathway[,-2],by.x="ID",by.y="term", all=TRUE)
      clusterGoSummary$Count[is.na(clusterGoSummary$Count)] <- 0
      clusterGoSummary$Percentage_of_Pathway_detected <-round(((clusterGoSummary$Count/clusterGoSummary$Metabolites_in_Pathway)*100),digits=2)
      clusterGoSummary <- clusterGoSummary[!duplicated(clusterGoSummary$ID),]
      clusterGoSummary <- clusterGoSummary[order(clusterGoSummary$p.adjust),]
      clusterGoSummary <- clusterGoSummary[,c(1,9,2:8, 10:11)]%>%
        dplyr::rename("MetaboliteIDs_in_pathway"="geneID")
    }else{
      message("None of the Input_data Metabolites of the cluster ", g ," were present in any terms of the PathwayFile. Hence the ClusterGoSummary ouput will be empty for this cluster. Please check that the metabolite IDs match the pathway IDs.")
    }

      #Save file
      g_save <- gsub("/", "-", g)

      if (Save_as_Results == "csv"){
        if(PathwayName ==""){
          csvORA <-  paste(Results_folder_MC_ORA,"/", 'ClusterGoSummary_', g_save, '.csv', sep="")
          write_csv(clusterGoSummary,csvORA )#Export the ORA results as .csv
        } else{
          csvORA <-  paste(Results_folder_MC_ORA,"/", 'ClusterGoSummary_', PathwayName, '-', g_save, '.csv', sep="")#Export the ORA results as .csv
          write_csv(clusterGoSummary,csvORA )
        }
      }else if(Save_as_Results == "txt"){
      if(PathwayName ==""){
        txtORA <-  paste(Results_folder_MC_ORA,"/", 'ClusterGoSummary_', g_save, '.txt', sep="")
        write.table(clusterGoSummary,txtORA, col.names = TRUE, row.names = FALSE) #Export the ORA results as txt
      } else{
        txtORA <-  paste(Results_folder_MC_ORA,"/", 'ClusterGoSummary_', PathwayName, '-', g_save, '.txt', sep="")#Export the ORA results as .csv
        write.table(clusterGoSummary,txtORA , col.names = TRUE, row.names = FALSE) #Export the ORA results as txt
      }
      }

    # Make list of DFs
    if(!(dim(clusterGoSummary)[1] == 0)){ #Check if the data frame is not empty
      # Add it to the list
      df_list[[g]] <- clusterGoSummary
    }
  }

  #Save excel file with output sheets:
  if(Save_as_Results == "xlsx"){
    if(PathwayName ==""){
      writexl::write_xlsx(df_list, paste(Results_folder_MC_ORA, "/ORASummary.xlsx", sep = ""), col_names = TRUE)#Export the ORA results as .xlxs
      } else{
        writexl::write_xlsx(df_list, paste(Results_folder_MC_ORA, "/ORASummary_", PathwayName, '.xlsx', sep = ""), col_names = TRUE)#Export the ORA results as .xlxs
      }
  }

  #return <- clusterGoSummary
  invisible(return(df_list))
}





#################################
### ### ### DM_ORA ### ### ###
#################################
#' DM_ORA
#'
#' Uses enricher to run ORA on the differential metabolites (DM) using a pathway list
#'
#' @param Input_data Input DF with column "t.val" and column "Metabolite", which needs to match the identifier type in column "Metabolite" in the PathwayFile.
#' @param PathwayFile DF that must include column "term" with the pathway name, column "Metabolite" with the Metabolite name or ID and column "Description" with pathway description that will be depicted on the plots.
#' @param PathwayName \emph{Optional: } Name of the pathway list used \strong{default: ""}
#' @param minGSSize \emph{Optional: } minimum group size in ORA \strong{default: 10}
#' @param maxGSSize \emph{Optional: } maximum group size in ORA \strong{default: 1000}
#' @param Save_as_Results \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt" \strong{default: "csv"}
#' @param pCutoff \emph{Optional: } p-adjusted value cutoff from ORA results. \strong{default: 0.05}
#' @param PercentageCutoff \emph{Optional: } Percentage cutoff of metabolites that should be considered for ORA. \strong{default: 10}
#' @param Folder_Name {Optional:} String which is added to the resulting folder name \strong(Default = NULL)
#'
#' @return Saves results as individual .csv files.
#' @export

DM_ORA <- function(Input_data,
                   PathwayFile,
                   PathwayName="",
                   minGSSize=10,
                   maxGSSize=1000 ,
                   pCutoff=0.05,
                   PercentageCutoff=10,
                   Save_as_Results="csv",
                   Folder_Name = NULL
){
  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse","clusterProfiler", "ggupset")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  suppressMessages(library(tidyverse))
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0){
    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(new.packages)
  }

  ## ------------ Check Input files ----------- ##
  # 1. The input data:
  if("Metabolite" %in% names(Input_data)==FALSE){
    stop("Check input. Input_data must contain a column named `Metabolite` including the metabolite names.")
  }
  if(length(Input_data[duplicated(Input_data$Metabolite), "Metabolite"]) > 0){
    doublons <- as.character(Input_data[duplicated(Input_data$Metabolite), "Metabolite"])#number of duplications
    Input_data <-Input_data[!duplicated(Input_data$Metabolite),]#remove duplications
    warning("Input_data contained duplicates based on Metabolite! Dropping duplicate IDs and kept only the first entry. You had ", length(doublons), " duplicates. Note that you should do this before running VizVolcano.")
  }

  # 2. PathwayFile and name
  if("Metabolite" %in% names(PathwayFile)==FALSE){
    stop("Check input. PathwayFile must contain a column named `Metabolite` including the metabolite names or IDs.")
  } else if("Metabolite" %in% names(PathwayFile)==TRUE){
    PathwayFile<-PathwayFile%>%
      dplyr::rename("gene"="Metabolite")
  }
  if("term" %in% names(PathwayFile)==FALSE){
    stop("Check input. PathwayFile must contain a column named `term` with pathway names.")
  }
  if("Description" %in% names(PathwayFile)==FALSE){
    stop("Check input. PathwayFile must contain a column named `Description` with description of the pathway.")
  }
  if(is.character(PathwayName)==FALSE){
    stop("Check input. PathwayName must be a character of syntax 'example'.")
  }

  # 3. General Settings
  if(is.numeric(minGSSize)== FALSE){
    stop("Check input. The selected minGSSize value should be numeric.")
  }
  if(is.numeric(maxGSSize)== FALSE){
    stop("Check input. The selected maxGSSize value should be numeric.")
  }
  Save_as_Results_options <- c("txt","csv", "xlsx" )
  if(Save_as_Results %in% Save_as_Results_options == FALSE){
    stop("Check input. The selected Save_as_Results option is not valid. Please select one of the folowwing: ",paste(Save_as_Results_options,collapse = ", "),"." )
  }
  if( is.numeric(pCutoff)== FALSE | pCutoff > 1 |  pCutoff < 0){
    stop("Check input. The selected Plot_pCutoff value should be numeric and between 0 and 1.")
  }
  if( is.numeric(PercentageCutoff)== FALSE  | PercentageCutoff > 100 | PercentageCutoff < 0){
    stop("Check input. The selected PercentageCutoff value should be numeric and between 0 and 100.")
  }

  ## ------------ Create Output folders ----------- ##
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
  Results_folder_DM_ORA = file.path(Results_folder, "DM_ORA")  # This searches for a folder within the "Results" folder in the current working directory and if its not found it creates one
  if (!dir.exists( Results_folder_DM_ORA)) {dir.create( Results_folder_DM_ORA)}  # check and create folder

  ############################################################################################################
  ## ------------ Load the data and check ----------- ##
  #Select universe
  allMetabolites <- as.character(Input_data$Metabolite)

  #select top changed metabolites (Up and down together)
  #check if the metabolites are significantly changed.
  value <- PercentageCutoff/100

  allMetabolites_DF <- Input_data[order(Input_data$t.val),]# rank by t.val
  selectMetabolites_DF <- allMetabolites_DF[c(1:(ceiling(value * nrow(allMetabolites_DF))),(nrow(allMetabolites_DF)-(ceiling(value * nrow(allMetabolites_DF)))):(nrow(allMetabolites_DF))),]
  selectMetabolites_DF$`Top/Bottom`<- "TRUE"
  selectMetabolites_DF <-merge(allMetabolites_DF,selectMetabolites_DF[,c("Metabolite", "Top/Bottom")], by="Metabolite", all.x=TRUE)

  InputSelection <- selectMetabolites_DF%>%
    mutate(`Top/Bottom_Percentage` = case_when(`Top/Bottom`==TRUE ~ 'TRUE',
                                 TRUE ~ 'FALSE'))%>%
    mutate(Significant = case_when(p.adj <= pCutoff ~ 'TRUE',
                                  TRUE ~ 'FALSE'))%>%
    mutate(Cluster_ChangedMetabolites = case_when(Significant==TRUE & `Top/Bottom_Percentage`==TRUE ~ 'TRUE',
                                   TRUE ~ 'FALSE'))
  InputSelection$`Top/Bottom` <- NULL #remove column as its not needed for output

  selectMetabolites <- InputSelection%>%
    subset(Cluster_ChangedMetabolites==TRUE)
  selectMetabolites <-as.character(selectMetabolites$Metabolite)

  #Load Pathways
  Pathway <- PathwayFile
  Term2gene <- Pathway[,c("term", "gene")]# term and MetaboliteID (MetaboliteID= gene as syntax required for enricher)
  term2name <- Pathway[,c("term", "Description")]# term and description

  #Add the number of genes present in each pathway
  Pathway$Count <- 1
  Pathway_Mean <- aggregate(Pathway$Count, by=list(term=Pathway$term), FUN=sum)
  names(Pathway_Mean)[names(Pathway_Mean) == "x"] <- "Metabolites_in_Pathway"
  Pathway <- merge(x= Pathway[,-4], y=Pathway_Mean,by="term", all.x=TRUE)

  ## ------------ Run ----------- ##
  #Run ORA
  clusterGo <- clusterProfiler::enricher(gene=selectMetabolites,
                                           pvalueCutoff = 1,
                                           pAdjustMethod = "BH",
                                           universe = allMetabolites,
                                           minGSSize=minGSSize,
                                           maxGSSize=maxGSSize,
                                           qvalueCutoff = 1,
                                           TERM2GENE=Term2gene ,
                                           TERM2NAME = term2name)
  clusterGoSummary <- data.frame(clusterGo)

  #Save file and plots
  if(!(dim(clusterGoSummary)[1] == 0)){
      #Add pathway information % of genes in pathway detected)
      clusterGoSummary <- merge(x= clusterGoSummary[,-2], y=Pathway[,-2],by.x="ID",by.y="term", all=TRUE)
      clusterGoSummary$Count[is.na(clusterGoSummary$Count)] <- 0
      clusterGoSummary$Percentage_of_Pathway_detected <-round(((clusterGoSummary$Count/clusterGoSummary$Metabolites_in_Pathway)*100),digits=2)
      clusterGoSummary <- clusterGoSummary[!duplicated(clusterGoSummary$ID),]
      clusterGoSummary <- clusterGoSummary[order(clusterGoSummary$p.adjust),]
      clusterGoSummary <- clusterGoSummary[,c(1,9,2:8, 10:11)]%>%
        dplyr::rename("MetaboliteIDs_in_pathway"="geneID")
  }else{
    stop("None of the Input_data Metabolites were present in any terms of the PathwayFile. Hence the ClusterGoSummary ouput will be empty. Please check that the metabolite IDs match the pathway IDs.")
    }

  #Save file
  if (Save_as_Results == "xlsx"){
    if(PathwayName ==""){
      ORA_output_list <- list(InputSelection = InputSelection , ClusterGoSummary = clusterGoSummary)
      writexl::write_xlsx(ORA_output_list, paste(Results_folder_DM_ORA, "/ORASummary.xlsx", sep = ""), col_names = TRUE)#Export the ORA results as .xlxs
      }else{
        ORA_output_list <- list(InputSelection = InputSelection , ClusterGoSummary = clusterGoSummary)
        writexl::write_xlsx(ORA_output_list, paste(Results_folder_DM_ORA, "/ORASummary_", PathwayName, '.xlsx', sep = ""), col_names = TRUE)#Export the ORA results as .xlxs
        }
    }else if (Save_as_Results == "csv"){
      if(PathwayName ==""){
        csvORA <- paste(Results_folder_DM_ORA,"/", 'ClusterGoSummary',',csv', sep="")
        csvInput <- paste(Results_folder_DM_ORA,"/", 'InputSelection',',csv', sep="")

        write_csv(clusterGoSummary,csvORA)#Export the ORA results as .csv
        write_csv(InputSelection,csvInput)#Export Input with selection columns as .csv
        }else{
          csvORA <- paste(Results_folder_DM_ORA,"/", 'ClusterGoSummary_', PathwayName, '.csv', sep="")#Export the ORA results as .csv
          csvInput <- paste(Results_folder_DM_ORA,"/", 'InputSelection_', PathwayName, '.csv', sep="")#Export the ORA results as .csv

          write_csv(clusterGoSummary,csvORA)
          write_csv(InputSelection,csvInput)
          }
      }else if (Save_as_Results == "txt"){
        if(PathwayName ==""){
          txtORA <-  paste(Results_folder_DM_ORA,"/", 'ClusterGoSummary','.txt', sep="")
          txtInput <-  paste(Results_folder_DM_ORA,"/", 'InputSelection','.txt', sep="")

          write.table(clusterGoSummary,txtORA, col.names = TRUE, row.names = FALSE) #Export the ORA results as txt
          write.table(InputSelection,txtInput, col.names = TRUE, row.names = FALSE) #Export the ORA results as txt
          }else{
            txtORA <-  paste(Results_folder_DM_ORA,"/", 'ClusterGoSummary_', PathwayName, '.txt', sep="")#Export the ORA results as .csv
            txtInput <-  paste(Results_folder_DM_ORA,"/", 'InputSelection_', PathwayName, '.txt', sep="")#Export the ORA results as .csv

            write.table(clusterGoSummary,txtORA , col.names = TRUE, row.names = FALSE) #Export the ORA results as txt
            write.table(InputSelection,txtInput , col.names = TRUE, row.names = FALSE) #Export the ORA results as txt
          }
      }
  #return list of DFs
  ORA_output_list <- list(InputSelection = InputSelection , ClusterGoSummary = clusterGoSummary)
  invisible(return(ORA_output_list))
  }


