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
### ### ### ClusterORA ### ### ###
#################################
#' ClusterORA
#'
#' Uses enricher to run ORA on each of the metabolite cluster from any of the MCA functions using a pathway list
#'
#' @param InputData DF with metabolite names/metabolite IDs as row names. Metabolite names/IDs need to match the identifier type (e.g. HMDB IDs) in the PathwayFile.
#' @param SettingsInfo \emph{Optional: } Pass ColumnName of the column including the cluster names that ORA should be performed on (=ClusterColumn). BackgroundColumn passes the column name needed if RemoveBackground=TRUE. Also pass ColumnName for PathwayFile including term and feature names. (ClusterColumn= ColumnName InputData, BackgroundColumn = ColumnName InputData, PathwayTerm= ColumnName PathwayFile, PathwayFeature= ColumnName PathwayFile) \strong{c(FeatureName="Metabolite", ClusterColumn="RG2_Significant", BackgroundColumn="BG_Method", PathwayTerm= "term", PathwayFeature= "Metabolite")}
#' @param RemoveBackground \emph{Optional: } If TRUE, column BackgroundColumn  name needs to be in SettingsInfo, which includes TRUE/FALSE for each metabolite to fall into background based on the chosen Background method for e.g. MCA_2Cond are removed from the universe. \strong{default: TRUE}
#' @param PathwayFile DF that must include column "term" with the pathway name, column "Feature" with the Metabolite name or ID and column "Description" with pathway description.
#' @param PathwayName \emph{Optional: } Name of the pathway list used \strong{default: ""}
#' @param minGSSize \emph{Optional: } minimum group size in ORA \strong{default: 10}
#' @param maxGSSize \emph{Optional: } maximum group size in ORA \strong{default: 1000}
#' @param SaveAs_Table \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt" \strong{default: "csv"}
#' @param FolderPath \emph{Optional:} Path to the folder the results should be saved at. \strong{default: NULL}
#'
#' @return Saves results as individual .csv files.
#' @export

ClusterORA <- function(InputData,
                       SettingsInfo=c(ClusterColumn="RG2_Significant", BackgroundColumn="BG_Method", PathwayTerm= "term", PathwayFeature= "Metabolite"),
                       RemoveBackground=TRUE,
                       PathwayFile,
                       PathwayName="",
                       minGSSize=10,
                       maxGSSize=1000 ,
                       SaveAs_Table= "csv",
                       FolderPath = NULL){
   ## ------------ Check Input files ----------- ##
  Pathways <- MetaProViz:::CheckInput_ORA(InputData=InputData,
                                          SettingsInfo=SettingsInfo,
                                          RemoveBackground=RemoveBackground,
                                          PathwayFile=PathwayFile,
                                          PathwayName=PathwayName,
                                          minGSSize=minGSSize,
                                          maxGSSize=maxGSSize,
                                          SaveAs_Table=SaveAs_Table,
                                          pCutoff=NULL,
                                          PercentageCutoff=NULL)

  ## ------------ Create Results output folder ----------- ##
  if(is.null(SaveAs_Table)==FALSE){
    Folder <- MetaProViz:::SavePath(FolderName= "ClusterORA",
                                    FolderPath=FolderPath)
    }

  ############################################################################################################
  ## ------------ Prepare Data ----------- ##
  # open the data
  if(RemoveBackground==TRUE){
    df <- subset(InputData, !InputData[[SettingsInfo[["BackgroundColumn"]]]] == "FALSE")%>%
      rownames_to_column("Metabolite")
  } else{
    df <- InputData%>%
      rownames_to_column("Metabolite")
  }

  #Select universe
  allMetabolites <- as.character(df$Metabolite)

  #Select clusters
  grps_labels <- unlist(unique(df[[SettingsInfo[["ClusterColumn"]]]]))

  #Load Pathways
  Pathway <- Pathways
  Term2gene <- Pathway[,c("term", "gene")]# term and MetaboliteID (MetaboliteID= gene as syntax required for enricher)
  term2name <- Pathway[,c("term", "Description")]# term and description

  #Add the number of genes present in each pathway
  Pathway$Count <- 1
  Pathway_Mean <- aggregate(Pathway$Count, by=list(term=Pathway$term), FUN=sum)
  names(Pathway_Mean)[names(Pathway_Mean) == "x"] <- "Metabolites_in_Pathway"
  Pathway <- merge(x= Pathway, y=Pathway_Mean,by="term", all.x=TRUE)
  Pathway$Count <- NULL

  ## ------------ Run ----------- ##
  df_list = list()# Make an empty list to store the created DFs
  clusterGo_list = list()
  #Run ORA
  for(g in grps_labels){
    grpMetabolites <- subset(df, df[[SettingsInfo[["ClusterColumn"]]]] == g)
    message("Number of metabolites in cluster `",g, "`: ", nrow(grpMetabolites), sep="")

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
    clusterGo_list[[g]]<- clusterGo
    if(!(dim(clusterGoSummary)[1] == 0)){
      #Add pathway information (% of genes in pathway detected)
      clusterGoSummary <- merge(x= clusterGoSummary%>% select(-Description), y=Pathway%>% select(term, Metabolites_in_Pathway), by.x="ID",by.y="term", all=TRUE)
      clusterGoSummary$Count[is.na(clusterGoSummary$Count)] <- 0
      clusterGoSummary$Percentage_of_Pathway_detected <-round(((clusterGoSummary$Count/clusterGoSummary$Metabolites_in_Pathway)*100),digits=2)
      clusterGoSummary <- clusterGoSummary[!duplicated(clusterGoSummary$ID),]
      clusterGoSummary <- clusterGoSummary[order(clusterGoSummary$p.adjust),]
      clusterGoSummary <- clusterGoSummary%>%
        dplyr::rename("Metabolites_in_pathway"="geneID")

      g_save <- gsub("/", "-", g)
      df_list[[g_save]] <- clusterGoSummary
    }else{
      message("None of the Input_data Metabolites of the cluster ", g ," were present in any terms of the PathwayFile. Hence the ClusterGoSummary ouput will be empty for this cluster. Please check that the metabolite IDs match the pathway IDs.")
    }
  }
  #Save files
  suppressMessages(suppressWarnings(
    MetaProViz:::SaveRes(InputList_DF=df_list,
                         InputList_Plot= NULL,
                         SaveAs_Table=SaveAs_Table,
                         SaveAs_Plot=NULL,
                         FolderPath= Folder,
                         FileName= paste("ClusterGoSummary",PathwayName, sep="_"),
                         CoRe=FALSE,
                         PrintPlot=FALSE)))

  #return <- clusterGoSummary
  ORA_Output <- list("DF"= df_list, "ClusterGo"=clusterGo_list)

  invisible(return(ORA_Output))
}




#################################
### ### ### StandardORA ### ### ###
#################################
#' StandardORA
#'
#' Uses enricher to run ORA on the differential metabolites (DM) using a pathway list
#'
#' @param InputData DF with metabolite names/metabolite IDs as row names. Metabolite names/IDs need to match the identifier type (e.g. HMDB IDs) in the PathwayFile.
#' @param SettingsInfo \emph{Optional: } Pass ColumnName of the column including parameters to use for pCutoff and PercentageCutoff. Also pass ColumnName for PathwayFile including term and feature names. (pvalColumn = ColumnName InputData, PercentageColumn= ColumnName InputData, PathwayTerm= ColumnName PathwayFile, PathwayFeature= ColumnName PathwayFile) \strong{c(pvalColumn="p.adj", PercentageColumn="t.val", PathwayTerm= "term", PathwayFeature= "Metabolite")}
#' @param pCutoff \emph{Optional: } p-adjusted value cutoff from ORA results. Must be a numeric value. \strong{default: 0.05}
#' @param PercentageCutoff \emph{Optional: } Percentage cutoff of metabolites that should be considered for ORA. Selects Top/Bottom % of selected PercentageColumn, usually t.val or Log2FC \strong{default: 10}
#' @param PathwayFile DF that must include column "term" with the pathway name, column "Metabolite" with the Metabolite name or ID and column "Description" with pathway description that will be depicted on the plots.
#' @param PathwayName \emph{Optional: } Name of the PathwayFile used \strong{default: ""}
#' @param minGSSize \emph{Optional: } minimum group size in ORA \strong{default: 10}
#' @param maxGSSize \emph{Optional: } maximum group size in ORA \strong{default: 1000}
#' @param Save_as_Results \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt" \strong{default: "csv"}
#' @param FolderPath \emph{Optional:} Path to the folder the results should be saved at. \strong{default: NULL}
#'
#' @return Saves results as individual .csv files.
#' @export

StandardORA <- function(InputData,
                        SettingsInfo=c(pvalColumn="p.adj", PercentageColumn="t.val", PathwayTerm= "term", PathwayFeature= "Metabolite"),
                        pCutoff=0.05,
                        PercentageCutoff=10,
                        PathwayFile,
                        PathwayName="",
                        minGSSize=10,
                        maxGSSize=1000 ,
                        SaveAs_Table="csv",
                        FolderPath = NULL

){
 ## ------------ Check Input files ----------- ##
  Pathways <- MetaProViz:::CheckInput_ORA(InputData=InputData,
                                          SettingsInfo=SettingsInfo,
                                          RemoveBackground=FALSE,
                                          PathwayFile=PathwayFile,
                                          PathwayName=PathwayName,
                                          minGSSize=minGSSize,
                                          maxGSSize=maxGSSize,
                                          SaveAs_Table=SaveAs_Table,
                                          pCutoff=pCutoff,
                                          PercentageCutoff=PercentageCutoff)

  ## ------------ Create Results output folder ----------- ##
  if(is.null(SaveAs_Table)==FALSE){
    Folder <- MetaProViz:::SavePath(FolderName= "ORA",
                                    FolderPath=FolderPath)
  }

  ############################################################################################################
  ## ------------ Load the data and check ----------- ##
  InputData<- InputData %>%
    rownames_to_column("Metabolite")

  #Select universe
  allMetabolites <- as.character(InputData$Metabolite)

  #select top changed metabolites (Up and down together)
  #check if the metabolites are significantly changed.
  value <- PercentageCutoff/100

  allMetabolites_DF <- InputData[order(InputData[[SettingsInfo[["PercentageColumn"]]]]),]# rank by t.val
  selectMetabolites_DF <- allMetabolites_DF[c(1:(ceiling(value * nrow(allMetabolites_DF))),(nrow(allMetabolites_DF)-(ceiling(value * nrow(allMetabolites_DF)))):(nrow(allMetabolites_DF))),]
  selectMetabolites_DF$`Top/Bottom`<- "TRUE"
  selectMetabolites_DF <-merge(allMetabolites_DF,selectMetabolites_DF[,c("Metabolite", "Top/Bottom")], by="Metabolite", all.x=TRUE)

  InputSelection <- selectMetabolites_DF%>%
    mutate(`Top/Bottom_Percentage` = case_when(`Top/Bottom`==TRUE ~ 'TRUE',
                                 TRUE ~ 'FALSE'))%>%
    mutate(Significant = case_when(get(SettingsInfo[["pvalColumn"]]) <= pCutoff ~ 'TRUE',
                                  TRUE ~ 'FALSE'))%>%
    mutate(Cluster_ChangedMetabolites = case_when(Significant==TRUE & `Top/Bottom_Percentage`==TRUE ~ 'TRUE',
                                   TRUE ~ 'FALSE'))
  InputSelection$`Top/Bottom` <- NULL #remove column as its not needed for output

  selectMetabolites <- InputSelection%>%
    subset(Cluster_ChangedMetabolites==TRUE)
  selectMetabolites <-as.character(selectMetabolites$Metabolite)

  #Load Pathways
  Pathway <- Pathways
  Term2gene <- Pathway[,c("term", "gene")]# term and MetaboliteID (MetaboliteID= gene as syntax required for enricher)
  term2name <- Pathway[,c("term", "Description")]# term and description

  #Add the number of genes present in each pathway
  Pathway$Count <- 1
  Pathway_Mean <- aggregate(Pathway$Count, by=list(term=Pathway$term), FUN=sum)
  names(Pathway_Mean)[names(Pathway_Mean) == "x"] <- "Metabolites_in_Pathway"
  Pathway <- merge(x= Pathway, y=Pathway_Mean,by="term", all.x=TRUE)
  Pathway$Count <- NULL

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

  #Make DF:
  if(!(dim(clusterGoSummary)[1] == 0)){
      #Add pathway information % of genes in pathway detected)
      clusterGoSummary <- merge(x= clusterGoSummary%>% select(-Description), y=Pathway%>% select(term, Metabolites_in_Pathway),by.x="ID",by.y="term", all=TRUE)
      clusterGoSummary$Count[is.na(clusterGoSummary$Count)] <- 0
      clusterGoSummary$Percentage_of_Pathway_detected <-round(((clusterGoSummary$Count/clusterGoSummary$Metabolites_in_Pathway)*100),digits=2)
      clusterGoSummary <- clusterGoSummary[!duplicated(clusterGoSummary$ID),]
      clusterGoSummary <- clusterGoSummary[order(clusterGoSummary$p.adjust),]
      clusterGoSummary <- clusterGoSummary%>%
        dplyr::rename("Metabolites_in_pathway"="geneID")
  }else{
    stop("None of the Input_data Metabolites were present in any terms of the PathwayFile. Hence the ClusterGoSummary ouput will be empty. Please check that the metabolite IDs match the pathway IDs.")
    }

  #Return and save list of DFs
  ORA_output_list <- list("InputSelection" = InputSelection , "ClusterGoSummary" = clusterGoSummary)

  #save:
  suppressMessages(suppressWarnings(
    MetaProViz:::SaveRes(InputList_DF=ORA_output_list,
                         InputList_Plot= NULL,
                         SaveAs_Table=SaveAs_Table,
                         SaveAs_Plot=NULL,
                         FolderPath= Folder,
                         FileName= paste(PathwayName),
                         CoRe=FALSE,
                         PrintPlot=FALSE)))

  #Return
  ORA_output_list <- c( ORA_output_list, list("ClusterGo"=clusterGo))

  invisible(return(ORA_output_list))
  }


################################################################################################
### ### ### ORA helper function: Internal Function to check function input ### ### ###
################################################################################################

#' Check input parameters
#'
#' @param InputData Passed to main function MetaProViz::PreProcessing()
#' @param SettingsInfo Passed to main function MetaProViz::PreProcessing()
#'
#' @keywords Input check
#' @noRd
#'
#'

CheckInput_ORA <- function(InputData,
                           SettingsInfo,
                           RemoveBackground,
                           PathwayFile,
                           PathwayName,
                           minGSSize,
                           maxGSSize,
                           SaveAs_Table,
                           pCutoff,
                           PercentageCutoff
                           ){
  # 1. The input data:
  if(class(InputData) != "data.frame"){
    stop("InputData should be a data.frame. It's currently a ", paste(class(InputData), ".",sep = ""))
  }
  if(any(duplicated(row.names(InputData)))==TRUE){
    stop("Duplicated row.names of InputData, whilst row.names must be unique")
  }

  # 2. Settings Columns:
  if(is.vector(SettingsInfo)==FALSE & is.null(SettingsInfo)==FALSE){
    stop("SettingsInfo should be NULL or a vector. It's currently a ", paste(class(SettingsInfo), ".", sep = ""))
  }

  if(is.null(SettingsInfo)==FALSE){
    #"ClusterColumn"
    if("ClusterColumn" %in% names(SettingsInfo)){
      if(SettingsInfo[["ClusterColumn"]] %in% colnames(InputData)== FALSE){
        stop("The ", SettingsInfo[["ClusterColumn"]], " column selected as ClusterColumn in SettingsInfo was not found in InputData. Please check your input.")
      }
    }

    #"BackgroundColumn"
    if("BackgroundColumn" %in% names(SettingsInfo)){
      if(SettingsInfo[["BackgroundColumn"]] %in% colnames(InputData)== FALSE){
        stop("The ", SettingsInfo[["BackgroundColumn"]], " column selected as BackgroundColumn in SettingsInfo was not found in InputData. Please check your input.")
      }
    }

    #"pvalColumn"
    if("pvalColumn" %in% names(SettingsInfo)){
      if(SettingsInfo[["pvalColumn"]] %in% colnames(InputData)== FALSE){
        stop("The ", SettingsInfo[["pvalColumn"]], " column selected as pvalColumn in SettingsInfo was not found in InputData. Please check your input.")
      }
    }

    #"PercentageColumn"
    if("PercentageColumn" %in% names(SettingsInfo)){
      if(SettingsInfo[["PercentageColumn"]] %in% colnames(InputData)== FALSE){
        stop("The ", SettingsInfo[["PercentageColumn"]], " column selected as PercentageColumn in SettingsInfo was not found in InputData. Please check your input.")
      }
    }


    #"PathwayTerm"
    if("PathwayTerm" %in% names(SettingsInfo)){
      if(SettingsInfo[["PathwayTerm"]] %in% colnames(PathwayFile)== FALSE){
        stop("The ", SettingsInfo[["PathwayTerm"]], " column selected as PathwayTerm in SettingsInfo was not found in PathwayFile. Please check your input.")
      }else{
        PathwayFile <- PathwayFile%>%
          dplyr::rename("term"=SettingsInfo[["PathwayTerm"]])
        PathwayFile$Description <- PathwayFile$term
      }
    }else{
      stop("SettingsInfo must provide the column name for PathwayTerm in PathwayFile")
    }

    # PathwayFeature
    if("PathwayFeature" %in% names(SettingsInfo)){
      if(SettingsInfo[["PathwayFeature"]] %in% colnames(PathwayFile)== FALSE){
        stop("The ", SettingsInfo[["PathwayFeature"]], " column selected as PathwayFeature in SettingsInfo was not found in PathwayFile. Please check your input.")
      }else{
        PathwayFile <- PathwayFile%>%
          dplyr::rename("gene"=SettingsInfo[["PathwayFeature"]])
      }
    }else{
      stop("SettingsInfo must provide the column name for PathwayFeature in PathwayFile")
    }

  }else{
    stop("you must provide SettingsInfo.")
  }

  # 3. General Settings
  if(is.character(PathwayName)==FALSE){
    stop("Check input. PathwayName must be a character of syntax 'example'.")
  }

  if(is.logical(RemoveBackground) == FALSE){
    stop("Check input. RemoveBackground value should be either =TRUE or = FALSE.")
  }

  if(is.numeric(minGSSize)== FALSE){
    stop("Check input. The selected minGSSize value should be numeric.")
  }

  if(is.numeric(maxGSSize)== FALSE){
    stop("Check input. The selected maxGSSize value should be numeric.")
  }

  SaveAs_Table_options <- c("txt","csv", "xlsx", "RData")#RData = SummarizedExperiment (?)
  if(is.null(SaveAs_Table)==FALSE){
    if((SaveAs_Table %in% SaveAs_Table_options == FALSE)| (is.null(SaveAs_Table)==TRUE)){
      stop("Check input. The selected SaveAs_Table option is not valid. Please select one of the folowwing: ",paste(SaveAs_Table_options,collapse = ", "),"." )
    }
  }

  if(is.null(pCutoff)== FALSE){
    if(is.numeric(pCutoff)== FALSE | pCutoff > 1 |  pCutoff < 0){
    stop("Check input. The selected Plot_pCutoff value should be numeric and between 0 and 1.")
  }
  }

  if(is.null(PercentageCutoff)== FALSE){
    if( is.numeric(PercentageCutoff)== FALSE  | PercentageCutoff > 100 | PercentageCutoff < 0){
    stop("Check input. The selected PercentageCutoff value should be numeric and between 0 and 100.")
  }
  }


  ## -------- Return Pathways ---------##
  return(invisible(PathwayFile))
}
