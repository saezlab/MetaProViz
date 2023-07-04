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
#' @param RemoveBackgroundGenes\emph{Optional: } If TRUE, column "BG_Method" needs to be in Input_data, which includes TRUE/FALSE for each metabolite to fall into background based on the chosen Background method for e.g. MCA_2Cond are removed from the universe. \strong{default: "TRUE"}
#' @param PathwayFile DF that must include column "term" with the pathway name, column "Metabolite" with the Metabolite name or ID and column "Description" with pathway description that will be depicted on the plots.
#' @param PathwayName \emph{Optional: } Name of the pathway list used \strong{default: ""}
#' @param minGSSize \emph{Optional: } minimum group size in ORA \strong{default: 10}
#' @param maxGSSize \emph{Optional: } maximum group size in ORA \strong{default: 1000}
#' @param Save_as \emph{Optional: } If Save_as=NULL no plots will be saved. Otherwise, file types for the figures are: "svg", "pdf", "png" \strong{default: "pdf"}
#' @param pCutoff \emph{Optional: } p-adjusted value cutoff from ORA results that should be plotted. If Save_as=NULL, this is ignored. \strong{default: 0.2}
#' @param PercentageCutoff \emph{Optional: } Percentage cutoff of metabolites that are detected of a pathway from ORA results that should be plotted. If Save_as=NULL, this is ignored. \strong{default: 10}
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
                   Save_as="pdf",
                   pCutoff=0.2,
                   PercentageCutoff=10
                   ){
  ## ------------ Setup and installs ----------- ##
  packages <- c("tidyverse","clusterProfiler", "enrichplot", "ggupset")
  install.packages(setdiff(packages, rownames(installed.packages())))
  suppressMessages(library(tidyverse))

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
  Save_as_options <- c("svg","pdf")
  if(Save_as %in% Save_as_options == FALSE | is.null(Save_as)==FALSE){
    stop("Check input. The selected Save_as option is not valid. Please set Save_as=NULL or select one of the following: ",paste(Save_as_options,collapse = ", "),"." )
  }
  if( is.numeric(pCutoff)== FALSE | pCutoff > 1 |  pCutoff < 0){
    stop("Check input. The selected Plot_pCutoff value should be numeric and between 0 and 1.")
  }
  if( is.numeric(PercentageCutoff)== FALSE  | PercentageCutoff > 100 | PercentageCutoff < 0){
    stop("Check input. The selected PercentageCutoff value should be numeric and between 0 and 100.")
  }

  ## ------------ Create Output folders ----------- ##
  name <- paste0("MetaProViz_Results_",Sys.Date())
  WorkD <- getwd()
  Results_folder <- file.path(WorkD, name)
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)} # Make Results folder
  Results_folder_MC_ORA = file.path(Results_folder, "MC_ORA")  # This searches for a folder within the "Results" folder in the current working directory and if its not found it creates one
  if (!dir.exists( Results_folder_MC_ORA)) {dir.create( Results_folder_MC_ORA)}  # check and create folder

  ############################################################################################################
  ## ------------ Run ----------- ##
  # open the data
  if(RemoveBackground=="TRUE"){
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

  #Run ORA
  for(g in grps_labels) {
    grpMetabolites <- subset(df, df[[MetabCluster_lab]] == g)
    print(g)
    print(dim(grpMetabolites))

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
    if (!(dim(clusterGoSummary)[1] == 0)){
      #Add pathway information (% of genes in pathway detected)
      clusterGoSummary <- merge(x= clusterGoSummary[,-2], y=Pathway[,-2],by.x="ID",by.y="term", all=TRUE)
      clusterGoSummary$Count[is.na(clusterGoSummary$Count)] <- 0
      clusterGoSummary$Percentage_of_Pathway_detected <-round(((clusterGoSummary$Count/clusterGoSummary$Metabolites_in_Pathway)*100),digits=2)
      clusterGoSummary <- clusterGoSummary[!duplicated(clusterGoSummary$ID),]
      clusterGoSummary <- clusterGoSummary[order(clusterGoSummary$p.adjust),]
      clusterGoSummary <- clusterGoSummary[,c(1,9,2:8, 10:11)]%>%
        dplyr::rename("MetaboliteIDs_in_pathway"="geneID")

      #Safe file
      if(PathwayName ==""){
        write_csv(clusterGoSummary, paste(Results_folder_MC_ORA,"/", 'ClusterGoSummary_', g, '.csv', sep=""))#Export the ORA results as .csv
        } else{
          write_csv(clusterGoSummary, paste(Results_folder_MC_ORA,"/", 'ClusterGoSummary_', PathwayName, '-', g, '.csv', sep=""))#Export the ORA results as .csv
        }

      #Make Selection of terms that should be displayed on the plots
      clusterGoSummary_Select <- clusterGoSummary %>%
        subset(p.adjust <= pCutoff & Percentage_of_Pathway_detected >= PercentageCutoff)
      rownames(clusterGoSummary_Select)<-clusterGoSummary_Select$ID
      clusterGoSummary_Select<-  clusterGoSummary_Select%>%
        dplyr::rename("geneID"="MetaboliteIDs_in_pathway")

      ## ------------ Plots ----------- ##
      if(is.null(Save_as)==FALSE){
        if (!(dim(clusterGoSummary_Select)[1] == 0)) {#exclude df's that have no observations
        clusterGo@result <- clusterGoSummary_Select[,1:9]
        #1. Dotplot:
        Dotplot <-  enrichplot::dotplot(clusterGo, showCategory=nrow(clusterGoSummary_Select)) +
          ggtitle(paste("Dotplot:", g, PathwayName, sep=" "))
        if(PathwayName ==""){
          ggsave(file=paste(Results_folder_MC_ORA,"/", "Dotplot_", g, ".", Save_as, sep=""), plot=Dotplot, width=10, height=8)
        } else{
          ggsave(file=paste(Results_folder_MC_ORA,"/", "Dotplot_", g,"_" ,PathwayName, ".", Save_as, sep=""), plot=Dotplot, width=10, height=8)
        }
        plot(Dotplot)

        #2. Emapplot
        x2 <- enrichplot::pairwise_termsim(clusterGo)
        Emapplot <-  enrichplot::emapplot(x2, pie_scale=1.5, layout = "nicely", showCategory=nrow(clusterGoSummary_Select))+
          ggtitle(paste("Emapplot:", g, PathwayName, sep=" "))
        if(PathwayName ==""){
          ggsave(file=paste(Results_folder_MC_ORA,"/", "Emapplot_", g, ".", Save_as, sep=""), plot=Emapplot, width=10, height=8)
        } else{
          ggsave(file=paste(Results_folder_MC_ORA,"/", "Emapplot_", g,"_" ,PathwayName, ".", Save_as, sep=""), plot=Emapplot, width=10, height=8)
        }
        plot(Emapplot)

        #3. Upsetplot:
        UpsetPlot <-  enrichplot::upsetplot(clusterGo, showCategory=nrow(clusterGoSummary_Select))+
          ggtitle(paste("UpsetPlot", g, PathwayName, sep=" "))
        if(PathwayName ==""){
          ggsave(file=paste(Results_folder_MC_ORA,"/", "UpsetPlot_", g, ".", Save_as, sep=""), plot=UpsetPlot, width=10, height=8)
        } else{
          ggsave(file=paste(Results_folder_MC_ORA,"/", "UpsetPlot_", g,"_" ,PathwayName, ".", Save_as, sep=""), plot=UpsetPlot, width=10, height=8)
        }
        plot(UpsetPlot)
        }
      }
    }
  }
}

