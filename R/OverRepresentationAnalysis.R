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


#################################
### ### ### ClusterORA ### ### ###
#################################
#' ClusterORA
#'
#' Uses enricher to run ORA on each of the metabolite cluster from any of the MCA functions using a pathway list
#'
#' @param data DF with metabolite names/metabolite IDs as row names. Metabolite names/IDs need to match the identifier type (e.g. HMDB IDs) in the input_pathway.
#' @param metadata_info \emph{Optional: } Pass ColumnName of the column including the cluster names that ORA should be performed on (=ClusterColumn). BackgroundColumn passes the column name needed if remove_background=TRUE. Also pass ColumnName for input_pathway including term and feature names. (ClusterColumn= ColumnName data, BackgroundColumn = ColumnName data, PathwayTerm= ColumnName input_pathway, PathwayFeature= ColumnName input_pathway) \strong{c(FeatureName="Metabolite", ClusterColumn="RG2_Significant", BackgroundColumn="BG_method", PathwayTerm= "term", PathwayFeature= "Metabolite")}
#' @param remove_background \emph{Optional: } If TRUE, column BackgroundColumn  name needs to be in metadata_info, which includes TRUE/FALSE for each metabolite to fall into background based on the chosen Background method for e.g. MCA_2Cond are removed from the universe. \strong{default: TRUE}
#' @param input_pathway DF that must include column "term" with the pathway name, column "Feature" with the Metabolite name or ID and column "Description" with pathway description.
#' @param pathway_name \emph{Optional: } Name of the pathway list used \strong{default: ""}
#' @param min_gssize \emph{Optional: } minimum group size in ORA \strong{default: 10}
#' @param max_gssize \emph{Optional: } maximum group size in ORA \strong{default: 1000}
#' @param save_table \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt" \strong{default: "csv"}
#' @param path \emph{Optional:} Path to the folder the results should be saved at. \strong{default: NULL}
#'
#' @return Saves results as individual .csv files.
#'
#' @importFrom logger log_info
#' @importFrom dplyr rename select
#' @export
ClusterORA <- function(data,
                       metadata_info=c(ClusterColumn="RG2_Significant", BackgroundColumn="BG_method", PathwayTerm= "term", PathwayFeature= "Metabolite"),
                       remove_background=TRUE,
                       input_pathway,
                       pathway_name="",
                       min_gssize=10,
                       max_gssize=1000 ,
                       save_table= "csv",
                       path = NULL){
  ## ------------ Create log file ----------- ##
  MetaProViz_Init()


   ## ------------ Check Input files ----------- ##
  Pathways <- CheckInput_ORA(data=data,
                                          metadata_info=metadata_info,
                                          remove_background=remove_background,
                                          input_pathway=input_pathway,
                                          pathway_name=pathway_name,
                                          min_gssize=min_gssize,
                                          max_gssize=max_gssize,
                                          save_table=save_table,
                                          cutoff_stat=NULL,
                                          cutoff_percentage=NULL)

  ## ------------ Create Results output folder ----------- ##
  if(is.null(save_table)==FALSE){
    Folder <- SavePath(folder_name= "ClusterORA",
                                    path=path)
    }

  ############################################################################################################
  ## ------------ Prepare data ----------- ##
  # open the data
  if(remove_background==TRUE){
    df <- subset(data, !data[[metadata_info[["BackgroundColumn"]]]] == "FALSE")%>%
      rownames_to_column("Metabolite")
  } else{
    df <- data%>%
      rownames_to_column("Metabolite")
  }

  #Select universe
  allMetabolites <- as.character(df$Metabolite)

  #Select clusters
  grps_labels <- unlist(unique(df[[metadata_info[["ClusterColumn"]]]]))

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
    grpMetabolites <- subset(df, df[[metadata_info[["ClusterColumn"]]]] == g)
    log_info("Number of metabolites in cluster `", g, "`: ", nrow(grpMetabolites), sep="")

    clusterGo <- enricher_internal(
        gene=as.character(grpMetabolites$Metabolite),
        pvalueCutoff = 1L,
        pAdjustmethod = "BH",
        universe = allMetabolites,
        minGSSize=min_gssize,
        maxGSSize=max_gssize,
        qvalueCutoff = 1L,
        TERM2GENE=Term2gene ,
        TERM2NAME = term2name
    )
    clusterGosummary <- data.frame(clusterGo)
    clusterGo_list[[g]]<- clusterGo
    if(!(dim(clusterGosummary)[1] == 0)){
      #Add pathway information (% of genes in pathway detected)
      clusterGosummary <- merge(x= clusterGosummary%>% select(-Description), y=Pathway%>% select(term, Metabolites_in_Pathway), by.x="ID",by.y="term", all=TRUE)
      clusterGosummary$Count[is.na(clusterGosummary$Count)] <- 0
      clusterGosummary$percentage_of_Pathway_detected <-round(((clusterGosummary$Count/clusterGosummary$Metabolites_in_Pathway)*100),digits=2)
      clusterGosummary <- clusterGosummary[!duplicated(clusterGosummary$ID),]
      clusterGosummary <- clusterGosummary[order(clusterGosummary$p.adjust),]
      clusterGosummary <- clusterGosummary%>%
        dplyr::rename("Metabolites_in_pathway"="geneID")

      g_save <- gsub("/", "-", g)
      df_list[[g_save]] <- clusterGosummary
    }else{
      log_info("None of the Input_data Metabolites of the cluster ", g ," were present in any terms of the input_pathway. Hence the ClusterGosummary ouput will be empty for this cluster. Please check that the metabolite IDs match the pathway IDs.")
    }
  }
  #Save files
  suppressMessages(suppressWarnings(
    SaveRes(inputlist_df=df_list,
                         inputlist_plot= NULL,
                         save_table=save_table,
                         save_plot=NULL,
                         path= Folder,
                         FileName= paste("ClusterGosummary",pathway_name, sep="_"),
                         core=FALSE,
                         print_plot=FALSE)))

  #return <- clusterGosummary
  ORA_Output <- list("DF"= df_list, "ClusterGo"=clusterGo_list)

  invisible(return(ORA_Output))
}


###################################
### ### ### StandardORA ### ### ###
###################################

#' Uses enricher to run ORA on the differential metabolites (DM) using a pathway list
#'
#' @param data DF with metabolite names/metabolite IDs as row names. Metabolite names/IDs need to match the identifier type (e.g. HMDB IDs) in the input_pathway.
#' @param metadata_info \emph{Optional: } Pass ColumnName of the column including parameters to use for cutoff_stat and cutoff_percentage. Also pass ColumnName for input_pathway including term and feature names. (pvalColumn = ColumnName data, percentageColumn= ColumnName data, PathwayTerm= ColumnName input_pathway, PathwayFeature= ColumnName input_pathway) \strong{c(pvalColumn="p.adj", percentageColumn="t.val", PathwayTerm= "term", PathwayFeature= "Metabolite")}
#' @param cutoff_stat \emph{Optional: } p-adjusted value cutoff from ORA results. Must be a numeric value. \strong{default: 0.05}
#' @param cutoff_percentage \emph{Optional: } percentage cutoff of metabolites that should be considered for ORA. Selects top/Bottom % of selected percentageColumn, usually t.val or Log2FC \strong{default: 10}
#' @param input_pathway DF that must include column "term" with the pathway name, column "Metabolite" with the Metabolite name or ID and column "Description" with pathway description that will be depicted on the plots.
#' @param pathway_name \emph{Optional: } Name of the input_pathway used \strong{default: ""}
#' @param min_gssize \emph{Optional: } minimum group size in ORA \strong{default: 10}
#' @param max_gssize \emph{Optional: } maximum group size in ORA \strong{default: 1000}
#' @param save_table \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt" \strong{default: "csv"}
#' @param path \emph{Optional:} Path to the folder the results should be saved at. \strong{default: NULL}
#'
#' @return Saves results as individual .csv files.
#'
#' @export
#'
StandardORA <- function(data,
                        metadata_info=c(pvalColumn="p.adj", percentageColumn="t.val", PathwayTerm= "term", PathwayFeature= "Metabolite"),
                        cutoff_stat=0.05,
                        cutoff_percentage=10,
                        input_pathway,
                        pathway_name="",
                        min_gssize=10,
                        max_gssize=1000 ,
                        save_table="csv",
                        path = NULL

){
  ## ------------ Create log file ----------- ##
  MetaProViz_Init()

 ## ------------ Check Input files ----------- ##
  Pathways <- CheckInput_ORA(data=data,
                                          metadata_info=metadata_info,
                                          remove_background=FALSE,
                                          input_pathway=input_pathway,
                                          pathway_name=pathway_name,
                                          min_gssize=min_gssize,
                                          max_gssize=max_gssize,
                                          save_table=save_table,
                                          cutoff_stat=cutoff_stat,
                                          cutoff_percentage=cutoff_percentage)

  ## ------------ Create Results output folder ----------- ##
  if(is.null(save_table)==FALSE){
    Folder <- SavePath(folder_name= "ORA",
                                    path=path)
  }

  ############################################################################################################
  ## ------------ Load the data and check ----------- ##
  data<- data %>%
    rownames_to_column("Metabolite")

  #Select universe
  allMetabolites <- as.character(data$Metabolite)

  #select top changed metabolites (Up and down together)
  #check if the metabolites are significantly changed.
  value <- cutoff_percentage/100

  allMetabolites_DF <- data[order(data[[metadata_info[["percentageColumn"]]]]),]# rank by t.val
  selectMetabolites_DF <- allMetabolites_DF[c(1:(ceiling(value * nrow(allMetabolites_DF))),(nrow(allMetabolites_DF)-(ceiling(value * nrow(allMetabolites_DF)))):(nrow(allMetabolites_DF))),]
  selectMetabolites_DF$`top/Bottom`<- "TRUE"
  selectMetabolites_DF <-merge(allMetabolites_DF,selectMetabolites_DF[,c("Metabolite", "top/Bottom")], by="Metabolite", all.x=TRUE)

  InputSelection <- selectMetabolites_DF%>%
    mutate(`top/Bottom_percentage` = case_when(`top/Bottom`==TRUE ~ 'TRUE',
                                 TRUE ~ 'FALSE'))%>%
    mutate(Significant = case_when(get(metadata_info[["pvalColumn"]]) <= cutoff_stat ~ 'TRUE',
                                  TRUE ~ 'FALSE'))%>%
    mutate(Cluster_ChangedMetabolites = case_when(Significant==TRUE & `top/Bottom_percentage`==TRUE ~ 'TRUE',
                                   TRUE ~ 'FALSE'))
  InputSelection$`top/Bottom` <- NULL #remove column as its not needed for output

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
  clusterGo <- enricher_internal(
      gene=selectMetabolites,
      pvalueCutoff = 1L,
      pAdjustmethod = "BH",
      universe = allMetabolites,
      minGSSize=min_gssize,
      maxGSSize=max_gssize,
      qvalueCutoff = 1L,
      TERM2GENE=Term2gene ,
      TERM2NAME = term2name
  )
  clusterGosummary <- data.frame(clusterGo)

  #Make DF:
  if(!(dim(clusterGosummary)[1] == 0)){
      #Add pathway information % of genes in pathway detected)
      clusterGosummary <- merge(x= clusterGosummary%>% select(-Description), y=Pathway%>% select(term, Metabolites_in_Pathway),by.x="ID",by.y="term", all=TRUE)
      clusterGosummary$Count[is.na(clusterGosummary$Count)] <- 0
      clusterGosummary$percentage_of_Pathway_detected <-round(((clusterGosummary$Count/clusterGosummary$Metabolites_in_Pathway)*100),digits=2)
      clusterGosummary <- clusterGosummary[!duplicated(clusterGosummary$ID),]
      clusterGosummary <- clusterGosummary[order(clusterGosummary$p.adjust),]
      clusterGosummary <- clusterGosummary%>%
        dplyr::rename("Metabolites_in_pathway"="geneID")
  }else{
    stop("None of the Input_data Metabolites were present in any terms of the input_pathway. Hence the ClusterGosummary ouput will be empty. Please check that the metabolite IDs match the pathway IDs.")
    }

  #Return and save list of DFs
  ORA_output_list <- list("InputSelection" = InputSelection , "ClusterGosummary" = clusterGosummary)

  #save:
  suppressMessages(suppressWarnings(
    SaveRes(inputlist_df=ORA_output_list,
                         inputlist_plot= NULL,
                         save_table=save_table,
                         save_plot=NULL,
                         path= Folder,
                         FileName= paste(pathway_name),
                         core=FALSE,
                         print_plot=FALSE)))

  #Return
  ORA_output_list <- c( ORA_output_list, list("ClusterGo"=clusterGo))

  invisible(return(ORA_output_list))
  }


#################################
### ### ### Helper Fishers exact test ### ### ###
#################################
#'
#'

# perform test on n groups of features. DF with column feature and column group
# output list of DFs named after groups
