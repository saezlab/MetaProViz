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

ClusterORA <- function(se, #InputData,
    SettingsInfo = c(ClusterColumn = "RG2_Significant", 
        BackgroundColumn = "BG_Method", PathwayTerm = "term", 
        PathwayFeature = "Metabolite"),
    RemoveBackground = TRUE,
    PathwayFile,
    PathwayName = "",
    minGSSize = 10,
    maxGSSize = 1000 ,
    SaveAs_Table = "csv",
    FolderPath = NULL) {
  
    ## ------------ Create log file ----------- ##
    MetaProViz_Init()

    ## ------------ Check Input files ----------- ##
    Pathways <- CheckInput_ORA(se = se, ##InputData = InputData,
        SettingsInfo = SettingsInfo, RemoveBackground = RemoveBackground,
        PathwayFile = PathwayFile, PathwayName = PathwayName,
        minGSSize = minGSSize, maxGSSize = maxGSSize, 
        SaveAs_Table = SaveAs_Table, pCutoff = NULL, PercentageCutoff = NULL)

    ## ------------ Create Results output folder ----------- ##
    if (!is.null(SaveAs_Table)) {
        Folder <- SavePath(FolderName = "ClusterORA", FolderPath = FolderPath)
    }

    ############################################################################
    ## ------------ Prepare Data ----------- ##
    InputData <- assay(se)
    ## open the data
    if (RemoveBackground) {
        df <- subset(InputData, 
            !InputData[[SettingsInfo[["BackgroundColumn"]]]] == "FALSE") %>%
            rownames_to_column("Metabolite")
    } else{
        df <- InputData %>%
            rownames_to_column("Metabolite")
    }

    ## select universe
    allMetabolites <- as.character(df$Metabolite)

    ## select clusters
    grps_labels <- unlist(unique(df[[SettingsInfo[["ClusterColumn"]]]]))

    ## load Pathways
    Pathway <- Pathways
    ## term and MetaboliteID (MetaboliteID= gene as syntax required for enricher)
    Term2gene <- Pathway[, c("term", "gene")]
    ## term and description
    term2name <- Pathway[, c("term", "Description")]

    ## add the number of genes present in each pathway
    Pathway$Count <- 1
    Pathway_Mean <- aggregate(Pathway$Count, 
        by = list(term = Pathway$term), FUN = sum)
    names(Pathway_Mean)[names(Pathway_Mean) == "x"] <- "Metabolites_in_Pathway"
    Pathway <- merge(x = Pathway, y = Pathway_Mean, by = "term", all.x = TRUE)
    Pathway$Count <- NULL

    ## ------------ Run ----------- ##
    ## make an empty list to store the created DFs
    df_list = list()
    clusterGo_list = list()
    
    ## run ORA
    for (g in grps_labels) {
        grpMetabolites <- subset(df, df[[SettingsInfo[["ClusterColumn"]]]] == g)
        message("Number of metabolites in cluster `", g, "`: ",
            nrow(grpMetabolites), sep = "")

        clusterGo <- clusterProfiler::enricher(
            gene = as.character(grpMetabolites$Metabolite),
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            universe = allMetabolites,
            minGSSize = minGSSize,
            maxGSSize = maxGSSize,
            qvalueCutoff = 1,
            TERM2GENE = Term2gene ,
            TERM2NAME = term2name)
        clusterGoSummary <- data.frame(clusterGo)
        clusterGo_list[[g]]<- clusterGo
        if (!(nrow(clusterGoSummary) == 0)) {
            ## add pathway information (% of genes in pathway detected)
            clusterGoSummary <- merge(
                x = select(clusterGoSummary, -Description), 
                y = select(Pathway, term, Metabolites_in_Pathway), 
                by.x = "ID", by.y = "term", all = TRUE)
            clusterGoSummary$Count[is.na(clusterGoSummary$Count)] <- 0
            clusterGoSummary$Percentage_of_Pathway_detected <-round(
                clusterGoSummary$Count / clusterGoSummary$Metabolites_in_Pathway * 100,
                digits = 2)
            clusterGoSummary <- clusterGoSummary[!duplicated(clusterGoSummary$ID), ]
            clusterGoSummary <- clusterGoSummary[order(clusterGoSummary$p.adjust), ]
            clusterGoSummary <- clusterGoSummary %>%
                dplyr::rename("Metabolites_in_pathway"="geneID")

            g_save <- gsub("/", "-", g)
            df_list[[g_save]] <- clusterGoSummary
        } else {
            message("None of the Input_data Metabolites of the cluster ", g,
                " were present in any terms of the PathwayFile. Hence the ClusterGoSummary ouput will be empty for this cluster. Please check that the metabolite IDs match the pathway IDs.")
        }
    }
  
    ## save files
    suppressMessages(suppressWarnings(
        SaveRes(InputList_DF = df_list,
            InputList_Plot = NULL,
            SaveAs_Table = SaveAs_Table,
            SaveAs_Plot = NULL,
            FolderPath = Folder,
            FileName = paste("ClusterGoSummary", PathwayName, sep = "_"),
            CoRe = FALSE,
            PrintPlot = FALSE)))

    ## return <- clusterGoSummary
    ORA_Output <- list(
        "data" = list(
            "DF" = df_list, 
            "ClusterGo" = clusterGo_list)
    )

    ## return
    invisible(ORA_Output)
}


###################################
### ### ### StandardORA ### ### ###
###################################

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
#'
#' @export
#'
StandardORA <- function(se,
    SettingsInfo = c(pvalColumn = "p.adj", PercentageColumn = "t.val", 
        PathwayTerm = "term", PathwayFeature = "Metabolite"),
    pCutoff = 0.05,
    PercentageCutoff = 10,
    PathwayFile,
    PathwayName = "",
    minGSSize = 10,
    maxGSSize = 1000 ,
    SaveAs_Table = "csv",
    FolderPath = NULL) {
  
    ## ------------ Create log file ----------- ##
    MetaProViz_Init()

    ## ------------ Check Input files ----------- ##
    Pathways <- CheckInput_ORA(se = se,
        SettingsInfo = SettingsInfo,
        RemoveBackground = FALSE,
        PathwayFile = PathwayFile,
        PathwayName = PathwayName,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        SaveAs_Table = SaveAs_Table,
        pCutoff = pCutoff,
        PercentageCutoff = PercentageCutoff)

    ## ------------ Create Results output folder ----------- ##
    if(!is.null(SaveAs_Table)) {
        Folder <- SavePath(FolderName= "ORA", FolderPath = FolderPath)
    }

    ############################################################################
    ## ------------ Load the data and check ----------- ##
    #InputData <- se |> 
    #    assay() |> 
    #    rownames_to_column("Metabolite")

    ## select universe
    allMetabolites <- assay(se)$KEGGCompound

    ## select top changed metabolites (Up and down together)
    ## check if the metabolites are significantly changed.
    value <- PercentageCutoff / 100

    ## rank by t.val
    allMetabolites_DF <- assay(se)[
        order(assay(se)[, SettingsInfo[["PercentageColumn"]]]), ]
    selectMetabolites_DF <- allMetabolites_DF[
        c(seq_len(ceiling(value * nrow(se))),
            (nrow(se) - (ceiling(value * nrow(se)))):(nrow(se))), ]
    selectMetabolites_DF$`Top/Bottom` <- TRUE
    selectMetabolites_DF <- merge(allMetabolites_DF, 
        selectMetabolites_DF[, c("Top/Bottom"), drop = FALSE], 
        by = "row.names", all.x = TRUE)

    InputSelection <- selectMetabolites_DF %>%
        mutate(`Top/Bottom_Percentage` = case_when(
            `Top/Bottom` ~ TRUE,
            TRUE ~ FALSE)) %>%
        mutate(Significant = case_when(
            get(SettingsInfo[["pvalColumn"]]) <= pCutoff ~ TRUE,
            TRUE ~ FALSE)) %>%
        mutate(Cluster_ChangedMetabolites = case_when(
            Significant & `Top/Bottom_Percentage` ~ TRUE,
            TRUE ~ FALSE))
    ## remove column as its not needed for output
    InputSelection$`Top/Bottom` <- NULL 

    selectMetabolites <- InputSelection %>%
        subset(Cluster_ChangedMetabolites)
    selectMetabolites <- as.character(selectMetabolites$KEGGCompound)

    ## load Pathways
    Pathway <- Pathways
    ## term and MetaboliteID (MetaboliteID= gene as syntax required for enricher)
    Term2gene <- Pathway[, c("term", "gene")]
    ## term and description
    term2name <- Pathway[, c("term", "Description")]

    ## add the number of genes present in each pathway
    Pathway$Count <- 1
    Pathway_Mean <- aggregate(Pathway$Count, 
        by = list(term = Pathway$term), FUN = sum)
    names(Pathway_Mean)[names(Pathway_Mean) == "x"] <- "Metabolites_in_Pathway"
    Pathway <- merge(x = Pathway, y = Pathway_Mean, by = "term", all.x = TRUE)
    Pathway$Count <- NULL

    ## ------------ Run ----------- ##
    ## run ORA
    clusterGo <- clusterProfiler::enricher(gene = selectMetabolites,
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        universe = allMetabolites,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        qvalueCutoff = 1,
        TERM2GENE = Term2gene ,
        TERM2NAME = term2name)
    clusterGoSummary <- data.frame(clusterGo)

    ## make DF:
    if(!nrow(clusterGoSummary) == 0){
        ## add pathway information % of genes in pathway detected)
        clusterGoSummary <- merge(x = select(clusterGoSummary, -Description), 
            y = select(Pathway, term, Metabolites_in_Pathway), 
            by.x = "ID", by.y = "term", all = TRUE)
        clusterGoSummary$Count[is.na(clusterGoSummary$Count)] <- 0
        clusterGoSummary$Percentage_of_Pathway_detected <- round(
            clusterGoSummary$Count / clusterGoSummary$Metabolites_in_Pathway * 100,
            digits = 2)
        clusterGoSummary <- clusterGoSummary[!duplicated(clusterGoSummary$ID), ]
        clusterGoSummary <- clusterGoSummary[order(clusterGoSummary$p.adjust), ]
        clusterGoSummary <- clusterGoSummary %>%
            dplyr::rename("Metabolites_in_pathway" = "geneID")
    } else {
        stop("None of the Input_data Metabolites were present in any terms of the PathwayFile. Hence the ClusterGoSummary ouput will be empty. Please check that the metabolite IDs match the pathway IDs.")
    }
    
    ## create SummarizedExperiment output
    cols_a <- c("RichFactor", "FoldEnrichment", "zScore", "pvalue", "p.adjust", "qvalue", "Count", "Metabolites_in_Pathway", "Percentage_of_Pathway_detected")
    a <- as.matrix(clusterGoSummary[, cols_a])
    cD <- DataFrame(feature = cols_a)
    rownames(cD) <- cols_a
    rD <- DataFrame(clusterGoSummary[, !colnames(clusterGoSummary) %in% cols_a])
    rownames(rD) <- rownames(a)
    se_go <- SummarizedExperiment(assays = a, rowData = rD, colData = cD)
    
    ## return and save list of DFs
    ORA_output_list <- list(
        "data" = list(
            "se" = se_go, 
            "InputSelection" = InputSelection , 
            "ClusterGoSummary" = clusterGoSummary)
    )

    ## save
    suppressMessages(suppressWarnings(
        SaveRes(data = ORA_output_list,
            plot = NULL, SaveAs_Table = SaveAs_Table,
            SaveAs_Plot = NULL, FolderPath = Folder,
            FileName = paste(PathwayName), CoRe = FALSE,
            PrintPlot = FALSE)))

    ## return
    ORA_output_list[["data"]][["ClusterGo"]] <- clusterGo
    invisible(ORA_output_list)
}


#################################
### ### ### Helper Fishers exact test ### ### ###
#################################
#'
#'

# perform test on n groups of features. DF with column feature and column group
# output list of DFs named after groups
