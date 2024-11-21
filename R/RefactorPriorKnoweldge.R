## ---------------------------
##
## Script name: GetPriorKnowledge
##
## Purpose of script: Create gene-metabolite sets for pathway enrichment analysis.
##
## Author: Christina Schmidt, Denes Turei and Macabe Daley
##
## Date Created: 2024-01-21
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

##########################################################################################
### ### ### Translate IDs to/from KEGG, PubChem, Chebi, HMDB ### ### ###
##########################################################################################

#' Translate IDs to/from KEGG, PubChem, Chebi, HMDB
#'
#' @param InputData Dataframe with at least one column with the target (e.g. metabolite), you can add other columns such as source (e.g. term)
#' @param SettingsInfo \emph{Optional: } Column name of Target in Input_GeneSet. \strong{Default = list(InputID="MetaboliteID" , GroupingVariable="term")}
#' @param From ID type your
#' @param To One or multiple ID types
#' @param SaveAs_Table \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt". \strong{Default = "csv"}
#' @param FolderPath {Optional:} String which is added to the resulting folder name \strong{Default = NULL}
#'
#' @return List with three DFs: 1) Original data and the new column of translated ids. 2) Mapping summary from Original ID to Translated. 3) Mapping summary from Translated to Original.
#'
#' @examples
#' KEGG_Pathways <- MetaProViz::LoadKEGG()
#' Res <- MetaProViz::TranslateID(InputData= KEGG_Pathways, SettingsInfo = c(InputID="MetaboliteID", GroupingVariable="term"), From = c("kegg"), To = c("pubchem","chebi","hmdb"), SaveAs_Table= "csv", FolderPath=NULL)
#'
#'
#' @keywords Translate IDs
#'
#' @importFrom dplyr mutate
#' @importFrom rlang !!! !! := sym syms
#' @importFrom tidyselect everything starts_with
#' @importFrom dplyr across summarize first ungroup group_by select
#' @importFrom OmnipathR id_types translate_ids
#' @importFrom logger log_warn
#' @importFrom stringr str_to_lower
#'
#' @export
#'
TranslateID <- function(
    InputData,
    SettingsInfo = c(InputID="MetaboliteID", GroupingVariable="term"),
    From = "kegg",
    To = c("pubchem","chebi","hmdb"),
    SaveAs_Table= "csv",
    FolderPath=NULL
  ){

   MetaProViz_Init()

  ## ------------------  Check Input ------------------- ##

  # Specific checks:
  unknown_types <-
    OmnipathR::id_types() %>%
    dplyr::select(tidyselect::starts_with('in_')) %>%
    unlist %>%
    unique %>%
    str_to_lower %>%
    setdiff(union(From, To), .)

  if (length(unknown_types) > 0L) {

    msg <- sprintf(
      'The following ID types are not recognized: %s',
      paste(unknown_types, collapse = ', ')
    )
    logger::log_warn(msg)
    warning(msg)

  }

  # Check that SettingsInfo[['InputID']] has no duplications within one group --> should not be the case --> remove duplications and inform the user/ ask if they forget to set groupings column
  doublons <- InputData %>%
    dplyr::group_by(!!sym(SettingsInfo[['InputID']]), !!sym(SettingsInfo[['GroupingVariable']]))%>%
    dplyr::filter(n() > 1) %>%
    dplyr::ungroup()

  if(nrow(doublons) > 0){
    message <- sprintf(
      'The following ID types are duplicated within one group: %s',
      paste(doublons, collapse = ', ')
    )
    logger::log_warn(message)
    warning(message)
  }


  ## ------------------  Create output folders and path ------------------- ##
  Folder <- MetaProViz:::SavePath(FolderName = "TranslateID", FolderPath = NULL)

  ## ------------------ Translate To-From for each pair ------------------- ##
  TranslatedDF <- OmnipathR::translate_ids(
      InputData,
      !!sym(SettingsInfo[['InputID']]) :=  !!sym(From),
      !!!syms(To),#list of symbols, hence three !!!
      ramp = TRUE,
      expand = FALSE,
      quantify_ambiguity = TRUE,
      qualify_ambiguity = TRUE,
      ambiguity_groups =  SettingsInfo[['GroupingVariable']]#Checks within the groups, without it checks across groups
    )

  #replace 0 with NA!

  ## ------------------ Create Tables for each translation ------------------- ##

  ResList <- list()
  for(item in To){
    DF <- TranslatedDF %>%
      dplyr::select(any_of(names(InputData)), item)%>%
      tidyr::unnest(cols = item)



  }







    # Task 1: Check that SettingsInfo[['InputID']] has the same items in to across the different entries (would be in different Groupings, otherwise there should not be any duplications)


    #Create Mapping column
    item_columns <- TranslatedDF %>%
      dplyr::select(item) %>%
      apply(1, function(row) paste(row, collapse = " + "))

    Mapping <- TranslatedDF %>%
      dplyr::mutate(Mapping = paste(!!sym(SettingsInfo[['InputID']]), "=", item_columns))%>%
      dplyr::select(any_of(names(InputData)), dplyr::contains(item), "Mapping")

    ExpandID <- TranslatedDF %>%
      dplyr::select(any_of(names(InputData)), dplyr::contains(item))  %>%
      tidyr::unnest(cols = all_of(dplyr::contains(item)))

    #return results
    ResList[[item]] <- ExpandID







  ## ------------------
  ambi <-
  To %>%
  rlang::set_names() %>%
  purrr::map(
    ~MappingAmbiguity(
      TranslatedDF %>% select(all_of(c(names(InputData), .x))),
      SettingsInfo[['InputID']],
      .x,
      SettingsInfo[['GroupingVariable']]
    )
  )






  # dplyr::select(any_of(names(InputData)), dplyr::contains(item))  %>%
  #   tidyr::unnest(cols = all_of(dplyr::contains(item)))
  #
  #
  #
  # ## ------------------ Save the results ------------------- ##
  # res <- list(
  #   InputDF = InputData,
  #   TranslatedDF = TranslatedDF)


  #save function to folder

  #return


}


##########################################################################################
### ### ### NEW ### ### ###
##########################################################################################

#this is now part of OmnipathR --> denes said it will fit better there!
# Problem: We need to be able to use this function here too. So we will need to extract it from Omnipath to also be present here! (either wrapper or copy!)
# needs to be @export

#' Translate IDs to/from KEGG, PubChem, Chebi, HMDB
#'
#' @param InputData Translated DF from MetaProViz::TranslateID reults or Dataframe with at least one column with the target ID (e.g. metabolite KEGG IDs) and another MetaboliteID type (e.g. KEGG and HMDB). Optional: add other columns such as source (e.g. term) or more metabolite IDs.
#' @param To Column name of original metabolite identifier in InputData. Here should only be one ID per row
#' @param From Column name of the secondary or translated metabolite identifier in InputData. Here can be multiple IDs per row.
#' @param GroupingVariable \emph{Optional: } If NULL no groups are used. If TRUE provide column name in InputData containing the GroupingVariable; features are grouped. \strong{Default = NULL}
#' @param SaveAs_Table \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt". \strong{Default = "csv"}
#' @param FolderPath {Optional:} String which is added to the resulting folder name \strong{Default = NULL}
#'
#' @return List with three DFs: 1) Original data and the new column of translated ids. 2) Mapping summary from Original ID to Translated. 3) Mapping summary from Translated to Original.
#'
#' @examples
#' KEGG_Pathways <- MetaProViz::LoadKEGG()
#' Res <- MetaProViz::TranslateID(InputData= KEGG_Pathways, SettingsInfo = c(InputID="MetaboliteID", GroupingVariable="term"), From = c("kegg"), To = c("pubchem","chebi","hmdb"), SaveAs_Table= "csv", FolderPath=NULL)
#'
#'
#' @keywords Mapping ambiguity
#'
#' @importFrom dplyr mutate
#' @importFrom rlang !!! !! := sym syms
#' @importFrom OmnipathR ambiguity
#'
#' @export
#'

MappingAmbiguity <- function( InputData=TranslatedDF,
                              #SettingsInfo = c(InputID="MetaboliteID", GroupingVariable="term"),
                              From = "MetaboliteID",
                              To = "hmdb",
                              GroupingVariable = NULL,
                              SaveAs_Table= "csv",
                              FolderPath=NULL
) {

  ## ------------------  Check Input ------------------- ##


  ## ------------------  Create output folders and path ------------------- ##




  ## ------------------  General checks of wrong occurences ------------------- ##
  # Task 1: Check that SettingsInfo[['InputID']] has no duplications within one group --> should not be the case --> remove duplications and inform the user/ ask if they forget to set groupings column
  # Task 2: Check that SettingsInfo[['InputID']] has the same items in to across the different entries (would be in different Groupings, otherwise there should not be any duplications) --> List of Miss-Mappings across terms

  # FYI: The above can not happen if our translateID function was used, but may be the case when the user has done something manually before



  ## Create input


  ids <- c(From, To)
  ResList <- c(from_to_to = identity, to_to_from = rev) %>%
    map(
      function(direction) {
        cols <- ids %>% direction

        InputData %>%
          dplyr::select(all_of(c(cols, GroupingVariable))) %>%
          tidyr::unnest(cols) %>%# unlist the columns in case they are not expaned
          OmnipathR::ambiguity(
          from_col = !!sym(cols[1L]),
          to_col = !!sym(cols[2L]),
          groups = GroupingVariable,
          quantify = 'Group',
          qualify = 'Group',
          #global = TRUE,#across groups will be done additionally --> suffix _AcrossGroup
          #summary=TRUE, #summary of the mapping column
          expand = TRUE
        )
      }
    )


  #2. --> use different column names:  "MetaboliteID_hmdb_to_ambiguity" = "To_hmdb" , "MetaboliteID_hmdb_from_ambiguity" = "From_MetaboliteID"
  #3. --> different order of columns: InputData columns, "hmdb", "From_MetaboliteID", "To_MetaboliteID"



  #---- Summary

  for(df in names(ResList)){
    #Flank problematic cases by adding further columns to df


    #Create Summary
    Summary <- ResList[[df]] %>%
      count(paste0())%>%#Column we will need to distinguish
      pivot_wider(names_from = Relationship, values_from = n, values_fill = 0) #n=count --> column names are one-to-one, etc. and entries are the occurrences



  }

  ## ------------------  Perform ambiguity mapping ------------------- ##
  #1. From-to-To: OriginalID-to-TranslatedID
  #2. From-to-To: TranslatedID-to-OriginalID
  Comp <- list(
    list(From = From, To = To),
    list(From = To, To = From)
  )

  ResList <- list()
  for(comp in seq_along(Comp)){
  # Change OmnipathR::ambiguity
    #1. --> create a column "From-to-To" (= hmdb-to-kegg) only with one-to-one, one-to-none, one-to-many, create a column "From-to-To_Bidirectional" (= hmdb-to-kegg_Bidirectional) one-to-one, one-to-none, one-to-many AND many-to-many
    #2. --> use different column names:  "MetaboliteID_hmdb_to_ambiguity" = "To_hmdb" , "MetaboliteID_hmdb_from_ambiguity" = "From_MetaboliteID"
    #3. --> different order of columns: InputData columns, "hmdb", "From_MetaboliteID", "To_MetaboliteID"
    if(is.null(GroupingVariable)){
      ResList[[paste0(Comp[[comp]]$From, "-to-", Comp[[comp]]$To, "_OneGroup", sep="")]] <- InputData %>%
        OmnipathR::ambiguity(
        from_col = !!sym(Comp[[comp]]$From),
        to_col = !!sym(Comp[[comp]]$To),
        groups = NULL,
        quantify = 'OneGroup', # maybe name NoGrouping?
        qualify = 'OneGroup',
        expand=TRUE)
    }else{
      ResList[[paste0(Comp[[comp]]$From, "-to-", Comp[[comp]]$To, "_WithinGroup", sep="")]] <- InputData %>%
        OmnipathR::ambiguity(
        from_col = !!sym(Comp[[comp]]$From),
        to_col = !!sym(Comp[[comp]]$To),
        groups = GroupingVariable,
        quantify = 'WithinGroup',
        qualify = 'WithinGroup',
        expand=TRUE)

      ResList[[paste0(Comp[[comp]]$From, "-to-", Comp[[comp]]$To,"_AcrossGroup", sep="")]] <- InputData %>%
        OmnipathR::ambiguity(
        from_col = !!sym(Comp[[comp]]$From),
        to_col = !!sym(Comp[[comp]]$To),
        groups = NULL,
        quantify = 'AcrossGroup',
        qualify = 'AcrossGroup',
        expand=TRUE)
    }



    }


  ## ------------------  summaries the results and log messages ------------------- ##

  for(df in names(ResList)){
    #Flank problematic cases by adding further columns to df


    #Create Summary
    Summary <- ResList[[df]] %>%
      count(Mapping)%>%#Column we will need to distinguish
      pivot_wider(names_from = Relationship, values_from = n, values_fill = 0) #n=count --> column names are one-to-one, etc. and entries are the occurrences





  }






  # if(SettingsInfo[["GroupingVariable"]] %in% colnames(ExpandID)){
  #   ExpandID <- ExpandID %>% #many-to-many = within or across pathways? --> add column with this information
  #      group_by(MetaboliteID, term) %>%
  #     mutate(GroupingVariable = case_when(
  #       n_distinct(hmdb) > 1 & MetaboliteID_hmdb_to_ambiguity > 1 & MetaboliteID_hmdb_ambiguity== "one-to-many" & n_distinct(term) >=2 & duplicated(term)==TRUE ~ "one-to-many_Within-and-AcrossGroups",  # Multiple KEGG IDs, multiple terms --> should not happen!
  #       n_distinct(hmdb) > 1 & MetaboliteID_hmdb_to_ambiguity > 1 & MetaboliteID_hmdb_ambiguity== "one-to-many" & n_distinct(term) >=2 & duplicated(term)==FALSE ~ "one-to-many_AcrossGroups",  # Multiple KEGG IDs, multiple terms --> should not happen!
  #       n_distinct(hmdb) > 1 & MetaboliteID_hmdb_to_ambiguity > 1 & MetaboliteID_hmdb_ambiguity == "one-to-many" & n_distinct(term) <= 1 ~ "one-to-many_WithinGroups",  # Multiple KEGG IDs, same term
  #       TRUE ~ NA_character_  #
  #     )) %>%
  #     ungroup()%>%
  #     group_by(hmdb) %>%
  #     mutate(GroupingVariable = case_when(
  #       n_distinct(MetaboliteID) > 1 & MetaboliteID_hmdb_to_ambiguity > 1 & MetaboliteID_hmdb_ambiguity== "many-to-many" & n_distinct(term) == 1 ~ "many-to-many_WithinGroups",  # Multiple KEGG IDs, same term
  #       n_distinct(MetaboliteID) > 1 & MetaboliteID_hmdb_to_ambiguity > 1 & MetaboliteID_hmdb_ambiguity== "many-to-many" & n_distinct(term) > 1 ~ "many-to-many_AcrossGroups",  # Multiple KEGG IDs, multiple terms --> should not happen!
  #       TRUE ~  paste(GroupingVariable) #
  #     )) %>%
  #     ungroup()
  # }


  #Create a summary file about the instances of one-to-many etc. also include a descriptive column that verbalizes issues
  # --> e.g. pathway inflation/deflation




  ## ------------------  Save the results ------------------- ##




  #was inspectID
  #Summary and translation one-to-many, many-to-one

  #Step 1: +/-term --> One Group needed
  #Step 2: Omnipath function to get numeric column summary of 1-to-9 map
  #Step 3: Case_when --> column ( do things before and not within case_when)


  #Step: create summary for the specific problem of metabolism


}


##########################################################################################
### ### ### Reduce mapping ambiguities by using detected IDs  ### ### ###
##########################################################################################

#' Reduce mapping ambiguities by using detected IDs
#'
#' @param InputData Dataframe with at least one column with the target (e.g. metabolite), you can add other columns such as source (e.g. term)
#' @param SettingsInfo
#'
#' @export
#'
CleanMapping <- function(InputData,
                       SettingsInfo = c(InputID="MetaboliteID", GroupingVariable="term")
){

  ## ------------------ Check Input ------------------- ##
  # HelperFunction `CheckInput`

  # Function that can use Prior knowledge that has multiple IDs "X1, X2, X3" and use the measured features to remove IDs based on measured data.
  # This is to ensure that not two detected metabolites map to the same entry and if the original PK was translated to  ensure the one-to-many, many-to-one issues are taken care of (within and across DBs)

  # --> Cleans translated ID in prior knowledge based on measured features


}


##########################################################################################
### ### ### Check Measured ID's in prior knowledge ### ### ###
##########################################################################################

#' Check and summarize PriorKnowledge-to-MeasuredFeatures relationship
#'
#' @param InputData Dataframe with at least one column with the target (e.g. metabolite), you can add other columns such as source (e.g. term)
#' @param SettingsInfo
#'
#' @return
#'
#' @examples
#'
#' @keywords
#'
#' @importFrom dplyr mutate
#' @importFrom rlang !!! !! := sym syms
#'
#' @export
#'

CheckMatchID <- function(InputData #Maybe name DetectedID

){
  ## ------------ Create log file ----------- ##
  MetaProViz_Init()

  ## ------------ Check Input files ----------- ##

  ## ------------ Create Results output folder ----------- ##
  if(is.null(SaveAs_Table)==FALSE){
    Folder <- SavePath(FolderName= "PriorKnowledgeChecks",
                       FolderPath=FolderPath)
  }
  ################################################################################################################################################################################################
  ## ------------ Prepare the Input -------- ##





}


##########################################################################################
### ### ### Helper (?) ### ### ###
##########################################################################################

#' Deal with detected input features. If we only have one ID, could it also be another one?
#'
#' @param InputData Dataframe with at least one column with the target (e.g. metabolite), you can add other columns such as source (e.g. term)
#' @param SettingsInfo
#'
#' @export
#'
AssignID <- function(InputData,
                       SettingsInfo = c(InputID="MetaboliteID", GroupingVariable="term")
){

  ## ------------------ Check Input ------------------- ##
  # HelperFunction `CheckInput`

  # A user has one HMDB IDs for their measured metabolites (one ID per measured peak) --> this is often the case as the user either gets a trivial name and they have searched for the ID themselves or because the facility only provides one ID at random
  # We have mapped the HMDB IDs with the pathways and 20 do not map
  # We want to check if it is because the pathways don't include them, or because the user just gave the wrong ID by chance (i.e. They picked D-Alanine, but the prior knowledge includes L-Alanine)

  # Do this by using structural information via  accessing the structural DB in OmniPath!
  # Output is DF with the original ID column and a new column with additional possible IDs based on structure

  #Is it possile to do this at the moment without structures, but by using other pior knowledge?

}


##########################################################################################
### ### ### Cluster Prior Knowledge ### ### ###
##########################################################################################

#' Deal with pathway overlap in prior knowledge
#'
#' @param InputData Dataframe with at least one column with the target (e.g. metabolite) and a column source (e.g. term).
#' @param SettingsInfo = c(InputID="MetaboliteID", GroupingVariable="term"),
#' @examples
#' KEGG_Pathways <- MetaProViz::LoadKEGG()
#' InputData = KEGG_Pathways
#'
#'
#'



ClusterPK <- function(InputData, # This can be either the original PK (e.g. KEGG pathways), but it can also be the output of enrichment results (--> meaning here we would cluster based on detection!)
                      SettingsInfo= c(InputID="MetaboliteID", GroupingVariable="term"),
                      Clust = "Graph", # Options: "Graph", "Hierarchical",
                      matrix ="percentage", # Choose "pearson", "spearman", "kendall", or "percentage"
                      min= 2 # minimum pathways per cluster

){

  # Cluster PK before running enrichment analysis --> add another column that groups the data based on the pathway overlap:
  # provide different options for clustering (e.g. % of overlap, semantics similarity) --> Ramp uses % of overlap, semnatics similarity: https://yulab-smu.top/biomedical-knowledge-mining-book/GOSemSim.html


  ## ------------------ Check Input ------------------- ##


  ## ------------------ Create output folders and path ------------------- ##




  ## ------------------ Cluster the data ------------------- ##
  # 1. Create a list of unique MetaboliteIDs for each term
  term_metabolites <- InputData %>%
    dplyr::group_by(!!sym(SettingsInfo[["GroupingVariable"]])) %>%
    dplyr::summarize(MetaboliteIDs = list(unique(!!sym(SettingsInfo[["InputID"]])))) %>%
    dplyr::ungroup()

  #2. Create the overlap matrix based on different methods:
  if (matrix == "percentage") {# Compute pairwise overlaps
    term_overlap <- combn(term_metabolites[[SettingsInfo[["GroupingVariable"]]]], 2, function(terms) {
      term1_ids <- term_metabolites$MetaboliteIDs[term_metabolites[[SettingsInfo[["GroupingVariable"]]]] == terms[1]][[1]]
      term2_ids <- term_metabolites$MetaboliteIDs[term_metabolites[[SettingsInfo[["GroupingVariable"]]]] == terms[2]][[1]]

      overlap <- length(intersect(term1_ids, term2_ids)) / length(union(term1_ids, term2_ids))
      data.frame(Term1 = terms[1], Term2 = terms[2], Overlap = overlap)
    }, simplify = FALSE) %>%
      dplyr::bind_rows()

    # Create overlap matrix: An overlap matrix is typically used to quantify the degree of overlap between two sets or groups.
    # overlap coefficient (or Jaccard Index) Overlap(A,B)= ∣A∪B∣ / ∣A∩B
    #The overlap matrix measures the similarity between sets or groups based on common elements.
    terms <- unique(c(term_overlap$Term1, term_overlap$Term2))
    overlap_matrix <- matrix(1, nrow = length(terms), ncol = length(terms), dimnames = list(terms, terms))
    for (i in seq_len(nrow(term_overlap))) {
      t1 <- term_overlap$Term1[i]
      t2 <- term_overlap$Term2[i]
      overlap_matrix[t1, t2] <- 1 - term_overlap$Overlap[i]
      overlap_matrix[t2, t1] <- 1 - term_overlap$Overlap[i]
    }
  } else {
    # Create a binary matrix for correlation methods
    terms <- term_metabolites[[SettingsInfo[["GroupingVariable"]]]]
    metabolites <- unique(unlist(term_metabolites$MetaboliteIDs)) #[[SettingsInfo[["InputID"]]]]

    binary_matrix <- matrix(0, nrow = length(terms), ncol = length(metabolites), dimnames = list(terms, metabolites))
    for (i in seq_along(terms)) {
      metabolites_for_term <- term_metabolites$MetaboliteIDs[[i]] #[[SettingsInfo[["InputID"]]]]
      binary_matrix[i, colnames(binary_matrix) %in% metabolites_for_term] <- 1
    }

    # Compute correlation matrix: square matrix used to represent the pairwise correlation coefficients between variables or terms
    # correlation matrix 𝐶 C is an 𝑛 × 𝑛 n×n matrix where each element 𝐶 𝑖 𝑗 C ij  is the correlation coefficient between the variables 𝑋𝑖 Xi  and 𝑋 𝑗 X j
    #The correlation matrix measures the strength and direction of linear relationships between variables.
    correlation_matrix <- cor(t(binary_matrix), method = matrix)

    # Convert to distance matrix
    overlap_matrix <- 1 - correlation_matrix
  }
  # 3. Cluster terms based on overlap threshold
  threshold <- 0.7 # Define similarity threshold
  term_clusters <- term_overlap %>%
    dplyr::filter(Overlap >= threshold) %>%
    dplyr::select(Term1, Term2)

  # 4. Clustering
  if (Clust == "Graph") { #Use Graph-based clustering
  # Here we need the distance matrix:
  overlap_matrix <- 1 - correlation_matrix

  # An adjacency matrix represents a graph structure and encodes the relationships between nodes (vertices)
  # Add weight (can also represent unweighted graphs)
  adjacency_matrix <- exp(-overlap_matrix^2)  # Applying Gaussian kernel to convert distance into similarity

  # Create a graph from the adjacency matrix
  g <- igraph::graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected", weighted = TRUE)
  initial_clusters <- igraph::components(g)$membership
  term_metabolites$Cluster <- initial_clusters[match(term_metabolites[[SettingsInfo[["GroupingVariable"]]]], names(initial_clusters))]
  } else if (Clust == "Hierarchical") { # Hierarchical clustering
    hclust_result <- hclust(as.dist(distance_matrix), method = "average") # make methods into parameters!
    num_clusters <- 4
    term_clusters_hclust <- cutree(hclust_result, k = num_clusters)

    term_metabolites$Cluster <- paste0("Cluster", term_clusters_hclust[match(terms, names(term_clusters_hclust))])
    #term_metabolites$Cluster <- clusters[match(term_metabolites[[SettingsInfo[["GroupingVariable"]]]], names(clusters))]
  } else {
    stop("Invalid clustering method specified in Clust parameter.")
  }

  # 5. Merge cluster group information back to the original data
  df <- InputData %>%
    dplyr::left_join(term_metabolites %>% select(!!sym(SettingsInfo[["GroupingVariable"]]), Cluster), by = SettingsInfo[["GroupingVariable"]])%>%
    dplyr::mutate(Cluster = ifelse(
      is.na(Cluster),
        "None", # Assign "None" to NAs
        paste0("Cluster", Cluster) # Convert numeric IDs to descriptive labels
      )
    )

  # 6. Summarize the clustering results





  ## ------------------ Save and return ------------------- ##


}




##########################################################################################
### ### ### Helper function to add term information to Enrichment Results ### ### ###
##########################################################################################

#' Adds extra columns to enrichment output that inform about 1. The amount of genes associated with term in prior knowledge, 2. The amount of genes detected in input data associated with term in prior knowledge, and 3. The percentage of genes detected in input data associated with term in prior knowledge.
#'
#' @param mat Data matrix used as input for enrichment analysis
#' @param net Prior Knowledge used as input for enrichment analysis
#' @param res Results returned from the enrichment analysis
#' @param .source used as input for enrichment analysis
#' @param .target used as input for enrichment analysis
#' @param complete TRUE or FALSE, weather only .source with results should be returned or all .source in net.
#'
#' @export

# Better function Name and parameter names needed
# Use in ORA functions and showcase in vignette with decoupleR output

AddInfo <- function(mat,
                    net,
                    res,
                    .source,
                    .target,
                    complete=FALSE){

  ## ------------------ Check Input ------------------- ##


  ## ------------------ Create output folders and path ------------------- ##


  ## ------------------ Add information to enrichment results ------------------- ##

  # add number of Genes_targeted_by_TF_num
  net$Count <- 1
  net_Mean <- aggregate(net$Count, by=list(source=net[[.source]]), FUN=sum)%>%
    rename("targets_num" = 2)

  if(complete==TRUE){
    res_Add<- merge(x= res, y=net_Mean, by="source", all=TRUE)
  }else{
    res_Add<- merge(x= res, y=net_Mean, by="source", all.x =TRUE)
  }

  # add list of Genes_targeted_by_TF_chr
  net_List <- aggregate(net[[.target]]~net[[.source]], FUN=toString)%>%
    rename("source" = 1,
           "targets_chr"=2)
  res_Add<- merge(x= res_Add, y=net_List, by="source", all.x=TRUE)

  # add number of Genes_targeted_by_TF_detected_num
  mat <- as.data.frame(mat)%>% #Are these the normalised counts?
    tibble::rownames_to_column("Symbol")

  Detected <- merge(x= mat , y=net[,c(.source, .target)], by.x="Symbol", by.y=.target, all.x=TRUE)%>%
    filter(!is.na(across(all_of(.source))))
  Detected$Count <-1
  Detected_Mean <- aggregate(Detected$Count, by=list(source=Detected[[.source]]), FUN=sum)%>%
    rename("targets_detected_num" = 2)

  res_Add<- merge(x= res_Add, y=Detected_Mean, by="source", all.x=TRUE)%>%
    mutate(targets_detected_num = replace_na(targets_detected_num, 0))

  # add list of Genes_targeted_by_TF_detected_chr
  Detected_List <- aggregate(Detected$Symbol~Detected[[.source]], FUN=toString)%>%
    rename("source"=1,
           "targets_detected_chr" = 2)

  res_Add<- merge(x= res_Add, y=Detected_List, by="source", all.x=TRUE)

  #add percentage of Percentage_of_Genes_detected
  res_Add$targets_detected_percentage <-round(((res_Add$targets_detected_num/res_Add$targets_num)*100),digits=2)

  #sort by score
  res_Add<-res_Add%>%
    arrange(desc(as.numeric(as.character(score))))

  ## ------------------ Save and return ------------------- ##
  Output<-res_Add
}




