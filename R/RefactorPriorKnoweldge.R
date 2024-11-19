## ---------------------------
##
## Script name: GetPriorKnowledge
##
## Purpose of script: Create gene-metabolite sets for pathway enrichment analysis.
##
## Author: Christina Schmidt and Macabe Daley
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
#' @param SaveAs_Table \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt". \strong{Default = "csv"}
#' @param FolderPath {Optional:} String which is added to the resulting folder name \strong{Default = NULL}
#'
#' @return List with three DFs: 1) Original data and the new column of translated ids. 2) Mapping summary from Original ID to Translated. 3) Mapping summary from Translated to Original.
#'
#' @examples
#' KEGG_Pathways <- MetaProViz::LoadKEGG()
#' Res <- MetaProViz::TranslateID(InputData= KEGG_Pathways, SettingsInfo = c(InputID="MetaboliteID", GroupingVariable="term"), From = c("kegg"), To = c("pubchem","chebi","hmdb"), SaveAs_Table= "csv", FolderPath=NULL)

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

  # Refactoring: make arguments similar to OmnipathR::translate_ids
  # i.e. instead of From / To dynamic dots
  # hence we dont need to rename the column that includes the "From"
  MetaProViz_Init()

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

  ## ------------------  Create output folders and path ------------------- ##
  # Export to this location was not part of the original function, we should
  # set it up later - Denes
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
      ambiguity_groups =  SettingsInfo[['GroupingVariable']]
    )

  ## ------------------ Add information to the results and Create Summary------------------- ##
  ResList <- list()
  for(item in  To){
    #Extract and prepare table for each metabolite ID:
    ExpandID <-  TranslatedDF %>%
      dplyr::select(any_of(names(InputData)), dplyr::contains(item))  %>%
      tidyr::unnest(cols = all_of(dplyr::contains(item)))

    # Add information about instances across or within pathways!
    if(SettingsInfo[["GroupingVariable"]] %in% colnames(ExpandID)){
      ExpandID <- ExpandID %>% #many-to-many = within or across pathways? --> add column with this information
        group_by(MetaboliteID, term) %>%
        mutate(GroupingVariable = case_when(
          n_distinct(hmdb) > 1 & MetaboliteID_hmdb_to_ambiguity > 1 & MetaboliteID_hmdb_ambiguity== "one-to-many" & n_distinct(term) >=2 & duplicated(term)==TRUE ~ "one-to-many_Within-and-AcrossGroups",  # Multiple KEGG IDs, multiple terms --> should not happen!
          n_distinct(hmdb) > 1 & MetaboliteID_hmdb_to_ambiguity > 1 & MetaboliteID_hmdb_ambiguity== "one-to-many" & n_distinct(term) >=2 & duplicated(term)==FALSE ~ "one-to-many_AcrossGroups",  # Multiple KEGG IDs, multiple terms --> should not happen!
          n_distinct(hmdb) > 1 & MetaboliteID_hmdb_to_ambiguity > 1 & MetaboliteID_hmdb_ambiguity == "one-to-many" & n_distinct(term) <= 1 ~ "one-to-many_WithinGroups",  # Multiple KEGG IDs, same term
          TRUE ~ NA_character_  #
        )) %>%
        ungroup()%>%
        group_by(hmdb) %>%
        mutate(GroupingVariable = case_when(
          n_distinct(MetaboliteID) > 1 & MetaboliteID_hmdb_to_ambiguity > 1 & MetaboliteID_hmdb_ambiguity== "many-to-many" & n_distinct(term) == 1 ~ "many-to-many_WithinGroups",  # Multiple KEGG IDs, same term
          n_distinct(MetaboliteID) > 1 & MetaboliteID_hmdb_to_ambiguity > 1 & MetaboliteID_hmdb_ambiguity== "many-to-many" & n_distinct(term) > 1 ~ "many-to-many_AcrossGroups",  # Multiple KEGG IDs, multiple terms --> should not happen!
          TRUE ~  paste(GroupingVariable) #
        )) %>%
        ungroup()
      }

    #Create a summary file about the instances of one-to-many etc. also include a descriptive column that verbalizes issues
    # --> e.g. pathway inflation/deflation


    #return results
    ResList[[item]] <- ExpandID
  }

  #many-to-many = within or across pathways? --> add column with this information


  ## ------------------ Save the results ------------------- ##
  res <- list(
    InputDF = InputData,
    TranslatedDF = TranslatedDF)
}


##########################################################################################
### ### ### NEW ### ### ###
##########################################################################################

#

MappingAmbiguity <- function(
    InputData,
    SettingsInfo = c(InputID="MetaboliteID", GroupingVariable="term"),
    From = "kegg",
    To = c("pubchem","chebi","hmdb"),
    SaveAs_Table= "csv",
    FolderPath=NULL
){
  #was inspectID
  #Summary and translation one-to-many, many-to-one

  #Step 1: +/-term --> One Group needed
  #Step 2: Omnipath function to get numeric column summary of 1-to-9 map
  #Step 3: Case_when --> column ( do things before and not within case_when)


  #Step: create summary for the specific problem of metabolism


}


##########################################################################################
### ### ### Helper ### ### ### --> Filtering should be done here!
##########################################################################################

#' Clean translated ID in prior knowledge based on measured features
#'
#' @param InputData Dataframe with at least one column with the target (e.g. metabolite), you can add other columns such as source (e.g. term)
#' @param SettingsInfo
#'
#' @export
#'
DetectedID <- function(InputData,
                       SettingsInfo = c(InputID="MetaboliteID", GroupingVariable="term")
){

  ## ------------------ Check Input ------------------- ##
  # HelperFunction `CheckInput`

  # Function that can use Prior knowledge that has multiple IDs "X1, X2, X3" and use the measured features to remove IDs based on measured data.
  # This is to enxure that not two detected metabolites map to the same entry and if the original PK was translated to  ensure the one-to-many, many-to-one issues are taken care of (within and across DBs)




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

CheckMatchID <- function(InputData

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

#' Deal with detected input features.
#'
#' @param InputData Dataframe with at least one column with the target (e.g. metabolite), you can add other columns such as source (e.g. term)
#' @param SettingsInfo
#'
#' @export
#'
PossibleID <- function(InputData,
                       SettingsInfo = c(InputID="MetaboliteID", GroupingVariable="term")
){

  ## ------------------ Check Input ------------------- ##
  # HelperFunction `CheckInput`

  # A user has one HMDB IDs for their measured metabolites (one ID per measured peak) --> this is often the case as the user either gets a trivial name and they have searched for the ID themselves or because the facility only provides one ID at random
  # We have mapped the HMDB IDs with the pathways and 20 do not map
  # We want to check if it is because the pathways don't include them, or because the user just gave the wrong ID by chance (i.e. They picked D-Alanine, but the prior knowledge includes L-Alanine)

  # Do this by using structural information via  accessing the structural DB in OmniPath!
  # Output is DF with the original ID column and a new column with additional possible IDs based on structure

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



ClusterPK <- function(InputData,
                      SettingsInfo= c(InputID="MetaboliteID", GroupingVariable="term")

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

  # 2. Compute pairwise overlaps
  term_overlap <- combn(term_metabolites[[SettingsInfo[["GroupingVariable"]]]], 2, function(terms) {
    term1_ids <- term_metabolites$MetaboliteIDs[term_metabolites[[SettingsInfo[["GroupingVariable"]]]] == terms[1]][[1]]
    term2_ids <- term_metabolites$MetaboliteIDs[term_metabolites[[SettingsInfo[["GroupingVariable"]]]] == terms[2]][[1]]

    overlap <- length(intersect(term1_ids, term2_ids)) / length(union(term1_ids, term2_ids))
    data.frame(Term1 = terms[1], Term2 = terms[2], Overlap = overlap)
    }, simplify = FALSE) %>%
    dplyr::bind_rows()

  # 3. Cluster terms based on overlap threshold
  threshold <- 0.7 # Define similarity threshold
  term_clusters <- term_overlap %>%
    dplyr::filter(Overlap >= threshold) %>%
    dplyr::select(Term1, Term2)

  # 4. Merge cluster group information back to the original data
  df <- InputData %>%
    dplyr::left_join(term_metabolites %>% select(!!sym(SettingsInfo[["GroupingVariable"]]), Group), by = SettingsInfo[["GroupingVariable"]])




  #Maybe use igraph?
  #g <- igraph::graph_from_data_frame(term_clusters, directed = FALSE)
  #clusters <- igraph::clusters(g)$membership
  #term_metabolites$Group <- clusters[match(term_metabolites[[SettingsInfo[["GroupingVariable"]]]], names(clusters))]

}












