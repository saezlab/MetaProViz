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
#' res <- MetaProViz::TranslateID(InputData= KEGG_Pathways, SettingsInfo = c(InputID="MetaboliteID", GroupingVariable="term"), From = c("kegg"), To = c("pubchem","chebi","hmdb"), Method="GetAll", SaveAs_Table= "csv", FolderPath=NULL)
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
    From = c("kegg"),
    To = c("pubchem","chebi","hmdb"),
    SaveAs_Table= "csv",
    FolderPath=NULL
  ){

  MetaProViz_Init()

  From %<>% (stringr::str_to_lower)
  To %<>% (stringr::str_to_lower)

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
  list(
    InputDF = InputData,
    TranslatedDF = OmnipathR::translate_ids(
      InputData,
      !!rlang::sym(From),
      !!!rlang::syms(To),
      ramp = TRUE,
      expand = FALSE,
      inspect = TRUE,
      inspect_grp = SettingsInfo[['GroupingVariable']]
    )
  )

}


##########################################################################################
### ### ### Helper ### ### ###
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
