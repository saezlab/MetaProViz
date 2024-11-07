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
#' @importFrom tidyselect everything
#' @importFrom dplyr across summarize first ungroup group_by
#'
#' @export
#'
TranslateID <- function(InputData,
                        SettingsInfo = c(InputID="MetaboliteID", GroupingVariable="term"),
                        From = c("kegg"),
                        To = c("pubchem","chebi","hmdb"),
                        SaveAs_Table= "csv",
                        FolderPath=NULL
){

  ## ------------------ Check Input ------------------- ##
  # HelperFunction `CheckInput`


  # Specific checks:
  VALID_ID_TYPES <- c("kegg","pubchem","chebi","hmdb")
  From %<>% (stringr::str_to_lower)
  To %<>% (stringr::str_to_lower)

  if (!all(c(From, To) %in% VALID_ID_TYPES)) {
    m <- paste("We currently 'only' support ID translation between: ", paste(VALID_ID_TYPES, collapse = ", ") , ". You entered:",
               paste((c(From, To))[!(c(From, To)) %in% ID_options], collapse = ", "))
    logger::log_error(m)
    stop(m)
  }

  ## ------------------  Create output folders and path ------------------- ##
  Folder <- MetaProViz:::SavePath(FolderName= "TranslateID",
                                  FolderPath=NULL)

  ## ------------------ Messages (log later) ------------------- ##
  logger::log_trace(paste('Using Method "', Method, '"', sep=""))

  ## ------------------ Translate To-From for each pair ------------------- ##
  DF_List <- list()
  DF_List[["InputDF"]] <- InputData
  DF_List[["TranslatedDF"]] <- InputData

  for (MetaboliteID in To) {# For each ID type in the ToFormat list, we will translate the IDs and then collapse them into a single row for each group (e.g. by path term and metaboliteID)
    #----- 1. Translate
    df_translated <-
      DF_List[['TranslatedDF']] %>%
      OmnipathR::translate_ids(
        !!rlang::sym(SettingsInfo[['InputID']]) := !!rlang::sym(From),
        !!rlang::sym(MetaboliteID),
        ramp = TRUE
      )  # Use OmnipathR to translate the ids. Note that the returned object (df_translated) will most likely have multiple mappings- note, the results will be in long format likely with double ups.

    #----- 2. Collapse
    # Now collapse the desired translated column rows into a single row for each group: "16680, 57856, 181457"
    grouping_vars <- c(SettingsInfo[['GroupingVariable']], SettingsInfo[['InputID']])

    df_translated <- df_translated %>%
      dplyr::group_by(!!!rlang::syms(grouping_vars)) %>%
      dplyr::summarize(!!rlang::sym(MetaboliteID) := paste(!!rlang::sym(MetaboliteID), collapse = ', '), dplyr::across(tidyselect::everything(), dplyr::first)) %>%
      dplyr::ungroup()


    logger::log_success(paste('Converting from ', From, " to ", MetaboliteID))

    DF_List[["TranslatedDF"]] <- df_translated


    #----- 5. Mapping Inspection: One-to-Many and Many-to-One relationshipts
    # Now Inspect the IDs to identify potential problems with mapping - let's do this even if the user chose the GetFirst method, because it is very important to show the limitations...
    if("GroupingVariable" %in% names(SettingsInfo)){
      mapping_inspection = MetaProViz:::InspectID(InputData= df_translated,
                                                  SettingsInfo = c(OriginalID= From, TranslatedID=(MetaboliteID), GroupingVariable="term"))
    } else {
      mapping_inspection = MetaProViz:::InspectID(InputData= df_translated,
                                                SettingsInfo = c(OriginalID= From, TranslatedID=paste0(MetaboliteID, '_AllIDs')))
    }

    DF_List <- c(DF_List, setNames(list(mapping_inspection$Orig2Trans), paste0("Mapping_Orig2Trans_", From, '_to_', MetaboliteID)))
    DF_List <- c(DF_List, setNames(list(mapping_inspection$Trans2Orig), paste0("Mapping_Trans2Orig_", MetaboliteID, '_to_', From)))
  }




  ## ------------------ Create a final summary table of the mapping issues ------------------- ##
  # Create a function to get the counts
  get_counts <- function(table) {
    table %>%
      count(Relationship) %>%
      rename(n = n)  # Rename to n for consistency
  }

  # Extract counts for each table and add the table name
  summary_table <- map_dfr(
    seq_along(DF_List)[-1],
    function(i) {
      table <- DF_List[[i]]
      table_name <- names(DF_List)[i]
      counts <- get_counts(table)
      counts <- counts %>% dplyr::mutate(Table = table_name)  # Add table name as a new column
      return(counts)
    }
  )

  # Pivot to a wider format
  final_table_summary <- summary_table %>%
    pivot_wider(names_from = Relationship, values_from = n, values_fill = 0)

  # Reorder safely
  desired_order <- c("Table", "One-to-None","One-to-One","One-to-Many")
  existing_columns <- names(final_table_summary) # Get the existing columns in final_table
  safe_order <- intersect(desired_order, existing_columns) # Create a safe order by intersecting desired_order with existing_columns

  # Check if "Table" is in the safe order and ensure it's at the beginning
  if ("Table" %in% safe_order) {
    safe_order <- c("Table", setdiff(safe_order, "Table"))
  }

  final_table_summary <- final_table_summary %>%
    dplyr::select(all_of(safe_order))

  ## ------------------ Save and return DF ------------------- ##
  DF_List <- c(DF_List, setNames(list(final_table_summary), "TranslationSummary"))

  return(DF_List)
}


##########################################################################################
### ### ### Helper to check One-to-Many and Many-to-One Relationships ### ### ###
##########################################################################################
#'
#' Inspect ID
#'
#' @param InputData Dataframe input data frame should be what we have after the translation has occurred, and include a column with collapsed values in the case of multi-mapping (e.g. C1, C2, C3), but could be all 1-to-1 in rare cases.
#' @param SettingsInfo column name in InputData of the OriginalID and column name of TranslatedID.
#'
#' @description Inspect how well IDs map from the translated format (e.g. PubChem) to the original data format (e.g. KEGG), in terms of direct mapping, or one-to-many relationships.
#' @return Two data frames, the first of which is a summary of the mapping from Original to Translated, and the second of which is the reverse, from Translated to Original, with counts per unique ID and pathway.
#' @noRd
#'
#'


InspectID <- function(InputData,
                      SettingsInfo
                      ) {

  ## ------------------ Check Input  ------------------- ##

  #If the TranslatedID and OriginalID, both do not have multiple IDs (ID1, ID2, etc), return a message that we stop cause we do not have any occurences!




  ## ------------------ Prepare Input  ------------------- ##
  # Check if the GroupingVariable is in the InputData, if not add an GroupingVariable for all of them
  if("GroupingVariable" %in% names(SettingsInfo)==FALSE){
    SettingsInfo[['GroupingVariable']] <- "GroupingVariable"
    InputData$GroupingVariable <- "OneGroup"
    message("No GroupingVariable provided, hence assuming all metabolite IDs are in the same group.")
  }

  ## ------------------ Map IDs  ------------------- ##
  mapping_list <- list()

  #----- Iterate over rows of the dataframe
  for (i in 1:nrow(InputData)) {
    # Process the collapsed field (e.g. "C1, C2, C3, C4")
    origin_info <- rlang::list2(!!SettingsInfo[['GroupingVariable']] := InputData[[SettingsInfo[['GroupingVariable']]]][i], !!SettingsInfo[['OriginalID']] := InputData[[SettingsInfo[['OriginalID']]]][i])

      if (length(InputData[[SettingsInfo[['TranslatedID']]]][i]) == 0 || is.null(InputData[[SettingsInfo[['TranslatedID']]]][i]) || is.na(InputData[[SettingsInfo[['TranslatedID']]]][i])) {
        processed_data <- (list(None = origin_info))
      } else {
        # Split the string by commas and return a named list
        ids <- strsplit(InputData[[SettingsInfo[['TranslatedID']]]][i], ",\\s*")[[1]]
        processed_data <-(setNames(rep(list(origin_info), length(ids)), ids))
      }

     # Append data to the mapping list
     for (item in names(processed_data)) {
        if (is.null(mapping_list[[item]])) {
          mapping_list[[item]] <- list(processed_data[[item]])
        } else {
          mapping_list[[item]] <- append(mapping_list[[item]], list(processed_data[[item]]))
        }
      }
  }

  #----- Convert the mapping list into a data frame
  rows <- list()  # Initialize a list to store the rows

  # Iterate through the mapping list and extract data
  for (key in names(mapping_list)) {
    for (mapping in mapping_list[[key]]) {
      # Ensure that mapping is not NULL and contains term/id
      if (!is.null(mapping) && !is.null(mapping[[SettingsInfo[['GroupingVariable']]]]) && !is.null(mapping[[SettingsInfo[['OriginalID']]]])) {
        # Append each row as a list to the rows list
        rows[[length(rows) + 1]] <- rlang::list2(!!SettingsInfo[['TranslatedID']] := key, term = mapping[[SettingsInfo[['GroupingVariable']]]], !!SettingsInfo[['OriginalID']] := mapping[[SettingsInfo[['OriginalID']]]])
      }
    }
  }

  # Convert the list of rows into a data frame
  mapping_df <- dplyr::bind_rows(rows)

  ## ------------------ Count occurences of TranslatedID-to-OriginalID ------------------- ##
  #----- Count how many unique terms and ids are associated with each ID
  mapping_df_summary <- mapping_df %>%
    dplyr::group_by(!!rlang::sym(SettingsInfo[['TranslatedID']])) %>%
    dplyr::summarize(!!paste0(SettingsInfo[['GroupingVariable']], "_count", sep="") := as.numeric(dplyr::n_distinct(!!rlang::sym(SettingsInfo[['GroupingVariable']]))), !!paste0(SettingsInfo[['OriginalID']], "_count", sep="") := as.numeric(dplyr::n_distinct(!!rlang::sym(SettingsInfo[['OriginalID']])))) %>%
    dplyr::filter(paste0(SettingsInfo[['GroupingVariable']], "_count", sep="") > 0 | paste0(SettingsInfo[['OriginalID']], "_count", sep="") > 0) # was previously filter(GroupingVariable_count > 1 | ID_original_count > 1)

  ### Nested results
  mapping_df_nested <- mapping_df %>%
    dplyr::group_by(!!rlang::sym(SettingsInfo[['TranslatedID']])) %>%
    dplyr::summarize(!!SettingsInfo[['GroupingVariable']] := toString(unique(!!sym(SettingsInfo[['GroupingVariable']]))), !!SettingsInfo[['OriginalID']] := toString(unique(!!sym(SettingsInfo[['OriginalID']]))))

  # Join the summary with the nested table for easier user inspection
  mapping_df_joined <- mapping_df_nested %>%
    dplyr::left_join(mapping_df_summary, by = SettingsInfo[['TranslatedID']])%>%
    dplyr::mutate(Relationship = dplyr::case_when(!!rlang::sym(paste0(SettingsInfo[['OriginalID']], "_count", sep="")) == 0 ~ "One-to-None",
                                           !!rlang::sym(paste0(SettingsInfo[['OriginalID']], "_count", sep="")) == 1 ~ "One-to-One",
                                           !!rlang::sym(paste0(SettingsInfo[['OriginalID']], "_count", sep="")) > 1 ~ "One-to-Many",
                                           TRUE ~ "FALSE"))

  #Deal with NAs:
  #It is to be expected to have NA's occuring between the original and translated IDs, as not all IDs will be able to be translated. We will remove these from the count.


  # Now rename the columns so the user knows what exactly the Original and Translated IDs are for each DF...
  mapping_inspection$Orig2Trans <- mapping_inspection$Orig2Trans %>%
    dyplyr::rename(!!paste0("Original_ID_", From) := "Original_ID", !!paste0("Translated_IDs_", MetaboliteID) := "Translated_IDs")
  mapping_inspection$Trans2Orig <- mapping_inspection$Trans2Orig %>%
    dyplyr::rename(!!paste0("Original_IDs_", From) := "Original_IDs", !!paste0("Translated_ID_", MetaboliteID) := "Translated_ID")


  ## ------------------ Output and Return ------------------- ##
  ### Output
  Mapping_Trans2Orig  <- list(mapping = mapping_df,
                                        mapping_summary = mapping_df_summary,
                                        mapping_nested = mapping_df_nested,
                                        mapping_joined= mapping_df_joined)

  Res <- c(
    list("Orig2Trans" = NULL),
    setNames(list(mapping_df_joined), paste0(SettingsInfo[['TranslatedID']], "-to-", SettingsInfo[['OriginalID']]))
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
#'



