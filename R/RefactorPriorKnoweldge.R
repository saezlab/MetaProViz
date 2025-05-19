## ---------------------------
##
## Script name: Getprior_knowledge
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
#' @param data dataframe with at least one column with the target (e.g. metabolite), you can add other columns such as source (e.g. term). Must be "long" DF, meaning one ID per row.
#' @param metadata_info \emph{Optional: } Column name of Target in input_pk. \strong{Default = list(InputID="MetaboliteID" , grouping_variable="term")}
#' @param from ID type that is present in your data. Choose between "kegg", "pubchem", "chebi", "hmdb". \strong{Default = "kegg"}
#' @param to One or multiple ID types to which you want to translate your data. Choose between "kegg", "pubchem", "chebi", "hmdb". \strong{Default = c("pubchem","chebi","hmdb")}
#' @param summary \emph{Optional: } If TRUE a long summary tables are created. \strong{Default = FALSE}
#' @param save_table \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt". \strong{Default = "csv"}
#' @param path {Optional:} Path to the folder the results should be saved at. \strong{Default = NULL}
#'
#' @return List with at least three DFs: 1) Original data and the new column of translated ids spearated by comma. 2) Mapping information between Original ID to Translated ID. 3) Mapping summary between Original ID to Translated ID.
#'
#' @examples
#' KEGG_Pathways <- MetaProViz::metsigdb_kegg()
#' Res <- MetaProViz::translate_id(data= KEGG_Pathways, metadata_info = c(InputID="MetaboliteID", grouping_variable="term"), from = c("kegg"), to = c("pubchem", "hmdb"))
#'
#' @keywords Translate metabolite IDs
#'
#' @importFrom dplyr mutate select group_by ungroup distinct filter across n
#' @importFrom tidyselect all_of starts_with
#' @importFrom rlang !!! !! := sym syms
#' @importFrom OmnipathR id_types translate_ids
#' @importFrom logger log_warn
#' @importFrom stringr str_to_lower
#'
#' @export
#'
translate_id <- function(data,
                        metadata_info = c(InputID="MetaboliteID", grouping_variable="term"),
                        from = "kegg",
                        to = c("pubchem","chebi","hmdb"),
                        summary=FALSE,
                        save_table= "csv",
                        path=NULL,
                        plot=FALSE #toDO
  ){# Add ability to also get metabolite names that are human readable from an ID type!

  metaproviz_init()

  ## ------------------  Check Input ------------------- ##
  # HelperFunction `check_param`
  check_param(data=data,
                          data_num=FALSE,
                          save_table=save_table)

  # Specific checks:
  if("InputID" %in% names(metadata_info)){
    if(metadata_info[["InputID"]] %in% colnames(data)== FALSE){
      message <- paste0("The ", metadata_info[["InputID"]], " column selected as InputID in metadata_info was not found in data. Please check your input.")
      logger::log_trace(paste("Error ", message, sep=""))
      stop(message)
    }
  }

  if("grouping_variable" %in% names(metadata_info)){
    if(metadata_info[["grouping_variable"]] %in% colnames(data)== FALSE){
      message <- paste0("The ", metadata_info[["grouping_variable"]], " column selected as grouping_variable in metadata_info was not found in data. Please check your input.")
      logger::log_trace(paste("Error ", message, sep=""))
      stop(message)
    }
  }

  if(is.logical(summary) == FALSE){
    message <- paste0("Check input. The summary parameter should be either =TRUE or =FALSE.")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }

  unknown_types <- OmnipathR::id_types() %>%
    dplyr::select(tidyselect::starts_with('in_')) %>%
    unlist %>%
    unique %>%
    str_to_lower %>%
    setdiff(union(from, to), .)

  if (length(unknown_types) > 0L) {
    msg <- sprintf('The following ID types are not recognized: %s', paste(unknown_types, collapse = ', '))
    logger::log_warn(msg)
    warning(msg)
  }

  # Check that metadata_info[['InputID']] has no duplications within one group --> should not be the case --> remove duplications and inform the user/ ask if they forget to set groupings column
  doublons <- data %>%
    dplyr::filter(!is.na(!!sym(metadata_info[['InputID']]))) %>%
    dplyr::group_by(!!sym(metadata_info[['InputID']]), !!sym(metadata_info[['grouping_variable']]))%>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::summarize()

  if(nrow(doublons) > 0){
    message <- sprintf(
        'The following IDs are duplicated within one group: %s',
        paste(doublons %>% dplyr::pull(metadata_info[['InputID']]), collapse = ', ')
    )
    logger::log_warn(message)
    warning(message)
  }

  ## ------------------  Create output folders and path ------------------- ##
  if(is.null(save_table)==FALSE ){
    folder <- save_path(folder_name= "prior_knowledge",
                                    path=path)

    Subfolder <- file.path(folder, "ID_Translation")
    if (!dir.exists(Subfolder)) {dir.create(Subfolder)}
  }

  ########################################################################################################################################################
  ## ------------------ Translate to-from for each pair ------------------- ##
  TranslatedDF <- OmnipathR::translate_ids(
      data,
      !!sym(metadata_info[['InputID']]) :=  !!sym(from),
      !!!syms(to),#list of symbols, hence three !!!
      ramp = TRUE,
      expand = FALSE,
      quantify_ambiguity = TRUE,
      qualify_ambiguity = TRUE,
      ambiguity_groups =  metadata_info[['grouping_variable']],#Checks within the groups, without it checks across groups
      ambiguity_summary = TRUE
    )
  #TranslatedDF %>% attributes %>% names
  #TranslatedDF%>% attr('ambiguity_MetaboliteID_hmdb')

  ## --------------- Create output DF -------------------- ##
  ResList <- list()

  ## Create DF for TranslatedIDs only with the original data and the translatedID columns
  DF_subset <- TranslatedDF %>%
    dplyr::select(tidyselect::all_of(intersect(names(.), names(data))), tidyselect::all_of(to)) %>%
    dplyr::mutate(across(all_of(to), ~ map_chr(., ~ paste(unique(.), collapse = ", ")))) %>%
    dplyr::group_by(!!sym(metadata_info[['InputID']]), !!sym(metadata_info[['grouping_variable']])) %>%
    dplyr::mutate(across(tidyselect::all_of(to), ~ paste(unique(.), collapse = ", "), .names = "{.col}")) %>%
    dplyr::ungroup() %>%
    dplyr::distinct() %>%
    dplyr::mutate(dplyr::across(tidyselect::all_of(to), ~ ifelse(. == "0", NA, .)))

  ResList[["TranslatedDF"]] <- DF_subset

  ## Add DF with mapping information
  ResList[["TranslatedDF_MappingInfo"]] <- TranslatedDF

  ## Also save the different mapping summaries!
  for(item in to){
    summaryDF <- TranslatedDF%>% attr(paste0("ambiguity_", metadata_info[['InputID']], "_", item, sep=""))
    ResList[[paste0("Mappingsummary_", item, sep="")]] <-  summaryDF
  }

  ## Create the long DF summary if summary =TRUE
  if(summary==TRUE){
    for(item in to){
      summary <- MetaProViz::mapping_ambiguity(data= TranslatedDF,
                                            from = metadata_info[['InputID']],
                                            to = item,
                                            grouping_variable = metadata_info[['grouping_variable']],
                                            summary=TRUE)[["summary"]]
      ResList[[paste0("Mappingsummary_Long_", from, "-to-", item, sep="")]] <- summary
    }
  }


  ## ------------------ Save the results ------------------- ##
  suppressMessages(suppressWarnings(
    save_res(inputlist_df=ResList,
                         inputlist_plot= NULL,
                         save_table=save_table,
                         save_plot=NULL,
                         path= Subfolder,
                         file_name= "translate_id",
                         core=FALSE,
                         print_plot=FALSE)))

  ## ------------------ Show Upset plot of the results ------------------- ##
  ## toDO: current issue with Lists vs doubles
  # if (plot==TRUE) {
  #   pk_list <- list(translated = TranslatedDF)
  #   metadata_info <- list(translated = to)
  #   print(pk_list)
  #   print(to)
  #   pk_comp_res <- MetaProViz:::compare_pk(data = pk_list,
  #                                         metadata_info = metadata_info,
  #                                         plot_title = "IDs available after ID Translation")
  #   ## Add Upset plot
  #   ResList[["Translated_UpsetPlot"]] <- pk_comp_res$upset_plot
  #}

  #Return
  invisible(return(ResList))
}



##########################################################################################
### ### ### Find additional potential IDs  ### ### ###
##########################################################################################

#' Find additional potential IDs for  "kegg", "pubchem", "chebi", "hmdb"
#'
#' @param data dataframe with at least one column with the detected metabolite IDs (one ID per row).
#' @param metadata_info \emph{Optional: } Column name of metabolite IDs. \strong{Default = list(InputID="MetaboliteID")}
#' @param from ID type that is present in your data. Choose between "kegg", "pubchem", "chebi", "hmdb". \strong{Default = "hmdb"}
#' @param save_table \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt". \strong{Default = "csv"}
#' @param path {Optional:} Path to the folder the results should be saved at. \strong{Default = NULL}
#'
#' @return Input DF with additional column including potential additional IDs.
#'
#' @examples
#' DetectedIDs <- MetaProViz::ToyData(data="Cells_Metadata")%>% tibble::rownames_to_column("TrivialName")%>%tidyr::drop_na()
#' Res <- MetaProViz::equivalent_ids(data= DetectedIDs, metadata_info = c(InputID="HMDB"), from = "hmdb")
#'
#' @keywords Find potential additional IDs for one metabolite identifier
#'
#' @importFrom dplyr mutate select group_by ungroup distinct filter across rowwise
#' @importFrom tidyr separate_rows unnest
#' @importFrom purrr map_chr
#' @importFrom tidyselect all_of starts_with
#' @importFrom rlang !!! !! := sym syms
#' @importFrom OmnipathR id_types translate_ids
#' @importFrom logger log_warn log_trace
#' @importFrom stringr str_to_lower str_split
#' @export
equivalent_ids <- function(data,
                          metadata_info = c(InputID="MetaboliteID"),
                          from = "hmdb",
                          save_table= "csv",
                          path=NULL){
  # FUTURE: Once we have the structural similarity tool available in OmniPath, we can start creating this function!

  ### 1)
  #Check Measured ID's in prior knowledge


  ### 2)
  # A user has one HMDB IDs for their measured metabolites (one ID per measured peak) --> this is often the case as the user either gets a trivial name and they have searched for the ID themselves or because the facility only provides one ID at random
  # We have mapped the HMDB IDs with the pathways and 20 do not map
  # We want to check if it is because the pathways don't include them, or because the user just gave the wrong ID by chance (i.e. They picked D-Alanine, but the prior knowledge includes L-Alanine)

  # Do this by using structural information via  accessing the structural DB in OmniPath!
  # Output is DF with the original ID column and a new column with additional possible IDs based on structure

  #Is it possible to do this at the moment without structures, but by using other pior knowledge?

  metaproviz_init()

  ## ------------------  Check Input ------------------- ##
  # HelperFunction `check_param`
  check_param(data=data,
                          data_num=FALSE,
                          save_table=save_table)

  # Specific checks:
  if("InputID" %in% names(metadata_info)){
    if(metadata_info[["InputID"]] %in% colnames(data)== FALSE){
      message <- paste0("The ", metadata_info[["InputID"]], " column selected as InputID in metadata_info was not found in data. Please check your input.")
      logger::log_trace(paste("Error ", message, sep=""))
      stop(message)
    }
  }

  unknown_types <- OmnipathR::id_types() %>%
    dplyr::select(tidyselect::starts_with('in_')) %>%
    unlist %>%
    unique %>%
    str_to_lower %>%
    setdiff(from, .)

  if (length(unknown_types) > 0L) {
    msg <- sprintf('The following ID types are not recognized: %s', paste(unknown_types, collapse = ', '))
    logger::log_warn(msg)
    warning(msg)
  }

  # Check that metadata_info[['InputID']] has no duplications within one group --> should not be the case --> remove duplications and inform the user/ ask if they forget to set groupings column
  doublons <- data[duplicated(data[[metadata_info[['InputID']]]]), ]

  if(nrow(doublons) > 0){
    data <- data %>%
      dplyr::distinct(!!sym(metadata_info[['InputID']]), .keep_all = TRUE)

    message <- sprintf('The following IDs are duplicated and removed: %s',paste(doublons[[metadata_info[['InputID']]]], collapse = ', '))
    logger::log_warn(message)
    warning(message)
  }

  ## ------------------  Create output folders and path ------------------- ##
  if(is.null(save_table)==FALSE ){
    folder <- save_path(folder_name= "prior_knowledge",
                                    path=path)

    Subfolder <- file.path(folder, "equivalent_ids")
    if (!dir.exists(Subfolder)) {dir.create(Subfolder)}
  }

  ## ------------------ Set the ID type for to ----------------- ##
  to <- case_when(
    from == "chebi" ~ "pubchem",  # If to is "pubchem", choose "chebi"
    TRUE ~ "chebi"              # For other cases, don't use a secondary column
  )

  message <- paste0(to, " is used to find additional potential IDs for ", from, ".", sep="")
  logger::log_trace(message)
  message(message)

  ## ------------------ Load manual table ----------------- ##
  if((from == "kegg") == FALSE){
    EquivalentFeatures <- equivalent_features%>%
      dplyr::select(from)
  }

  ## ------------------ Translate from-to-to ------------------- ##
  TranslatedDF <- OmnipathR::translate_ids(
    data,
    !!sym(metadata_info[['InputID']]) :=  !!sym(from),
    !!!syms(to),#list of symbols, hence three !!!
    ramp = TRUE,
    expand = FALSE,
    quantify_ambiguity =FALSE,
    qualify_ambiguity =  TRUE, # Can not be set to FALSE!
    ambiguity_groups =  NULL,#Checks within the groups, without it checks across groups
    ambiguity_summary =  FALSE
  )%>%
    dplyr::select(tidyselect::all_of(intersect(names(.), names(data))), tidyselect::all_of(to)) %>%
    dplyr::mutate(across(all_of(to), ~ purrr::map_chr(., ~ paste(unique(.), collapse = ", ")))) %>%
    dplyr::group_by(!!sym(metadata_info[['InputID']])) %>%
    dplyr::mutate(across(tidyselect::all_of(to), ~ paste(unique(.), collapse = ", "), .names = "{.col}")) %>%
    dplyr::ungroup() %>%
    dplyr::distinct() %>%
    dplyr::mutate(dplyr::across(tidyselect::all_of(to), ~ ifelse(. == "0", NA, .)))


  ## ------------------ Translate to-to-from ------------------- ##
  TranslatedDF_Long <- TranslatedDF%>%
    dplyr::select(!!sym(metadata_info[['InputID']]), !!sym(to))%>%
    dplyr::rename("InputID" = !!sym(metadata_info[['InputID']]))%>%
    tidyr::separate_rows(!!sym(to), sep = ", ") %>%
    dplyr::mutate(across(all_of(to), ~trimws(.))) %>%  # Remove extra spaces
    dplyr::filter(!!sym(to) != "")  # Remove empty entries

  OtherIDs <- OmnipathR::translate_ids(
    TranslatedDF_Long ,
    !!sym(to),
    !!sym(from),#list of symbols, hence three !!!
    ramp = TRUE,
    expand = FALSE,
    quantify_ambiguity =FALSE,
    qualify_ambiguity =  TRUE, # Can not be set to FALSE!
    ambiguity_groups =  NULL,#Checks within the groups, without it checks across groups
    ambiguity_summary =  FALSE
  )%>%
    dplyr::select("InputID", !!sym(to), !!sym(from))%>%
    dplyr::distinct(InputID, !!sym(from), .keep_all = TRUE) %>%  # Remove duplicates based on InputID and from
    dplyr::mutate(AdditionalID = dplyr::if_else(InputID == !!sym(from), FALSE, TRUE)) %>%
    dplyr::select("InputID",!!sym(from), "AdditionalID")%>%
    dplyr::filter(AdditionalID == TRUE) %>%
    dplyr::mutate(across(all_of(from), ~ purrr::map_chr(., ~ paste(unique(.), collapse = ", "))))%>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      fromList = list(stringr::str_split(!!sym(from), ",\\s*")[[1]]),  # Wrap in list
      SameAsInput = ifelse(any(fromList == InputID), InputID, NA_character_),  # Match InputID
      PotentialAdditionalIDs = paste(fromList[fromList != InputID], collapse = ", ")  # Combine other IDs
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(InputID, PotentialAdditionalIDs, hmdb)%>%  # Final selection
    dplyr::rename("AllIDs"= "hmdb")

  ## ------------------ Merge to Input ------------------- ##
  OtherIDs <- merge(data, OtherIDs, by.x= metadata_info[['InputID']] , by.y= "InputID", all.x=TRUE)

  ##------------------- Add additional IDs -------------- ##

  if (exists("EquivalentFeatures")) {
   EquivalentFeatures$AllIDs <- EquivalentFeatures[[from]]
    EquivalentFeatures_Long <- EquivalentFeatures  %>%
      separate_rows(!!sym(from), sep = ",")

    OtherIDs <- merge(OtherIDs, EquivalentFeatures_Long, by.x= metadata_info[['InputID']] , by.y= "hmdb", all.x=TRUE)%>%
      rowwise() %>%
      mutate(AllIDs = paste(unique(na.omit(unlist(stringr::str_split(paste(na.omit(c(AllIDs.x, AllIDs.y)), collapse = ","), ",\\s*")))), collapse = ",")) %>%
      ungroup()%>%
      rowwise() %>%
      mutate(
        PotentialAdditionalIDs = paste(
          setdiff(
            unlist(stringr::str_split(AllIDs, ",\\s*")),  # Split merged_column into individual IDs
            as.character(!!sym(metadata_info[['InputID']]))  # Split hmdb into individual IDs
          ),
          collapse = ", "  # Combine the remaining IDs back into a comma-separated string
        )
      ) %>%
      ungroup()%>%
      select(-AllIDs.x, -AllIDs.y)
  }

  ##------------------- Fill empty cells -------------- ##
  OtherIDs <- OtherIDs%>%
    mutate(PotentialAdditionalIDs = ifelse(PotentialAdditionalIDs == "", NA, PotentialAdditionalIDs))%>%
    mutate(AllIDs = ifelse(AllIDs == "", !!sym(metadata_info[['InputID']]), AllIDs))

  ## ------------------ Create count_ids plot ------------------- ##
  #QC plot of before and after
  Before <-MetaProViz:::count_ids(data=data,
                                 column=metadata_info[["InputID"]],
                                 save_plot = NULL,
                                 save_table = NULL,
                                 print_plot = FALSE,
                                 path = NULL)

  After <-MetaProViz:::count_ids(data=OtherIDs,
                                  column="AllIDs",
                                  save_plot = NULL,
                                  save_table = NULL,
                                  print_plot = FALSE,
                                  path = NULL)


  ## ------------------ Create Output ------------------- ##
  OutputDF <- OtherIDs

  ## ------------------ Save the results ------------------- ##
  ResList <- list("equivalent_ids" = OutputDF)

  suppressMessages(suppressWarnings(
    save_res(inputlist_df=ResList,
                         inputlist_plot= NULL,
                         save_table=save_table,
                         save_plot=NULL,
                         path= Subfolder,
                         file_name= "equivalent_ids",
                         core=FALSE,
                         print_plot=FALSE)))

  return(invisible(OutputDF))
}

##########################################################################################
### ### ### Mapping Ambiguity ### ### ###
##########################################################################################

#' Create Mapping Ambiguities between two ID types
#'
#' @param data Translated DF from MetaProViz::translate_id reults or dataframe with at least one column with the target metabolite ID and another MetaboliteID type. One of the IDs can only have one ID per row, the other ID can be either separated by comma or a list. Optional: add other columns such as source (e.g. term).
#' @param to Column name of original metabolite identifier in data. Here should only have one ID per row.
#' @param from Column name of the secondary or translated metabolite identifier in data. Here can be multiple IDs per row either separated by comma " ," or a list of IDs.
#' @param grouping_variable \emph{Optional: } If NULL no groups are used. If TRUE provide column name in data containing the grouping_variable and features are grouped. \strong{Default = NULL}
#' @param summary \emph{Optional: } If TRUE a long summary tables are created. \strong{Default = FALSE}
#' @param save_table \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt". \strong{Default = "csv"}
#' @param path {Optional:} Path to the folder the results should be saved at. \strong{Default = NULL}
#'
#' @return List with at least 4 DFs: 1-3) from-to-to: 1. MappingIssues, 2. MappingIssues summary, 3. Long summary (If summary=TRUE) & 4-6) to-to-from: 4. MappingIssues, 5. MappingIssues summary, 6. Long summary (If summary=TRUE) & 7) Combined summary table (If summary=TRUE)
#'
#' @examples
#' KEGG_Pathways <- MetaProViz::metsigdb_kegg()
#' InputDF <- MetaProViz::translate_id(data= KEGG_Pathways, metadata_info = c(InputID="MetaboliteID", grouping_variable="term"), from = c("kegg"), to = c("pubchem"))[["TranslatedDF"]]
#' Res <- MetaProViz::mapping_ambiguity(data= InputDF, from = "MetaboliteID", to = "pubchem", grouping_variable = "term", summary=TRUE)
#'
#' @keywords Mapping ambiguity
#'
#' @importFrom dplyr mutate bind_cols bind_rows
#' @importFrom rlang !!! !! := sym syms
#' @importFrom OmnipathR ambiguity
#'
#' @export
#'
mapping_ambiguity <- function(data,
                             from,
                             to,
                             grouping_variable = NULL,
                             summary=FALSE,
                             save_table= "csv",
                             path=NULL
) {

  metaproviz_init()
  ## ------------------  Check Input ------------------- ##
  # HelperFunction `check_param`
  check_param(data=data,
                          data_num=FALSE,
                          save_table=save_table)

  # Specific checks:
  if(from %in% colnames(data)== FALSE){
      message <- paste0(from, " column was not found in data. Please check your input.")
      logger::log_trace(paste("Error ", message, sep=""))
      stop(message)
      }

  if(to %in% colnames(data)== FALSE){
    message <- paste0(to, " column was not found in data. Please check your input.")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }

  if(is.null(grouping_variable)==FALSE){
    if(grouping_variable %in% colnames(data)== FALSE){
      message <- paste0(grouping_variable, " column was not found in data. Please check your input.")
      logger::log_trace(paste("Error ", message, sep=""))
      stop(message)
    }
  }

  if(is.logical(summary) == FALSE){
    message <- paste0("Check input. The summary parameter should be either =TRUE or =FALSE.")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }

  ## ------------------  General checks of wrong occurences ------------------- ##
  # Task 1: Check that from has no duplications within one group --> should not be the case --> remove duplications and inform the user/ ask if they forget to set groupings column
  # Task 2: Check that from has the same items in to across the different entries (would be in different Groupings, otherwise there should not be any duplications) --> List of Miss-Mappings across terms

  # FYI: The above can not happen if our translateID function was used, but may be the case when the user has done something manually before


  ## ------------------  Create output folders and path ------------------- ##
  if(is.null(save_table)==FALSE ){
    folder <- save_path(folder_name= "prior_knowledge",
                                    path=path)

    Subfolder <- file.path(folder, "MappingAmbiguities")
    if (!dir.exists(Subfolder)) {dir.create(Subfolder)}
  }

  #####################################################################################################################################################################################
  ## ------------------  Prepare Input data ------------------- ##
  #If the user provides a DF where the to column is a list of IDs, then we can use it right away
  #If the to column is not a list of IDs, but a character column, we need to convert it into a list of IDs
  if(is.character(data[[to]])==TRUE){
    data[[to]] <- data[[to]]%>%
      strsplit(", ")%>%
      lapply(as.character)
  }

  ## ------------------  Perform ambiguity mapping ------------------- ##
  #1. from-to-to: OriginalID-to-TranslatedID
  #2. from-to-to: TranslatedID-to-OriginalID
  Comp <- list(
    list(from = from, to = to),
    list(from = to, to = from)
  )

  ResList <- list()
  for(comp in seq_along(Comp)){
    #Run Omnipath ambiguity
    ResList[[paste0(Comp[[comp]]$from, "-to-", Comp[[comp]]$to , sep="")]] <- data %>%
      tidyr::unnest(cols = all_of(Comp[[comp]]$from))%>% # unlist the columns in case they are not expaned
      filter(!is.na(!!sym(Comp[[comp]]$from)))%>%#Remove NA values, otherwise they are counted as column is character
      OmnipathR::ambiguity(
          from_col = !!sym(Comp[[comp]]$from),
          to_col = !!sym(Comp[[comp]]$to),
          groups = grouping_variable,
          quantify = TRUE,
          qualify = TRUE,
          global = TRUE,#across groups will be done additionally --> suffix _AcrossGroup
          summary=TRUE, #summary of the mapping column
          expand=TRUE)

    #Extract summary table:
    ResList[[paste0(Comp[[comp]]$from, "-to-", Comp[[comp]]$to, "_summary", sep="")]] <-
        ResList[[paste0(Comp[[comp]]$from, "-to-", Comp[[comp]]$to , sep="")]]%>%
        attr(paste0("ambiguity_", Comp[[comp]]$from , "_",Comp[[comp]]$to, sep=""))

    ############################################################################################################
    if(summary==TRUE){
      if(is.null(grouping_variable)==FALSE){
        # Add further information we need to summarise the table and combine Original-to-Translated and Translated-to-Original
        # If we have a grouping_variable we need to combine it with the MetaboliteID before merging
        ResList[[paste0(Comp[[comp]]$from, "-to-", Comp[[comp]]$to, "_Long", sep="")]] <- ResList[[paste0(Comp[[comp]]$from, "-to-", Comp[[comp]]$to , sep="")]]%>%
          tidyr::unnest(cols = all_of(Comp[[comp]]$from))%>%
          mutate(!!sym(paste0("AcrossGroupMappingIssue(", Comp[[comp]]$from, "_to_", Comp[[comp]]$to, ")", sep="")) := case_when(
            !!sym(paste0(Comp[[comp]]$from, "_", Comp[[comp]]$to, "_ambiguity_bygroup", sep="")) != !!sym(paste0(Comp[[comp]]$from, "_", Comp[[comp]]$to, "_ambiguity", sep=""))  ~ "TRUE",
            TRUE ~ "FALSE" ))%>%
          group_by(!!sym(Comp[[comp]]$from), !!sym(grouping_variable))%>%
          mutate(!!sym(Comp[[comp]]$to) := ifelse(!!sym(Comp[[comp]]$from) == 0, NA,  # Or another placeholder
                                                  paste(unique(!!sym(Comp[[comp]]$to)), collapse = ", ")
          )) %>%
          mutate( !!sym(paste0("Count(", Comp[[comp]]$from, "_to_", Comp[[comp]]$to, ")")) := ifelse(all(!!sym(Comp[[comp]]$to) == 0), 0, n()))%>%
          ungroup()%>%
          distinct() %>%
          unite(!!sym(paste0(Comp[[comp]]$from, "_to_", Comp[[comp]]$to)), c(Comp[[comp]]$from, Comp[[comp]]$to), sep=" --> ", remove=FALSE)%>%
          separate_rows(!!sym(Comp[[comp]]$to), sep = ", ") %>%
          unite(UniqueID, c(from, to, grouping_variable), sep="_", remove=FALSE)%>%
          distinct()
        }else{
          ResList[[paste0(Comp[[comp]]$from, "-to-", Comp[[comp]]$to, "_Long", sep="")]] <- ResList[[paste0(Comp[[comp]]$from, "-to-", Comp[[comp]]$to , sep="")]]%>%
            tidyr::unnest(cols = all_of(Comp[[comp]]$from))%>%
            group_by(!!sym(Comp[[comp]]$from))%>%
            mutate(!!sym(Comp[[comp]]$to) := ifelse(!!sym(Comp[[comp]]$from) == 0, NA,  # Or another placeholder
                                                    paste(unique(!!sym(Comp[[comp]]$to)), collapse = ", ")
            )) %>%
            mutate( !!sym(paste0("Count(", Comp[[comp]]$from, "_to_", Comp[[comp]]$to, ")")) := ifelse(all(!!sym(Comp[[comp]]$to) == 0), 0, n()))%>%
            ungroup()%>%
            distinct() %>%
            unite(!!sym(paste0(Comp[[comp]]$from, "_to_", Comp[[comp]]$to)), c(Comp[[comp]]$from, Comp[[comp]]$to), sep=" --> ", remove=FALSE)%>%
            separate_rows(!!sym(Comp[[comp]]$to), sep = ", ") %>%
            unite(UniqueID, c(from, to), sep="_", remove=FALSE)%>%
            distinct()%>%
            mutate(!!sym(paste0("AcrossGroupMappingIssue(", from, "_to_", to, ")", sep="")) := NA)
        }
    }

    # Add NA metabolite maps back if they do exist:
    Removed <- data %>%
      tidyr::unnest(cols = all_of(Comp[[comp]]$from))%>% # unlist the columns in case they are not expaned
      filter(is.na(!!sym(Comp[[comp]]$from)))
    if(nrow(Removed)>0){
      ResList[[paste0(Comp[[comp]]$from, "-to-", Comp[[comp]]$to , sep="")]] <- dplyr::bind_rows(ResList[[paste0(Comp[[comp]]$from, "-to-", Comp[[comp]]$to , sep="")]],
                                                                                         test<- Removed%>%
                                                                                            dplyr::bind_cols(setNames(as.list(rep(NA, length(setdiff(names(ResList[[paste0(Comp[[comp]]$from, "-to-", Comp[[comp]]$to , sep="")]]), names(Removed))))),
                                                                                                               setdiff(names(ResList[[paste0(Comp[[comp]]$from, "-to-", Comp[[comp]]$to , sep="")]]), names(Removed))))
      )
    }
  }

  ## ------------------ Create summaryTable ------------------- ##
  if(summary==TRUE){
    # Combine the two tables
    summary <- merge(x= ResList[[paste0(from, "-to-", to, "_Long", sep="")]][,c("UniqueID", paste0(from, "_to_", to), paste0("Count(", from, "_to_", to, ")"), paste0("AcrossGroupMappingIssue(", from, "_to_", to, ")", sep=""))],
                     y= ResList[[paste0(to, "-to-", from, "_Long", sep="")]][,c("UniqueID", paste0(to, "_to_", from), paste0("Count(", to, "_to_", from, ")"), paste0("AcrossGroupMappingIssue(", to, "_to_", from, ")", sep=""))],
                     by = "UniqueID",
                     all = TRUE)%>%
      separate(UniqueID, into = c(from, to, grouping_variable), sep="_", remove=FALSE)%>%
      distinct()

    # Add relevant mapping information
    summary <- summary %>%
      mutate(Mapping = case_when(
        !!sym(paste0("Count(", from, "_to_", to, ")")) == 1 & !!sym(paste0("Count(", to, "_to_", from, ")")) == 1  ~ "one-to-one",
        !!sym(paste0("Count(", from, "_to_", to, ")")) > 1 & !!sym(paste0("Count(", to, "_to_", from, ")")) == 1  ~ "one-to-many",
        !!sym(paste0("Count(", from, "_to_", to, ")")) > 1 & !!sym(paste0("Count(", to, "_to_", from, ")")) > 1  ~ "many-to-many",
        !!sym(paste0("Count(", from, "_to_", to, ")")) == 1 & !!sym(paste0("Count(", to, "_to_", from, ")")) > 1  ~ "many-to-one",
        !!sym(paste0("Count(", from, "_to_", to, ")")) >= 1 & !!sym(paste0("Count(", to, "_to_", from, ")")) == NA  ~ "one-to-none",
        !!sym(paste0("Count(", from, "_to_", to, ")")) >= 1 & is.na(!!sym(paste0("Count(", to, "_to_", from, ")")))  ~ "one-to-none",
        !!sym(paste0("Count(", from, "_to_", to, ")")) == NA & !!sym(paste0("Count(", to, "_to_", from, ")")) >= 1  ~ "none-to-one",
        is.na(!!sym(paste0("Count(", from, "_to_", to, ")"))) & !!sym(paste0("Count(", to, "_to_", from, ")")) >= 1  ~ "none-to-one",
        TRUE ~ NA )) %>%
      mutate( !!sym(paste0("Count(", from, "_to_", to, ")")) := replace_na( !!sym(paste0("Count(", from, "_to_", to, ")")), 0)) %>%
      mutate( !!sym(paste0("Count(", to, "_to_", from, ")")) := replace_na( !!sym(paste0("Count(", to, "_to_", from, ")")), 0))

    ResList[["summary"]] <- summary
  }

  ## ------------------ Save the results ------------------- ##
  suppressMessages(suppressWarnings(
    save_res(inputlist_df=ResList,
                         inputlist_plot= NULL,
                         save_table=save_table,
                         save_plot=NULL,
                         path= Subfolder,
                         file_name= "mapping_ambiguity",
                         core=FALSE,
                         print_plot=FALSE)))

  #Return
  invisible(return(ResList))
}

##########################################################################################
### ### ### Check Measured ID's in prior knowledge ### ### ###
##########################################################################################

#' Check and summarize prior_knowledge-to-MeasuredFeatures relationship
#'
#' @param data dataframe with at least one column with the detected metabolite IDs (e.g. HMDB). If there are multiple IDs per detected peak, please separate them by comma ("," or ", " or chr list). If there is a main ID and additional IDs, please provide them in separate columns.
#' @param prior_knowledge dataframe with at least one column with the metabolite ID (e.g. HMDB) that need to match data metabolite IDs "source" (e.g. term). If there are multiple IDs, as the original pathway IDs (e.g. KEGG) where translated (e.g. to HMDB), please separate them by comma ("," or ", " or chr list).
#' @param metadata_info Colum name of Metabolite IDs in data and prior_knowledge as well as column name of grouping_variable in prior_knowledge. \strong{Default = c(InputID="HMDB", PriorID="HMDB", grouping_variable="term")}
#' @param save_table \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt". \strong{Default = "csv"}
#' @param path {Optional:} Path to the folder the results should be saved at. \strong{Default = NULL}
#'
#' @importFrom dplyr mutate
#' @importFrom rlang !!! !! := sym syms
#'
#' @examples
#' DetectedIDs <-  MetaProViz::ToyData(data="Cells_Metadata")%>% rownames_to_column("Metabolite") %>%dplyr::select("Metabolite", "HMDB")%>%tidyr::drop_na()
#' input_pathway <- MetaProViz::translate_id(data= MetaProViz::metsigdb_kegg(), metadata_info = c(InputID="MetaboliteID", grouping_variable="term"), from = c("kegg"), to = c("hmdb"))[["TranslatedDF"]]%>%tidyr::drop_na()
#' Res <- MetaProViz:::check_match_id(data= DetectedIDs, prior_knowledge= input_pathway, metadata_info = c(InputID="HMDB", PriorID="hmdb", grouping_variable="term"))
#'
#' @noRd
#'

check_match_id <- function(data,
                         prior_knowledge,
                         metadata_info = c(InputID="HMDB", PriorID="HMDB", grouping_variable="term"),
                         save_table= "csv",
                         path=NULL
){

  ## ------------ Create log file ----------- ##
  metaproviz_init()

  ## ------------ Check Input files ----------- ##

  ## data:
  if("InputID" %in% names(metadata_info)){
    if(metadata_info[["InputID"]] %in% colnames(data)== FALSE){
      message <- paste0("The ", metadata_info[["InputID"]], " column selected as InpuID in metadata_info was not found in data. Please check your input.")
      logger::log_trace(paste("Error ", message, sep=""))
      stop(message)
    }
  }else{
    message <- paste0("No ", metadata_info[["InputID"]], " provided. Please check your input.")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }

  ### This is after the main input checks (before NA removal), so we will save original df here for later merging to get the Null and duplicates back.
  data_Original <- data

  if(sum(is.na(data[[metadata_info[["InputID"]]]])) >=1){#remove NAs:
     message <- paste0(sum(is.na(data[[metadata_info[["InputID"]]]])), " NA values were removed from column", metadata_info[["InputID"]])
     logger::log_trace(paste("Warning: ", message, sep=""))

     data <- data %>%
      filter(!is.na(.data[[metadata_info[["InputID"]]]]))

    warning(message)
  }

  if(nrow(data) - nrow(distinct(data, .data[[metadata_info[["InputID"]]]])) >= 1){# Remove duplicate IDs
    message <- paste0(nrow(data) - nrow(distinct(data, .data[[metadata_info[["InputID"]]]])), " duplicated IDs were removed from column", metadata_info[["InputID"]])
    logger::log_trace(paste("Warning: ", message, sep=""))

    data <- data %>%
      distinct(.data[[metadata_info[["InputID"]]]], .keep_all = TRUE)

    warning(message)
  }

  data_MultipleIDs <- any(
     grepl(",\\s*", data[[metadata_info[["InputID"]]]]) |  # Comma-separated
       sapply(data[[metadata_info[["InputID"]]]] , function(x) {
         if (grepl("^c\\(|^list\\(", x)) {
           parsed <- tryCatch(eval(parse(text = x)), error = function(e) NULL)
           return(is.list(parsed) && length(parsed) > 1 || is.vector(parsed) && length(parsed) > 1)
         }
         FALSE
       })
   )

  ## prior_knowledge:
  if("PriorID" %in% names(metadata_info)){
    if(metadata_info[["PriorID"]] %in% colnames(prior_knowledge)== FALSE){
      message <- paste0("The ", metadata_info[["PriorID"]], " column selected as InpuID in metadata_info was not found in prior_knowledge. Please check your input.")
      logger::log_trace(paste("Error ", message, sep=""))
      stop(message)
    }
  }else{
    message <- paste0("No ", metadata_info[["PriorID"]], " provided. Please check your input.")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }

  ### This is after the main input checks (before NA removal), so we will save original df here for later merging to get the Null and duplicates back.
  prior_knowledge_Original <- prior_knowledge

  if(sum(is.na(prior_knowledge[[metadata_info[["PriorID"]]]])) >=1){#remove NAs:
    message <- paste0(sum(is.na(prior_knowledge[[metadata_info[["PriorID"]]]])), " NA values were removed from column", metadata_info[["PriorID"]])
    logger::log_trace(paste("Warning: ", message, sep=""))

    prior_knowledge <- prior_knowledge %>%
      filter(!is.na(.data[[metadata_info[["PriorID"]]]]))

    warning(message)
  }

   if("grouping_variable" %in% names(metadata_info)){#Add grouping_variable
    if(metadata_info[["grouping_variable"]] %in% colnames(prior_knowledge)== FALSE){
      message <- paste0("The ", metadata_info[["grouping_variable"]], " column selected as InpuID in metadata_info was not found in prior_knowledge. Please check your input.")
      logger::log_trace(paste("Error ", message, sep=""))
      stop(message)
    }
  }else{
    #Add grouping_variable
    metadata_info["grouping_variable"] <- "grouping_variable"
    prior_knowledge["grouping_variable"] <- "None"

    message <- paste0("No ", metadata_info[["PriorID"]], " provided. If this was not intentional, please check your input.")
    logger::log_trace(message)
    message(message)
  }

  if(nrow(prior_knowledge) - nrow(distinct(prior_knowledge, .data[[metadata_info[["PriorID"]]]], .data[[metadata_info[["grouping_variable"]]]])) >= 1){# Remove duplicate IDs
    message <- paste0(nrow(prior_knowledge) - nrow(distinct(prior_knowledge, .data[[metadata_info[["PriorID"]]]], .data[[metadata_info[["grouping_variable"]]]])) , " duplicated IDs were removed from column", metadata_info[["PriorID"]])
    logger::log_trace(paste("Warning: ", message, sep=""))

    prior_knowledge <- prior_knowledge %>%
      distinct(.data[[metadata_info[["PriorID"]]]], !!sym(metadata_info[["grouping_variable"]]), .keep_all = TRUE)%>%
      group_by(!!sym(metadata_info[["PriorID"]])) %>%
      mutate(across(everything(), ~ if (is.character(.)) paste(unique(.), collapse = ", ")))%>%
      ungroup()%>%
      distinct(.data[[metadata_info[["PriorID"]]]], .keep_all = TRUE)

    warning(message)
  }

  PK_MultipleIDs <- any(# Check if multiple IDs are present:
    grepl(",\\s*", prior_knowledge[[metadata_info[["PriorID"]]]]) |  # Comma-separated
      sapply(prior_knowledge[[metadata_info[["PriorID"]]]] , function(x) {
        if (grepl("^c\\(|^list\\(", x)) {
          parsed <- tryCatch(eval(parse(text = x)), error = function(e) NULL)
          return(is.list(parsed) && length(parsed) > 1 || is.vector(parsed) && length(parsed) > 1)
        }
        FALSE
      })
  )

  ## ------------ Create Results output folder ----------- ##
  if(is.null(save_table)==FALSE){
    folder <- save_path(folder_name= "prior_knowledgeChecks",
                       path=path)
    Subfolder <- file.path(folder, "CheckMatchID_Detected-to-PK")
    if (!dir.exists(Subfolder)) {dir.create(Subfolder)}
  }

  ######################################################################################################################################
  ## ------------ Check how IDs match and if needed remove unmatched IDs ----------- ##

  # 1.  Create long DF
  create_long_df <- function(df, id_col, df_name) {
    df %>%
      mutate(row_id = dplyr::row_number()) %>%
      mutate(!!paste0("OriginalEntry_", df_name, sep="") := !!sym(id_col)) %>%  # Store original values
      separate_rows(!!sym(id_col), sep = ",\\s*") %>%
      group_by(row_id) %>%
      mutate(!!(paste0("OriginalGroup_", df_name, sep="")) := paste0(df_name, "_", dplyr::cur_group_id())) %>%
      ungroup()
  }

  if(data_MultipleIDs){
    data_long <- create_long_df(data, metadata_info[["InputID"]], "data")%>%
      select(metadata_info[["InputID"]],"OriginalEntry_data", OriginalGroup_data)
  }else{
    data_long <- data %>%
      mutate(OriginalGroup_data := paste0("data_", dplyr::row_number()))%>%
      select(metadata_info[["InputID"]], OriginalGroup_data)
  }

  if(PK_MultipleIDs){
    PK_long <- create_long_df(prior_knowledge, metadata_info[["PriorID"]], "PK")%>%
      select(metadata_info[["PriorID"]], "OriginalEntry_PK", OriginalGroup_PK, metadata_info[["grouping_variable"]])
  }else{
    PK_long <- prior_knowledge %>%
      mutate(OriginalGroup_PK := paste0("PK_", dplyr::row_number()))%>%
      select(metadata_info[["PriorID"]],OriginalGroup_PK, metadata_info[["grouping_variable"]])
  }

  # 2. Merge DF
  merged_df <- merge(PK_long, data_long, by.x= metadata_info[["PriorID"]],  by.y= metadata_info[["InputID"]], all=TRUE)%>%
    distinct(!!sym(metadata_info[["PriorID"]]), OriginalGroup_data, .keep_all = TRUE)

  #3. Add information to summarize and describe problems
  merged_df <- merged_df %>%
    # num_PK_entries
    group_by(OriginalGroup_PK, !!sym(metadata_info[["grouping_variable"]])) %>%
    mutate(
      num_PK_entries = sum(!is.na(OriginalGroup_PK)),
      num_PK_entries_groups = dplyr::n_distinct(OriginalGroup_PK, na.rm = TRUE)) %>% # count the times we have the same PK_entry match with multiple data entries --> extend below!
    ungroup()%>%
    # num_Input_entries
    group_by(OriginalGroup_data, !!sym(metadata_info[["grouping_variable"]])) %>%
    mutate(
      num_Input_entries = sum(!is.na(OriginalGroup_data)),
      num_Input_entries_groups = dplyr::n_distinct(OriginalGroup_data,, na.rm = TRUE))%>%
    ungroup()%>%
    mutate(
      ActionRequired = case_when(
        num_Input_entries ==1 & num_Input_entries_groups == 1 & num_PK_entries_groups == 1 ~ "None",
        num_Input_entries ==1 & num_Input_entries_groups == 1 & num_PK_entries_groups >= 2  ~ "Check",
        num_Input_entries ==1 & num_Input_entries_groups >= 2 & num_PK_entries_groups == 1  ~ "Check",
        num_Input_entries > 1 & num_Input_entries_groups == 1 & num_PK_entries_groups >= 2 ~ "Check",
        num_Input_entries > 1 & num_Input_entries_groups >= 2 & num_PK_entries_groups >= 2 ~ "Check",
        num_Input_entries == 0 ~ "None", # ID(s) of PK not measured
        TRUE ~ NA_character_
      )
    )%>%
    mutate(
      Detection = case_when(
        num_Input_entries ==1 & num_Input_entries_groups == 1 & num_PK_entries_groups == 1 ~ "One input ID of the same group maps to at least ONE PK ID of ONE group",
        num_Input_entries ==1 & num_Input_entries_groups == 1 & num_PK_entries_groups >= 2  ~ "One input ID of the same group maps to at least ONE PK ID of MANY groups",
        num_Input_entries ==1 & num_Input_entries_groups >= 2 & num_PK_entries_groups == 1  ~ "One input ID of MANY groups maps to at least ONE PK ID of ONE groups",
        num_Input_entries > 1 & num_Input_entries_groups == 1 & num_PK_entries_groups >= 2 ~ "MANY input IDs of the same group map to at least ONE PK ID of MANY PK groups",
        num_Input_entries > 1 & num_Input_entries_groups >= 2 & num_PK_entries_groups >= 2 ~ "MANY input IDs of the MANY groups map to at least ONE PK ID of MANY PK groups",
        num_Input_entries == 0 ~ "Not Detected", # ID(s) of PK not measured
        TRUE ~ NA_character_
      )
    )
    #  Handle "Detected-to-PK" (When PK has multiple IDs)
    #group_by(OriginalGroup_PK, !!sym(metadata_info[["grouping_variable"]])) %>%
    #mutate(
    #  `Detected-to-PK` = case_when(
    #    num_Input_entries == 0 & num_PK_entries == 1 ~ "none-to-one", # No match & OriginalGroup_PK appears once
    #    num_Input_entries == 0 & num_PK_entries > 1 ~ "none-to-many", # No match & OriginalGroup_PK appears multiple times
    #    num_Input_entries == 1 & num_PK_entries == 1 ~ "one-to-one", # One unique match in data, one PK
    #    num_Input_entries == 1 & num_PK_entries > 1 ~ "one-to-many", # One unique match in data, multiple PKs
    #    num_Input_entries > 1 & num_PK_entries == 1 ~ "many-to-one", # Multiple matches in data, one PK
    #    num_Input_entries > 1 & num_PK_entries > 1 ~ "many-to-many", # Multiple matches in data, multiple PKs
    #    num_Input_entries == 1 & num_PK_entries == 0 ~ "one-to-none",
    #    num_Input_entries > 1 & num_PK_entries == 0 ~ "many-to-none",
    #    TRUE ~ NA_character_
    #  )
    #) %>%
    #ungroup() %>%
    # Handle "PK-to-Detected" (When data has multiple IDs)
    #group_by(OriginalGroup_data, !!sym(metadata_info[["grouping_variable"]])) %>%
    #mutate(
    #  `PK-to-Detected` = case_when(
    #    num_PK_entries == 0 & num_Input_entries == 1 ~ "none-from-one",
    #    num_PK_entries == 0 & num_Input_entries > 1 ~ "none-from-many",
    #    num_PK_entries  == 1 & num_Input_entries == 1 ~ "one-from-one",
    #    num_PK_entries  == 1 & num_Input_entries > 1 ~ "one-from-many",
    #    num_PK_entries  > 1 & num_Input_entries == 1 ~ "many-from-one",
    #    num_PK_entries  > 1 & num_Input_entries > 1 ~ "many-from-many",
    #    num_PK_entries > 1 & num_Input_entries == 0 ~ "many-from-none",
    #    num_PK_entries == 1 & num_Input_entries == 0 ~ "one-from-none",
    #    TRUE ~ NA_character_
    #  )
    #) %>%
    #ungroup() %>%
    # Assign ActionRequired (If many-to-many in either case)
    #mutate(
    #  ActionRequired = case_when(
    #    `Detected-to-PK` == "many-to-many" | `PK-to-Detected` == "many-from-many" ~ "Check",
    #    TRUE ~ "None"
     # )
    #)

  # 4. Create summary table
  Values_data <- unique(data[[metadata_info[["InputID"]]]])
  Values_PK <- unique(PK_long[[metadata_info[["PriorID"]]]])

  summary_df <- tibble::tibble(
    !!sym(metadata_info[["InputID"]]) := Values_data,
    found_match_in_PK = NA,
    matches = NA_character_,
    match_overlap_percentage = NA_real_,
    original_count = NA_integer_,
    matches_count = NA_integer_
  )

  # Populate the summary data frame
  for(i in seq_along(Values_data)) {
    # Handle NA case explicitly
    if (is.na(Values_data[i])) {
      summary_df$original_count[i] <- 0
      summary_df$matches_count[i] <- 0
      summary_df$match_overlap_percentage[i] <- NA
      summary_df$found_match_in_PK[i] <- NA # could also set it to FALSE but making NA for now for plotting
      summary_df$matches[i] <- NA
    } else {
      # Split each cell into individual entries and trim whitespace
      entries <- trimws(unlist(strsplit(as.character(Values_data[i]), ",\\s*"))) # delimiter = "," or ", "

      # Identify which entries are in the lookup set
      matched <- entries[entries %in% Values_PK]

      # Determine if any match was found
      summary_df$found_match_in_PK[i] <- length(matched) > 0

      # Concatenate matched entries into a single string
      summary_df$matches[i] <- paste(matched, collapse = ", ")

      # Calculate and store counts
      summary_df$original_count[i] <- length(entries)
      summary_df$matches_count[i] <- length(matched)

      # Calculate fraction: matched entries / total entries
      if(length(entries) > 0) {
        summary_df$match_overlap_percentage[i] <- (length(matched) / length(entries))*100
      } else {
        summary_df$match_overlap_percentage[i] <- NA
      }
    }
  }

  summary_df <- merge(x= summary_df,
                      y= merged_df%>%
                        dplyr::select(-c(OriginalGroup_PK, OriginalGroup_data))%>%
                        distinct(!!sym(metadata_info[["PriorID"]]), .keep_all = TRUE),
                      by.x= metadata_info[["InputID"]] ,
                      by.y= metadata_info[["PriorID"]],
                      all.x=TRUE)

  if(PK_MultipleIDs){
    summary_df <- merge(x= summary_df, y= data, by=metadata_info[["InputID"]], all.x=TRUE)%>%
      distinct(!!sym(metadata_info[["InputID"]]), OriginalEntry_PK, .keep_all = TRUE)
  }else{
    summary_df <- merge(x= summary_df, y= data, by=metadata_info[["InputID"]], all.x=TRUE)%>%
      distinct(!!sym(metadata_info[["InputID"]]), .keep_all = TRUE)
  }

  # 5. Merge back on input data to retain Nulls and duplications in case the user wants this (e.g. for plotting or inspecting further)

  # Function: add_NA_to_table
  #
  # Description:
  #   This function takes two data frames:
  #     - table_with_NA: the original table that may contain rows where the key column is NA.
  #     - table_without_NA: a processed table (e.g. from a join) that contains extra columns and excludes rows where the key is NA.
  #
  #   The function extracts rows from table_with_NA where the key is NA, extends these rows by adding any extra columns
  #   (present in table_without_NA but not in table_with_NA) with NA values, and then binds these extended rows to table_without_NA.
  #
  # Parameters:
  #   table_with_NA: data frame containing the original data (e.g. FeatureMetadata_Biocrates).
  #   table_without_NA: data frame containing the processed data with extra columns (e.g. tempnew).
  #   key: The column name (as a string) used as the key for matching (e.g. "HMDB").
  #
  # Returns:
  #   A combined data frame that includes the rows from table_without_NA along with the extended NA rows from table_with_NA.
  add_NA_to_table <- function(table_with_NA, table_without_NA, key) {

    # Subset rows from the original table where the key column is NA
    na_rows <- dplyr::filter(table_with_NA, is.na(.data[[key]]))

    # Identify extra columns present in table_without_NA that are not in the original table
    extra_cols <- setdiff(names(table_without_NA), names(table_with_NA))

    # Extend the NA rows by adding the extra columns, setting their values to NA
    na_rows_extended <- dplyr::mutate(na_rows, !!!setNames(rep(list(NA), length(extra_cols)), extra_cols))

    # Reorder columns to match the structure of table_without_NA
    na_rows_extended <- dplyr::select(na_rows_extended, dplyr::all_of(names(table_without_NA)))

    # Combine the processed table with the extended NA rows
    combined_table <- dplyr::bind_rows(table_without_NA, na_rows_extended)

    return(combined_table)
  }


  # Function: create_duplicates_table
  #
  # Description:
  #   This function takes two data frames:
  #     - table_with_duplicates: the original table that may contain duplicate rows based on the key.
  #     - table_without_duplicates: a deduplicated table (e.g. from a join) that includes extra columns.
  #
  #   The function identifies duplicate rows (non-NA keys that appear more than once, excluding the first occurrence)
  #   in table_with_duplicates, then left joins these duplicate rows with the extra columns from table_without_duplicates.
  #   This ensures that all duplicate rows receive the same extra column values as the first occurrence.
  #
  # Parameters:
  #   table_with_duplicates: data frame containing the original data that may include duplicate keys.
  #   table_without_duplicates: data frame with the processed data (first occurrence for each key and extra columns).
  #   key: The column name (as a string) used as the key for matching (e.g. "HMDB").
  #
  # Returns:
  #   A data frame containing the duplicate rows, extended with the extra columns from table_without_duplicates.
  create_duplicates_table <- function(table_with_duplicates, table_without_duplicates, key) {

    # Identify extra columns present in table_without_duplicates that are not in the original table
    extra_cols <- setdiff(names(table_without_duplicates), names(table_with_duplicates))

    # Extract duplicate rows: for each non-NA key, keep rows beyond the first occurrence
    dup_rows <- table_with_duplicates %>%
      dplyr::filter(!is.na(.data[[key]])) %>%
      dplyr::group_by(.data[[key]]) %>%
      dplyr::filter(dplyr::row_number() > 1) %>%
      dplyr::ungroup()

    # For each duplicate row, join the extra columns from table_without_duplicates so that all duplicates
    # receive the same extra values as the first occurrence
    dup_rows_extended <- dplyr::left_join(
      dup_rows,
      dplyr::select(table_without_duplicates, dplyr::all_of(c(key, extra_cols))),
      by = key
    ) %>%
      dplyr::select(dplyr::all_of(names(table_without_duplicates)))

    return(dup_rows_extended)
  }

  # Create the table with NA rows added
  temp_results_NAs_added <- add_NA_to_table(data_Original, summary_df, metadata_info[["InputID"]])
  # Create the table with duplicate key (metadata_info[["InputID"]]) rows extended
  temp_results_of_duplicates <- create_duplicates_table(data_Original, summary_df, metadata_info[["InputID"]])
  # Combine these to get a summary table that includes both NA and duplicate rows
  summary_df_with_NA_and_duplicates <- dplyr::bind_rows(temp_results_NAs_added, temp_results_of_duplicates)

  # Now for the user let's also create separate dfs with just the NA values and just the duplicates, in case they want to inspect this easier
  summary_df_only_NA <- summary_df_with_NA_and_duplicates %>%
    dplyr::filter(is.na(.data[[metadata_info[["InputID"]]]]))

  summary_df_only_duplicates <- summary_df_with_NA_and_duplicates %>%
    dplyr::filter(!is.na(.data[[metadata_info[["InputID"]]]])) %>%  # Exclude NA values
    dplyr::group_by(.data[[metadata_info[["InputID"]]]]) %>%
    dplyr::filter(dplyr::n() > 1) %>%             # Keep groups with duplicates
    dplyr::ungroup()

  # 6. Messages and summarise
  message <- paste0("data has multiple IDs per measurement = ", data_MultipleIDs, ". prior_knowledge has multiple IDs per entry = ", PK_MultipleIDs, ".", sep="")
  message1 <- paste0("data has ", dplyr::n_distinct(unique(data[[metadata_info[["InputID"]]]])), " unique entries with " ,dplyr::n_distinct(unique(data_long[[metadata_info[["InputID"]]]])) ," unique ", metadata_info[["InputID"]], " IDs. Of those IDs, ", nrow(summary_df%>% dplyr::filter(matches_count == 1)), " match, which is ", (nrow(summary_df%>% dplyr::filter(matches_count == 1)) / dplyr::n_distinct(unique(data_long[[metadata_info[["InputID"]]]])))*100, "%." , sep="")
  message2 <- paste0("prior_knowledge has ", dplyr::n_distinct(prior_knowledge[[metadata_info[["PriorID"]]]]), " unique entries with " ,dplyr::n_distinct(PK_long[[metadata_info[["PriorID"]]]]) ," unique ", metadata_info[["PriorID"]], " IDs. Of those IDs, ", nrow(summary_df%>% dplyr::filter(matches_count == 1)), " are detected in the data, which is ", (nrow(summary_df%>% dplyr::filter(matches_count == 1)) / dplyr::n_distinct(PK_long[[metadata_info[["PriorID"]]]]))*100, "%.")

  if(nrow(summary_df%>% dplyr::filter(ActionRequired == "Check"))>=1){
    #warning <- paste0("There are cases where multiple detected IDs match to multiple prior knowledge IDs of the same category") # "Check"

  }

  ## ------------------ Plot summary ----------------------##
  # x = "Class" and y = Frequency. Match Status can be colour of if no class provided class = Match status.
  # Check Biocrates code.


  ## ------------------ Save Results ----------------------##
  ResList <- list("data_Matched" = summary_df,
                  "data_Matched_NA_and_duplicates" = summary_df_with_NA_and_duplicates,
                  "data_Matched_only_NA" = summary_df_only_NA,
                  "data_Matched_only_duplicates" = summary_df_only_duplicates)

  suppressMessages(suppressWarnings(
  save_res(inputlist_df=ResList,
                       inputlist_plot= NULL,
                       save_table=save_table,
                       save_plot=NULL,
                       path= Subfolder,
                       file_name= "CheckMatchID_Detected-to-PK",
                       core=FALSE,
                       print_plot=FALSE)))

   #Return
   invisible(return(ResList))
}


##########################################################################################
### ### ### cluster Prior Knowledge ### ### ###
##########################################################################################

#' Deal with pathway overlap in prior knowledge
#'
#' @param data dataframe with at least one column with the target (e.g. metabolite) and a column source (e.g. term).
#' @param metadata_info = c(InputID="MetaboliteID", grouping_variable="term"),
#'
#' @examples
#' KEGG_Pathways <- MetaProViz::metsigdb_kegg()
#' data = KEGG_Pathways
#'
#'
#' @importFrom igraph graph_from_adjacency_matrix components
#' @noRd
cluster_pk <- function(data, # This can be either the original PK (e.g. KEGG pathways), but it can also be the output of enrichment results (--> meaning here we would cluster based on detection!)
                      metadata_info= c(InputID="MetaboliteID", grouping_variable="term"),
                      clust = "Graph", # Options: "Graph", "Hierarchical",
                      matrix ="percentage", # Choose "pearson", "spearman", "kendall", or "percentage"
                      min= 2 # minimum pathways per cluster

){

  # cluster PK before running enrichment analysis --> add another column that groups the data based on the pathway overlap:
  # provide different options for clustering (e.g. % of overlap, semantics similarity) --> Ramp uses % of overlap, semnatics similarity: https://yulab-smu.top/biomedical-knowledge-mining-book/GOSemSim.html


  ## ------------------ Check Input ------------------- ##


  ## ------------------ Create output folders and path ------------------- ##



  ######################################################################################################################################
  ## ------------------ cluster the data ------------------- ##
  # 1. Create a list of unique MetaboliteIDs for each term
  term_metabolites <- data %>%
    dplyr::group_by(!!sym(metadata_info[["grouping_variable"]])) %>%
    dplyr::summarize(MetaboliteIDs = list(unique(!!sym(metadata_info[["InputID"]])))) %>%
    dplyr::ungroup()

  #2. Create the overlap matrix based on different methods:
  if (matrix == "percentage") {# Compute pairwise overlaps
    term_overlap <- combn(term_metabolites[[metadata_info[["grouping_variable"]]]], 2, function(terms) {
      term1_ids <- term_metabolites$MetaboliteIDs[term_metabolites[[metadata_info[["grouping_variable"]]]] == terms[1]][[1]]
      term2_ids <- term_metabolites$MetaboliteIDs[term_metabolites[[metadata_info[["grouping_variable"]]]] == terms[2]][[1]]

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
    terms <- term_metabolites[[metadata_info[["grouping_variable"]]]]
    metabolites <- unique(unlist(term_metabolites$MetaboliteIDs)) #[[metadata_info[["InputID"]]]]

    binary_matrix <- matrix(0, nrow = length(terms), ncol = length(metabolites), dimnames = list(terms, metabolites))
    for (i in seq_along(terms)) {
      metabolites_for_term <- term_metabolites$MetaboliteIDs[[i]] #[[metadata_info[["InputID"]]]]
      binary_matrix[i, colnames(binary_matrix) %in% metabolites_for_term] <- 1
    }

    # Compute correlation matrix: square matrix used to represent the pairwise correlation coefficients between variables or terms
    # correlation matrix 𝐶 C is an 𝑛 × 𝑛 n×n matrix where each element 𝐶 𝑖 𝑗 C ij  is the correlation coefficient between the variables 𝑋𝑖 Xi  and 𝑋 𝑗 X j
    #The correlation matrix measures the strength and direction of linear relationships between variables.
    correlation_matrix <- cor(t(binary_matrix), method = matrix)

    # Convert to distance matrix
    overlap_matrix <- 1 - correlation_matrix
  }
  # 3. cluster terms based on overlap threshold
  threshold <- 0.7 # Define similarity threshold
  term_clusters <- term_overlap %>%
    dplyr::filter(Overlap >= threshold) %>%
    dplyr::select(Term1, Term2)

  # 4. clustering
  if (clust == "Graph") { #Use Graph-based clustering
  # Here we need the distance matrix:
  overlap_matrix <- 1 - correlation_matrix

  # An adjacency matrix represents a graph structure and encodes the relationships between nodes (vertices)
  # Add weight (can also represent unweighted graphs)
  adjacency_matrix <- exp(-overlap_matrix^2)  # Applying Gaussian kernel to convert distance into similarity

  # Create a graph from the adjacency matrix
  g <- igraph::graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected", weighted = TRUE)
  initial_clusters <- igraph::components(g)$membership
  term_metabolites$cluster <- initial_clusters[match(term_metabolites[[metadata_info[["grouping_variable"]]]], names(initial_clusters))]
  } else if (clust == "Hierarchical") { # Hierarchical clustering
    hclust_result <- hclust(as.dist(distance_matrix), method = "average") # make methods into parameters!
    num_clusters <- 4
    term_clusters_hclust <- cutree(hclust_result, k = num_clusters)

    term_metabolites$cluster <- paste0("cluster", term_clusters_hclust[match(terms, names(term_clusters_hclust))])
    #term_metabolites$cluster <- clusters[match(term_metabolites[[metadata_info[["grouping_variable"]]]], names(clusters))]
  } else {
    stop("Invalid clustering method specified in clust parameter.")
  }

  # 5. Merge cluster group information back to the original data
  df <- data %>%
    dplyr::left_join(term_metabolites %>% select(!!sym(metadata_info[["grouping_variable"]]), cluster), by = metadata_info[["grouping_variable"]])%>%
    dplyr::mutate(cluster = ifelse(
      is.na(cluster),
        "None", # Assign "None" to NAs
        paste0("cluster", cluster) # Convert numeric IDs to descriptive labels
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
#' @param mat data matrix used as input for enrichment analysis
#' @param net Prior Knowledge used as input for enrichment analysis
#' @param res Results returned from the enrichment analysis
#' @param .source used as input for enrichment analysis
#' @param .target used as input for enrichment analysis
#' @param complete TRUE or FALSE, weather only .source with results should be returned or all .source in net.
#'
#' @noRd

# Better function Name and parameter names needed
# Use in ORA functions and showcase in vignette with decoupleR output

add_info <- function(mat,
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

  #add percentage of percentage_of_Genes_detected
  res_Add$targets_detected_percentage <-round(((res_Add$targets_detected_num/res_Add$targets_num)*100),digits=2)

  #sort by score
  res_Add<-res_Add%>%
    arrange(desc(as.numeric(as.character(score))))

  ## ------------------ Save and return ------------------- ##
  Output<-res_Add
}


##########################################################################################
### ### ### Compare Prior Knowledge resources against each other or themselves ### ### ###
##########################################################################################

#' Compare Prior Knowledge Resources and/or Columns within a Single Resource and Generate an UpSet Plot
#'
#' This function compares gene and/or metabolite features across multiple prior knowledge (PK) resources or,
#' if a single resource is provided with a vector of column names in \code{metadata_info}, compares columns within that resource.
#'
#' In the multi-resource mode, each element in \code{data} represents a PK resource (either as a data frame or a recognized resource name)
#' from which a set of features is extracted. A binary summary table is then constructed and used to create an UpSet plot.
#'
#' In the within-resource mode, a single data frame is provided (with \code{data} containing one element) and its \code{metadata_info} entry
#' is a vector of column names to compare (e.g., binary indicators for different annotations). In this case, the function expects the data frame
#' to have a grouping column named \code{"Class"} (or, alternatively, a column specified via the \code{class_col} attribute in \code{metadata_info})
#' that is used for grouping in the UpSet plot.
#'
#' @param data A named list where each element corresponds to a prior knowledge (PK) resource. Each element can be:
#'        \itemize{
#'          \item A data frame containing gene/metabolite identifiers (and additional columns for within-resource comparison),
#'          \item A character string indicating the resource name. Recognized names include (but are not limited to): \code{"Hallmarks"},
#'                \code{"Gaude"}, \code{"MetalinksDB"}, and \code{"RAMP"} (or \code{"metsigdb_chemicalclass"}). In the latter case, the function
#'                will attempt to load the corresponding data automatically.
#'        }
#'
#' @param metadata_info A named list (with names matching those in \code{data}) where each element is either a character string or a
#'        character vector indicating the column name(s) to extract features. For multiple-resource comparisons, these refer to the columns
#'        containing feature identifiers. For within-resource comparisons, the vector should list the columns to compare (e.g., \code{c("CHEBI", "HMDB", "LIMID")}).
#'        In within-resource mode, the input data frame is expected to contain a column named \code{"Class"} (or a grouping column specified via the
#'        \code{class_col} attribute). \emph{If no grouping column is found, a default grouping column named \code{"Group"} (with all rows assigned the same value) is created.}
#'
#' @param filter_by Character. Optional filter for the resulting features when comparing multiple resources.
#'        Options are: \code{"both"} (default), \code{"gene"}, or \code{"metabolite"}. This parameter is ignored in within-resource mode.
#'
#' @param plot_title Character. Title for the UpSet plot. Default is \code{"Overlap of Prior Knowledge Resources"}.
#' @param palette_type Character. Color palette to be used in the plot. Default is \code{"polychrome"}.
#' @param output_file Character. Optional file path to save the generated plot; if \code{NULL}, the plot is not saved.
#'
#' @return A list containing two elements: \itemize{
#'     \item{summary_table: A data frame representing either:
#'         \itemize{
#'             \item{the binary summary matrix of feature presence/absence across
#'                 multiple resources, or}
#'             \item{the original data frame (augmented with binary columns and a
#'                 \code{None} column) in within-resource mode.}
#'         }
#'     }
#'     \item{upset_plot: The UpSet plot object generated by the function.}
#' }
#'
#' @examples
#' ## Example 1: Multi-Resource Comparison
#'
#' # Using automatic data loading for multiple resources.
#' data <- list(Hallmarks = "Hallmarks", Gaude = "Gaude",
#'                 MetalinksDB = "MetalinksDB", RAMP = "metsigdb_chemicalclass")
#' res <- compare_pk(data = data)
#'
#' # Filtering to include only gene features:
#' res_genes <- compare_pk(data = data, filter_by = "gene")
#'
#' ## Example 2: Within-Resource Comparison (Comparing Columns Within a Single data Frame)
#'
#' # Assume FeatureMetadata_Biocrates is a data frame with columns: "TrivialName", "CHEBI", "HMDB", "LIMID", and "Class".
#' # Here the "Class" column is used as the grouping variable in the UpSet plot.
#' data_single <- list(Biocft = FeatureMetadata_Biocrates)
#' metadata_info_single <- list(Biocft = c("CHEBI", "HMDB", "LIMID"))
#'
#' res_single <- compare_pk(data = data_single, metadata_info = metadata_info_single,
#'                           plot_title = "Overlap of BioCrates Columns")
#'
#' ## Example 3: Custom data Frames with Custom Column Names
#'
#' # Example with preloaded data frames and custom column names:
#' hallmarks_df <- data.frame(feature = c("HMDB0001", "GENE1", "GENE2"), stringsAsFactors = FALSE)
#' gaude_df <- data.frame(feature = c("GENE2", "GENE3"), stringsAsFactors = FALSE)
#' metalinks_df <- data.frame(hmdb = c("HMDB0001", "HMDB0002"),
#'                            gene_symbol = c("GENE1", "GENE4"), stringsAsFactors = FALSE)
#' ramp_df <- data.frame(class_source_id = c("HMDB0001", "HMDB0003"), stringsAsFactors = FALSE)
#' data <- list(Hallmarks = hallmarks_df, Gaude = gaude_df,
#'                 MetalinksDB = metalinks_df, RAMP = ramp_df)
#' metadata_info <- list(Hallmarks = "feature", Gaude = "feature",
#'                      MetalinksDB = c("hmdb", "gene_symbol"), RAMP = "class_source_id")
#' res <- compare_pk(data = data, metadata_info = metadata_info, filter_by = "metabolite")
#'
#' @importFrom dplyr mutate select
#' @importFrom utils write.csv
#'
#' @export
compare_pk <- function(data,
                       metadata_info = NULL,
                       filter_by = c("both", "gene", "metabolite"),
                       plot_title = "Overlap of Prior Knowledge Resources",
                       name_col = "TrivialName",
                       palette_type = "polychrome",
                       output_file = NULL) {
  ###########################################################################
  ## ------------ Create log file ----------- ##
  metaproviz_init()

  ## ------------ Check Input files ----------- ##
  # Match filter argument
  filter_by <- match.arg(filter_by)

  # Validate data input
  if (!is.list(data) || length(data) < 1) {
    message <- paste0("data must be a non-empty list.")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }
  if (is.null(names(data)) || any(names(data) == "")) {
    message <- paste0("data must be a named list with resource names.")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }

  ###########################################################################
  ## ----------- Input ----------- ##
  # Define resource lookup table with information on how to retrieve and transform each resource.
  resource_definitions <- list(
    hallmarks = list(
      var = "Hallmark_Pathways",
      load_fun = hallmarks,
      transform_fun = function(x) {
        resource_object <- MetaProViz::make_gene_metab_set(input_pk = x,
                                                         metadata_info = c(Target = "gene"),
                                                         pk_name = "Hallmarks")
        if ("GeneMetabSet" %in% names(resource_object)) {
          resource_object$GeneMetabSet
        } else {
          stop("make_gene_metab_set for Hallmarks did not return 'GeneMetabSet'.")
        }
      },
      default_col = "feature"
    ),
    gaude = list(
      var = "Gaude_Pathways",
      load_fun = gaude_pathways,
      transform_fun = function(x) {
        resource_object <- MetaProViz::make_gene_metab_set(input_pk = x,
                                                         metadata_info = c(Target = "gene"),
                                                         pk_name = "Gaude")
        if ("GeneMetabSet" %in% names(resource_object)) {
          resource_object$GeneMetabSet
        } else {
          stop("make_gene_metab_set for Gaude did not return 'GeneMetabSet'.")
        }
      },
      default_col = "feature"
    ),
    metalinksdb = list(
      var = "MetalinksDB",
      load_fun = MetaProViz::metsigdb_metalinks,
      transform_fun = function(x) {
        if ("MetalinksDB" %in% names(x)) {
          x$MetalinksDB
        } else {
          stop("Loaded MetalinksDB does not contain a 'MetalinksDB' element.")
        }
      },
      default_col = c("hmdb", "gene_symbol")
    ),
    ramp = list(
      var = "ChemicalClass_MetabSet",
      load_fun = MetaProViz::metsigdb_chemicalclass,
      transform_fun = function(x) {
        x  # for RAMP, assume the global variable itself is the data frame.
      },
      default_col = "class_source_id"
    )
  )

  # Preprocess data: auto‑load resources if they are provided as strings.
  for (res in names(data)) {
    if (!inherits(data[[res]], "data.frame") && is.character(data[[res]])) {
      resource_id <- tolower(data[[res]])
      if (resource_id %in% names(resource_definitions)) {
        res_def <- resource_definitions[[resource_id]]
        if (exists(res_def$var, envir = .GlobalEnv)) {
          resource_object <- get(res_def$var, envir = .GlobalEnv)
        } else {
          resource_object <- res_def$load_fun()
        }
        data[[res]] <- res_def$transform_fun(resource_object)
        if (is.null(metadata_info[[res]])) {
          metadata_info[[res]] <- res_def$default_col
        }
      }
    }
  }

  # Initialize metadata_info if not provided.
  if (is.null(metadata_info)) {
    metadata_info <- list()
  }

  # Determine if we are in within-resource mode.
  # If only one resource is provided and its metadata_info entry has >1 column, assume within-resource comparison.
  single_resource <- (length(data) == 1)
  within_resource_mode <- FALSE
  if (single_resource) {
    resource_name <- names(data)[1]
    if (!is.null(metadata_info[[resource_name]]) &&
        length(metadata_info[[resource_name]]) > 1) {
      within_resource_mode <- TRUE
    }
  }

  if (within_resource_mode) {
    # ===== Within-Resource Comparison Mode =====
    # Retrieve the single data frame.
    resource_data <- data[[resource_name]]

    # Identify the intersection columns based on metadata_info.
    intersect_cols <- metadata_info[[resource_name]]
    missing_cols <- setdiff(intersect_cols, colnames(resource_data))
    if (length(missing_cols) > 0) {
      stop("The following intersection column(s) specified in metadata_info were not found in resource '",
           resource_name, "': ", paste(missing_cols, collapse = ", "))
    }

    # Identify a column for grouping. If none exists, create a default grouping column.
    if ("Class" %in% colnames(resource_data)) {
      class_col <- "Class"
    } else if (!is.null(attr(metadata_info[[resource_name]], "class_col"))) {
      class_col <- attr(metadata_info[[resource_name]], "class_col")
    } else {
      # No grouping column provided—create a default column named "Group" with the same value for all rows.
      resource_data$Group <- "All"
      class_col <- "Group"
    }

    # Convert the specified intersection columns to binary (0/1). Here non-NA and values != 0 are treated as present.
    binary_suffix <- "_bin"
    for (col in intersect_cols) {
      new_col <- paste0(col, binary_suffix)
      resource_data[[new_col]] <- as.integer(!is.na(resource_data[[col]]) & (resource_data[[col]] != 0) & (resource_data[[col]] != ''))
    }

    # Identify the binary columns based on the suffix
    bin_cols <- grep(paste0(binary_suffix, "$"), colnames(resource_data), value = TRUE)

    # Create the "None" column using the binary columns
    resource_data$None <- as.integer(rowSums(resource_data[, bin_cols, drop = FALSE]) == 0)

    # Create df_summary, potentially with the name column (default is TrivialName, if not found, it will not be included - only used to help interpret summary table)
    if (name_col %in% colnames(resource_data)) {
      summary_cols <- c(name_col, bin_cols, "None", class_col)
    }
    else {
      summary_cols <- c(bin_cols, "None", class_col)
    }
    df_summary <- resource_data[, summary_cols, drop = FALSE]
    # Rename the cols again
    # Find indices of columns ending in "_bin"
    bin_cols_idx <- grep("_bin$", names(df_summary))
    # Remove the suffix from these column names
    names(df_summary)[bin_cols_idx] <- sub("_bin$", "", names(df_summary)[bin_cols_idx])

    # Generate the UpSet plot.
    upset_plot <- MetaProViz:::viz_upset(
      df = df_summary,
      class_col = class_col,
      intersect_cols = c(intersect_cols, "None"),
      plot_title = plot_title,
      palette_type = palette_type,
      output_file = output_file
    )

    return(list(summary_table = df_summary, upset_plot = upset_plot))

  } else {
    # ===== Multi-Resource Comparison Mode =====
    # Process each resource in data.
    for (res in names(data)) {
      resource_val <- data[[res]]
      if (!inherits(resource_val, "data.frame")) {
        if (!is.character(resource_val)) {
          stop("Each element in data must be either a data frame or a character string indicating a resource name.")
        }
        resource_id <- tolower(resource_val)
        if (resource_id %in% c("metsigdb_chemicalclass")) {
          resource_id <- "ramp"
        }
        if (!resource_id %in% names(resource_definitions)) {
          stop("Unknown resource identifier: ", resource_val,
               ". Please provide a data frame or a valid resource name.")
        }
        res_def <- resource_definitions[[resource_id]]
        if (exists(res_def$var, envir = .GlobalEnv)) {
          resource_object <- get(res_def$var, envir = .GlobalEnv)
        } else {
          resource_object <- res_def$load_fun()
        }
        data[[res]] <- res_def$transform_fun(resource_object)
        if (is.null(metadata_info[[res]])) {
          metadata_info[[res]] <- res_def$default_col
        }
      } else {
        resource_id <- tolower(res)
        if (is.null(metadata_info[[res]]) && resource_id %in% names(resource_definitions)) {
          metadata_info[[res]] <- resource_definitions[[resource_id]]$default_col
        } else if (is.null(metadata_info[[res]])) {
          stop("metadata_info must be provided for resource: ", res)
        }
      }
    }

    # Extract features from each resource based on metadata_info.
    resource_features <- list()
    for (res in names(data)) {
      resource_data <- data[[res]]
      cols <- metadata_info[[res]]
      if (!all(cols %in% colnames(resource_data))) {
        stop(paste("Column(s)", paste(cols, collapse = ", "),
                   "not found in resource", res))
      }
      features <- if (length(cols) > 1) {
        unique(unlist(lapply(cols, function(col) na.omit(resource_data[[col]]))))
      } else {
        unique(na.omit(resource_data[[cols]]))
      }
      resource_features[[res]] <- as.character(features)
    }

    # Compile all unique features across resources.
    all_features <- unique(unlist(resource_features))

    # Create the binary summary table.
    df_binary <- data.frame(Feature = all_features, stringsAsFactors = FALSE)
    for (res in names(resource_features)) {
      df_binary[[res]] <- as.integer(all_features %in% resource_features[[res]])
    }
    df_binary$Type <- ifelse(grepl("^HMDB", df_binary$Feature), "metabolite (HMDB)", "gene")
    resource_cols <- names(resource_features)
    df_binary$None <- as.integer(rowSums(df_binary[, resource_cols, drop = FALSE]) == 0)

    # Optionally filter the summary table.
    if (filter_by == "gene") {
      df_binary <- subset(df_binary, Type == "gene")
    } else if (filter_by == "metabolite") {
      df_binary <- subset(df_binary, Type == "metabolite (HMDB)")
    }

    # Generate the UpSet plot.
    upset_plot <- MetaProViz:::viz_upset(
      df = df_binary,
      class_col = "Type",
      intersect_cols = resource_cols,
      plot_title = plot_title,
      palette_type = palette_type,
      output_file = output_file
    )

    return(list(summary_table = df_binary, upset_plot = upset_plot))
  }
}


##########################################################################################
### ### ### Helper function to count number of entries for an ID column value and plot ### ### ###
##########################################################################################

#' Count Entries and Generate a Histogram Plot for a Specified Column
#'
#' This function processes a data frame column by counting the number of entries within each cell.
#' It considers both \code{NA} values and empty strings as zero entries, and categorizes each cell as
#' "No ID", "Single ID", or "Multiple IDs" based on the count. A histogram is then generated to visualize
#' the distribution of entry counts.
#'
#' @param data A data frame containing the data to be analyzed.
#' @param column A string specifying the name of the column in \code{data} to analyze.
#' @param delimiter A string specifying the delimiter used to split cell values. Defaults to \code{","}.
#' @param fill_colors A named character vector providing colors for each category. Defaults to
#'   \code{c("No ID" = "#FB8072", "Single ID" = "#B3DE69", "Multiple IDs" = "#80B1D3")}.
#' @param binwidth Numeric value specifying the bin width for the histogram. Defaults to \code{1}.
#' @param title_prefix A string to use as the title of the plot. If \code{NULL} (default), the title
#'   will be generated as "Number of <column> IDs per Biocrates Cell".
#' @param save_plot \emph{Optional: } Select the file type of output plots. Options are svg, png, pdf. \strong{Default = svg}
#' @param save_table \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt". \strong{Default = "csv"}
#' @param print_plot \emph{Optional: } TRUE or FALSE, if TRUE Volcano plot is saved as an overview of the results. \strong{Default = TRUE}
#' @param path \emph{Optional:} Path to the folder the results should be saved at. \strong{Default = NULL}
#'
#' @return A list with two elements:
#'   \item{result}{A data frame that includes three additional columns: \code{was_na} (logical indicator
#'                  of missing or empty cells), \code{entry_count} (number of entries in each cell), and
#'                  \code{id_label} (a categorical label based on the entry count).}
#'   \item{plot}{A \code{ggplot} object representing the histogram of entry counts.}
#'
#' @export
count_ids <- function(data,
                      column,
                      delimiter = ",",
                      fill_colors = c("No ID" = "#FB8072",
                                                   "Single ID" = "#B3DE69",
                                                   "Multiple IDs" = "#80B1D3"),
                      binwidth = 1,
                      title_prefix = NULL,
                      save_plot = "svg",
                      save_table = "csv",
                      print_plot = TRUE,
                      path = NULL) {
  #@Macabe:
  #move named fill colors inside of the function. if the user provides other string of colors, we change them.
  #we need to specify that NA is counted as none
  #we need to check for duplications (i.e. is the trivialname duplicated in the data frame, remove this and give a warning)
  #give the user the change to pass multiple columns to analyse, which would mean we create a plot for each column and label the plot with the column name
  #add count_ids function into Equivalent IDs: make a plot before and after equivalent IDs and put them side by side to return as the QC plot of the function
  #add save_plot and save_table
  #return both standard plot and Plot_Sized (=nice version of the plot which is saved. the other plot is still returned, since this is the ggplot version and can be changed further)
  #create subtitle and not title prefix


  ## ------------------  logger initiation ------------------- ##


  ## ------------------  Checks ------------------- ##


  ## ------------------  Create output folders and path ------------------- ##
  if(is.null(save_table)==FALSE ){
    folder <- save_path(folder_name= "prior_knowledge",
                       path=path)

    Subfolder <- file.path(folder, "Count_MetaboliteIDs")
    if (!dir.exists(Subfolder)) {dir.create(Subfolder)}
  }

  ## ------------------  data table ------------------- ##
  # Process the data: count entries and label each cell based on the number of entries.
  processed_data <- dplyr::mutate(
    data,
    was_na = is.na(.data[[column]]) | .data[[column]] == "",
    entry_count = sapply(.data[[column]], function(cell) {
      if (is.na(cell) || cell == "") {
        0  # Treat NA or empty as 0 entries for counting
      } else {
        length(unlist(strsplit(as.character(cell), delimiter)))
      }
    }),
    id_label = dplyr::case_when(
      entry_count == 0 ~ "No ID",
      entry_count == 1 ~ "Single ID",
      entry_count >= 2 ~ "Multiple IDs"
    )
  )


  ## ------------------  plot ------------------- ##
  # Generate the plot title if not provided
  if (is.null(title_prefix)) {
    plot_title <- paste("Number of", column, "IDs per measured peak.")
  } else {
    plot_title <- title_prefix
  }

  # Create plot using ggplot2
  if(length(unique(processed_data$entry_count))>1){# Create the histogram plot using ggplot2
    plot_obj <- ggplot2::ggplot(processed_data, ggplot2::aes(x = entry_count, fill = id_label)) +
      ggplot2::geom_histogram(binwidth = binwidth, boundary = -0.5, color = "black")+# Set axis to start at 0
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.05)))+
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05)))
  }else{# Create bargraph plot using ggplot2
    plot_obj <- ggplot2::ggplot(processed_data, ggplot2::aes(x = as.factor(entry_count), fill = id_label)) +
      ggplot2::geom_bar(color = "black")+  # Set axis to start at 0
     ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05)))

  }

  plot_obj <- plot_obj +
    ggplot2::scale_fill_manual(values = fill_colors) +
    ggplot2::labs(title = plot_title,
                  x = "Number of IDs",
                  y = "Frequency",
                  fill = paste0(column, " IDs", sep="")) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 12),
      legend.position.inside = c(0.8, 0.8),
      legend.justification = c("right", "top"),
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 12)
    )

  # Make the nice plot:
  Plot_Sized <-  plot_grob_superplot(input_plot=plot_obj, metadata_info= c(Conditions="id_label", Superplot = TRUE), metadata_sample= processed_data%>%dplyr::rename("Conditions"="entry_count") , plot_name = plot_title, subtitle = "", plot_type="Bar")
  plot_height <- grid::convertUnit(Plot_Sized$height, 'cm', valueOnly = TRUE)
  plot_width <- grid::convertUnit(Plot_Sized$width, 'cm', valueOnly = TRUE)
  Plot_Sized %<>%
    {ggplot2::ggplot() + annotation_custom(.)} %>%
    add(theme(panel.background = ggplot2::element_rect(fill = "transparent")))


  ## ------------------  save and return ------------------- ##
  suppressMessages(suppressWarnings(
    save_res(inputlist_df=list("Table"=processed_data),#This needs to be a list, also for single comparisons
            inputlist_plot=list("Plot_Sized"=Plot_Sized) ,
            save_table= save_table,
            save_plot=save_plot,
            path= Subfolder,
            file_name= "Count_MetaboliteIDs",
            core=FALSE,
            print_plot=print_plot)))

  OutputList <- list()
  OutputList <- list("Table"=processed_data, "Plot"=plot_obj, "Plot_Sized"=Plot_Sized)

  # Return the processed data and the plot object as a list
  return(invisible(OutputList))
}

