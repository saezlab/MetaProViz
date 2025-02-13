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
#' @param InputData Dataframe with at least one column with the target (e.g. metabolite), you can add other columns such as source (e.g. term). Must be "long" DF, meaning one ID per row.
#' @param SettingsInfo \emph{Optional: } Column name of Target in Input_GeneSet. \strong{Default = list(InputID="MetaboliteID" , GroupingVariable="term")}
#' @param From ID type that is present in your data. Choose between "kegg", "pubchem", "chebi", "hmdb". \strong{Default = "kegg"}
#' @param To One or multiple ID types to which you want to translate your data. Choose between "kegg", "pubchem", "chebi", "hmdb". \strong{Default = c("pubchem","chebi","hmdb")}
#' @param Summary \emph{Optional: } If TRUE a long summary tables are created. \strong{Default = FALSE}
#' @param SaveAs_Table \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt". \strong{Default = "csv"}
#' @param FolderPath {Optional:} Path to the folder the results should be saved at. \strong{Default = NULL}
#'
#' @return List with at least three DFs: 1) Original data and the new column of translated ids spearated by comma. 2) Mapping information between Original ID to Translated ID. 3) Mapping summary between Original ID to Translated ID.
#'
#' @examples
#' KEGG_Pathways <- MetaProViz::LoadKEGG()
#' Res <- MetaProViz::TranslateID(InputData= KEGG_Pathways, SettingsInfo = c(InputID="MetaboliteID", GroupingVariable="term"), From = c("kegg"), To = c("pubchem", "hmdb"))
#'
#' @keywords Translate metabolite IDs
#'
#' @importFrom dplyr mutate select group_by ungroup distinct filter across n
#' @importFrom tidyselect all_of
#' @importFrom rlang !!! !! := sym syms
#' @importFrom OmnipathR id_types translate_ids
#' @importFrom logger log_warn
#' @importFrom stringr str_to_lower
#'
#' @export
#'
TranslateID <- function(InputData,
                        SettingsInfo = c(InputID="MetaboliteID", GroupingVariable="term"),
                        From = "kegg",
                        To = c("pubchem","chebi","hmdb"),
                        Summary=FALSE,
                        SaveAs_Table= "csv",
                        FolderPath=NULL
  ){# Add ability to also get metabolite names that are human readable from an ID type!

  MetaProViz_Init()

  ## ------------------  Check Input ------------------- ##
  # HelperFunction `CheckInput`
  MetaProViz:::CheckInput(InputData=InputData,
                          InputData_Num=FALSE,
                          SaveAs_Table=SaveAs_Table)

  # Specific checks:
  if("InputID" %in% names(SettingsInfo)){
    if(SettingsInfo[["InputID"]] %in% colnames(InputData)== FALSE){
      message <- paste0("The ", SettingsInfo[["InputID"]], " column selected as InputID in SettingsInfo was not found in InputData. Please check your input.")
      logger::log_trace(paste("Error ", message, sep=""))
      stop(message)
    }
  }

  if("GroupingVariable" %in% names(SettingsInfo)){
    if(SettingsInfo[["GroupingVariable"]] %in% colnames(InputData)== FALSE){
      message <- paste0("The ", SettingsInfo[["GroupingVariable"]], " column selected as GroupingVariable in SettingsInfo was not found in InputData. Please check your input.")
      logger::log_trace(paste("Error ", message, sep=""))
      stop(message)
    }
  }

  if(is.logical(Summary) == FALSE){
    message <- paste0("Check input. The Summary parameter should be either =TRUE or =FALSE.")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }

  unknown_types <- OmnipathR::id_types() %>%
    dplyr::select(tidyselect::starts_with('in_')) %>%
    unlist %>%
    unique %>%
    str_to_lower %>%
    setdiff(union(From, To), .)

  if (length(unknown_types) > 0L) {
    msg <- sprintf('The following ID types are not recognized: %s', paste(unknown_types, collapse = ', '))
    logger::log_warn(msg)
    warning(msg)
  }

  # Check that SettingsInfo[['InputID']] has no duplications within one group --> should not be the case --> remove duplications and inform the user/ ask if they forget to set groupings column
  doublons <- InputData %>%
    dplyr::group_by(!!sym(SettingsInfo[['InputID']]), !!sym(SettingsInfo[['GroupingVariable']]))%>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::ungroup()

  if(nrow(doublons) > 0){
    message <- sprintf('The following ID types are duplicated within one group: %s',paste(doublons, collapse = ', '))
    logger::log_warn(message)
    warning(message)
  }

  ## ------------------  Create output folders and path ------------------- ##
  if(is.null(SaveAs_Table)==FALSE ){
    Folder <- MetaProViz:::SavePath(FolderName= "PriorKnowledge",
                                    FolderPath=FolderPath)

    SubFolder <- file.path(Folder, "ID_Translation")
    if (!dir.exists(SubFolder)) {dir.create(SubFolder)}
  }

  ########################################################################################################################################################
  ## ------------------ Translate To-From for each pair ------------------- ##
  TranslatedDF <- OmnipathR::translate_ids(
      InputData,
      !!sym(SettingsInfo[['InputID']]) :=  !!sym(From),
      !!!syms(To),#list of symbols, hence three !!!
      ramp = TRUE,
      expand = FALSE,
      quantify_ambiguity = TRUE,
      qualify_ambiguity = TRUE,
      ambiguity_groups =  SettingsInfo[['GroupingVariable']],#Checks within the groups, without it checks across groups
      ambiguity_summary = TRUE
    )
  #TranslatedDF %>% attributes %>% names
  #TranslatedDF%>% attr('ambiguity_MetaboliteID_hmdb')

  ## --------------- Create output DF -------------------- ##
  ResList <- list()

  ## Create DF for TranslatedIDs only with the original data and the translatedID columns
  DF_subset <- TranslatedDF %>%
    dplyr::select(tidyselect::all_of(intersect(names(.), names(InputData))), tidyselect::all_of(To)) %>%
    dplyr::mutate(across(all_of(To), ~ map_chr(., ~ paste(unique(.), collapse = ", ")))) %>%
    dplyr::group_by(!!sym(SettingsInfo[['InputID']]), !!sym(SettingsInfo[['GroupingVariable']])) %>%
    dplyr::mutate(across(tidyselect::all_of(To), ~ paste(unique(.), collapse = ", "), .names = "{.col}")) %>%
    dplyr::ungroup() %>%
    dplyr::distinct() %>%
    dplyr::mutate(dplyr::across(tidyselect::all_of(To), ~ ifelse(. == "0", NA, .)))

  ResList[["TranslatedDF"]] <- DF_subset

  ## Add DF with mapping information
  ResList[["TranslatedDF_MappingInfo"]] <- TranslatedDF

  ## Also save the different mapping summaries!
  for(item in To){
    SummaryDF <- TranslatedDF%>% attr(paste0("ambiguity_", SettingsInfo[['InputID']], "_", item, sep=""))
    ResList[[paste0("MappingSummary_", item, sep="")]] <-  SummaryDF
  }

  ## Create the long DF summary if Summary =TRUE
  if(Summary==TRUE){
    for(item in To){
      Summary <- MetaProViz::MappingAmbiguity(InputData= TranslatedDF,
                                            From = SettingsInfo[['InputID']],
                                            To = item,
                                            GroupingVariable = SettingsInfo[['GroupingVariable']],
                                            Summary=TRUE)[["Summary"]]
      ResList[[paste0("MappingSummary_Long_", From, "-to-", item, sep="")]] <- Summary
    }
  }


  ## ------------------ Save the results ------------------- ##
  suppressMessages(suppressWarnings(
    MetaProViz:::SaveRes(InputList_DF=ResList,
                         InputList_Plot= NULL,
                         SaveAs_Table=SaveAs_Table,
                         SaveAs_Plot=NULL,
                         FolderPath= SubFolder,
                         FileName= "TranslateID",
                         CoRe=FALSE,
                         PrintPlot=FALSE)))

  #Return
  invisible(return(ResList))
}



##########################################################################################
### ### ### Find additional potential IDs  ### ### ###
##########################################################################################

#' Find additional potential IDs for  "kegg", "pubchem", "chebi", "hmdb"
#'
#' @param InputData Dataframe with at least one column with the detected metabolite IDs (one ID per row).
#' @param SettingsInfo \emph{Optional: } Column name of metabolite IDs. \strong{Default = list(InputID="MetaboliteID")}
#' @param From ID type that is present in your data. Choose between "kegg", "pubchem", "chebi", "hmdb". \strong{Default = "kegg"}
#' @param SaveAs_Table \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt". \strong{Default = "csv"}
#' @param FolderPath {Optional:} Path to the folder the results should be saved at. \strong{Default = NULL}
#'
#' @return Input DF with additional column including potential additional IDs.
#'
#' @examples
#' DetectedIDs <- MetaProViz::ToyData(Data="Cells_MetaData")%>% tibble::rownames_to_column("TrivialName")%>%tidyr::drop_na()
#' Res <- MetaProViz::EquivalentIDs(InputData= DetectedIDs, SettingsInfo = c(InputID="KEGG.ID"), From = "kegg")
#'
#' @keywords Find potential additional IDs for one metabolite identifier
#'
#' @importFrom dplyr mutate select group_by ungroup distinct filter across rowwise
#' @importFrom tidyr separate_rows
#' @importFrom purrr map_chr
#' @importFrom tidyselect all_of
#' @importFrom rlang !!! !! := sym syms
#' @importFrom OmnipathR id_types translate_ids
#' @importFrom logger log_warn log_trace
#' @importFrom stringr str_to_lower
#'
#' @export
#'
EquivalentIDs <- function(InputData,
                          SettingsInfo = c(InputID="MetaboliteID"),
                          From = "kegg",
                          SaveAs_Table= "csv",
                          FolderPath=NULL){

  MetaProViz_Init()

  ## ------------------  Check Input ------------------- ##
  # HelperFunction `CheckInput`
  MetaProViz:::CheckInput(InputData=InputData,
                          InputData_Num=FALSE,
                          SaveAs_Table=SaveAs_Table)

  # Specific checks:
  if("InputID" %in% names(SettingsInfo)){
    if(SettingsInfo[["InputID"]] %in% colnames(InputData)== FALSE){
      message <- paste0("The ", SettingsInfo[["InputID"]], " column selected as InputID in SettingsInfo was not found in InputData. Please check your input.")
      logger::log_trace(paste("Error ", message, sep=""))
      stop(message)
    }
  }

  unknown_types <- OmnipathR::id_types() %>%
    dplyr::select(tidyselect::starts_with('in_')) %>%
    unlist %>%
    unique %>%
    str_to_lower %>%
    setdiff(From, .)

  if (length(unknown_types) > 0L) {
    msg <- sprintf('The following ID types are not recognized: %s', paste(unknown_types, collapse = ', '))
    logger::log_warn(msg)
    warning(msg)
  }

  # Check that SettingsInfo[['InputID']] has no duplications within one group --> should not be the case --> remove duplications and inform the user/ ask if they forget to set groupings column
  doublons <- InputData[duplicated(InputData[[SettingsInfo[['InputID']]]]), ]

  if(nrow(doublons) > 0){
    InputData <- InputData %>%
      dplyr::distinct(!!sym(SettingsInfo[['InputID']]), .keep_all = TRUE)

    message <- sprintf('The following IDs are duplicated and removed: %s',paste(doublons[[SettingsInfo[['InputID']]]], collapse = ', '))
    logger::log_warn(message)
    warning(message)
  }

  ## ------------------  Create output folders and path ------------------- ##
  if(is.null(SaveAs_Table)==FALSE ){
    Folder <- MetaProViz:::SavePath(FolderName= "PriorKnowledge",
                                    FolderPath=FolderPath)

    SubFolder <- file.path(Folder, "EquivalentIDs")
    if (!dir.exists(SubFolder)) {dir.create(SubFolder)}
  }

  ## ------------------ Set the ID type for To ----------------- ##
  To <- case_when(
    From == "pubchem" ~ "chebi",  # If To is "pubchem", choose "chebi"
    TRUE ~ "pubchem"              # For other cases, don't use a secondary column
  )

  message <- paste0(To, " is used to find additional potential IDs for ", From, ".", sep="")
  logger::log_trace(message)
  message(message)

  ## ------------------ Translate From-to-To ------------------- ##
  TranslatedDF <- OmnipathR::translate_ids(
    InputData,
    !!sym(SettingsInfo[['InputID']]) :=  !!sym(From),
    !!!syms(To),#list of symbols, hence three !!!
    ramp = TRUE,
    expand = FALSE,
    quantify_ambiguity =FALSE,
    qualify_ambiguity =  TRUE, # Can not be set to FALSE!
    ambiguity_groups =  NULL,#Checks within the groups, without it checks across groups
    ambiguity_summary =  FALSE
  )%>%
    dplyr::select(tidyselect::all_of(intersect(names(.), names(InputData))), tidyselect::all_of(To)) %>%
    dplyr::mutate(across(all_of(To), ~ purrr::map_chr(., ~ paste(unique(.), collapse = ", ")))) %>%
    dplyr::group_by(!!sym(SettingsInfo[['InputID']])) %>%
    dplyr::mutate(across(tidyselect::all_of(To), ~ paste(unique(.), collapse = ", "), .names = "{.col}")) %>%
    dplyr::ungroup() %>%
    dplyr::distinct() %>%
    dplyr::mutate(dplyr::across(tidyselect::all_of(To), ~ ifelse(. == "0", NA, .)))


  ## ------------------ Translate To-to-From ------------------- ##
  TranslatedDF_Long <- TranslatedDF%>%
    dplyr::select(!!sym(SettingsInfo[['InputID']]), !!sym(To))%>%
    dplyr::rename("InputID" = !!sym(SettingsInfo[['InputID']]))%>%
    tidyr::separate_rows(!!sym(To), sep = ", ") %>%
    dplyr::mutate(across(all_of(To), ~trimws(.))) %>%  # Remove extra spaces
    dplyr::filter(!!sym(To) != "")  # Remove empty entries

  OtherIDs <- OmnipathR::translate_ids(
    TranslatedDF_Long ,
    !!sym(To),
    !!sym(From),#list of symbols, hence three !!!
    ramp = TRUE,
    expand = FALSE,
    quantify_ambiguity =FALSE,
    qualify_ambiguity =  TRUE, # Can not be set to FALSE!
    ambiguity_groups =  NULL,#Checks within the groups, without it checks across groups
    ambiguity_summary =  FALSE
  )%>%
    dplyr::select("InputID", !!sym(To), !!sym(From))%>%
    dplyr::distinct(InputID, !!sym(From), .keep_all = TRUE) %>%  # Remove duplicates based on InputID and From
    dplyr::mutate(AdditionalID = dplyr::if_else(InputID == !!sym(From), FALSE, TRUE)) %>%
    dplyr::select("InputID",!!sym(From), "AdditionalID")%>%
    dplyr::filter(AdditionalID == TRUE) %>%
    dplyr::mutate(across(all_of(From), ~ purrr::map_chr(., ~ paste(unique(.), collapse = ", "))))%>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      FromList = list(stringr::str_split(!!sym(From), ",\\s*")[[1]]),  # Wrap in list
      SameAsInput = ifelse(any(FromList == InputID), InputID, NA_character_),  # Match InputID
      PotentialAdditionalIDs = paste(FromList[FromList != InputID], collapse = ", ")  # Combine other IDs
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(InputID, PotentialAdditionalIDs)  # Final selection

  ## ------------------ Create Output ------------------- ##
  OutputDF <- merge(InputData, OtherIDs, by.x= SettingsInfo[['InputID']] , by.y= "InputID", all.x=TRUE)

  ## ------------------ Save the results ------------------- ##
  ResList <- list("EquivalentIDs" = OutputDF)

  suppressMessages(suppressWarnings(
    MetaProViz:::SaveRes(InputList_DF=ResList,
                         InputList_Plot= NULL,
                         SaveAs_Table=SaveAs_Table,
                         SaveAs_Plot=NULL,
                         FolderPath= SubFolder,
                         FileName= "EquivalentIDs",
                         CoRe=FALSE,
                         PrintPlot=FALSE)))

  return(invisible(OutputDF))
}

##########################################################################################
### ### ### Mapping Ambiguity ### ### ###
##########################################################################################

#' Create Mapping Ambiguities between two ID types
#'
#' @param InputData Translated DF from MetaProViz::TranslateID reults or Dataframe with at least one column with the target metabolite ID and another MetaboliteID type. One of the IDs can only have one ID per row, the other ID can be either separated by comma or a list. Optional: add other columns such as source (e.g. term).
#' @param To Column name of original metabolite identifier in InputData. Here should only have one ID per row.
#' @param From Column name of the secondary or translated metabolite identifier in InputData. Here can be multiple IDs per row either separated by comma " ," or a list of IDs.
#' @param GroupingVariable \emph{Optional: } If NULL no groups are used. If TRUE provide column name in InputData containing the GroupingVariable and features are grouped. \strong{Default = NULL}
#' @param Summary \emph{Optional: } If TRUE a long summary tables are created. \strong{Default = FALSE}
#' @param SaveAs_Table \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt". \strong{Default = "csv"}
#' @param FolderPath {Optional:} Path to the folder the results should be saved at. \strong{Default = NULL}
#'
#' @return List with at least 4 DFs: 1-3) From-to-To: 1. MappingIssues, 2. MappingIssues Summary, 3. Long summary (If Summary=TRUE) & 4-6) To-to-From: 4. MappingIssues, 5. MappingIssues Summary, 6. Long summary (If Summary=TRUE) & 7) Combined summary table (If Summary=TRUE)
#'
#' @examples
#' KEGG_Pathways <- MetaProViz::LoadKEGG()
#' InputDF <- MetaProViz::TranslateID(InputData= KEGG_Pathways, SettingsInfo = c(InputID="MetaboliteID", GroupingVariable="term"), From = c("kegg"), To = c("pubchem"))[["TranslatedDF"]]
#' Res <- MetaProViz::MappingAmbiguity(InputData= InputDF, From = "MetaboliteID", To = "pubchem", GroupingVariable = "term", Summary=TRUE)
#'
#' @keywords Mapping ambiguity
#'
#' @importFrom dplyr mutate
#' @importFrom rlang !!! !! := sym syms
#' @importFrom OmnipathR ambiguity
#'
#' @export
#'
MappingAmbiguity <- function(InputData,
                             From,
                             To,
                             GroupingVariable = NULL,
                             Summary=FALSE,
                             SaveAs_Table= "csv",
                             FolderPath=NULL
) {

  MetaProViz_Init()
  ## ------------------  Check Input ------------------- ##
  # HelperFunction `CheckInput`
  MetaProViz:::CheckInput(InputData=InputData,
                          InputData_Num=FALSE,
                          SaveAs_Table=SaveAs_Table)

  # Specific checks:
  if(From %in% colnames(InputData)== FALSE){
      message <- paste0(From, " column was not found in InputData. Please check your input.")
      logger::log_trace(paste("Error ", message, sep=""))
      stop(message)
      }

  if(To %in% colnames(InputData)== FALSE){
    message <- paste0(To, " column was not found in InputData. Please check your input.")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }

  if(is.null(GroupingVariable)==FALSE){
    if(GroupingVariable %in% colnames(InputData)== FALSE){
      message <- paste0(GroupingVariable, " column was not found in InputData. Please check your input.")
      logger::log_trace(paste("Error ", message, sep=""))
      stop(message)
    }
  }

  if(is.logical(Summary) == FALSE){
    message <- paste0("Check input. The Summary parameter should be either =TRUE or =FALSE.")
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }

  ## ------------------  General checks of wrong occurences ------------------- ##
  # Task 1: Check that From has no duplications within one group --> should not be the case --> remove duplications and inform the user/ ask if they forget to set groupings column
  # Task 2: Check that From has the same items in to across the different entries (would be in different Groupings, otherwise there should not be any duplications) --> List of Miss-Mappings across terms

  # FYI: The above can not happen if our translateID function was used, but may be the case when the user has done something manually before


  ## ------------------  Create output folders and path ------------------- ##
  if(is.null(SaveAs_Table)==FALSE ){
    Folder <- MetaProViz:::SavePath(FolderName= "PriorKnowledge",
                                    FolderPath=FolderPath)

    SubFolder <- file.path(Folder, "MappingAmbiguities")
    if (!dir.exists(SubFolder)) {dir.create(SubFolder)}
  }

  #####################################################################################################################################################################################
  ## ------------------  Prepare Input data ------------------- ##
  #If the user provides a DF where the To column is a list of IDs, then we can use it right away
  #If the To column is not a list of IDs, but a character column, we need to convert it into a list of IDs
  if(is.character(InputData[[To]])==TRUE){
    InputData[[To]] <- InputData[[To]]%>%
      strsplit(", ")%>%
      lapply(as.character)
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
    #Run Omnipath ambiguity
    ResList[[paste0(Comp[[comp]]$From, "-to-", Comp[[comp]]$To , sep="")]] <- InputData %>%
      tidyr::unnest(cols = all_of(Comp[[comp]]$From))%>% # unlist the columns in case they are not expaned
      filter(!is.na(!!sym(Comp[[comp]]$From)))%>%#Remove NA values, otherwise they are counted as column is character
      OmnipathR::ambiguity(
          from_col = !!sym(Comp[[comp]]$From),
          to_col = !!sym(Comp[[comp]]$To),
          groups = GroupingVariable,
          quantify = TRUE,
          qualify = TRUE,
          global = TRUE,#across groups will be done additionally --> suffix _AcrossGroup
          summary=TRUE, #summary of the mapping column
          expand=TRUE)

    #Extract summary table:
    ResList[[paste0(Comp[[comp]]$From, "-to-", Comp[[comp]]$To, "_Summary", sep="")]] <-
        ResList[[paste0(Comp[[comp]]$From, "-to-", Comp[[comp]]$To , sep="")]]%>%
        attr(paste0("ambiguity_", Comp[[comp]]$From , "_",Comp[[comp]]$To, sep=""))

    ############################################################################################################
    if(Summary==TRUE){
      if(is.null(GroupingVariable)==FALSE){
        # Add further information we need to summarise the table and combine Original-to-Translated and Translated-to-Original
        # If we have a GroupingVariable we need to combine it with the MetaboliteID before merging
        ResList[[paste0(Comp[[comp]]$From, "-to-", Comp[[comp]]$To, "_Long", sep="")]] <- ResList[[paste0(Comp[[comp]]$From, "-to-", Comp[[comp]]$To , sep="")]]%>%
          tidyr::unnest(cols = all_of(Comp[[comp]]$From))%>%
          mutate(!!sym(paste0("AcrossGroupMappingIssue(", Comp[[comp]]$From, "_to_", Comp[[comp]]$To, ")", sep="")) := case_when(
            !!sym(paste0(Comp[[comp]]$From, "_", Comp[[comp]]$To, "_ambiguity_bygroup", sep="")) != !!sym(paste0(Comp[[comp]]$From, "_", Comp[[comp]]$To, "_ambiguity", sep=""))  ~ "TRUE",
            TRUE ~ "FALSE" ))%>%
          group_by(!!sym(Comp[[comp]]$From), !!sym(GroupingVariable))%>%
          mutate(!!sym(Comp[[comp]]$To) := ifelse(!!sym(Comp[[comp]]$From) == 0, NA,  # Or another placeholder
                                                  paste(unique(!!sym(Comp[[comp]]$To)), collapse = ", ")
          )) %>%
          mutate( !!sym(paste0("Count(", Comp[[comp]]$From, "_to_", Comp[[comp]]$To, ")")) := ifelse(all(!!sym(Comp[[comp]]$To) == 0), 0, n()))%>%
          ungroup()%>%
          distinct() %>%
          unite(!!sym(paste0(Comp[[comp]]$From, "_to_", Comp[[comp]]$To)), c(Comp[[comp]]$From, Comp[[comp]]$To), sep=" --> ", remove=FALSE)%>%
          separate_rows(!!sym(Comp[[comp]]$To), sep = ", ") %>%
          unite(UniqueID, c(From, To, GroupingVariable), sep="_", remove=FALSE)%>%
          distinct()
        }else{
          ResList[[paste0(Comp[[comp]]$From, "-to-", Comp[[comp]]$To, "_Long", sep="")]] <- ResList[[paste0(Comp[[comp]]$From, "-to-", Comp[[comp]]$To , sep="")]]%>%
            tidyr::unnest(cols = all_of(Comp[[comp]]$From))%>%
            group_by(!!sym(Comp[[comp]]$From))%>%
            mutate(!!sym(Comp[[comp]]$To) := ifelse(!!sym(Comp[[comp]]$From) == 0, NA,  # Or another placeholder
                                                    paste(unique(!!sym(Comp[[comp]]$To)), collapse = ", ")
            )) %>%
            mutate( !!sym(paste0("Count(", Comp[[comp]]$From, "_to_", Comp[[comp]]$To, ")")) := ifelse(all(!!sym(Comp[[comp]]$To) == 0), 0, n()))%>%
            ungroup()%>%
            distinct() %>%
            unite(!!sym(paste0(Comp[[comp]]$From, "_to_", Comp[[comp]]$To)), c(Comp[[comp]]$From, Comp[[comp]]$To), sep=" --> ", remove=FALSE)%>%
            separate_rows(!!sym(Comp[[comp]]$To), sep = ", ") %>%
            unite(UniqueID, c(From, To), sep="_", remove=FALSE)%>%
            distinct()%>%
            mutate(!!sym(paste0("AcrossGroupMappingIssue(", From, "_to_", To, ")", sep="")) := NA)
        }
    }

    # Add NA metabolite maps back if they do exist:
    Removed <- InputData %>%
      tidyr::unnest(cols = all_of(Comp[[comp]]$From))%>% # unlist the columns in case they are not expaned
      filter(is.na(!!sym(Comp[[comp]]$From)))
    if(nrow(Removed)>0){
      ResList[[paste0(Comp[[comp]]$From, "-to-", Comp[[comp]]$To , sep="")]] <- bind_rows(ResList[[paste0(Comp[[comp]]$From, "-to-", Comp[[comp]]$To , sep="")]],
                                                                                         test<- Removed%>%
                                                                                            bind_cols(setNames(as.list(rep(NA, length(setdiff(names(ResList[[paste0(Comp[[comp]]$From, "-to-", Comp[[comp]]$To , sep="")]]), names(Removed))))),
                                                                                                               setdiff(names(ResList[[paste0(Comp[[comp]]$From, "-to-", Comp[[comp]]$To , sep="")]]), names(Removed))))
      )
    }
  }

  ## ------------------ Create SummaryTable ------------------- ##
  if(Summary==TRUE){
    # Combine the two tables
    Summary <- merge(x= ResList[[paste0(From, "-to-", To, "_Long", sep="")]][,c("UniqueID", paste0(From, "_to_", To), paste0("Count(", From, "_to_", To, ")"), paste0("AcrossGroupMappingIssue(", From, "_to_", To, ")", sep=""))],
                     y= ResList[[paste0(To, "-to-", From, "_Long", sep="")]][,c("UniqueID", paste0(To, "_to_", From), paste0("Count(", To, "_to_", From, ")"), paste0("AcrossGroupMappingIssue(", To, "_to_", From, ")", sep=""))],
                     by = "UniqueID",
                     all = TRUE)%>%
      separate(UniqueID, into = c(From, To, GroupingVariable), sep="_", remove=FALSE)%>%
      distinct()

    # Add relevant mapping information
    Summary <- Summary %>%
      mutate(Mapping = case_when(
        !!sym(paste0("Count(", From, "_to_", To, ")")) == 1 & !!sym(paste0("Count(", To, "_to_", From, ")")) == 1  ~ "one-to-one",
        !!sym(paste0("Count(", From, "_to_", To, ")")) > 1 & !!sym(paste0("Count(", To, "_to_", From, ")")) == 1  ~ "one-to-many",
        !!sym(paste0("Count(", From, "_to_", To, ")")) > 1 & !!sym(paste0("Count(", To, "_to_", From, ")")) > 1  ~ "many-to-many",
        !!sym(paste0("Count(", From, "_to_", To, ")")) == 1 & !!sym(paste0("Count(", To, "_to_", From, ")")) > 1  ~ "many-to-one",
        !!sym(paste0("Count(", From, "_to_", To, ")")) >= 1 & !!sym(paste0("Count(", To, "_to_", From, ")")) == NA  ~ "one-to-none",
        !!sym(paste0("Count(", From, "_to_", To, ")")) >= 1 & is.na(!!sym(paste0("Count(", To, "_to_", From, ")")))  ~ "one-to-none",
        !!sym(paste0("Count(", From, "_to_", To, ")")) == NA & !!sym(paste0("Count(", To, "_to_", From, ")")) >= 1  ~ "none-to-one",
        is.na(!!sym(paste0("Count(", From, "_to_", To, ")"))) & !!sym(paste0("Count(", To, "_to_", From, ")")) >= 1  ~ "none-to-one",
        TRUE ~ NA )) %>%
      mutate( !!sym(paste0("Count(", From, "_to_", To, ")")) := replace_na( !!sym(paste0("Count(", From, "_to_", To, ")")), 0)) %>%
      mutate( !!sym(paste0("Count(", To, "_to_", From, ")")) := replace_na( !!sym(paste0("Count(", To, "_to_", From, ")")), 0))

    ResList[["Summary"]] <- Summary
  }

  ## ------------------ Save the results ------------------- ##
  suppressMessages(suppressWarnings(
    MetaProViz:::SaveRes(InputList_DF=ResList,
                         InputList_Plot= NULL,
                         SaveAs_Table=SaveAs_Table,
                         SaveAs_Plot=NULL,
                         FolderPath= SubFolder,
                         FileName= "MappingAmbiguity",
                         CoRe=FALSE,
                         PrintPlot=FALSE)))

  #Return
  invisible(return(ResList))
}


##########################################################################################
### ### ### Reduce mapping ambiguities by using detected IDs  ### ### ###
##########################################################################################

#' Reduce mapping ambiguities by using detected IDs
#'
#' @param InputData Dataframe with at least one column with the target (e.g. metabolite)
#' @param PriorKnowledge Dataframe with at least one column with the "original target metabolite ID", "source" (e.g. term) and the "translated target metabolite ID".
#' @param SettingsInfo
#'
#' @examples
#' DetectedIDs <-  MetaProViz::ToyData(Data="Cells_MetaData")%>% rownames_to_column("Metabolite") %>%select("Metabolite", "HMDB")
#' PathwayFile <- MetaProViz::TranslateID(InputData= MetaProViz::LoadKEGG(), SettingsInfo = c(InputID="MetaboliteID", GroupingVariable="term"), From = c("kegg"), To = c("hmdb"))[["TranslatedDF"]]
#' Res <- MetaProViz::CleanMapping(InputData= DetectedIDs, PriorKnowledge= PathwayFile, DetectedID="HMDB", GroupingVariable="term", From="MetaboliteID", To="hmdb")
#'
#' @return Prior Knowledge usable for enrichment analysis
#'
#' @noRd
#'
CleanMapping <- function(InputData,
                         PriorKnowledge,
                         DetectedID="HMDB", #enable that lists can be passed if multuiple IDs are assigned to one measurement!
                         GroupingVariable="term",
                         From="MetaboliteID",
                         To="hmdb",
                         SaveAs_Table= "csv",
                         FolderPath=NULL
){

  MetaProViz_Init()
  ## ------------------ Check Input ------------------- ##
  # HelperFunction `CheckInput`

  # Function that can use Prior knowledge that has multiple IDs "X1, X2, X3" and use the measured features to remove IDs based on measured data.
  # This is to ensure that not two detected metabolites map to the same entry and if the original PK was translated to  ensure the one-to-many, many-to-one issues are taken care of (within and across DBs)

  # --> Cleans translated ID in prior knowledge based on measured features

  ## ------------------  Create output folders and path ------------------- ##
  if(is.null(SaveAs_Table)==FALSE ){
    Folder <- MetaProViz:::SavePath(FolderName= "PriorKnowledge",
                                    FolderPath=FolderPath)
  }

}


##########################################################################################
### ### ### Check Measured ID's in prior knowledge ### ### ###
##########################################################################################

#' Check and summarize PriorKnowledge-to-MeasuredFeatures relationship
#'
#' @param InputData Dataframe with at least one column with the target (e.g. metabolite), you can add other columns such as source (e.g. term)
#' @param SettingsInfo
#'
#' @importFrom dplyr mutate
#' @importFrom rlang !!! !! := sym syms
#'
#' @noRd
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
    SubFolder <- file.path(Folder, "MetaboliteSet")
    if (!dir.exists(SubFolder)) {dir.create(SubFolder)}
  }

  ######################################################################################################################################
  ## ------------------ Prepare InputData -------------------##
  DetectedDF <- InputData%>%
    select(!!sym(DetectedID))%>%
    as.data.frame()%>%
    distinct()%>%
    filter(!is.na(!!sym(DetectedID)))


  ## ------------------ Use the ambiguity function to create the "Summary File" ------------------- ##
  Summary <- MetaProViz::MappingAmbiguity(InputData= PriorKnowledge,
                                          From = From,
                                          To = To,
                                          GroupingVariable = GroupingVariable,
                                          Summary=TRUE)[["Summary"]]
    #remove To=0 or NA


  ## ------------------ Extract IDs that "are not" versus "are" present in prior knowledge and return ----------------------##
  # Add detected IDs to the PriorKnowledge
  CombinedDF_Detected <- merge(x= Summary,
                y= DetectedDF,
                by.x= To,
                by.y=DetectedID,
                all.y=TRUE)

  DetectedIDs_NoPathway <- CombinedDF_Detected%>%
    filter(is.na(!!sym(From)))%>%
    select(!!sym(To))%>%
    distinct()

  DetectedIDs_InPathway <- CombinedDF_Detected%>%
    filter(!is.na(!!sym(From)))%>%
    select(!!sym(To))%>%
    distinct()

  message <- paste0("Of ", nrow(DetectedDF), " metabolite IDs that were detected in the measured data, ", nrow(DetectedIDs_InPathway), " were found in the prior knowledge and ", nrow(DetectedIDs_NoPathway), " were not found in the prior knowledge.")
  logger::log_info(message)
  message(message)

  ## ------------------ Clean Prior Knowledge based on detected IDs ----------------------##
  CombinedDF_Detected <- CombinedDF_Detected %>%
    filter(!is.na(!!sym(From)))%>%
    group_by(!!sym(From), !!sym(GroupingVariable))%>%
    mutate(Collapsed_To = paste(unique(!!sym(To)), collapse = ", "))%>%
    mutate(Collapsed_To_Number = n())%>%
    mutate(Action = case_when(
      (!!sym(To) != 0 | !is.na(!!sym(To))) & Mapping == "one-to-one" & dplyr::if_all(starts_with("AcrossGroupMappingIssue("), ~ . != TRUE, na.rm = TRUE) & Collapsed_To_Number == 1 ~ "None",# Detected ID is a one-to-one mapping within and across groups.
      (!!sym(To) != 0 | !is.na(!!sym(To))) & Mapping == "one-to-many" & dplyr::if_all(starts_with("AcrossGroupMappingIssue("), ~ . != TRUE, na.rm = TRUE) & Collapsed_To_Number == 1 ~ "None",
      (!!sym(To) != 0 | !is.na(!!sym(To))) & Mapping == "one-to-many" & dplyr::if_all(starts_with("AcrossGroupMappingIssue("), ~ . != TRUE, na.rm = TRUE) & Collapsed_To_Number > 1 ~ "one-to-many_detected",
      (!!sym(To) != 0 | !is.na(!!sym(To))) & Mapping == "many-to-many" & dplyr::if_all(starts_with("AcrossGroupMappingIssue("), ~ . != TRUE, na.rm = TRUE) & Collapsed_To_Number == 1 ~ "many-to-many_1detected",#Check that the detected ID is not the one mapping to many!
      (!!sym(To) != 0 | !is.na(!!sym(To))) & Mapping == "many-to-many" & dplyr::if_all(starts_with("AcrossGroupMappingIssue("), ~ . != TRUE, na.rm = TRUE) & Collapsed_To_Number > 1 ~ "many-to-many_Ndetected",#Check that the detected ID is not the one mapping to many!
      #if_any(starts_with("AcrossGroupMappingIssue("), ~ . == TRUE, na.rm = TRUE)~ "TRUE",
      TRUE ~ "FALSE" ))

  #Clean Prior knowledge:






  ## ------------------ Save Results ----------------------##

}


##########################################################################################
### ### ###  ### ### ###
##########################################################################################

#' Deal with detected input features. If we only have one ID, could it also be another one?
#'
#'
#AssignID <- function(){

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

#}


##########################################################################################
### ### ### Cluster Prior Knowledge ### ### ###
##########################################################################################

#' Deal with pathway overlap in prior knowledge
#'
#' @param InputData Dataframe with at least one column with the target (e.g. metabolite) and a column source (e.g. term).
#' @param SettingsInfo = c(InputID="MetaboliteID", GroupingVariable="term"),
#'
#' @examples
#' KEGG_Pathways <- MetaProViz::LoadKEGG()
#' InputData = KEGG_Pathways
#'
#'
#' @noRd



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



  ######################################################################################################################################
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
#' @noRd

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




