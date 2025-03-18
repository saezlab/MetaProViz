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
#' KEGG_Pathways <- LoadKEGG()
#' Res <- TranslateID(InputData= KEGG_Pathways, 
#'     SettingsInfo = c(InputID="MetaboliteID", GroupingVariable="term"), 
#'     From = c("kegg"), To = c("pubchem", "hmdb"))
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
    SettingsInfo = c(InputID = "MetaboliteID", GroupingVariable = "term"),
    From = "kegg", ## EDIT: name the options here and use match.arg
    To = c("pubchem", "chebi", "hmdb"), ## EDIT: are here multiple allowed?, should be checked with match.arg
    Summary = FALSE,
    SaveAs_Table = "csv", ## EDIT: name the options here and use match.arg
    FolderPath = NULL) {
    
    ## add ability to also get metabolite names that are human readable from an ID type!

    MetaProViz_Init()

    ## ------------------  Check Input ------------------- ##
    # HelperFunction `CheckInput`
    CheckInput(InputData = InputData, InputData_Num = FALSE, 
        SaveAs_Table = SaveAs_Table)

    ## Specific checks:
    if ("InputID" %in% names(SettingsInfo)) {
        if (!SettingsInfo[["InputID"]] %in% colnames(InputData)) {
            message <- paste0("The ", SettingsInfo[["InputID"]], 
                " column selected as InputID in SettingsInfo was not found in ",
                "InputData. Please check your input.")
            logger::log_trace(paste0("Error ", message))
            stop(message)
        }
    }

    if ("GroupingVariable" %in% names(SettingsInfo)) {
        if (!SettingsInfo[["GroupingVariable"]] %in% colnames(InputData)) {
            message <- paste0("The ", SettingsInfo[["GroupingVariable"]], 
                " column selected as GroupingVariable in SettingsInfo was not ",
                "found in InputData. Please check your input.")
            logger::log_trace(paste0("Error ", message))
            stop(message)
        }
    }

    if (!is.logical(Summary)) {
        message <- paste0("Check input. The Summary parameter should be either TRUE or FALSE.")
        logger::log_trace(paste0("Error ", message))
        stop(message)
    }

    unknown_types <- OmnipathR::id_types() %>%
        dplyr::select(tidyselect::starts_with('in_')) %>%
        unlist() %>%
        unique() %>%
        str_to_lower() %>%
        setdiff(union(From, To), .)

    if (length(unknown_types) > 0L) {
        msg <- sprintf('The following ID types are not recognized: %s', 
            paste(unknown_types, collapse = ', '))
        logger::log_warn(msg)
        warning(msg)
    }

    ## check that SettingsInfo[['InputID']] has no duplications within 
    ## one group --> should not be the case --> remove duplications and 
    ## inform the user/ ask if they forget to set groupings column
    doublons <- InputData %>%
        dplyr::group_by(!!sym(SettingsInfo[['InputID']]), 
            !!sym(SettingsInfo[['GroupingVariable']]))%>%
        dplyr::filter(dplyr::n() > 1) %>%
        dplyr::ungroup()

    if (nrow(doublons) > 0) {
        message <- sprintf(
            "The following ID types are duplicated within one group: %s",
            paste(doublons, collapse = ', '))
        logger::log_warn(message)
        warning(message)
    }

    ## ------------------  Create output folders and path ------------------- ##
    if (!is.null(SaveAs_Table)) {
        Folder <- SavePath(FolderName = "PriorKnowledge", 
            FolderPath = FolderPath)

        SubFolder <- file.path(Folder, "ID_Translation")
        if (!dir.exists(SubFolder)) {
            dir.create(SubFolder)
        }
    }

    ############################################################################
    ## ------------------ Translate To-From for each pair ------------------- ##
    TranslatedDF <- OmnipathR::translate_ids(
        InputData,
        !!sym(SettingsInfo[['InputID']]) :=  !!sym(From),
        !!!syms(To), ## list of symbols, hence three !!!
        ramp = TRUE,
        expand = FALSE,
        quantify_ambiguity = TRUE,
        qualify_ambiguity = TRUE,
        ## checks within the groups, without it checks across groups
        ambiguity_groups =  SettingsInfo[['GroupingVariable']],
        ambiguity_summary = TRUE
    )
    #TranslatedDF %>% attributes %>% names
    #TranslatedDF%>% attr('ambiguity_MetaboliteID_hmdb')

    ## --------------- Create output DF -------------------- ##
    ResList <- list()

    ## Create DF for TranslatedIDs only with the original data and the translatedID columns
    DF_subset <- TranslatedDF %>%
        dplyr::select(tidyselect::all_of(intersect(names(.), names(InputData))), 
            tidyselect::all_of(To)) %>%
        dplyr::mutate(across(all_of(To), ~ map_chr(., ~ paste(unique(.), 
            collapse = ", ")))) %>%
        dplyr::group_by(!!sym(SettingsInfo[['InputID']]), 
            !!sym(SettingsInfo[['GroupingVariable']])) %>%
        dplyr::mutate(across(
            tidyselect::all_of(To), ~ paste(unique(.), collapse = ", "), 
            .names = "{.col}")) %>%
        dplyr::ungroup() %>%
        dplyr::distinct() %>%
        dplyr::mutate(dplyr::across(tidyselect::all_of(To), ~ ifelse(. == "0", NA, .))) ## EDIT: this seems quite complicated and takes some time to understand, is it possible to simplify or add comments?

    ResList[["TranslatedDF"]] <- DF_subset

    ## add DF with mapping information
    ResList[["TranslatedDF_MappingInfo"]] <- TranslatedDF

    ## also save the different mapping summaries!
    for (item in To) {
        SummaryDF <- TranslatedDF %>% 
            attr(paste0("ambiguity_", SettingsInfo[['InputID']], "_", item))
        ResList[[paste0("MappingSummary_", item, sep)]] <-  SummaryDF
    }

    ## create the long DF summary if Summary =TRUE
    if (Summary) {
        for (item in To) {
            Summary <- MappingAmbiguity(InputData = TranslatedDF,
                From = SettingsInfo[['InputID']],
                To = item,
                GroupingVariable = SettingsInfo[['GroupingVariable']],
                Summary = TRUE)[["Summary"]]
            ResList[[paste0("MappingSummary_Long_", From, "-to-", item)]] <- Summary
        }
    }

    ## ------------------ Save the results ------------------- ##
    suppressMessages(suppressWarnings(
        SaveRes(
            InputList_DF = ResList,
            InputList_Plot = NULL,
            SaveAs_Table = SaveAs_Table,
            SaveAs_Plot = NULL,
            FolderPath = SubFolder,
            FileName = "TranslateID",
            CoRe = FALSE,
            PrintPlot = FALSE)))
    #Return
    invisible(ResList)
}

##########################################################################################
### ### ### Find additional potential IDs  ### ### ###
##########################################################################################

#' Find additional potential IDs for  "kegg", "pubchem", "chebi", "hmdb"
#'
#' @param InputData Dataframe with at least one column with the detected metabolite IDs (one ID per row).
#' @param SettingsInfo \emph{Optional: } Column name of metabolite IDs. \strong{Default = list(InputID="MetaboliteID")}
#' @param From ID type that is present in your data. Choose between "kegg", "pubchem", "chebi", "hmdb". \strong{Default = "hmdb"}
#' @param SaveAs_Table \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt". \strong{Default = "csv"}
#' @param FolderPath {Optional:} Path to the folder the results should be saved at. \strong{Default = NULL}
#'
#' @return Input DF with additional column including potential additional IDs.
#'
#' @examples
#' DetectedIDs <- ToyData(Data="Cells_MetaData") %>% 
#'     tibble::rownames_to_column("TrivialName") %>%
#'     tidyr::drop_na()
#' Res <- EquivalentIDs(InputData = DetectedIDs, 
#'     SettingsInfo = c(InputID = "HMDB"), From = "hmdb")
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
    SettingsInfo = c(InputID = "MetaboliteID"),
    From = "hmdb",
    SaveAs_Table= "csv",
    FolderPath = NULL) {
  
    # FUTURE: Once we have the structural similarity tool available in OmniPath, we can start creating this function!

    ### 1)
    ## check Measured ID's in prior knowledge

    ### 2)
    ## A user has one HMDB IDs for their measured metabolites 
    ## (one ID per measured peak) --> this is often the case as the user 
    ## either gets a trivial name and they have searched for the ID 
    ## themselves or because the facility only provides one ID at random
    # We have mapped the HMDB IDs with the pathways and 20 do not map
    # We want to check if it is because the pathways don't include them, 
    ## or because the user just gave the wrong ID by chance (i.e. They 
    ## picked D-Alanine, but the prior knowledge includes L-Alanine)

    ## Do this by using structural information via  accessing the 
    ## structural DB in OmniPath!
    ## Output is DF with the original ID column and a new column with 
    ## additional possible IDs based on structure

    ## Is it possible to do this at the moment without structures, 
    ## but by using other pior knowledge?

    MetaProViz_Init()

    ## ------------------  Check Input ------------------- ##
    ## HelperFunction `CheckInput`
    CheckInput(InputData = InputData, InputData_Num = FALSE,
        SaveAs_Table = SaveAs_Table)

    ## specific checks:
    if ("InputID" %in% names(SettingsInfo)) {
        if (!SettingsInfo[["InputID"]] %in% colnames(InputData)) {
            message <- paste0("The ", SettingsInfo[["InputID"]], 
                " column selected as InputID in SettingsInfo was not found ",
                "in InputData. Please check your input.")
            logger::log_trace(paste0("Error ", message))
            stop(message)
        }
    }

    unknown_types <- OmnipathR::id_types() %>%
        dplyr::select(tidyselect::starts_with('in_')) %>%
        unlist() %>%
        unique() %>%
        str_to_lower() %>%
        setdiff(From, .) ## EDIT: every term that is replicated should be written as a function and tested

    if (length(unknown_types) > 0L) {
        msg <- sprintf('The following ID types are not recognized: %s', 
            paste(unknown_types, collapse = ', '))
        logger::log_warn(msg)
        warning(msg)
    }

    ## check that SettingsInfo[['InputID']] has no duplications within 
    ## one group --> should not be the case --> remove duplications and inform 
    ## the user/ ask if they forget to set groupings column
    doublons <- InputData[duplicated(InputData[[SettingsInfo[['InputID']]]]), ]

    if (nrow(doublons) > 0) {
        InputData <- InputData %>%
            dplyr::distinct(!!sym(SettingsInfo[['InputID']]), .keep_all = TRUE)

        message <- sprintf("The following IDs are duplicated and removed: %s",
            paste(doublons[[SettingsInfo[['InputID']]]], collapse = ', '))
        logger::log_warn(message)
        warning(message)
    }

    ## ------------------  Create output folders and path ------------------- ##
    if (!is.null(SaveAs_Table)) {
        Folder <- SavePath(FolderName = "PriorKnowledge", FolderPath = FolderPath)

        SubFolder <- file.path(Folder, "EquivalentIDs")
        if (!dir.exists(SubFolder)) {
            dir.create(SubFolder)
        }
    }

    ## ------------------ Set the ID type for To ----------------- ##
    To <- case_when(
        ## if To is "pubchem", choose "chebi"
        From == "chebi" ~ "pubchem",
        ## for other cases, don't use a secondary column
        TRUE ~ "chebi"              
    )

    message <- paste0(To, " is used to find additional potential IDs for ", 
        From, ".")
    logger::log_trace(message)
    message(message)

    ## ------------------ Load manual table ----------------- ##
    if (From != "kegg") {
        EquivalentFeatures <- ToyData("EquivalentFeatures") %>%
            dplyr::select(From)
    }

    ## ------------------ Translate From-to-To ------------------- ##
    TranslatedDF <- OmnipathR::translate_ids(
            InputData,
            !!sym(SettingsInfo[['InputID']]) := !!sym(From),
            !!!syms(To), ## list of symbols, hence three !!!
            ramp = TRUE,
            expand = FALSE,
            quantify_ambiguity = FALSE,
            ## Can not be set to FALSE!
            qualify_ambiguity = TRUE,
            ## Checks within the groups, without it checks across groups
            ambiguity_groups = NULL,
            ambiguity_summary = FALSE) %>%
        dplyr::select(
            tidyselect::all_of(intersect(names(.), names(InputData))), 
            tidyselect::all_of(To)) %>%
        dplyr::mutate(
            across(all_of(To), ~ purrr::map_chr(., ~ paste(unique(.), collapse = ", ")))) %>%
        dplyr::group_by(!!sym(SettingsInfo[['InputID']])) %>%
        dplyr::mutate(
            across(tidyselect::all_of(To), ~ paste(unique(.), collapse = ", "), 
                .names = "{.col}")) %>%
        dplyr::ungroup() %>%
        dplyr::distinct() %>%
        dplyr::mutate(dplyr::across(tidyselect::all_of(To), ~ ifelse(. == "0", NA, .))) ## EDIT: looks quite complicated, could it be simplified or comments be added?


    ## ------------------ Translate To-to-From ------------------- ##
    TranslatedDF_Long <- TranslatedDF %>%
        dplyr::select(!!sym(SettingsInfo[['InputID']]), !!sym(To)) %>%
        dplyr::rename("InputID" = !!sym(SettingsInfo[['InputID']])) %>%
        tidyr::separate_rows(!!sym(To), sep = ", ") %>%
        ## remove extra spaces
        dplyr::mutate(across(all_of(To), ~ trimws(.))) %>%
        ## remove empty entries
        dplyr::filter(!!sym(To) != "")

    OtherIDs <- OmnipathR::translate_ids(
            TranslatedDF_Long ,
            !!sym(To),
            !!sym(From),#list of symbols, hence three !!!
            ramp = TRUE,
            expand = FALSE,
            quantify_ambiguity =FALSE,
            ## can not be set to FALSE!
            qualify_ambiguity = TRUE, 
            ## checks within the groups, without it checks across groups
            ambiguity_groups = NULL,
            ambiguity_summary = FALSE) %>%
        dplyr::select("InputID", !!sym(To), !!sym(From)) %>%
        ## remove duplicates based on InputID and From
        dplyr::distinct(InputID, !!sym(From), .keep_all = TRUE) %>%  
        dplyr::mutate(AdditionalID = dplyr::if_else(InputID == !!sym(From), FALSE, TRUE)) %>%
        dplyr::select("InputID",!!sym(From), "AdditionalID") %>%
        dplyr::filter(AdditionalID == TRUE) %>%
        dplyr::mutate(across(all_of(From), ~ purrr::map_chr(., ~ paste(unique(.), collapse = ", ")))) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
            ## wrap in list
            FromList = list(stringr::str_split(!!sym(From), ",\\s*")[[1]]),  
            ## match InputID
            SameAsInput = ifelse(any(FromList == InputID), InputID, NA_character_),  
            ## combine other IDs
            PotentialAdditionalIDs = paste(FromList[FromList != InputID], collapse = ", ")  
        ) %>%
        dplyr::ungroup() %>%
        ## final selection
        dplyr::select(InputID, PotentialAdditionalIDs, hmdb) %>%  
        dplyr::rename("AllIDs" = "hmdb") ## EDIT: looks quite complicated, could it be simplified or comments be added?

    ## ------------------ Merge to Input ------------------- ##
    OtherIDs <- merge(InputData, OtherIDs, by.x = SettingsInfo[['InputID']], 
        by.y = "InputID", all.x = TRUE)

    ##------------------- Add additional IDs -------------- ##
    if (exists("EquivalentFeatures")) {
        EquivalentFeatures$AllIDs <- EquivalentFeatures[[From]]
        EquivalentFeatures_Long <- EquivalentFeatures  %>%
            separate_rows(!!sym(From), sep = ",")

        OtherIDs <- merge(OtherIDs, EquivalentFeatures_Long, 
                by.x = SettingsInfo[['InputID']] , by.y = "hmdb", 
                all.x = TRUE) %>%
            rowwise() %>%
            mutate(AllIDs = paste(unique(
                na.omit(unlist(str_split(paste(na.omit(c(AllIDs.x, AllIDs.y)), collapse = ","), ",\\s*")))),
                collapse = ",")) %>%
            ungroup()%>%
            rowwise() %>%
            mutate(
            PotentialAdditionalIDs = paste(
                setdiff(
                    ## split merged_column into individual IDs
                    unlist(str_split(AllIDs, ",\\s*")),
                    ## split hmdb into individual IDs
                    as.character(!!sym(SettingsInfo[['InputID']]))
                ),
                ## Combine the remaining IDs back into a comma-separated string
                collapse = ", ")) %>%
            ungroup() %>%
            select(-AllIDs.x, -AllIDs.y) ## EDIT: looks quite complicated, could it be simplified or comments be added?
    }

    ## ------------------ Create Output ------------------- ##
    OutputDF <- OtherIDs

    ## ------------------ Save the results ------------------- ##
    ResList <- list("EquivalentIDs" = OutputDF)

    suppressMessages(suppressWarnings(
        SaveRes(
            InputList_DF = ResList,
            InputList_Plot = NULL,
            SaveAs_Table = SaveAs_Table,
            SaveAs_Plot = NULL,
            FolderPath = SubFolder,
            FileName = "EquivalentIDs",
            CoRe = FALSE,
            PrintPlot = FALSE)))

    invisible(OutputDF)
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
#' KEGG_Pathways <- LoadKEGG()
#' InputDF <- TranslateID(InputData= KEGG_Pathways, 
#'     SettingsInfo = c(InputID="MetaboliteID", GroupingVariable="term"), 
#'     From = c("kegg"), To = c("pubchem"))[["TranslatedDF"]]
#' Res <- MappingAmbiguity(InputData = InputDF, From = "MetaboliteID", 
#'     To = "pubchem", GroupingVariable = "term", Summary = TRUE)
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
    Summary = FALSE,
    SaveAs_Table = "csv",
    FolderPath = NULL) {

    MetaProViz_Init()
    
    ## ------------------  Check Input ------------------- ##
    ## HelperFunction `CheckInput`
    CheckInput(InputData = InputData, InputData_Num = FALSE, 
        SaveAs_Table = SaveAs_Table)

    # Specific checks:
    if (!From %in% colnames(InputData)) {
        message <- paste0(From, " column was not found in InputData. Please check your input.")
        logger::log_trace(paste0("Error ", message))
        stop(message)
    }

    if (!To %in% colnames(InputData)) {
        message <- paste0(To, 
            " column was not found in InputData. Please check your input.")
        logger::log_trace(paste0("Error ", message))
        stop(message)
    }

    if (!is.null(GroupingVariable)) {
        if (!GroupingVariable %in% colnames(InputData)) {
            message <- paste0(GroupingVariable, 
                " column was not found in InputData. Please check your input.")
            logger::log_trace(paste0("Error ", message))
            stop(message)
        }
    }

    if (!is.logical(Summary)) {
        message <- paste0(
            "Check input. The Summary parameter should be either TRUE or FALSE.")
        logger::log_trace(paste0("Error ", message))
        stop(message)
    }

    ## ------------------  General checks of wrong occurences ------------------- ##
    ## Task 1: Check that From has no duplications within one group --> should 
    ## not be the case --> remove duplications and inform the user/ ask if 
    ## they forget to set groupings column
    ## Task 2: Check that From has the same items in to across the different 
    ## entries (would be in different Groupings, otherwise there should not be 
    ## any duplications) --> List of Miss-Mappings across terms

    ## FYI: The above can not happen if our translateID function was used, 
    ## but may be the case when the user has done something manually before


    ## ------------------  Create output folders and path ------------------- ##
    if (!is.null(SaveAs_Table)) {
        Folder <- SavePath(FolderName = "PriorKnowledge", 
            FolderPath = FolderPath)

        SubFolder <- file.path(Folder, "MappingAmbiguities")
        if (!dir.exists(SubFolder)) {
            dir.create(SubFolder)
        }
    }

    ############################################################################
    ## ------------------  Prepare Input data ------------------- ##
    ## if the user provides a DF where the To column is a list of IDs, then we 
    ## can use it right away
    ## if the To column is not a list of IDs, but a character column, we need 
    ## to convert it into a list of IDs
    if (is.character(InputData[[To]])) {
        InputData[[To]] <- InputData[[To]] %>%
            strsplit(", ") %>%
            lapply(as.character)
    }

    ## ------------------  Perform ambiguity mapping ------------------- ##
    ## 1. From-to-To: OriginalID-to-TranslatedID
    ## 2. From-to-To: TranslatedID-to-OriginalID
    Comp <- list(
        list(From = From, To = To),
        list(From = To, To = From)
    )

    ResList <- list()
    for (comp in seq_along(Comp)) {
        
        ## run Omnipath ambiguity
        ResList[[paste0(Comp[[comp]]$From, "-to-", Comp[[comp]]$To)]] <- InputData %>%
            ## unlist the columns in case they are not expaned
            tidyr::unnest(cols = all_of(Comp[[comp]]$From)) %>% 
            ## Remove NA values, otherwise they are counted as column is character
            filter(!is.na(!!sym(Comp[[comp]]$From))) %>% 
            OmnipathR::ambiguity(
                from_col = !!sym(Comp[[comp]]$From),
                to_col = !!sym(Comp[[comp]]$To),
                groups = GroupingVariable,
                quantify = TRUE,
                qualify = TRUE,
                ## across groups will be done additionally --> suffix _AcrossGroup
                global = TRUE,
                ## summary of the mapping column
                summary = TRUE, 
                expand = TRUE)

        ## extract summary table:
        ResList[[paste0(Comp[[comp]]$From, "-to-", Comp[[comp]]$To, "_Summary")]] <-
            ResList[[paste0(Comp[[comp]]$From, "-to-", Comp[[comp]]$To)]] %>% ## EDIT: the index could be assigned to an object since it is used several times and then be reused
            attr(paste0("ambiguity_", Comp[[comp]]$From , "_", Comp[[comp]]$To))

        ########################################################################
        if (Summary) {
        if (!is.null(GroupingVariable)) {
            ## add further information we need to summarise the table and 
            ## combine Original-to-Translated and Translated-to-Original
            ## If we have a GroupingVariable we need to combine it 
            ## with the MetaboliteID before merging
            ResList[[paste0(Comp[[comp]]$From, "-to-", Comp[[comp]]$To, "_Long")]] <- ResList[[paste0(Comp[[comp]]$From, "-to-", Comp[[comp]]$To)]] %>%
                tidyr::unnest(cols = all_of(Comp[[comp]]$From)) %>%
                mutate(!!sym(paste0("AcrossGroupMappingIssue(", Comp[[comp]]$From, "_to_", Comp[[comp]]$To, ")")) := case_when(
                    !!sym(paste0(Comp[[comp]]$From, "_", Comp[[comp]]$To, "_ambiguity_bygroup")) != !!sym(paste0(Comp[[comp]]$From, "_", Comp[[comp]]$To, "_ambiguity"))  ~ "TRUE",
                TRUE ~ "FALSE" )) %>%
                group_by(!!sym(Comp[[comp]]$From), !!sym(GroupingVariable)) %>%
                mutate(!!sym(Comp[[comp]]$To) := ifelse(
                    !!sym(Comp[[comp]]$From) == 0, NA,  ## Or another placeholder
                    paste(unique(!!sym(Comp[[comp]]$To)), collapse = ", "))) %>%
                mutate(!!sym(paste0("Count(", Comp[[comp]]$From, "_to_", Comp[[comp]]$To, ")")) := ifelse(all(!!sym(Comp[[comp]]$To) == 0), 0, n())) %>%
                ungroup()%>%
                distinct() %>%
                unite(!!sym(paste0(Comp[[comp]]$From, "_to_", Comp[[comp]]$To)), c(Comp[[comp]]$From, Comp[[comp]]$To), sep = " --> ", remove = FALSE) %>%
                separate_rows(!!sym(Comp[[comp]]$To), sep = ", ") %>%
                unite(UniqueID, c(From, To, GroupingVariable), sep = "_", remove = FALSE)%>%
                distinct() ## EDIT: looks quite complicated, could it be simplified or comments be added?
            } else{
                ResList[[paste0(Comp[[comp]]$From, "-to-", Comp[[comp]]$To, "_Long")]] <- ResList[[paste0(Comp[[comp]]$From, "-to-", Comp[[comp]]$To)]] %>%
                    tidyr::unnest(cols = all_of(Comp[[comp]]$From))%>%
                    group_by(!!sym(Comp[[comp]]$From))%>%
                    mutate(!!sym(Comp[[comp]]$To) := ifelse(!!sym(Comp[[comp]]$From) == 0, NA,  # Or another placeholder
                        paste(unique(!!sym(Comp[[comp]]$To)), collapse = ", "))) %>%
                    mutate(!!sym(paste0("Count(", Comp[[comp]]$From, "_to_", Comp[[comp]]$To, ")")) := ifelse(all(!!sym(Comp[[comp]]$To) == 0), 0, n())) %>%
                    ungroup()%>%
                    distinct() %>%
                    unite(!!sym(paste0(Comp[[comp]]$From, "_to_", Comp[[comp]]$To)), c(Comp[[comp]]$From, Comp[[comp]]$To), sep=" --> ", remove = FALSE) %>%
                    separate_rows(!!sym(Comp[[comp]]$To), sep = ", ") %>%
                    unite(UniqueID, c(From, To), sep="_", remove = FALSE) %>%
                    distinct() %>%
                    mutate(!!sym(paste0("AcrossGroupMappingIssue(", From, "_to_", To, ")")) := NA) ## EDIT: looks quite complicated, could it be simplified or comments be added?
            }
        }

        ## add NA metabolite maps back if they do exist:
        Removed <- InputData %>%
            ## unlist the columns in case they are not expaned
            tidyr::unnest(cols = all_of(Comp[[comp]]$From)) %>% 
                filter(is.na(!!sym(Comp[[comp]]$From)))
        
            if (nrow(Removed)>0) {
            ResList[[paste0(Comp[[comp]]$From, "-to-", Comp[[comp]]$To)]] <- bind_rows(
                ResList[[paste0(Comp[[comp]]$From, "-to-", Comp[[comp]]$To)]],
                test <- Removed %>%
                    bind_cols(setNames(as.list(rep(NA, length(setdiff(names(ResList[[paste0(Comp[[comp]]$From, "-to-", Comp[[comp]]$To)]]), names(Removed))))),
                        setdiff(names(ResList[[paste0(Comp[[comp]]$From, "-to-", Comp[[comp]]$To)]]), names(Removed))))  ## EDIT: looks quite complicated, could it be simplified or comments be added?
                )
            }
        }

        ## ------------------ Create SummaryTable ------------------- ##
        if (Summary) {
        ## combine the two tables
        Summary <- merge(
            x = ResList[[paste0(From, "-to-", To, "_Long")]][, c("UniqueID", paste0(From, "_to_", To), paste0("Count(", From, "_to_", To, ")"), paste0("AcrossGroupMappingIssue(", From, "_to_", To, ")"))],
            y= ResList[[paste0(To, "-to-", From, "_Long")]][, c("UniqueID", paste0(To, "_to_", From), paste0("Count(", To, "_to_", From, ")"), paste0("AcrossGroupMappingIssue(", To, "_to_", From, ")"))], ## EDIT: looks quite complicated, could it be simplified or comments be added?
            by = "UniqueID",
            all = TRUE) %>%
            separate(UniqueID, into = c(From, To, GroupingVariable), sep = "_", remove = FALSE) %>%
            distinct()

        ## Add relevant mapping information
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
            mutate(!!sym(paste0("Count(", From, "_to_", To, ")")) := replace_na(!!sym(paste0("Count(", From, "_to_", To, ")")), 0)) %>%
            mutate(!!sym(paste0("Count(", To, "_to_", From, ")")) := replace_na(!!sym(paste0("Count(", To, "_to_", From, ")")), 0))

            ResList[["Summary"]] <- Summary
    }

    ## ------------------ Save the results ------------------- ##
    suppressMessages(suppressWarnings(
        SaveRes(InputList_DF = ResList,
            InputList_Plot = NULL,
            SaveAs_Table = SaveAs_Table,
            SaveAs_Plot = NULL,
            FolderPath = SubFolder,
            FileName = "MappingAmbiguity",
            CoRe = FALSE,
            PrintPlot = FALSE)))

    ## return
    invisible(ResList)
}

##########################################################################################
### ### ### Check Measured ID's in prior knowledge ### ### ###
##########################################################################################

#' Check and summarize PriorKnowledge-to-MeasuredFeatures relationship
#'
#' @param InputData Dataframe with at least one column with the detected metabolite IDs (e.g. HMDB). If there are multiple IDs per detected peak, please separate them by comma ("," or ", " or chr list). If there is a main ID and additional IDs, please provide them in separate columns.
#' @param PriorKnowledge Dataframe with at least one column with the metabolite ID (e.g. HMDB) that need to match InputData metabolite IDs "source" (e.g. term). If there are multiple IDs, as the original pathway IDs (e.g. KEGG) where translated (e.g. to HMDB), please separate them by comma ("," or ", " or chr list).
#' @param SettingsInfo Colum name of Metabolite IDs in InputData and PriorKnowledge as well as column name of GroupingVariable in PriorKnowledge. \strong{Default = c(InputID="HMDB", PriorID="HMDB", GroupingVariable="term")}
#'
#' @importFrom dplyr mutate
#' @importFrom rlang !!! !! := sym syms
#'
#' @examples
#' DetectedIDs <-  ToyData(Data="Cells_MetaData") %>% 
#'     rownames_to_column("Metabolite") %>% 
#'     dplyr::select("Metabolite", "HMDB") %>%
#'     tidyr::drop_na()
#' PathwayFile <- TranslateID(InputData = LoadKEGG(), 
#'         SettingsInfo = c(InputID="MetaboliteID", GroupingVariable="term"), 
#'         From = c("kegg"), To = c("hmdb"))[["TranslatedDF"]] %>%
#'     tidyr::drop_na()
#' Res <- CleanMapping(InputData= DetectedIDs, PriorKnowledge = PathwayFile, 
#'     SettingsInfo = c(InputID = "HMDB", PriorID = "hmdb", 
#'         GroupingVariable = "term"))
#'
#' @noRd
#'
CheckMatchID <- function(InputData,
    PriorKnowledge,
    SettingsInfo = c(InputID = "HMDB", PriorID ="HMDB", GroupingVariable = "term")) {

    ## ------------ Create log file ----------- ##
    MetaProViz_Init()

    ## ------------ Check Input files ----------- ##

    ## InputData:
    if ("InputID" %in% names(SettingsInfo)) {
        if (!SettingsInfo[["InputID"]] %in% colnames(InputData)) {
            message <- paste0("The ", SettingsInfo[["InputID"]], " column selected as InpuID in SettingsInfo was not found in InputData. Please check your input.")
            logger::log_trace(paste("Error ", message, sep=""))
            stop(message)
        }
    } else{
        message <- paste0("No ", SettingsInfo[["InputID"]], " provided. Please check your input.")
        logger::log_trace(paste("Error ", message, sep=""))
        stop(message)
    }

    if (sum(is.na(InputData[[SettingsInfo[["InputID"]]]])) >= 1) {
        ## remove NAs:
        message <- paste0(sum(is.na(InputData[[SettingsInfo[["InputID"]]]])), 
            " NA values were removed from column", SettingsInfo[["InputID"]])
        logger::log_trace(paste0("Warning: ", message))

        InputData <- InputData %>%
            filter(!is.na(.data[[SettingsInfo[["InputID"]]]]))
        warning(message)
    }

    if (nrow(InputData) - nrow(distinct(InputData, .data[[SettingsInfo[["InputID"]]]])) >= 1) {
        ## Remove duplicate IDs
        message <- paste0(
            nrow(InputData) -  nrow(distinct(InputData, .data[[SettingsInfo[["InputID"]]]])), 
            " duplicated IDs were removed from column", SettingsInfo[["InputID"]])
        logger::log_trace(paste("Warning: ", message, sep=""))

        InputData <- InputData %>%
            distinct(.data[[SettingsInfo[["InputID"]]]], .keep_all = TRUE)

        warning(message)
    }

    InputData_MultipleIDs <- any(
        grepl(",\\s*", InputData[[SettingsInfo[["InputID"]]]]) |  # Comma-separated
            sapply(InputData[[SettingsInfo[["InputID"]]]], function(x) {
            if (grepl("^c\\(|^list\\(", x)) {
                parsed <- tryCatch(eval(parse(text = x)), error = function(e) NULL)
                return(is.list(parsed) && length(parsed) > 1 || is.vector(parsed) && length(parsed) > 1)
            }
            FALSE
        })
    )

    ## PriorKnowledge:
    if ("PriorID" %in% names(SettingsInfo)) {
        if (!SettingsInfo[["PriorID"]] %in% colnames(PriorKnowledge)) {
            message <- paste0("The ", SettingsInfo[["PriorID"]], 
                " column selected as InpuID in SettingsInfo was not found in PriorKnowledge. Please check your input.")
            logger::log_trace(paste0("Error ", message))
            stop(message)
        }
    } else {
        message <- paste0("No ", SettingsInfo[["PriorID"]], " provided. Please check your input.")
        logger::log_trace(paste0("Error ", message))
        stop(message)
    }

    if (sum(is.na(PriorKnowledge[[SettingsInfo[["PriorID"]]]])) >= 1) {
        ## remove NAs:
        message <- paste0(sum(is.na(PriorKnowledge[[SettingsInfo[["PriorID"]]]])), 
            " NA values were removed from column", SettingsInfo[["PriorID"]])
        logger::log_trace(paste0("Warning: ", message))

        PriorKnowledge <- PriorKnowledge %>%
            filter(!is.na(.data[[SettingsInfo[["PriorID"]]]]))
        warning(message)
    }

    if ("GroupingVariable" %in% names(SettingsInfo)) {
        ## add GroupingVariable
        if (!SettingsInfo[["GroupingVariable"]] %in% colnames(PriorKnowledge)) {
            message <- paste0("The ", SettingsInfo[["GroupingVariable"]], 
                " column selected as InpuID in SettingsInfo was not found in ",
                "PriorKnowledge. Please check your input.")
            logger::log_trace(paste0("Error ", message))
            stop(message)
        }
    } else{
        ## add GroupingVariable
        SettingsInfo["GroupingVariable"] <- "GroupingVariable"
        PriorKnowledge["GroupingVariable"] <- "None"

        message <- paste0("No ", SettingsInfo[["PriorID"]], 
            " provided. If this was not intentional, please check your input.")
        logger::log_trace(message)
        message(message)
    }

    if (nrow(PriorKnowledge) - nrow(distinct(PriorKnowledge, .data[[SettingsInfo[["PriorID"]]]], .data[[SettingsInfo[["GroupingVariable"]]]])) >= 1) {
        # Remove duplicate IDs
        message <- paste0(
            nrow(PriorKnowledge) - nrow(distinct(PriorKnowledge, .data[[SettingsInfo[["PriorID"]]]], .data[[SettingsInfo[["GroupingVariable"]]]])), 
            " duplicated IDs were removed from column", SettingsInfo[["PriorID"]])
        logger::log_trace(paste("Warning: ", message, sep=""))

        PriorKnowledge <- PriorKnowledge %>%
            distinct(.data[[SettingsInfo[["PriorID"]]]], !!sym(SettingsInfo[["GroupingVariable"]]), .keep_all = TRUE) %>%
            group_by(!!sym(SettingsInfo[["PriorID"]])) %>%
            mutate(across(everything(), ~ if (is.character(.)) paste(unique(.), collapse = ", "))) %>%
            ungroup() %>%
            distinct(.data[[SettingsInfo[["PriorID"]]]], .keep_all = TRUE)

        warning(message)
    }

    PK_MultipleIDs <- any(
        ## check if multiple IDs are present:
        grepl(",\\s*", PriorKnowledge[[SettingsInfo[["PriorID"]]]]) |  # Comma-separated
            sapply(PriorKnowledge[[SettingsInfo[["PriorID"]]]] , 
            function(x) {
                if (grepl("^c\\(|^list\\(", x)) {
                    parsed <- tryCatch(eval(parse(text = x)), error = function(e) NULL)
                    return(is.list(parsed) && length(parsed) > 1 || is.vector(parsed) && length(parsed) > 1)
                }
                FALSE
            }) ## EDIT: is this the same as above, if so, it should be simplified (functionised?)
    )

    ## ------------ Create Results output folder ----------- ##
    if (!is.null(SaveAs_Table)) {
        Folder <- SavePath(FolderName = "PriorKnowledgeChecks",
            FolderPath = FolderPath)
        SubFolder <- file.path(Folder, "CheckMatchID_Detected-to-PK")
        if (!dir.exists(SubFolder)) {
            dir.create(SubFolder)
        }
    }

    ############################################################################
    ## ------------ Check how IDs match and if needed remove unmatched IDs ----------- ##

    # 1.  Create long DF
    create_long_df <- function(df, id_col, df_name) { ## EDIT: function should be not defined within other functions
        df %>%
            mutate(row_id = dplyr::row_number()) %>%
            mutate(!!paste0("OriginalEntry_", df_name, sep="") := !!sym(id_col)) %>%  
            ## Store original values
            separate_rows(!!sym(id_col), sep = ",\\s*") %>%
            group_by(row_id) %>%
            mutate(!!(paste0("OriginalGroup_", df_name)) := paste0(df_name, "_", dplyr::cur_group_id())) %>%
            ungroup()
    }

    if (InputData_MultipleIDs) {
        InputData_long <- create_long_df(InputData, SettingsInfo[["InputID"]], "InputData") %>%
            select(SettingsInfo[["InputID"]], "OriginalEntry_InputData", OriginalGroup_InputData)
    } else{
        InputData_long <- InputData %>%
            mutate(OriginalGroup_InputData := paste0("InputData_", dplyr::row_number()))%>%
            select(SettingsInfo[["InputID"]], OriginalGroup_InputData)
    }

    if (PK_MultipleIDs) {
        PK_long <- create_long_df(PriorKnowledge, SettingsInfo[["PriorID"]], "PK") %>%
            select(SettingsInfo[["PriorID"]], "OriginalEntry_PK", OriginalGroup_PK, SettingsInfo[["GroupingVariable"]])
    } else{
        PK_long <- PriorKnowledge %>%
            mutate(OriginalGroup_PK := paste0("PK_", dplyr::row_number()))%>%
            select(SettingsInfo[["PriorID"]],OriginalGroup_PK, SettingsInfo[["GroupingVariable"]])
    }

    ## 2. merge DF
    merged_df <- merge(PK_long, InputData_long, 
            by.x = SettingsInfo[["PriorID"]],  by.y = SettingsInfo[["InputID"]], 
            all = TRUE) %>%
        distinct(!!sym(SettingsInfo[["PriorID"]]), OriginalGroup_InputData, 
            .keep_all = TRUE)

    ## 3. add information to summarize and describe problems
    merged_df <- merged_df %>%
        ## num_PK_entries
        group_by(OriginalGroup_PK, !!sym(SettingsInfo[["GroupingVariable"]])) %>%
        mutate(
            num_PK_entries = sum(!is.na(OriginalGroup_PK)),
            num_PK_entries_groups = dplyr::n_distinct(OriginalGroup_PK, na.rm = TRUE)) %>% # count the times we have the same PK_entry match with multiple InputData entries --> extend below!
        ungroup()%>%
        # num_Input_entries
        group_by(OriginalGroup_InputData, !!sym(SettingsInfo[["GroupingVariable"]])) %>%
        mutate(
            num_Input_entries = sum(!is.na(OriginalGroup_InputData)),
            num_Input_entries_groups = dplyr::n_distinct(OriginalGroup_InputData,, na.rm = TRUE)) %>%
        ungroup() %>%
        mutate(
            ActionRequired = case_when(
                num_Input_entries == 1 & num_Input_entries_groups == 1 & num_PK_entries_groups == 1 ~ "None",
                num_Input_entries == 1 & num_Input_entries_groups == 1 & num_PK_entries_groups >= 2  ~ "Check",
                num_Input_entries == 1 & num_Input_entries_groups >= 2 & num_PK_entries_groups == 1  ~ "Check",
                num_Input_entries > 1 & num_Input_entries_groups == 1 & num_PK_entries_groups >= 2 ~ "Check",
                num_Input_entries > 1 & num_Input_entries_groups >= 2 & num_PK_entries_groups >= 2 ~ "Check",
                num_Input_entries == 0 ~ "None", # ID(s) of PK not measured
                TRUE ~ NA_character_
            )) %>%
        mutate(
            Detection = case_when(
                num_Input_entries == 1 & num_Input_entries_groups == 1 & num_PK_entries_groups == 1 ~ "One input ID of the same group maps to at least ONE PK ID of ONE group",
                num_Input_entries == 1 & num_Input_entries_groups == 1 & num_PK_entries_groups >= 2  ~ "One input ID of the same group maps to at least ONE PK ID of MANY groups",
                num_Input_entries == 1 & num_Input_entries_groups >= 2 & num_PK_entries_groups == 1  ~ "One input ID of MANY groups maps to at least ONE PK ID of ONE groups",
                num_Input_entries > 1 & num_Input_entries_groups == 1 & num_PK_entries_groups >= 2 ~ "MANY input IDs of the same group map to at least ONE PK ID of MANY PK groups",
                num_Input_entries > 1 & num_Input_entries_groups >= 2 & num_PK_entries_groups >= 2 ~ "MANY input IDs of the MANY groups map to at least ONE PK ID of MANY PK groups",
                num_Input_entries == 0 ~ "Not Detected", # ID(s) of PK not measured
                TRUE ~ NA_character_))
    
    ##  Handle "Detected-to-PK" (When PK has multiple IDs)
    #group_by(OriginalGroup_PK, !!sym(SettingsInfo[["GroupingVariable"]])) %>%
    #mutate(
    #  `Detected-to-PK` = case_when(
    #    num_Input_entries == 0 & num_PK_entries == 1 ~ "none-to-one", # No match & OriginalGroup_PK appears once
    #    num_Input_entries == 0 & num_PK_entries > 1 ~ "none-to-many", # No match & OriginalGroup_PK appears multiple times
    #    num_Input_entries == 1 & num_PK_entries == 1 ~ "one-to-one", # One unique match in InputData, one PK
    #    num_Input_entries == 1 & num_PK_entries > 1 ~ "one-to-many", # One unique match in InputData, multiple PKs
    #    num_Input_entries > 1 & num_PK_entries == 1 ~ "many-to-one", # Multiple matches in InputData, one PK
    #    num_Input_entries > 1 & num_PK_entries > 1 ~ "many-to-many", # Multiple matches in InputData, multiple PKs
    #    num_Input_entries == 1 & num_PK_entries == 0 ~ "one-to-none",
    #    num_Input_entries > 1 & num_PK_entries == 0 ~ "many-to-none",
    #    TRUE ~ NA_character_
    #  )
    #) %>%
    #ungroup() %>%
    # Handle "PK-to-Detected" (When InputData has multiple IDs)
    #group_by(OriginalGroup_InputData, !!sym(SettingsInfo[["GroupingVariable"]])) %>%
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

    ## 4. Create summary table
    Values_InputData <- unique(InputData[[SettingsInfo[["InputID"]]]])
    Values_PK <- unique(PK_long[[SettingsInfo[["PriorID"]]]])

    summary_df <- tibble::tibble(
        !!sym(SettingsInfo[["InputID"]]) := Values_InputData,
        found_match_in_PK = NA,
        matches = NA_character_,
        match_overlap_percentage = NA_real_,
        original_count = NA_integer_,
        matches_count = NA_integer_
    )

    ## populate the summary data frame
    for (i in seq_along(Values_InputData)) {
        # Handle NA case explicitly
        if (is.na(Values_InputData[i])) {
            summary_df$original_count[i] <- 0
            summary_df$matches_count[i] <- 0
            summary_df$match_overlap_percentage[i] <- NA
            ## could also set it to FALSE but making NA for now for plotting
            summary_df$found_match_in_PK[i] <- NA 
            summary_df$matches[i] <- NA
        } else {
            ## split each cell into individual entries and trim whitespace
            entries <- trimws(
                unlist(strsplit(as.character(Values_InputData[i]), ",\\s*"))) ## delimiter = "," or ", "

            ## identify which entries are in the lookup set
            matched <- entries[entries %in% Values_PK]

            ## determine if any match was found
            summary_df$found_match_in_PK[i] <- length(matched) > 0

            ## concatenate matched entries into a single string
            summary_df$matches[i] <- paste(matched, collapse = ", ")

            ## calculate and store counts
            summary_df$original_count[i] <- length(entries)
            summary_df$matches_count[i] <- length(matched)

            ## calculate fraction: matched entries / total entries
            if (length(entries) > 0) {
                summary_df$match_overlap_percentage[i] <- (length(matched) / length(entries))*100
            } else {
                summary_df$match_overlap_percentage[i] <- NA
            }
        }
    }

    summary_df <- merge(
        x = summary_df,
        y = merged_df %>%
            dplyr::select(-c(OriginalGroup_PK, OriginalGroup_InputData)) %>%
            distinct(!!sym(SettingsInfo[["PriorID"]]), .keep_all = TRUE),
        by.x = SettingsInfo[["InputID"]] ,
        by.y = SettingsInfo[["PriorID"]],
        all.x = TRUE)

    summary_df <- merge(x = summary_df, y = InputData, 
            by = SettingsInfo[["InputID"]], all.x = TRUE) %>%
        distinct(!!sym(SettingsInfo[["InputID"]]), OriginalEntry_PK, 
            .keep_all = TRUE)

    ## 5. Messages and summarise
    message <- paste0("InputData has multiple IDs per measurement = ", 
        InputData_MultipleIDs, ". PriorKnowledge has multiple IDs per entry = ", 
        PK_MultipleIDs, ".")
    message1 <- paste0("InputData has ", 
        dplyr::n_distinct(unique(InputData[[SettingsInfo[["InputID"]]]])), 
        " unique entries with ", 
        dplyr::n_distinct(unique(InputData_long[[SettingsInfo[["InputID"]]]])), 
        " unique ", SettingsInfo[["InputID"]], " IDs. Of those IDs, ", 
        nrow(dplyr::filter(summary_df, matches_count == 1)), 
        " match, which is ", 
        nrow(dplyr::filter(summary_df, matches_count == 1)) / dplyr::n_distinct(unique(InputData_long[[SettingsInfo[["InputID"]]]])) * 100, "%." , sep = "")
    message2 <- paste0("PriorKnowledge has ", 
        dplyr::n_distinct(PriorKnowledge[[SettingsInfo[["PriorID"]]]]), 
        " unique entries with ", 
        dplyr::n_distinct(PK_long[[SettingsInfo[["PriorID"]]]]), 
        " unique ", SettingsInfo[["PriorID"]], " IDs. Of those IDs, ", 
        nrow(dplyr::filter(summary_df, matches_count == 1)), 
        " are detected in the data, which is ", 
        nrow(dplyr::filter(summary_df, matches_count == 1)) / dplyr::n_distinct(PK_long[[SettingsInfo[["PriorID"]]]])) * 100, "%.") ## EDIT: consider precomputing repeated objects

    if (nrow(dplyr::filter(summary_df, ActionRequired == "Check")) >= 1) {
        #warning <- paste0("There are cases where multiple detected IDs match to multiple prior knowledge IDs of the same category") # "Check"

    }

    ## ------------------ Plot Summary ----------------------##
    # x = "Class" and y = Frequency. Match Status can be colour of if no class provided class = Match status.
    # Check Biocrates code.

    ## ------------------ Save Results ----------------------##
    ResList <- list("InputData_Matched" = summary_df)

    suppressMessages(
        suppressWarnings(
            SaveRes(InputList_DF = ResList,
                InputList_Plot = NULL,
                SaveAs_Table = SaveAs_Table,
                SaveAs_Plot = NULL,
                FolderPath = SubFolder,
                FileName = "CheckMatchID_Detected-to-PK",
                CoRe = FALSE,
                PrintPlot = FALSE)))

   ## return
   invisible(ResList[["InputData_Matched"]])
}


##########################################################################################
### ### ### Cluster Prior Knowledge ### ### ###
##########################################################################################

#' Deal with pathway overlap in prior knowledge
#'
#' @param InputData Dataframe with at least one column with the target (e.g. metabolite) and a column source (e.g. term).
#' @param SettingsInfo = c(InputID="MetaboliteID", GroupingVariable="term"),
#'
#' @examples
#' KEGG_Pathways <- LoadKEGG()
#' InputData = KEGG_Pathways
#'
#'
#' @noRd
ClusterPK <- function(
    InputData, # This can be either the original PK (e.g. KEGG pathways), but it can also be the output of enrichment results (--> meaning here we would cluster based on detection!)
    SettingsInfo = c(InputID = "MetaboliteID", GroupingVariable = "term"),
    Clust = c("Graph", "Hierarchical"), # Options: "Graph", "Hierarchical",
    matrix = c("percentage", "pearson", "spearman", "kendall"), # Choose "pearson", "spearman", "kendall", or "percentage"
    min = 2 # minimum pathways per cluster) {
    
    # Cluster PK before running enrichment analysis --> add another column that groups the data based on the pathway overlap:
    # provide different options for clustering (e.g. % of overlap, semantics similarity) --> Ramp uses % of overlap, semnatics similarity: https://yulab-smu.top/biomedical-knowledge-mining-book/GOSemSim.html


    ## ------------------ Check Input ------------------- ##
    ## match arguments
    Clust <- match.arg(Clust)
    matrix <- match.arg(matrix)

    ## ------------------ Create output folders and path ------------------- ##



    ############################################################################
    ## ------------------ Cluster the data ------------------- ##
    ## 1. create a list of unique MetaboliteIDs for each term
    term_metabolites <- InputData %>%
        dplyr::group_by(!!sym(SettingsInfo[["GroupingVariable"]])) %>%
        dplyr::summarize(MetaboliteIDs = list(unique(!!sym(SettingsInfo[["InputID"]])))) %>%
        dplyr::ungroup()

    ## 2. Create the overlap matrix based on different methods:
    if (matrix == "percentage") {
        # Compute pairwise overlaps
        term_overlap <- combn(
            term_metabolites[[SettingsInfo[["GroupingVariable"]]]], 2, 
            function(terms) {
                term1_ids <- term_metabolites$MetaboliteIDs[term_metabolites[[SettingsInfo[["GroupingVariable"]]]] == terms[1]][[1]]
                term2_ids <- term_metabolites$MetaboliteIDs[term_metabolites[[SettingsInfo[["GroupingVariable"]]]] == terms[2]][[1]]

                overlap <- length(intersect(term1_ids, term2_ids)) / length(union(term1_ids, term2_ids))
                data.frame(Term1 = terms[1], Term2 = terms[2], Overlap = overlap)
            }, simplify = FALSE) %>%
            dplyr::bind_rows()

        ## create overlap matrix: An overlap matrix is typically used to 
        ## quantify the degree of overlap between two sets or groups,
        ## overlap coefficient (or Jaccard Index) Overlap(A,B)= ∣A∪B∣ / ∣A∩B
        ## the overlap matrix measures the similarity between sets or groups based on common elements.
        terms <- unique(c(term_overlap$Term1, term_overlap$Term2))
        overlap_matrix <- matrix(1, nrow = length(terms), 
            ncol = length(terms), dimnames = list(terms, terms))
        for (i in seq_len(nrow(term_overlap))) {
            t1 <- term_overlap$Term1[i]
            t2 <- term_overlap$Term2[i]
            overlap_matrix[t1, t2] <- 1 - term_overlap$Overlap[i]
            overlap_matrix[t2, t1] <- 1 - term_overlap$Overlap[i]
        }
    } else {
        ## create a binary matrix for correlation methods
        terms <- term_metabolites[[SettingsInfo[["GroupingVariable"]]]]
        metabolites <- unique(unlist(term_metabolites$MetaboliteIDs)) #[[SettingsInfo[["InputID"]]]]

        binary_matrix <- matrix(0, nrow = length(terms), 
            ncol = length(metabolites), dimnames = list(terms, metabolites))
        
        for (i in seq_along(terms)) {
            metabolites_for_term <- term_metabolites$MetaboliteIDs[[i]] #[[SettingsInfo[["InputID"]]]]
            binary_matrix[i, colnames(binary_matrix) %in% metabolites_for_term] <- 1
        }

        ## Compute correlation matrix: square matrix used to represent the 
        ## pairwise correlation coefficients between variables or terms
        ## correlation matrix 𝐶 C is an 𝑛 × 𝑛 n×n matrix where each 
        ## element 𝐶 𝑖 𝑗 C ij  is the correlation coefficient between the
        ## variables 𝑋𝑖 Xi  and 𝑋 𝑗 X j
        ## the correlation matrix measures the strength and direction of 
        ## linear relationships between variables.
        correlation_matrix <- cor(t(binary_matrix), method = matrix)

        ## Convert to distance matrix
        overlap_matrix <- 1 - correlation_matrix
    }
  
    ## 3. Cluster terms based on overlap threshold
    ## define similarity threshold
    threshold <- 0.7
    term_clusters <- term_overlap %>%
        dplyr::filter(Overlap >= threshold) %>%
        dplyr::select(Term1, Term2)

    ## 4. Clustering
    if (Clust == "Graph") { 
        ## use Graph-based clustering
        ## calculate the distance matrix:
        overlap_matrix <- 1 - correlation_matrix

        ## an adjacency matrix represents a graph structure and encodes the 
        ## relationships between nodes (vertices)
        ## add weight (can also represent unweighted graphs)
        ## Applying Gaussian kernel to convert distance into similarity
        adjacency_matrix <- exp(-overlap_matrix ^ 2)  

        ## create a graph from the adjacency matrix
        g <- igraph::graph_from_adjacency_matrix(adjacency_matrix, 
            mode = "undirected", weighted = TRUE)
        initial_clusters <- igraph::components(g)$membership
        term_metabolites$Cluster <- initial_clusters[
            match(term_metabolites[[SettingsInfo[["GroupingVariable"]]]], 
                names(initial_clusters))]
    } else if (Clust == "Hierarchical") {
        ## hierarchical clustering
        ## make methods into parameters!
        hclust_result <- hclust(as.dist(distance_matrix), method = "average") 
        num_clusters <- 4
        term_clusters_hclust <- cutree(hclust_result, k = num_clusters)

        term_metabolites$Cluster <- paste0("Cluster", 
            term_clusters_hclust[match(terms, names(term_clusters_hclust))])
        #term_metabolites$Cluster <- clusters[match(term_metabolites[[SettingsInfo[["GroupingVariable"]]]], names(clusters))]
    } else { ## EDIT: not needed with match.arg
        stop("Invalid clustering method specified in Clust parameter.")
    }

    ## 5. Merge cluster group information back to the original data
    df <- InputData %>%
        dplyr::left_join(
            select(term_metabolites, !!sym(SettingsInfo[["GroupingVariable"]]), Cluster), 
            by = SettingsInfo[["GroupingVariable"]]) %>%
        dplyr::mutate(Cluster = ifelse(
            is.na(Cluster),
            ## assign "None" to NAs
            "None", 
            ## convert numeric IDs to descriptive labels
            paste0("Cluster", Cluster)))

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
    complete = FALSE) {

    ## ------------------ check input ------------------- ##


    ## ------------------ create output folders and path ------------------- ##


    ## ------------------ add information to enrichment results ------------------- ##

    ## add number of Genes_targeted_by_TF_num
    net$Count <- 1
    net_Mean <- aggregate(net$Count, 
            by = list(source = net[[.source]]), FUN = sum) %>%
        rename("targets_num" = 2)

    if (complete) {
        res_Add <- merge(x = res, y = net_Mean, by = "source", all = TRUE)
    } else{
        res_Add <- merge(x = res, y = net_Mean, by = "source", all.x = TRUE)
    } ## EDIT: would a left_join do the same without the if/else? 

    ## add list of Genes_targeted_by_TF_chr
    net_List <- aggregate(net[[.target]] ~ net[[.source]], FUN = toString) %>%
        rename("source" = 1, "targets_chr" = 2)
    res_Add <- merge(x = res_Add, y = net_List, by = "source", all.x = TRUE)

    ## add number of Genes_targeted_by_TF_detected_num
    mat <- as.data.frame(mat) %>% ## Are these the normalised counts?
        tibble::rownames_to_column("Symbol")

    Detected <- merge(x = mat , y = net[,c(.source, .target)], 
            by.x = "Symbol", by.y = .target, all.x = TRUE) %>%
        filter(!is.na(across(all_of(.source))))
    Detected$Count <-1
    Detected_Mean <- aggregate(Detected$Count, 
            by = list(source = Detected[[.source]]), FUN = sum) %>%
        rename("targets_detected_num" = 2)

    res_Add <- merge(x= res_Add, y = Detected_Mean, by = "source", 
        all.x = TRUE) %>%
    mutate(targets_detected_num = replace_na(targets_detected_num, 0))

    ## add list of Genes_targeted_by_TF_detected_chr
    Detected_List <- aggregate(Detected$Symbol ~ Detected[[.source]], 
            FUN = toString) %>%
        rename("source" = 1, "targets_detected_chr" = 2)

    res_Add <- merge(x = res_Add, y = Detected_List, by = "source", 
        all.x = TRUE)

    ## add percentage of Percentage_of_Genes_detected
    res_Add$targets_detected_percentage <- round(
        res_Add$targets_detected_num / res_Add$targets_num * 100, digits = 2)

    ## sort by score
    res_Add <- res_Add %>%
        arrange(desc(as.numeric(as.character(score))))

    ## ------------------ Save and return ------------------- ##
    res_Add
}
