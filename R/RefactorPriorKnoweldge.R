#!/usr/bin/env Rscript

#
#  This file is part of the `MetaProViz` R package
#
#  Copyright 2022-2025
#  Saez Lab, Heidelberg University
#
#  Authors: see the file `README.md`
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file `LICENSE` or copy at
#      https://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: https://saezlab.github.io/MetaProViz
#  Git repo: https://github.com/saezlab/MetaProViz
#


#
# Translate IDs to/from KEGG, PubChem, Chebi, HMDB
#

#' Translate IDs to/from KEGG, PubChem, Chebi, HMDB
#'
#' @param data dataframe with at least one column with the target (e.g. metabolite),
#'     you can add other columns such as source (e.g. term). Must be "long" DF,
#'     meaning one ID per row.
#' @param metadata_info \emph{Optional: } Column name of Target in input_pk. \strong{Default =
#'     list(InputID="MetaboliteID" , grouping_variable="term")}
#' @param from ID type that is present in your data. Choose between "kegg", "pubchem",
#'     "chebi", "hmdb". \strong{Default = "kegg"}
#' @param to One or multiple ID types to which you want to translate your data.
#'     Choose between "kegg", "pubchem", "chebi", "hmdb". \strong{Default =
#'     c("pubchem","chebi","hmdb")}
#' @param summary \emph{Optional: } If TRUE a long summary tables are created.
#'     \strong{Default = FALSE}
#' @param save_table \emph{Optional: } File types for the analysis results are: "csv",
#'     "xlsx", "txt". \strong{Default = "csv"}
#' @param path {Optional:} Path to the folder the results should be saved at.
#'     \strong{Default = NULL}
#'
#' @return List with at least three DFs: 1) Original data and the new column of
#'     translated ids spearated by comma. 2) Mapping information between
#'     Original ID to Translated ID. 3) Mapping summary between Original ID to
#'     Translated ID.
#'
#' @examples
#' \dontrun{
#' KEGG_Pathways <- metsigdb_kegg()
#' Res <- translate_id(
#'     data = KEGG_Pathways,
#'     metadata_info = c(
#'         InputID = "MetaboliteID",
#'         grouping_variable = "term"
#'     ),
#'     from = c("kegg"),
#'     to = c("pubchem", "hmdb")
#' )
#' }
#'
#' @importFrom dplyr mutate select group_by ungroup distinct filter across n
#' @importFrom tidyselect all_of starts_with
#' @importFrom rlang !!! !! := sym syms
#' @importFrom OmnipathR id_types translate_ids
#' @importFrom logger log_warn log_trace
#' @importFrom stringr str_to_lower
#' @export
translate_id <- function(
    data,
    metadata_info = c(
        InputID = "MetaboliteID",
        grouping_variable = "term"
    ),
    from = "kegg",
    to = c("pubchem", "chebi", "hmdb"),
    summary = FALSE,
    save_table = "csv",
    path = NULL
) {
    # Add ability to also get metabolite names that are human readable
    # from an ID type!

    metaproviz_init()

    # # ------------------  Check Input ------------------- ##
    # HelperFunction `check_param`
    check_param(
        data = data,
        data_num = FALSE,
        save_table = save_table
    )

    # Specific checks:
    if ("InputID" %in% names(metadata_info)) {
        if (!(metadata_info[["InputID"]] %in% colnames(data))) {
            message <-
                paste0(
                    "The ",
                    metadata_info[["InputID"]],
                    " column selected as InputID in metadata_info was ",
                    "not found in data. ",
                    "Please check your input."
                )
            log_trace(paste("Error ", message, sep = ""))
            stop(message)
        }
    }

    if ("grouping_variable" %in% names(metadata_info)) {
        if (!(metadata_info[["grouping_variable"]] %in% colnames(data))) {
            message <-
                paste0(
                    "The ",
                    metadata_info[["grouping_variable"]],
                    " column selected as grouping_variable in metadata_info ",
                    "was not found in data. ",
                    "Please check your input."
                )
            log_trace(paste("Error ", message, sep = ""))
            stop(message)
        }
    }

    if (!is.logical(summary)) {
        message <-
            paste0(
                "Check input. The summary parameter should be either =TRUE or ",
                "=FALSE."
            )
        log_trace(paste("Error ", message, sep = ""))
        stop(message)
    }

    unknown_types <-
        id_types() %>%
        select(
            starts_with("in_")
        ) %>%
        unlist() %>%
        unique() %>%
        str_to_lower() %>%
        setdiff(
            union(from, to), .
        )

    if (length(unknown_types) > 0L) {
        msg <-
            sprintf(
                "The following ID types are not recognized: %s",
                paste(unknown_types, collapse = ", ")
            )
        log_warn(msg)
        warning(msg)
    }

    # Check that metadata_info[['InputID']] has no duplications within one group
    # --> should not be the case --> remove duplications and inform the user/
    # ask if they forget to set groupings column
    doublons <-
        data %>%
        filter(
            !is.na(
                !!sym(
                    metadata_info[["InputID"]]
                )
            )
        ) %>%
        group_by(
            !!sym(metadata_info[["InputID"]]),
            !!sym(metadata_info[["grouping_variable"]])
        ) %>%
        filter(
            n() > 1L
        ) %>%
        summarize()

    if (nrow(doublons) > 0L) {
        message <-
            sprintf(
                "The following IDs are duplicated within one group: %s",
                paste(
                    doublons %>% pull(metadata_info[["InputID"]]),
                    collapse = ", "
                )
            )
        log_warn(message)
        warning(message)
    }

    # # ------------------  Create output folders and path ------------------- ##
    if (!is.null(save_table)) {
        folder <-
            save_path(
                folder_name = "PK",
                path = path
            )

        Subfolder <- file.path(folder, "TranslateIDs")
        if (!dir.exists(Subfolder)) {
            dir.create(Subfolder)
        }
    }

    # ###########################################################################
    # # ------------------ Translate to-from for each pair ------------------- ##
    TranslatedDF <-
        translate_ids(
            data,
            !!sym(metadata_info[["InputID"]]) := !!sym(from),
            !!!syms(to),  # list of symbols, hence three !!!
            ramp = TRUE,
            expand = FALSE,
            quantify_ambiguity = TRUE,
            qualify_ambiguity = TRUE,
            # Checks within the groups, without it checks across groups
            ambiguity_groups = metadata_info[["grouping_variable"]],
            ambiguity_summary = TRUE
        )
    # TranslatedDF %>% attributes %>% names
    # TranslatedDF%>% attr('ambiguity_MetaboliteID_hmdb')

    # # --------------- Create output DF -------------------- ##
    ResList <- list()

    # # Create DF for TranslatedIDs only with the original data and the
    # # translatedID columns
    DF_subset <-
        TranslatedDF %>%
        select(
            all_of(intersect(names(.), names(data))), all_of(to)
        ) %>%
        mutate(
            across(
                all_of(to),
                ~ map_chr(., ~ paste(unique(.), collapse = ", "))
            )
        ) %>%
        group_by(
            !!sym(metadata_info[["InputID"]]),
            !!sym(metadata_info[["grouping_variable"]])
        ) %>%
        mutate(
            across(all_of(to),
                ~ paste(unique(.),
                    collapse = ", "
                ),
                .names = "{.col}"
            )
        ) %>%
        ungroup() %>%
        distinct() %>%
        mutate(
            across(
                all_of(to),
                ~ ifelse(. == "0", NA, .)
            )
        )

    ResList[["TranslatedDF"]] <- DF_subset

    # # Add DF with mapping information
    ResList[["TranslatedDF_MappingInfo"]] <- TranslatedDF

    # # Also save the different mapping summaries!
    for (item in to) {
        summaryDF <-
            TranslatedDF %>%
            attr(
                paste0(
                    "ambiguity_",
                    metadata_info[["InputID"]],
                    "_",
                    item,
                    sep = ""
                )
            )
        ResList[[paste0("Mappingsummary_", item, sep = "")]] <- summaryDF
    }

    # # Create the long DF summary if summary =TRUE
    if (summary) {
        for (item in to) {
            summary <-
                mapping_ambiguity(
                    data = TranslatedDF,
                    from = metadata_info[["InputID"]],
                    to = item,
                    grouping_variable = metadata_info[["grouping_variable"]],
                    summary = TRUE
                )[["summary"]]
            ResList[[paste0(
                "Mappingsummary_Long_",
                from,
                "-to-",
                item,
                sep = ""
            )]] <- summary
        }
    }


    # # ------------------ Save the results ------------------- ##


    save_res(
        inputlist_df = ResList,
        inputlist_plot = NULL,
        save_table = save_table,
        save_plot = NULL,
        path = Subfolder,
        file_name = "translate_id",
        core = FALSE,
        print_plot = FALSE
    )


    # # ------------------ Show Upset plot of the results ------------------- ##
    # # toDO: current issue with Lists vs doubles
    # if (plot==TRUE) {
    #   pk_list <- list(translated = TranslatedDF)
    #   metadata_info <- list(translated = to)
    #   print(pk_list)
    #   print(to)
    #   pk_comp_res <-
    #       compare_pk(
    #           data = pk_list,
    #           metadata_info = metadata_info,
    #           plot_name = "IDs available after
    #           ID Translation"
    #       )
    #   ## Add Upset plot
    #   ResList[["Translated_UpsetPlot"]] <- pk_comp_res$upset_plot
    # }

    # Return
    invisible(
        return(ResList)
    )
}


#
# Find additional potential IDs
#


#' @importFrom rlang exec !!!
#' @importFrom purrr map pmap_chr
#' @importFrom logger log_warn
#' @noRd
.log_df <- function(d) {

    match.call() %>% {deparse(.$d)} %>% log_warn('DF: %s', .)

    d %>% dim %>% {exec(log_warn, 'DF of size %i x %i', !!!.)}

    d %>%
    map(class) %>%
    sprintf('%s <%s>', names(.), .) %>%
    paste0(collapse = ', ') %>%
    log_warn('Columns: %s', .)

    log_warn('First 3 rows:')

    d %>% head(3L) %>% pmap_chr(~paste0(c(...), collapse = '|')) %>% log_warn()

}


#' Find additional potential IDs for  "kegg", "pubchem", "chebi", "hmdb"
#'
#' @param data dataframe with at least one column with the detected metabolite IDs (one
#'     ID per row).
#' @param metadata_info \emph{Optional: } Column name of metabolite IDs. \strong{Default =
#'     list(InputID="MetaboliteID")}
#' @param from ID type that is present in your data. Choose between "kegg", "pubchem",
#'     "chebi", "hmdb". \strong{Default = "hmdb"}
#' @param save_table \emph{Optional: } File types for the analysis results are: "csv",
#'     "xlsx", "txt". \strong{Default = "csv"}
#' @param path {Optional:} Path to the folder the results should be saved at.
#'     \strong{Default = NULL}
#'
#' @return Input DF with additional column including potential additional IDs.
#'
#' @examples
#' data(cellular_meta)
#' DetectedIDs <- cellular_meta %>% tidyr::drop_na()
#' Res <- equivalent_id(
#'     data = DetectedIDs,
#'     metadata_info = c(InputID = "HMDB"),
#'     from = "hmdb"
#' )
#'
#' @importFrom dplyr mutate select group_by ungroup distinct filter across
#' @importFrom dplyr rowwise if_else
#' @importFrom tidyr separate_rows unnest
#' @importFrom purrr map_chr map_lgl map_int
#' @importFrom tidyselect all_of starts_with
#' @importFrom rlang !!! !! := sym syms
#' @importFrom OmnipathR id_types translate_ids ensembl_organisms
#' @importFrom OmnipathR uniprot_organisms oma_organisms
#' @importFrom logger log_warn log_trace
#' @importFrom stringr str_to_lower str_split
#' @importFrom utils data
#' @export
equivalent_id <- function(
        data,
        metadata_info = c(InputID = "MetaboliteID"),
        from = "hmdb",
        save_table = "csv",
        path = NULL
    ) {
    # FUTURE: Once we have the structural similarity tool available in OmniPath,
    # we can start creating this function!

    # ## 1)
    # Check Measured ID's in prior knowledge


    # ## 2)
    # A user has one HMDB IDs for their measured metabolites (one ID per
    # measured peak) --> this is often the case as the user either gets a
    # trivial name and they have searched for the ID themselves or because the
    # facility only provides one ID at random
    # We have mapped the HMDB IDs with the pathways and 20 do not map
    # We want to check if it is because the pathways don't include them, or
    # because the user just gave the wrong ID by chance (i.e. They picked
    # D-Alanine, but the prior knowledge includes L-Alanine)

    # Do this by using structural information via  accessing the structural DB
    # in OmniPath!
    # Output is DF with the original ID column and a new column with additional
    # possible IDs based on structure

    # Is it possible to do this at the moment without structures, but by using
    # other prior knowledge?

    # NSE vs. R CMD check workaround
    InputID <- AdditionalID <- fromList <-
        PotentialAdditionalIDs <- AllIDs <- AllIDs.x <- AllIDs.y <- NULL

    metaproviz_init()

    # # ------------------  Check Input ------------------- ##
    # HelperFunction `check_param`
    check_param(
        data = data,
        data_num = FALSE,
        save_table = save_table
    )

    # Specific checks:
    if ("InputID" %in% names(metadata_info)) {
        if (!(metadata_info[["InputID"]] %in% colnames(data))) {
            message <-
                paste0(
                    "The ",
                    metadata_info[["InputID"]],
                    " column selected as InputID in metadata_info was ",
                    "not found in data. ",
                    "Please check your input."
                )
            log_trace(
                paste("Error ", message, sep = "")
            )
            stop(message)
        }
    }

    unknown_types <-
        id_types() %>%
        select(
            starts_with("in_")
        ) %>%
        unlist() %>%
        unique() %>%
        str_to_lower() %>%
        setdiff(from, .)

    if (length(unknown_types) > 0L) {
        msg <-
            sprintf(
                "The following ID types are not recognized: %s",
                paste(unknown_types, collapse = ", ")
            )
        log_warn(msg)
        warning(msg)
    }

    # Check that metadata_info[['InputID']] has no duplications within one group
    # --> should not be the case --> remove duplications and inform the user/
    # ask if they forget to set groupings column
    doublons <- data[duplicated(data[[metadata_info[["InputID"]]]]), ]

    if (nrow(doublons) > 0L) {
        data <-
            data %>%
            distinct(
                !!sym(metadata_info[["InputID"]]),
                .keep_all = TRUE
            )

        message <-
            sprintf(
                "The following IDs are duplicated and removed: %s",
                paste(
                    doublons[[metadata_info[["InputID"]]]],
                    collapse = ", "
                )
            )
        log_warn(message)
        warning(message)
    }

    # # ------------------  Create output folders and path ------------------- ##
    if (!is.null(save_table)) {
        folder <-
            save_path(
                folder_name = "PK",
                path = path
            )

        Subfolder <- file.path(folder, "EquivalentIDs")
        if (!dir.exists(Subfolder)) {
            dir.create(Subfolder)
        }
    }

    # # ------------------ Set the ID type for to ----------------- ##
    to <-
        case_when(
            from == "chebi" ~ "pubchem",
            # If to is "pubchem", choose "chebi"
            TRUE ~ "chebi"
        # For other cases, don't use a secondary column
        )

    message <-
        paste0(
            to,
            " is used to find additional potential IDs for ",
            from,
            ".",
            sep = ""
        )
    log_trace(message)
    message(message)

    # # ------------------ Load manual table ----------------- ##
    if (!(from == "kegg")) {
        data(
            list = "equivalent_features",
            package = "MetaProViz",
            envir = environment()
        )
        EquivalentFeatures <-
            "equivalent_features" %>%
            get(envir = environment()) %>%
            select(from)
    }

    # # ------------------ Translate from-to-to ------------------- ##
    TranslatedDF <-
        translate_ids(
            data,
            !!sym(metadata_info[["InputID"]]) := !!sym(from),
            !!!syms(to),
            # list of symbols, hence three !!!
            ramp = TRUE,
            expand = FALSE,
            quantify_ambiguity = FALSE,
            qualify_ambiguity = TRUE,
            ambiguity_groups = NULL,  # Can not be set to FALSE!
            # Checks within the groups, without it checks across groups
            ambiguity_summary = FALSE
        ) %>%
        select(
            all_of(intersect(names(.), names(data))), all_of(to)
        ) %>%
        mutate(
            across(
                all_of(to),
                ~ map_chr(., ~ paste(unique(.), collapse = ", "))
            )
        ) %>%
        group_by(
            !!sym(metadata_info[["InputID"]])
        ) %>%
        mutate(
            across(
                all_of(to),
                ~ paste(unique(.), collapse = ", "),
                .names = "{.col}"
            )
        ) %>%
        ungroup() %>%
        distinct() %>%
        mutate(
            across(all_of(to), ~ ifelse(. == "0", NA, .))
        )


    # # ------------------ Translate to-to-from ------------------- ##
    TranslatedDF_Long <-
        TranslatedDF %>%
        select(
            !!sym(metadata_info[["InputID"]]), !!sym(to)
        ) %>%
        rename(
            "InputID" = !!sym(metadata_info[["InputID"]])
        ) %>%
        separate_rows(
            !!sym(to),
            sep = ", "
        ) %>%
        mutate(
            across(all_of(to), ~ trimws(.))
        ) %>%
        filter(  # Remove extra spaces
            !!sym(to) != ""
        )
    # Remove empty entries

    OtherIDs <-
        translate_ids(
            TranslatedDF_Long,
            !!sym(to),
            !!sym(from),
            ramp = TRUE,
            expand = FALSE,
            quantify_ambiguity = FALSE,
            qualify_ambiguity = TRUE,
            ambiguity_groups = NULL,  # Can not be set to FALSE!
            # Checks within the groups, without it checks across groups
            ambiguity_summary = FALSE
        ) %>%
        select(
            "InputID", !!sym(to), !!sym(from)
        ) %>%
        distinct(
            InputID, !!sym(from),
            .keep_all = TRUE
        ) %>%
        # Remove duplicates based on InputID and from
        mutate(
            AdditionalID = if_else(InputID == !!sym(from), FALSE, TRUE)
        ) %>%
        select(
            "InputID", !!sym(from), "AdditionalID"
        ) %>%
        filter(
            AdditionalID
        ) %>%
        mutate(
            across(
                all_of(from),
                ~ map_chr(., ~ paste(unique(.), collapse = ", "))
            )
        ) %>%
        rowwise() %>%
        mutate(
            fromList = list(str_split(!!sym(from), ", \\s*")[[1]]),
            # Wrap in list
            SameAsInput = ifelse(
                any(fromList == InputID),
                InputID,
                NA_character_
            ),
            # Match InputID
            PotentialAdditionalIDs = paste(
                fromList[fromList != InputID],
                collapse = ", "
            )
        # Combine other IDs
        ) %>%
        ungroup() %>%
        select(
            InputID, PotentialAdditionalIDs, !!sym(from)
        ) %>%
        # Final selection
        rename(
            "AllIDs" = from
        )

    # # ------------------ Merge to Input ------------------- ##
    OtherIDs <-
        merge(
            data,
            OtherIDs,
            by.x = metadata_info[["InputID"]],
            by.y = "InputID",
            all.x = TRUE
        )

    # # ------------------- Add additional IDs -------------- ##

    if (exists("EquivalentFeatures")) {
        EquivalentFeatures$AllIDs <- EquivalentFeatures[[from]]
        EquivalentFeatures_Long <-
            EquivalentFeatures %>%
            separate_rows(
                !!sym(from),
                sep = ", "
            )

        OtherIDs <-
            merge(
                OtherIDs,
                EquivalentFeatures_Long,
                by.x = metadata_info[["InputID"]],
                by.y = "hmdb",
                all.x = TRUE
            ) %>%
            rowwise() %>%
            mutate(
                AllIDs = paste(
                    unique(
                        na.omit(
                            unlist(
                                str_split(
                                    paste(
                                        na.omit(
                                            c(
                                                AllIDs.x,
                                                AllIDs.y
                                            )
                                        ),
                                        collapse = ", "
                                    ),
                                    ", \\s*"
                                )
                            )
                        )
                    ),
                    collapse = ", "
                )
            ) %>%
            ungroup() %>%
            rowwise() %>%
            mutate(
                PotentialAdditionalIDs = paste(
                    setdiff(
                        unlist(
                            str_split(
                                AllIDs, ", \\s*"
                            )
                        ),
                        # Split merged_column into individual IDs
                        as.character(
                            !!sym(
                                # Split hmdb into individual IDs
                                metadata_info[["InputID"]]
                            )
                        )
                    ),
                    # Combine the remaining IDs back into a comma-separated
                    # string
                    collapse = ", "
                )
            ) %>%
            ungroup() %>%
            select(
                -AllIDs.x,
                -AllIDs.y
            )
    }

    # # ------------------- Fill empty cells -------------- ##
    OtherIDs <-
        OtherIDs %>%
        mutate(
            PotentialAdditionalIDs = ifelse(
                is.na(PotentialAdditionalIDs) | PotentialAdditionalIDs == "",
                NA,
                PotentialAdditionalIDs
            )
        ) %>%
        mutate(
            AllIDs = ifelse(
                is.na(AllIDs) | AllIDs == "",
                !!sym(metadata_info[["InputID"]]),
                AllIDs
            )
        )

    # # ------------------ Create count_id plot ------------------- ##
    # QC plot of before and after
    Before <-
        count_id(
            data = data,
            column = metadata_info[["InputID"]],
            save_plot = NULL,
            save_table = NULL,
            print_plot = FALSE,
            path = NULL
        )

    After <-
        count_id(
            data = OtherIDs,
            column = "AllIDs",
            save_plot = NULL,
            save_table = NULL,
            print_plot = FALSE,
            path = NULL
        )


    # # ------------------ Create Output ------------------- ##
    OutputDF <- OtherIDs

    # # ------------------ Save the results ------------------- ##
    ResList <- list("equivalent_id" = OutputDF)


    save_res(
        inputlist_df = ResList,
        inputlist_plot = NULL,
        save_table = save_table,
        save_plot = NULL,
        path = Subfolder,
        file_name = "equivalent_id",
        core = FALSE,
        print_plot = FALSE
    )

    return(invisible(OutputDF))
}

#
# Mapping Ambiguity
#

#' Create Mapping Ambiguities between two ID types
#'
#' @param data Translated DF from translate_id reults or dataframe with at least one
#'     column with the target metabolite ID and another MetaboliteID type. One
#'     of the IDs can only have one ID per row, the other ID can be either
#'
# separated


#' by comma or a list. Optional: add other columns such as source (e.g. term).
#'
#' @param to Column name of original metabolite identifier in data. Here should only
#'     have one ID per row.
#' @param from Column name of the secondary or translated metabolite identifier in
#'     data. Here can be multiple IDs per row either separated by comma " ," or
#'     a list of IDs.
#' @param grouping_variable \emph{Optional: } If NULL no groups are used. If TRUE provide column
#'     name in data containing the grouping_variable and features are grouped.
#'     \strong{Default = NULL}
#' @param summary \emph{Optional: } If TRUE a long summary tables are created.
#'     \strong{Default = FALSE}
#' @param save_table \emph{Optional: } File types for the analysis results are: "csv",
#'     "xlsx", "txt". \strong{Default = "csv"}
#' @param path {Optional:} Path to the folder the results should be saved at.
#'     \strong{Default = NULL}
#'
#' @return List with at least 4 DFs: 1-3) from-to-to: 1. MappingIssues, 2.
#'     MappingIssues summary, 3. Long summary (If summary=TRUE) & 4-6)
#'     to-to-from: 4. MappingIssues, 5. MappingIssues summary, 6. Long summary
#'     (If summary=TRUE) & 7) Combined summary table (If summary=TRUE)
#'
#' @examples
#' \dontrun{
#' KEGG_Pathways <- metsigdb_kegg()
#' InputDF <- translate_id(
#'     data = KEGG_Pathways,
#'     metadata_info = c(
#'         InputID = "MetaboliteID",
#'         grouping_variable = "term"
#'     ),
#'     from = c("kegg"),
#'     to = c("pubchem")
#' )[["TranslatedDF"]]
#' Res <- mapping_ambiguity(
#'     data = InputDF,
#'     from = "MetaboliteID",
#'     to = "pubchem",
#'     grouping_variable = "term",
#'     summary = TRUE
#' )
#' }
#'
#' @importFrom dplyr mutate bind_cols bind_rows
#' @importFrom rlang !!! !! := sym syms
#' @importFrom OmnipathR ambiguity
#' @importFrom logger log_trace
#' @importFrom tidyr unnest
#' @importFrom purrr map_chr map_int map_lgl
#' @importFrom dplyr first
#' @export
mapping_ambiguity <- function(
data,
from,
to,
grouping_variable = NULL,
summary = FALSE,
save_table = "csv",
path = NULL
) {
    # NSE vs. R CMD check workaround
    UniqueID <- NULL

    metaproviz_init()
    # # ------------------  Check Input ------------------- ##
    # HelperFunction `check_param`
    check_param(
        data = data,
        data_num = FALSE,
        save_table = save_table
    )

    # Specific checks:
    if (!(from %in% colnames(data))) {
        message <-
            paste0(
                from,
                " column was not found in data. Please check your input."
            )
        log_trace(
            paste(
                "Error ", message,
                sep = ""
            )
        )
        stop(message)
    }

    if (!(to %in% colnames(data))) {
        message <-
            paste0(
                to,
                " column was not found in data. Please check your input."
            )
        log_trace(
            paste(
                "Error ", message,
                sep = ""
            )
        )
        stop(message)
    }

    if (!is.null(grouping_variable)) {
        if (!(grouping_variable %in% colnames(data))) {
            message <-
                paste0(
                    grouping_variable,
                    " column was not found in data. Please check your input."
                )
            log_trace(
                paste(
                    "Error ", message,
                    sep = ""
                )
            )
            stop(message)
        }
    }

    if (!is.logical(summary)) {
        message <-
            paste0(
                "Check input. The summary parameter should be either =TRUE or ",
                "=FALSE."
            )
        log_trace(
            paste(
                "Error ", message,
                sep = ""
            )
        )
        stop(message)
    }

    # # ------------------  General checks of wrong occurences
    # # ------------------- ##
    # Task 1: Check that from has no duplications within one group --> should
    # not be the case --> remove duplications and inform the user/ ask if they
    # forget to set groupings column
    # Task 2: Check that from has the same items in to across the different
    # entries (would be in different Groupings, otherwise there should not be
    # any duplications) --> List of Miss-Mappings across terms

    # FYI: The above can not happen if our translateID function was used, but
    # may be the case when the user has done something manually before


    # # ------------------  Create output folders and path ------------------- ##
    if (!is.null(save_table)) {
        folder <-
            save_path(
                folder_name = "PK",
                path = path
            )

        Subfolder <- file.path(folder, "MappingAmbiguities")
        if (!dir.exists(Subfolder)) {
            dir.create(Subfolder)
        }
    }

    # ###########################################################################
    # # ------------------  Prepare Input data ------------------- ##
    # If the user provides a DF where the to column is a list of IDs, then we
    # can use it right away
    # If the to column is not a list of IDs, but a character column, we need to
    # convert it into a list of IDs
    if (is.character(data[[to]])) {
        data[[to]] %<>% strsplit(", ") %>% lapply(as.character)
    }

    # # ------------------  Perform ambiguity mapping ------------------- ##
    # 1. from-to-to: OriginalID-to-TranslatedID
    # 2. from-to-to: TranslatedID-to-OriginalID
    Comp <- list(
        list(from = from, to = to),
        list(from = to, to = from)
    )

    ResList <- list()
    for (comp in seq_along(Comp)) {
        # Run Omnipath ambiguity
        ResList[[
            paste0(
                Comp[[comp]]$from,
                "-to-",
                Comp[[comp]]$to,
                sep = ""
            )
        ]] <-
            data %>%
            unnest(  # unlist the columns in case they are not expaned
                cols = all_of(Comp[[comp]]$from)
            ) %>%
            # Remove NA values, otherwise they are counted as column is
            # character
            filter(
                !is.na(
                    !!sym(Comp[[comp]]$from)
                )
            ) %>%
            ambiguity(
                from_col = !!sym(Comp[[comp]]$from),
                to_col = !!sym(Comp[[comp]]$to),
                groups = grouping_variable,
                quantify = TRUE,
                qualify = TRUE,
                global = TRUE,
                # across groups will be done additionally --> suffix
                # _AcrossGroup
                summary = TRUE,
                # summary of the mapping column
                expand = TRUE
            )

        # Extract summary table:
        ResList[[
            paste0(
                Comp[[comp]]$from,
                "-to-",
                Comp[[comp]]$to,
                "_summary",
                sep = ""
            )
        ]] <-
            ResList[[
                paste0(
                    Comp[[comp]]$from,
                    "-to-",
                    Comp[[comp]]$to,
                    sep = ""
                )
            ]] %>%
            attr(
                paste0(
                    "ambiguity_",
                    Comp[[comp]]$from,
                    "_",
                    Comp[[comp]]$to,
                    sep = ""
                )
            )

        # #######################################################################
        if (summary) {
            if (!is.null(grouping_variable)) {
                from_to_v1 <-
                    sprintf(
                        '%s_%s',
                        Comp[[comp]]$from,
                        Comp[[comp]]$to
                    )
                from_to_v2 <-
                    sprintf(
                        '%s-to-%s',
                        Comp[[comp]]$from,
                        Comp[[comp]]$to
                    )
                from_to_v3 <-
                    sprintf(
                        '%s_to_%s',
                        Comp[[comp]]$from,
                        Comp[[comp]]$to
                    )
                from_to_amb <- sprintf('%s_ambiguity_bygroup', from_to_v1)
                from_to_amb_bg <- sprintf('%s_ambiguity', from_to_v1)
                from_to_long <- sprintf('%s_Long', from_to_v2)
                from_to_acr_grp <-
                    sprintf(
                        'AcrossGroupMappingIssue(%s)',
                        from_to_v3
                    )
                from_to_count <- sprintf('Count(%s)', from_to_v3)
                # Add further information we need to summarise the table and
                # combine Original-to-Translated and Translated-to-Original
                # If we have a grouping_variable we need to combine it with the
                # MetaboliteID before merging
                ResList[[from_to_long]] <-
                    ResList[[from_to_v2]] %>%
                    unnest(
                        cols = all_of(Comp[[comp]]$from)
                    ) %>%
                    mutate(
                        !!sym(from_to_acr_grp) := case_when(
                            (
                                !!sym(from_to_amb_bg) !=
                                !!sym(from_to_amb)
                            ) ~ "TRUE",
                            TRUE ~ "FALSE"
                        )
                    ) %>%
                    group_by(
                        !!sym(Comp[[comp]]$from),
                        !!sym(grouping_variable)
                    ) %>%
                    mutate(
                        !!sym(
                            Comp[[comp]]$to
                        ) := ifelse(
                            !!sym(Comp[[comp]]$from) == 0L,
                            NA,  # Or another placeholder
                            paste(
                                unique(!!sym(Comp[[comp]]$to)),
                                collapse = ", "
                            )
                        )
                    ) %>%
                    mutate(
                        !!sym(from_to_count) := ifelse(
                            all(!!sym(Comp[[comp]]$to) == 0L),
                            0L,
                            n()
                        )
                    ) %>%
                    ungroup() %>%
                    distinct() %>%
                    unite(
                        !!sym(from_to_v3),
                        c(
                            Comp[[comp]]$from,
                            Comp[[comp]]$to
                        ),
                        sep = " --> ",
                        remove = FALSE
                    ) %>%
                    separate_rows(
                        !!sym(Comp[[comp]]$to),
                        sep = ", "
                    ) %>%
                    unite(
                        UniqueID,
                        c(from, to, grouping_variable),
                        sep = "_",
                        remove = FALSE
                    ) %>%
                    distinct()
            } else {
                ResList[[paste0(
                    Comp[[comp]]$from,
                    "-to-",
                    Comp[[comp]]$to,
                    "_Long",
                    sep = ""
                )]] <-
                    ResList[[paste0(
                        Comp[[comp]]$from,
                        "-to-",
                        Comp[[comp]]$to,
                        sep = ""
                    )]] %>%
                    unnest(
                        cols = all_of(Comp[[comp]]$from)
                    ) %>%
                    group_by(
                        !!sym(Comp[[comp]]$from)
                    ) %>%
                    mutate(
                        !!sym(
                            Comp[[comp]]$to
                        ) := ifelse(
                            !!sym(Comp[[comp]]$from) == 0L,
                            NA,  # Or another placeholder
                            paste(
                                unique(!!sym(Comp[[comp]]$to)),
                                collapse = ", "
                            )
                        )
                    ) %>%
                    mutate(
                        !!sym(
                            paste0(
                                "Count(",
                                Comp[[comp]]$from,
                                "_to_",
                                Comp[[comp]]$to,
                                ")"
                            )
                        ) := ifelse(
                            all(!!sym(Comp[[comp]]$to) == 0L),
                            0,
                            n()
                        )
                    ) %>%
                    ungroup() %>%
                    distinct() %>%
                    unite(
                        !!sym(
                            paste0(
                                Comp[[comp]]$from,
                                "_to_",
                                Comp[[comp]]$to
                            )
                        ),
                        c(
                            Comp[[comp]]$from,
                            Comp[[comp]]$to
                        ),
                        sep = " --> ",
                        remove = FALSE
                    ) %>%
                    separate_rows(
                        !!sym(Comp[[comp]]$to),
                        sep = ", "
                    ) %>%
                    unite(
                        UniqueID,
                        c(from, to),
                        sep = "_",
                        remove = FALSE
                    ) %>%
                    distinct() %>%
                    mutate(
                        !!sym(
                            paste0(
                                "AcrossGroupMappingIssue(",
                                from,
                                "_to_",
                                to,
                                ")",
                                sep = ""
                            )
                        ) := NA
                    )
            }
        }

        # Add NA metabolite maps back if they do exist:
        Removed <-
            data %>%
            unnest(  # unlist the columns in case they are not expaned
                cols = all_of(Comp[[comp]]$from)
            ) %>%
            filter(
                is.na(!!sym(Comp[[comp]]$from))
            )

        if (nrow(Removed) > 0L) {
            ResList[[paste0(
                Comp[[comp]]$from,
                "-to-",
                Comp[[comp]]$to,
                sep = ""
            )]] <-
                bind_rows(
                    ResList[[paste0(
                        Comp[[comp]]$from,
                        "-to-",
                        Comp[[comp]]$to,
                        sep = ""
                    )]],
                    test <-
                        Removed %>%
                        bind_cols(
                            setNames(
                                as.list(
                                    rep(
                                        NA,
                                        length(
                                            setdiff(
                                                names(
                                                    ResList[[paste0(
                                                        Comp[[comp]]$from,
                                                        "-to-",
                                                        Comp[[comp]]$to,
                                                        sep = ""
                                                    )]]
                                                ),
                                                names(Removed)
                                            )
                                        )
                                    )
                                ),
                                setdiff(
                                    names(
                                        ResList[[paste0(
                                            Comp[[comp]]$from,
                                            "-to-",
                                            Comp[[comp]]$to,
                                            sep = ""
                                        )]]
                                    ),
                                    names(Removed)
                                )
                            )
                        )
                )
        }
    }

    # # ------------------ Create summaryTable ------------------- ##
    if (summary) {
        # Combine the two tables
        summary <-
            merge(
                x = ResList[[
                    paste0(
                        from, "-to-", to, "_Long",
                        sep = ""
                    )
                ]][
,
                    c(
                        "UniqueID",
                        paste0(from, "_to_", to),
                        paste0("Count(", from, "_to_", to, ")"),
                        paste0(
                            "AcrossGroupMappingIssue(",
                            from,
                            "_to_",
                            to,
                            ")",
                            sep = ""
                        )
                    )
                ],
                y = ResList[[
                    paste0(
                        to, "-to-", from, "_Long",
                        sep = ""
                    )
                ]][
,
                    c(
                        "UniqueID",
                        paste0(to, "_to_", from),
                        paste0("Count(", to, "_to_", from, ")"),
                        paste0(
                            "AcrossGroupMappingIssue(",
                            to,
                            "_to_",
                            from,
                            ")",
                            sep = ""
                        )
                    )
                ],
                by = "UniqueID",
                all = TRUE
            ) %>%
            separate(
                UniqueID,
                into = c(from, to, grouping_variable),
                sep = "_",
                remove = FALSE
            ) %>%
            distinct()

        # Add relevant mapping information
        summary <-
            summary %>%
            mutate(
                Mapping = case_when(
                    !!sym(paste0("Count(", from, "_to_", to, ")")) == 1L &
                    !!sym(paste0("Count(", to, "_to_", from, ")")) == 1L
                        ~"one-to-one",
                    !!sym(paste0("Count(", from, "_to_", to, ")")) > 1L &
                    !!sym(paste0("Count(", to, "_to_", from, ")")) == 1L
                        ~"one-to-many",
                    !!sym(paste0("Count(", from, "_to_", to, ")")) > 1L &
                    !!sym(paste0("Count(", to, "_to_", from, ")")) > 1L
                        ~"many-to-many",
                    !!sym(paste0("Count(", from, "_to_", to, ")")) == 1L &
                    !!sym(paste0("Count(", to, "_to_", from, ")")) > 1L
                        ~"many-to-one",
                    !!sym(paste0("Count(", from, "_to_", to, ")")) >= 1L &
                    !!sym(paste0("Count(", to, "_to_", from, ")")) == NA
                        ~"one-to-none",
                    !!sym(paste0("Count(", from, "_to_", to, ")")) >= 1L &
                    is.na(!!sym(paste0("Count(", to, "_to_", from, ")")))
                        ~"one-to-none",
                    !!sym(paste0("Count(", from, "_to_", to, ")")) == NA &
                    !!sym(paste0("Count(", to, "_to_", from, ")")) >= 1L
                        ~"none-to-one",
                    is.na(!!sym(paste0("Count(", from, "_to_", to, ")"))) &
                    !!sym(paste0("Count(", to, "_to_", from, ")")) >= 1L
                        ~"none-to-one",
                    TRUE
                        ~NA
                )
            ) %>%
            mutate(
                !!sym(
                    paste0("Count(", from, "_to_", to, ")")
                ) := replace_na(
                    !!sym(
                        paste0("Count(", from, "_to_", to, ")")
                    ),
                    0L
                )
            ) %>%
            mutate(
                !!sym(
                    paste0("Count(", to, "_to_", from, ")")
                ) := replace_na(
                    !!sym(
                        paste0("Count(", to, "_to_", from, ")")
                    ),
                    0L
                )
            )

        ResList[["summary"]] <- summary
    }

    # # ------------------ Save the results ------------------- ##


    save_res(
        inputlist_df = ResList,
        inputlist_plot = NULL,
        save_table = save_table,
        save_plot = NULL,
        path = Subfolder,
        file_name = "mapping_ambiguity",
        core = FALSE,
        print_plot = FALSE
    )


    # Return
    invisible(
        return(ResList)
    )
}

#
# Check Measured ID's in prior knowledge
#

#' Check and summarize relationship between prixor knowledge to measured
#'
#' features
#'
#' @param data dataframe with at least one column with the detected metabolite IDs
#'     (e.g. HMDB). If there are multiple IDs per detected peak, please
#'     separate them by comma ("," or ", " or chr list). If there is a main ID
#'     and additional IDs, please provide them in separate columns.
#' @param input_pk dataframe with at least one column with the metabolite ID (e.g. HMDB)
#'     that need to match data metabolite IDs "source" (e.g. term). If there
#'     are multiple IDs, as the original pathway IDs (e.g. KEGG) where
#'     translated (e.g. to HMDB), please separate them by comma ("," or ", " or
#'     chr list).
#' @param metadata_info Colum name of Metabolite IDs in data and input_pk as well as column name
#'     of grouping_variable in input_pk. \strong{Default = c(InputID="HMDB",
#'     PriorID="HMDB", grouping_variable="term")}
#' @param save_table \emph{Optional: } File types for the analysis results are: "csv",
#'     "xlsx", "txt". \strong{Default = "csv"}
#' @param path {Optional:} Path to the folder the results should be saved at.
#'     \strong{Default = NULL}
#'
#' @return A \code{list} with three elements: \itemize{ \item \code{data_summary} —
#'     a data frame summarising matching results per input ID, including
#'     counts, conflicts, and recommended actions. \item
#'     \code{GroupingVariable_summary} — a detailed data frame showing matches
#'     grouped by the specified variable, with conflict annotations. \item
#'     \code{data_long} — a merged data frame of prior knowledge IDs and
#'     detected IDs in long format. }
#'
#' @examples
#' \dontrun{
#' data(cellular_meta)
#' DetectedIDs <- cellular_meta %>%
#'     dplyr::select("Metabolite", "HMDB") %>%
#'     tidyr::drop_na()
#' input_pathway <- translate_id(
#'     data = metsigdb_kegg(),
#'     metadata_info = c(
#'         InputID = "MetaboliteID",
#'         grouping_variable = "term"
#'     ),
#'     from = c("kegg"),
#'     to = c("hmdb")
#' )[["TranslatedDF"]] %>% tidyr::drop_na()
#' Res <- checkmatch_pk_to_data(
#'     data = DetectedIDs,
#'     input_pk = input_pathway,
#'     metadata_info = c(
#'         InputID = "HMDB",
#'         PriorID = "hmdb",
#'         grouping_variable = "term"
#'     )
#' )
#' }
#'
#' @importFrom dplyr cur_group_id filter mutate n_distinct row_number select
#' @importFrom logger log_trace
#' @importFrom tibble tibble
#' @importFrom tidyr separate_rows
#' @importFrom rlang !!! !! := sym syms
#' @importFrom purrr map_chr
#' @importFrom stringr str_split
#' @importFrom dplyr first
#' @export
checkmatch_pk_to_data <- function(
    data,
    input_pk,
    metadata_info = c(
        InputID = "HMDB",
        PriorID = "HMDB",
        grouping_variable = "term"
        ),
    save_table = "csv",
    path = NULL
) {
    # NSE vs. R CMD check workaround
    .data <- OriginalGroup_PK <- OriginalGroup_data <-
        Count_FeatureIDs_to_GroupingVariable <- GroupConflict_Notes <-
        ActionRequired <- matches <- matches_count <- Group_Conflict_Notes <-
        row_id <- NULL

    # # ------------ Create log file ----------- ##
    metaproviz_init()

    # # ------------ Check Input files ----------- ##

    # # data:
    if ("InputID" %in% names(metadata_info)) {
        if (!(metadata_info[["InputID"]] %in% colnames(data))) {
            message <-
                paste0(
                    "The ",
                    metadata_info[["InputID"]],
                    " column selected as InpuID in metadata_info was not found in ",
                    "data. Please check your input."
                )
            log_trace(
                paste("Error ", message, sep = "")
            )
            stop(message)
        }
    } else {
        message <-
            paste0(
                "No ",
                metadata_info[["InputID"]],
                " provided. Please check your input."
            )
        log_trace(
            paste("Error ", message, sep = "")
        )
        stop(message)
    }

    # ## This is after the main input checks (before NA removal), so we will save
    # ## original df here for later merging to get the Null and duplicates back.
    data_Original <- data

    if (sum(is.na(data[[metadata_info[["InputID"]]]])) >= 1L) {
        # remove NAs:
        message <-
            paste0(
                sum(is.na(data[[metadata_info[["InputID"]]]])),
                " NA values were removed from column ",
                metadata_info[["InputID"]]
            )
        log_trace(
            paste("Warning: ", message, sep = "")
        )

        data <-
            data %>%
            filter(
                !is.na(.data[[metadata_info[["InputID"]]]])
            )

        warning(message)
    }

    if (nrow(data) - nrow(distinct(data, .data[[metadata_info[["InputID"]]]])) >= 1L) {
        # Remove duplicate IDs
        message <-
            paste0(
                nrow(data) - nrow(distinct(data, .data[[metadata_info[["InputID"]]]])),
                " duplicated IDs were removed from column ",
                metadata_info[["InputID"]]
            )
        log_trace(
            paste("Warning: ", message, sep = "")
        )

        data <-
            data %>%
            distinct(
                .data[[metadata_info[["InputID"]]]],
                .keep_all = TRUE
            )

        warning(message)
    }

    data_MultipleIDs <-
        any(
            grepl(
                ", \\s*",
                data[[metadata_info[["InputID"]]]]
            ) |  # Comma-separated
                map_lgl(
                    data[[metadata_info[["InputID"]]]],
                    function(x) {
                        if (grepl("^c\\(|^list\\(", x)) {
                            parsed <- tryCatch(
                                eval(parse(text = x)),
                                error = function(e) NULL
                            )
                            return(
                                is.list(parsed) &&
                                length(parsed) > 1L ||
                                is.vector(parsed) &&
                                length(parsed) > 1L
                            )
                        }
                        FALSE
                    }
                )
        )

    # # input_pk:
    if ("PriorID" %in% names(metadata_info)) {
        if (!(metadata_info[["PriorID"]] %in% colnames(input_pk))) {
            message <-
                paste0(
                    "The ",
                    metadata_info[["PriorID"]],
                    " column selected as InpuID in metadata_info was not found in ",
                    "input_pk. Please check your input."
                )
            log_trace(
                paste("Error ", message, sep = "")
            )
            stop(message)
        }
    } else {
        message <-
            paste0(
                "No ",
                metadata_info[["PriorID"]],
                " provided. Please check your input."
            )
        log_trace(
            paste("Error ", message, sep = "")
        )
        stop(message)
    }

    # ## This is after the main input checks (before NA removal), so we will save
    # ## original df here for later merging to get the Null and duplicates back.
    prior_knowledge_Original <- input_pk

    if (sum(is.na(input_pk[[metadata_info[["PriorID"]]]])) >= 1L) {
        # remove NAs:
        message <-
            paste0(
                sum(is.na(input_pk[[metadata_info[["PriorID"]]]])),
                " NA values were removed from column ",
                metadata_info[["PriorID"]]
            )
        log_trace(
            paste("Warning: ", message, sep = "")
        )

        input_pk <-
            input_pk %>%
            filter(
                !is.na(.data[[metadata_info[["PriorID"]]]])
            )

        warning(message)
    }

    if ("grouping_variable" %in% names(metadata_info)) {
        # Add grouping_variable
        if (!(metadata_info[["grouping_variable"]] %in% colnames(input_pk))) {
            message <-
                paste0(
                    "The ",
                    metadata_info[["grouping_variable"]],
                    " column selected as InpuID in metadata_info was not found in ",
                    "input_pk. Please check your input."
                )
            log_trace(
                paste("Error ", message, sep = "")
            )
            stop(message)
        }
    } else {
        # Add grouping_variable
        metadata_info["grouping_variable"] <- "GroupingVariable"
        input_pk["GroupingVariable"] <- "OneGroup"

        message <-
            paste0(
                "No metadata_info grouping_variable provided. ",
                "If this was not intentional, please check your input."
            )
        log_trace(message)
        message(message)
    }

    n_removed <-
        nrow(input_pk) -
        nrow(
            distinct(
                input_pk,
                .data[[metadata_info[["PriorID"]]]],
                .data[[metadata_info[["grouping_variable"]]]]
            )
        )
    if (n_removed >= 1L) {
        # Remove duplicate IDs
        message <-
            paste0(
                n_removed,
                " duplicated IDs were removed from PK column ",
                metadata_info[["PriorID"]]
            )
        log_trace(
            paste("Warning: ", message, sep = "")
        )

        input_pk <-
            input_pk %>%
            distinct(
                .data[[metadata_info[["PriorID"]]]],
                !!sym(metadata_info[["grouping_variable"]]),
                .keep_all = TRUE
            ) %>%
            group_by(
                !!sym(metadata_info[["PriorID"]])
            ) %>%
            mutate(
                across(
                    everything(),
                    ~ if (is.character(.)) paste(unique(.), collapse = ", ")
                )
            ) %>%
            ungroup() %>%
            distinct(
                .data[[metadata_info[["PriorID"]]]],
                .keep_all = TRUE
            )

        warning(message)
    }

    # Check if multiple IDs are present:
    PK_MultipleIDs <-
        any(
            grepl(
                ", \\s*",
                input_pk[[metadata_info[["PriorID"]]]]
            ) |  # Comma-separated
            map_lgl(
                input_pk[[metadata_info[["PriorID"]]]],
                function(x) {
                    if (grepl("^c\\(|^list\\(", x)) {
                        parsed <- tryCatch(
                            eval(parse(text = x)),
                            error = function(e) NULL
                        )
                        return(
                            is.list(parsed) &&
                            length(parsed) > 1L ||
                            is.vector(parsed) &&
                            length(parsed) > 1L
                        )
                    }
                    FALSE
                }
            )
        )

    # # ------------ Create Results output folder ----------- ##
    if (!is.null(save_table)) {
        folder <-
            save_path(
                folder_name = "PK",
                path = path
            )
        Subfolder <- file.path(folder, "CheckMatchID_Detected-to-PK")
        if (!dir.exists(Subfolder)) {
            dir.create(Subfolder)
        }
    }

    # ###########################################################################
    # # ------------ Check how IDs match and if needed remove unmatched IDs
    # # ----------- ##

    # 1.  Create long DF
    create_long_df <-
    function(
        df,
        id_col,
        df_name
    ) {
            df %>%
            mutate(
                row_id = row_number()
            ) %>%
            mutate(  # Store original values
                !!paste0("OriginalEntry_", df_name, sep = "") :=
                    !!sym(id_col)
            ) %>%
            separate_rows(
                !!sym(id_col),
                sep = ", \\s*"
            ) %>%
            group_by(
                row_id
            ) %>%
            mutate(
                !!(paste0("OriginalGroup_", df_name, sep = "")) :=
                    paste0(df_name, "_", cur_group_id())
            ) %>%
            ungroup()
        }

    if (data_MultipleIDs) {
        data_long <-
            create_long_df(
                data,
                metadata_info[["InputID"]],
                "data"
            ) %>%
            select(
                metadata_info[["InputID"]],
                "OriginalEntry_data",
                OriginalGroup_data
            )
    } else {
        data_long <-
            data %>%
            mutate(
                OriginalGroup_data := paste0("data_", row_number())
            ) %>%
            select(
                metadata_info[["InputID"]],
                OriginalGroup_data
            )
        data_long$`OriginalEntry_data` <-
            data_long[[metadata_info[["InputID"]]]]
    }

    if (PK_MultipleIDs) {
        PK_long <-
            create_long_df(
                input_pk,
                metadata_info[["PriorID"]],
                "PK"
            ) %>%
            select(
                metadata_info[["PriorID"]],
                "OriginalEntry_PK",
                OriginalGroup_PK,
                metadata_info[["grouping_variable"]]
            )
    } else {
        PK_long <-
            input_pk %>%
            mutate(
                OriginalGroup_PK := paste0("PK_", row_number())
            ) %>%
            select(
                metadata_info[["PriorID"]],
                OriginalGroup_PK,
                metadata_info[["grouping_variable"]]
            )
    # PK_long$`OriginalEntry_PK` <- PK_long[[metadata_info[["PriorID"]]]]
    }

    # 2. Merge DF
    merged_df <-
        merge(
            PK_long,
            data_long,
            by.x = metadata_info[["PriorID"]],
            by.y = metadata_info[["InputID"]],
            all = TRUE
        ) %>%
        distinct(
            !!sym(metadata_info[["PriorID"]]),
            OriginalGroup_data,
            !!sym(metadata_info[["grouping_variable"]]),
            .keep_all = TRUE
        )

    # 3. Create summary table
    Values_data <-
        unique(
            data[[metadata_info[["InputID"]]]]
        )
    Values_PK <-
        unique(
            PK_long[[metadata_info[["PriorID"]]]]
        )

    summary_df <-
        tibble(
            !!sym(metadata_info[["InputID"]]) := Values_data,
            found_match_in_PK = NA,
            matches = NA_character_,
            match_overlap_percentage = NA_real_,
            original_count = NA_integer_,
            matches_count = NA_integer_
        )

    # Populate the summary data frame
    for (i in seq_along(Values_data)) {
        # Handle NA case explicitly
        if (is.na(Values_data[i])) {
            summary_df$original_count[i] <- 0
            summary_df$matches_count[i] <- 0
            summary_df$match_overlap_percentage[i] <- NA
            # could also set it to FALSE but making NA for now for plotting
            summary_df$found_match_in_PK[i] <- NA
            summary_df$matches[i] <- NA
        } else {
            # Split each cell into individual entries and trim whitespace
            entries <-
                trimws(
                    unlist(
                        strsplit(
                            as.character(Values_data[i]),
                            ", \\s*"
                        )
                    )
                )
            # delimiter = "," or ", "

            # Identify which entries are in the lookup set
            matched <- entries[entries %in% Values_PK]

            # Determine if any match was found
            summary_df$found_match_in_PK[i] <- length(matched) > 0L

            # Concatenate matched entries into a single string
            summary_df$matches[i] <- paste(matched, collapse = ", ")

            # Calculate and store counts
            summary_df$original_count[i] <- length(entries)
            summary_df$matches_count[i] <- length(matched)

            # Calculate fraction: matched entries / total entries
            if (length(entries) > 0L) {
                summary_df$match_overlap_percentage[i] <-
                    (length(matched) / length(entries)) * 100
            } else {
                summary_df$match_overlap_percentage[i] <- NA
            }
        }
    }

    summary_df <-
        merge(
            x = summary_df,
            y = merged_df %>%
            select(
                    -c(OriginalGroup_PK, OriginalGroup_data)
                ),
            by.x = metadata_info[["InputID"]],
            by.y = "OriginalEntry_data",
            all.x = TRUE
        ) %>%
        group_by(
            !!sym(metadata_info[["InputID"]]),
            !!sym(metadata_info[["grouping_variable"]])
        ) %>%
        mutate(
            Count_FeatureIDs_to_GroupingVariable = case_when(
                is.na(!!sym(metadata_info[["grouping_variable"]]))
                    ~NA_integer_,
                TRUE
                    ~n_distinct(!!sym(metadata_info[["PriorID"]]), na.rm = TRUE)
            )
        ) %>%
        mutate(
            Group_Conflict_Notes = case_when(
                any(
                    Count_FeatureIDs_to_GroupingVariable > 1L,
                    na.rm = TRUE
                ) ~ ifelse(
                    all(
                        is.na(!!sym(metadata_info[["grouping_variable"]]))
                    ),
                    "None",
                    paste(
                        na.omit(
                            unique(
                                !!sym(metadata_info[["grouping_variable"]])
                            )
                        ),
                        collapse = " || "
                    )
                ),
                TRUE ~ "None"
            )
        ) %>%
        ungroup() %>%
        mutate(
            matches = ifelse(matches == "", NA, matches)
        )

    summary_df_short <-
        summary_df %>%
        group_by(
            !!sym(metadata_info[["InputID"]])
        ) %>%
        mutate(
            Unique_GroupingVariable_count = n_distinct(
                !!sym(metadata_info["grouping_variable"]),
                na.rm = TRUE
            ),
            !!sym(metadata_info[["grouping_variable"]]) := list(!!sym(metadata_info[["grouping_variable"]])),
            Count_FeatureIDs_to_GroupingVariable = list(Count_FeatureIDs_to_GroupingVariable),
            Group_Conflict_Notes = paste(unique(Group_Conflict_Notes), collapse = " || ")
        ) %>%
        ungroup() %>%
        distinct(
            !!sym(metadata_info[["InputID"]]),
            !!sym(metadata_info[["grouping_variable"]]),
            .keep_all = TRUE
        ) %>%
        mutate(
            ActionRequired = case_when(
                original_count >= 1L &
                matches_count == 1L &
                Unique_GroupingVariable_count >= 1L
                    ~"None",
                original_count >= 1L &
                matches_count == 0L &
                Unique_GroupingVariable_count >= 0L
                    ~"None",
                original_count > 1L &
                matches_count > 1L &
                Unique_GroupingVariable_count == 1L
                    ~"Check",
                original_count > 1L &
                matches_count > 1L &
                Unique_GroupingVariable_count > 1L
                    ~"Check",
                TRUE
                    ~NA_character_
            )
        ) %>%
        select(
            -!!sym(metadata_info[["PriorID"]])
        ) %>%
        mutate(
            matches = ifelse(matches == "", NA, matches),
            !!sym(metadata_info[["grouping_variable"]]) :=
                ifelse(
                    !!sym(metadata_info[["grouping_variable"]]) == "",
                    NA,
                    !!sym(metadata_info[["grouping_variable"]])
                )
        ) %>%
        mutate(
            InputID_select = case_when(
                original_count == 1L &
                matches_count <= 1L
                    ~metadata_info[["InputID"]],
                original_count >= 2L &
                matches_count == 0L
                    ~str_split(metadata_info[["InputID"]], ", \\s*") %>%
                    map_chr(first),
                original_count >= 2L &
                matches_count == 1L
                    ~matches,
                TRUE
                    ~NA_character_
            ),
            Action_Specific = case_when(
                matches_count >= 2L &
                Group_Conflict_Notes == "None"
                    ~"KeepEachID",
                matches_count >= 2L &
                Group_Conflict_Notes != "None"
                    ~"KeepOneID",
                TRUE
                    ~"None"
            )
        )

    # 4. Messages and summarise
    message0 <-
        paste0(
            "data has multiple IDs per measurement = ",
            data_MultipleIDs,
            ". input_pk has multiple IDs per entry = ",
            PK_MultipleIDs,
            ".",
            sep = ""
        )
    message1 <-
        paste0(
            "data has ",
            n_distinct(unique(data[[metadata_info[["InputID"]]]])),
            " unique entries with ",
            n_distinct(unique(data_long[[metadata_info[["InputID"]]]])),
            " unique ",
            metadata_info[["InputID"]],
            " IDs. Of those IDs, ",
            nrow(summary_df_short %>% filter(matches_count >= 1L)),
            " match, which is ",
            (
                nrow(summary_df_short %>% filter(matches_count >= 1L)) /
                n_distinct(unique(data_long[[metadata_info[["InputID"]]]])) *
                100
            ),
            "%.",
            sep = ""
        )
    message2 <-
        paste0(
            "input_pk has ",
            n_distinct(input_pk[[metadata_info[["PriorID"]]]]),
            " unique entries with ",
            n_distinct(PK_long[[metadata_info[["PriorID"]]]]),
            " unique ",
            metadata_info[["PriorID"]],
            " IDs. Of those IDs, ",
            nrow(summary_df_short %>% filter(matches_count >= 1L)),
            " are detected in the data, which is ",
            (
                nrow(summary_df_short %>% filter(matches_count >= 1L)) /
                n_distinct(PK_long[[metadata_info[["PriorID"]]]]) *
                100
            ),
            "%."
        )

    message(message0)
    message(message1)
    message(message2)

    if (nrow(summary_df_short %>% filter(ActionRequired == "Check")) >= 1L) {
        warning1 <-  # "Check"
            paste0(
                "There are cases where multiple detected IDs match to ",
                "multiple ",
                "prior knowledge ",
                "IDs of the same category"
            )
        log_warn(warning1)
        warning(warning1)
    }

    # # ------------------ Plot summary ----------------------##
    # x = "Class" and y = Frequency. Match Status can be colour of if no class
    # provided class = Match status.
    # Check Biocrates code.


    # # ------------------ Save Results ----------------------##
    ResList <-
        list(
            "data_summary" = summary_df_short,
            "GroupingVariable_summary" = summary_df,
            "data_long" = merged_df
        )


    save_res(
        inputlist_df = ResList,
        inputlist_plot = NULL,
        save_table = save_table,
        save_plot = NULL,
        path = Subfolder,
        file_name = "CheckMatchID_Detected-to-PK",
        core = FALSE,
        print_plot = FALSE
    )


    # Return
    invisible(return(ResList))
}

##
## Cluster prior knowledge by pathway overlap
##

#' Cluster terms in prior knowledge by set overlap
#'
#' @param data Long data frame with one ID per row, or enrichment-style table
#'     with a delimited metabolite list per term (see input_format).
#' @param metadata_info List with entries `metabolite_column` (metabolite ID
#'     column or delimited list column) and `pathway_column` (term column).
#'     Defaults to `c(metabolite_column = "MetaboliteID", pathway_column = "term")`.
#' @param similarity Similarity measure between term ID sets. Options:
#'     "jaccard" (default), "overlap_coefficient", or "correlation".
#'     Jaccard similarity is |A ∩ B| / |A ∪ B|. Overlap coefficient is
#'     |A ∩ B| / min(|A|, |B|). Jaccard is stricter for large sets, while
#'     overlap_coefficient is more permissive for nested sets.
#' @param correlation_method Correlation method when `similarity = "correlation"`.
#'     One of "pearson", "spearman", "kendall". Ignored otherwise.
#' @param input_format Input format of `data`. Use "long" for one ID per row
#'     (default) or "enrichment" for one term per row with a delimited ID list.
#'     The `metabolite_column` entry in metadata_info is interpreted accordingly.
#' @param delimiter Delimiter for metabolite ID lists when input_format =
#'     "enrichment". Ignored for input_format = "long". Default = "/".
#' @param threshold Similarity cutoff for keeping edges (applies to all
#'     clustering modes). Default = 0.5.
#' @param clust Clustering strategy: "components" (connected components on
#'     thresholded unweighted graph), "community" (Louvain on thresholded
#'     weighted graph), or "hierarchical" (hclust on distance = 1 - similarity).
#' @param hclust_method Linkage method for hierarchical clustering. One of
#'     "average" (default), "single", "complete", "ward.D", "ward.D2",
#'     "mcquitty", "median", "centroid". Used only when clust = "hierarchical".
#' @param min Minimum cluster size; smaller clusters are relabeled to "None".
#'     Default = 2.
#' @param plot_name \emph{Optional: } String added to output files of the plot.
#'     Default = "ClusterGraph".
#' @param max_nodes \emph{Optional: } Maximum nodes for plotting. If set,
#'     keeps nodes from the largest component up to this limit (by degree).
#'     Used only for the graph plot. Default = 10000.
#' @param min_degree \emph{Optional: } Minimum degree filter for graph plotting.
#'     Used only for the graph plot. Default = 1.
#' @param node_size_column \emph{Optional: } Numeric column name from `data`
#'     used to scale node sizes in the graph. Aggregated per term (mean) when
#'     multiple rows map to the same term. Default = NULL.
#' @param show_density \emph{Optional: } If TRUE, add a hull background per
#'     cluster to the graph. Default = FALSE.
#' @param seed \emph{Optional: } Random seed for graph layout reproducibility.
#'     Default = NULL.
#' @param save_plot \emph{Optional: } Select the file type of output plots.
#'     Options are svg, pdf, png or NULL. \strong{Default = "svg"}
#' @param print_plot \emph{Optional: } If TRUE prints an overview of resulting
#'     plots. \strong{Default = FALSE}
#' @param path {Optional:} String which is added to the resulting folder name.
#'     \strong{default: NULL}
#'
#' @return A list with:
#'     \item{data}{Input data with a `cluster` column added.}
#'     \item{cluster_summary}{Summary of cluster sizes and percentages.}
#'     \item{clusters}{Named vector of term -> cluster assignment.}
#'     \item{similarity_matrix}{Term-by-term similarity matrix.}
#'     \item{distance_matrix}{Term-by-term distance matrix (1 - similarity).}
#'     \item{graph_plot}{Graph plot returned by viz_graph.}
#' 
#' @examples
#' 
#' # Load example data
#' d <- metsigdb_kegg()
#' 
#' # Run clustering with graph plotting
#' r <- cluster_pk(
#'     d,
#'     metadata_info = c(
#'         metabolite_column = "MetaboliteID",
#'         pathway_column = "term"
#'     ),
#'     input_format = "long",
#'     similarity = "jaccard",
#'     threshold = 0.2,
#'     clust = "community",
#'     min = 2,
#'     plot_name = "GraphExample_long_format",
#'     save_plot = "png",
#'     min_degree = 1,
#'     seed = 123,
#'     show_density = TRUE,
#'     max_nodes = 1000
#' ) 
#' 
#' print(head(r$cluster_summary))
#' print(r$graph_plot)
#' 
#' ## add an example for an enrichment format result
#'
#' @importFrom dplyr group_by summarize ungroup mutate select left_join
#' @importFrom dplyr across n distinct filter tibble arrange
#' @importFrom igraph graph_from_adjacency_matrix components cluster_louvain
#' @importFrom stats cor as.dist hclust cutree
#' @importFrom rlang sym !! .data
#' @importFrom logger log_trace log_warn
#' @export
cluster_pk <- function(
    data,
    metadata_info = c(
        metabolite_column = "MetaboliteID",
        pathway_column = "term"
    ),
    similarity = c("jaccard", "overlap_coefficient", "correlation"),
    correlation_method = "pearson",
    input_format = c("long", "enrichment"),
    delimiter = "/",
    threshold = 0.5,
    clust = c("components", "community", "hierarchical"),
    hclust_method = "average",
    min = 2,
    plot_name = "ClusterGraph",
    max_nodes = 10000,
    min_degree = 1,
    node_size_column = NULL,
    show_density = FALSE,
    seed = NULL,
    save_plot = "png",
    print_plot = FALSE,
    path = NULL
) {

    # NSE vs. R CMD check workaround
    cluster <- .data <- NULL

    # ---- Input checks ----------------------------------------------------
    similarity <- match.arg(similarity)
    input_format <- match.arg(input_format)
    clust <- match.arg(clust)

    if (clust == "hierarchical") {
        hclust_method <- match.arg(
            hclust_method,
            c(
                "average", "single", "complete",
                "ward.D", "ward.D2",
                "mcquitty", "median", "centroid"
            )
        )
    }

    check_param_cluster_pk(
        data = data,
        metadata_info = metadata_info,
        input_format = input_format,
        delimiter = delimiter,
        threshold = threshold,
        min = min,
        node_size_column = node_size_column,
        show_density = show_density,
        seed = seed
    )

    id_col <- metadata_info[["metabolite_column"]]
    term_col <- metadata_info[["pathway_column"]]
    min <- as.integer(min)

    if (similarity == "correlation") {
        correlation_method <-
            match.arg(correlation_method, c("pearson", "spearman", "kendall"))
    }

    # ---- Build term -> ID list ------------------------------------------
    if (input_format == "long") {
        # Drop duplicated term-ID rows to avoid inflating overlaps
        data <- dplyr::distinct(data, !!sym(term_col), !!sym(id_col), .keep_all = TRUE)

        # Warn if terms have no IDs
        empty_terms <-
            data %>%
            dplyr::group_by(!!sym(term_col)) %>%
            dplyr::summarize(
                n_ids = sum(!is.na(!!sym(id_col)) & !!sym(id_col) != ""),
                .groups = "drop"
            ) %>%
            dplyr::filter(n_ids == 0L)
        if (nrow(empty_terms) > 0L) {
            log_warn(
                "Terms with no IDs will be dropped: %s",
                paste(empty_terms[[term_col]], collapse = ", ")
            )
        }

        term_metabolites <-
            data %>%
            dplyr::filter(!is.na(!!sym(id_col)), !!sym(id_col) != "") %>%
            dplyr::group_by(!!sym(term_col)) %>%
            dplyr::summarize(
                MetaboliteIDs = list(unique(!!sym(id_col))),
                .groups = "drop"
            )
    } else {
        # Enrichment input: split delimited metabolite IDs per term
        term_metabolites <-
            data %>%
            dplyr::mutate(
                MetaboliteIDs = strsplit(
                    as.character(.data[[id_col]]),
                    delimiter,
                    fixed = TRUE
                )
            ) %>%
            dplyr::mutate(
                MetaboliteIDs = lapply(
                    MetaboliteIDs,
                    function(x) {
                        x <- trimws(x)
                        x[x != ""]
                    }
                )
            ) %>%
            dplyr::group_by(!!sym(term_col)) %>%
            dplyr::summarize(
                MetaboliteIDs = list(unique(unlist(MetaboliteIDs))),
                .groups = "drop"
            )

        empty_terms <- term_metabolites[lengths(term_metabolites$MetaboliteIDs) == 0L, , drop = FALSE]
        if (nrow(empty_terms) > 0L) {
            log_warn(
                "Terms with no IDs will be dropped: %s",
                paste(empty_terms[[term_col]], collapse = ", ")
            )
        }
        term_metabolites <- term_metabolites[lengths(term_metabolites$MetaboliteIDs) > 0L, , drop = FALSE]
    }

    terms <- term_metabolites[[term_col]]
    n_terms <- length(terms)
    if (n_terms == 0L) {
        stop("No terms with non-missing IDs after filtering.")
    }

    # ---- Similarity matrix ----------------------------------------------
    similarity_matrix <-
        matrix(
            1,
            nrow = n_terms,
            ncol = n_terms,
            dimnames = list(terms, terms)
        )

    if (similarity %in% c("jaccard", "overlap_coefficient")) {
        combs <- combn(seq_len(n_terms), 2L)
        for (j in seq_len(ncol(combs))) {
            i1 <- combs[1, j]
            i2 <- combs[2, j]
            set1 <- term_metabolites$MetaboliteIDs[[i1]]
            set2 <- term_metabolites$MetaboliteIDs[[i2]]
            inter <- length(intersect(set1, set2))
            if (similarity == "jaccard") {
                denom <- length(union(set1, set2))
            } else { # overlap coefficient
                denom <- min(length(set1), length(set2))
            }
            sim <- if (denom == 0) 0 else inter / denom
            t1 <- terms[i1]
            t2 <- terms[i2]
            similarity_matrix[t1, t2] <- sim
            similarity_matrix[t2, t1] <- sim
        }
    } else { # correlation
        metabolites <- unique(unlist(term_metabolites$MetaboliteIDs))
        binary_matrix <-
            matrix(
                0,
                nrow = n_terms,
                ncol = length(metabolites),
                dimnames = list(terms, metabolites)
            )
        for (i in seq_len(n_terms)) {
            ids <- term_metabolites$MetaboliteIDs[[i]]
            binary_matrix[i, colnames(binary_matrix) %in% ids] <- 1
        }
        corr <- cor(t(binary_matrix), method = correlation_method)
        corr[corr < 0] <- 0
        diag(corr) <- 1
        similarity_matrix <- corr
    }

    # Distance matrix
    distance_matrix <- 1 - similarity_matrix
    diag(distance_matrix) <- 0

    # Thresholded versions
    similarity_thr <- similarity_matrix
    similarity_thr[similarity_thr < threshold] <- 0
    diag(similarity_thr) <- 0

    distance_thr <- distance_matrix
    distance_thr[similarity_thr == 0] <- 1
    diag(distance_thr) <- 0

    # ---- Clustering ------------------------------------------------------
    clusters <- rep(NA_integer_, n_terms)
    names(clusters) <- terms

    if (clust == "components") {
        adj <- similarity_thr
        adj[adj > 0] <- 1
        g <- igraph::graph_from_adjacency_matrix(
            adj,
            mode = "undirected",
            weighted = NULL
        )
        mem <- igraph::components(g)$membership
        if (is.null(names(mem))) {
            names(mem) <- igraph::V(g)$name
        }
        clusters[names(mem)] <- mem
    } else if (clust == "community") {
        adj <- similarity_thr
        g <- igraph::graph_from_adjacency_matrix(
            adj,
            mode = "undirected",
            weighted = TRUE
        )
        mem <- igraph::cluster_louvain(g)$membership
        if (is.null(names(mem))) {
            names(mem) <- igraph::V(g)$name
        }
        clusters[names(mem)] <- mem
    } else if (clust == "hierarchical") {
        hc <- hclust(as.dist(distance_thr), method = hclust_method)
        cut_height <- 1 - threshold
        mem <- cutree(hc, h = cut_height)
        if (is.null(names(mem))) {
            names(mem) <- terms
        }
        clusters[names(mem)] <- mem
    }

    # ---- Apply minimum size filter --------------------------------------
    cluster_sizes <- table(clusters, useNA = "no")
    small <- names(cluster_sizes[cluster_sizes < min])
    clusters[as.character(clusters) %in% small] <- NA_integer_

    # Label clusters
    cluster_labels <- ifelse(
        is.na(clusters),
        "None",
        paste0("cluster", clusters)
    )
    names(cluster_labels) <- terms

    # Merge back to data
    term_metabolites$cluster <- cluster_labels[term_metabolites[[term_col]]]
    df <-
        data %>%
        dplyr::left_join(
            term_metabolites %>% dplyr::select(!!sym(term_col), cluster),
            by = term_col
        )

    # ---- Summary ---------------------------------------------------------
    cluster_summary <-
        term_metabolites %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarize(n_terms = dplyr::n(), .groups = "drop") %>%
        dplyr::mutate(pct_terms = 100 * n_terms / sum(n_terms))

    # ---- Node sizes ------------------------------------------------------
    node_sizes <- NULL
    if (!is.null(node_size_column)) {
        node_size_df <-
            data %>%
            dplyr::group_by(!!sym(term_col)) %>%
            dplyr::summarize(
                node_size = mean(.data[[node_size_column]], na.rm = TRUE),
                .groups = "drop"
            )
        node_size_df$node_size[is.nan(node_size_df$node_size)] <- NA_real_
        node_sizes <- node_size_df$node_size
        names(node_sizes) <- node_size_df[[term_col]]
        node_sizes <- node_sizes[terms]
        attr(node_sizes, "label") <- node_size_column
    }

    # ---- Graph plot ------------------------------------------------------
    graph_plot <- viz_graph(
        similarity_matrix = similarity_matrix,
        clusters = cluster_labels,
        threshold = threshold,
        plot_name = plot_name,
        max_nodes = max_nodes,
        min_degree = min_degree,
        node_sizes = node_sizes,
        show_density = show_density,
        seed = seed,
        save_plot = save_plot,
        print_plot = print_plot,
        path = path
    )

    out <- list(
        data = df,
        cluster_summary = cluster_summary,
        clusters = cluster_labels,
        similarity_matrix = similarity_matrix,
        distance_matrix = distance_matrix,
        graph_plot = graph_plot
    )

    log_trace(paste0(
        "cluster_pk completed with ",
        clust,
        " clustering; ",
        length(unique(cluster_labels)),
        " clusters (including None)."
    ))

    return(out)
}



#
# Helper function to add term information to Enrichment Results
#

# Better function Name and parameter names needed
# Use in ORA functions and showcase in vignette with decoupleR output
#' Adds extra columns to enrichment output
#'
#' These columns inform about 1. The amount of genes associated with term in
#' prior knowledge, 2. The amount of genes detected in input data associated
#' with term in prior knowledge, and 3. The percentage of genes detected in
#' input data associated with term in prior knowledge.
#'
#' @param mat data matrix used as input for enrichment analysis
#' @param net Prior Knowledge used as input for enrichment analysis
#' @param res Results returned from the enrichment analysis
#' @param .source used as input for enrichment analysis
#' @param .target used as input for enrichment analysis
#' @param complete TRUE or FALSE, weather only .source with results should be returned or
#'     all .source in net.
#'
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr desc
#' @noRd
add_info <- function(
    mat,
    net,
    res,
    .source,
    .target,
    complete = FALSE
) {

    # NSE vs. R CMD check workaround
    targets_detected_num <- score <- NULL

    # add number of Genes_targeted_by_TF_num
    net$Count <- 1
    net_Mean <-
        aggregate(
            net$Count,
            by = list(source = net[[.source]]),
            FUN = sum
        ) %>%
        rename("targets_num" = 2)

    if (complete) {
        res_Add <-
            merge(
                x = res,
                y = net_Mean,
                by = "source",
                all = TRUE
            )
    } else {
        res_Add <-
            merge(
                x = res,
                y = net_Mean,
                by = "source",
                all.x = TRUE
            )
    }

    # add list of Genes_targeted_by_TF_chr
    net_List <-
        aggregate(
            net[[.target]] ~ net[[.source]],
            FUN = toString
        ) %>%
        rename(
            "source" = 1,
            "targets_chr" = 2
        )
    res_Add <-
        merge(
            x = res_Add,
            y = net_List,
            by = "source",
            all.x = TRUE
        )

    # add number of Genes_targeted_by_TF_detected_num
    mat <-
        as.data.frame(mat) %>%  # Are these the normalised counts?
        rownames_to_column("Symbol")

    Detected <-
        merge(
            x = mat,
            y = net[, c(.source, .target)],
            by.x = "Symbol",
            by.y = .target,
            all.x = TRUE
        ) %>%
        filter(
            !is.na(
                across(
                    all_of(.source)
                )
            )
        )
    Detected$Count <- 1
    Detected_Mean <-
        aggregate(
            Detected$Count,
            by = list(source = Detected[[.source]]),
            FUN = sum
        ) %>%
        rename("targets_detected_num" = 2)

    res_Add <-
        merge(
            x = res_Add,
            y = Detected_Mean,
            by = "source",
            all.x = TRUE
        ) %>%
        mutate(
            targets_detected_num = replace_na(targets_detected_num, 0)
        )

    # add list of Genes_targeted_by_TF_detected_chr
    Detected_List <-
        aggregate(
            Detected$Symbol ~ Detected[[.source]],
            FUN = toString
        ) %>%
        rename(
            "source" = 1,
            "targets_detected_chr" = 2
        )

    res_Add <-
        merge(
            x = res_Add,
            y = Detected_List,
            by = "source",
            all.x = TRUE
        )

    # add percentage of percentage_of_Genes_detected
    res_Add$targets_detected_percentage <-
        round(
            ((res_Add$targets_detected_num / res_Add$targets_num) * 100),
            digits = 2
        )

    # sort by score
    res_Add <-
        res_Add %>%
        arrange(
            desc(
                as.numeric(
                    as.character(score)
                )
            )
        )

    # # ------------------ Save and return ------------------- ##
    Output <- res_Add
}


#
# Compare Prior Knowledge resources against each other or themselves
#

#' Compare Prior Knowledge Resources and/or Columns within a Single Resource
#'
#' and Generate an UpSet Plot This function compares gene and/or metabolite
#' features across multiple prior knowledge (PK) resources or, if a single
#' resource is provided with a vector of column names in \code{metadata_info},
#' compares columns within that resource. In the multi-resource mode, each
#' element in \code{data} represents a PK resource (either as a data frame or a
#' recognized resource name) from which a set of features is extracted. A
#' binary summary table is then constructed and used to create an UpSet plot.
#' In the within-resource mode, a single data frame is provided (with
#' \code{data} containing one element) and its \code{metadata_info} entry is a
#' vector of column names to compare (e.g., binary indicators for different
#' annotations). In this case, the function expects the data frame to have a
#' grouping column named \code{"Class"} (or, alternatively, a column specified
#' via the \code{class_col} attribute in \code{metadata_info}) that is used for
#' grouping in the UpSet plot.
#'
#' @param data A named list where each element corresponds to a prior knowledge (PK)
#'     resource. Each element can be: \itemize{ \item   A data frame containing
#'     gene/metabolite identifiers (and additional columns for within-resource
#'     comparison), \item   A character string indicating the resource name.
#'     Recognized names include (but are not limited to): \code{"Hallmarks"},
#'     \code{"Gaude"}, \code{"MetalinksDB"}, and \code{"RAMP"} (or
#'     \code{"metsigdb_chemicalclass"}). In the latter case, the function will
#'     attempt to load the corresponding data automatically. }
#' @param metadata_info A named list (with names matching those in \code{data}) where each
#'     element is either a character string or a character vector indicating
#'     the column name(s) to extract features. For multiple-resource
#'     comparisons, these refer to the columns containing feature identifiers.
#'     For within-resource comparisons, the vector should list the columns to
#'     compare (e.g., \code{c("CHEBI", "HMDB", "LIMID")}). In within-resource
#'     mode, the input data frame is expected to contain a column named
#'     \code{"Class"} (or a grouping column specified via the \code{class_col}
#'     attribute). \emph{If no grouping column is found, a default grouping
#'     column named \code{"Group"} (with all rows assigned the same value) is
#'     created.}
#' @param filter_by Character. Optional filter for the resulting features when comparing
#'     multiple resources. Options are: \code{"both"} (default), \code{"gene"},
#'     or \code{"metabolite"}. This parameter is ignored in within-resource
#'     mode.
#' @param plot_name \emph{Optional: } String which is added to the output files of the
#'     Upsetplot \strong{Default = ""}
#' @param name_col \emph{Optional: } column name including the feature names. Default is
#'     \code{"TrivialName"}.
#' @param palette_type Character. Color palette to be used in the plot. Default is
#'     \code{"polychrome"}.
#' @param save_plot \emph{Optional: } Select the file type of output plots. Options are svg,
#'     png, pdf. \strong{Default = svg}
#' @param save_table \emph{Optional: } File types for the analysis results are: "csv",
#'     "xlsx", "txt". \strong{Default = "csv"}
#' @param print_plot \emph{Optional: } TRUE or FALSE, if TRUE Volcano plot is saved as an
#'     overview of the results. \strong{Default = TRUE}
#' @param path \emph{Optional:} Path to the folder the results should be saved at.
#'     \strong{Default = NULL}
#'
#' @return A list containing two elements: \itemize{ \item{summary_table: A data
#'     frame representing either: \itemize{ \item{the binary summary matrix of
#'     feature presence/absence across multiple resources, or} \item{the
#'     original data frame (augmented with binary columns and a \code{None}
#'     column) in within-resource mode.} } } \item{upset_plot: The UpSet plot
#'     object generated by the function.} }
#'
#' @examples
#' \dontrun{
#' ## Example 1: Within-Resource Comparison
#' ## (Comparing Columns Within a Single data Frame)
#'
#' # biocrates_features is a data frame with columns:
#' # "TrivialName", "CHEBI", "HMDB", "LIMID", and "Class".
#' # Here the "Class" column is used as the grouping variable
#' # in the UpSet plot.
#' data(biocrates_features)
#' data_single <- list(Biocft = biocrates_features)
#' metadata_info_single <- list(Biocft = c("CHEBI", "HMDB", "LIMID"))
#'
#' res_single <-
#' compare_pk(
#' data = data_single, metadata_info = metadata_info_single,
#' plot_name = "Overlap of BioCrates Columns"
#' )
#'
#' ## Example 2: Custom data Frames with Custom Column Names
#'
#' # Example with preloaded data frames and custom column names:
#' hallmarks_df <- data.frame(
#' feature = c("HMDB0001", "GENE1", "GENE2"),
#' stringsAsFactors = FALSE
#' )
#' gaude_df <- data.frame(
#' feature = c("GENE2", "GENE3"),
#' stringsAsFactors = FALSE
#' )
#' metalinks_df <- data.frame(
#' hmdb = c("HMDB0001", "HMDB0002"),
#' gene_symbol = c("GENE1", "GENE4"),
#' stringsAsFactors = FALSE
#' )
#' ramp_df <- data.frame(
#' class_source_id = c("HMDB0001", "HMDB0003"),
#' stringsAsFactors = FALSE
#' )
#' data <- list(
#' Hallmarks = hallmarks_df, Gaude = gaude_df,
#' MetalinksDB = metalinks_df, RAMP = ramp_df
#' )
#' metadata_info <- list(
#' Hallmarks = "feature", Gaude = "feature",
#' MetalinksDB = c("hmdb", "gene_symbol"),
#' RAMP = "class_source_id"
#' )
#' res <- compare_pk(
#' data = data, metadata_info = metadata_info,
#' filter_by = "metabolite"
#' )
#' }
#'
#' @importFrom dplyr mutate select
#' @importFrom logger log_trace
#' @export
compare_pk <- function(
    data,
    metadata_info = NULL,
    filter_by = c("both", "gene", "metabolite"),
    plot_name = "Overlap of Prior Knowledge Resources",
    name_col = "TrivialName",
    palette_type = "polychrome",
    save_plot = "svg",
    save_table = "csv",
    print_plot = TRUE,
    path = NULL
) {

    # NSE vs. R CMD check workaround
    Type <- NULL

    # # ------------ Create log file ----------- ##
    metaproviz_init()

    # # ------------ Check Input files ----------- ##
    # Match filter argument
    filter_by <-
        match.arg(filter_by)

    # Validate data input
    if (!is.list(data) || length(data) < 1L) {
        message <-
            paste0("data must be a non - empty list.")
        log_trace(
            paste("Error ", message, sep = "")
        )
        stop(message)
    }
    if (is.null(names(data)) || any(names(data) == "")) {
        message <-
            paste0("data must be a named list with resource names.")
        log_trace(
            paste("Error ", message, sep = "")
        )
        stop(message)
    }

    # # ------------ Create folders ----------- ##
    if (!is.null(save_plot) | !is.null(save_table)) {
        folder <-
            save_path(
                folder_name = "PK",
                path = path
            )

        subfolder <- file.path(folder, "ComparePK")
        if (!dir.exists(subfolder)) {
            dir.create(subfolder)
        }
    }

    # ##########################################################################
    # # ----------- Input ----------- ##
    # Define resource lookup table with information on how to retrieve and
    # transform each resource.
    default_cols <-
        list(
            hallmarks = "feature",
            gaude = "feature",
            metalinksdb = c("hmdb", "gene_symbol"),
            ramp = "class_source_id"
        )

    # Initialize metadata_info if not provided.
    if (is.null(metadata_info)) {
        metadata_info <- list()
    }

    # Determine if we are in within-resource mode.
    # If only one resource is provided and its metadata_info entry has >1
    # column, assume within-resource comparison.
    single_resource <- (length(data) == 1L)
    within_resource_mode <- FALSE
    if (single_resource) {
        resource_name <- names(data)[1]
        if (!is.null(metadata_info[[resource_name]]) &&
        length(metadata_info[[resource_name]]) > 1L) {
            within_resource_mode <- TRUE
        }
    }

    if (within_resource_mode) {
        # ===== Within-Resource Comparison Mode =====
        # Retrieve the single data frame.
        resource_data <- data[[resource_name]]

        # Identify the intersection columns based on metadata_info.
        intersect_cols <- metadata_info[[resource_name]]
        missing_cols <-
            setdiff(
                intersect_cols,
                colnames(resource_data)
            )
        if (length(missing_cols) > 0L) {
            stop(
                "The following intersection column(s) specified in",
                "metadata_info were not found in resource '",
                resource_name, "': ", paste(missing_cols, collapse = ", ")
            )
        }

        # Identify a column for grouping. If none exists, create a default
        # grouping column.
        if ("Class" %in% colnames(resource_data)) {
            class_col <- "Class"
        } else if (!is.null(attr(metadata_info[[resource_name]], "class_col"))) {
            class_col <- attr(metadata_info[[resource_name]], "class_col")
        } else {
            # No grouping column provided—create a default column named "Group"
            # with the same value for all rows.
            resource_data$Group <- "All"
            class_col <- "Group"
        }

        # Convert the specified intersection columns to binary (0/1). Here
        # non-NA and values != 0 are treated as present.
        binary_suffix <- "_bin"
        for (col in intersect_cols) {
            new_col <- paste0(col, binary_suffix)
            resource_data[[new_col]] <-
                as.integer(
                    !is.na(resource_data[[col]]) &
                    (resource_data[[col]] != 0L) &
                        (resource_data[[col]] != "")
                )
        }

        # Identify the binary columns based on the suffix
        bin_cols <-
            grep(
                paste0(
                    binary_suffix,
                    "$"
                ),
                colnames(resource_data),
                value = TRUE
            )

        # Create the "None" column using the binary columns
        resource_data$None <-
            as.integer(
                rowSums(
                    resource_data[, bin_cols, drop = FALSE]
                ) == 0L
            )

        # Create df_summary, potentially with the name column (default is
        # TrivialName, if not found, it will not be included-only used to help
        # interpret summary table)
        if (name_col %in% colnames(resource_data)) {
            summary_cols <-
                c(
                    name_col,
                    bin_cols,
                    "None",
                    class_col
                )
        } else {
            summary_cols <-
                c(
                    bin_cols,
                    "None",
                    class_col
                )
        }
        df_summary <-
            resource_data[, summary_cols, drop = FALSE]
        # Rename the cols again
        # Find indices of columns ending in "_bin"
        bin_cols_idx <-
            grep(
                "_bin$",
                names(df_summary)
            )
        # Remove the suffix from these column names
        names(df_summary)[bin_cols_idx] <-
            sub(
                "_bin$",
                "",
                names(df_summary)[bin_cols_idx]
            )

        # Generate the UpSet plot.
        upset_plot <-
            tryCatch(
                viz_upset(
                    df = df_summary,
                    class_col = class_col,
                    intersect_cols = c(intersect_cols, "None"),
                    plot_name = plot_name,
                    palette_type = palette_type,
                    save_plot = NULL,
                    print_plot = FALSE
                ),
                error = function(e) {
                    log_warn(
                        paste0(
                            "viz_upset failed inside compare_pk(): ",
                            conditionMessage(e)
                        )
                    )
                    NULL   ## viz_upset failed -> set upset_plot variable to NULL
                }
            )

        summary_table <- df_summary
    } else {
        # ===== Multi-Resource Comparison Mode =====
        # Process each resource in data.
        for (res in names(data)) {
            resource_val <- data[[res]]
            resource_id <- tolower(res)
            if (
                is.null(metadata_info[[res]]) &&
                resource_id %in% names(default_cols)
            ) {
                metadata_info[[res]] <-
                    default_cols[[resource_id]]
            } else if (is.null(metadata_info[[res]])) {
                stop(
                    "metadata_info must be provided for resource: ",
                    res
                )
            }
        }

        # Extract features from each resource based on metadata_info.
        resource_features <- list()
        for (res in names(data)) {
            resource_data <- data[[res]]
            cols <- metadata_info[[res]]
            if (!all(cols %in% colnames(resource_data))) {
                stop(
                    paste(
                        "Column(s)",
                        paste(cols, collapse = ", "),
                        "not found in resource",
                        res
                    )
                )
            }
            features <-
                if (length(cols) > 1L) {
                    unique(
                        unlist(
                            lapply(
                                cols,
                                function(col) na.omit(resource_data[[col]])
                            )
                        )
                    )
                } else {
                    unique(
                        na.omit(resource_data[[cols]])
                    )
                }
            resource_features[[res]] <-
                as.character(features)
        }

        # Compile all unique features across resources.
        all_features <-
            unique(
                unlist(resource_features)
            )

        # Create the binary summary table.
        df_binary <-
            data.frame(
                Feature = all_features,
                stringsAsFactors = FALSE
            )

        for (res in names(resource_features)) {
            df_binary[[res]] <-
                as.integer(
                    all_features %in% resource_features[[res]]
                )
        }

        df_binary$Type <-
            ifelse(
                grepl(
                    "^HMDB",
                    df_binary$Feature
                ),
                "metabolite (HMDB)",
                "gene"
            )

        resource_cols <-
            names(resource_features)

        df_binary$None <-
            as.integer(
                rowSums(df_binary[, resource_cols, drop = FALSE]) == 0L
            )

        # Optionally filter the summary table.
        if (filter_by == "gene") {
            df_binary <-
                subset(
                    df_binary,
                    Type == "gene"
            )
        } else if (filter_by == "metabolite") {
            df_binary <-
                subset(
                    df_binary,
                    Type == "metabolite (HMDB)"
                )
        }

        # Generate the UpSet plot.
        upset_plot <-
            tryCatch(
                viz_upset(
                    df = df_binary,
                    class_col = "Type",
                    intersect_cols = resource_cols,
                    plot_name = plot_name,
                    palette_type = palette_type,
                    save_plot = NULL,
                    print_plot = FALSE
                ),
                error = function(e) {
                    log_warn(
                        paste0(
                            "viz_upset failed inside compare_pk(): ",
                            conditionMessage(e)
                        )
                    )
                    NULL     ## viz_upset failed -> set upset_plot variable to NULL
                }
            )

        summary_table <- df_binary
    }

    # ##########################################################################
    # ## -------- Save and return ----------###


    save_res(
        inputlist_df = list(summary_table = summary_table),
        inputlist_plot = list(upset_plot = upset_plot),
        save_table = save_table,
        save_plot = save_plot,
        path = subfolder,
        file_name = "compare_pk",
        core = FALSE,
        print_plot = print_plot
    )


    return(
        list(
            summary_table = summary_table,
            upset_plot = upset_plot
        )
    )
}


################################################################################
### ### Helper function to count number of entries for an ID column value    ###
### ### and plot                                                             ###
#' Count Entries and Generate a Histogram Plot for a Specified Column
#'
#' This function processes a data frame column by counting the number of
#' entries within each cell. It considers both \code{NA} values and empty
#' strings as zero entries, and categorizes each cell as "No ID", "Single ID",
#' or "Multiple IDs" based on the count. A histogram is then generated to
#' visualize the distribution of entry counts. scale_x_continuous
#'
#' @param data A data frame containing the data to be analyzed.
#' @param column A string specifying the name of the column in \code{data} to analyze.
#' @param delimiter A string specifying the delimiter used to split cell values. Defaults to
#'     \code{","}.
#' @param fill_colors A named character vector providing colors for each category. Defaults to
#'     \code{c("No ID" = "#FB8072", "Single ID" = "#B3DE69", "Multiple IDs" =
#'     "#80B1D3")}.
#' @param binwidth Numeric value specifying the bin width for the histogram. Defaults to
#'     \code{1}.
#' @param title_prefix A string to use as the title of the plot. If \code{NULL} (default), the
#'     title will be generated as "Number of <column> IDs per Biocrates Cell".
#' @param save_plot \emph{Optional: } Select the file type of output plots. Options are svg,
#'     png, pdf. \strong{Default = svg}
#' @param save_table \emph{Optional: } File types for the analysis results are: "csv",
#'     "xlsx", "txt". \strong{Default = "csv"}
#' @param print_plot \emph{Optional: } TRUE or FALSE, if TRUE Volcano plot is saved as an
#'     overview of the results. \strong{Default = TRUE}
#' @param path \emph{Optional:} Path to the folder the results should be saved at.
#'     \strong{Default = NULL}
#'
#' @return A list with two elements: \item{result}{A data frame that includes three
#'     additional columns: \code{was_na} (logical indicator of missing or empty
#'     cells), \code{entry_count} (number of entries in each cell), and
#'     \code{id_label} (a categorical label based on the entry count).}
#'     \item{plot}{A \code{ggplot} object representing the histogram of entry
#'     counts.}
#'
#' @examples
#' data(biocrates_features)
#' count_id(biocrates_features, "HMDB")
#'
#' @importFrom dplyr case_when mutate rename
#' @importFrom ggplot2 aes element_rect element_text expansion geom_bar
#' @importFrom ggplot2 geom_histogram ggplot labs scale_fill_manual
#' @importFrom ggplot2 scale_y_continuous theme theme_classic
#' @importFrom grid convertUnit
#' @export
count_id <- function(
    data,
    column,
    # *sobbing*
    delimiter = ",
    ",
    fill_colors = c(
        "No ID" = "#FB8072",
        "Single ID" = "#B3DE69",
        "Multiple IDs" = "#80B1D3"
    ),
    binwidth = 1,
    title_prefix = NULL,
    save_plot = "svg",
    save_table = "csv",
    print_plot = TRUE,
    path = NULL
) {
    # @Macabe:
    # move named fill colors inside of the function. if the user provides other
    # string of colors, we change them.
    # we need to specify that NA is counted as none
    # we need to check for duplications (i.e. is the trivialname duplicated in
    # the data frame, remove this and give a warning)
    # give the user the change to pass multiple columns to analyse, which would
    # mean we create a plot for each column and label the plot with the column
    # name
    # add count_id function into Equivalent IDs: make a plot before and after
    # equivalent IDs and put them side by side to return as the QC plot of the
    # function
    # add save_plot and save_table
    # return both standard plot and Plot_Sized (=nice version of the plot which
    # is saved. the other plot is still returned, since this is the ggplot
    # version and can be changed further)
    # create subtitle and not title prefix

    # NSE vs. R CMD check workaround
    .data <- entry_count <- id_label <- NULL

    # # ------------------  Create output folders and path ------------------- ##
    if (!is.null(save_table)) {
        folder <- save_path(
            folder_name = "PK",
            path = path
        )

        Subfolder <- file.path(folder, "CountIDs")
        if (!dir.exists(Subfolder)) {
            dir.create(Subfolder)
        }
    }

    # # ------------------  data table ------------------- ##
    # Process the data: count entries and label each cell based on the number of
    # entries.
    processed_data <-
        mutate(
            data,
            was_na = is.na(.data[[column]]) | .data[[column]] == "",
            entry_count = map_int(
                .data[[column]],
                function(cell) {
                    if (is.na(cell) || cell == "") {
                        0L
                    # Treat NA or empty as 0 entries for counting
                    } else {
                        as.integer(length(unlist(strsplit(
                            as.character(cell),
                            delimiter
                        ))))
                    }
                }
            ),
            id_label = case_when(
                entry_count == 0L ~ "No ID",
                entry_count == 1L ~ "Single ID",
                entry_count >= 2L ~ "Multiple IDs"
            )
        )


    # # ------------------  plot ------------------- ##
    # Generate the plot title if not provided
    if (is.null(title_prefix)) {
        plot_name <- paste("Number of", column, "IDs per measured peak.")
    } else {
        plot_name <- title_prefix
    }

    # Create the histogram plot using ggplot2
    if (length(unique(processed_data$entry_count)) > 1L) {
        plot_obj <-
            ggplot(
                processed_data,
                aes(
                    x = entry_count,
                    fill = id_label
                )
            ) +
            geom_histogram(
                binwidth = binwidth,
                boundary = -0.5,
                color = "black"
            ) +
            # Set axis to start at 0
            scale_x_continuous(
                expand = expansion(mult = c(0, 0.05))
            ) +
            scale_y_continuous(
                expand = expansion(mult = c(0, 0.05))
            )
    } else {
        # Create bargraph plot using ggplot2
        plot_obj <-
            ggplot(
                processed_data,
                aes(
                    x = as.factor(entry_count),
                    fill = id_label
                )
            ) +
            geom_bar(
                color = "black"
            ) +
            # Set axis to start at 0
            scale_y_continuous(
                expand = expansion(mult = c(0, 0.05))
            )
    }

    plot_obj <-
        plot_obj +
        scale_fill_manual(
            values = fill_colors
        ) +
        labs(
            title = plot_name,
            x = "Number of IDs",
            y = "Frequency",
            fill = paste0(column, " IDs", sep = "")
        ) +
        theme_classic() +
        theme(
            plot.title = element_text(hjust = 0.5, size = 12),
            legend.position.inside = c(0.8, 0.8),
            legend.justification = c("right", "top"),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 12)
        )

    # Make the nice plot:
    Plot_Sized <-
        plot_grob_superplot(
            input_plot = plot_obj,
            metadata_info = c(
                Conditions = "id_label",
                Superplot = TRUE
            ),
            metadata_sample = processed_data %>%
            dplyr::rename("Conditions" = "entry_count"),
            plot_name = plot_name,
            subtitle = "",
            plot_type = "Bar"
        )

    plot_height <-
        convertUnit(
            Plot_Sized$height,
            "cm",
            valueOnly = TRUE
        )

    plot_width <-
        convertUnit(
            Plot_Sized$width,
            "cm",
            valueOnly = TRUE
        )

    Plot_Sized %<>%
        {
            ggplot() +
            annotation_custom(.)
        } %>%
        add(
            theme(
                panel.background = element_rect(fill = "transparent")
            )
        )

    # # ------------------  save and return ------------------- ##


    save_res(
        #  inputlist_df needs to be a list, also for single comparisons
        inputlist_df = list("Table" = processed_data),
        inputlist_plot = list("Plot_Sized" = Plot_Sized),
        save_table = save_table,
        save_plot = save_plot,
        path = Subfolder,
        file_name = "Count_MetaboliteIDs",
        core = FALSE,
        print_plot = print_plot
    )


    OutputList <- list()
    OutputList <-
        list(
            "Table" = processed_data,
            "Plot" = plot_obj,
            "Plot_Sized" = Plot_Sized
        )

    # Return the processed data and the plot object as a list
    return(
        invisible(OutputList)
    )
}
