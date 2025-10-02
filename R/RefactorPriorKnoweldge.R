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


################################################################################
### ### ### Translate IDs to/from KEGG, PubChem, Chebi, HMDB ### ### ###
################################################################################

#' Translate IDs to/from KEGG, PubChem, Chebi, HMDB
#'
#' @param data dataframe with at least one column with the target (e.g.
#' metabolite), you can add other columns such as source (e.g. term). Must be
#' "long" DF, meaning one ID per row.
#' @param metadata_info \emph{Optional: } Column name of Target in input_pk.
#' \strong{Default = list(InputID="MetaboliteID" , grouping_variable="term")}
#' @param from ID type that is present in your data. Choose between "kegg",
#' "pubchem", "chebi", "hmdb". \strong{Default = "kegg"}
#' @param to One or multiple ID types to which you want to translate your data.
#' Choose between "kegg", "pubchem", "chebi", "hmdb". \strong{Default =
#' c("pubchem","chebi","hmdb")}
#' @param summary \emph{Optional: } If TRUE a long summary tables are created.
#' \strong{Default = FALSE}
#' @param save_table \emph{Optional: } File types for the analysis results are:
#' "csv", "xlsx", "txt". \strong{Default = "csv"}
#' @param path {Optional:} Path to the folder the results should be saved at.
#' \strong{Default = NULL}
#'
#' @return List with at least three DFs: 1) Original data and the new column of
#' translated ids spearated by comma. 2) Mapping information between Original
#' ID
#' to Translated ID. 3) Mapping summary between Original ID to Translated ID.
#'
#' @examples
#' KEGG_Pathways <- metsigdb_kegg()
#' Res <- translate_id(data = KEGG_Pathways, metadata_info = c(InputID =
#' "MetaboliteID", grouping_variable = "term"), from = c("kegg"), to =
#' c("pubchem", "hmdb"))
#'
#' @keywords Translate metabolite IDs
#'
#' @importFrom dplyr mutate select group_by ungroup distinct filter across n
#' @importFrom tidyselect all_of starts_with
#' @importFrom rlang !!! !! := sym syms
#' @importFrom OmnipathR id_types translate_ids
#' @importFrom logger log_warn log_trace
#' @importFrom stringr str_to_lower
#' @export
translate_id <- 
    function(
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

  ## ------------------  Check Input ------------------- ##
  # HelperFunction `check_param`
    check_param(
        data = data,
        data_num = FALSE,
        save_table = save_table
    )

  # Specific checks:
    if ("InputID" %in% names(metadata_info)) {
        if (metadata_info[["InputID"]] %in% colnames(data) == FALSE) {
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
        if (metadata_info[["grouping_variable"]] %in% colnames(data) == FALSE) {
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

    if (is.logical(summary) == FALSE) {
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
  doublons <- data %>%
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
            n() > 1
        ) %>%
    summarize()

    if (nrow(doublons) > 0) {
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

  ## ------------------  Create output folders and path ------------------- ##
    if (is.null(save_table) == FALSE) {
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

    ############################################################################
  ## ------------------ Translate to-from for each pair ------------------- ##
    TranslatedDF <- 
        translate_ids(
      data,
            !!sym(metadata_info[["InputID"]]) := !!sym(from),
            !!!syms(to), # list of symbols, hence three !!!
      ramp = TRUE,
      expand = FALSE,
      quantify_ambiguity = TRUE,
      qualify_ambiguity = TRUE,
            ambiguity_groups = metadata_info[["grouping_variable"]], # Checks within the groups, without it checks across groups
      ambiguity_summary = TRUE
    )
    # TranslatedDF %>% attributes %>% names
    # TranslatedDF%>% attr('ambiguity_MetaboliteID_hmdb')

  ## --------------- Create output DF -------------------- ##
  ResList <- list()

    ## Create DF for TranslatedIDs only with the original data and the
    ## translatedID columns
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
            collapse = ", "),
            .names = "{.col}")
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

  ## Add DF with mapping information
  ResList[["TranslatedDF_MappingInfo"]] <- TranslatedDF

  ## Also save the different mapping summaries!
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

  ## Create the long DF summary if summary =TRUE
    if (summary == TRUE) {
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


  ## ------------------ Save the results ------------------- ##
    suppressMessages(
        suppressWarnings(
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
        )
    )

  ## ------------------ Show Upset plot of the results ------------------- ##
  ## toDO: current issue with Lists vs doubles
  # if (plot==TRUE) {
  #   pk_list <- list(translated = TranslatedDF)
  #   metadata_info <- list(translated = to)
  #   print(pk_list)
  #   print(to)
    #   pk_comp_res <- 
    #       MetaProViz:::compare_pk(
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



################################################################################
### ### ###             Find additional potential IDs                ### ### ###
################################################################################

#' Find additional potential IDs for  "kegg", "pubchem", "chebi", "hmdb"
#'
#' @param data dataframe with at least one column with the detected metabolite
#' IDs (one ID per row).
#' @param metadata_info \emph{Optional: } Column name of metabolite IDs.
#' \strong{Default = list(InputID="MetaboliteID")}
#' @param from ID type that is present in your data. Choose between "kegg",
#' "pubchem", "chebi", "hmdb". \strong{Default = "hmdb"}
#' @param save_table \emph{Optional: } File types for the analysis results are:
#' "csv", "xlsx", "txt". \strong{Default = "csv"}
#' @param path {Optional:} Path to the folder the results should be saved at.
#' \strong{Default = NULL}
#'
#' @return Input DF with additional column including potential additional IDs.
#'
#' @examples
#' DetectedIDs <- cellular_meta %>% tidyr::drop_na()
#' Res <- equivalent_id(data = DetectedIDs, metadata_info = c(InputID =
# "HMDB"),
#' from = "hmdb")
#'
#' @keywords Find potential additional IDs for one metabolite identifier
#'
#' @importFrom dplyr mutate select group_by ungroup distinct filter across
#' rowwise if_else
#' @importFrom tidyr separate_rows unnest
#' @importFrom purrr map_chr map_lgl map_int 
#' @importFrom tidyselect all_of starts_with
#' @importFrom rlang !!! !! := sym syms
#' @importFrom OmnipathR id_types translate_ids
#' @importFrom logger log_warn log_trace
#' @importFrom stringr str_to_lower str_split
#' @export
equivalent_id <- 
    function(
        data,
        metadata_info = c(InputID = "MetaboliteID"),
                          from = "hmdb",
        save_table = "csv",
        path = NULL
        ) {
    # FUTURE: Once we have the structural similarity tool available in OmniPath,
    # we can start creating this function!

  ### 1)
    # Check Measured ID's in prior knowledge


  ### 2)
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

  metaproviz_init()

  ## ------------------  Check Input ------------------- ##
  # HelperFunction `check_param`
    check_param(
        data = data,
        data_num = FALSE,
        save_table = save_table
    )

  # Specific checks:
    if ("InputID" %in% names(metadata_info)) {
        if (metadata_info[["InputID"]] %in% colnames(data) == FALSE) {
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

    if (nrow(doublons) > 0) {
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

  ## ------------------  Create output folders and path ------------------- ##
    if (is.null(save_table) == FALSE) {
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

  ## ------------------ Set the ID type for to ----------------- ##
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

  ## ------------------ Load manual table ----------------- ##
    if ((from == "kegg") == FALSE) {
        EquivalentFeatures <- 
            equivalent_features %>%
      select(from)
  }

  ## ------------------ Translate from-to-to ------------------- ##
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
            ambiguity_groups = NULL, # Can not be set to FALSE!
            ambiguity_summary = FALSE # Checks within the groups, without it checks across groups
        ) %>%
        select(
            all_of(intersect(names(.), names(data))), all_of(to)
        ) %>%
        mutate(
            across(all_of(to), ~ map_chr(., ~ paste(unique(.), collapse = ", "))
            )
        ) %>%
        group_by(
            !!sym(metadata_info[["InputID"]])
        ) %>%
        mutate(
            across(all_of(to), ~ paste(unique(.), collapse = ", "), .names = "{.col}")
        ) %>%
    ungroup() %>%
    distinct() %>%
        mutate(
            across(all_of(to), ~ ifelse(. == "0", NA, .))
        )


  ## ------------------ Translate to-to-from ------------------- ##
    TranslatedDF_Long <- 
        TranslatedDF %>%
            select(
                !!sym(metadata_info[["InputID"]]), !!sym(to)
            ) %>%
            rename(
                "InputID" = !!sym(metadata_info[["InputID"]])
            ) %>%
            separate_rows(
                !!sym(to), sep = ", "
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
            ambiguity_summary = FALSE  # Checks within the groups, without it checks across groups
        ) %>%
        select(
            "InputID", !!sym(to), !!sym(from)
        ) %>%
        distinct(
            InputID, !!sym(from), .keep_all = TRUE
        ) %>%
        # Remove duplicates based on InputID and from
        mutate(
            AdditionalID = if_else(InputID == !!sym(from), FALSE, TRUE)
        ) %>%
        select(
            "InputID", !!sym(from), "AdditionalID"
        ) %>%
        filter(
            AdditionalID == TRUE
        ) %>%
        mutate(
            across(
                all_of(from),
                ~ map_chr(., ~ paste(unique(.), collapse = ", "))
            )
        ) %>%
    rowwise() %>%
    mutate(
            fromList = list(str_split(!!sym(from), ",\\s*")[[1]]),
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

  ## ------------------ Merge to Input ------------------- ##
    OtherIDs <- 
        merge(
            data,
            OtherIDs,
            by.x = metadata_info[["InputID"]],
            by.y = "InputID",
            all.x = TRUE
        )

    ## ------------------- Add additional IDs -------------- ##

  if (exists("EquivalentFeatures")) {
   EquivalentFeatures$AllIDs <- EquivalentFeatures[[from]]
        EquivalentFeatures_Long <- 
            EquivalentFeatures %>%
                separate_rows(
                    !!sym(from), sep = ","
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
                                        collapse = ","
                                    ),
                                ",\\s*"
                                )
                            )
                        )
                    ),
                    collapse = ",")
            ) %>%
            ungroup() %>%
      rowwise() %>%
      mutate(
        PotentialAdditionalIDs = paste(
          setdiff(
                        unlist(
                            str_split(
                                AllIDs, ",\\s*"
                            )
                        ),
                        # Split merged_column into individual IDs
                        as.character(
                            !!sym(
                                metadata_info[["InputID"]] # Split hmdb into individual IDs
                            )
                        )
          ),
                    collapse = ", " # Combine the remaining IDs back into a comma-separated string
        )
      ) %>%
            ungroup() %>%
            select(
                -AllIDs.x, 
                -AllIDs.y
            )
  }

    ## ------------------- Fill empty cells -------------- ##
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

  ## ------------------ Create count_id plot ------------------- ##
    # QC plot of before and after
    Before <- 
        MetaProViz:::count_id(
            data = data,
            column = metadata_info[["InputID"]],
                                 save_plot = NULL,
                                 save_table = NULL,
                                 print_plot = FALSE,
            path = NULL
        )

    After <- 
        MetaProViz:::count_id(
            data = OtherIDs,
            column = "AllIDs",
                                  save_plot = NULL,
                                  save_table = NULL,
                                  print_plot = FALSE,
            path = NULL
        )


  ## ------------------ Create Output ------------------- ##
  OutputDF <- OtherIDs

  ## ------------------ Save the results ------------------- ##
  ResList <- list("equivalent_id" = OutputDF)

    suppressMessages(
        suppressWarnings(
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
        )
    )

  return(invisible(OutputDF))
}

################################################################################
### ### ###                    Mapping Ambiguity                     ### ### ###
################################################################################

#' Create Mapping Ambiguities between two ID types
#'
#' @param data Translated DF from translate_id reults or dataframe with at
#' least
#' one column with the target metabolite ID and another MetaboliteID type. One
#' of the IDs can only have one ID per row, the other ID can be either
# separated
#' by comma or a list. Optional: add other columns such as source (e.g. term).
#' @param to Column name of original metabolite identifier in data. Here should
#' only have one ID per row.
#' @param from Column name of the secondary or translated metabolite identifier
#' in data. Here can be multiple IDs per row either separated by comma " ," or
#' a
#' list of IDs.
#' @param grouping_variable \emph{Optional: } If NULL no groups are used. If
#' TRUE provide column name in data containing the grouping_variable and
#' features are grouped. \strong{Default = NULL}
#' @param summary \emph{Optional: } If TRUE a long summary tables are created.
#' \strong{Default = FALSE}
#' @param save_table \emph{Optional: } File types for the analysis results are:
#' "csv", "xlsx", "txt". \strong{Default = "csv"}
#' @param path {Optional:} Path to the folder the results should be saved at.
#' \strong{Default = NULL}
#'
#' @return List with at least 4 DFs: 1-3) from-to-to: 1. MappingIssues, 2.
#' MappingIssues summary, 3. Long summary (If summary=TRUE) & 4-6) to-to-from:
#' 4. MappingIssues, 5. MappingIssues summary, 6. Long summary (If
#' summary=TRUE)
#' & 7) Combined summary table (If summary=TRUE)
#'
#' @examples
#' KEGG_Pathways <- metsigdb_kegg()
#' InputDF <- translate_id(data = KEGG_Pathways, metadata_info = c(InputID =
#' "MetaboliteID", grouping_variable = "term"), from = c("kegg"), to =
#' c("pubchem"))[["TranslatedDF"]]
#' Res <- mapping_ambiguity(data = InputDF, from = "MetaboliteID", to =
#' "pubchem", grouping_variable = "term", summary = TRUE)
#'
#' @keywords Mapping ambiguity
#'
#' @importFrom dplyr mutate bind_cols bind_rows
#' @importFrom rlang !!! !! := sym syms
#' @importFrom OmnipathR ambiguity
#' @importFrom logger log_trace
#' @importFrom tidyr unnest
#' @importFrom purrr map_chr map_int map_lgl
#' @importFrom dplyr first
#' 
#' @export
mapping_ambiguity <- 
    function(
        data,
                             from,
                             to,
                             grouping_variable = NULL,
        summary = FALSE,
        save_table = "csv",
        path = NULL
) {
            
  metaproviz_init()
  ## ------------------  Check Input ------------------- ##
  # HelperFunction `check_param`
    check_param(
        data = data,
        data_num = FALSE,
        save_table = save_table
    )

  # Specific checks:
    if (from %in% colnames(data) == FALSE) {
        message <- 
            paste0(
                from,
                " column was not found in data. Please check your input."
            )
        log_trace(
            paste(
                "Error ", message, sep = ""
            )
        )
      stop(message)
      }

    if (to %in% colnames(data) == FALSE) {
        message <- 
            paste0(
                to,
                " column was not found in data. Please check your input."
            )
        log_trace(
            paste(
                "Error ", message, sep = ""
            )
        )
    stop(message)
  }

    if (is.null(grouping_variable) == FALSE) {
        if (grouping_variable %in% colnames(data) == FALSE) {
            message <- 
                paste0(
                    grouping_variable,
                    " column was not found in data. Please check your input."
                )
            log_trace(
                paste(
                    "Error ", message, sep = ""
                )
            )
      stop(message)
    }
  }

    if (is.logical(summary) == FALSE) {
        message <- 
            paste0(
                "Check input. The summary parameter should be either =TRUE or ",
                "=FALSE."
            )
        log_trace(
            paste(
                "Error ", message, sep = ""
            )
        )
    stop(message)
  }

    ## ------------------  General checks of wrong occurences
    ## ------------------- ##
    # Task 1: Check that from has no duplications within one group --> should
    # not be the case --> remove duplications and inform the user/ ask if they
    # forget to set groupings column
    # Task 2: Check that from has the same items in to across the different
    # entries (would be in different Groupings, otherwise there should not be
    # any duplications) --> List of Miss-Mappings across terms

    # FYI: The above can not happen if our translateID function was used, but
    # may be the case when the user has done something manually before


  ## ------------------  Create output folders and path ------------------- ##
    if (is.null(save_table) == FALSE) {
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

    ############################################################################
  ## ------------------  Prepare Input data ------------------- ##
    # If the user provides a DF where the to column is a list of IDs, then we
    # can use it right away
    # If the to column is not a list of IDs, but a character column, we need to
    # convert it into a list of IDs
    if (is.character(data[[to]]) == TRUE) {
        data[[to]] <- 
            data[[to]] %>%
            strsplit(", ") %>%
      lapply(as.character)
  }

  ## ------------------  Perform ambiguity mapping ------------------- ##
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
            )]] <- data %>%
            unnest(  # unlist the columns in case they are not expaned
                cols = all_of(Comp[[comp]]$from)
            ) %>%
            filter(  # Remove NA values, otherwise they are counted as column is character
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
            )]] <- ResList[[
                paste0(
                    Comp[[comp]]$from,
                    "-to-",
                    Comp[[comp]]$to,
                    sep = ""
                )]] %>%
                attr(
                    paste0(
                        "ambiguity_",
                        Comp[[comp]]$from,
                        "_",
                        Comp[[comp]]$to,
                        sep = ""
                    )
                )

        ########################################################################
        if (summary == TRUE) {
            if (is.null(grouping_variable) == FALSE) {
                # Add further information we need to summarise the table and
                # combine Original-to-Translated and Translated-to-Original
                # If we have a grouping_variable we need to combine it with the
                # MetaboliteID before merging
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
                    mutate(
                        !!sym(
                            paste0(
                                "AcrossGroupMappingIssue(",
                                Comp[[comp]]$from,
                                "_to_",
                                Comp[[comp]]$to,
                                ")",
                                sep = ""
                            )
                        ) := case_when(
                            !!sym(
                                paste0(
                                    Comp[[comp]]$from,
                                    "_",
                                    Comp[[comp]]$to,
                                    "_ambiguity_bygroup",
                                    sep = ""
                                )
                            ) != !!sym(
                                paste0(
                                    Comp[[comp]]$from,
                                    "_",
                                    Comp[[comp]]$to,
                                    "_ambiguity",
                                    sep = ""
                                )
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
                            !!sym(Comp[[comp]]$from) == 0,
                            NA, # Or another placeholder
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
                            all(!!sym(Comp[[comp]]$to) == 0),
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
                        c(Comp[[comp]]$from,
                        Comp[[comp]]$to),
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
                            !!sym(Comp[[comp]]$from) == 0,
                            NA, # Or another placeholder
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
                            all(!!sym(Comp[[comp]]$to) == 0),
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
                        !!sym(Comp[[comp]]$to), sep = ", "
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
            unnest(   # unlist the columns in case they are not expaned
               cols = all_of(Comp[[comp]]$from)
            ) %>%
            filter(
                is.na(!!sym(Comp[[comp]]$from))
            )
        if (nrow(Removed) > 0) {
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

  ## ------------------ Create summaryTable ------------------- ##
    if (summary == TRUE) {
    # Combine the two tables
        summary <- 
            merge(
                x = ResList[[
                    paste0(
                        from, "-to-", to, "_Long", 
                        sep = ""
                    )]][,
                        c(
                            "UniqueID", 
                            paste0(from, "_to_", to),
                            paste0("Count(", from, "_to_", to, ")"),
                            paste0("AcrossGroupMappingIssue(", from, "_to_", to, ")", sep = "")
                        )
                    ],
                y = ResList[[
                    paste0(
                        to, "-to-", from, "_Long", 
                        sep = ""
                    )]][,
                        c(
                            "UniqueID", 
                            paste0(to, "_to_", from), 
                            paste0("Count(", to, "_to_", from, ")"),
                            paste0("AcrossGroupMappingIssue(", to, "_to_", from, ")", sep = "")
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
    summary <- summary %>%
            mutate(
                Mapping = case_when(
                    !!sym( paste0("Count(", from, "_to_", to, ")") ) == 1   &   !!sym( paste0("Count(", to, "_to_", from, ")") ) == 1   ~ "one-to-one",
                    !!sym( paste0("Count(", from, "_to_", to, ")") ) > 1    &   !!sym( paste0("Count(", to, "_to_", from, ")") ) == 1   ~ "one-to-many", 
                    !!sym( paste0("Count(", from, "_to_", to, ")") ) > 1    &   !!sym( paste0("Count(", to, "_to_", from, ")") ) > 1    ~ "many-to-many",
                    !!sym( paste0("Count(", from, "_to_", to, ")") ) == 1   &   !!sym( paste0("Count(", to, "_to_", from, ")") ) > 1    ~ "many-to-one",
                    !!sym( paste0("Count(", from, "_to_", to, ")") ) >= 1   &   !!sym( paste0("Count(", to, "_to_", from, ")") ) == NA  ~ "one-to-none",
                    !!sym( paste0("Count(", from, "_to_", to, ")") ) >= 1   &   is.na( !!sym(paste0("Count(", to, "_to_", from, ")")))  ~ "one-to-none",
                    !!sym( paste0("Count(", from, "_to_", to, ")") ) == NA  &   !!sym( paste0("Count(", to, "_to_", from, ")") ) >= 1   ~ "none-to-one",
                    is.na( !!sym(paste0("Count(", from, "_to_", to, ")")) ) &   !!sym( paste0("Count(", to, "_to_", from, ")") ) >= 1   ~ "none-to-one",
                    TRUE ~ NA
                )) %>%
                mutate(
                    !!sym(
                        paste0("Count(", from, "_to_", to, ")")
                    ) := replace_na(
                        !!sym(
                            paste0("Count(", from, "_to_", to, ")")
                        ),
                        0
                    )
                ) %>%
                mutate(
                    !!sym(
                        paste0("Count(", to, "_to_", from, ")")
                    ) := replace_na(
                        !!sym(
                            paste0("Count(", to, "_to_", from, ")")
                        ),
                        0
                    )
                )

    ResList[["summary"]] <- summary
  }

  ## ------------------ Save the results ------------------- ##
    suppressMessages(
            suppressWarnings(
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
            )
        )

    # Return
    invisible(
        return(ResList)
    )
}

################################################################################
### ### ### Check Measured ID's in prior knowledge ### ### ###
################################################################################

#' Check and summarize relationship between prixor knowledge to measured
#' features
#'
#' @param data dataframe with at least one column with the detected metabolite
#' IDs (e.g. HMDB). If there are multiple IDs per detected peak, please
#' separate
#' them by comma ("," or ", " or chr list). If there is a main ID and
#' additional
#' IDs, please provide them in separate columns.
#' @param input_pk dataframe with at least one column with the metabolite ID
#' (e.g. HMDB) that need to match data metabolite IDs "source" (e.g. term). If
#' there are multiple IDs, as the original pathway IDs (e.g. KEGG) where
#' translated (e.g. to HMDB), please separate them by comma ("," or ", " or chr
#' list).
#' @param metadata_info Colum name of Metabolite IDs in data and input_pk as
#' well as column name of grouping_variable in input_pk. \strong{Default =
#' c(InputID="HMDB", PriorID="HMDB", grouping_variable="term")}
#' @param save_table \emph{Optional: } File types for the analysis results are:
#' "csv", "xlsx", "txt". \strong{Default = "csv"}
#' @param path {Optional:} Path to the folder the results should be saved at.
#' \strong{Default = NULL}
#'
#' @return A \code{list} with three elements:
#' \itemize{
#' \item \code{data_summary} — a data frame summarising matching results per
#' input ID, including counts, conflicts, and recommended actions.
#' \item \code{GroupingVariable_summary} — a detailed data frame showing
# matches
#' grouped by the specified variable, with conflict annotations.
#' \item \code{data_long} — a merged data frame of prior knowledge IDs and
#' detected IDs in long format.
#' }
#' 
#' @examples
#' DetectedIDs <- cellular_meta %>%
#' dplyr::select("Metabolite", "HMDB") %>%
#' tidyr::drop_na()
#' input_pathway <- translate_id(data = metsigdb_kegg(), metadata_info =
#' c(InputID = "MetaboliteID", grouping_variable = "term"), from = c("kegg"),
#' to
#' = c("hmdb"))[["TranslatedDF"]] %>% tidyr::drop_na()
#' Res <- checkmatch_pk_to_data(data = DetectedIDs, input_pk = input_pathway,
#' metadata_info = c(InputID = "HMDB", PriorID = "hmdb", grouping_variable =
#' "term"))
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
checkmatch_pk_to_data <- 
    function(
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

  ## ------------ Create log file ----------- ##
  metaproviz_init()

  ## ------------ Check Input files ----------- ##

  ## data:
    if ("InputID" %in% names(metadata_info)) {
        if (metadata_info[["InputID"]] %in% colnames(data) == FALSE) {
            message <- 
                paste0(
                    "The ",
                    metadata_info[["InputID"]],
                    " column selected as InpuID in metadata_info was not found in",
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

    ### This is after the main input checks (before NA removal), so we will save
    ### original df here for later merging to get the Null and duplicates back.
  data_Original <- data

    if (sum(is.na(data[[metadata_info[["InputID"]]]])) >= 1) {
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

    if (nrow(data) - nrow(distinct(data, .data[[metadata_info[["InputID"]]]])) >= 1) {
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
                    .data[[metadata_info[["InputID"]]]], .keep_all = TRUE
                )

    warning(message)
  }

    data_MultipleIDs <- 
        any(
            grepl(",\\s*", data[[metadata_info[["InputID"]]]]) |   # Comma-separated
       map_lgl(
         data[[metadata_info[["InputID"]]]],
         function(x) {
           if (grepl("^c\\(|^list\\(", x)) {
                        parsed <- tryCatch(
                            eval(parse(text = x)),
                            error = function(e) NULL
                        )
                        return(
                            is.list(parsed) && length(parsed) > 1 || is.vector(parsed) && length(parsed) > 1
                        )
           }
           FALSE
         }
       )
   )

  ## input_pk:
    if ("PriorID" %in% names(metadata_info)) {
        if (metadata_info[["PriorID"]] %in% colnames(input_pk) == FALSE) {
            message <- 
                paste0(
                    "The ",
                    metadata_info[["PriorID"]],
                    " column selected as InpuID in metadata_info was not found in",
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

    ### This is after the main input checks (before NA removal), so we will save
    ### original df here for later merging to get the Null and duplicates back.
  prior_knowledge_Original <- input_pk

    if (sum(is.na(input_pk[[metadata_info[["PriorID"]]]])) >= 1) {
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
        if (metadata_info[["grouping_variable"]] %in% colnames(input_pk) == FALSE) {
            message <- 
                paste0(
                    "The ",
                    metadata_info[["grouping_variable"]],
                    " column selected as InpuID in metadata_info was not found in",
                    "input_pk."
                    "Please check your input."
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

    if (nrow(input_pk) - nrow(distinct(input_pk, .data[[metadata_info[["PriorID"]]]], .data[[metadata_info[["grouping_variable"]]]])) >= 1) {
    # Remove duplicate IDs
        message <- 
            paste0(
                nrow(input_pk) - nrow(distinct(input_pk, .data[[metadata_info[["PriorID"]]]], .data[[metadata_info[["grouping_variable"]]]])),
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
                ",\\s*", 
                input_pk[[metadata_info[["PriorID"]]]]) |   # Comma-separated
      map_lgl(
        input_pk[[metadata_info[["PriorID"]]]],
        function(x) {
          if (grepl("^c\\(|^list\\(", x)) {
                            parsed <- tryCatch(
                                eval(parse(text = x)),
                                error = function(e) NULL
                            )
                            return(
                                is.list(parsed) && length(parsed) > 1 || is.vector(parsed) && length(parsed) > 1
                            )
          }
          FALSE
        }
      )
  )

  ## ------------ Create Results output folder ----------- ##
    if (is.null(save_table) == FALSE) {
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

    ############################################################################
    ## ------------ Check how IDs match and if needed remove unmatched IDs
    ## ----------- ##

  # 1.  Create long DF
    create_long_df <- 
        function(df, id_col, df_name) {
    df %>%
            mutate(
                row_id = row_number()
            ) %>%
            mutate(   # Store original values
                !!paste0("OriginalEntry_", df_name, sep = "") := 
                    !!sym(id_col)
            ) %>%
            separate_rows(
                !!sym(id_col), sep = ",\\s*"
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
            summary_df$found_match_in_PK[i] <- NA  # could also set it to FALSE but making NA for now for plotting
      summary_df$matches[i] <- NA
    } else {
      # Split each cell into individual entries and trim whitespace
            entries <- 
                trimws(
                    unlist(
                        strsplit(
                            as.character(Values_data[i]),
                            ",\\s*"
                        )
                    )
                )
            # delimiter = "," or ", "

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
            if (length(entries) > 0) {
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
      is.na(!!sym(metadata_info[["grouping_variable"]])) ~ NA_integer_,
      TRUE ~ n_distinct(!!sym(metadata_info[["PriorID"]]), na.rm = TRUE)
    )
            ) %>%
            mutate(
      Group_Conflict_Notes = case_when(
                    any(
                        Count_FeatureIDs_to_GroupingVariable > 1, na.rm = TRUE
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
                    original_count >= 1 & matches_count == 1 & Unique_GroupingVariable_count >= 1 ~ "None",
                    original_count >= 1 & matches_count == 0 & Unique_GroupingVariable_count >= 0 ~ "None",
                    original_count > 1  & matches_count > 1  & Unique_GroupingVariable_count == 1 ~ "Check",
                    original_count > 1  & matches_count > 1  & Unique_GroupingVariable_count > 1  ~ "Check",
        TRUE ~ NA_character_
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
                    original_count == 1 & matches_count <= 1 ~ metadata_info[["InputID"]], 
                    original_count >= 2 & matches_count == 0 ~ str_split(metadata_info[["InputID"]], ",\\s*") %>% map_chr(first),
                    original_count >= 2 & matches_count == 1 ~ matches,
        TRUE ~ NA_character_
      ),
      Action_Specific = case_when(
        matches_count >= 2 & Group_Conflict_Notes == "None" ~ "KeepEachID",
        matches_count >= 2 & Group_Conflict_Notes != "None" ~ "KeepOneID",
        TRUE ~ "None"
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
            nrow(summary_df_short %>% filter(matches_count >= 1)),
            " match, which is ",
            ( nrow(summary_df_short %>% filter(matches_count >= 1)) / n_distinct(unique(data_long[[metadata_info[["InputID"]]]])) ) * 100,
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
            nrow(summary_df_short %>% filter(matches_count >= 1)),
            " are detected in the data, which is ",
            ( nrow(summary_df_short %>% filter(matches_count >= 1)) / n_distinct(PK_long[[metadata_info[["PriorID"]]]]) ) * 100,
            "%."
        )

  message(message0)
  message(message1)
  message(message2)

    if (nrow(summary_df_short %>% filter(ActionRequired == "Check")) >= 1) {
        warning1 <-    # "Check"
            paste0(   
                "There are cases where multiple detected IDs match to ",
                "multiple ",
                "prior knowledge ",
                "IDs of the same category"
            )
    warning(warning1)
  }

  ## ------------------ Plot summary ----------------------##
    # x = "Class" and y = Frequency. Match Status can be colour of if no class
    # provided class = Match status.
  # Check Biocrates code.


  ## ------------------ Save Results ----------------------##
    ResList <- 
        list(
            "data_summary" = summary_df_short,
                  "GroupingVariable_summary" = summary_df,
            "data_long" = merged_df
        )

    suppressMessages(
        suppressWarnings(
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
        )
    )

    # Return
   invisible(return(ResList))
}


################################################################################
### ### ### cluster Prior Knowledge ### ### ###
################################################################################

#' Deal with pathway overlap in prior knowledge
#'
#' @param data dataframe with at least one column with the target (e.g.
#' metabolite) and a column source (e.g. term).
#' @param metadata_info = c(InputID="MetaboliteID", grouping_variable="term"),
#'
#' @examples
#' metsigdb_kegg()
#'
#' @importFrom dplyr bind_rows filter group_by left_join mutate
#' @importFrom dplyr select summarize ungroup
#' @importFrom igraph graph_from_adjacency_matrix components
#' @noRd
cluster_pk <- 
    function(
        data,
        metadata_info = c(
            InputID = "MetaboliteID",
            grouping_variable = "term"
        ),
        clust = "Graph", 
        matrix = "percentage",
        min = 2  
    ) {
    # Notes about function arguments
    ## data:    This can be either the original PK (e.g. KEGG pathways), but it 
    ##          can also be the output of enrichment results (--> meaning here 
    ##          we would cluster based on detection!)
    ## clust:   Options: "Graph", "Hierarchical",
    ## matrix:  Choose "pearson", "spearman", "kendall", or "percentage"
    ## min:     # minimum pathways per cluster


    # cluster PK before running enrichment analysis --> add another column that
    # groups the data based on the pathway overlap:
    # provide different options for clustering (e.g. % of overlap, semantics
    # similarity) --> Ramp uses % of overlap, semnatics similarity:
    # https://yulab-smu.top/biomedical-knowledge-mining-book/GOSemSim.html


  ## ------------------ Check Input ------------------- ##

  ## ------------------ Create output folders and path ------------------- ##

    ############################################################################
  ## ------------------ cluster the data ------------------- ##
  # 1. Create a list of unique MetaboliteIDs for each term
    term_metabolites <- 
        data %>%
            group_by(
                !!sym(metadata_info[["grouping_variable"]])
            ) %>%
            summarize(
                MetaboliteIDs = list(
                    unique(
                        !!sym(metadata_info[["InputID"]])
                    )
                )
            ) %>%
    ungroup()

    # 2. Create the overlap matrix based on different methods:
    if (matrix == "percentage") {
    # Compute pairwise overlaps
        term_overlap <-
            combn(
                term_metabolites[[metadata_info[["grouping_variable"]]]], 
                2,
                function(terms) {
                    term1_ids <-
                        term_metabolites$MetaboliteIDs[term_metabolites[[metadata_info[["grouping_variable"]]]] == terms[1]][[1]]
                    term2_ids <-
                        term_metabolites$MetaboliteIDs[term_metabolites[[metadata_info[["grouping_variable"]]]] == terms[2]][[1]]

                    overlap <- 
                        length(
                            intersect(term1_ids,term2_ids)) / length(union(term1_ids, term2_ids)
                        )
                    data.frame(
                        Term1 = terms[1], 
                        Term2 = terms[2], 
                        Overlap = overlap
                    )
                }, 
                simplify = FALSE
            ) %>%
      bind_rows()

        # Create overlap matrix: An overlap matrix is typically used to quantify
        # the degree of overlap between two sets or groups.
    # overlap coefficient (or Jaccard Index) Overlap(A,B)= ∣A∪B∣ / ∣A∩B
        # The overlap matrix measures the similarity between sets or groups
        # based on common elements.
        terms <- 
            unique(
                c(
                    term_overlap$Term1, 
                    term_overlap$Term2
                )
            )
        overlap_matrix <- 
            matrix(
                1,
                nrow = length(terms),
                ncol = length(terms),
                dimnames = list(terms, terms)
            )
    for (i in seq_len(nrow(term_overlap))) {
      t1 <- term_overlap$Term1[i]
      t2 <- term_overlap$Term2[i]
      overlap_matrix[t1, t2] <- 1 - term_overlap$Overlap[i]
      overlap_matrix[t2, t1] <- 1 - term_overlap$Overlap[i]
    }
  } else {
    # Create a binary matrix for correlation methods
        terms <- 
            term_metabolites[[metadata_info[["grouping_variable"]]]]
        metabolites <- 
            unique(
                unlist(
                    term_metabolites$MetaboliteIDs
                )
            )
        # [[metadata_info[["InputID"]]]]

        binary_matrix <- 
            matrix(
                0,
                nrow = length(terms),
                ncol = length(metabolites),
                dimnames = list(terms, metabolites)
            )
    for (i in seq_along(terms)) {
            metabolites_for_term <- 
                term_metabolites$MetaboliteIDs[[i]]
            # [[metadata_info[["InputID"]]]]
      binary_matrix[i, colnames(binary_matrix) %in% metabolites_for_term] <- 1
    }

        # Compute correlation matrix: square matrix used to represent the
        # pairwise correlation coefficients between variables or terms
        # correlation matrix 𝐶 C is an 𝑛 × 𝑛 n×n matrix where each element 𝐶 𝑖 𝑗
        # C ij  is the correlation coefficient between the variables 𝑋𝑖 Xi  and
        # 𝑋 𝑗 X j
        # The correlation matrix measures the strength and direction of linear
        # relationships between variables.
        correlation_matrix <- 
            cor(
                t(binary_matrix), 
                method = matrix
            )

    # Convert to distance matrix
        overlap_matrix <- 
            1 - correlation_matrix
  }
  # 3. cluster terms based on overlap threshold
    threshold <- 0.7   # Define similarity threshold
    
    term_clusters <- 
        term_overlap %>%
            filter(
                Overlap >= threshold
            ) %>%
            select(
                Term1, 
                Term2
            )

  # 4. clustering
    if (clust == "Graph") { # Use Graph-based clustering
  # Here we need the distance matrix:
        overlap_matrix <- 
            1 - correlation_matrix

        # An adjacency matrix represents a graph structure and encodes the
        # relationships between nodes (vertices)
  # Add weight (can also represent unweighted graphs)

        # Applying Gaussian kernel to convert distance into similarity
        adjacency_matrix <- 
            exp(-overlap_matrix^2)

  # Create a graph from the adjacency matrix
        g <- 
            graph_from_adjacency_matrix(
                adjacency_matrix,
                mode = "undirected",
                weighted = TRUE
            )
        initial_clusters <- 
            components(g)$membership
        term_metabolites$cluster <-
            initial_clusters[match(term_metabolites[[metadata_info[["grouping_variable"]]]], names(initial_clusters))]
    } else if (clust == "Hierarchical") {
    # Hierarchical clustering
        hclust_result <- 
            hclust(
                as.dist(distance_matrix), 
                method = "average"
            )
        # make methods into parameters!
    num_clusters <- 4
        term_clusters_hclust <- 
            cutree(
                hclust_result, 
                k = num_clusters
            )

        term_metabolites$cluster <- 
            paste0(
                "cluster",
                term_clusters_hclust[match(terms, names(term_clusters_hclust))]
            )
        # term_metabolites$cluster <- clusters[...] (refer to previous block) ... ) names(clusters))]
  } else {
    stop("Invalid clustering method specified in clust parameter.")
  }

  # 5. Merge cluster group information back to the original data
    df <- 
        data %>%
            left_join(
                term_metabolites %>%
                    select(
                        !!sym(metadata_info[["grouping_variable"]]), 
                        cluster
                    ),
                by = metadata_info[["grouping_variable"]]
            ) %>%
            mutate(
                cluster = ifelse(
      is.na(cluster),
                    "None",
                    paste0("cluster", cluster)  # Convert numeric IDs to descriptive labels
      )
    )

  # 6. Summarize the clustering results

  ## ------------------ Save and return ------------------- ##
}

################################################################################
### ###     Helper function to add term information to Enrichment Results    ###
################################################################################

# Better function Name and parameter names needed
# Use in ORA functions and showcase in vignette with decoupleR output

#' Adds extra columns to enrichment output
#'
#' These columns inform about 1. The amount of genes associated with term in prior knowledge, 2. The amount of genes detected in input data associated with term in prior knowledge, and 3. The percentage of genes detected in input data associated with term in prior knowledge.
#'
#' @param mat data matrix used as input for enrichment analysis
#' @param net Prior Knowledge used as input for enrichment analysis
#' @param res Results returned from the enrichment analysis
#' @param .source used as input for enrichment analysis
#' @param .target used as input for enrichment analysis
#' @param complete TRUE or FALSE, weather only .source with results should be returned or all .source in net.
#'
#' @importFrom tibble rownames_to_column
#' @noRd
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
    rownames_to_column("Symbol")

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
#' @param metadata_info A named list (with names matching those in \code{data}) where each element is either a character string or a
#'        character vector indicating the column name(s) to extract features. For multiple-resource comparisons, these refer to the columns
#'        containing feature identifiers. For within-resource comparisons, the vector should list the columns to compare (e.g., \code{c("CHEBI", "HMDB", "LIMID")}).
#'        In within-resource mode, the input data frame is expected to contain a column named \code{"Class"} (or a grouping column specified via the
#'        \code{class_col} attribute). \emph{If no grouping column is found, a default grouping column named \code{"Group"} (with all rows assigned the same value) is created.}
#' @param filter_by Character. Optional filter for the resulting features when comparing multiple resources.
#'        Options are: \code{"both"} (default), \code{"gene"}, or \code{"metabolite"}. This parameter is ignored in within-resource mode.
#' @param plot_name \emph{Optional: } String which is added to the output files of the Upsetplot \strong{Default = ""}
#' @param name_col \emph{Optional: } column name including the feature names. Default is \code{"TrivialName"}.
#' @param palette_type Character. Color palette to be used in the plot. Default is \code{"polychrome"}.
#' @param save_plot \emph{Optional: } Select the file type of output plots. Options are svg, png, pdf. \strong{Default = svg}
#' @param save_table \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt". \strong{Default = "csv"}
#' @param print_plot \emph{Optional: } TRUE or FALSE, if TRUE Volcano plot is saved as an overview of the results. \strong{Default = TRUE}
#' @param path \emph{Optional:} Path to the folder the results should be saved at. \strong{Default = NULL}
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
#' ## Example 1: Within-Resource Comparison (Comparing Columns Within a Single data Frame)
#'
#' # biocrates_features is a data frame with columns: "TrivialName", "CHEBI", "HMDB", "LIMID", and "Class".
#' # Here the "Class" column is used as the grouping variable in the UpSet plot.
#' data_single <- list(Biocft = biocrates_features)
#' metadata_info_single <- list(Biocft = c("CHEBI", "HMDB", "LIMID"))
#'
#' res_single <- compare_pk(data = data_single, metadata_info = metadata_info_single,
#'                           plot_name = "Overlap of BioCrates Columns")
#'
#' ## Example 2: Custom data Frames with Custom Column Names
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
#' @importFrom logger log_trace
#' @export
compare_pk <- function(data,
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
  ###########################################################################
  ## ------------ Create log file ----------- ##
  metaproviz_init()

  ## ------------ Check Input files ----------- ##
  # Match filter argument
  filter_by <- match.arg(filter_by)

  # Validate data input
  if (!is.list(data) || length(data) < 1) {
    message <- paste0("data must be a non-empty list.")
    log_trace(paste("Error ", message, sep=""))
    stop(message)
  }
  if (is.null(names(data)) || any(names(data) == "")) {
    message <- paste0("data must be a named list with resource names.")
    log_trace(paste("Error ", message, sep=""))
    stop(message)
  }

  ###########################################################################
  ## ------------ Create folders ----------- ##
  if(is.null(save_plot)==FALSE |is.null(save_table)==FALSE){
    folder <- save_path(folder_name= "PK",
                        path=path)

    subfolder <- file.path(folder, "ComparePK")
    if (!dir.exists(subfolder)) {dir.create(subfolder)}
  }

  ###########################################################################
  ## ----------- Input ----------- ##
  # Define resource lookup table with information on how to retrieve and transform each resource.
  default_cols <- list(
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
      plot_name = plot_name,
      palette_type = palette_type,
      save_plot = NULL,
      print_plot = FALSE
    )

    summary_table <- df_summary
  } else {
    # ===== Multi-Resource Comparison Mode =====
    # Process each resource in data.
    for (res in names(data)) {
      resource_val <- data[[res]]
			resource_id <- tolower(res)
			if (is.null(metadata_info[[res]]) && resource_id %in% names(default_cols)) {
				metadata_info[[res]] <- default_cols[[resource_id]]
			} else if (is.null(metadata_info[[res]])) {
				stop("metadata_info must be provided for resource: ", res)
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
      plot_name = plot_name,
      palette_type = palette_type,
      save_plot = NULL,
      print_plot = FALSE
    )

    summary_table <- df_binary
  }

  ###########################################################################
  ###-------- Save and return ----------###
  suppressMessages(suppressWarnings(
    save_res(inputlist_df=list(summary_table = summary_table),
             inputlist_plot= list(upset_plot = upset_plot),
             save_table=save_table,
             save_plot=save_plot,
             path= subfolder,
             file_name= "compare_pk",
             core=FALSE,
             print_plot=print_plot)))

  return(list(summary_table = summary_table, upset_plot = upset_plot))
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
#' @examples
#' count_id(biocrates_features, "HMDB")
#'
#' @return A list with two elements:
#'   \item{result}{A data frame that includes three additional columns: \code{was_na} (logical indicator
#'                  of missing or empty cells), \code{entry_count} (number of entries in each cell), and
#'                  \code{id_label} (a categorical label based on the entry count).}
#'   \item{plot}{A \code{ggplot} object representing the histogram of entry counts.}
#'
#' @importFrom dplyr case_when mutate rename
#' @importFrom ggplot2 aes element_rect element_text expansion geom_bar
#' @importFrom ggplot2 geom_histogram ggplot labs scale_fill_manual scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous theme theme_classic
#' @importFrom grid convertUnit
#' @export
count_id <- function(data,
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
  #add count_id function into Equivalent IDs: make a plot before and after equivalent IDs and put them side by side to return as the QC plot of the function
  #add save_plot and save_table
  #return both standard plot and Plot_Sized (=nice version of the plot which is saved. the other plot is still returned, since this is the ggplot version and can be changed further)
  #create subtitle and not title prefix


  ## ------------------  logger initiation ------------------- ##


  ## ------------------  Checks ------------------- ##


  ## ------------------  Create output folders and path ------------------- ##
  if(is.null(save_table)==FALSE ){
    folder <- save_path(folder_name= "PK",
                       path=path)

    Subfolder <- file.path(folder, "CountIDs")
    if (!dir.exists(Subfolder)) {dir.create(Subfolder)}
  }

  ## ------------------  data table ------------------- ##
  # Process the data: count entries and label each cell based on the number of entries.
  processed_data <- mutate(
    data,
    was_na = is.na(.data[[column]]) | .data[[column]] == "",
    entry_count = map_int(
      .data[[column]],
      function(cell) {
        if (is.na(cell) || cell == "") {
          0L  # Treat NA or empty as 0 entries for counting
        } else {
          as.integer(length(unlist(strsplit(as.character(cell), delimiter))))
        }
      }
    ),
    id_label = case_when(
      entry_count == 0 ~ "No ID",
      entry_count == 1 ~ "Single ID",
      entry_count >= 2 ~ "Multiple IDs"
    )
  )


  ## ------------------  plot ------------------- ##
  # Generate the plot title if not provided
  if (is.null(title_prefix)) {
    plot_name <- paste("Number of", column, "IDs per measured peak.")
  } else {
    plot_name <- title_prefix
  }

  # Create plot using ggplot2
  if(length(unique(processed_data$entry_count))>1){# Create the histogram plot using ggplot2
    plot_obj <- ggplot(processed_data, aes(x = entry_count, fill = id_label)) +
      geom_histogram(binwidth = binwidth, boundary = -0.5, color = "black")+# Set axis to start at 0
      scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+
      scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
  }else{# Create bargraph plot using ggplot2
    plot_obj <- ggplot(processed_data, aes(x = as.factor(entry_count), fill = id_label)) +
      geom_bar(color = "black")+  # Set axis to start at 0
     scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

  }

  plot_obj <- plot_obj +
    scale_fill_manual(values = fill_colors) +
    labs(title = plot_name,
                  x = "Number of IDs",
                  y = "Frequency",
                  fill = paste0(column, " IDs", sep="")) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      legend.position.inside = c(0.8, 0.8),
      legend.justification = c("right", "top"),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12)
    )

  # Make the nice plot:
  Plot_Sized <-  plot_grob_superplot(input_plot=plot_obj, metadata_info= c(Conditions="id_label", Superplot = TRUE), metadata_sample= processed_data%>%dplyr::rename("Conditions"="entry_count") , plot_name = plot_name, subtitle = "", plot_type="Bar")
  plot_height <- convertUnit(Plot_Sized$height, 'cm', valueOnly = TRUE)
  plot_width <- convertUnit(Plot_Sized$width, 'cm', valueOnly = TRUE)
  Plot_Sized %<>%
    {ggplot() + annotation_custom(.)} %>%
    add(theme(panel.background = element_rect(fill = "transparent")))

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

