# !/usr/bin/env Rscript

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
# Get KEGG prior knowledge
#

#' Metabolites excluded from prior knowledge resources
#'
#' Builds an exclusion table from hard-coded seed identifiers and translates
#' them across ID systems (HMDB, KEGG, CHEBI, PUBCHEM).
#'
#' @return A data frame with columns `metabolite_id`, `class`, and `id_type`.
#'
#' @examples
#' \dontrun{
#' get_exclusion_metabolites()
#' }
#'
#' @export
get_exclusion_metabolites <- function() {

    # example partial dataframes to exclude more metabolite sets:
    ## add into the seed_df below:
    
    #     data.frame(  ## imaginary test set
    #         seed_id = c(
    #             "HMDB0000094", "HMDB0001165", "HMDB0001405", "HMDB0001272"
    #         ),
    #         class = "ions",
    #         seed_type = "HMDB",
    #         stringsAsFactors = FALSE
    #     ),
    # data.frame(  ## imaginary test set
    #     seed_id = c(
    #         "C00533", "C00027", "C00704", "C00001", "C00011", "C00544",
    #         "C00014", "C00288", "C00007", "C00013", "C00088", "C00244",
    #         "C00282", "C19909"
    #     ),
    #     class = "small_molecules",
    #     seed_type = "KEGG",
    #     stringsAsFactors = FALSE
    # ),
    # data.frame(   ## imaginary test set
    #     seed_id = c(
    #         "CHEBI:16236", "CHEBI:16856", "CHEBI:17754", "CHEBI:16716"
    #     ),
    #     class = "xenobiotics",
    #     seed_type = "CHEBI",
    #     stringsAsFactors = FALSE
    # )
        
    seed_df <- rbind(
        data.frame(  ## KEGG ions
            seed_id = c(
                "C00076", "C00238", "C01330", "C00080", "C00698", 
                "C00162", "C00704", "C00001", "C00011", "C00070", 
                "C14818", "C00038", "C00291", "C01342"
            ),
            class = "ions",
            seed_type = "KEGG",
            stringsAsFactors = FALSE
        ),
        data.frame(  ## KEGG small molecules
            seed_id = c(
                "C00533", "C00027", "C00704", "C00001", "C00011",
                "C16844", "C00014", "C00288", "C00007", "C00013", 
                "C22381", "C00088", "C00244", "C00282", "C01322",
                "C00703"
            ),
            class = "small_molecules",
            seed_type = "KEGG",
            stringsAsFactors = FALSE
        ),
        data.frame(   ## CHEBI monoatomic ions (children nodes of root CHEBI:24867)
            seed_id = c(
                "CHEBI:149691", "CHEBI:149692", "CHEBI:15378", "CHEBI:15858", "CHEBI:16042",
                "CHEBI:16382", "CHEBI:16793", "CHEBI:17051", "CHEBI:17996", "CHEBI:18420",
                "CHEBI:23336", "CHEBI:23378", "CHEBI:23905", "CHEBI:23906", "CHEBI:24636",
                "CHEBI:24867", "CHEBI:24875", "CHEBI:25155", "CHEBI:25158", "CHEBI:25197",
                "CHEBI:25198", "CHEBI:25213", "CHEBI:25414", "CHEBI:25430", "CHEBI:25516",
                "CHEBI:26937", "CHEBI:27153", "CHEBI:27365", "CHEBI:29033", "CHEBI:29034",
                "CHEBI:29035", "CHEBI:29036", "CHEBI:29037", "CHEBI:29041", "CHEBI:29101",
                "CHEBI:29102", "CHEBI:29103", "CHEBI:29104", "CHEBI:29105", "CHEBI:29108",
                "CHEBI:29233", "CHEBI:29234", "CHEBI:29352", "CHEBI:29435", "CHEBI:29436",
                "CHEBI:29832", "CHEBI:30030", "CHEBI:30120", "CHEBI:30144", "CHEBI:30150",
                "CHEBI:30151", "CHEBI:30165", "CHEBI:30166", "CHEBI:30168", "CHEBI:30180",
                "CHEBI:30216", "CHEBI:30220", "CHEBI:30240", "CHEBI:30399", "CHEBI:30412",
                "CHEBI:30417", "CHEBI:30433", "CHEBI:30439", "CHEBI:30475", "CHEBI:30476",
                "CHEBI:30502", "CHEBI:30503", "CHEBI:30511", "CHEBI:30516", "CHEBI:30517",
                "CHEBI:30549", "CHEBI:30550", "CHEBI:30556", "CHEBI:30582", "CHEBI:30584",
                "CHEBI:30585", "CHEBI:30685", "CHEBI:30686", "CHEBI:30691", "CHEBI:32991",
                "CHEBI:32992", "CHEBI:32993", "CHEBI:32994", "CHEBI:32995", "CHEBI:33002",
                "CHEBI:33003", "CHEBI:33004", "CHEBI:33005", "CHEBI:33006", "CHEBI:33007",
                "CHEBI:33008", "CHEBI:33009", "CHEBI:33010", "CHEBI:33116", "CHEBI:33316",
                "CHEBI:33398", "CHEBI:33422", "CHEBI:33423", "CHEBI:33429", "CHEBI:33467",
                "CHEBI:33469", "CHEBI:33496", "CHEBI:33500", "CHEBI:33502", "CHEBI:33503",
                "CHEBI:33504", "CHEBI:33513", "CHEBI:33515", "CHEBI:33516", "CHEBI:33974",
                "CHEBI:35104", "CHEBI:35172", "CHEBI:37130", "CHEBI:37132", "CHEBI:37136",
                "CHEBI:37137", "CHEBI:37239", "CHEBI:37254", "CHEBI:37255", "CHEBI:37264",
                "CHEBI:37286", "CHEBI:37294", "CHEBI:37317", "CHEBI:39099", "CHEBI:39123",
                "CHEBI:39124", "CHEBI:39125", "CHEBI:39127", "CHEBI:39129", "CHEBI:39132",
                "CHEBI:48775", "CHEBI:48782", "CHEBI:48828", "CHEBI:49414", "CHEBI:49415",
                "CHEBI:49423", "CHEBI:49446", "CHEBI:49468", "CHEBI:49470", "CHEBI:49482",
                "CHEBI:49496", "CHEBI:49544", "CHEBI:49547", "CHEBI:49552", "CHEBI:49588",
                "CHEBI:49591", "CHEBI:49618", "CHEBI:49650", "CHEBI:49664", "CHEBI:49701",
                "CHEBI:49704", "CHEBI:49713", "CHEBI:49746", "CHEBI:49786", "CHEBI:49789",
                "CHEBI:49807", "CHEBI:49812", "CHEBI:49832", "CHEBI:49836", "CHEBI:49847",
                "CHEBI:49862", "CHEBI:49867", "CHEBI:49890", "CHEBI:49902", "CHEBI:49920",
                "CHEBI:49948", "CHEBI:49955", "CHEBI:49962", "CHEBI:49978", "CHEBI:49980",
                "CHEBI:50086", "CHEBI:50236", "CHEBI:50237", "CHEBI:60240", "CHEBI:60252",
                "CHEBI:60253", "CHEBI:60401", "CHEBI:60871", "CHEBI:63056", "CHEBI:63062",
                "CHEBI:63063", "CHEBI:63526", "CHEBI:84043", "CHEBI:85033", "CHEBI:85544",
                "CHEBI:85545", "CHEBI:88186"
            ),
            class = "ions",
            seed_type = "CHEBI",
            stringsAsFactors = FALSE
        ),
        data.frame(   ## CHEBI atoms (children nodes of root CHEBI:33250)
            seed_id = c(
                "CHEBI:137980", "CHEBI:155827", "CHEBI:176566", "CHEBI:176570", "CHEBI:176571",
                "CHEBI:176572", "CHEBI:176573", "CHEBI:176574", "CHEBI:176575", "CHEBI:176578",
                "CHEBI:176583", "CHEBI:176584", "CHEBI:18248", "CHEBI:18291", "CHEBI:194531",
                "CHEBI:194533", "CHEBI:194535", "CHEBI:194537", "CHEBI:194539", "CHEBI:194541",
                "CHEBI:196958", "CHEBI:196959", "CHEBI:22313", "CHEBI:22314", "CHEBI:22927",
                "CHEBI:22977", "CHEBI:22984", "CHEBI:23116", "CHEBI:233500", "CHEBI:24061",
                "CHEBI:24473", "CHEBI:24859", "CHEBI:25016", "CHEBI:25107", "CHEBI:25195",
                "CHEBI:25555", "CHEBI:25585", "CHEBI:25805", "CHEBI:26216", "CHEBI:26708",
                "CHEBI:26833", "CHEBI:27007", "CHEBI:27081", "CHEBI:27214", "CHEBI:27363",
                "CHEBI:27560", "CHEBI:27563", "CHEBI:27568", "CHEBI:27573", "CHEBI:27594",
                "CHEBI:27638", "CHEBI:27698", "CHEBI:27998", "CHEBI:28073", "CHEBI:28112",
                "CHEBI:28659", "CHEBI:28685", "CHEBI:28694", "CHEBI:28984", "CHEBI:29236",
                "CHEBI:29237", "CHEBI:29238", "CHEBI:29287", "CHEBI:30145", "CHEBI:30217",
                "CHEBI:30218", "CHEBI:30219", "CHEBI:30415", "CHEBI:30430", "CHEBI:30440",
                "CHEBI:30441", "CHEBI:30452", "CHEBI:30501", "CHEBI:30512", "CHEBI:30513",
                "CHEBI:30514", "CHEBI:30682", "CHEBI:30687", "CHEBI:32594", "CHEBI:32999",
                "CHEBI:33250", "CHEBI:33300", "CHEBI:33301", "CHEBI:33303", "CHEBI:33306",
                "CHEBI:33309", "CHEBI:33310", "CHEBI:33313", "CHEBI:33314", "CHEBI:33317",
                "CHEBI:33318", "CHEBI:33319", "CHEBI:33320", "CHEBI:33321", "CHEBI:33322",
                "CHEBI:33323", "CHEBI:33324", "CHEBI:33325", "CHEBI:33330", "CHEBI:33331",
                "CHEBI:33335", "CHEBI:33336", "CHEBI:33337", "CHEBI:33340", "CHEBI:33341",
                "CHEBI:33342", "CHEBI:33343", "CHEBI:33344", "CHEBI:33345", "CHEBI:33346",
                "CHEBI:33347", "CHEBI:33348", "CHEBI:33349", "CHEBI:33350", "CHEBI:33351",
                "CHEBI:33352", "CHEBI:33353", "CHEBI:33355", "CHEBI:33356", "CHEBI:33357",
                "CHEBI:33358", "CHEBI:33359", "CHEBI:33361", "CHEBI:33362", "CHEBI:33363",
                "CHEBI:33364", "CHEBI:33365", "CHEBI:33366", "CHEBI:33367", "CHEBI:33368",
                "CHEBI:33369", "CHEBI:33371", "CHEBI:33372", "CHEBI:33373", "CHEBI:33374",
                "CHEBI:33375", "CHEBI:33376", "CHEBI:33377", "CHEBI:33379", "CHEBI:33380",
                "CHEBI:33381", "CHEBI:33382", "CHEBI:33385", "CHEBI:33386", "CHEBI:33387",
                "CHEBI:33388", "CHEBI:33389", "CHEBI:33390", "CHEBI:33391", "CHEBI:33392",
                "CHEBI:33393", "CHEBI:33394", "CHEBI:33395", "CHEBI:33396", "CHEBI:33397",
                "CHEBI:33491", "CHEBI:33492", "CHEBI:33493", "CHEBI:33517", "CHEBI:33521",
                "CHEBI:33559", "CHEBI:33560", "CHEBI:33561", "CHEBI:33562", "CHEBI:33815",
                "CHEBI:33818", "CHEBI:33819", "CHEBI:36927", "CHEBI:36928", "CHEBI:36929",
                "CHEBI:36930", "CHEBI:36931", "CHEBI:36932", "CHEBI:36933", "CHEBI:36934",
                "CHEBI:36935", "CHEBI:36936", "CHEBI:36937", "CHEBI:36938", "CHEBI:36939",
                "CHEBI:36940", "CHEBI:37003", "CHEBI:37004", "CHEBI:37340", "CHEBI:37341",
                "CHEBI:37342", "CHEBI:37343", "CHEBI:37344", "CHEBI:37345", "CHEBI:37346",
                "CHEBI:37347", "CHEBI:37348", "CHEBI:37350", "CHEBI:37351", "CHEBI:37352",
                "CHEBI:37353", "CHEBI:37354", "CHEBI:37355", "CHEBI:37356", "CHEBI:37357",
                "CHEBI:37358", "CHEBI:37359", "CHEBI:37360", "CHEBI:37361", "CHEBI:37362",
                "CHEBI:37363", "CHEBI:37364", "CHEBI:37365", "CHEBI:37366", "CHEBI:37367",
                "CHEBI:37368", "CHEBI:37369", "CHEBI:37802", "CHEBI:37803", "CHEBI:37804",
                "CHEBI:37805", "CHEBI:37968", "CHEBI:37969", "CHEBI:37970", "CHEBI:37971",
                "CHEBI:37972", "CHEBI:37973", "CHEBI:37974", "CHEBI:37975", "CHEBI:37976",
                "CHEBI:37977", "CHEBI:37978", "CHEBI:37979", "CHEBI:37980", "CHEBI:37981",
                "CHEBI:37982", "CHEBI:37983", "CHEBI:37984", "CHEBI:37985", "CHEBI:49475",
                "CHEBI:49631", "CHEBI:49637", "CHEBI:49648", "CHEBI:49666", "CHEBI:49696",
                "CHEBI:49828", "CHEBI:49882", "CHEBI:49957", "CHEBI:50076", "CHEBI:52230",
                "CHEBI:52231", "CHEBI:52232", "CHEBI:52233", "CHEBI:52234", "CHEBI:52235",
                "CHEBI:52451", "CHEBI:52452", "CHEBI:52453", "CHEBI:52457", "CHEBI:52458",
                "CHEBI:52459", "CHEBI:52460", "CHEBI:52462", "CHEBI:52619", "CHEBI:52620",
                "CHEBI:52621", "CHEBI:52622", "CHEBI:52623", "CHEBI:52624", "CHEBI:52626",
                "CHEBI:52627", "CHEBI:52631", "CHEBI:52632", "CHEBI:52633", "CHEBI:52634",
                "CHEBI:52635", "CHEBI:52636", "CHEBI:52637", "CHEBI:52743", "CHEBI:52758",
                "CHEBI:52763", "CHEBI:5631", "CHEBI:77014", "CHEBI:82623", "CHEBI:88184"
            ),
            class = "atoms",
            seed_type = "CHEBI",
            stringsAsFactors = FALSE
        )
    )

    target_types <- c("hmdb", "kegg", "chebi", "pubchem")

    extract_ids <- function(x) {
        if (is.list(x)) {
            x <- unlist(x, use.names = FALSE)
        }
        x <- as.character(x)
        unique(x[!is.na(x) & x != "" & x != "0"])
    }

    translated_rows <- list()

    for (source_type in unique(seed_df$seed_type)) {
        source_type_lower <- tolower(source_type)
        source_targets <- setdiff(target_types, source_type_lower)

        input_df <- seed_df[
            seed_df$seed_type == source_type,
            c("seed_id", "class", "seed_type"),
            drop = FALSE
        ]

        translated <- tryCatch(
            translate_ids(
                input_df,
                !!sym("seed_id") := !!sym(source_type_lower),
                !!!syms(source_targets),
                ramp = TRUE,
                expand = FALSE,
                quantify_ambiguity = FALSE,
                qualify_ambiguity = FALSE
            ),
            error = function(e) input_df
        )

        for (target in source_targets) {
            if (!(target %in% colnames(translated))) {
                next
            }

            out <- lapply(seq_len(nrow(translated)), function(i) {
                ids <- extract_ids(translated[[target]][i])
                if (length(ids) == 0L) {
                    return(NULL)
                }
                data.frame(
                    metabolite_id = ids,
                    class = translated$class[i],
                    id_type = toupper(target),
                    stringsAsFactors = FALSE
                )
            })

            out <- out[!vapply(out, is.null, logical(1))]
            if (length(out) > 0L) {
                translated_rows[[paste(source_type, target, sep = "_to_")]] <-
                    do.call(rbind, out)
            }
        }
    }

    seed_rows <- data.frame(
        metabolite_id = seed_df$seed_id,
        class = seed_df$class,
        id_type = seed_df$seed_type,
        stringsAsFactors = FALSE
    )

    if (length(translated_rows) > 0L) {
        translation_df <- do.call(rbind, translated_rows)
        exclusion_df <- rbind(seed_rows, translation_df)
    } else {
        exclusion_df <- seed_rows
    }

    exclusion_df <- unique(
        exclusion_df[, c("metabolite_id", "class", "id_type"), drop = FALSE]
    )
    rownames(exclusion_df) <- NULL
    exclusion_df
}


# Helper: normalize and validate exclusion classes.
.resolve_exclusion_classes <- function(exclude_metabolites) {
    valid_classes <- c("ions", "small_molecules", "xenobiotics", "atoms")

    if (is.null(exclude_metabolites)) {
        return(character(0))
    }

    if (!is.character(exclude_metabolites)) {
        stop(
            "`exclude_metabolites` must be NULL, \"all\", or a character vector ",
            "with values from: ", paste(valid_classes, collapse = ", "), "."
        )
    }

    if (length(exclude_metabolites) == 1L && identical(exclude_metabolites, "all")) {
        return(valid_classes)
    }

    invalid_values <- setdiff(exclude_metabolites, valid_classes)
    if (length(invalid_values) > 0L) {
        stop(
            "Invalid `exclude_metabolites` value(s): ",
            paste(invalid_values, collapse = ", "),
            ". Allowed values are: NULL, \"all\", or any of: ",
            paste(valid_classes, collapse = ", "), "."
        )
    }

    unique(exclude_metabolites)
}


# Helper: apply exclusion table by ID column and ID type.
.apply_metabolite_exclusion <- function(
    df,
    id_column,
    id_type,
    exclude_metabolites
) {
    selected_classes <- .resolve_exclusion_classes(exclude_metabolites)
    if (length(selected_classes) == 0L) {
        return(df)
    }

    if (!(id_column %in% colnames(df))) {
        stop("Column `", id_column, "` not found in data.")
    }

    exclusion_df <- get_exclusion_metabolites()
    exclude_ids <- unique(exclusion_df$metabolite_id[
        exclusion_df$class %in% selected_classes &
        exclusion_df$id_type == id_type
    ])

    if (length(exclude_ids) == 0L) {
        return(df)
    }

    df[!(as.character(df[[id_column]]) %in% exclude_ids), , drop = FALSE]
}

#' KEGG pathways
#'
#' @param exclude_metabolites Optional metabolite classes to exclude:
#'     NULL (exclude nothing), "all" (default), or any combination of
#'     c("ions", "small_molecules", "xenobiotics", "atoms").
#'
#' @return A data frame containing the KEGG pathways suitable for ORA.
#'
#' @examples
#' metsigdb_kegg()
#'
#' @importFrom OmnipathR kegg_conv kegg_link kegg_list
#' @importFrom dplyr filter inner_join mutate rename
#' @importFrom logger log_info
#' @importFrom purrr map_chr map_lgl set_names
#' @importFrom stringr str_replace str_split
#' @importFrom tidyselect everything
#' @export
metsigdb_kegg <- function(
    exclude_metabolites = "all"
) {

    # NSE vs. R CMD check workaround
    name.x <- name.y <- id_a <- compound_name <- compound_names <- pathway_name <- compound <- pathway <- NULL
    ## ------------ Create log file ----------- ##
    metaproviz_init()

    log_info("Load KEGG.")

    # install_github("saezlab/OmnipathR@devel")
    KEGG_H <- kegg_link('compound', 'pathway') %>%
    mutate(across(everything(), ~str_replace(., '^\\w+:', ''))) %>%
    set_names(c('pathway', 'compound')) %>%
    inner_join(kegg_list('pathway'), by = c(pathway = 'id')) %>%
    inner_join(kegg_list('compound'), by = c(compound = 'id')) %>%
    inner_join(
        kegg_conv('compound', 'pubchem') %>%
        mutate(across(everything(), ~str_replace(., '^\\w+:', ''))),
        by = c(compound = 'id_b')
    ) %>%
    rename(pathway_name = name.x, compound_name = name.y, pubchem = id_a) %>%
    mutate(
        compound_names = str_split(compound_name, '; '),
        compound_name = map_chr(compound_names, ~extract(.x, 1L))
    ) %>%
    rename(
        term = pathway_name,
        Metabolite = compound_name,
        MetaboliteID = compound,
        Description = pathway
    )

    KEGG_H <- .apply_metabolite_exclusion(
        df = KEGG_H,
        id_column = "MetaboliteID",
        id_type = "KEGG",
        exclude_metabolites = exclude_metabolites
    )

    return(KEGG_H)

}


#
# Load RaMP prior knowledge
#

#' Metabolite chemical classes from RaMP DB
#'
#' @param version \emph{Optional: } Version of the RaMP database loaded from OmniPathR.
#'     \strong{default: "2.5.4"}
#' @param save_table \emph{Optional: } File types for the analysis results are: "csv",
#'     "xlsx", "txt". \strong{Default = "csv"}
#' @param path {Optional:} String which is added to the resulting folder name
#'     \strong{default: NULL}
#' @param exclude_metabolites Optional metabolite classes to exclude:
#'     NULL (exclude nothing), "all" (default), or any combination of
#'     c("ions", "small_molecules", "xenobiotics", "atoms").
#'
#' @return A data frame containing the Prior Knowledge.
#'
#' @examples
#' ChemicalClass <- metsigdb_chemicalclass()
#'
#' @importFrom OmnipathR ramp_table
#' @importFrom dplyr filter group_by mutate select summarise
#' @importFrom rappdirs user_cache_dir
#' @importFrom stringr str_remove str_starts
#' @importFrom tidyr pivot_wider
#' @export
metsigdb_chemicalclass <- function(
    version = "2.5.4",
    save_table = "csv",
    path = NULL,
    exclude_metabolites = "all"
) {

    # NSE vs. R CMD check workaround
    class_source_id <- chem_source_id <- class_level_name <- class_name <- common_name <- ClassyFire_class <- ClassyFire_super_class <- ClassyFire_sub_class <- NULL
    ## ------------ Create log file ----------- ##
    metaproviz_init()

    ## ------------ folder ----------- ##
    if (!is.null(save_table)) {
        folder <- save_path(
            folder_name = "PK",
            path = path
        )

        SubfolderPK <- file.path(folder, "LoadedPK")
    if (!dir.exists(SubfolderPK)) {dir.create(SubfolderPK)}

        Subfolder <- file.path(SubfolderPK, "MetaboliteSet")
    if (!dir.exists(Subfolder)) {dir.create(Subfolder)}
    }


    ##########################################################################
    # Get the directory and filepath of cache results of R
    directory <- user_cache_dir()  # get chache directory
    File_path <- paste(directory, "/RaMP-ChemicalClass_Metabolite.rds", sep = "")

    # First we will check the users chache directory and weather there are rds
    # files with KEGG_pathways already:
    if (file.exists(File_path)) {
        HMDB_ChemicalClass <- readRDS(File_path)
        message("Cached file loaded from: ", File_path)
    } else {  # load from OmniPath
    # Get RaMP via OmnipathR and extract ClassyFire classes
    Structure <- ramp_table( "metabolite_class", version = version)
    Class <- ramp_table( "chem_props", version = version)

    HMDB_ChemicalClass <-
        merge(
            Structure,
            Class[, c(seq_len(3), 10)],
            by = "ramp_id",
            all.x = TRUE
        ) %>%
    filter(str_starts(class_source_id, "hmdb:")) %>%  # Select HMDB only!
    filter(str_starts(chem_source_id, "hmdb:")) %>%  # Select HMDB only!
    select(-c("chem_data_source", "chem_source_id")) %>%
    pivot_wider(
        names_from = class_level_name,
        # Use class_level_name as the new column names
        values_from = class_name,
        # Use class_name as the values for the new columns
        # Combine duplicate values
        values_fn = list(class_name = ~paste(unique(.), collapse = ", "))
    ) %>%
    group_by(across(-common_name)) %>%
    summarise(
        common_name = paste(unique(common_name), collapse = "; "),
        # Combine all common names into one
        .groups = "drop"  # Ungroup after summarising
    ) %>%
    # Remove 'hmdb:' prefix
    mutate(class_source_id = str_remove(class_source_id, "^hmdb:")) %>%
    select(
        class_source_id,
        common_name,
        ClassyFire_class,
        ClassyFire_super_class,
        ClassyFire_sub_class
    )  # Reorder columns

    # Save the results as an RDS file in the Cache directory of R
    if (!dir.exists(directory)) {dir.create(directory)}
        saveRDS(
            HMDB_ChemicalClass,
            file = paste(directory, "/RaMP-ChemicalClass_Metabolite.rds", sep = "")
        )

    }

    HMDB_ChemicalClass <- .apply_metabolite_exclusion(
        df = HMDB_ChemicalClass,
        id_column = "class_source_id",
        id_type = "HMDB",
        exclude_metabolites = exclude_metabolites
    )

    # # -------------- Save and return
    DF_List <- list("ChemicalClass_MetabSet" = HMDB_ChemicalClass)
    # This needs to be a list, also for single comparisons
    save_res(inputlist_df = DF_List,
        inputlist_plot = NULL,
        save_table = save_table,
        save_plot = NULL,
        path = Subfolder,
        file_name = "ChemicalClass",
        core = FALSE,
        print_plot = FALSE
    )

    return(HMDB_ChemicalClass)

}


#
# Get Metabolite Pathways using Cosmos prior knowledge
#

#' Create metabolite sets from existing genesets
#'
#' Gene to metabolite translation is based on mappings in Recon-3D (cosmosR).
#'
#' @param input_pk dataframe with two columns for source (=term) and Target (=gene), e.g.
#'     Hallmarks.
#' @param metadata_info \emph{Optional: }  Column name of Target in input_pk. \strong{Default =
#'     c(Target="gene")}
#' @param pk_name \emph{Optional: } Name of the prior knowledge resource. \strong{default:
#'     NULL}
#' @param save_table \emph{Optional: } File types for the analysis results are: "csv",
#'     "xlsx", "txt". \strong{Default = "csv"}
#' @param path {Optional:} String which is added to the resulting folder name
#'     \strong{default: NULL}
#' @param exclude_metabolites Optional metabolite classes to exclude:
#'     NULL (exclude nothing), "all" (default), or any combination of
#'     c("ions", "small_molecules", "xenobiotics", "atoms").
#'
#' @return List of two data frames: "GeneMetabSet" and "MetabSet".
#'
#' @examples
#' data(hallmarks)
#' make_gene_metab_set(hallmarks)
#'
#' @importFrom dplyr rename
#' @importFrom logger log_info
#' @importFrom utils data
#' @importFrom cosmosR default_CARNIVAL_options
#' @export
make_gene_metab_set <- function(
    input_pk,
    metadata_info = c(Target = "gene"),
    pk_name = NULL,
    save_table = "csv",
    path = NULL,
    exclude_metabolites = "all"
) {

    # NSE vs. R CMD check workaround
    feature <- NULL

    ## ------------ Create log file ----------- ##
    metaproviz_init()

    log_info("make_gene_metab_set.")

    ## ------------ Check Input files ----------- ##
    # 1. The input data:
    if (!is.data.frame(input_pk)) {
        msg <- paste0(
            "`input_pk` must be of class data.frame with columns for source (= ",
            "term) and Target (= gene). Please check your input"
        )
        log_error(msg)
        stop(msg)
    }
    # 2. Target:
    if ("Target" %in% names(metadata_info)) {
        if (!(metadata_info[["Target"]] %in% colnames(input_pk))) {
            stop(
                "The ",
                metadata_info[["Target"]],
                " column selected as Conditions in metadata_info was not found in input_pk. Please check your input."
            )
        }
    } else {
        stop("Please provide a column name for the Target in metadata_info.")
    }

    if (is.null(pk_name)) {
        pk_name <- "GeneMetabSet"
    }

    ## ------------ folder ----------- ##
    if (!is.null(save_table)) {
        ## in case the user wants to save the results -> save_table is not NULL, folder is created
        folder <- save_path(
            folder_name = "PK",
            path = path
        )
    } else if (is.null(path)) {
        ## in case the user does NOT to save, i.e. save_table == NULL
        ## but no path is specified, i.e. path == NULL (default)
        folder <- file.path(getwd(), "MetaProViz_Results", "PK")
    } else {
        ## in case the user does NOT to save, i.e. save_table == NULL
        ## but a path IS specified, i.e. path is not NULL
        folder <- file.path(path, "PK")
    }

    SubfolderPK <- file.path(folder, "LoadedPK")
    if (!dir.exists(SubfolderPK)) dir.create(SubfolderPK, recursive = TRUE)

    Subfolder <- file.path(SubfolderPK, "MetaboliteSet")
    if (!dir.exists(Subfolder)) dir.create(Subfolder, recursive = TRUE)



    ##########################################################################
    # # -------------- Cosmos PKN

    # R CMD check workaround:
    nothing <- default_CARNIVAL_options("cbc")

    # load the network from cosmos
    data("meta_network", package = "cosmosR", envir = environment())
    meta_network <- get("meta_network", envir = environment())
    meta_network <- meta_network[which(meta_network$source != meta_network$target), ]

    # adapt to our needs extracting the metabolites:
    # extract entries with metabolites in source or Target
    meta_network_metabs <- meta_network[grepl("Metab__", meta_network$source) |
    grepl("Metab__HMDB", meta_network$target), -2]
    # extract entries with genes in source or Target
    meta_network_metabs <- meta_network_metabs[grepl("Gene", meta_network_metabs$source) |
    grepl("Gene", meta_network_metabs$target), ]

    # Get reactant and product
    meta_network_metabs_reactant <- meta_network_metabs[grepl("Metab__HMDB", meta_network_metabs$source), ] %>% rename("metab" = 1, "gene" = 2)
    meta_network_metabs_products <- meta_network_metabs[grepl("Metab__HMDB", meta_network_metabs$target), ] %>% rename("gene" = 1, "metab" = 2)

    meta_network_metabs <- as.data.frame(rbind(meta_network_metabs_reactant, meta_network_metabs_products))
    meta_network_metabs$gene <- gsub("Gene.*__", "", meta_network_metabs$gene)
    meta_network_metabs$metab <- gsub("_[a-z]$", "", meta_network_metabs$metab)
    meta_network_metabs$metab <- gsub("Metab__", "", meta_network_metabs$metab)

    # # --------------metalinks transporters
    # Add metalinks transporters to Cosmos PKN

    # # -------------- Combine with input_pk
    # add pathway names --> File that can be used for metabolite pathway analysis
    MetabSet <-
        merge(
            meta_network_metabs,
            input_pk,
            by.x = "gene",
            by.y = metadata_info[["Target"]]
        )

    # combine with pathways --> File that can be used for combined pathway analysis (metabolites and gene t.vals)
    GeneMetabSet <- unique(as.data.frame(rbind(input_pk %>% dplyr::rename("feature" = metadata_info[["Target"]]), MetabSet[, -1] %>% dplyr::rename("feature" = 1))))

    # Apply metabolite exclusion only on metabolite rows (genes are kept unchanged).
    gene_rows <- GeneMetabSet[!grepl("^HMDB", GeneMetabSet$feature), , drop = FALSE]
    metab_rows <- GeneMetabSet[grepl("^HMDB", GeneMetabSet$feature), , drop = FALSE]

    if (nrow(metab_rows) > 0L) {
        metab_rows <- .apply_metabolite_exclusion(
            df = metab_rows,
            id_column = "feature",
            id_type = "HMDB",
            exclude_metabolites = exclude_metabolites
        )
    }

    GeneMetabSet <- unique(rbind(gene_rows, metab_rows))

    # # ------------ Select metabolites only
    MetabSet <- GeneMetabSet %>%
    filter(grepl("^HMDB", feature))

    # # -------------- Save and return
    DF_List <- list(
        "GeneMetabSet" = GeneMetabSet,
        "MetabSet" = MetabSet
    )
    # This needs to be a list, also for single comparisons
    save_res(inputlist_df = DF_List,
        inputlist_plot = NULL,
        save_table = save_table,
        save_plot = NULL,
        path = Subfolder,
        file_name = pk_name,
        core = FALSE,
        print_plot = FALSE
    )

    return(invisible(DF_List))
}


#
# Load MetaLinksDB prior knowledge
#

#' Annotated metabolite-protein interactions from MetalinksDB
#'
#' @param types Desired edge types. Options are: "lr", "pd", where 'lr' stands for
#'     'ligand-receptor' and 'pd' stands for
#'     'production-degradation'.\strong{default: NULL}
#' @param cell_location Desired metabolite cell locations. Pass selection using c("Select1",
#'     "Select2", "Selectn"). View options setting "?". Options are:
#'     "Cytoplasm", "Endoplasmic reticulum", "Extracellular", "Lysosome" ,
#'     "Mitochondria", "Peroxisome", "Membrane", "Nucleus", "Golgi apparatus" ,
#'     "Inner mitochondrial membrane". \strong{default: NULL}
#' @param tissue_location Desired metabolite tissue locations. Pass selection using c("Select1",
#'     "Select2", "Selectn"). View options setting "?". Options are:
#'     "Placenta", "Adipose Tissue","Bladder", "Brain", "Epidermis","Kidney",
#'     "Liver", "Neuron", "Pancreas", "Prostate", "Skeletal Muscle", "Spleen",
#'     "Testis", "Thyroid Gland", "Adrenal Medulla",
#'     "Erythrocyte","Fibroblasts", "Intestine", "Ovary", "Platelet", "All
#'     Tissues", "Semen", "Adrenal Gland", "Adrenal Cortex", "Heart", "Lung",
#'     "Hair", "Eye Lens", "Leukocyte", Retina", "Smooth Muscle", "Gall
#'     Bladder", "Bile",  "Bone Marrow", "Blood", "Basal Ganglia", "Cartilage".
#'     \strong{default: NULL}
#' @param biospecimen_location Desired metabolite biospecimen locations.Pass selection using
#'     c("Select1", "Select2", "Selectn").View options setting "?".  "Blood",
#'     "Feces", "Saliva", "Sweat", "Urine", "Breast Milk", "Cellular
#'     Cytoplasm", "Cerebrospinal Fluid (CSF)", "Amniotic Fluid" , "Aqueous
#'     Humour", "Ascites Fluid", "Lymph", "Tears", "Breath", "Bile", "Semen",
#'     "Pericardial Effusion".\strong{default: NULL}
#' @param disease Desired metabolite diseases.Pass selection using c("Select1", "Select2",
#'     "Selectn"). View options setting "?". \strong{default: NULL}
#' @param pathway Desired metabolite pathways.Pass selection using c("Select1", "Select2",
#'     "Selectn"). View options setting "?".\strong{default: NULL}
#' @param hmdb_ids Desired HMDB IDs.Pass selection using c("Select1", "Select2",
#'     "Selectn"). View options setting "?".\strong{default: NULL}
#' @param uniprot_ids Desired UniProt IDs.Pass selection using c("Select1", "Select2",
#'     "Selectn"). View options setting "?".\strong{default: NULL}
#' @param save_table \emph{Optional: } File types for the analysis results are: "csv",
#'     "xlsx", "txt". \strong{Default = "csv"}
#' @param path \emph{Optional:} Path to the folder the results should be saved at.
#'     \strong{default: NULL}
#' @param exclude_metabolites Optional metabolite classes to exclude:
#'     NULL (exclude nothing), "all" (default), or any combination of
#'     c("ions", "small_molecules", "xenobiotics", "atoms").
#'
#' @return A data frame of metabolite-protein interactions from MetalinksDB.
#'
#' @examples
#' metsigdb_metalinks()
#'
#' @importFrom DBI dbListTables dbDisconnect dbGetQuery
#' @importFrom OmnipathR metalinksdb_sqlite
#' @importFrom logger log_info
#' @importFrom dplyr mutate case_when mutate
#' @importFrom stringr str_replace_all
#' @export
metsigdb_metalinks <- function(
    types = NULL,
    cell_location = NULL,
    tissue_location = NULL,
    biospecimen_location = NULL,
    disease = NULL,
    pathway = NULL,
    hmdb_ids = NULL,
    uniprot_ids = NULL,
    save_table = "csv",
    path = NULL,
    exclude_metabolites = "all"
) {

    # NSE vs. R CMD check workaround
    combined_score <- experiment_score <- mor <- protein_type <- source <- transport_direction <- type <- NULL
    ## ------------ Create log file ----------- ##
    metaproviz_init()

    log_info("MetaLinksDB.")


    # ------------------------------------------------------------------
    if (!(any(c(types, cell_location, tissue_location, biospecimen_location, disease, pathway, hmdb_ids, uniprot_ids) == "?"))) {
    # Check Input parameters


    }

    # Python version enables the user to add their own link to the database dump (probably to obtain a specific version. Lets check how the link was generated and see if it would make sense for us to do the same.)
    # --> At the moment arbitrary!
    # We could provide the user the ability to point to their own path were they already dumpled/stored QA version of metalinks they like to use!

    ## ------------ folder ----------- ##
    if (!is.null(save_table)) {
        folder <- save_path(
            folder_name = "PK",
            path = path
        )

        SubfolderPK <- file.path(folder, "LoadedPK")
    if (!dir.exists(SubfolderPK)) {dir.create(SubfolderPK)}

        Subfolder <- file.path(SubfolderPK, "MetaboliteSet")
    if (!dir.exists(Subfolder)) {dir.create(Subfolder)}
    }

    con <- metalinksdb_sqlite()
    on.exit(dbDisconnect(con))
    # ------------------------------------------------------------------
    # Query the database for a specific tables
    tables <- dbListTables(con)

    TablesList <- list()
    for (table in tables) {
    query <- paste("SELECT * FROM", table)
    data <- dbGetQuery(con, query)
    TablesList[[table]] <- data
    }

    # Close the connection

    MetalinksDB <- TablesList[["edges"]]  # extract the edges table
    # ------------------------------------------------------------------
    # Answer questions about the database
    # if any parameter is ? then return the data
    if (any(c(types, cell_location, tissue_location, biospecimen_location, disease, pathway, hmdb_ids, uniprot_ids) == "?")) {
        Questions <- which(c(types, cell_location, tissue_location, biospecimen_location, disease, pathway, hmdb_ids, uniprot_ids) == "?")
    # Check tables where the user has questions
    if (length(Questions)>0) {
        for (i in Questions) {
        if (i == 1) {
            log_info("Types:")
            log_info(unique(MetalinksDB$type))
        }
        if (i == 2) {
            log_info("Cell Location:")
            CellLocation <- TablesList[["cell_location"]]
            log_info(unique(CellLocation$cell_location))
        }
        if (i == 3) {
            log_info("Tissue Location:")
            TissueLocation <- TablesList[["tissue_location"]]
            log_info(unique(TissueLocation$tissue_location))
        }
        if (i == 4) {
            log_info("Biospecimen Location:")
            BiospecimenLocation <- TablesList[["biospecimen_location"]]
            log_info(unique(BiospecimenLocation$biospecimen_location))
        }
        if (i == 5) {
            log_info("Disease:")
            Disease <- TablesList[["disease"]]
            log_info(unique(Disease$disease))
        }
        if (i == 6) {
            log_info("Pathway:")
            Pathway <- TablesList[["pathway"]]
            log_info(unique(Pathway$pathway))
        }
        if (i == 7) {
            log_info("HMDB IDs:")
            log_info(unique(MetalinksDB$hmdb))
        }
        if (i == 8) {
            log_info("UniProt IDs:")
            log_info(unique(MetalinksDB$uniprot))
        }
        }
    }
    msg <- paste0(
        "No result is returned unless correct options for your selections are ",
        "used. `?` is not a valid option, but only returns you the list of ",
        "options."
    )
    log_info(msg)
    message(msg)
    return()
    }

    # ------------------------------------------------------------------
    # Extract specific connections based on parameter settings. If any parameter is not NULL, filter the data:
    ## types
    if (!is.null(types)) {
        MetalinksDB <- MetalinksDB[MetalinksDB$type %in% types, ]
    }

    ## cell_location
    if (!is.null(cell_location)) {
        CellLocation <- TablesList[["cell_location"]]
        # Filter the cell location
        CellLocation <- CellLocation[CellLocation$cell_location %in% cell_location, ]

    # Get unique HMDB IDs
        CellLocation_HMDB <- unique(CellLocation$hmdb)

    # Only keep selected HMDB IDs
        MetalinksDB <- MetalinksDB[MetalinksDB$hmdb %in%  CellLocation_HMDB, ]
    }

    ## tissue_location
    if (!is.null(tissue_location)) {  # "All Tissues"?
        TissueLocation <- TablesList[["tissue_location"]]
        # Filter the tissue location
        TissueLocation <- TissueLocation[TissueLocation$tissue_location %in% tissue_location, ]

    # Get unique HMDB IDs
        TissueLocation_HMDB <- unique(TissueLocation$hmdb)

    # Only keep selected HMDB IDs
        MetalinksDB <- MetalinksDB[MetalinksDB$hmdb %in%  TissueLocation_HMDB, ]
    }

    ## biospecimen_location
    if (!is.null(biospecimen_location)) {
        BiospecimenLocation <- TablesList[["biospecimen_location"]]
        # Filter the biospecimen location
        BiospecimenLocation <- BiospecimenLocation[BiospecimenLocation$biospecimen_location %in% biospecimen_location, ]

    # Get unique HMDB IDs
        BiospecimenLocation_HMDB <- unique(BiospecimenLocation$hmdb)

    # Only keep selected HMDB IDs
        MetalinksDB <- MetalinksDB[MetalinksDB$hmdb %in%  BiospecimenLocation_HMDB, ]
    }

    ## disease
    if (!is.null(disease)) {
        Disease <- TablesList[["disease"]]
        Disease <- Disease[Disease$disease %in% disease, ]  # Filter the disease

    # Get unique HMDB IDs
        Disease_HMDB <- unique(Disease$hmdb)

    # Only keep selected HMDB IDs
        MetalinksDB <- MetalinksDB[MetalinksDB$hmdb %in%  Disease_HMDB, ]
    }

    ## pathway
    if (!is.null(pathway)) {
        Pathway <- TablesList[["pathway"]]
        Pathway <- Pathway[Pathway$pathway %in% pathway, ]  # Filter the pathway

    # Get unique HMDB IDs
        Pathway_HMDB <- unique(Pathway$hmdb)

    # Only keep selected HMDB IDs
        MetalinksDB <- MetalinksDB[MetalinksDB$hmdb %in%  Pathway_HMDB, ]
    }

    ## hmdb_ids
    if (!is.null(hmdb_ids)) {
    # Only keep selected HMDB IDs
        MetalinksDB <- MetalinksDB[MetalinksDB$hmdb %in%  hmdb_ids, ]
    }

    ## uniprot_ids
    if (!is.null(uniprot_ids)) {
    # Only keep selected UniProt IDs
        MetalinksDB <- MetalinksDB[MetalinksDB$uniprot %in%  uniprot_ids, ]
    }

    # ------------------------------------------------------------------
    # Add other ID types:
    ## Metabolite Name
    MetalinksDB <-
        merge(
            MetalinksDB,
            TablesList[["metabolites"]],
            by = "hmdb",
            all.x = TRUE
        )

    ## Gene Name
    MetalinksDB <-
        merge(
            MetalinksDB,
            TablesList[["proteins"]],
            by = "uniprot",
            all.x = TRUE
        )


    ## Rearrange columns:
    MetalinksDB <- MetalinksDB[, c(2, 10:12, 1, 13:14, 3:9)] %>%
    mutate(
        type = case_when(
        type == "lr" ~ "Ligand-Receptor",
        type == "pd" ~ "Production-Degradation",
        # this keeps the original value if it doesn't match any condition
        TRUE ~ type
    )
    ) %>%
    mutate(
        mode_of_regulation = case_when(
        mor == -1 ~ "Inhibiting",
        mor == 1 ~ "Activating",
        mor == 0 ~ "Binding",
        # this keeps the original value if it doesn't match any condition
        TRUE ~ as.character(mor)
    )
    ) %>%
    mutate(
        protein_type_clean = str_replace_all(protein_type, '"', ''),
        receptor_class = case_when(
            protein_type_clean == "gpcr" ~ "GPCR",
            protein_type_clean == "lgic" ~ "Ligand-gated ion channel",
            protein_type_clean == "vgic" ~ "Voltage-gated ion channel",
            protein_type_clean == "catalytic_receptor" ~ "Catalytic receptor",
            protein_type_clean == "nhr" ~ "Nuclear hormone receptor",
            TRUE ~ NA_character_
        ),
        interaction_family = case_when(
            protein_type_clean == "transporter" | !is.na(transport_direction) ~ "Transporter-metabolite",
            type == "Ligand-Receptor" | !is.na(receptor_class) ~ "Receptor-metabolite",
            protein_type_clean == "enzyme" ~ "Enzyme-metabolite",
            TRUE ~ "Other protein-metabolite"
        ),
        interaction_mechanism = case_when(
            interaction_family == "Transporter-metabolite" ~ "Transport",
            interaction_family == "Receptor-metabolite" ~ "Ligand-receptor signaling",
            interaction_family == "Enzyme-metabolite" | type == "Production-Degradation" ~ "Metabolic conversion/turnover",
            TRUE ~ "Other/unspecified"
        ),
        regulation_polarity = case_when(
            mor == 1 ~ "Positive (activating)",
            mor == -1 ~ "Negative (inhibiting)",
            mor == 0 ~ "Neutral (binding)",
            TRUE ~ "Unknown"
        ),
        transport_direction_label = case_when(
            transport_direction == "in" ~ "Import/Uptake",
            transport_direction == "out" ~ "Export/Secretion",
            is.na(transport_direction) ~ "Not specified",
            TRUE ~ "Other"
        ),
        transport_mode = case_when(
            interaction_family == "Transporter-metabolite" & transport_direction == "in" ~ "Uptake transporter",
            interaction_family == "Transporter-metabolite" & transport_direction == "out" ~ "Efflux transporter",
            interaction_family == "Transporter-metabolite" ~ "Transporter (direction unspecified)",
            TRUE ~ NA_character_
        ),
        evidence_class = case_when(
            source %in% c("rhea", "hmr", "recon") ~ "Metabolic knowledgebase",
            source %in% c("CellPhoneDB", "scConnect", "Cellinker", "NeuronChat") ~ "Cell-cell communication resource",
            source == "Stitch" ~ "Chemical-protein interaction resource",
            is.na(source) ~ "Unknown",
            TRUE ~ "Other"
        ),
        experiment_evidence_present = ifelse(!is.na(experiment_score) & experiment_score > 0, "Yes", "No/Unknown"),
        combined_confidence_tier = case_when(
            is.na(combined_score) ~ "Unknown",
            combined_score >= 900 ~ "Very high",
            combined_score >= 700 ~ "High",
            combined_score >= 400 ~ "Medium",
            TRUE ~ "Low"
        ),
        interaction_detail = case_when(
            interaction_family == "Transporter-metabolite" ~ paste0(transport_mode, "; ", regulation_polarity),
            interaction_family == "Receptor-metabolite" ~ paste0(ifelse(is.na(receptor_class), "Receptor", receptor_class), "; ", regulation_polarity),
            interaction_family == "Enzyme-metabolite" ~ paste0("Enzymatic link; ", regulation_polarity),
            TRUE ~ paste0("Other link; ", regulation_polarity)
        )
    )
    # --------------------------------------------------------
    # Remove metabolites that are not detectable by mass spectrometry


    # ------------------------------------------------------------------
    # Decide on useful selections term-metabolite for MetaProViz.
    # MetalinksDB_Pathways <- merge(MetalinksDB[, c(1:3)], TablesList[["pathway"]], by="hmdb", all.x=TRUE)

    MetalinksDB %<>%
        mutate(
            term_specific = ifelse(
            is.na(protein_type),
            NA,
            sprintf(
            '%s_%s',
            str_replace_all(protein_type, '"', ''),
            type
        )
        )
        )
    MetalinksDB <- .apply_metabolite_exclusion(
        df = MetalinksDB,
        id_column = "hmdb",
        id_type = "HMDB",
        exclude_metabolites = exclude_metabolites
    )

    # ------------------------------------------------------------------
    # Save results in folder
    # # -------------- Save and return
    # This needs to be a list, also for single comparisons
    save_res(inputlist_df = list(MetalinksDB),
        inputlist_plot = NULL,
        save_table = save_table,
        save_plot = NULL,
        path = Subfolder,
        file_name = "MetaLinksDB",
        core = FALSE,
        print_plot = FALSE
    )

    return(MetalinksDB)

}



#' Retrieve Reactome metabolite sets suitable for ORA.
#' 
#' Queries the OmniPath resource through OmniPathR to obtail Reactome pathway
#' level metabolite sets.
#'
#' @param species String. Optionally specify pathways to query from a species
#' via full name or three letter code. Default = "Homo sapiens". NULL for all species.
#' @param pathway_ids String vector. Optionally provide pathway_ids to query. Default NULL to query all pathways.
#' @param out_path String. Optionally save results as csv into out_path. Default NULL.
#' @param exclude_metabolites Optional metabolite classes to exclude:
#'     NULL (exclude nothing), "all" (default), or any combination of
#'     c("ions", "small_molecules", "xenobiotics", "atoms").
#' 
#' @return 
#' A tibble in long format containing one row per metabolite for the Reactome pathways.
#' 
#' @examples 
#' \dontrun{
#' df <- metsigdb_reactome()
#' head(df)
#' }
#' 
#' @importFrom OmnipathR reactome_chebi
#' @importFrom magrittr %>%
#' @importFrom tidyr separate_rows
#' @importFrom dplyr mutate rename
#' 
#' @export
metsigdb_reactome <- function(
    species = "Homo sapiens",
    pathway_ids = NULL,
    out_path = NULL,
    exclude_metabolites = "all"
){

    # NSE vs. R CMD check workaround
    chebi_ids <- pathway_df <- pathway_df_long <- NULL

    pathway_df <- OmnipathR::reactome_chebi(
        organism = species,
        pathway_ids = pathway_ids,
        out_path = out_path
    )

    pathway_df_long <- pathway_df %>%
        separate_rows(chebi_ids, sep = ", ") %>%
        mutate(chebi_ids = trimws(chebi_ids)) %>%
        rename(chebi_id = chebi_ids)

    pathway_df_long <- .apply_metabolite_exclusion(
        df = pathway_df_long,
        id_column = "chebi_id",
        id_type = "CHEBI",
        exclude_metabolites = exclude_metabolites
    )

    pathway_df_long

}


#' Retrieve WikiPathways metabolite mapping suitable for ORA.
#'
#' Retrieves pathway to metabolite mappings from WikiPathways (via
#' `wikipathways_metabolites_sparql()`) via OmnipathR and returns a long-format table with
#' one metabolite identifier per row.
#'
#' @param species Character. Species name. Default is `"Homo sapiens"`.
#' @param exclude_metabolites Optional metabolite classes to exclude:
#'     NULL (exclude nothing), "all" (default), or any combination of
#'     c("ions", "small_molecules", "xenobiotics", "atoms").
#'
#' @return A tibble in long format with columns `pathway_id`, `pathway_name`,
#'   `pathway_url`, `n_metabolites_in_pathway`, and `metabolite_id`.
#'
#' @importFrom dplyr mutate rename
#' @importFrom tidyr separate_rows
#' @importFrom stringr str_trim
#' @importFrom magrittr %>%
#' @importFrom OmnipathR wikipathways_metabolites_sparql
#'
#' @export
metsigdb_wikipathways <- function(
    species = "Homo sapiens",
    exclude_metabolites = "all"
) {

    # NSE vs. R CMD check workaround
    metabolites <- pathway_df <- pathway_df_long <- NULL

    pathway_df <- wikipathways_metabolites_sparql(organism = species)

    pathway_df_long <- pathway_df %>%
        tidyr::separate_rows(metabolites, sep = ";\\s*") %>%
        dplyr::mutate(metabolites = stringr::str_trim(metabolites)) %>%
        dplyr::rename(metabolite_id = metabolites)

    pathway_df_long <- .apply_metabolite_exclusion(
        df = pathway_df_long,
        id_column = "metabolite_id",
        id_type = "CHEBI",
        exclude_metabolites = exclude_metabolites
    )

    pathway_df_long

}





#' Retrieve MACDB metabolite-cancer associations.
#'
#' Retrieves metabolite-cancer associations from MACDB via OmnipathR and
#' summarizes them to unique metabolite-cancer associations (by `term` and
#' `Metabolite_PubchemID`).
#'
#' @param exclude_metabolites Optional metabolite classes to exclude:
#'     NULL (exclude nothing), "all" (default), or any combination of
#'     c("ions", "small_molecules", "xenobiotics", "atoms").
#'
#' @return A data frame with one row per unique metabolite-cancer association,
#'   collapsed metadata columns, and summary metrics (`evidence_count`,
#'   `significance_count`, `association_score`).
#'
#' @importFrom OmnipathR macdb_metabolite_cancer_associations
#' @importFrom dplyr mutate
#'
#' @export
metsigdb_macdb <- function(
    exclude_metabolites = "all"
) {

    # NSE vs. R CMD check workaround
    metabolite_pubchem_cid <- metabolite_name <- study_cancer_subtype <- study_cancer_type <- NULL

    macdb_df <- OmnipathR::macdb_metabolite_cancer_associations()

    required_columns <- c(
        "metabolite_pubchem_cid",
        "metabolite_name",
        "study_cancer_type",
        "study_cancer_subtype",
        "study_pubmed_id",
        "association_pvalue"
    )
    missing_columns <- setdiff(required_columns, colnames(macdb_df))
    if (length(missing_columns) > 0L) {
        stop(
            "Unexpected MACDB format. Missing required columns: ",
            paste(missing_columns, collapse = ", ")
        )
    }

    macdb_df <- macdb_df %>%
        dplyr::mutate(
            term = ifelse(
                is.na(study_cancer_subtype) | study_cancer_subtype == "",
                as.character(study_cancer_type),
                as.character(study_cancer_subtype)
            ),
            Metabolite = metabolite_name,
            Metabolite_PubchemID = as.character(metabolite_pubchem_cid),
            Description = study_cancer_type
        )

    macdb_df <- .apply_metabolite_exclusion(
        df = macdb_df,
        id_column = "Metabolite_PubchemID",
        id_type = "PUBCHEM",
        exclude_metabolites = exclude_metabolites
    )

    keep_columns <- c(
        "term",
        "Description",
        "Metabolite",
        "Metabolite_PubchemID",
        "Cohort_id",
        "study_trait_onto_id",
        "study_cancer_doid",
        "study_tissue",
        "study_pubmed_id",
        "association_pvalue",
        "association_log2fc",
        "association_case_concentration",
        "association_control_concentration",
        "has_pValue",
        "has_log2FC",
        "has_caseConcentration",
        "has_controlConcentration",
        "has_bothConcentrations"
    )

    keep_columns <- keep_columns[keep_columns %in% colnames(macdb_df)]
    macdb_df <- macdb_df[, keep_columns, drop = FALSE]

    .summarize_macdb_associations(macdb_df)

}


.summarize_macdb_associations <- function(macdb_df) {

    # NSE vs. R CMD check workaround
    term <- Metabolite_PubchemID <- study_pubmed_id <- Cohort_id <- association_pvalue <-
        evidence_count <- significance_count <- association_score <-
        study_significant <- NULL

    key_columns <- c("term", "Metabolite_PubchemID")
    missing_key_columns <- setdiff(key_columns, colnames(macdb_df))
    if (length(missing_key_columns) > 0L) {
        stop(
            "Missing key column(s) for MACDB summarization: ",
            paste(missing_key_columns, collapse = ", ")
        )
    }
    collapse_columns <- setdiff(colnames(macdb_df), key_columns)
    summarized_df <- macdb_df %>%
        dplyr::group_by(term, Metabolite_PubchemID) %>%
        dplyr::summarise(
            dplyr::across(
                dplyr::all_of(collapse_columns),
                function(x) {
                    x <- as.character(x)
                    x <- x[!is.na(x) & nzchar(x)]
                    if (length(x) == 0L) {
                        return(NA_character_)
                    }
                    x <- x[!duplicated(x)]
                    paste(x, collapse = "; ")
                }
            ),
            .groups = "drop"
        )

    per_study <- macdb_df %>%
        dplyr::filter(!is.na(study_pubmed_id)) %>%
        dplyr::group_by(term, Metabolite_PubchemID, study_pubmed_id) %>%
        dplyr::summarise(
            study_significant = any(!is.na(association_pvalue) & association_pvalue <= 0.05),
            .groups = "drop"
        )

    counts_df <- per_study %>%
        dplyr::group_by(term, Metabolite_PubchemID) %>%
        dplyr::summarise(
            evidence_count = dplyr::n_distinct(study_pubmed_id),
            significance_count = sum(study_significant),
            .groups = "drop"
        )

    counts_df <- summarized_df %>%
        dplyr::select(term, Metabolite_PubchemID) %>%
        dplyr::left_join(counts_df, by = c("term", "Metabolite_PubchemID")) %>%
        dplyr::mutate(
            evidence_count = ifelse(is.na(evidence_count), 0L, as.integer(evidence_count)),
            significance_count = ifelse(is.na(significance_count), 0L, as.integer(significance_count)),
            association_score = evidence_count + significance_count
        )

    summarized_df %>%
        dplyr::left_join(
            counts_df,
            by = c("term", "Metabolite_PubchemID")
        )

}
