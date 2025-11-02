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

#' KEGG pathways
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
metsigdb_kegg <- function() {

    # NSE vs. R CMD check workaround
    name.x <- name.y <- id_a <- compound_name <- compound_names <- pathway_name <- compound <- pathway <- NULL
    ## ------------ Create log file ----------- ##
    metaproviz_init()

    log_info("Load KEGG.")

    # ------------------------------------------------------------------
    # Remove Metabolites
    to_remove <-
    list(
        # # # 1. Ions should be removed
        Remove_Ions = c("Calcium cation", "Potassium cation", "Sodium cation", "H+", "Cl-", "Fatty acid", "Superoxide", "H2O", "CO2", "Copper", "Fe2+", "Magnesium cation", "Fe3+", "Zinc cation", "Nickel", "NH4+"),
        # # # 2. Unspecific small molecules
        Remove_Small = c("Nitric oxide", "Hydrogen peroxide", "Superoxide", "H2O", "CO2", "Hydroxyl radical", "Ammonia", "HCO3-", "Oxygen", "Diphosphate", "Reactive oxygen species", "Nitrite", "Nitrate", "Hydrogen", "RX", "Hg")
    ) %>%
    unlist()


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
    filter(
        map_lgl(
        compound_names,
        ~intersect(.x, to_remove) %>% length() == 0
    )
    ) %>%
    rename(
        term = pathway_name,
        Metabolite = compound_name,
        MetaboliteID = compound,
        Description = pathway
    )  # Update vignettes and remove rename


    # Return into environment
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
    path = NULL
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
    path = NULL
) {

    # NSE vs. R CMD check workaround
    feature <- NULL

    ## ------------ Create log file ----------- ##
    metaproviz_init()

    log_info("make_gene_metab_set.")

    ## ------------ Check Input files ----------- ##
    # 1. The input data:
    if (!is.data.frame(input_pk)) {
        stop("`input_pk` must be of class data.frame with columns for source (= term) and Target (= gene). Please check your input")
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
        folder <- save_path(
            folder_name = "PK",
            path = path
        )
        SubfolderPK <- file.path(folder, "LoadedPK")
    } else if (!dir.exists(SubfolderPK)) {
        dir.create(SubfolderPK)
        Subfolder <- file.path(SubfolderPK, "MetaboliteSet")
    } else if (!dir.exists(Subfolder)) {
        dir.create(Subfolder)
    }


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


    # # ------------ Select metabolites only
    MetabSet <- GeneMetabSet %>%
    filter(grepl("HMDB", feature))


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
    path = NULL
) {

    # NSE vs. R CMD check workaround
    protein_type <- type <- NULL
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
    message("No result is returned unless correct options for your selections are used. `?` is not a valid option, but only returns you the list of options.")
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

#
# Helper: Load compound lists of ions, xenobiotics, cofactors (not detectable
# by LC-MS), etc
#


# This is needed to remove ions, xenobiotics, cofactors, etc. from the metabolite list for any of the functions above.

