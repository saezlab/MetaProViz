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
# Heatmap
#

#' Heatmap visualization
#'
#' @param data SummarizedExperiment (se) file including assay and colData. If se file
#'     is provided, metadata_sample is extracted from the colData of the se
#'     object. metadata_feature, if available, are extracted from the rowData.
#'     Alternatively provide a DF with unique sample identifiers as row names
#'     and metabolite numerical values in columns with metabolite identifiers
#'     as column names. Use NA for metabolites that were not detected.
#' @param metadata_info \emph{Optional: } NULL or Named vector where you can include vectors or
#'     lists for annotation c(individual_Metab=
#'     "ColumnName_metadata_feature",individual_Sample=
#'     "ColumnName_metadata_sample", color_Metab="ColumnName_metadata_feature",
#'     color_Sample= list("ColumnName_metadata_sample",
#'     "ColumnName_metadata_sample",...)).\strong{Default = NULL}
#' @param metadata_sample \emph{Optional: } Only required if you did not provide se file in
#'     parameter data. Provide DF which contains metadata information about the
#'     samples, which will be combined with your input data based on the unique
#'     sample identifiers used as rownames. \strong{Default = NULL}
#' @param metadata_feature \emph{Optional: } To provide metadata information for each metabolite.
#'     Only used if you did not provide se file in parameter data. Provide DF
#'     where the row names must match the metabolite names in the columns of
#'     the data. \strong{Default = NULL}
#' @param plot_name \emph{Optional: } String which is added to the output files of the plot
#' @param scale \emph{Optional: } String with the information for scale row, column or
#'     none. \strong{Default = row}
#' @param save_plot \emph{Optional: } Select the file type of output plots. Options are svg,
#'     pdf, png or NULL. \strong{Default = "svg"}
#' @param enforce_featurenames \emph{Optional: } If there are more than 100 features no rownames will
#'     be shown, which is due to readability. You can Enforce this by setting
#'     this parameter to TRUE. \strong{Default = FALSE}
#' @param enforce_samplenames \emph{Optional: } If there are more than 50 sampless no colnames will be
#'     shown, which is due to readability. You can Enforce this by setting this
#'     parameter to TRUE. \strong{Default = FALSE}
#' @param print_plot \emph{Optional: } print the plots to the active graphic device.
#' @param path {Optional:} String which is added to the resulting folder name
#'     \strong{default: NULL}
#'
#' @return List with two elements: Plot and Plot_Sized
#'
#' @examples
#' data(intracell_raw_se)
#' Res <- viz_heatmap(data = intracell_raw_se)
#'
#' data(intracell_raw)
#' Intra <- intracell_raw %>% tibble::column_to_rownames("Code")
#' Res <- viz_heatmap(data = Intra[, -c(1:3)])
#'
#' @importFrom ggplot2 ggplot theme element_rect
#' @importFrom grid convertUnit
#' @importFrom dplyr rename select
#' @importFrom magrittr %>% %<>% is_in not
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom logger log_trace
#' @importFrom pheatmap pheatmap
#' @export
viz_heatmap <- function(
    data,
    metadata_info = NULL,
    metadata_sample = NULL,
    metadata_feature = NULL,
    plot_name = "",
    scale = "row",
    save_plot = "svg",
    enforce_featurenames = FALSE,
    enforce_samplenames = FALSE,
    print_plot = TRUE,
    path = NULL
) {
    # # ------------ Create log file ----------- ##
    metaproviz_init()

    ## ------------- Check SummarizedExperiment file ---------- ##
    input_data <- data
    if (inherits(data, "SummarizedExperiment")) {
        log_info('Processing input SummarizedExperiment object.')
        se_list <- process_se(data)
        data <- se_list$data
        metadata_sample <- se_list$metadata_sample
        metadata_feature <- se_list$metadata_feature
    }

    # # ------------ Check Input files ----------- ##
    # HelperFunction `check_param`
    check_param(
        data = data,
        metadata_sample = metadata_sample,
        metadata_feature = metadata_feature,
        metadata_info = metadata_info,
        save_plot = save_plot,
        save_table = NULL,
        core = FALSE,
        print_plot = print_plot
    )

    # check_param` Specific
    if (!is.logical(enforce_featurenames) | !is.logical(enforce_samplenames)) {
        message <- paste0("Check input. The enforce_featurenames and enforce_samplenames value should be either =TRUE or = FALSE.")
        log_trace(paste("Error ", message, sep = ""))
        stop(message)
    }

    scale_options <- c("row", "column", "none")
    if (!(scale %in% scale_options)) {
        message <-
            paste0(
                "Check input. The selected scale option is not valid. Please select one of the folowwing: ",
                paste(scale_options, collapse = ", "),
                "."
            )
        log_trace(paste("Error ", message, sep = ""))
        stop(message)
    }

    # # ------------ Create Results output folder ----------- ##
    if (!is.null(save_plot)) {
        folder <- save_path(
            folder_name = "Heatmap",
            path = path
        )
    }

    # ####################################################
    # # -------------- Load data --------------- ##
    # what is this??
    data <- data

    # removes information about metabolites that are not included in the data
    if (!is.null(metadata_feature)) {
        metadata_feature <-
            merge(
                x = metadata_feature,
                y = as.data.frame(t(data)),
                by = 0,
                all.y = TRUE
            ) %>%
            column_to_rownames("Row.names")
        metadata_feature <- metadata_feature[, -c((ncol(metadata_feature) - nrow(data) + 1):ncol(metadata_feature))]
    }

    # removes information about samples that are not included in the data
    if (!is.null(metadata_sample)) {
        metadata_sample <-
            merge(x = metadata_sample, y = data, by = 0, all.y = TRUE) %>%
            column_to_rownames("Row.names")
        metadata_sample <- metadata_sample[, -c((ncol(metadata_sample) - ncol(data) + 1):ncol(metadata_sample))]
    }

    # # -------------- Plot --------------- ##
    if (
        !("individual_Metab" %in% names(metadata_info)) &
        "individual_Sample" %in% names(metadata_info)
    ) {
        # Ensure that groups that are assigned NAs do not cause problems:
        metadata_feature[[metadata_info[["individual_Metab"]]]] <- ifelse(is.na(metadata_feature[[metadata_info[["individual_Metab"]]]]), "NA", metadata_feature[[metadata_info[["individual_Metab"]]]])
        unique_paths <- unique(metadata_feature[[metadata_info[["individual_Metab"]]]])

        for (i in unique_paths) {  # Check pathways with 1 metabolite
            selected_path <-
                metadata_feature %>%
                filter(get(metadata_info[["individual_Metab"]]) == i)
            selected_path_metabs <- colnames(data)[colnames(data) %in% row.names(selected_path)]
            if (length(selected_path_metabs) == 1L) {
                message <-
                    paste0(
                        "The metadata group ",
                        i,
                        " includes only 1 metabolite. Heatmap cannot be made for 1 metabolite, thus it will be ignored."
                    )
                log_trace(paste("Warning ", message, sep = ""))
                warning(message)
                # Remove the pathway
                unique_paths <- unique_paths[!unique_paths %in% i]
            }
        }

        IndividualPlots <- unique_paths

        PlotList <- list()  # Empty list to store all the plots
        PlotList_adaptedGrid <- list()  # Empty list to store all the plots

        for (i in IndividualPlots) {
            selected_path <-
                metadata_feature %>%
                filter(get(metadata_info[["individual_Metab"]]) == i)
            selected_path_metabs <-
                colnames(data)[colnames(data) %in% row.names(selected_path)]
            data_path <- data %>% select(all_of(selected_path_metabs))

            # Column annotation
            col_annot_vars <- metadata_info[grepl("color_Sample", names(metadata_info))]
            col_annot <- NULL
            if (length(col_annot_vars) > 0L) {
                for (x in seq_along(col_annot_vars)) {
                    annot_sel <- col_annot_vars[[x]]
                    col_annot[x] <- metadata_sample %>%
                    select(annot_sel) %>%
                    as.data.frame()
                    names(col_annot)[x] <- annot_sel
                }
                col_annot <- as.data.frame(col_annot)
                rownames(col_annot) <- rownames(data_path)
            }

            # Row annotation
            row_annot_vars <- metadata_info[grepl("color_Metab", names(metadata_info))]
            row_annot <- NULL
            if (length(row_annot_vars) > 0L) {
                for (y in seq_along(row_annot_vars)) {
                    annot_sel <-
                        row_annot_vars[[y]]
                    row_annot[y] <- metadata_feature %>% select(all_of(annot_sel))
                    row_annot <-
                        row_annot %>%
                        as.data.frame()
                    names(row_annot)[y] <- annot_sel
                }
                row_annot <- as.data.frame(row_annot)
                rownames(row_annot) <- rownames(metadata_feature)
            }

            # Check number of features:
            Features <- as.data.frame(t(data_path))
            if (enforce_featurenames) {
                show_rownames <- TRUE
                cellheight_Feature <- 9
            } else if (nrow(Features) > 100L) {
                show_rownames <- FALSE
                cellheight_Feature <- 1
            } else {
                show_rownames <- TRUE
                cellheight_Feature <- 9
            }

            # Check number of samples
            if (enforce_samplenames) {
                show_colnames <- TRUE
                cellwidth_Sample <- 9
            } else if (nrow(data_path) > 50L) {
                show_colnames <- FALSE
                cellwidth_Sample <- 1
            } else {
                show_colnames <- TRUE
                cellwidth_Sample <- 9
            }

            # Make the plot
            if (nrow(t(data_path)) >= 2L) {
                # set.seed(1234)

                heatmap <- pheatmap(
                    t(data_path),
                    show_rownames = as.logical(show_rownames),
                    show_colnames = as.logical(show_colnames),
                    clustering_method = "complete",
                    scale = scale,
                    clustering_distance_rows = "correlation",
                    annotation_col = col_annot,
                    annotation_row = row_annot,
                    legend = TRUE,
                    cellwidth = cellwidth_Sample,
                    cellheight = cellheight_Feature,
                    fontsize_row = 10,
                    fontsize_col = 10,
                    fontsize = 9,
                    main = paste(plot_name, " Metabolites: ", i, sep = " "),
                    silent = TRUE
                )

                # # Store the plot in the 'plots' list
                cleaned_i <-
                    gsub(
                        "[[:space:], /\\\\]",
                        "-",
                        i
                    )  # removes empty spaces and replaces /, \ with -
                PlotList[[cleaned_i]] <- heatmap

                # Width and height according to Sample and metabolite number
                Plot_Sized <-
                    plot_grob_heatmap(
                        input_plot = heatmap,
                        metadata_info = metadata_info,
                        metadata_sample = metadata_sample,
                        metadata_feature = metadata_feature,
                        plot_name = cleaned_i
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
                    add(theme(panel.background = element_rect(fill = "transparent")))

                PlotList_adaptedGrid[[cleaned_i]] <- Plot_Sized

                # ----- Save
                save_res(
                    inputlist_df = NULL,
                    inputlist_plot = PlotList_adaptedGrid,
                    save_table = NULL,
                    save_plot = save_plot,
                    path = folder,
                    file_name = paste("Heatmap_", plot_name, sep = ""),
                    core = FALSE,
                    print_plot = print_plot,
                    plot_height = plot_height,
                    plot_width = plot_width,
                    plot_unit = "cm"
                )
            } else {
                message <-
                    paste0(
                        i,
                        " includes <= 2L objects and is hence not plotted."
                    )
                log_trace(paste("Message ", message, sep = ""))
                message(message)
            }
        }

        # Return if assigned:
        return(invisible(list("Plot" = PlotList, "Plot_Sized" = PlotList_adaptedGrid)))

    } else if (
        !("individual_Metab" %in% names(metadata_info)) &
        "individual_Sample" %in% names(metadata_info)
    ) {

        # Ensure that groups that are assigned NAs do not cause problems:
        metadata_sample[[metadata_info[["individual_Sample"]]]] <- ifelse(
            is.na(metadata_sample[[metadata_info[["individual_Sample"]]]]),
            "NA",
            metadata_sample[[metadata_info[["individual_Sample"]]]]
        )

        unique_paths_Sample <- unique(metadata_sample[[metadata_info[["individual_Sample"]]]])

        for (i in unique_paths_Sample) {  # Check pathways with 1 metabolite
            selected_path <-
                metadata_sample %>%
                filter(get(metadata_info[["individual_Sample"]]) == i)
            selected_path_metabs <- colnames(data)[colnames(data) %in% row.names(selected_path)]
            if (length(selected_path_metabs) == 1L) {
                message <-
                    paste0(
                        "The metadata group ",
                        i,
                        " includes only 1 metabolite. Heatmap cannot be made for 1 metabolite, thus it will be ignored."
                    )
                log_trace(paste("Warning ", message, sep = ""))
                warning(message)
                # Remove the pathway
                unique_paths_Sample <- unique_paths_Sample[!unique_paths_Sample %in% i]
            }
        }

        IndividualPlots <- unique_paths_Sample
        PlotList <- list()  # Empty list to store all the plots
        PlotList_adaptedGrid <- list()  # Empty list to store all the plots

        for (i in IndividualPlots) {
            # Select the data:
            selected_path <-
                metadata_sample %>%
                filter(get(metadata_info[["individual_Sample"]]) == i) %>%
                rownames_to_column("UniqueID")
            selected_path <-
                as.data.frame(selected_path[, 1]) %>%
                rename("UniqueID" = 1)
            data_path <-
                merge(
                    selected_path,
                    data %>% rownames_to_column("UniqueID"),
                    by = "UniqueID",
                    all.x = TRUE
                )
            data_path <- data_path %>%
            column_to_rownames("UniqueID")

            # Column annotation
            selected_metadata_sample <-
                merge(
                    selected_path,
                    metadata_sample %>% rownames_to_column("UniqueID"),
                    by = "UniqueID",
                    all.x = TRUE
                )

            col_annot_vars <- metadata_info[grepl("color_Sample", names(metadata_info))]
            col_annot <- NULL
            if (length(col_annot_vars) > 0L) {
                for (x in seq_along(col_annot_vars)) {
                    annot_sel <-
                        col_annot_vars[[x]]
                    col_annot[x] <- selected_metadata_sample %>%
                    select(annot_sel) %>%
                    as.data.frame()
                    names(col_annot)[x] <- annot_sel
                }
                col_annot <- as.data.frame(col_annot)
                rownames(col_annot) <- rownames(data_path)
            }

            # Row annotation
            row_annot_vars <- metadata_info[grepl("color_Metab", names(metadata_info))]
            row_annot <- NULL
            if (length(row_annot_vars) > 0L) {
                for (y in seq_along(row_annot_vars)) {
                    annot_sel <-
                        row_annot_vars[[y]]
                    row_annot[y] <- metadata_feature %>% select(all_of(annot_sel))
                    row_annot <-
                        row_annot %>%
                        as.data.frame()
                    names(row_annot)[y] <- annot_sel
                }
                row_annot <- as.data.frame(row_annot)
                rownames(row_annot) <- rownames(metadata_feature)
            }

            # Check number of features:
            Features <- as.data.frame(t(data_path))
            if (enforce_featurenames) {
                show_rownames <- TRUE
                cellheight_Feature <- 9
            } else if (nrow(Features) > 100L) {
                show_rownames <- FALSE
                cellheight_Feature <- 1
            } else {
                show_rownames <- TRUE
                cellheight_Feature <- 9
            }

            # Check number of samples
            if (enforce_samplenames) {
                show_colnames <- TRUE
                cellwidth_Sample <- 9
            } else if (nrow(data_path) > 50L) {
                show_colnames <- FALSE
                cellwidth_Sample <- 1
            } else {
                show_colnames <- TRUE
                cellwidth_Sample <- 9
            }

            # Make the plot
            if (nrow(t(data_path)) >= 2L) {
                # set.seed(1234)

                heatmap <- pheatmap(t(data_path),
                    show_rownames = as.logical(show_rownames),
                    show_colnames = as.logical(show_colnames),
                    clustering_method = "complete",
                    scale = scale,
                    clustering_distance_rows = "correlation",
                    annotation_col = col_annot,
                    annotation_row = row_annot,
                    legend = TRUE,
                    cellwidth = cellwidth_Sample,
                    cellheight = cellheight_Feature,
                    fontsize_row = 10,
                    fontsize_col = 10,
                    fontsize = 9,
                    main = paste(plot_name, " Samples: ", i, sep = " "),
                    silent = TRUE
                )

                # ----- Save
                cleaned_i <-
                    gsub(
                        "[[:space:], /\\\\]",
                        "-",
                        i
                    )  # removes empty spaces and replaces /, \ with -
                PlotList[[cleaned_i]] <- heatmap

                # Width and height according to Sample and metabolite number
                Plot_Sized <-
                    plot_grob_heatmap(
                        input_plot = heatmap,
                        metadata_info = metadata_info,
                        metadata_sample = metadata_sample,
                        metadata_feature = metadata_feature,
                        plot_name = cleaned_i
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
                    add(theme(panel.background = element_rect(fill = "transparent")))

                PlotList_adaptedGrid[[cleaned_i]] <- Plot_Sized

                # ----- Save
                save_res(
                    inputlist_df = NULL,
                    inputlist_plot = PlotList_adaptedGrid,
                    save_table = NULL,
                    save_plot = save_plot,
                    path = folder,
                    file_name = paste("Heatmap_", plot_name, sep = ""),
                    core = FALSE,
                    print_plot = print_plot,
                    plot_height = plot_height,
                    plot_width = plot_width,
                    plot_unit = "cm"
                )
            } else {
                message <-
                    paste0(
                        i,
                        " includes <= 2L objects and is hence not plotted."
                    )
                log_trace(paste("Message ", message, sep = ""))
                message(message)
            }
        }

        # Return if assigned:
        return(invisible(list("Plot" = PlotList, "Plot_Sized" = PlotList_adaptedGrid)))

    } else if (
        "individual_Metab" %in% names(metadata_info) &
        "individual_Sample" %in% names(metadata_info)
    ) {
        # Ensure that groups that are assigned NAs do not cause problems:
        metadata_feature[[metadata_info[["individual_Metab"]]]] <- ifelse(is.na(metadata_feature[[metadata_info[["individual_Metab"]]]]), "NA", metadata_feature[[metadata_info[["individual_Metab"]]]])

        unique_paths <- unique(metadata_feature[[metadata_info[["individual_Metab"]]]])

        for (i in unique_paths) {  # Check pathways with 1 metabolite
            selected_path <-
                metadata_feature %>%
                filter(get(metadata_info[["individual_Metab"]]) == i)
            selected_path_metabs <- colnames(data)[colnames(data) %in% row.names(selected_path)]
            if (length(selected_path_metabs) == 1L) {
                message <-
                    paste0(
                        "The metadata group ",
                        i,
                        " includes only 1 metabolite. Heatmap cannot be made for 1 metabolite, thus it will be ignored."
                    )
                log_trace(paste("Warning ", message, sep = ""))
                warning(message)
                # Remove the pathway
                unique_paths <- unique_paths[!unique_paths %in% i]
            }
        }

        # Ensure that groups that are assigned NAs do not cause problems:
        metadata_sample[[metadata_info[["individual_Sample"]]]] <- ifelse(is.na(metadata_sample[[metadata_info[["individual_Sample"]]]]), "NA", metadata_sample[[metadata_info[["individual_Sample"]]]])

        unique_paths_Sample <- unique(metadata_sample[[metadata_info[["individual_Sample"]]]])

        for (i in unique_paths_Sample) {  # Check pathways with 1 metabolite
            selected_path <-
                metadata_sample %>%
                filter(get(metadata_info[["individual_Sample"]]) == i)
            selected_path_metabs <- colnames(data)[colnames(data) %in% row.names(selected_path)]
            if (length(selected_path_metabs) == 1L) {
                message <-
                    paste0(
                        "The metadata group ",
                        i,
                        " includes only 1 metabolite. Heatmap cannot be made for 1 metabolite, thus it will be ignored."
                    )
                log_trace(paste("Warning ", message, sep = ""))
                warning(message)
                # Remove the pathway
                unique_paths_Sample <- unique_paths_Sample[!unique_paths_Sample %in% i]
            }
        }

        IndividualPlots_Metab <- unique_paths
        IndividualPlots_Sample <- unique_paths_Sample

        PlotList <- list()  # Empty list to store all the plots
        PlotList_adaptedGrid <- list()  # Empty list to store all the plots

        for (i in IndividualPlots_Metab) {
            selected_path <-
                metadata_feature %>%
                filter(get(metadata_info[["individual_Metab"]]) == i)
            selected_path_metabs <-
                colnames(data)[colnames(data) %in% row.names(selected_path)]
            data_path_metab <- data %>% select(all_of(selected_path_metabs))

            # Row annotation
            row_annot_vars <- metadata_info[grepl("color_Metab", names(metadata_info))]
            row_annot <- NULL
            if (length(row_annot_vars) > 0L) {
                for (y in seq_along(row_annot_vars)) {
                    annot_sel <-
                        row_annot_vars[[y]]
                    row_annot[y] <- metadata_feature %>% select(all_of(annot_sel))
                    row_annot <-
                        row_annot %>%
                        as.data.frame()
                    names(row_annot)[y] <- annot_sel
                }
                row_annot <- as.data.frame(row_annot)
                rownames(row_annot) <- rownames(metadata_feature)
            }

            # Col annotation:
            for (s in IndividualPlots_Sample) {
                # Select the data:
                selected_path <-
                    metadata_sample %>%
                    filter(get(metadata_info[["individual_Sample"]]) == s) %>%
                    rownames_to_column("UniqueID")
                selected_path <-
                    as.data.frame(selected_path[, 1]) %>%
                    rename("UniqueID" = 1)
                data_path <-
                    merge(
                        selected_path,
                        data_path_metab %>% rownames_to_column("UniqueID"),
                        by = "UniqueID",
                        all.x = TRUE
                    ) %>%
                    column_to_rownames("UniqueID")

                # Column annotation
                selected_metadata_sample <-
                    merge(
                        selected_path,
                        metadata_sample %>% rownames_to_column("UniqueID"),
                        by = "UniqueID",
                        all.x = TRUE
                    )

                col_annot_vars <- metadata_info[grepl("color_Sample", names(metadata_info))]
                col_annot <- NULL
                if (length(col_annot_vars) > 0L) {
                    for (x in seq_along(col_annot_vars)) {
                        annot_sel <-
                            col_annot_vars[[x]]
                        col_annot[x] <- selected_metadata_sample %>%
                        select(annot_sel) %>%
                        as.data.frame()
                        names(col_annot)[x] <- annot_sel
                    }
                    col_annot <- as.data.frame(col_annot)
                    rownames(col_annot) <- rownames(data_path)
                }

                # Check number of features:
                Features <- as.data.frame(t(data_path))
                if (enforce_featurenames) {
                    show_rownames <- TRUE
                    cellheight_Feature <- 9
                } else if (nrow(Features) > 100L) {
                    show_rownames <- FALSE
                    cellheight_Feature <- 1
                } else {
                    show_rownames <- TRUE
                    cellheight_Feature <- 9
                }

                # Check number of samples
                if (enforce_samplenames) {
                    show_colnames <- TRUE
                    cellwidth_Sample <- 9
                } else if (nrow(data_path) > 50L) {
                    show_colnames <- FALSE
                    cellwidth_Sample <- 1
                } else {
                    show_colnames <- TRUE
                    cellwidth_Sample <- 9
                }

                # Make the plot
                if (nrow(t(data_path)) >= 2L) {
                    # set.seed(1234)

                    heatmap <- pheatmap(
                        t(data_path),
                        show_rownames = as.logical(show_rownames),
                        show_colnames = as.logical(show_colnames),
                        clustering_method = "complete",
                        scale = scale,
                        clustering_distance_rows = "correlation",
                        annotation_col = col_annot,
                        annotation_row = row_annot,
                        legend = TRUE,
                        cellwidth = cellwidth_Sample,
                        cellheight = cellheight_Feature,
                        fontsize_row = 10,
                        fontsize_col = 10,
                        fontsize = 9,
                        main = paste(plot_name, " Metabolites: ", i, " Sample:", s, sep = ""),
                        silent = TRUE
                    )

                    # # Store the plot in the 'plots' list
                    cleaned_i <-
                        gsub(
                            "[[:space:], /\\\\]",
                            "-",
                            i
                        )  # removes empty spaces and replaces /, \ with -
                    cleaned_s <-
                        gsub(
                            "[[:space:], /\\\\]",
                            "-",
                            s
                        )  # removes empty spaces and replaces /, \ with -
                    PlotList[[paste(cleaned_i, cleaned_s, sep = "_")]] <- heatmap

                    # -------- Plot width and heights
                    # Width and height according to Sample and metabolite number
                    plot_name <- paste(cleaned_i, cleaned_s, sep = "_")
                    Plot_Sized <-
                        plot_grob_heatmap(
                            input_plot = heatmap,
                            metadata_info = metadata_info,
                            metadata_sample = metadata_sample,
                            metadata_feature = metadata_feature,
                            plot_name = plot_name
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
                        add(theme(panel.background = element_rect(fill = "transparent")))

                    PlotList_adaptedGrid[[paste(cleaned_i, cleaned_s, sep = "_")]] <- Plot_Sized

                    # ----- Save
                    save_res(
                        inputlist_df = NULL,
                        inputlist_plot = PlotList_adaptedGrid,
                        save_table = NULL,
                        save_plot = save_plot,
                        path = folder,
                        file_name = paste("Heatmap_", plot_name, sep = ""),
                        core = FALSE,
                        print_plot = print_plot,
                        plot_height = plot_height,
                        plot_width = plot_width,
                        plot_unit = "cm"
                    )
                } else {
                    message(
                        i,
                        " includes <= 2L objects and is hence not plotted."
                    )
                }
            }
        }

        return(invisible(list("Plot" = PlotList, "Plot_Sized" = PlotList_adaptedGrid)))

    } else if (c("individual_Metab", "individual_Sample") %>% is_in(names(metadata_info)) %>% any %>% not) {

        PlotList <- list()  # Empty list to store all the plots
        PlotList_adaptedGrid <- list()  # Empty list to store all the plots

        # Column annotation
        col_annot_vars <- metadata_info[grepl("color_Sample", names(metadata_info))]
        col_annot <- NULL
        if (length(col_annot_vars) > 0L) {
            for (i in seq_along(col_annot_vars)) {
                annot_sel <-
                    col_annot_vars[[i]]
                col_annot[i] <- metadata_sample %>%
                select(annot_sel) %>%
                as.data.frame()
                names(col_annot)[i] <- annot_sel
            }
            col_annot <- as.data.frame(col_annot)
            rownames(col_annot) <- rownames(data)
        }

        # Row annotation
        row_annot_vars <- metadata_info[grepl("color_Metab", names(metadata_info))]
        row_annot <- NULL
        if (length(row_annot_vars) > 0L) {
            for (i in seq_along(row_annot_vars)) {
                annot_sel <-
                    row_annot_vars[[i]]
                row_annot[i] <- metadata_feature %>% select(all_of(annot_sel))
                row_annot <-
                    row_annot %>%
                    as.data.frame()
                names(row_annot)[i] <- annot_sel
            }
            row_annot <- as.data.frame(row_annot)
            rownames(row_annot) <- rownames(metadata_feature)
        }

        # Check number of features:
        Features <- as.data.frame(t(data))
        if (enforce_featurenames) {
            show_rownames <- TRUE
            cellheight_Feature <- 9
        } else if (nrow(Features) > 100L) {
            show_rownames <- FALSE
            cellheight_Feature <- 1
        } else {
            show_rownames <- TRUE
            cellheight_Feature <- 9
        }

        # Check number of samples
        if (enforce_samplenames) {
            show_colnames <- TRUE
            cellwidth_Sample <- 9
        } else if (nrow(data) > 50L) {
            show_colnames <- FALSE
            cellwidth_Sample <- 1
        } else {
            show_colnames <- TRUE
            cellwidth_Sample <- 9
        }

        # Make the plot:
        if (nrow(t(data)) >= 2L) {
            # set.seed(1234)

            heatmap <- pheatmap(t(data),
                show_rownames = as.logical(show_rownames),
                show_colnames = as.logical(show_colnames),
                clustering_method = "complete",
                scale = scale,
                clustering_distance_rows = "correlation",
                annotation_col = col_annot,
                annotation_row = row_annot,
                legend = TRUE,
                cellwidth = cellwidth_Sample,
                cellheight = cellheight_Feature,
                fontsize_row = 10,
                fontsize_col = 10,
                fontsize = 9,
                main = plot_name,
                silent = TRUE
            )

            # # Store the plot in the 'plots' list
            PlotList[[plot_name]] <- heatmap

            # -------- Plot width and heights
            # Width and height according to Sample and metabolite number
            Plot_Sized <-
                plot_grob_heatmap(
                    input_plot = heatmap,
                    metadata_info = metadata_info,
                    metadata_sample = metadata_sample,
                    metadata_feature = metadata_feature,
                    plot_name = plot_name
                )
            plot_height <-
                convertUnit(
                    Plot_Sized$height,
                    "cm",
                    valueOnly = TRUE
                )
            plot_width <- convertUnit(Plot_Sized$width, "cm", valueOnly = TRUE)
            Plot_Sized %<>%
                {
                    ggplot() +
                        annotation_custom(.)
                } %>%
                add(theme(panel.background = element_rect(fill = "transparent")))

            PlotList_adaptedGrid[[paste("Heatmap_", plot_name, sep = "")]] <- Plot_Sized

            # ----- Save
            save_res(
                inputlist_df = NULL,
                inputlist_plot = PlotList_adaptedGrid,
                save_table = NULL,
                save_plot = save_plot,
                path = folder,
                file_name = paste("Heatmap_", plot_name, sep = ""),
                core = FALSE,
                print_plot = print_plot,
                plot_height = plot_height,
                plot_width = plot_width,
                plot_unit = "cm"
            )
        } else {
            message <-
                paste0(
                    plot_name,
                    " includes <= 2L objects and is hence not plotted."
                )
            log_trace(paste("Message ", message, sep = ""))
            message(message)
        }

        return(invisible(list("Plot" = PlotList, "Plot_Sized" = PlotList_adaptedGrid)))

    }

}
