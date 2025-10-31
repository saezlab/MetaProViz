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
# General helper function: Save folder path
#

#' save_path is the helper function to create the folder structure and path
#'
#' Create folder and path
#'
#' @param folder_name Name of the folder, which can not contain any special characters.
#'     Created within  the individual MetaProViz functions and can not be
#'     changed by the user.
#' @param path Passed to main function by the user
#'
#' @noRd
save_path <- function(
    folder_name,
    path
) {
    # Check if folder_name includes special characters that are not allowed
    cleaned_folder_name <- gsub("[^a-zA-Z0-9 ]", "", folder_name)
    if (folder_name != cleaned_folder_name) {
        message("Special characters were removed from `folder_name`.")
    }

    # Check if path exist
    if (is.null(path)) {
        path <- getwd()
        path <- file.path(path, "MetaProViz_Results")
        if (!dir.exists(path)) {
            dir.create(path)
        }
    } else {
        if (dir.exists(path) == FALSE) {
            path <- getwd()
            message(
                "Provided `path` does not exist and hence results are saved here: ",
                path,
                sep = ""
            )
        }
    }

    # Create the folder name
    Results_folder <- file.path(path, cleaned_folder_name)
    if (!dir.exists(Results_folder)) {
        dir.create(Results_folder)
    }

    # Return the folder path:
    return(invisible(Results_folder))
}


#' Make sure the results directory exists
#'
#' @importFrom magrittr %>%
#' @noRd
results_dir <- function(
    path = "MetaProViz_Results"
) {
    # toDO: options?
    path %>%
    {
            `if`(
                !dir.exists(.),
                {
                    dir.create(.)
                    .
                },
                .
            )
        }
}


#
# General helper function: Save tables and plot
#

#' save_res is the helper function to save the plots and tables
#'
#' @param inputlist_df \emph{Optional: } Generated within the MetaProViz function. Contains
#'     named DFs. If not avalailable can be set to NULL.\strong{Default = NULL}
#' @param inputlist_plot \emph{Optional: } Generated within the MetaProViz function. Contains
#'     named Plots. If not avalailable can be set to NULL.\strong{Default =
#'     NULL}
#' @param save_table \emph{Optional: } Passed to main function by the user. If not
#'     avalailable can be set to NULL.\strong{Default = NULL}
#' @param save_plot \emph{Optional: } Passed to main function by the user. If not
#'     avalailable can be set to NULL. \strong{Default = NULL}
#' @param path Passed to main function by the user.
#' @param file_name Passed to main function by the user.
#' @param core \emph{Optional: } Passed to main function by the user. If not
#'     avalailable can be set to NULL.\strong{Default = FALSE}
#' @param print_plot \emph{Optional: } Passed to main function by the user. If not
#'     avalailable can be set to NULL.\strong{Default = TRUE}
#' @param plot_height \emph{Optional: } Parameter for ggsave.\strong{Default = NULL}
#' @param plot_width \emph{Optional: } Parameter for ggsave. \strong{Default = NULL}
#' @param plot_unit \emph{Optional: } Parameter for ggsave. \strong{Default = NULL}
#'
#' @importFrom ggplot2 ggsave
#' @importFrom readr write_csv write_delim
#' @importFrom writexl write_xlsx
#' @importFrom tidyselect where
#' @importFrom grDevices svg png pdf
#' @importFrom graphics plot.new text
#' @noRd
save_res <- function(
    inputlist_df = NULL,
    inputlist_plot = NULL,
    save_table = NULL,
    save_plot = NULL,
    path,
    file_name,
    core = FALSE,
    print_plot = TRUE,
    plot_height = NULL,
    plot_width = NULL,
    plot_unit = NULL
) {
    # ############### Save Tables:
    if (is.null(save_table) == FALSE) {
        # Excel File: One file with multiple sheets:
        if (save_table == "xlsx") {
            # Make file_name
            if (core == FALSE | is.null(core) == TRUE) {
                file_name <-
                    paste0(
                        path,
                        "/",
                        file_name,
                        "_",
                        Sys.Date(),
                        sep = ""
                    )
            } else {
                file_name <-
                    paste0(
                        path,
                        "/core_",
                        file_name,
                        "_",
                        Sys.Date(),
                        sep = ""
                    )
            }
            # Save Excel
            write_xlsx(
                inputlist_df,
                paste0(file_name, ".xlsx", sep = ""),
                col_names = TRUE
            )
        } else {
            for (DF in names(inputlist_df)) {
                # Make file_name
                if (core == FALSE | is.null(core) == TRUE) {
                    file_name_Save <-
                        paste0(
                            path,
                            "/",
                            file_name,
                            "_",
                            DF,
                            "_",
                            Sys.Date(),
                            sep = ""
                        )
                } else {
                    file_name_Save <-
                        paste0(
                            path,
                            "/core_",
                            file_name,
                            "_",
                            DF,
                            "_",
                            Sys.Date(),
                            sep = ""
                        )
                }

                # unlist DF columns if needed
                inputlist_df[[DF]] <- inputlist_df[[DF]] %>%
                mutate(
                        across(
                            where(is.list),
                            ~ map_chr(.x, ~ paste(sort(unique(.x)), collapse = "; "))
                        )
                    )
                # Save table
                if (save_table == "csv") {
                    inputlist_df[[DF]] %>%
                    write_csv(paste0(file_name_Save, ".csv", sep = ""))
                } else if (save_table == "txt") {
                    inputlist_df[[DF]] %>%
                    write_delim(paste0(file_name_Save, ".csv", sep = ""))
                }
            }
        }
    }

    # ############### Save Plots:
    if (is.null(save_plot) == FALSE) {
        for (Plot in names(inputlist_plot)) {
            # Make file_name
            if (core == FALSE | is.null(core) == TRUE) {
                file_name_Save <-
                    paste0(
                        path,
                        "/",
                        file_name,
                        "_",
                        Plot,
                        "_",
                        Sys.Date(),
                        sep = ""
                    )
            } else {
                file_name_Save <-
                    paste0(
                        path,
                        "/core_",
                        file_name,
                        "_",
                        Plot,
                        "_",
                        Sys.Date(),
                        sep = ""
                    )
            }

            # Save
            if (is.null(plot_height)) {
                plot_height <- 12
            }
            if (is.null(plot_width)) {
                plot_width <- 16
            }
            if (is.null(plot_unit)) {
                plot_unit <- "cm"
            }

            # Check if this is a ComplexUpset plot (has S7 class from ComplexUpset package)
            # ComplexUpset plots have theme conflicts with ggsave, so use device-based saving
            plot_obj <- inputlist_plot[[Plot]]
            is_upset_plot <- inherits(plot_obj, "upset_plot") ||
                ("package:ComplexUpset" %in% search() && inherits(plot_obj, "ggplot"))

            if (is_upset_plot || grepl("upset", Plot, ignore.case = TRUE)) {
                # Use device-based saving for upset plots
                file_name_full <-
                    paste0(
                        file_name_Save,
                        ".",
                        save_plot,
                        sep = ""
                    )

                # Convert to inches for device functions
                width_in <- grid::convertUnit(grid::unit(plot_width, plot_unit), "inches", valueOnly = TRUE)
                height_in <- grid::convertUnit(grid::unit(plot_height, plot_unit), "inches", valueOnly = TRUE)

                # Open device based on save_plot type
                if (save_plot == "svg") {
                    svg(file_name_full, width = width_in, height = height_in)
                } else if (save_plot == "pdf") {
                    pdf(file_name_full, width = width_in, height = height_in)
                } else if (save_plot == "png") {
                    png(
                        file_name_full,
                        width = width_in,
                        height = height_in,
                        units = "in",
                        res = 300
                    )
                }

                # Plot and close device
                # ComplexUpset plots may have theme issues with some ggplot2 versions
                # Use gridExtra to render the plot which handles patchwork objects better
                tryCatch(
                    {
                        grid::grid.draw(patchwork::patchworkGrob(plot_obj))
                    },
                    error = function(e) {
                        # Fallback: try direct print which sometimes works
                        tryCatch(
                            {
                                log_info(plot_obj)
                            },
                            error = function(e2) {
                                # Last resort: save a blank plot with error message
                                plot.new()
                                text(
                                    0.5,
                                    0.5,
                                    paste("Error rendering ComplexUpset plot:\n", e2$message
                                ),
                                    cex = 0.8, col = "red"
                                )
                            }
                        )
                    }
                )
                dev.off()
            } else {
                # Use ggsave for regular ggplot2 plots
                ggsave(
                    filename = paste0(file_name_Save, ".", save_plot, sep = ""),
                    plot = plot_obj,
                    width = plot_width,
                    height = plot_height,
                    units = plot_unit
                )
            }

            if (print_plot == TRUE) {
                # Special handling for upset plots-they may have theme issues
                if (is_upset_plot || grepl("upset", Plot, ignore.case = TRUE)) {
                    tryCatch(
                        {
                            grid::grid.draw(patchwork::patchworkGrob(plot_obj))
                        },
                        error = function(e) {
                            message("Note: Could not render upset plot to screen due to theme compatibility issues. Plot was saved to file successfully.")
                        }
                    )
                } else {
                    plot(inputlist_plot[[Plot]])
                }
            }
        }
    }
}
