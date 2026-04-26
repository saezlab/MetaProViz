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
        path %<>% file.path("MetaProViz_Results")
        if (!dir.exists(path)) {
            dir.create(path)
        }
    } else {
        if (!dir.exists(path)) {
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
#' @importFrom grid convertUnit unit grid.draw
#' @importFrom logger log_warn
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
    if (!is.null(save_table)) {
        # Excel File: One file with multiple sheets:
        if (save_table == "xlsx") {
            # Make file_name
            if (!core | is.null(core)) {
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
                if (!core | is.null(core)) {
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
    if (!is.null(save_plot)) {
        for (Plot in names(inputlist_plot)) {
            # Make file_name
            if (!core | is.null(core)) {
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

            # Plots that don't go through ggsave cleanly:
            #   - ComplexUpset has theme conflicts with ggsave.
            #   - Plots wrapped with `with_canvas_size` (heatmap, PCA,
            #     volcano, superplot, …) carry their own canvas
            #     dimensions and must be drawn on a device that matches.
            #     Otherwise grid's viewport math goes non-finite and
            #     dies with "non-finite location and/or size for
            #     viewport".
            plot_obj <- inputlist_plot[[Plot]]
            is_upset_plot <- inherits(plot_obj, "upset_plot") ||
                (
                    "package:ComplexUpset" %in% search() &&
                    inherits(plot_obj, "ggplot")
                )
            is_canvas_sized <- inherits(plot_obj, "with_canvas_size")

            if (is_upset_plot || grepl("upset", Plot, ignore.case = TRUE) ||
                is_canvas_sized) {

                file_name_full <- paste0(file_name_Save, ".", save_plot, sep = "")

                # If the plot carries a canvas size, prefer that over
                # the caller-passed plot_width/plot_height — using the
                # wrong dimensions is exactly what trips grid's viewport
                # check.
                if (is_canvas_sized && !is.null(plot_obj$width) &&
                        !is.null(plot_obj$height)) {
                    width_in <- convertUnit(plot_obj$width, "inches", valueOnly = TRUE)
                    height_in <- convertUnit(plot_obj$height, "inches", valueOnly = TRUE)
                } else {
                    width_in <- convertUnit(unit(plot_width, plot_unit), "inches", valueOnly = TRUE)
                    height_in <- convertUnit(unit(plot_height, plot_unit), "inches", valueOnly = TRUE)
                }

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

                tryCatch(
                    {
                        if (is_upset_plot) {
                            grid.draw(patchwork::patchworkGrob(plot_obj))
                        } else {
                            # with_canvas_size gtables and other grobs
                            # render directly with grid.draw.
                            grid.draw(plot_obj)
                        }
                    },
                    error = function(e) {
                        log_warn(
                            "Could not render plot %s to %s: %s",
                            Plot, save_plot, conditionMessage(e)
                        )
                        # Last resort: blank page with error label so
                        # the file exists and the loop continues.
                        plot.new()
                        text(
                            0.5, 0.5,
                            paste("Error rendering plot:\n", conditionMessage(e)),
                            cex = 0.8, col = "red"
                        )
                    }
                )
                dev.off()
            } else {
                # Use ggsave for regular ggplot2 plots. Wrap in
                # tryCatch: heatmaps wrapped via ggplot()+annotation_custom
                # over a pheatmap-derived gtable can fail viewport math
                # in batch contexts ("non-finite location and/or size for
                # viewport"). Emit a clear warning rather than halting
                # the entire enclosing analysis.
                tryCatch(
                    ggsave(
                        filename = paste0(file_name_Save, ".", save_plot, sep = ""),
                        plot = plot_obj,
                        width = plot_width,
                        height = plot_height,
                        units = plot_unit
                    ),
                    error = function(e) {
                        log_warn(
                            "ggsave failed for plot %s: %s",
                            Plot, conditionMessage(e)
                        )
                        warning(
                            "MetaProViz could not save plot ", Plot, " (",
                            save_plot, "): ", conditionMessage(e),
                            "\nThe analysis result is still returned in memory; ",
                            "you can save it manually with `ggsave()`.",
                            call. = FALSE
                        )
                    }
                )
            }

            if (print_plot) {
                # Special handling for upset plots-they may have theme issues
                if (is_upset_plot || grepl("upset", Plot, ignore.case = TRUE)) {
                    tryCatch(
                        {
                            grid.draw(patchwork::patchworkGrob(plot_obj))
                        },
                        error = function(e) {
                            msg <- paste0(
                                "Note: Could not render upset plot to screen due to theme compatibility ",
                                "issues. Plot was saved to file successfully."
                            )
                            log_info(msg)
                            message(msg)
                        }
                    )
                } else {
                    plot(inputlist_plot[[Plot]])
                }
            }
        }
    }
}
