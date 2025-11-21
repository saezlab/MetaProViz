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
# Function to create complex upset plots to visualise PK coverage
#

#' Generate upset plot for set intersections
#'
#' Visualizes Prior Knowledge (PK) coverage by displaying intersections
#' among a set of metabolite ID columns. If a class column is provided,
#' data is grouped by that column and a color palette ("viridis" or
#' "polychrome") is used to represent class levels, with a corresponding
#' legend (hidden if too many unique classes). If no class column is
#' provided, a basic upset plot is generated.
#'
#' @param df Data frame containing the data to be plotted.
#' @param class_col Character (optional). Name of the column in df that
#'      represents the class of each observation. This column is coerced
#'      to a factor if provided. Default: NULL.
#' @param intersect_cols Character vector (optional). Names of the
#'      columns in df to be used for generating intersections. Default:
#'      c("LIMID", "HMDB", "CHEBI", "None").
#' @param plot_name Character (optional). String added to output files
#'      of the upset plot. Default: "".
#' @param palette_type Character (optional). Color palette to use for
#'      fill aesthetic when class_col is provided. Options are "viridis"
#'      and "polychrome". Default: c("viridis", "polychrome").
#' @param max_legend_terms Numeric (optional). Maximum number of unique
#'      terms in class_col for which legend should be displayed. If
#'      number of levels exceeds this value, legend will be hidden.
#'      Ignored if class_col is NULL. Default: 20.
#' @param save_plot Character (optional). File type of output plots:
#'      "svg", "png", "pdf", or NULL. Default: "svg".
#' @param print_plot Logical (optional). Whether volcano plot is saved
#'      as overview of results. Default: TRUE.
#' @param path Character (optional). Path to folder where results
#'      should be saved. Default: NULL.
#'
#' @return ggplot object representing the generated upset plot.
#'
#' @importFrom ggplot2 scale_fill_manual scale_fill_viridis_d element_text theme
#' @importFrom ggplot2 aes_string theme_minimal margin labs
#' @importFrom Polychrome palette36.colors
#' @importFrom ComplexUpset intersection_size upset upset_set_size
#' @importFrom logger log_info log_trace
#' @importFrom stats setNames
#' @importFrom grid convertUnit
#' @noRd
viz_upset <- function(
    df,
    class_col = NULL,
    intersect_cols = c("LIMID", "HMDB", "CHEBI", "None"),
    plot_name = "Metabolite IDs",
    palette_type = c("viridis", "polychrome"),
    max_legend_terms = 20,
    save_plot = "svg",
    print_plot = TRUE,
    path = NULL) {
    # NSE vs. R CMD check workaround
    upset_plot <- NULL

    # ##########################################################################
    # # ------------ Create log file ----------- ##
    metaproviz_init()

    log_info("viz_upset: Upset plot visualization")
    # # ------------ Check Input files ----------- ##
    #
    palette_type <- match.arg(
        palette_type,
        c("viridis", "polychrome")
    )


    # # ------------ Create Results output folder ----------- ##
    if (!is.null(save_plot)) {
        folder <- save_path(
            folder_name = "UpsetPlots",
            path = path
        )
        log_info("viz_upset results saved at ", folder)
    }


    # ##########################################################################
    # # ----------- Check input data frame ----------- ##
    # If either a list of DFs (different PK resources, same IDs) or a single DF but multiple ID columns or a combination of both!
    # currently there is also the compare_pk function, but I think we could maybe use it in here too!


    # # ----------- Set the plot parameters: ------------ ##
    # # --- Prepare colour palette

    # If a class column is provided, process it for fill aesthetics
    if (!is.null(class_col)) {
        df[[class_col]] <- as.factor(df[[class_col]])
        fill_scale <- scale_fill_viridis_d(option = "viridis")
        if (palette_type == "viridis") {
            fill_scale <- scale_fill_viridis_d(option = "viridis")
        } else if (palette_type == "polychrome") {
            class_levels <- levels(df[[class_col]])
            my_palette <- palette36.colors(n = 36)
            if (length(my_palette) < length(class_levels)) {
                fill_scale <- scale_fill_viridis_d(option = "viridis")
                message <- paste0("Not enough colors in the Polychrome palette for the number of classes. Hence viridis palette was used instead.")
                warning("Not enough colors in the Polychrome palette for the number of classes!")

                log_trace(paste("Warning ", message, sep = ""))
            }
            my_palette_named <- setNames(my_palette[seq_along(class_levels)], class_levels)
            fill_scale <- scale_fill_manual(values = my_palette_named)
        }

        # Build the base annotation with a mapping for fill based on the class column.
        base_annotation <- list(
            "Intersection size" = intersection_size(
                mapping = aes_string(fill = class_col),
                counts = TRUE
            ) + fill_scale
        )
    } else {
        # No class column provided: use default annotation without fill mapping.
        base_annotation <- list(
            "Intersection size" = intersection_size(counts = TRUE)
        )
    }

    # # ----------- Make the  plot based on the choosen parameters ------------ ##
    # Create the upset plot
    p <- upset(
        data = df,
        intersect = intersect_cols,
        name = plot_name,
        base_annotations = base_annotation,
        set_sizes = upset_set_size()
    )


    # # ----------- Save and return -------------#
    save_res(
        inputlist_df = NULL,
        inputlist_plot = list(upset_plot = p),
        save_table = NULL,
        save_plot = save_plot,
        path = folder,
        file_name = "UpsetPlot",
        core = FALSE,
        print_plot = print_plot
    )

    return(p)
}
