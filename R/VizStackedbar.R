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



##########################################################################################
### ### ### Function to create stacked bar plots to visualise PK coverage ### ### ###
##########################################################################################

#' Generate Stacked Bar Plot, e.g. for PK Coverage
#'
#' This function creates a stacked bar plot, for example to visualize Prior Knowledge (PK) coverage.
#' The plot groups data by a specified column and fills the bars according to another column,
#' allowing you to inspect the distribution of match statuses (or any categorical variable).
#'
#' @param data A data frame containing the data to be plotted.
#' @param group_col A string specifying the name of the column to group the data by.
#' @param fill_col A string specifying the name of the column to use for fill aesthetics.
#' @param fill_values A vector of color values to be used for the fill aesthetic.
#' @param fill_labels A vector of labels corresponding to the fill levels for the legend.
#' @param plot_name A string specifying the title of the plot.
#' @param x_label A string for the x-axis label. Defaults to "Frequency".
#' @param y_label A string for the y-axis label. If \code{NULL}, the value of \code{group_col} is used.
#' @param legend_position A numeric vector of length 2 specifying the (x, y) position of the legend.
#'                        Defaults to \code{c(0.95, 0.05)}.
#'
#' @return A \code{ggplot} object representing the stacked bar plot.
#'
#' @importFrom dplyr arrange group_by mutate n pull
#' @importFrom dplyr summarise
#' @importFrom ggplot2 aes_string element_text geom_bar ggplot labs
#' @importFrom ggplot2 scale_fill_manual theme theme_minimal
#' @importFrom rlang sym
#' @noRd
viz_stackedbar <- function(data,
                          group_col,
                          fill_col,
                          fill_values,
                          fill_labels,
                          plot_name,
                          x_label = "Frequency",
                          y_label = NULL,
                          legend_position = c(0.95, 0.05)) {
  # Convert column names to symbols for tidy evaluation
  group_sym <- sym(group_col)
  fill_sym  <- sym(fill_col)

  # Determine order of groups by overall frequency (ascending)
  group_order <- data %>%
    group_by(!!group_sym) %>%
    summarise(total = n(), .groups = 'drop') %>%
    arrange(total) %>%
    pull(!!group_sym)

  # Summarize data by group and fill status, then reorder the group factor
  summary_data <- data %>%
    group_by(!!group_sym, !!fill_sym) %>%
    summarise(count = n(), .groups = 'drop') %>%
    mutate(!!group_sym := factor(!!group_sym, levels = group_order))

  # If y_label is not provided, use the grouping column name
  if (is.null(y_label)) {
    y_label <- group_col
  }

  # Create the plot
  plot <- ggplot(summary_data, aes_string(y = group_col,
                                                    x = "count",
                                                    fill = paste0("as.factor(", fill_col, ")"))) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = fill_values,
                               labels = fill_labels,
                               name = "Match Status") +
    labs(title = plot_name,
                  x = x_label,
                  y = y_label) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 1),
                   legend.position = legend_position,
                   legend.justification = c("right", "bottom"),
                   plot.title = element_text(hjust = 0.4))

  ## ----------- Save and return -------------#
  # Only needed if exported!


  return(plot)
}
