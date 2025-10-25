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
### ### ### Function to create complex upset plots to visualise PK coverage ### ### ###
##########################################################################################

#' Generate Complex Upset Plot for PK Coverage (Optional Class Grouping)
#'
#' This function creates a complex upset plot.
#' It is mostly intended to visualize Prior Knowledge (PK) coverage by displaying the intersections among
#' a set of metabolite ID columns, but its use could be extended too. If a class column is provided,
#' the data is grouped by that column and a color palette ("viridis" or "polychrome") is used to represent
#' the class levels, along with a corresponding legend (which can be hidden if there are too many unique classes).
#' If no class column is provided (\code{class_col = NULL}), a basic upset plot is generated.
#'
#' @param df A data frame containing the data to be plotted.
#' @param class_col \emph{Optional: } An optional string specifying the name of the column in \code{df} that represents the class
#'                  of each observation. This column is coerced to a factor if provided. \strong{Default = NULL}
#' @param intersect_cols \emph{Optional: } A character vector specifying the names of the columns in \code{df} to be used for generating intersections.
#'                       \strong{Default = c("LIMID", "HMDB", "CHEBI", "None")}.
#' @param plot_name \emph{Optional: } String which is added to the output files of the Upsetplot\strong{Default = ""}
#' @param palette_type \emph{Optional: } A string specifying the color palette to use for the fill aesthetic when \code{class_col} is provided.
#'                     Options are \code{"viridis"} and \code{"polychrome"}. \strong{Default = c("viridis", "polychrome")}
#' @param max_legend_terms \emph{Optional: } Numeric value specifying the maximum number of unique terms in \code{class_col}
#'                         for which the legend should be displayed. If the number of levels exceeds this value,
#'                         the legend will be hidden. Ignored if \code{class_col} is \code{NULL}. \strong{Default = 20}.
#' @param save_plot \emph{Optional: } Select the file type of output plots. Options are svg, png, pdf or NULL. \strong{Default = svg}
#' @param print_plot \emph{Optional: } TRUE or FALSE, if TRUE Volcano plot is saved as an overview of the results. \strong{Default = TRUE}
#' @param path \emph{Optional:} Path to the folder the results should be saved at. \strong{Default= NULL}

#'
#' @return A \code{ggplot} object representing the generated upset plot.
#'
#'
#' @keywords upset plot, complex upset
#'
#' @importFrom ggplot2 scale_fill_manual scale_fill_viridis_d element_text theme
#' @importFrom ggplot2 aes_string theme_minimal margin labs
#' @importFrom Polychrome palette36.colors
#' @importFrom ComplexUpset intersection_size upset upset_set_size
#' @importFrom logger log_info log_trace
#' @importFrom stats setNames
#' @importFrom grid convertUnit
#'
#' @noRd
viz_upset <- function(df,
                     class_col = NULL,
                     intersect_cols = c("LIMID", "HMDB", "CHEBI", "None"),
                     plot_name = "Metabolite IDs",
                     palette_type = c("viridis", "polychrome"),
                     max_legend_terms = 20,
                     save_plot = "svg",
                     print_plot=TRUE,
                     path = NULL) {

  ###########################################################################
  ## ------------ Create log file ----------- ##
  metaproviz_init()

  log_info("viz_upset: Upset plot visualization")
  ## ------------ Check Input files ----------- ##
  #
  palette_type <- match.arg(palette_type,
                            c("viridis", "polychrome"))


  ## ------------ Create Results output folder ----------- ##
  if(is.null(save_plot)==FALSE){
    folder <- save_path(folder_name= "UpsetPlots",
                       path=path)
    log_info("viz_upset results saved at ", folder)
  }


  ###########################################################################
  ## ----------- Check input data frame ----------- ##
  # If either a list of DFs (different PK resources, same IDs) or a single DF but multiple ID columns or a combination of both!
  # currently there is also the compare_pk function, but I think we could maybe use it in here too!



  ## ----------- Set the plot parameters: ------------ ##
  ##--- Prepare colour palette

  # If a class column is provided, process it for fill aesthetics
  if (!is.null(class_col)) {
    df[[class_col]] <- as.factor(df[[class_col]])
    fill_scale <- scale_fill_viridis_d(option = "viridis")
    if(palette_type == "viridis"){
      fill_scale <- scale_fill_viridis_d(option = "viridis")
    } else if(palette_type == "polychrome"){
      class_levels <- levels(df[[class_col]])
      my_palette <- palette36.colors(n = 36)
      if(length(my_palette) < length(class_levels)) {
        fill_scale <- scale_fill_viridis_d(option = "viridis")
        message <- paste0("Not enough colors in the Polychrome palette for the number of classes. Hence viridis palette was used instead.")
        warning("Not enough colors in the Polychrome palette for the number of classes!")

        log_trace(paste("Warning ", message, sep=""))
      }
      my_palette_named <- setNames(my_palette[1:length(class_levels)], class_levels)
      fill_scale <- scale_fill_manual(values = my_palette_named)
    }
    # Build the base annotation with a mapping for fill based on the class column.
    base_annotation <- list(
      "Intersection size" = intersection_size(
        mapping = aes_string(fill = class_col),
        counts = TRUE
      ) + fill_scale +
        theme(
          legend.position = "right"
        )
    )
  } else {
    # No class column provided: use default annotation without fill mapping.
    base_annotation <- list(
      "Intersection size" = intersection_size(counts = TRUE)
    )
  }

  ## ----------- Make the  plot based on the choosen parameters ------------ ##
  # Create the upset plot
  p <- upset(
    data = df,
    intersect = intersect_cols,
    name = plot_name,
    base_annotations = base_annotation,
    set_sizes = (
      upset_set_size() +
        theme(axis.text.y = element_text(size = 10))
    )
  ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.margin = margin(1, 1, 1, 1, "cm")
    )

  # If a class column was provided, hide the legend if there are too many unique terms
  if (!is.null(class_col) && length(levels(df[[class_col]])) > max_legend_terms) {
    p <- p + theme(legend.position = "none")
  }


  ## ----------- Save and return -------------#
  suppressMessages(suppressWarnings(
    save_res(inputlist_df=NULL,
             inputlist_plot= list(upset_plot = upset_plot),
             save_table=NULL,
             save_plot=save_plot,
             path= folder,
             file_name= "UpsetPlot",
             core=FALSE,
             print_plot=print_plot)))

  return(p)
}
