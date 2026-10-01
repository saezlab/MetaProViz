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

#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:\%>\%]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling \code{rhs(lhs)}.
#'
#' @examples
#' c(1, 4, 9) %>% sqrt()
#'
#' @importFrom magrittr %>%
#' @export
NULL
