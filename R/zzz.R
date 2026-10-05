#!/usr/bin/env Rscript

#
#  This file is part of the `MetaProViz` R package
#
#  Copyright 2023-2026
#  Saez Lab, Heidelberg University, University Hospital Heidelberg,
#  European Bioinformatics Institute (EMBL-EBI), University of Cologne
#
#  Authors: see the file `README.md`
#
#  Distributed under the BSD 3-Clause License.
#  See accompanying file `LICENSE.md` or copy at
#      https://opensource.org/license/bsd-3-clause
#
#  Website: https://saezlab.github.io/MetaProViz
#  Git repo: https://github.com/saezlab/MetaProViz
#


#' Package initialization hook
#'
#' @importFrom magrittr %>%
#' @importFrom rlang !!
#' @noRd
.onLoad <- function(
    libname,
    pkgname
) {
    opr <- "OmnipathR"
    ddb <- "disable_doctest_bypass"

    if (exists(ddb, where = asNamespace(opr), mode = "function")) {
        ((!!opr) %:::% (!!ddb))()
    }
}
