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
