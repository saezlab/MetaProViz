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

#' Makes sure we have a string even if the argument was passed by NSE
#'
#' @importFrom magrittr %>%
#' @importFrom rlang enquo quo_get_expr quo_text is_symbol
#' @importFrom stringr str_remove
#' @noRd
.nse_ensure_str <- function(arg) {
    enquo(arg) %>%
    {
            `if`(
                is_symbol(quo_get_expr(.)),
                quo_text(.),
                quo_get_expr(.)
            )
        } %>%
        str_remove("`$") %>%
        str_remove("^`")
}


#' Workaround against R CMD check notes about using `:::`
#'
#' @importFrom rlang enquo !!
#' @noRd
`%:::%` <- function(pkg, fun) {
    pkg <- .nse_ensure_str(!!enquo(pkg))
    fun <- .nse_ensure_str(!!enquo(fun))

    get(fun, envir = asNamespace(pkg), inherits = FALSE)
}


#' Hotelling's T2 control chart for individual multivariate observations
#'
#' Base-R replacement for `qcc::mqcc(type = "T2.single")`; returns the same
#' fields `outlier_detection()` reads (statistics, limits, type,
#' confidence.level, violations$beyond.limits). 
#'
#' `stats.T2.single()` (cov() + mahalanobis()) computes the T2 statistic
#' itself: Hotelling, H. (1931), Annals of Mathematical Statistics. 2 (3),
#' 360-378, doi:https://doi.org/10.1214/aoms/1177732979.
#'
#' `limits.T2.single()` (qbeta() control limit): Tracy, Young & Young
#' (1992); see also Mason & Young (2002), Montgomery (2013).
#'
#' @importFrom stats mahalanobis qbeta cov
#' @noRd
.hotelling_t2_single <- function(data, confidence.level) {
    m <- nrow(data)
    p <- ncol(data)

    center <- colMeans(data)
    cov <- cov(data)

    statistics <- mahalanobis(data, center, cov)
    names(statistics) <- rownames(data)

    ucl <- (m - 1)^2 / m * qbeta(confidence.level, p / 2, (m - p - 1) / 2)
    limits <- matrix(c(0, ucl), ncol = 2, dimnames = list("", c("LCL", "UCL")))

    list(
        statistics = statistics,
        limits = limits,
        type = "T2.single",
        confidence.level = confidence.level,
        violations = list(beyond.limits = which(statistics > ucl))
    )
}
