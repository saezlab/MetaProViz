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


#' Path to the current MetaProViz log file
#'
#' @return Character: path to the current logfile, or \code{NULL} if no logfile is
#'     available.
#'
#' @examples
#' metaproviz_logfile()
#' # [1] "path/metaproviz/metaproviz-log/metaproviz-20210309-1642.log"
#'
#'
#' @seealso
#' \code{\link{metaproviz_log}}
#'
#' @importFrom OmnipathR logfile
#' @export
metaproviz_logfile <- function() {
    # Creates the path for the log file
    logfile("MetaProViz")
}


#' Browse the current MetaProViz log file
#'
#' @return Returns `NULL`.
#'
#' @examples
#' \dontrun{
#' metaproviz_log()
#' # then you can browse the log file, and exit with `q`
#' }
#'
#'
#' @seealso
#' \code{\link{metaproviz_logfile}}
#'
#' @importFrom OmnipathR read_log
#' @export
metaproviz_log <- function() {
    read_log("MetaProViz")
}


#' Sets the log level for the package logger
#'
#' @param level Character or class `loglevel`. The desired log level.
#' @param target Character, either 'logfile' or 'console'
#'
#' @return Returns `NULL`.
#'
#' @examples
#' metaproviz_set_loglevel(logger::FATAL, target = "console")
#'
#'
#' @importFrom OmnipathR set_loglevel
#' @export
metaproviz_set_loglevel <- function(
    level,
    target = "logfile"
) {
    # to change log-level e.g. to see all messages being printed
    set_loglevel(level, target = target, pkg = "MetaProViz")
}


#' Set the console log level to "trace"
#'
#' @noRd
mpv_trace <- function() {  # MPV=MetaProViz

    # Useful for debugging
    metaproviz_set_loglevel("trace", target = "console")
}


#' Set the console log level to "trace"
#'
#' @importFrom OmnipathR set_loglevel
#' @importFrom logger log_info
#' @importFrom utils packageVersion
#' @noRd
metaproviz_init <- function() {
    # NSE vs. R CMD check workaround
    OmnipathR <- init_config <- init_log <- .on_buildserver <- NULL
    # Only run the first time a MetaProViz function is used in an environment --> Creates log file
    if (is.null(metaproviz.env$init)) {
        # Initial log
        (OmnipathR %:::% init_config)(pkg = "MetaProViz")
        (OmnipathR %:::% init_log)(pkg = "MetaProViz")
        log_info("Welcome to MetaProViz!")
        log_info("MetaProViz version: %s", packageVersion("MetaProViz"))

        # Do we run in a build server like bioconductor?
        buildserver <- (OmnipathR %:::% .on_buildserver)()

        if (buildserver) {
            set_loglevel("trace", target = "console", pkg = "MetaProViz")
        }

        metaproviz.env$init <- TRUE
    }
}
