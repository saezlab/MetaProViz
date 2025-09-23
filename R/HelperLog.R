## ---------------------------
##
## Script name: HelperFunctions
##
## Purpose of script: Log files
##
## Author: Denes Turei
##
## Date Created: 2024-06-14
##
## Copyright (c) Denes Turei
## Email:
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------


#' Path to the current MetaProViz log file
#'
#' @return Character: path to the current logfile, or \code{NULL} if no
#'     logfile is available.
#'
#' @examples
#' metaproviz_logfile()
#' # [1] "path/metaproviz/metaproviz-log/metaproviz-20210309-1642.log"
#'
#' @importFrom OmnipathR logfile
#'
#' @export
#'
#' @seealso \code{\link{metaproviz_log}}
metaproviz_logfile <- function(){
  #Creates the path for the log file
  OmnipathR::logfile('MetaProViz')
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
#' @importFrom OmnipathR read_log
#'
#' @export
#' @seealso \code{\link{metaproviz_logfile}}
metaproviz_log <- function(){

  OmnipathR::read_log('MetaProViz')

}


#' Sets the log level for the package logger
#'
#' @param level Character or class `loglevel`. The desired log level.
#' @param target Character, either 'logfile' or 'console'
#'
#' @return Returns `NULL`.
#'
#' @examples
#' metaproviz_set_loglevel(logger::FATAL, target = 'console')
#'
#' @importFrom OmnipathR set_loglevel
#' @export
metaproviz_set_loglevel <- function(level, target = 'logfile'){
  #to change log-level e.g. to see all messages being printed
  OmnipathR::set_loglevel(level, target = target, pkg = 'MetaProViz')
}


#' Set the console log level to "trace"
#'
#' @noRd
mpv_trace <- function() {#MPV=MetaProViz
  # Useful for debugging
  metaproviz_set_loglevel('trace', target = 'console')
}

#' Set the console log level to "trace"
#'
#' @importFrom OmnipathR set_loglevel
#' @importFrom logger log_info
#' @importFrom utils packageVersion
#'
#' @noRd
#'
metaproviz_init <- function(){
  # Only run the first time a MetaProViz function is used in an environment --> Creates log file
  if(is.null(metaproviz.env$init)){
    # Initial log
    (OmnipathR%:::%init_config)(pkg = "MetaProViz")
    (OmnipathR%:::%init_log)(pkg = "MetaProViz")
    log_info('Welcome to MetaProViz!')
    log_info('MetaProViz version: %s', packageVersion("MetaProViz"))

    # Do we run in a build server like bioconductor?
    buildserver <- (OmnipathR%:::%.on_buildserver)()

    if(buildserver){
      set_loglevel('trace', target = 'console', pkg = "MetaProViz")
    }

    metaproviz.env$init <- TRUE
  }
}
