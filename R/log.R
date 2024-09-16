
#' Path to the current MetaProViz log file
#'
#' @return Character: path to the current logfile, or \code{NULL} if no
#'     logfile is available.
#'
#' @examples
#' metaproviz_logfile()
#' # [1] "/home/denes/metaproviz/metaproviz-log/metaproviz-20210309-1642.log"
#'
#' @importFrom OmnipathR logfile
#' @export
#' @seealso \code{\link{metaproviz_log}}
metaproviz_logfile <- function(){

    logfile('MetaProViz')

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
#' @export
#' @seealso \code{\link{metaproviz_logfile}}
metaproviz_log <- function(){

    read_log('MetaProViz')

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
#' @export
metaproviz_set_loglevel <- function(level, target = 'logfile'){

    set_loglevel(level, target = target, pkg = 'MetaProViz')

}


#' Set the console log level to "trace"
#'
#' @noRd
.mpvtrace <- function() {

    metaproviz_set_loglevel('trace', target = 'console')

}
