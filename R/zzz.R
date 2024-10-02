
#' @importFrom OmnipathR set_loglevel
#' @importFrom logger log_info
#' @importFrom utils packageVersion
#' @noRd
.onLoad <- function(libname, pkgname){

    OmnipathR:::init_config(pkg = pkgname)
    OmnipathR:::init_log(pkg = pkgname)
    log_info('Welcome to MetaProViz!')
    log_info('MetaProViz version: %s', packageVersion(pkgname))

    buildserver <- OmnipathR:::.on_buildserver()

    if(buildserver){

        set_loglevel('trace', target = 'console', pkg = pkgname)

    }

}
