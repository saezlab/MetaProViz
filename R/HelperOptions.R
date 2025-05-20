
.metaproviz_options_defaults <- list(
  metaproviz.loglevel = 'trace',
  metaproviz.console_loglevel = 'success'
)


#' Save the current package configuration
#'
#' @param path Path to the config file. Directories and the file will be
#'     created if don't exist.
#' @param title Save the config under this title. One config file might
#'     contain multiple configurations, each identified by a title.
#' @param local Save into a config file in the current directory instead of
#'     a user level config file. When loading, the config in the current
#'     directory has priority over the user level config.
#'
#' @return Returns `NULL`.
#'
#' @examples
#' \dontrun{
#' # after this, all downloads will default to commercial licenses
#' # i.e. the resources that allow only academic use will be excluded:
#' options(metaproviz.console_loglevel = 'trace')
#' metaproviz_save_config()
#' }
#'
#' @importFrom OmnipathR save_config
#' @export
metaproviz_save_config <- function(
        path = NULL,
        title = 'default',
        local = FALSE
    ){

    metaproviz_save_config(path = path, title = title, local = local, pkg = 'MetaProViz')

}


#' Load the package configuration from a config file
#'
#' @param path Path to the config file.
#' @param title Load the config under this title. One config file might
#'     contain multple configurations, each identified by a title. If the
#'     title is not available the first section of the config file will be
#'     used.
#' @param user Force to use the user level config even if a config file
#'     exists in the current directory. By default, the local config files
#'     have prioroty over the user level config.
#' @param ... Passed to \code{yaml::yaml.load_file}.
#'
#' @return Invisibly returns the config as a list.
#'
#' @examples
#' \dontrun{
#' # load the config from a custom config file:
#' metaproviz_load_config(path = 'my_custom_metaproviz_config.yml')
#' }
#'
#' @importFrom OmnipathR load_config
#' @export
metaproviz_load_config <- function(
        path = NULL,
        title = 'default',
        user = FALSE,
        ...
    ){

    metaproviz_load_config(
        path = path,
        title = title,
        user = user,
        pkg = 'MetaProViz',
        ...
    )

}


#' Restore the built-in default values of all config parameters of MetaProViz
#'
#' @param save If a path, the restored config will be also saved
#'     to this file. If TRUE, the config will be saved to the current default
#'     config path (see \code{\link{metaproviz_config_path}}).
#' @param reset_all Reset to their defaults also the options already set in
#'     the R options.
#'
#' @examples
#' \dontrun{
#' # restore the defaults and write them to the default config file:
#' metaproviz_reset_config()
#' metaproviz_save_config()
#' }
#'
#' @return The config as a list.
#'
#' @export
#' @importFrom OmnipathR reset_config
#' @seealso \code{\link{metaproviz_load_config}, \link{metaproviz_save_config}}
metaproviz_reset_config <- function(save = NULL, reset_all = FALSE) {

   reset_config(save = save, reset_all = reset_all, pkg = 'MetaProViz')

}


#' Current config file path of MetaProViz
#'
#' @param user Logical: prioritize the user level config even if a config in
#'     the current working directory is available.
#'
#' @return Character: path to the config file.
#'
#' @examples
#' metaproviz_config_path()
#'
#' @importFrom OmnipathR config_path
#' @export
metaproviz_config_path <- function(user = FALSE){

    config_path(user = user, pkg = 'MetaProViz')

}
