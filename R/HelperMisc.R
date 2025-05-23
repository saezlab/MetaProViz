#' Makes sure we have a string even if the argument was passed by NSE
#'
#' @importFrom magrittr %>%
#' @importFrom rlang enquo quo_get_expr quo_text is_symbol
#' @importFrom stringr str_remove
#' @noRd
.nse_ensure_str <- function(arg){

    enquo(arg) %>%
    {`if`(
        is_symbol(quo_get_expr(.)),
        quo_text(.),
        quo_get_expr(.)
    )} %>%
    str_remove('`$') %>%
    str_remove('^`')

}


#' Workaround against R CMD check notes about using `:::`
#'
#' @importFrom rlang enquo !!
#'
#' @noRd
`%:::%` <- function(pkg, fun){

    pkg <- .nse_ensure_str(!!enquo(pkg))
    fun <- .nse_ensure_str(!!enquo(fun))

    get(fun, envir = asNamespace(pkg), inherits = FALSE)

}
