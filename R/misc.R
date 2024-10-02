
#' Use the second value if the first one is NULL
#'
#' @importFrom rlang %||%
#' @noRd
if_null <- function(value1, value2) {

    value1 %||% value2

}


#' Convert string lengths to centimeters
#'
#' @noRd
char2cm <- function(char) {

  nchar(char) * .25 + .8

}


#' Convert numeric to cm unit
#'
#' @importFrom grid unit
#' @noRd
cm <- function(value) {

  unit(value, 'cm')

}
