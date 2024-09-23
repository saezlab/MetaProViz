
#' A built in palette from MetaProViz
#'
#' @param name Character: name of the palette.
#'
#' @return Character vector of colors.
#'
#' @importFrom yaml read_yaml
#' @importFrom magrittr %>%
#' @noRd
metaproviz_palette <- function(name = "default") {

    system.file(
        'palettes',
        sprintf('%s.yaml', name),
        package = 'MetaProViz',
        mustWork = TRUE
    ) %>%
    read_yaml %>%
    unlist %>%
    unname

}
