#' Set the height or width of a row or column in a TableGrob (gtable)
#'
#' @param gtbl A TableGrob (gtable) object.
#' @param name Character or integer: Name or index of the element or row.
#' @param height Unit: the new height with its unit, e.g. ``unit(2.3, "cm")``.
#' @param dim Character: either ``'heights'`` or ``'widths'``.
#' @param offset Integer: offset from ``name``.
#' @param ifempty Logical: set the size only if the index looked up by ``name``
#'     and ``offset`` can not be found in the ``l`` or ``t`` column of the
#'     layout table, which means, no grob starts in that column or row.
#' @param callback Function: called with the new and the original value as
#'     arguments, and its value will override the new value. By default it is
#'     ``partial(switch, TRUE)``, which ignores the original value.
#'
#' @importFrom magrittr %>% extract2 add %<>%
#' @importFrom dplyr filter
#' @importFrom logger log_warn log_trace
#' @importFrom purrr partial
#' @noRd
set_size <- function(
        gtbl,
        name,
        size,
        dim,
        offset = 0L,
        ifempty = offset != 0L,
        callback = partial(switch, TRUE)
    ) {

    size %<>% parse_unit
    col <- list(heights = 't', widths = 'l') %>% extract2(dim)

    idx <-
        name %>%
        {`if`(
            is.character(.),
            filter(gtbl$layout, name == .) %>% extract2(col),
            .
        )} %>%
        add(offset)

    name_miss <- length(idx) == 0L

    outof_range <-
        idx < min(gtbl$layout[[col]]) ||
        idx > max(gtbl$layout[[col]])

    info <-
        sprintf(
            '[name=%s,offset=%i,empty=%s,original=%s]',
            name,
            offset,
            !(idx %in% gtbl$layout[[col]]),
            `if`(name_miss || outof_range, 'NA', gtbl[[dim]][idx])
        )

    if (name_miss) {

        log_warn(
            'No such name in gtable: %s; names available: %s',
            name,
            paste0(gtbl$layout$name, collapse = ', ')
        )

    } else if (outof_range) {

        log_warn(
            'Index %i is out of range for %s[%i-%i] %s',
            idx,
            dim,
            min(gtbl$layout[[col]]),
            max(gtbl$layout[[col]]),
            info
        )

    } else if (!ifempty || !(idx %in% gtbl$layout[[col]])) {

        size %<>% callback(gtbl[[dim]][idx])
        log_trace('Setting %s[%i] to %s %s', dim, idx, size, info)
        gtbl[[dim]][idx] <- size

    } else {

        log_trace('Not setting %s[%i] to %s %s', dim, idx, size, info)

    }

    invisible(gtbl)

}


#' Set the height of a row in a TableGrob (gtable)
#'
#' @importFrom purrr partial
#' @noRd
set_height <- partial(set_size, dim = 'heights')


#' Set the width of a column in a TableGrob (gtable)
#'
#' @importFrom purrr partial
#' @noRd
set_width <- partial(set_size, dim = 'widths')

#' Parse a unit from a string
#'
#' @param u Character: unit as a string, e.g. "1in".
#'
#' @importFrom magrittr %>%
#' @importFrom stringr str_trim str_match
#' @importFrom rlang !!! exec
#' @importFrom grid unit
#' @noRd
parse_unit <- function(u) {

    u %>%
    {`if`(
        is.character(.),
        str_trim(.) %>%
        str_match('^([0-9.]+)([a-z]+)$') %>%
        {`if`(
            !any(is.na(.)),
            exec(unit, !!!.[-1L]),
            {
                log_warn('Could not parse unit %s', u)
                unit(1, 'null')
            }
        )},
        .
    )}

}
