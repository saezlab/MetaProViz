#' Set the height or width of a row or column in a TableGrob (gtable)
#'
#' @param gtbl A TableGrob (gtable) object.
#' @param name Character or integer: Name or index of the element or row.
#' @param size Character, numeric or unit: the new size optinally with its
#'     unit, e.g. ``2.3`` or ``2.3cm`` ``unit(2.3, "cm")``.
#' @param dim Character: either ``'heights'`` or ``'widths'``.
#' @param offset Integer: offset from ``name``.
#' @param ifempty Logical: set the size only if the index looked up by ``name``
#'     and ``offset`` can not be found in the ``l`` or ``t`` column of the
#'     layout table, which means, no grob starts in that column or row.
#' @param callback Function: called with the new and the original value as
#'     arguments, and its value will override the new value. By default it is
#'     ``partial(switch, TRUE)``, which ignores the original value.
#' @param grow Logical or numeric: if the total width or height is stored in
#'     the plottable in the slots added by ``withCanvasSize``, increase the
#'     relevant dimension. By default, the increment is identical to the
#'     increase of the size adjusted here; if a number is provided, it will be
#'     added to the existing value.
#'
#' @importFrom magrittr %>% extract2 add %<>%
#' @importFrom dplyr filter
#' @importFrom logger log_warn log_trace
#' @importFrom purrr partial
#' @importFrom stringr str_sub
#' @noRd
set_size <- function(
        gtbl,
        name,
        size,
        dim,
        offset = 0L,
        ifempty = offset != 0L,
        callback = partial(switch, TRUE),
        grow = FALSE
    ) {

    callback %<>% {`if`(is.character(.), get(.), .)}
    size %<>% parse_unit
    col <- dim %>% gtable_col
    tdim <- dim %>% str_sub(end = -2L)

    idx <- gtable_idx(gtbl, name, offset, dim)

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

        original <- gtbl[[dim]][idx]
        size %<>% callback(original)

        if (grow && tdim %in% names(gtbl) && size != original) {
            grow %<>% {`if`(is.numeric(.), ., size - original)}
            log_trace(
                paste0(
                    'Adding %s to the total %s of the gtable, ',
                    'growing it from %s to %s.'
                ),
                grow, tdim, gtbl[[tdim]], gtbl[[tdim]] + grow
            )
            gtbl[[tdim]] %<>% add(grow)
        }

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
#' @importFrom grid unit
#' @importFrom logger log_warn
#' @noRd
parse_unit <- function(u) {

    u %>%
    {`if`(
        is.character(.),
        str_trim(.) %>%
        str_match('^(-?[0-9.]+)([a-z]+)$') %>%
        {`if`(
            !any(is.na(.)),
            unit(.[,2L], .[,3L]),
            {
                log_warn('Could not parse unit %s', paste0(u, collapse = ', '))
                unit(1, 'null')
            }
        )},
       `if`(is.numeric(.), unit(., 'cm'), .)
    )}

}


#' Little subclass of gtable that carries the canvas width and height
#'
#' @importFrom magrittr %>%
#' @noRd
withCanvasSize <- function(gtable, width = NULL, height = NULL) {

    cls <- class(gtable) %>% c('withCanvasSize')

    gtable %>%
    c(
        list(
            width = parse_unit(width),
            height = parse_unit(height)
        )
    ) %>%
    `class<-`(cls)

}


#' Numeric row or column index in gtable from name
#'
#' @importFrom dplyr filter
#' @importFrom magrittr %>% add
#' @noRd
gtable_idx <- function(gtbl, name, dim, offset = 0L) {

    dim %<>% gtable_col

    name %>%
    {`if`(
        is.character(.),
        filter(gtbl$layout, name == .) %>% extract2(dim),
        .
    )} %>%
    add(offset)

}


#' Relevant column in gtable$layout for a dimension ("widths" or "heights")
#'
#' @importFrom magrittr %>% extract2
#' @importFrom rlang %||%
#' @noRd
gtable_col <- function(dim) {

    list(heights = 't', widths = 'l') %>% extract2(dim) %||% dim

}
