## ---------------------------
##
## Script name: Helper for figures
##
## Purpose of script: Make ggplot figures nice
##
## Author:
##
## Date Created: 2024-09-28
##
## Copyright (c) Saez Lab
## Email:
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------
#'

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

    idx <-
        name %>%
        in_gtable(gtbl) %>%
        gtable_idx(gtbl, ., dim, offset = offset)

    name_miss <- length(idx) == 0L

    outof_range <-
        idx < min(gtbl$layout[[col]]) ||
        idx > max(gtbl$layout[[col]])

    info <-
        sprintf(
            '[name=%s,offset=%i,empty=%s,original=%s]',
            paste0(name, collapse = ','),
            offset,
            !(idx %in% gtbl$layout[[col]]),
            `if`(name_miss || outof_range, 'NA', gtbl[[dim]][idx])
        )

    if (name_miss) {

        log_warn(
            'No such name in gtable: %s; names available: %s',
            paste0(name, collapse = ', '),
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

        if (grow && tdim %in% names(gtbl)) {
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


#' Apply simple width and height adjustments to gtable layout
#'
#' @importFrom magrittr %>%
#' @importFrom purrr reduce
#' @noRd
adjust_layout <- function(gtbl, param) {

    c('widths', 'heights') %>%
    reduce(
        ~set_sizes(.x, .y, param[[.y]]),
        .init = gtbl
    )

}


#' Apply simple adjustments along one dimension of a gtable layout
#'
#' @importFrom magrittr %>%
#' @importFrom purrr reduce
#' @importFrom rlang !!!
#' @noRd
set_sizes <- function(gtbl, dim, param) {

    param %>%
    reduce(
        ~exec(set_size, .x, !!!.y, dim = dim),
        .init = gtbl
    )

}


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

    log_trace('Parsing unit from `%s`', u)

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
       `if`(is.numeric(.), cm(.), .)
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
        filter(gtbl$layout, name %in% .) %>% extract2(dim),
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


#' Adjust gtable paremeters to accommodate title(s)
#'
#' @param gtbl A TableGrob (gtable) object.
#' @param titles Character: titles to be added to the plot. Most commonly we
#'     have one title, sometimes there is also a subtitle.
#'
#' @importFrom magrittr %<>% %>% add
#' @importFrom logger log_trace
#' @noRd
adjust_title <- function(gtbl, titles) {

    if (titles %>% nchar %>% as.logical %>% any) {

        log_trace('The plot has title, adjusting layout to accommodate it.')

        gtbl %<>% set_height(c('title', 'main'), '1.5cm')#controls margins --> PlotName and subtitle
        # Sum up total heights:
        gtbl$height %<>% add(cm(.5))

        #------- Width: Check how much width is needed for the figure title/subtitle
        title_width <- titles %>% char2cm %>% max %>% cm
        gtbl <<- gtbl
        gtbl %<>% set_width(
            c('guide-box-right', 'legend'),
            sprintf('%.02fcm', title_width - gtbl$width),
            callback = max
        )
        gtbl$width %<>% max(title_width)

    }

    invisible(gtbl)

}


#' Adjust gtable paremeters to accommodate legend(s)
#'
#' @importFrom ggpubr get_legend
#' @importFrom magrittr %<>% %>% extract add multiply_by
#' @importFrom logger log_trace
#' @importFrom purrr map_int
#' @noRd
adjust_legend <- function(
        gtbl,
        InputPlot,
        sections = FALSE,
        SettingsInfo = NULL
    ) {

    if(
        any(sections %in% names(SettingsInfo)) ||
        identical(sections, TRUE)
    ) {

        log_trace('The plot has legend, adjusting layout to accommodate it.')
        log_trace('Sections: %s', paste0(sections, collapse = ', '))

        Legend <- get_legend(InputPlot) # Extract legend to adjust separately
        leg <<- Legend

        #------- Legend widths
        ## Legend titles:
        legend_nchar <-
            SettingsInfo %>%
            as.list %>%
            extract(sections) %>%
            unlist %>%
            map_int(nchar) %>%
            max(0L) %>%
            multiply_by(.25)

        legend_width <-
            Legend$widths[3L] %>%
            as.numeric %>%
            round(1L) %>%
            max(legend_nchar)

        log_trace('Legend nchar: %.02fcm, Legend width: %.02fcm', legend_nchar, legend_width)

        ## Legend space:
        gtbl %<>%
            set_width(
                'guide-box-right',
                sprintf('%.02fcm', legend_width),
                callback = max,
                # here we have a bug, this should be TRUE
                grow = FALSE
            )

        #------- Legend heights
        Legend_heights <-
            Legend$heights %>%
            extract(c(3L, 5L)) %>%
            as.numeric %>%
            sum(2) #+2 to secure space above and below plot

        if(as.numeric(gtbl$height) < Legend_heights){

            Add <- (Legend_heights - as.numeric(gtbl$height)) / 2

            gtbl %<>%
                #controls margins --> Can be increased if Figure legend needs more space on the top
                set_height('background', cm(Add)) %>%
                #controls margins --> Can be increased if Figure legend needs more space on the bottom
                set_height('xlab-b', cm(Add), offset = 1L)

            gtbl$height <- cm(Legend_heights)

        }

    }

    invisible(gtbl)

}


#' From a legend gtable extract the titles into a character vector
#'
#' @importFrom magrittr %>% extract2
#' @importFrom purrr map map_chr
#' @noRd
titles_from_legend <- function(leg) {

    leg$grobs %>%
    map(~{(.x$layout$name == 'title') %>% which %>% extract2(.x$grobs, .)}) %>%
    map(~{.x$children %>% map_chr(~.x$label)}) %>%
    unlist %>%
    unname

}


#' @importFrom utils head
#' @importFrom magrittr %>%
#' @noRd
in_gtable <- function(name, gtbl) {

    name %>%
    {`if`(is.character(.), intersect(., gtbl$layout$name) %>% head(1L), .)}

}
