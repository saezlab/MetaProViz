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
        ~rlang::exec(set_size, .x, !!!.y, dim = dim),
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
#'     have one title and a subtitle
#'
#' @importFrom magrittr %<>% %>% add
#' @importFrom logger log_trace
#' @noRd
adjust_title <- function(gtbl, titles) {


    if (titles %>% nchar %>% as.logical %>% any) {

        log_trace('The plot has title, adjusting layout to accommodate it.')

        gtbl %<>% set_height(c('title', 'main'), '0.5cm')#controls margins -->plot_name
        gtbl %<>% set_height(c('subtitle', 'main'), '0.5cm')#controls margins -->plot_name


        # Sum up total heights:
        gtbl$height %<>% add(cm(1))

        #------- Width: Check how much width is needed for the figure title/subtitle
        title_width <- titles %>% char2cm %>% max %>% cm

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
        metadata_info = NULL
    ) {

    if(
        any(sections %in% names(metadata_info)) ||
        identical(sections, TRUE)
    ) {

        log_trace('The plot has legend, adjusting layout to accommodate it.')
        log_trace('Sections: %s', paste0(sections, collapse = ', '))

        Legend <- get_legend(InputPlot) # Extract legend to adjust separately

        #------- Legend widths
        ## Legend titles:
        legend_nchar <-
            metadata_info %>%
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
                grow = TRUE
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


##############################################################
### ### ### Plot helper function: Processing ### ### ###
##############################################################

#' @param InputPlot This is the ggplot object generated within the in any of the processing functions function.
#' @param plot_name Generated within the processing functions.
#' @param plot_type Generated within the processing functions.
#'
#' @keywords Plot helper function
#' @importFrom ggplot2 annotation_custom theme
#' @importFrom magrittr %>% %<>% add
#' @noRd
#'

plotGrob_Processing <- function(InputPlot,plot_name, plot_type){

  if(plot_type == "Scree"){
    UNIT <- unit(12, "cm")
  }else if(plot_type == "Hotellings"){
    UNIT <- unit(12, "cm")#0.25*Sample
  }else{#CV and Hist
    UNIT <- unit(8, "cm")
  }
  # Make plot into nice format:
  SUPER_PARAM <- list(widths = list(
    list("axis-b", UNIT),
    list("ylab-l", "0cm", offset = -4L, ifempty = FALSE),
    list("axis-l", "1cm"),
    list("ylab-l", "1cm"),
    list("guide-box-left", "0cm"),
    list("axis-r", "0cm"),
    list("ylab-r", "0cm"),
    list("ylab-l", "1cm", offset = -1L),
    list("guide-box-right", "1cm")
  ),
  heights = list(
    list("axis-l", "8cm"),
    list("axis-b", "0.5cm"),#This is adjusted for the x-axis ticks!
    list("xlab-b", "0.75cm"),#This gives us the distance of the caption to the x-axis label
    list("title", "0cm", offset = -2L, ifempty = FALSE),
    list("title", "0cm", offset = -1L),
    list("title", "0.25cm"),# how much space is between title and y-axis label
    list("subtitle", "0cm"),
    list("caption", "0.5cm"), #plots statistics information, space to bottom
    list("guide-box-top", "0cm"),
    list("xlab-t", "0cm", offset = -1L)
  )
  )

  #Adjust the parameters:
  suppressWarnings(suppressMessages(
    Plot_Sized <- InputPlot %>%
      ggplotGrob %>%
      withCanvasSize(width = 12, height = 11) %>%
      adjust_layout(SUPER_PARAM) %>%
      adjust_title(c(PlotName))
  ))

  Plot_Sized %<>%
    {ggplot2::ggplot() + annotation_custom(.)} %>%
    add(theme(panel.background = ggplot2::element_rect(fill = "transparent")))

  return(Plot_Sized)
}

##############################################################
### ### ### Plot helper function: PCA ### ### ###
##############################################################

#' PCA helper function: Internal Function
#'
#' @param InputPlot This is the ggplot object generated within the VizPCA function.
#' @param metadata_info Passed to VizPCA
#' @param plot_name Passed to VizPCA
#'
#' @keywords PCA helper function
#' @importFrom ggplot2 ggplotGrob
#' @importFrom magrittr %>%
#' @noRd
PlotGrob_PCA <- function(InputPlot, metadata_info,plot_name){

  PCA_PARAM <- list(
    widths = list(
      list("axis-b", "8cm"),
      list("ylab-l", "0cm", offset = -4L, ifempty = FALSE),
      list("axis-l", "1cm"),
      list("ylab-l", "1cm"),
      list("guide-box-left", "0cm"),
      list("axis-r", "0cm"),
      list("ylab-r", "0cm"),
      list("ylab-l", "1cm", offset = -1L),
      list("guide-box-right", "1cm")
    ),
    heights = list(
      list("axis-l", "8cm"),
      list("axis-b", "1cm"),
      list("xlab-b", ".5cm"),
      list("xlab-b", "1cm", offset = 1L),
      list("title", "0cm", offset = -2L, ifempty = FALSE),
      list("title", "0cm", offset = -1L),
      list("title", "0.25cm"),# how much space is between title and y-axis label
      list("subtitle", "0cm"),
      list("guide-box-top", "0cm"),
      list("xlab-t", "0cm", offset = -1L)
    )
  )

  Plot_Sized <- InputPlot %>%
    ggplotGrob %>%
    withCanvasSize(width = 12, height = 11) %>%
    adjust_layout(PCA_PARAM) %>%
    adjust_title(PlotName) %>%
    adjust_legend(
      InputPlot,
      sections = c("color", "shape"),
      metadata_info = metadata_info
    )

  log_trace(
    'Sum of heights: %.02f, sum of widths: %.02f',
    grid::convertUnit(sum(Plot_Sized$height), 'cm', valueOnly = TRUE),
    grid::convertUnit(sum(Plot_Sized$width), 'cm', valueOnly = TRUE)
  )

  #Return
  Output <- Plot_Sized
}


##############################################################
### ### ### Plot helper function: Heatmap ### ### ###
##############################################################

#' @param InputPlot This is the ggplot object generated within the VizHeatmap function.
#' @param metadata_info Passed to VizHeatmap
#' @param metadata_sample Passed to VizHeatmap
#' @param metadata_feature Passed to VizHeatmap
#' @param plot_name Passed to VizHeatmap
#'
#' @keywords Heatmap helper function
#' @noRd

PlotGrob_Heatmap <- function(InputPlot, metadata_info, metadata_sample, metadata_feature,plot_name){

  # Set the parameters for the plot we would like to use as a basis, before we start adjusting it:
  HEAT_PARAM <- list(
    widths = list(
      list("legend", "2cm")
    ),
    heights = list(
      list("main", "1cm")
    )
  )

  #If we plot feature names on the x-axis, we need to adjust the height of the plot:
  #if(as.logical(show_rownames)==TRUE){
  #  Rows <- nrow(t(data))
  #}

  #Adjust the parameters:
  Input <- InputPlot$gtable

  Plot_Sized <- Input  %>%
    withCanvasSize(width = 12, height = 11) %>%
    adjust_layout(HEAT_PARAM) %>%
    adjust_title(c(PlotName))

  #Extract legend information and adjust:
  color_entries <- grep("^color", names(metadata_info), value = TRUE)
  if(length(color_entries)>0){#We need to adapt the plot Hights and widths
    if(sum(grepl("color_Sample", names(metadata_info)))>0){
      names <- metadata_info[grepl("color_Sample", names(metadata_info))]
      colour_names <- NULL
      legend_names <- NULL
      for (x in 1:length(names)){
        names_sel <- names[[x]]
        legend_names[x] <- names_sel
        colour_names[x] <- metadata_sample[names[[x]]]
      }
    }else{
      colour_names <- NULL
      legend_names <- NULL
    }

    if(sum(grepl("color_Metab", names(metadata_info)))>0){
      names <- metadata_info[grepl("color_Metab", names(metadata_info))]
      colour_names_M <- NULL
      legend_names_M <- NULL
      for (x in 1:length(names)){
        names_sel <- names[[x]]
        legend_names_M[x] <- names_sel
        colour_names_M[x] <- metadata_feature[names[[x]]]
      }
    }else{
      colour_names_M <- NULL
      legend_names_M <- NULL
    }

    legend_head <- c(legend_names, legend_names_M)
    longest_name <- legend_head[which.max(nchar(legend_head[[1]]))]
    character_count_head <- nchar(longest_name)+4#This is the length of the legend title name

    legend_names <- c(unlist(colour_names), unlist(colour_names_M))
    longest_name <- legend_names[which.max(nchar(legend_names[[1]]))]
    character_count <- nchar(longest_name)#This is the length of the legend colour names

    legendWidth <-  unit(((max(character_count_head, character_count))*0.3), "cm")#legend space

    # Sum up total heights:
    Plot_Sized$width %<>% add(legendWidth)

    legendHeights <- unit((sum(length(unique(legend_names_M))+length(unique(colour_names_M))+length(unique(legend_names))+length(unique(colour_names)))), "cm")

    Plot_Sized$width %<>% add(legendWidth)
    if((grid::convertUnit(legendHeights, 'cm', valueOnly = TRUE))>(grid::convertUnit(Plot_Sized$height, 'cm', valueOnly = TRUE))){
      Plot_Sized$height <- legendHeights
    }

  }

  Output <- Plot_Sized

}


##############################################################
### ### ### Plot helper function: Volcano   ### ### ###
##############################################################

#' @param InputPlot This is the ggplot object generated within the VizVolcano function.
#' @param metadata_info Passed to VizVolcano
#' @param plot_name Passed to VizVolcano
#' @param Subtitle
#'
#' @keywords Volcano helper function
#' @noRd
#'

plotGrob_Volcano <- function(InputPlot, metadata_info,plot_name, Subtitle){
  # Set the parameters for the plot we would like to use as a basis, before we start adjusting it:
  VOL_PARAM <- list(
    widths = list(
      list("axis-b", "6cm"),
      list("ylab-l", "0cm", offset = -4L, ifempty = FALSE),
      list("axis-l", "1cm"),
      list("ylab-l", "1cm"),
      list("guide-box-left", "0cm"),
      list("axis-r", "0cm"),
      list("ylab-r", "0cm"),
      list("ylab-l", "1cm", offset = -1L),
      list("guide-box-right", "1cm")
    ),
    heights = list(
      list("axis-l", "8cm"),
      list("axis-b", "0.75cm"),#This is the distance to x-axis!
      list("xlab-b", "0.75cm"),#This gives us the distance of the caption to the x-axis label
      #list("xlab-b", "1cm", offset = 1L),
      list("title", "0cm", offset = -2L, ifempty = FALSE),
      list("title", "0cm", offset = -1L),#Space above title
      list("title", "0.25cm"),# how much space is between title and y-axis label
      list("subtitle", "0cm"),
      list("guide-box-top", "0cm"),
      list("xlab-t", "0cm", offset = -1L)
    )
  )

  #Adjust the parameters:
  Plot_Sized <- InputPlot %>%
    ggplotGrob %>%
    withCanvasSize(width = 12, height = 11) %>%
    adjust_layout(VOL_PARAM) %>%
    adjust_title(c(PlotName, Subtitle)) %>%#Fix this (if there is no Subtitle!)
    adjust_legend(
      InputPlot,
      sections = c("color", "shape"),
      metadata_info = metadata_info
    )

  log_trace(
    'Sum of heights: %.02f, sum of widths: %.02f',
    grid::convertUnit(sum(Plot_Sized$height), 'cm', valueOnly = TRUE),
    grid::convertUnit(sum(Plot_Sized$width), 'cm', valueOnly = TRUE)
  )

  #Return
  Output <- Plot_Sized
}


#####################################################################
### ### ### Plot helper function: Superplots  ### ### ###
#####################################################################

#' @param InputPlot This is the ggplot object generated within the VizSuperplots function.
#' @param metadata_info Passed to VizSuperplots
#' @param metadata_sample Passed to VizSuperplots
#' @param Subtitle Passed to VizSuperplots
#' @param plot_name Passed to VizSuperplots
#' @param plot_type Passed to VizSuperplots
#'
#' @keywords PCA helper function
#' @noRd

plotGrob_Superplot <- function(InputPlot,
                               metadata_info,
                               metadata_sample,
                               Subtitle,
                              plot_name,
                               plot_type){
  # Set the parameters for the plot we would like to use as a basis, before we start adjusting it:
  X_Con <- metadata_sample%>%
    dplyr::distinct(Conditions)

  X_Tick <- unit(X_Con[[1]] %>% char2cm %>% max * 0.6, "cm")

  if(plot_type == "Bar"){
    UNIT <- unit(X_Con%>%nrow() * 0.5, "cm")
  }else{
    UNIT <- unit(X_Con%>%nrow() * 1, "cm")
  }

  SUPER_PARAM <- list(
    widths = list(
      list("axis-b", paste(UNIT)),
      list("ylab-l", "0cm", offset = -4L, ifempty = FALSE),
      list("axis-l", "1cm"),
      list("ylab-l", "1cm"),
      list("guide-box-left", "0cm"),
      list("axis-r", "0cm"),
      list("ylab-r", "0cm"),
      list("ylab-l", "1cm", offset = -1L),
      list("guide-box-right", "1cm")
    ),
    heights = list(
      list("axis-l", "8cm"),
      list("axis-b", X_Tick),#This is adjusted for the x-axis ticks!
      list("xlab-b", "0.75cm"),#This gives us the distance of the caption to the x-axis label
      list("title", "0cm", offset = -2L, ifempty = FALSE),
      list("title", "0cm", offset = -1L),
      list("title", "0.25cm"),# how much space is between title and y-axis label
      list("subtitle", "0cm"),
      list("caption", "0.5cm"), #plots statistics information, space to bottom
      list("guide-box-top", "0cm"),
      list("xlab-t", "0cm", offset = -1L)
    )
  )

  #Adjust the parameters:
  Plot_Sized <- InputPlot %>%
    ggplotGrob %>%
    withCanvasSize(width = 12, height = 11) %>%
    adjust_layout(SUPER_PARAM) %>%
    adjust_title(c(PlotName, Subtitle)) %>%
    adjust_legend(
      InputPlot,
      sections = c("Superplot"),#here we do not have colour and shape, but other parameters
      metadata_info = metadata_info
    )

  #Return
  Output <- Plot_Sized
}



