# !/usr/bin/env Rscript

#
#  This file is part of the `MetaProViz` R package
#
#  Copyright 2022-2025
#  Saez Lab, Heidelberg University
#
#  Authors: see the file `README.md`
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file `LICENSE` or copy at
#      https://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: https://saezlab.github.io/MetaProViz
#  Git repo: https://github.com/saezlab/MetaProViz
#

#' This script allows you to perform different visualizations (bar, box, violin
#' plots) using the results of the MetaProViz analysis
#'


#' Bar, Box or Violin plot in Superplot style visualization
#'
#' @param data SummarizedExperiment (se) file including assay and colData. If se file
#'     is provided, metadata_sample is extracted from the colData of the se
#'     object. metadata_feature, if available, are extracted from the rowData.
#'     Alternatively provide a DF with unique sample identifiers as row names
#'     and metabolite numerical values in columns with metabolite identifiers
#'     as column names. Use NA for metabolites that were not detected.
#' @param metadata_sample \emph{Optional: } Only required if you did not provide se file in
#'     parameter data. Provide DF which contains metadata information about the
#'     samples, which will be combined with your input data based on the unique
#'     sample identifiers used as rownames. \strong{Default = NULL}
#' @param metadata_info Named vector including at least information on the conditions column:
#'     c(Conditions="ColumnName_metadata_sample"). Additionally Superplots can
#'     be made by adding Superplot ="ColumnName_metadata_sample", which are
#'     usually biological replicates or patient IDs. \strong{Default =
#'     c(Conditions="Conditions", Superplot = NULL)}
#' @param plot_type String with the information of the Graph style. Available options are
#'     Bar. Box and Violin  \strong{Default = Box}
#' @param plot_name \emph{Optional: } String which is added to the output files of the plot.
#' @param plot_conditions Vector with names of selected Conditions for the plot. Can also be used
#'     to order the Conditions in the way they should be displayed on the
#'     x-axis of the plot. \strong{Default = NULL}
#' @param stat_comparison List of numeric vectors containing Condition pairs to
#'     compare based on the order of the plot_conditions vector. \strong{Default = NULL}
#' @param pval \emph{Optional: } String which contains an abbreviation of the selected
#'     test to calculate p.value. For one-vs-one comparisons choose t.test or
#'     wilcox.test , for one-vs-all or all-vs-all comparison choose aov
#'     (=anova) or kruskal.test \strong{Default = NULL}
#' @param padj \emph{Optional: } String which contains an abbreviation of the selected
#'     p.adjusted test for p.value correction for multiple Hypothesis testing.
#'     Search: ?p.adjust for more methods:"BH", "fdr", "bonferroni", "holm",
#'     etc.\strong{Default = NULL}
#' @param xlab \emph{Optional: } String to replace x-axis label in plot.
#'     \strong{Default = NULL}
#' @param ylab \emph{Optional: } String to replace y-axis label in plot.
#'     \strong{Default = NULL}
#' @param theme \emph{Optional: } Selection of theme for plot, e.g. theme_grey(). You
#'     can check for complete themes here:
#'     https://ggplot2.tidyverse.org/reference/ggtheme.html. \strong{Default =
#'     NULL}
#' @param color_palette \emph{Optional: } Provide customized color_palette in vector format.
#'     \strong{Default = NULL}
#' @param color_palette_dot \emph{Optional: } Provide customized color_palette in vector format.
#'     \strong{Default = NULL}
#' @param save_plot \emph{Optional: } Select the file type of output plots.
#'     Options are svg, pdf, png or NULL. \strong{Default = svg}
#' @param print_plot \emph{Optional: } TRUE or FALSE, if TRUE plots are saved
#'     as an overview of the results. \strong{Default = TRUE}
#' @param path \emph{Optional:} Path to the folder the results should be saved at.
#'     \strong{Default = NULL}
#'
#' @return List with two elements: Plot and Plot_Sized
#'
#' @examples
#' data(intracell_raw_se)
#' # only plot the first 2 metabolites
#' Res <- viz_superplot(data = intracell_raw_se[1:2, , drop = FALSE])
#'
#' data(intracell_raw)
#' Intra <- intracell_raw[, c(1:6)] %>% tibble::column_to_rownames("Code")
#' Res <- viz_superplot(
#'     data = Intra[, -c(1:3)],
#'     metadata_sample = Intra[, c(1:3)]
#' )
#'
#' @importFrom ggplot2 ggplot theme geom_violin stat_summary geom_boxplot
#' @importFrom ggplot2 position_dodge element_text theme_classic
#' @importFrom ggplot2 geom_bar labs scale_color_manual theme xlab ylab
#' @importFrom ggpubr stat_pvalue_manual stat_compare_means
#' @importFrom ggbeeswarm geom_beeswarm
#' @importFrom grid convertUnit
#' @importFrom dplyr rename select group_by summarise filter mutate n across
#' @importFrom tidyr separate unite
#' @importFrom magrittr %>% %<>%
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom logger log_trace log_info
#' @importFrom tidyselect all_of
#' @importFrom purrr map_dbl
#' @importFrom stats sd
#' @export
viz_superplot <- function(
    data,
    metadata_sample = NULL,
    metadata_info = c(Conditions = "Conditions", Superplot = NULL),
    plot_type = "Box",  # Bar, Box, Violin,
    plot_name = "",
    plot_conditions = NULL,
    stat_comparison = NULL,
    pval = NULL,
    padj = NULL,
    xlab = NULL,
    ylab = NULL,
    theme = NULL,
    color_palette = NULL,
    color_palette_dot = NULL,
    save_plot = "svg",
    print_plot = TRUE,
    path = NULL
) {

    # NSE vs. R CMD check workaround
    Conditions <- Superplot <- Intensity <- comparisons_rev <- NULL

    ## ------------ Create log file ----------- ##
    metaproviz_init()

    log_info("viz_superplot: Superplot visualization")

    ## ------------- Check SummarizedExperiment file ---------- ##
    input_data <- data
    if (inherits(data, "SummarizedExperiment")) {
        log_info('Processing input SummarizedExperiment object.')
        se_list <- process_se(data)
        data <- se_list$data
        metadata_sample <- se_list$metadata_sample
    }

    ## ------------ Check Input files ----------- ##
    # HelperFunction `check_param`
    check_param(
        data = data,
        metadata_sample = metadata_sample,
        metadata_feature = NULL,
        metadata_info = metadata_info,
        save_plot = save_plot,
        save_table = NULL,
        core = FALSE,
        print_plot = print_plot
    )

    # check_param` Specific
    if (is.null(metadata_info)) {
        message <- paste0("You must provide the column name for Conditions via metadata_info = c(Conditions = ColumnName) in order to plot the x-axis conditions.")
        log_trace(paste("Error ", message, sep = ""))
        stop(message)
    }

    if (!(plot_type %in% c("Box", "Bar", "Violin"))) {
        message <- paste0("plot_type must be either Box, Bar or Violin.")
        log_trace(paste("Error ", message, sep = ""))
        stop(message)
    }

    if (!is.null(plot_conditions)) {
        for (Condition in plot_conditions) {
            if (!(Condition %in% metadata_sample[[metadata_info[["Conditions"]]]])) {
                message <-
                    paste0(
                        "Check Input. The plot_conditions ",
                        Condition,
                        " were not found in the Conditions Column."
                    )
                log_trace(paste("Error ", message, sep = ""))
                stop(message)
            }
        }
    }

    if (!is.null(stat_comparison)) {
        for (Comp in stat_comparison) {
            # find the differences between the following two blocks
            # :'''(
            if (!is.null(plot_conditions)) {
                if (!(plot_conditions[Comp[1]] %in% metadata_sample[[metadata_info[["Conditions"]]]])) {
                    message <-
                        paste0(
                            "Check Input. The stat_comparison condition ",
                            Comp[1],
                            " is not found in the Conditions Column of the metadata_sample."
                        )
                    log_trace(paste("Error ", message, sep = ""))
                    stop(message)
                }
                if (!(plot_conditions[Comp[2]] %in% metadata_sample[[metadata_info[["Conditions"]]]])) {
                    message <-
                        paste0(
                            "Check Input. The stat_comparison condition ",
                            Comp[2],
                            " is not found in the Conditions Column of the metadata_sample."
                        )
                    log_trace(paste("Error ", message, sep = ""))
                    stop(message)
                }
            }
        }
    }

    if (is.null(color_palette)) {
        color_palette <- "grey"
    }

    ## ------------ Check Input metadata_info ----------- ##
    # 7. Check stat_comparison & plot_conditions
    if (is.null(plot_conditions)) {
        Number_Cond <- length(unique(tolower(metadata_sample[["Conditions"]])))
    if (Number_Cond<=2) {
        MultipleComparison <- FALSE
    } else {
        MultipleComparison <- TRUE
    }
    } else if (length(plot_conditions)>2) {
    MultipleComparison <- TRUE
    } else if (length(plot_conditions)<=2) {
    Number_Cond <- length(unique(tolower(metadata_sample[["Conditions"]])))
    if (Number_Cond<=2) {
        MultipleComparison <- FALSE
    } else {
        MultipleComparison <- TRUE
    }
    }

    if (!is.null(pval)) {
        if (MultipleComparison & (pval == "t.test" | pval == "wilcox.test")) {
        message <-
            paste0(
                "Check input. The selected pval option for Hypothesis testing,",
                pval,
                " is for multiple comparison, but you have only 2 conditions. Hence aov is performed."
            )
        log_trace(paste("Warning ", message, sep = ""))
        warning(message)
        pval <- "aov"
    } else if (!MultipleComparison & (pval == "aov" | pval == "kruskal.test")) {
        message <-
            paste0(
                "Check input. The selected pval option for Hypothesis testing,",
                pval,
                " is for multiple comparison, but you have only 2 conditions. Hence t.test is performed."
            )
        log_trace(paste("Warning ", message, sep = ""))
        warning(message)
        pval <- "t.test"
        }
    }

    if (is.null(pval) & !MultipleComparison) {
        pval <- "t.test"
    }

    if (is.null(pval) & MultipleComparison) {
        pval <- "aov"
    }

    STAT_padj_options <-
        c(
            "holm",
            "hochberg",
            "hommel",
            "bonferroni",
            "BH",
            "BY",
            "fdr",
            "none"
        )
    if (!is.null(padj)) {
        if (!(padj %in% STAT_padj_options)) {
            message <-
                paste0(
                    "Check input. The selected padj option for multiple Hypothesis testing correction is not valid. Please select NULL or one of the folowing: ",
                    paste(STAT_padj_options, collapse = ", "),
                    "."
                )
            log_trace(paste("Error ", message, sep = ""))
            stop(message)
        }
    }

    if (is.null(padj)) {
        padj <- "fdr"
    }

    ## ------------ Create Results output folder ----------- ##
    if (!is.null(save_plot)) {
        folder <- save_path(folder_name = paste(plot_type, "Plots", sep = ""),
                                    path = path)
    }
    log_info("viz_superplot results saved at ", folder)

    ##########################################################################
    ## ------------ Prepare Input ----------- ##
    metadata_sample %<>%
    rename("Conditions" = paste(metadata_info[["Conditions"]]) )

    if ("Superplot" %in% names(metadata_info)) {
        metadata_sample %<>%
        rename("Superplot" = paste(metadata_info[["Superplot"]]) )

        data_merge <-
            merge(
                metadata_sample[c("Conditions", "Superplot")],
                data,
                by = 0L
            )
        data_merge %<>% column_to_rownames("Row.names")
    } else {
        data_merge <- merge(metadata_sample[c("Conditions")], data, by = 0)
        data_merge %<>% column_to_rownames("Row.names")
    }

    # Rename the x and y lab if the information has been passed:
    if (is.null(xlab)) {  # use column name of x provided by user
        xlab <- bquote(.(as.symbol(metadata_info[["Conditions"]])))
    } else if (!is.null(xlab)) {
    xlab <- bquote(.(as.symbol(xlab)))
    }

    if (is.null(ylab)) {  # use column name of x provided by user
        ylab <- bquote(.(as.symbol("Intensity")))
    } else if (!is.null(ylab)) {
    ylab <- bquote(.(as.symbol(ylab)))
    }

    # Set the theme:
    if (is.null(theme)) {
        theme <- theme_classic()
    }

    ## ------------ Create plots ----------- ##
    # make a list for plotting all plots together
    PlotList <- list()  # Empty list to store all the plots
    PlotList_adaptedGrid <- list()  # Empty list to store all the plots

    for (i in colnames(data)) {
    # Prepare the dfs:

        dataMeans <-
        data_merge %>%
        select(i, Conditions) %>%
        group_by(Conditions) %>%
        summarise(
            across(
            all_of(i),
            list(mean = mean, sd = sd)
        )
        ) %>%
        as.data.frame()


    names(dataMeans)[2] <- "Intensity"

    if ("Superplot" %in% names(metadata_info)) {
        plotdata <- data_merge %>%
                        select(i, Conditions, Superplot) %>%
                        group_by(Conditions) %>%
                        as.data.frame()
    } else {
        plotdata <- data_merge %>%
                        select(i, Conditions) %>%
                        group_by(Conditions) %>%
                        as.data.frame()
    }
    names(plotdata)[1] <- c("Intensity")
    # Change conditions to factor
    plotdata$Conditions <- factor(plotdata$Conditions)

    # Take only selected conditions
    if (!is.null(plot_conditions)) {
        dataMeans %<>% filter(Conditions %in% plot_conditions)
        plotdata %<>% filter(Conditions %in% plot_conditions)
        plotdata$Conditions <- factor(plotdata$Conditions, levels = plot_conditions)
    }

    # Make the Plot
    Plot <- ggplot(plotdata, aes(x = Conditions, y = Intensity))

    # Add graph style and error bar
    data_summary <- function(
        x
    ) {
        m <- mean(x)
        ymin <- m-sd(x)
        ymax <- m+sd(x)
        return(c(y = m, ymin = ymin, ymax = ymax))
    }

    if (plot_type == "Bar") {
        Plot <- Plot+  geom_bar(stat = "summary", fun = "mean", fill = color_palette)+ stat_summary(fun.data = data_summary,
            geom = "errorbar", color = "black", width = 0.2)
    } else if (plot_type == "Violin") {
        Plot <- Plot+ geom_violin(fill = color_palette)+ stat_summary(fun.data = data_summary,
            geom = "errorbar", color = "black", width = 0.2)
    } else if (plot_type == "Box") {
        Plot <- Plot +  geom_boxplot(fill = color_palette, width = 0.5, position = position_dodge(width = 0.5))
    }

    # Add Superplot
    if ("Superplot" %in% names(metadata_info)) {
        if (!is.null(color_palette_dot)) {
            Plot <- Plot+ geom_beeswarm(aes(x = Conditions, y = Intensity, color = as.factor(Superplot)), size = 3)+
            labs(
                color = metadata_info[["Superplot"]],
                fill = metadata_info[["Superplot"]]
            )+
            scale_color_manual(values = color_palette_dot)
        } else {
        Plot <- Plot+ geom_beeswarm(aes(x = Conditions, y = Intensity, color = as.factor(Superplot)), size = 3)+
            labs(
                color = metadata_info[["Superplot"]],
                fill = metadata_info[["Superplot"]]
            )
        }
    } else {
        Plot <- Plot+ geom_beeswarm(aes(x = Conditions, y = Intensity), size = 2)
    }

    # # ##---- Add stats:
    if (pval == "t.test" | pval == "wilcox.test") {
        # One vs. One comparison: t-test
        if (!is.null(stat_comparison)) {
            Plot <- Plot+ stat_compare_means(comparisons = stat_comparison,
            label = "p.format", method = pval, hide.ns = TRUE,
            position = position_dodge(0.9), vjust = 0.25, show.legend = FALSE)
        } else {
        comparison <- unique(plotdata$Conditions)
        Plot <- Plot+ stat_compare_means(comparisons = comparison,
            label = "p.format", method = pval, hide.ns = TRUE,
            position = position_dodge(0.9), vjust = 0.25, show.legend = FALSE)

        }
        Plot <- Plot +labs(caption = paste("p.val using pairwise ", pval))
        } else {
        # All-vs-All comparisons table:
        conditions <- metadata_sample$Conditions
        denominator <- unique(metadata_sample$Conditions)
        numerator <- unique(metadata_sample$Conditions)
        comparisons <- combn(unique(conditions), 2) %>% as.matrix()

        # Prepare Stat results using dma STAT helper functions
        if (pval == "aov") {
            STAT_C1vC2 <- mpv_aov(data = data.frame("Intensity" = plotdata[, -c(2:3)]),
                            metadata_info = c(Conditions = "Conditions", Numerator = unique(metadata_sample$Conditions), Denominator = unique(metadata_sample$Conditions)),
                            metadata_sample = metadata_sample,
                            log2fc_table = NULL)
        } else if (pval == "kruskal.test") {
            STAT_C1vC2 <- mpv_kruskal(data = data.frame("Intensity" = plotdata[, -c(2:3)]),
                                            metadata_info = c(Conditions = "Conditions", Numerator = unique(metadata_sample$Conditions), Denominator = unique(metadata_sample$Conditions)),
                                            metadata_sample = metadata_sample,
                                            log2fc_table = NULL)
        }

        # Prepare df to add stats to plot
        df <-
            data.frame(
                comparisons = names(STAT_C1vC2),
                stringsAsFactors = FALSE
            ) %>%
            separate(
                comparisons,
                into = c("group1", "group2"),
                sep = "_vs_",
                remove = FALSE
            ) %>%
            unite(
                comparisons_rev,
                c("group2", "group1"),
                sep = "_vs_",
                remove = FALSE
            )
        df$p.adj <- round(map_dbl(STAT_C1vC2, ~ .x$p.adj), 5)

        # Add the 'res' column by repeating 'position' to match the number of rows
        position <- c(max(dataMeans$Intensity + 2*dataMeans$sd),
                        max(dataMeans$Intensity + 2*dataMeans$sd)+0.04* max(dataMeans$Intensity + 2*dataMeans$sd),
                        max(dataMeans$Intensity + 2*dataMeans$sd)+0.08* max(dataMeans$Intensity + 2*dataMeans$sd))

        df %<>%
        mutate(y.position = rep(position, length.out = n()))

        # select stats based on comparison_table
        if (!is.null(stat_comparison)) {
            # Generate the comparisons
            df_select <- data.frame()
            for (comp in stat_comparison) {
            entry <-
                paste0(
                    plot_conditions[comp[1]],
                    "_vs_",
                    plot_conditions[comp[2]]
                )
            df_select %<>% rbind(data.frame(entry))
            }

            df_merge <-
                merge(
                    df_select,
                    df,
                    by.x = "entry",
                    by.y = "comparisons",
                    all.x = TRUE
                ) %>%
            column_to_rownames("entry")

            # in case the reverse comparisons are needed
            if (all(is.na(df_merge))) {
                df_merge <-
                    merge(
                        df_select,
                        df,
                        by.x = "entry",
                        by.y = "comparisons_rev",
                        all.x = TRUE
                    ) %>%
                column_to_rownames("entry")
            }
        } else {
            df_merge <- df[, -2] %>%
            column_to_rownames("comparisons")
            }


        # add stats to plot
        if (plot_type == "Bar") {
            # http://rpkgs.datanovia.com/ggpubr/reference/stat_pvalue_manual.html
            Plot <- Plot +stat_pvalue_manual(df_merge, hide.ns = FALSE, size = 3, tip.length = 0.01, step.increase = 0.05)
        } else {
            # http://rpkgs.datanovia.com/ggpubr/reference/stat_pvalue_manual.html
            Plot <- Plot +stat_pvalue_manual(df_merge, hide.ns = FALSE, size = 3, tip.length = 0.01, step.increase = 0.01)
        }
        Plot <- Plot +labs(caption = paste("p.adj using ", pval, "and", padj))
    }

    Plot <- Plot + theme+ labs(title = plot_name,
        subtitle = i)  # ggtitle(paste(i))
    Plot <- Plot + theme(legend.position = "right", plot.title = element_text(size = 12, face = "bold"), axis.text.x = element_text(angle = 90, hjust = 1))+ xlab(xlab)+ ylab(ylab)

    ## Store the plot in the 'plots' list
    PlotList[[i]] <- Plot

    # Make plot into nice format:
    Plot_Sized <-  plot_grob_superplot(input_plot = Plot, metadata_info = metadata_info, metadata_sample = metadata_sample, plot_name = plot_name, subtitle = i, plot_type = plot_type)
    plot_height <- convertUnit(Plot_Sized$height, 'cm', valueOnly = TRUE)
    plot_width <- convertUnit(Plot_Sized$width, 'cm', valueOnly = TRUE)
    Plot_Sized %<>%
        {ggplot() + annotation_custom(.)} %>%
        add(theme(panel.background = element_rect(fill = "transparent")))

    ##########################################################################
    ## --------------- save ----------------- # #
    cleaned_i <-
        gsub(
            "[[:space:],/\\\\:*?\"<> |]",
            "-",
            i
        )  # removes empty spaces and replaces /,\ with -
    PlotList_adaptedGrid[[cleaned_i]] <- Plot_Sized

    SaveList <- list()
    SaveList[[cleaned_i]] <- Plot_Sized
    # ----- Save
    save_res(inputlist_df = NULL,
        inputlist_plot = SaveList,
        save_table = NULL,
        save_plot = save_plot,
        path = folder,
        file_name = paste(plot_type, "Plots_", plot_name, sep = ""),
        core = FALSE,
        print_plot = print_plot,
        plot_height = plot_height,
        plot_width = plot_width,
        plot_unit = "cm"
    )
    }
    return(invisible(list("Plot" = PlotList, "Plot_Sized" = PlotList_adaptedGrid)))
}

