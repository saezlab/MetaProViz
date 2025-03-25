## ---------------------------
##
## Script name: Visualization  PCA plot
##
## Purpose of script: Data Visualisation of the MetaProViz analysis to aid biological interpretation
##
## Author: Dimitrios Prymidis and Christina Schmidt
##
## Date Created: 2022-10-28
##
## Copyright (c) Dimitrios Prymidis and Christina Schmidt
## Email:
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------
#'
#' This script allows you to perform PCA plot visualization using the results of the MetaProViz analysis


#################################
### ### ### PCA Plots ### ### ###
#################################

#' PCA plot visualization
#'
#' @param Input_data DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected. includes experimental design and outlier column.
#' @param SettingsInfo \emph{Optional: } NULL or Named vector including at least one of those three information : c(color="ColumnName_Plot_SettingsFile", shape= "ColumnName_Plot_SettingsFile"). \strong{Default = NULL}
#' @param SettingsFile_Sample \emph{Optional: } DF which contains information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames. Column "Conditions" with information about the sample conditions (e.g. "N" and "T" or "Normal" and "Tumor"), can be used for feature filtering and colour coding in the PCA. Column "AnalyticalReplicate" including numerical values, defines technical repetitions of measurements, which will be summarised. Column "BiologicalReplicates" including numerical values. Please use the following names: "Conditions", "Biological_Replicates", "Analytical_Replicates".\strong{Default = NULL}
#' @param ColorPalette \emph{Optional: } Provide customiced color-palette in vector format. For continuous scale use e.g. scale_color_gradient(low = "#88CCEE", high = "red") and for discrete scale c("#88CCEE",  "#DDCC77","#661100",  "#332288")\strong{Default = NULL}
#' @param ColorScale \emph{Optional: } Either "continuous" or "discrete" colour scale. For numeric or integer you can choose either, for character you have to choose discrete. \strong{Default = NULL}
#' @param ShapePalette \emph{Optional: } Provide customiced shape-palette in vector format. \strong{Default = NULL}
#' @param ShowLoadings \emph{Optional: } TRUE or FALSE for whether PCA loadings are also plotted on the PCA (biplot) \strong{Default = FALSE}
#' @param Scaling \emph{Optional: } TRUE or FALSE for whether a data scaling is used \strong{Default = TRUE}
#' @param PCx \emph{Optional: } Numeric value of the PC that should be plotted on the x-axis \strong{Default = 1}
#' @param PCy \emph{Optional: } Numeric value of the PC that should be plotted on the y-axis \strong{Default = 2}
#' @param Theme \emph{Optional: } Selection of theme for plot, e.g. theme_grey(). You can check for complete themes here: https://ggplot2.tidyverse.org/reference/ggtheme.html. If default=NULL we use theme_classic(). \strong{Default = "discrete"}
#' @param PlotName \emph{Optional: } String which is added to the output files of the PCA \strong{Default = ""}
#' @param SaveAs_Plot \emph{Optional: } Select the file type of output plots. Options are svg, png, pdf or NULL. \strong{Default = svg}
#' @param PrintPlot \emph{Optional: } TRUE or FALSE, if TRUE Volcano plot is saved as an overview of the results. \strong{Default = TRUE}
#' @param FolderPath \emph{Optional:} Path to the folder the results should be saved at. \strong{default: NULL}
#'
#' @return List with two elements: Plot and Plot_Sized
#'
#' @examples
#' Intra <- ToyData("IntraCells_Raw")[,-c(1:3)]
#' Res <- VizPCA(Intra)
#'
#' @keywords PCA
#'
#' @importFrom ggplot2 ggplot theme element_rect autoplot scale_shape_manual geom_hline geom_vline ggtitle
#' @import ggfortify
#' @importFrom dplyr rename
#' @importFrom magrittr %>% %<>%
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom rlang !! :=
#' @importFrom logger log_info log_trace
#' @importFrom S4Vectors DataFrame
#'
#' @export
#'
VizPCA <- function(se, #InputData,
                   SettingsInfo = NULL,
                   #SettingsFile_Sample = NULL,
                   ColorPalette = NULL,
                   ColorScale = "discrete", ## EDIT: would list here the option and use match.arg
                   ShapePalette = NULL,
                   ShowLoadings = FALSE,
                   Scaling = TRUE,
                   PCx = 1,
                   PCy = 2,
                   Theme = NULL, ##theme_classic() ## EDIT: why not preset it, would simplify some things downstream?
                   PlotName = '',
                   SaveAs_Plot = "svg", ## EDIT: would list here the option and use match.arg
                   PrintPlot = TRUE,
                   FolderPath = NULL) {

    ###########################################################################
    ## ------------ Create log file ----------- ##
    MetaProViz_Init()

    logger::log_info("VizPCA: PCA plot visualization")
    ## ------------ Check Input files ----------- ##
    ## HelperFunction `CheckInput`
    CheckInput(se = se, ##InputData = InputData,
        ##SettingsFile_Sample = SettingsFile_Sample,
        ##SettingsFile_Metab = NULL,
        SettingsInfo = SettingsInfo,
        SaveAs_Plot = SaveAs_Plot,
        SaveAs_Table = NULL,
        CoRe = FALSE,
        PrintPlot = PrintPlot)

    # CheckInput` Specific
    if (!is.logical(ShowLoadings)) {
        message <- paste("The Show_Loadings value should be either TRUE if loadings are to be shown on the PCA plot or FALSE if not.")
        logger::log_trace(paste0("Error ", message))
        stop(message)
    }
    if (!is.logical(Scaling)) {
        message <- paste("The Scaling value should be either TRUE if data scaling is to be performed prior to the PCA or FALSE if not.")
        logger::log_trace(paste0("Error ", message))
        stop(message)
    }

    if (any(is.na(assay(se)))) {
        assay(se)[is.na(assay(se))] <- 0 #replace NA with 0
        message <- paste("NA values are included in InputData that were set to 0 prior to performing PCA.")
        logger::log_info(message)
        message(message)
    }

    ## ------------ Create Results output folder ----------- ##
    Folder <- NULL
    if (!is.null(SaveAs_Plot)) {
        Folder <- SavePath(FolderName = "PCAPlots", FolderPath = FolderPath)
    }
    logger::log_info("VizPCA results saved at ", Folder)

    ###########################################################################
    ## ----------- Set the plot parameters: ------------ ##
    ##--- Prepare colour and shape palette
    if (is.null(ColorPalette)) {
        if ((ColorScale == "discrete")) {
            safe_colorblind_palette <- c("#88CCEE",  "#DDCC77","#661100", ## EDIT: could this be defined outside of the function?
                "#332288", "#AA4499","#999933", "#44AA99", "#882215", "#6699CC", 
                "#117733", "#888888","#CC6677", "black", "gold1", "darkorchid4",
                "red", "orange", "blue")
        } else if(ColorScale == "continuous") {
            safe_colorblind_palette <- NULL
        }
    } else {
        safe_colorblind_palette <-ColorPalette
    }
    
    if (is.null(ShapePalette)) {
        safe_shape_palette <- c(15, 17, 16, 18, 6, 7, 8, 11, 12)
    } else {
        safe_shape_palette <-ShapePalette
    }

    logger::log_info(paste("VizPCA colour:", paste(safe_colorblind_palette, collapse = ", ")))
    logger::log_info(paste("VizPCA shape:", paste(safe_shape_palette, collapse = ", ")))

    ##--- Prepare the color scheme:
    if ("color" %in% names(SettingsInfo) & "shape" %in% names(SettingsInfo)) {
        if((SettingsInfo[["shape"]] == SettingsInfo[["color"]])){
            colData(se)[, "shape"] <- colData(se)[, SettingsInfo[["color"]]]
            se@colData <- colData(se) |>
                as.data.frame() |>
                dplyr::rename("color" = SettingsInfo[["color"]]) |>
                DataFrame()
        } else {
            se@colData <- colData(se) |>
                as.data.frame() |>
                dplyr::rename(
                    "color" = SettingsInfo[["color"]],
                    "shape" = SettingsInfo[["shape"]]) |>
                DataFrame()
        }
    } else if("color" %in% names(SettingsInfo) & !"shape" %in% names(SettingsInfo)) {
        if ("color" %in% names(SettingsInfo)) {
            se@colData <- colData(se) %>%
                as.data.frame() |>
                dplyr::rename("color" = SettingsInfo[["color"]]) |>
                DataFrame()
        }
        if ("shape" %in% names(SettingsInfo)) {
            se@colData <- colData(se) %>%
                as.data.frame() |>
                dplyr::rename("shape" = SettingsInfo[["shape"]]) |>
                DataFrame()
        }
    }

    ##--- Prepare Input Data:
    ##if (!is.null(SettingsFile_Sample)){
    InputPCA <- merge(colData(se), t(assay(se)), by = "row.names", all.y = TRUE) 
        ##InputPCA  <- merge(
       ##     x = tibble::rownames_to_column(as.data.frame(colData(se)), "UniqueID"), 
       ##     y = tibble::rownames_to_column(InputData, "UniqueID"), 
       ##     by = "UniqueID", all.y = TRUE) %>%
       ##     tibble::column_to_rownames("UniqueID")
   ## } else {
##        InputPCA  <- InputData
  ##  }

    ##--- Prepare the color and shape settings:
    if ("color" %in% names(colData(se))) {
        if (ColorScale == "discrete") {
            InputPCA$color <- as.factor(InputPCA$color)
            color_select <- safe_colorblind_palette[1:length(unique(InputPCA$color))]
        } else if (ColorScale=="continuous") {
            if (is.numeric(InputPCA$color) | is.integer(InputPCA$color)) {
                InputPCA$color <- as.numeric(InputPCA$color)
                color_select <- safe_colorblind_palette
            } else {
                InputPCA$color <- as.factor(InputPCA$color)
                ## Overwrite color pallette
                safe_colorblind_palette <- metaproviz_palette()
                ## color that will be used for distinct
                color_select <- safe_colorblind_palette[1:length(unique(InputPCA$color))]
                ## Overwrite color_scale
                ColorScale <- "discrete"
                logger::log_info("Warning: ColorScale=continuous, but is.numeric or is.integer is FALSE, hence colour scale is set to discrete.")
                warning("ColorScale=continuous, but is.numeric or is.integer is FALSE, hence colour scale is set to discrete.")
            }
        }
    }

    logger::log_info("VizPCA ColorScale: ", ColorScale)

    if ("shape" %in% names(colData(se))) {
        shape_select <- safe_shape_palette[1:length(unique(InputPCA$shape))]

        if (!is.character(InputPCA$shape)) {
            ## Convert the column to character
            InputPCA$shape <- as.character(InputPCA$shape)
        }
    }

    ##---  #assign column and legend name
    if ("color" %in% names(colData(se))) {
        InputPCA  <- InputPCA %>%
            dplyr::rename(!!SettingsInfo[["color"]] :="color")
        Param_Col <- SettingsInfo[["color"]]
    } else{
        color_select <- NULL
        Param_Col <- NULL
    }

    if ("shape" %in% names(colData(se))) {
        InputPCA  <- InputPCA %>%
            dplyr::rename(!!SettingsInfo[["shape"]] :="shape")
        Param_Sha <- SettingsInfo[["shape"]]
    } else {
        shape_select <-NULL
        Param_Sha <-NULL
    }

    ## ----------- Make the  plot based on the choosen parameters ------------ ##
    PlotList <- list() #Empty list to store all the plots
    PlotList_adaptedGrid <- list()

    ## Make the plot:
    PCA <- ggplot2::autoplot(
            object = stats::prcomp(as.matrix(InputPCA[, rownames(se)]), 
                scale. = as.logical(Scaling)), 
            data = InputPCA,
            x = PCx, y = PCy, colour = Param_Col, fill =  Param_Col, 
            shape = Param_Sha, size = 3, alpha = 0.8, label = TRUE,
            label.size = 2.5, label.repel = TRUE,
            loadings = as.logical(ShowLoadings), #draws Eigenvectors 
            loadings.label = as.logical(ShowLoadings), 
            loadings.label.vjust = 1.2, loadings.label.size=2.5,
            loadings.colour = "grey10", loadings.label.colour = "grey10") +
        ggplot2::scale_shape_manual(values = shape_select) +
        ggplot2::ggtitle(paste(PlotName)) +
        ggplot2::geom_hline(yintercept = 0,  color = "black", linewidth = 0.1) +
        ggplot2::geom_vline(xintercept = 0,  color = "black", linewidth = 0.1)

    if (ColorScale == "discrete") {
        PCA <- PCA + 
            ggplot2::scale_color_manual(values = color_select)
    } else if (ColorScale == "continuous" & is.null(ColorPalette)) {
        PCA <- PCA + 
            color_select
    }

    ## Add the theme
    if(!is.null(Theme)){
        PCA <- PCA + 
            Theme
    } else {
        PCA <- PCA + ggplot2::theme_classic()
    }

    ## Store the plot in the 'plots' list
    PlotList[["Plot"]] <- PCA

    ## Set the total heights and widths
    PCA %<>% 
        PlotGrob_PCA(SettingsInfo = SettingsInfo, PlotName = PlotName)
    PlotHeight <- grid::convertUnit(PCA$height, 'cm', valueOnly = TRUE)
    PlotWidth <- grid::convertUnit(PCA$width, 'cm', valueOnly = TRUE)
    PCA %<>%
        {ggplot2::ggplot() + annotation_custom(.)} %>%
        add(theme(panel.background = element_rect(fill = "transparent")))

    PlotList_adaptedGrid[["Plot_Sized"]] <- PCA

    ###########################################################################
    ##----- Save and Return
    #Here we make a list in which we will save the outputs:
    FileName <- PlotName %>% 
        {`if`(nchar(.), sprintf('PCA_%s', .), 'PCA')}

    suppressMessages(suppressWarnings(
        SaveRes(
            data = NULL,
            plot = PlotList_adaptedGrid,
            SaveAs_Table = NULL,
            SaveAs_Plot = SaveAs_Plot,
            FolderPath = Folder,
            FileName = FileName,
            CoRe = FALSE,
            PrintPlot = PrintPlot,
            PlotHeight = PlotHeight,
            PlotWidth = PlotWidth,
            PlotUnit = "cm")
    ))

    invisible(list(Plot = PlotList, Plot_Sized = PlotList_adaptedGrid))
}
