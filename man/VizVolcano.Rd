% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VizVolcano.R
\name{VizVolcano}
\alias{VizVolcano}
\title{This script allows you to perform different data visualizations using the results of the MetaProViz analysis}
\usage{
VizVolcano(
  Plot_Settings = "Standard",
  Plot_SettingsInfo = NULL,
  Plot_SettingsFile = NULL,
  Input_data,
  y = "p.adj",
  x = "Log2FC",
  AdditionalInput_data = NULL,
  OutputPlotName = "",
  Comparison_name = c(Input_data = "Cond1", AdditionalInput_data = "Cond2"),
  xlab = NULL,
  ylab = NULL,
  pCutoff = 0.05,
  FCcutoff = 0.5,
  color_palette = NULL,
  shape_palette = NULL,
  SelectLab = "",
  Connectors = FALSE,
  Subtitle = "",
  Theme = NULL,
  Save_as_Plot = "svg",
  Folder_Name = NULL
)
}
\arguments{
\item{Plot_Settings}{\emph{Optional: } Choose between "Standard" (Input_data), "Compare" (plot two comparisons together Input_data and Input_data2) or "PEA" (Pathway Enrichment Analysis) \strong{Default = "Standard"}}

\item{Plot_SettingsInfo}{\emph{Optional: } NULL or Named vector including at least one of those three information for Plot_Settings="Standard" or "Compare": c(color="ColumnName_Plot_SettingsFile", shape= "ColumnName_Plot_SettingsFile", individual="ColumnName_Plot_SettingsFile"). For Plot_Settings="PEA" a named vector with c(PEA_Pathway="ColumnNameAdditionalInput_data", PEA_score="ColumnNameAdditionalInput_data", PEA_stat= "ColumnNameAdditionalInput_data", individual="Plot_SettingsFile), optionally you can additionally include c(color="ColumnName_Plot_SettingsFile", shape= "ColumnName_Plot_SettingsFile").\strong{Default = NULL}}

\item{Plot_SettingsFile}{\emph{Optional: } DF with column "Metabolite" including the Metabolite names (needs to match Metabolite names of Input_data) and other columns with required PlotSettingInfo. \strong{Default = NULL}}

\item{Input_data}{DF with column "Metabolite" including the Metabolite names, Log2FC, pvalue/padjusted values. Can also include additional columns with metadata usable for Plot_Setting_Info.}

\item{y}{\emph{Optional: } Column name including the values that should be used for y-axis. Usually this would include the p.adjusted value. \strong{Default = "p.adj"}}

\item{x}{\emph{Optional: } Column name including the values that should be used for x-axis. Usually this would include the Log2FC value. \strong{Default = "Log2FC"}}

\item{AdditionalInput_data}{\emph{Optional: } DF to compare to main Input_data with the same column names x and y (Plot_Settings="Compare") or Pathway enrichment analysis results (Plot_Settings="PEA"). \strong{Default = NULL}}

\item{OutputPlotName}{\emph{Optional: } String which is added to the output files of the plot. \strong{Default = ""}}

\item{Comparison_name}{\emph{Optional: } Named vector including those information about the two datasets that are compared on the plots when choosing Plot_Settings= "Compare". \strong{Default = c(Input_data="Cond1", AdditionalInput_data= "Cond2")}}

\item{xlab}{\emph{Optional: } String to replace x-axis label in plot. \strong{Default = NULL}}

\item{ylab}{\emph{Optional: } String to replace y-axis label in plot. \strong{Default = NULL}}

\item{pCutoff}{\emph{Optional: } Number of the desired p value cutoff for assessing significance. \strong{Default = 0.05}}

\item{FCcutoff}{\emph{Optional: } Number of the desired log fold change cutoff for assessing significance. \strong{Default = 0.5}}

\item{color_palette}{\emph{Optional: } Provide customiced color-palette in vector format. \strong{Default = NULL}}

\item{shape_palette}{\emph{Optional: } Provide customiced shape-palette in vector format. \strong{Default = NULL}}

\item{SelectLab}{\emph{Optional: } If set to NULL, feature labels will be plotted randomly. If vector is provided, e.g. c("MetaboliteName1", "MetaboliteName2"), selected names will be plotted. If set to default "", no feature names will be plotted. \strong{Default = ""}}

\item{Connectors}{\emph{Optional: } TRUE or FALSE for whether Connectors from names to points are to be added to the plot. \strong{Default =  FALSE}}

\item{Subtitle}{\emph{Optional: } \strong{Default = ""}}

\item{Theme}{\emph{Optional: } Selection of theme for plot, e.g. theme_grey(). You can check for complete themes here: https://ggplot2.tidyverse.org/reference/ggtheme.html. \strong{Default = NULL}}

\item{Save_as_Plot}{\emph{Optional: } Select the file type of output plots. Options are svg, pdf, png or NULL. \strong{Default = "svg"}}

\item{Folder_Name}{{Optional:} String which is added to the resulting folder name \strong(Default = NULL)}
}
\description{
This script allows you to perform different data visualizations using the results of the MetaProViz analysis
}
\keyword{Volcano}
\keyword{pathways}
\keyword{plot,}
