% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VizHeatmap.R
\name{VizHeatmap}
\alias{VizHeatmap}
\title{This script allows you to perform different data visualizations using the results of the MetaProViz analysis}
\usage{
VizHeatmap(
  InputData,
  SettingsInfo = NULL,
  SettingsFile_Sample = NULL,
  SettingsFile_Metab = NULL,
  PlotName = "",
  Scale = "row",
  SaveAs_Plot = "svg",
  Enforce_FeatureNames = FALSE,
  Enforce_SampleNames = FALSE,
  PrintPlot = TRUE,
  FolderPath = NULL
)
}
\arguments{
\item{InputData}{DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Includes experimental design and outlier column.}

\item{SettingsInfo}{\emph{Optional: } NULL or Named vector  where you can include vectors or lists for annotation c(individual_Metab= "ColumnName_SettingsFile_Metab",individual_Sample= "ColumnName_SettingsFile_Sample", color_Metab="ColumnName_SettingsFile_Metab", color_Sample= list("ColumnName_SettingsFile_Sample", "ColumnName_SettingsFile_Sample",...)).\strong{Default = NULL}}

\item{SettingsFile_Sample}{DF which contains information about the samples, which will be combined with your input data based on the unique sample identifiers. and other columns with required PlotSettingInfo.\strong{Default = NULL}}

\item{SettingsFile_Metab}{\emph{Optional: } DF with column "Metabolite" including the Metabolite names (needs to match Metabolite names of Input_data) and other columns with required PlotSettingInfo. \strong{Default = NULL}}

\item{PlotName}{\emph{Optional: } String which is added to the output files of the plot}

\item{Scale}{\emph{Optional: } String with the information for Scale row, column or none. \strong{Default = row}}

\item{SaveAs_Plot}{\emph{Optional: } Select the file type of output plots. Options are svg, pdf, png or NULL. \strong{Default = "svg"}}

\item{Enforce_FeatureNames}{\emph{Optional: } If there are more than 100 features no rownames will be shown, which is due to readability. You can Enforce this by setting this parameter to TRUE. \strong{Default = FALSE}}

\item{Enforce_SampleNames}{\emph{Optional: } If there are more than 50 sampless no colnames will be shown, which is due to readability. You can Enforce this by setting this parameter to TRUE. \strong{Default = FALSE}}

\item{Folder_Name}{{Optional:} String which is added to the resulting folder name \strong{default: NULL}}
}
\description{
This script allows you to perform different data visualizations using the results of the MetaProViz analysis
}
\keyword{Volcano}
\keyword{pathways}
\keyword{plot,}
