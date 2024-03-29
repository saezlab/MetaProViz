% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OverRepresentationAnalysis.R
\name{MC_ORA}
\alias{MC_ORA}
\title{\subsection{---------------------------}{
}}
\usage{
MC_ORA(
  Input_data,
  MetabCluster_lab = "RG2_Significant",
  RemoveBackground = TRUE,
  PathwayFile,
  PathwayName = "",
  minGSSize = 10,
  maxGSSize = 1000,
  Save_as_Results = "csv",
  Folder_Name = NULL
)
}
\arguments{
\item{Input_data}{Input DF with column "Metabolite", which needs to match the identifier type in column ""Metabolite" in the PathwayFile.}

\item{MetabCluster_lab}{\emph{Optional: } Provide column name for the metabolite cluster labels in the Input_data.\strong{default: "RG2_Significant"}}

\item{PathwayFile}{DF that must include column "term" with the pathway name, column "Metabolite" with the Metabolite name or ID and column "Description" with pathway description that will be depicted on the plots.}

\item{PathwayName}{\emph{Optional: } Name of the pathway list used \strong{default: ""}}

\item{minGSSize}{\emph{Optional: } minimum group size in ORA \strong{default: 10}}

\item{maxGSSize}{\emph{Optional: } maximum group size in ORA \strong{default: 1000}}

\item{Save_as_Results}{\emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt" \strong{default: "csv"}}

\item{Folder_Name}{{Optional:} String which is added to the resulting folder name \strong(Default = NULL)}

\item{RemoveBackground\emph{Optional:}{} If TRUE, column "BG_Method" needs to be in Input_data, which includes TRUE/FALSE for each metabolite to fall into background based on the chosen Background method for e.g. MCA_2Cond are removed from the universe. \strong{default: TRUE}}

\item{Save_as_Plot}{\emph{Optional: } If Save_as_Plot=NULL no plots will be saved. Otherwise, file types for the figures are: "svg", "pdf", "png" \strong{default: "pdf"}}
}
\value{
Saves results as individual .csv files.
}
\description{
This script allows you
MC_ORA
}
\details{
Uses enricher to run ORA on each of the metabolite cluster from any of the MCA functions using a pathway list
}
