% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GetPriorKnoweldge.R
\name{TranslateID}
\alias{TranslateID}
\title{Translate IDs}
\usage{
TranslateID(
  Input_DataFrame,
  SettingsInfo = list(IdColumn = "MetaboliteID", FromFormat = c("kegg"), ToFormat =
    c("pubchem", "chebi", "hmdb"), Method = "GetAll", GroupingVariable = "term")
)
}
\arguments{
\item{Input_DataFrame}{Dataframe with two columns for source (=term) and Target (=gene), e.g. Hallmarks.}

\item{SettingsInfo}{\emph{Optional: }  Column name of Target in Input_GeneSet. \strong{Default = list(IdColumn="MetaboliteID", FromFormat=c("kegg"), ToFormat=c("pubchem","chebi","hmdb"), Method="GetAll", GroupingVariable="term")}}

\item{SaveAs_Table}{\emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt". \strong{Default = "csv"}}

\item{FolderPath}{{Optional:} String which is added to the resulting folder name \strong{Default = NULL}}
}
\value{
3 data frames: 1) Original data and the new column of translated ids. 2) Mapping summary from Original ID to Translated. 3) Mapping summary from Translated to Original.
}
\description{
Translate IDs to and from KEGG, PubChem, Chebi.
}
