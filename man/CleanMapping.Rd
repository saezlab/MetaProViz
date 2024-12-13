% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RefactorPriorKnoweldge.R
\name{CleanMapping}
\alias{CleanMapping}
\title{Reduce mapping ambiguities by using detected IDs}
\usage{
CleanMapping(
  InputData,
  PriorKnowledge,
  DetectedID = "HMDB",
  GroupingVariable = "term",
  From = "MetaboliteID",
  To = "hmdb",
  SaveAs_Table = "csv",
  FolderPath = NULL
)
}
\arguments{
\item{InputData}{Dataframe with at least one column with the target (e.g. metabolite)}

\item{PriorKnowledge}{Dataframe with at least one column with the "original target metabolite ID", "source" (e.g. term) and the "translated target metabolite ID".}

\item{SettingsInfo}{}
}
\value{
Prior Knowledge usable for enrichment analysis
}
\description{
Reduce mapping ambiguities by using detected IDs
}
\examples{
DetectedIDs <-  MetaProViz::ToyData(Data="Cells_MetaData")\%>\% rownames_to_column("Metabolite") \%>\%select("Metabolite", "HMDB")
PathwayFile <- MetaProViz::TranslateID(InputData= MetaProViz::LoadKEGG(), SettingsInfo = c(InputID="MetaboliteID", GroupingVariable="term"), From = c("kegg"), To = c("hmdb"))[["TranslatedDF"]]
Res <- MetaProViz::CleanMapping(InputData= DetectedIDs, PriorKnowledge= PathwayFile, DetectedID="HMDB", GroupingVariable="term", From="MetaboliteID", To="hmdb")

}
