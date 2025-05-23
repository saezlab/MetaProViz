% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Processing.R
\name{processing}
\alias{processing}
\title{Modularised Normalization: 80\%-filtering rule, total-ion count normalization, missing value imputation and Outlier Detection: HotellingT2.}
\usage{
processing(
  data,
  metadata_sample,
  metadata_info,
  featurefilt = "Modified",
  cutoff_featurefilt = 0.8,
  tic = TRUE,
  mvi = TRUE,
  mvi_percentage = 50,
  hotellins_confidence = 0.99,
  core = FALSE,
  save_plot = "svg",
  save_table = "csv",
  print_plot = TRUE,
  path = NULL
)
}
\arguments{
\item{data}{DF which contains unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected.}

\item{metadata_sample}{DF which contains information about the samples, which will be combined with the input data based on the unique sample identifiers used as rownames.}

\item{metadata_info}{Named vector containing the information about the names of the experimental parameters. c(Conditions="ColumnName_Plot_SettingsFile", Biological_Replicates="ColumnName_Plot_SettingsFile"). Column "Conditions" with information about the sample conditions (e.g. "N" and "T" or "Normal" and "Tumor"), can be used for feature filtering and colour coding in the PCA. Column "BiologicalReplicates" including numerical values. For core = TRUE a core_norm_factor = "Columnname_Input_SettingsFile" and core_media = "Columnname_Input_SettingsFile", have to also be added. Column core_norm_factor is used for normalization and core_media is used to specify the name of the media controls in the Conditions.}

\item{featurefilt}{\emph{Optional: }If NULL, no feature filtering is performed. If set to "Standard" then it applies the 80\%-filtering rule (Bijlsma S. et al., 2006) on the metabolite features on the whole dataset. If is set to "Modified",filtering is done based on the different conditions, thus a column named "Conditions" must be provided in the Input_SettingsFile input file including the individual conditions you want to apply the filtering to (Yang, J et al., 2015). \strong{Default = "Standard"}}

\item{cutoff_featurefilt}{\emph{Optional: } percentage of feature filtering. \strong{Default = 0.8}}

\item{tic}{\emph{Optional: } If TRUE, total Ion Count normalization is performed. \strong{Default = TRUE}}

\item{mvi}{\emph{Optional: } If TRUE, Missing Value Imputation (mvi) based on half minimum is performed \strong{Default = TRUE}}

\item{mvi_percentage}{\emph{Optional: } percentage 0-100 of imputed value based on the minimum value. \strong{Default = 50}}

\item{hotellins_confidence}{\emph{Optional: } Defines the Confidence of Outlier identification in HotellingT2 test. Must be numeric.\strong{Default = 0.99}}

\item{core}{\emph{Optional: } If TRUE, a consumption-release experiment has been performed and the core value will be calculated. Please consider providing a Normalisation factor column called "core_norm_factor" in your "Input_SettingsFile" DF, where the column "Conditions" matches. The normalisation factor must be a numerical value obtained from growth rate that has been obtained from a growth curve or growth factor that was obtained by the ratio of cell count/protein quantification at the start point to cell count/protein quantification at the end point.. Additionally control media samples have to be available in the "Input" DF and defined as "core_media" samples in the "Conditions" column in the "Input_SettingsFile" DF. \strong{Default = FALSE}}

\item{save_plot}{\emph{Optional: } Select the file type of output plots. Options are svg, png, pdf. If set to NULL, plots are not saved. \strong{Default = svg}}

\item{save_table}{\emph{Optional: } Select the file type of output table. Options are "csv", "xlsx", "txt". If set to NULL, plots are not saved. \strong{Default = "csv"}}

\item{print_plot}{\emph{Optional: } If TRUE prints an overview of resulting plots. \strong{Default = TRUE}}

\item{path}{\emph{Optional:} Path to the folder the results should be saved at. \strong{default: NULL}}
}
\value{
List with two elements: DF (including all output tables generated) and Plot (including all plots generated)
}
\description{
Modularised Normalization: 80\%-filtering rule, total-ion count normalization, missing value imputation and Outlier Detection: HotellingT2.
}
\examples{
Intra <- intracell_raw \%>\%tibble::column_to_rownames("Code")
ResI <- MetaProViz::processing(data=Intra[-c(49:58) ,-c(1:3)],
                                 metadata_sample=Intra[-c(49:58) , c(1:3)],
                                 metadata_info = c(Conditions = "Conditions", Biological_Replicates = "Biological_Replicates"))

Media <- medium_raw \%>\%tibble::column_to_rownames("Code")
ResM <- MetaProViz::processing(data = Media[-c(40:45) ,-c(1:3)],
                                  metadata_sample = Media[-c(40:45) ,c(1:3)] ,
                                  metadata_info = c(Conditions = "Conditions", Biological_Replicates = "Biological_Replicates", core_norm_factor = "GrowthFactor", core_media = "blank"),
                                  core=TRUE)

}
\keyword{80}
\keyword{Count}
\keyword{HotellingT2,}
\keyword{Imputation,}
\keyword{Ion}
\keyword{Missing}
\keyword{PCA,}
\keyword{Value}
\keyword{charts}
\keyword{control}
\keyword{filtering}
\keyword{multivariate}
\keyword{normalization,}
\keyword{percent}
\keyword{quality}
\keyword{rule,}
\keyword{total}
