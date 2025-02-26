% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HelperOptions.R
\name{metaproviz_save_config}
\alias{metaproviz_save_config}
\title{Save the current package configuration}
\usage{
metaproviz_save_config(path = NULL, title = "default", local = FALSE)
}
\arguments{
\item{path}{Path to the config file. Directories and the file will be
created if don't exist.}

\item{title}{Save the config under this title. One config file might
contain multiple configurations, each identified by a title.}

\item{local}{Save into a config file in the current directory instead of
a user level config file. When loading, the config in the current
directory has priority over the user level config.}
}
\value{
Returns \code{NULL}.
}
\description{
Save the current package configuration
}
\examples{
\dontrun{
# after this, all downloads will default to commercial licenses
# i.e. the resources that allow only academic use will be excluded:
options(metaproviz.console_loglevel = 'trace')
metaproviz_save_config()
}

}
