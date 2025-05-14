#' The MetaProViz package
#'
#' MetaProViz (Metabolomics Processing, functional analysis and Visualization),
#' a free open-source R-package that provides mechanistic hypotheses from
#' metabolomics data by integrating prior knowledge from literature with
#' metabolomics. MetaProViz offers an interactive framework consisting of five
#' modules: Processing, differential analysis, prior knwoledge access and
#' refactoring, functional analysis and visualization of both intracellular and
#' exometabolomics (=consumption-release (core) data).
#'
#' @author Christina Schmidt <\email{christina.schmidt@uni-heidelberg.de}> and
#' Denes Turei <\email{turei.denes@@gmail.com}> and Dimitrios Prymidis
#' and Macabe Daley and Julio Saez-Rodriguez and Christian Frezza
#'
#' @name MetaProViz
'_PACKAGE'

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
