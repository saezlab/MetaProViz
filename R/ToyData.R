## ---------------------------
##
## Script name: HelperFunctions
##
## Purpose of script: General helper functions to check function input and save results
##
## Author: Christina Schmidt
##
## Date Created: 2023-06-14
##
## Copyright (c) Christina Schmidt
## Email:
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

################################################################################################
### ### ### Example data ### ### ###
################################################################################################

#' Access built-in example data
#'
#' @param dataset Character: name of a built-in dataset:
#'     \itemize{
#'         \item{\code{"IntraCells_Raw"}: }
#'         \item{\code{"IntraCells_DMA"}: }
#'         \item{\code{"CultureMedia_Raw"}: }
#'         \item{\code{"Cells_Metadata"}: }
#'         \item{\code{"Tissue_Norm"}: }
#'         \item{\code{"Tissue_Metadata"}: }
#'         \item{\code{"Tissue_DMA"}: }
#'         \item{\code{"Tissue_DMA_Old"}: }
#'         \item{\code{"Tissue_DMA_Young"}: }
#'         \item{\code{"Tissue_TvN_Proteomics"}: }
#'         \item{\code{"Tissue_TvN_RNAseq"}: }
#'         \item{\code{"EquivalentFeatures"}: }
#'         \item{\code{"BiocratesFeatureTable"}: }
#'     }
#'
#' @return A data frame containing the toy data.
#'
#' @description Import and process .csv file to create toy data DF.
#'
#' @examples
#' Intra <- MetaProViz::ToyData("IntraCells_Raw")
#'
#' @importFrom readr read_csv cols
#' @importFrom magrittr %>% extract2
#' @importFrom tibble column_to_rownames
#' @importFrom logger log_trace
#'
#' @export
#'
ToyData <- function(dataset) {
  ## ------------ Create log file ----------- ##
  MetaProViz_Init()

  #Available datasets:
  datasets <- list(
    IntraCells_Raw = "MS55_RawPeakData.csv.gz",
    IntraCells_DMA = "MS55_DMA_786M1A_vs_HK2.csv.gz",
    CultureMedia_Raw = "MS51_RawPeakData.csv.gz",
    Cells_Metadata = "MappingTable_SelectPathways.csv.gz",
    Tissue_Norm = "Hakimi_ccRCC-Tissue_Data.csv.gz",
    Tissue_Metadata = "Hakimi_ccRCC-Tissue_FeatureMetaData.csv.gz",
    Tissue_DMA = "Hakimi_ccRCC-Tissue_DMA_TvsN.csv.gz",
    Tissue_DMA_Old ="Hakimi_ccRCC-Tissue_DMA_TvsN-Old.csv.gz",
    Tissue_DMA_Young ="Hakimi_ccRCC-Tissue_DMA_TvsN-Young.csv.gz",
    Tissue_TvN_Proteomics ="ccRCC-Tissue_TvN_Proteomics.csv.gz",
    Tissue_TvN_RNAseq = "ccRCC-Tissue_TvN_RNAseq.csv.gz",
    AlaninePathways = "AlaninePathways.csv.gz",
    EquivalentFeatures = "EquivalentFeatureTable.csv.gz",
    BiocratesFeatureTable = "BiocratesFeatureTable.csv.gz"
  )
  newnames <- list(
    IntraCells_Raw = "intracell_raw",
    IntraCells_DMA = "intracell_dma",
    CultureMedia_Raw = "medium_raw",
    Cells_Metadata = "cellular_meta",
    Tissue_Norm = "tissue_norm",
    Tissue_Metadata = "tissue_meta",
    Tissue_DMA = "tissue_dma",
    Tissue_DMA_Old ="tissue_dma_old",
    Tissue_DMA_Young ="tissue_dma_young",
    Tissue_TvN_Proteomics ="tissue_tvn_proteomics",
    Tissue_TvN_RNAseq = "tissue_tvn_rnaseq",
    AlaninePathways = "alanine_pathways",
    EquivalentFeatures = "equivalent_features",
    BiocratesFeatureTable = "biocrates_freatures"
  )
  for (key in names(datasets)) {
    newname <- newnames[[key]]
    path <- sprintf('inst/extdata/%s', datasets[key])
    rda_path <- sprintf('data/%s.rda', newname)
    the_data <- read_csv(path, col_types = cols())
    assign(newname, the_data, envir = .GlobalEnv)
    print(newname)
    print(is.symbol(sym(newname)))
    save(newname, file = rda_path)
  }

  rncols <- c("Code", "Metabolite")

  #Load dataset:
  if (!dataset %in% names(datasets)) {
    message <- sprintf("No such dataset: `%s`. Available datasets: %s", dataset, paste(names(datasets), collapse = ", "))
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }

  datasets %>%
  magrittr::extract2(dataset) %>%
    system.file("extdata", ., package = "MetaProViz") %>%
    readr::read_csv(col_types = readr::cols()) %>%
    {`if`(
    (rncol <- names(.) %>% intersect(rncols)) %>% length,
    tibble::column_to_rownames(., rncol),
    .
    )}
}


################################################################################################
### ### ### Information about the Example data ### ### ###
################################################################################################

####################################################
#' intracell_raw
#'
#' Metabolomics workbench project PR001418, study ST002224 where we exported integrated raw
#' peak values of intracellular metabolomics of HK2 and ccRCC cell lines 786-O, 786-M1A and 786-M2A.
#'
#' @format A data frame with multiple rows and columns:
#' \describe{
#'   \item{Conditions}{Character vector indicating cell line identity}
#'   \item{Analytical_Replicate}{Integer replicate number for analytical replicates}
#'   \item{Biological_Replicate}{Integer replicate number for biological replicates}
#'   \item{...}{Numeric columns for each measured metabolite (raw peak values)}
#' }
#'
#' @source Sciacovelli & Dugourd et. al., Dynamic partitioning of branched-chain amino acids-derived
#' nitrogen supports renal cancer progression, Nature Communications 2022, DOI:10.1038/s41467-022-35036-4.
"intracell_raw"

####################################################
#' intracell_dma
#'
#' Metabolomics workbench project PR001418, study ST002224 where we performed differential metabolite
#' analysis comparing intracellular metabolomics of 786-M1A versus HK2 cells.
#'
#' @format Columns include Log2FC, stats, metabolite identifiers, metabolite pathways and normalised
#' metabolite values used as input with row names being metabolitetrivial names.
#'
#' @source  Sciacovelli & Dugourd et. al., Dynamic partitioning of branched-chain amino acids-derived
#' nitrogen supports renal cancer progression , Nature Communications 2022, \doi{10.1038/s41467-022-35036-4}
"intracell_dma"


####################################################
#' medium_raw
#'
#' Metabolomics workbench project PR001418, study ST002226 where we exported integrated raw
#' peak values of intracellular metabolomics of HK2 and cccRCC cell lines 786-O, 786-M1A,
#' 786-M2A, OS-RC-2, OS-LM1 and RFX-631.
#'
#' @format Columns include Conditions (=Cell lines, blanks),Biological_Replicate, GrowthFactor and
#' a numeric column for each measured metabolite (raw data)
#'
#' @source  Sciacovelli & Dugourd et. al., Dynamic partitioning of branched-chain amino acids-derived
#' nitrogen supports renal cancer progression , Nature Communications 2022, \doi{10.1038/s41467-022-35036-4}
"medium_raw"

####################################################
#' cellular_meta
#'
#' Metabolomics workbench project PR001418, study ST002226 and ST002224 measured metabolites were assigned
#' HMDB and KEGG IDs as well as one main metabolic pathway.
#'
#' @format Columns include Metabolite IDs (HMDB, KEGG), Main metabolic pathway with row names being
#' metabolite trivial names.
#'
#' @source  Sciacovelli & Dugourd et. al., Dynamic partitioning of branched-chain amino acids-derived
#' nitrogen supports renal cancer progression , Nature Communications 2022, \doi{10.1038/s41467-022-35036-4}
"cellular_meta"

####################################################
#' tissue_norm
#'
#' This is median normalised data from the supplementary table 2 of Hakimi et al with metabolomic
#' profiling on 138 matched clear cell renal cell carcinoma (ccRCC)/normal tissue pairs.
#'
#' @format Columns include patient metadata (e.g. age, gender, sage, etc.) a numeric column for each
#' measured metabolite (normalised data)
#'
#' @source Hakimi et. al, An integrated metabolic atlas of clear cell renal cell carcinoma, Cancer Cell 2016,
#' \doi{10.1016/j.ccell.2015.12.004}
"tissue_norm"

####################################################
#' Tissue_Metadata
#'
#' In Hakimi et. al. metabolites were assigned to metabolite IDs, pathways, platform, mass and other
#'  fetaure metainformation.
#'
#' @format Columns include Metabolite IDs (HMDB, KEGG, etc), platform, mass, metabolic pathway with
#' row names being metabolite trivial names.
#'
#' @source Hakimi et. al, An integrated metabolic atlas of clear cell renal cell carcinoma, Cancer Cell 2016,
#' \doi{10.1016/j.ccell.2015.12.004}
"tissue_meta"

####################################################
#' tissue_dma
#'
#' We performed differential metabolite analysis comparing ccRCC tissue versus adjacent normal tissue using
#' median normalised data from the supplementary table 2 of Hakimi et. al.(="Tissue_Norm").
#'
#' @format Columns include Log2FC, stats, metabolite identifiers, metabolite pathways and normalised
#' metabolite values used as input with row names being metabolite trivial names.
#'
#' @source Hakimi et. al, An integrated metabolic atlas of clear cell renal cell carcinoma, Cancer Cell 2016,
#' \doi{10.1016/j.ccell.2015.12.004}
"tissue_dma"

####################################################
#' tissue_dma_old
#'
#' We performed differential metabolite analysis comparing ccRCC tissue versus adjacent normal tissue of
#' the patient's subset of old patient's (age > 58 years) using median normalised data from the supplementary
#' table 2 of Hakimi et. al.(="Tissue_Norm").
#'
#' @format Columns include Log2FC, stats, metabolite identifiers, metabolite pathways and normalised
#' metabolite values used as input with row names being metabolite trivial names.
#'
#' @source Hakimi et. al, An integrated metabolic atlas of clear cell renal cell carcinoma, Cancer Cell 2016,
#' \doi{10.1016/j.ccell.2015.12.004}
"tissue_dma_old"

####################################################
#' tissue_dma_young
#'
#' We performed differential metabolite analysis comparing ccRCC tissue versus adjacent normal tissue of
#' the patient's subset of young patient's (age <42 years) using median normalised data from the supplementary
#' table 2 of Hakimi et. al.(="Tissue_Norm").
#'
#' @format Columns include Log2FC, stats, metabolite identifiers, metabolite pathways and normalised
#' metabolite values used as input with row names being metabolite trivial names.
#'
#' @source Hakimi et. al, An integrated metabolic atlas of clear cell renal cell carcinoma, Cancer Cell 2016,
#' \doi{10.1016/j.ccell.2015.12.004}
"tissue_dma_young"

####################################################
#' tissue_tvn_proteomics
#'
#' The processed proteomics data was downloaded from the supplementary table 3 of Mora & Schmidt et. al., which
#' used the study from Clark et. al. under Proteomics data Commons PDC000127.
#'
#' @format Columns include Log2FC, stats, gene name and SiRCle cluster information that summarises genes based
#' on their regulation
#'
#' @source Mora & Schmidt, SiRCle (Signature Regulatory Clustering) model integration reveals mechanisms of phenotype
#' regulation in renal cancer, Genome Medicine 2024, \doi{10.1186/s13073-024-01415-3} Clark et. al, Integrated
#' proteogenomic characterization of clear cell renal cell carcinoma, Cell 2019, \doi{10.1016/j.cell.2019.10.007}
"tissue_tvn_proteomics"

####################################################
#' tissue_tvn_rnaseq
#'
#' The processed transcriptomics data was downloaded from the supplementary table 3 of Mora & Schmidt et. al., which
#' used the study from Clark et. al. under Proteomics data Commons PDC000127.
#'
#' @format Columns include Log2FC, stats, gene name and SiRCle cluster information that summarises genes based
#' on their regulation
#'
#' @source Mora & Schmidt, SiRCle (Signature Regulatory Clustering) model integration reveals mechanisms of phenotype
#' regulation in renal cancer, Genome Medicine 2024, \doi{10.1186/s13073-024-01415-3} Clark et. al, Integrated
#' proteogenomic characterization of clear cell renal cell carcinoma, Cell 2019, \doi{10.1016/j.cell.2019.10.007}
"tissue_tvn_rnaseq"

####################################################
#' equivalent_features
#'
#' Manually curated list of aminoacids and aminoacid-related metabolites with corresponding metabolite
#' identifiers (HMDB, KEGG, etc.) irrespective of chirality.
#'
#' @format Columns include metabolite trivial names, metabolite IDs (HMDB, KEGG, etc.), metabolite structural information (=INCHI).
#'
#' @source Schmidt et al, MetaProViz: METabolomics pre-PRocessing, functiOnal analysis and VIZualisation version 2.1.7, GitHub 2025.
"equivalent_features"

####################################################
#' biocrates_features
#'
#' Biocrates kit feature information of the "MxP Quant 500 XL kit" that covers more than 1,000 metabolites with biochemical
#' class information and the exported different metabolite IDs (HMDB, KEGG, etc.).
#'
#' @format Columns include metabolite trivial name, metabolite class, metabolite IDs (HMDB, KEGG, etc.), metabolite structural
#' information (INCHI, Key, etc.).
#'
#' @source Biocrates MxP® Quant 500 XL kit, https://biocrates.com/mxp-quant-500-xl/
"biocrates_features"

####################################################
#' mca_2cond
#'
#' Manually curated table defining the flow of information of the two condition biological regulatory clusters Regulatory labels
#'from the different grouping methods.
#'
#' @format Columns include Intra, core, core_Detection including state entries (e.g. up, down, etc.) and the Regulator Clustering
#' columns (RG1-RG3)
#'
#' @source Schmidt et al, MetaProViz: METabolomics pre-PRocessing, functiOnal analysis and VIZualisation version 2.1.7, GitHub 2025.
"mca_2cond"

####################################################
#' mca_core
#'
#' Manually curated table defining the flow of information of the Conusuption-Release and Intracellular metabolomics biological
#' regulatory clusters Regulatory labels from the different grouping methods.
#'
#' @format Columns include Intra, core, core_Detection including state entries (e.g. up, down, etc.) and the Regulator Clustering
#' columns (RG1-RG3)
#'
#' @source Schmidt et al, MetaProViz: METabolomics pre-PRocessing, functiOnal analysis and VIZualisation version 2.1.7, GitHub 2025.
"mca_core"

####################################################
#' alanine_pathways
#'
#' Manually curated table for the amino acid alanine toshowcase pathways (wiki, reactome, etc.) and alanine IDs (chebi, hmdb, etc.)
#' included in those pathways
#'
#' @format Columns include pathway_name, pathwayId and pathwaySource as well as inputID  and commonName.
#'
#' @source Schmidt et al, MetaProViz: METabolomics pre-PRocessing, functiOnal analysis and VIZualisation version 2.1.7, GitHub 2025.
"alanine_pathways"
