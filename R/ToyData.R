#!/usr/bin/env Rscript

#
#  This file is part of the `MetaProViz` R package
#
#  Copyright 2022-2025
#  Saez Lab, Heidelberg University
#
#  Authors: see the file `README.md`
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file `LICENSE` or copy at
#      https://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: https://saezlab.github.io/MetaProViz
#  Git repo: https://github.com/saezlab/MetaProViz
#


################################################################################################
### ### ### Built-in example data (docs) ### ### ###
################################################################################################

####################################################
#' intracell_raw
#'
#' Metabolomics workbench project PR001418, study ST002224 where we exported integrated raw
#' peak values of intracellular metabolomics of HK2 and ccRCC cell lines 786-O, 786-M1A and 786-M2A.
#'
#' @format A data frame with the following columns:
#'  -`Conditions`: Character vector indicating cell line identity
#'   - `Analytical_Replicate`: Integer replicate number for analytical replicates
#'  -`Biological_Replicate`: Integer replicate number for biological replicates
#'   - Additional numeric columns (183 in total) containing raw metabolite intensities
#'
#' @source Sciacovelli & Dugourd et. al., Dynamic partitioning of branched-chain amino acids-derived
#' nitrogen supports renal cancer progression, Nature Communications 2022, DOI:10.1038/s41467-022-35036-4.
#'
#' @examples
#' data(intracell_raw)
#' head(intracell_raw)
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
#'
#' @examples
#' data(intracell_dma)
#' head(intracell_dma)
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
#'
#' @examples
#' data(medium_raw)
#' head(medium_raw)
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
#'
#' @examples
#' data(cellular_meta)
#' head(cellular_meta)
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
#'
#' @examples
#' data(tissue_norm)
#' head(tissue_norm)
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
#' @examples
#' data(tissue_meta)
#' head(tissue_meta)
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
#'
#' @examples
#' data(tissue_dma)
#' head(tissue_dma)
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
#' @examples
#' data(tissue_dma_old)
#' head(tissue_dma_old)
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
#' @examples
#' data(tissue_dma_young)
#' head(tissue_dma_young)
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
#' @examples
#' data(tissue_tvn_proteomics)
#' head(tissue_tvn_proteomics)
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
#' @examples
#' data(tissue_tvn_rnaseq)
#' head(tissue_tvn_rnaseq)
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
#' @examples
#' data(equivalent_features)
#' head(equivalent_features)
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
#' @examples
#' data(biocrates_features)
#' head(biocrates_features)
"biocrates_features"

####################################################
#' mca_twocond_rules
#'
#' Manually curated table defining the flow of information of the two condition biological regulatory clusters Regulatory labels
#' from the different grouping methods.
#'
#' @format Columns include Intra, core, core_Detection including state entries (e.g. up, down, etc.) and the Regulator Clustering
#' columns (RG1-RG3)
#'
#' @source Schmidt et al, MetaProViz: METabolomics pre-PRocessing, functiOnal analysis and VIZualisation version 2.1.7, GitHub 2025.
#' @examples
#' data(mca_twocond_rules)
#' head(mca_twocond_rules)
"mca_twocond_rules"

####################################################
#' mca_core_rules
#'
#' Manually curated table defining the flow of information of the Conusuption-Release and Intracellular metabolomics biological
#' regulatory clusters Regulatory labels from the different grouping methods.
#'
#' @format Columns include Intra, core, core_Detection including state entries (e.g. up, down, etc.) and the Regulator Clustering
#' columns (RG1-RG3)
#'
#' @source Schmidt et al, MetaProViz: METabolomics pre-PRocessing, functiOnal analysis and VIZualisation version 2.1.7, GitHub 2025.
#' @examples
#' data(mca_core_rules)
#' head(mca_core_rules)
"mca_core_rules"

####################################################
#' alanine_pathways
#'
#' Manually curated table for the amino acid alanine toshowcase pathways (wiki, reactome, etc.) and alanine IDs (chebi, hmdb, etc.)
#' included in those pathways
#'
#' @format Columns include pathway_name, pathwayId and pathwaySource as well as inputID  and commonName.
#'
#' @source Schmidt et al, MetaProViz: METabolomics pre-PRocessing, functiOnal analysis and VIZualisation version 2.1.7, GitHub 2025.
#' @examples
#' data(alanine_pathways)
#' head(alanine_pathways)
"alanine_pathways"


#' hallmarks
#'
#' @format Data frame with columns: term, gene
#'
#' @source Omnipath
#' @examples
#' data(hallmarks)
#' head(hallmarks)
"hallmarks"

#' gaude_pathways
#'
#' @format Data frame with columns: term, gene
#'
#' @source Omnipath
#' @examples
#' data(gaude_pathways)
#' head(gaude_pathways)
"gaude_pathways"
