# Tutorials

All tutorials use example data that ship with **MetaProViz**, so you can
run every step yourself. If you are new to MetaProViz, start with
*Getting started*. The other tutorials go into more detail on the
different data types and on metabolite prior knowledge.

## Get started

![Overview of MetaProViz modules](figures/thumbs/quick-start.png)

##### [Getting started](https://saezlab.github.io/MetaProViz/articles/quick-start.md)

A short tour through the main steps on intracellular cell line data:
pre-processing, differential analysis, enrichment analysis, metabolite
clustering and plots.

Bioconductor vignette

## Analysis workflows

![Heatmap of metabolites coloured by
pathway](figures/thumbs/standard-metabolomics.png)

##### [Standard Metabolomics](https://saezlab.github.io/MetaProViz/articles/pkgdown/standard-metabolomics.md)

Intracellular metabolomics of kidney cancer cell lines: quality control,
pre-processing, differential analysis, ORA, metabolite clustering and
all visualisations.

Website tutorial

![Volcano plot of consumption-release
data](figures/thumbs/core-metabolomics.png)

##### [CoRe Metabolomics](https://saezlab.github.io/MetaProViz/articles/pkgdown/core-metabolomics.md)

Consumption-release data from cell culture media: blank and growth
normalisation, differential analysis, clustering together with
intracellular data and metabolite-receptor sets.

Website tutorial

![Variance of principal components explained by patient
metadata](figures/thumbs/sample-metadata.png)

##### [Sample Metadata Analysis](https://saezlab.github.io/MetaProViz/articles/pkgdown/sample-metadata.md)

Tumour and normal tissue of kidney cancer patients: find the metabolites
that separate patient groups, compare patient subsets and improve
metabolite IDs before enrichment analysis.

Website tutorial

![Network of metabolites and the KEGG pathways they
share](figures/thumbs/pk-networks.png)

##### [Prior Knowledge Networks](https://saezlab.github.io/MetaProViz/articles/pkgdown/pk-networks.md)

Investigate differential analysis results with networks: which
transporters the changed metabolites share, which metabolites drive
enriched pathways and in which cancers they were reported before.

Website tutorial

## Prior knowledge and metabolite IDs

![Graph of KEGG pathways coloured by
coverage](figures/thumbs/prior-knowledge.png)

##### [Prior Knowledge - Access & Integration](https://saezlab.github.io/MetaProViz/articles/pkgdown/prior-knowledge.md)

Load metabolite sets from MetSigDB, link them to your measured data,
translate metabolite IDs and check how well your data cover each
pathway.

Website tutorial

![Scheme of the ID processing
workflow](figures/thumbs/id-processing-workflow.png)

##### [ID Processing Workflow](https://saezlab.github.io/MetaProViz/articles/pkgdown/id-processing-workflow.md)

Check the metabolite IDs of your features for consistency and expand
them across HMDB, KEGG, ChEBI and PubChem with
[`id_processing()`](https://saezlab.github.io/MetaProViz/reference/id_processing.md).

Website tutorial

![Overview of the MetSigDB resources](figures/thumbs/metsigdb.png)

##### [MetSigDB](https://saezlab.github.io/MetaProViz/articles/pkgdown/metsigdb.md)

The resources in MetSigDB compared: their size, overlap and redundancy,
and how their terms cluster.

Website tutorial
