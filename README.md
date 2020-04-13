# annoFuse
Using annoFuse, users can filter out fusions known to be artifactual and retained high-quality fusion calls using support of at least one junction read and remove false calls if there is disproportionate spanning fragment support of more than 10 reads compared to the junction read count. 
    For prioritization, users can capture known as well as putative driver fusions reported in TCGA, or fusions containing gene partners that are known oncogenes, tumor suppressor genes, or COSMIC genes. 
    Finally, users can also determine recurrent fusions across the cohort and recurrently-fused genes within each histology. By providing a standardized filtering and annotation method from multiple callers (STAR-Fusion and Arriba) users are able to merge, filter and prioritize putative oncogenic fusions across the PBTA. 

## Getting Started
These instructions will get you a copy of the project up and running on your local machine. 

## Install package
`devtools::install_github("d3b-center/annoFuse")`

### Prerequisites
 - merge calls from each caller for you cohort and add tumor_id to each file before merging to be able to differentiate between the calls and a column annots with additional annotation (eg from FusionAnnotator or caller specific annotation)
 - reference folder <link to box for example> with a gene genelistreference.txt and fusionreference.txt
 - expression matrix with GeneSymbol per row or should be able to collap to per GeneSymbol
 - normal expression matrix from GTEx or your own normal cohort if you require zscore annotation using a normal expression matrix
 
### Application
Please find detailed information regarding the functionalities [here](https://github.com/d3b-center/annoFuse/wiki)

** Single Sample run **

```
# Load required packages
suppressPackageStartupMessages(library("annoFuse"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("reshape2"))

# Run annoFuse for Single sample with default expression filter and FusionAnnotator red flag artifact filter 
standardFusioncalls<-annoFuseSingleSample(fusionfileArriba = "~/path/to/fusionfileArriba/Sample1_arriba_annotated.tsv",fusionfileStarFusion = "~/path/to/fusionfileStarFusion/Sample1_starfusion.tsv",expressionFile = "~/path/to/expressionFile/Sample1.rsem.genes.results.gz",tumorID = "Sample1")

# Add domain level information for fusion
bioMartDataPfam<-readRDS(system.file("extdata","pfamDataBioMart.RDS", package="annoFuse"))
annDomain<-annoFuse::getPfamDomain(standardFusioncalls  = res,bioMartDataPfam = bioMartDataPfam,keepPartialAnno = TRUE)

```


## Authors
Krutika S. Gaonkar, Komal S. Rathi, Jaclyn N. Taroni, Jo Lynne Rokita

## License

This project is licensed under the MIT License 

## Built With
R version 3.5.1 (2018-07-02) -- "Feather Spray"

## Version 
0.1.8

## Vignette
To browse vignettes use:
`browseVignettes("annoFuse")` 
