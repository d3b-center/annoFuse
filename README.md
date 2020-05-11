# annoFuse
Using annoFuse, users can filter out fusions known to be artifactual and retained high-quality fusion calls using support of at least one junction read and remove false calls if there is disproportionate spanning fragment support of more than 10 reads compared to the junction read count. 
    For prioritization, users can capture known as well as putative driver fusions reported in TCGA, or fusions containing gene partners that are known oncogenes, tumor suppressor genes, or COSMIC genes. 
    Finally, users can also determine recurrent fusions across the cohort and recurrently-fused genes within each histology. By providing a standardized filtering and annotation method from multiple callers (STAR-Fusion and Arriba) users are able to merge, filter and prioritize putative oncogenic fusions across the PBTA. 

## Getting Started
These instructions will get you a copy of the package up and running on your local machine. 

## Install package
`devtools::install_github("d3b-center/annoFuse")`

### Prerequisites for cohort level analysis
 - merge calls from each caller for you cohort and add tumor_id to each file before merging to be able to differentiate between the calls and a column annots with additional annotation (eg from FusionAnnotator or caller specific annotation)
 - reference folder  with a gene [genelistreference.txt](https://github.com/d3b-center/annoFuse/blob/master/inst/extdata/genelistreference.txt) and [fusionreference.txt](https://github.com/d3b-center/annoFuse/blob/master/inst/extdata/fusionreference.txt)
 - expression matrix with GeneSymbol per row and samples as columns
 
### Prerequisites for single sample analysis
 - [STAR-Fusion star-fusion.fusion_predictions.tsv ](https://github.com/STAR-Fusion/STAR-Fusion/wiki#output-from-star-fusion)
 - [Arriba fusions.tsv](https://arriba.readthedocs.io/en/latest/output-files/)
 - RSEM genes.results.gz

## Vignette
To browse vignette to see example of ** Single Sample run ** 

`devtools::install_github("d3b-center/annoFuse", build_vignettes=TRUE)`

`browseVignettes("annoFuse")`

Please find detailed information regarding the functionalities [here](https://github.com/d3b-center/annoFuse/wiki)

## Authors
Krutika S. Gaonkar, Komal S. Rathi, Jaclyn N. Taroni, Jo Lynne Rokita

## License

This project is licensed under the MIT License 

## Built With
R version 3.5.1 (2018-07-02) -- "Feather Spray"

## Version 
0.2

