# annoFuse

Using annoFuse, users can filter out fusions known to be artifactual and retained high-quality fusion calls using support of at least one junction read and remove false calls if there is disproportionate spanning fragment support of more than 10 reads compared to the junction read count. 

For prioritization, users can capture known as well as putative driver fusions reported in TCGA, or fusions containing gene partners that are known oncogenes, tumor suppressor genes, or COSMIC genes. 

Finally, users can also determine recurrent fusions across the cohort and recurrently-fused genes within each histology. By providing a standardized filtering and annotation method from multiple callers (STAR-Fusion and Arriba) users are able to merge, filter and prioritize putative oncogenic fusions across the PBTA. 

## Getting Started

These instructions will get you a copy of the package up and running on your local machine. 

## Install package

```
devtools::install_github("d3b-center/annoFuse", dependencies = TRUE)
```

### Prerequisites for cohort level analysis

 - merge calls from each caller for you cohort and a column `annots` with additional annotation (eg from FusionAnnotator or caller specific annotation we have used FusionAnnotator in our vignettes)

 - reference folder  with a gene [genelistreference.txt](https://github.com/d3b-center/annoFuse/blob/master/inst/extdata/genelistreference.txt) and [fusionreference.txt](https://github.com/d3b-center/annoFuse/blob/master/inst/extdata/fusionreference.txt) inst/extdata has reference files we've used in our vignettes.
These files were created using publically available gene lists [kinases](http://kinase.com/human/kinome/tables/Kincat_Hsap.08.02.xls),[oncogenes](http://www.bushmanlab.org/assets/doc/allOnco_Feb2017.tsv),[tumor suppressors](https://bioinfo.uth.edu/TSGene/Human_TSGs.txt?csrt=5027697123997809089),curated transcription factors [@doi:10.1016/j.cell.2018.01.029],[COSMIC genes](https://cancer.sanger.ac.uk/census),[TCGA fusions](https://tumorfusions.org/PanCanFusV2/downloads/pancanfus.txt.gz) and _MYBL1_ [@doi:10.1073/pnas.1300252110], _SNCAIP_ [@doi:10.1038/nature11327], _FOXR2_ [@doi:10.1016/j.cell.2016.01.015], _TTYH1_ [@doi:10.1038/ng.2849], and _TERT_ [@doi:10.1038/ng.3438; @doi:10.1002/gcc.22110; @doi:10.1016/j.canlet.2014.11.057; @doi:10.1007/s11910-017-0722-5] were added to the oncogene list and _BCOR_ [@doi:10.1016/j.cell.2016.01.015] and _QKI_ [@doi:10.1038/ng.3500] were added to the tumor suppressor gene list based on pediatric cancer literature review.

 - expression matrix with GeneSymbol per row and samples as columns
 
### Prerequisites for single sample analysis

 - [STAR-Fusion star-fusion.fusion_predictions.tsv ](https://github.com/STAR-Fusion/STAR-Fusion/wiki#output-from-star-fusion)
 - [Arriba fusions.tsv](https://arriba.readthedocs.io/en/latest/output-files/)
 - RSEM genes.results.gz

## Overview of package

![](vignettes/Figure_1.png)


## Vignette

To browse vignettes

```
devtools::install_github("d3b-center/annoFuse", build_vignettes=TRUE, dependencies = TRUE)
browseVignettes("annoFuse")
```



## Authors
Krutika S. Gaonkar,Federico Marini, Komal S. Rathi, Jaclyn N. Taroni, Jo Lynne Rokita

## License

This project is licensed under the MIT License 
