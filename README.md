# annoFuse

<!-- badges: start -->
[![R build status](https://github.com/d3b-center/annoFuse/workflows/R-CMD-check/badge.svg)](https://github.com/d3b-center/annoFuse/actions)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4036789.svg)](https://doi.org/10.5281/zenodo.4036789)

<!-- badges: end -->

Using annoFuse, users can filter out fusions known to be artifactual and retained high-quality fusion calls using support of at least one junction read and remove false calls if there is disproportionate spanning fragment support of more than 100 reads compared to the junction read count. 

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
The fusion reference is a compilation of the annotations listed in the table below.
 
 
 
 Annotation | File | Source  
------ | ---------- | --------- 
| pfamID                        | http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/pfamDesc.txt.gz     | UCSC pfamID Description database |
| Domain Location               | http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ucscGenePfam.txt.gz | UCSC pfamID Description database |
| TCGA fusions                  | https://tumorfusions.org/PanCanFusV2/downloads/pancanfus.txt.gz             | TumorFusions: an integrative   resource for cancer-associated transcript fusions PMID: 29099951  |
| Oncogenes                     | http://www.bushmanlab.org/assets/doc/allOnco_Feb2017.tsv                    | www.bushmanlab.org |
| Tumor suppressor genes (TSGs) | https://bioinfo.uth.edu/TSGene/Human_TSGs.txt?csrt=5027697123997809089      | Tumor Suppressor Gene Database   2.0 PMIDs: 23066107, 26590405 |
| Kinases                       | http://kinase.com/human/kinome/tables/Kincat_Hsap.08.02.xls |      The protein kinase complement of the human genome PMID: 12471243 |
| COSMIC genes                  | https://cancer.sanger.ac.uk/census | Catalogue of Somatic Mutations   in Cancer |
| Pediatric-specific oncogenes  | _MYBL1, SNCAIP, FOXR2, TTYH1, TERT_ | doi:10.1073/pnas.1300252110,   doi:10.1038/nature11327, doi:10.1016/j.cell.2016.01.015, doi:10.1038/ng.2849,   doi:10.1038/ng.3438, doi:10.1002/gcc.22110, doi:10.1016/j.canlet.2014.11.057,   doi:10.1007/s11910-017-0722-5 |
| Pediatric-specific TSGs | _BCOR_, _QKI_  | doi:10.1016/j.cell.2016.01.015, doi:10.1038/ng.3500 |

 
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
Krutika S. Gaonkar, Federico Marini, Komal S. Rathi, Jaclyn N. Taroni, Jo Lynne Rokita

## How to cite

Krutika S. Gaonkar, Federico Marini, Komal S. Rathi, Payal Jain, Yuankun Zhu, Nicholas A. Chimicles, Miguel A. Brown, Ammar S. Naqvi, Bo Zhang, Phillip B. Storm, John M. Maris, Pichai Raman, Adam C. Resnick, Konstantin Strauch, Jaclyn N. Taroni & Jo Lynne Rokita (2020). annoFuse: an R Package to annotate, prioritize, and interactively explore putative oncogenic RNA fusions. BMC Bioinformatics, 21(1), 577. https://doi.org/10.1186/s12859-020-03922-7

## License

This project is licensed under the MIT License 
