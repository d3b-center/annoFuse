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
The fusion reference is made up of TCGA fusions (Hu et al. 2018), gene reference file was created from publically available reference files : oncogene (http://www.bushmanlab.org/assets/doc/allOnco_Feb2017.tsv), tumor suppressor (Zhao et al. 2013)-(Zhao et al. 2016), kinase from pfam location database (http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ucscGenePfam.txt.gz) and id description from pfamID Description database (http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/pfamDesc.txt.gz), as well as from Human Kinome database (Manning et al. 2002), transcription factor (Lambert et al. 2018), and  whether it has been reported in the Catalogue of Somatic Mutations in Cancer (COSMIC) Cancer Gene Census (Sondka et al. 2018). In addition,  we added the following genes specific to pediatric cancer from literature review; _MYBL1_ (Ramkissoon et al. 2013), _SNCAIP_ (Northcott et al. 2012), _FOXR2_ (Sturm et al. 2016), _TTYH1_ (Kleinman et al. 2014), and _TERT_ (Valentijn et al. 2015);(Lindner et al. 2015);(Karlsson et al. 2015);(Karsy et al. 2017)  were added to the oncogene list and _BCOR_ (Sturm et al. 2016) and _QKI_ (Bandopadhayay et al. 2016).

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
