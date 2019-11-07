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
 
 ## The package has the following 4 steps to filter and annotate fusion calls
 ### 1) Standardize calls from fusion callers to retain information regarding fused genes,breakpoints, reading frame information as well as annotation from FusionAnnotator. 
 
 Input : Merged fusion calls per caller with additional columns "annots" and "tumor_id"
 
 Output: Standardized fusion calls with following columns
 
   LeftBreakpoint : 5' gene breakpoint

   RightBreakpoint : 3' gene breakpoint

   FusionName : geneA--geneB

   Sample : tumor_id used by user in merged samples set

   Caller : eg StarFusion, Arriba etc

   Fusion_Type : reading frame information
    
   JunctionReadCount : junction supporting reads

   SpanningFragCount : fragments spanning the fusion 

   Confidence : Confidence provided from caller if not NA

   annots : Annotation provided by user; recommended FusionAnnotator

 ### 2) Filter standardized fusion calls to remove false positives with low read support, annotated as read-throughs, found in normal and gene homolog databases and if both fused genes are not expressed above the given threshold.
 
 Input : Standardized fusion calls from step1
 
 Output: Standardized fusion calls after filtering readthroughs, artifacts from annotation column and expression filter above threshold (default 1)
 
 ### 3) a. Annotate genes in standardized and filtered fusion calls with useful biological features of interest eg. Kinase, Tumor suppressor etc. 
 
 Input: Filtered standardized fusion calls from step2
 
 Ouptut: Standardized fusion calls with annotation per gene. Since callers like arriba also call intergenic fusions we have divided the fused genes as gene1A--gene1B geneic fusion between gene1A and gene2A; if fusion has intergenic 5' breakpoint then the fusion name would be gene1A/gene2A--gene1B and if 3' breakpoint is intergenic the fusion name would be gene1A--gene1B/gene2B 
 Additional columns are:
 
 Gene1A_anno : annotation per gene from reference gene list for gene1a
 
 Gene1B_anno : annotation per gene from reference gene list for gene1b
 
 Gene2A_anno : annotation per gene from reference gene list for gene2a
 
 Gene2B_anno : annotation per gene from reference gene list for gene2b
 
 Fusion_anno : annotation per gene from reference fusion list
 
 ###  b. (OPTIONAL) Annotation of zscore using a normal expression matrix
 
 Input : Filtered standardized fusion calls from step2
 
 Output : Annotation from zscore calculation using a normal expression matrix to identify if expression is differential compared to normal per gene in fusion 
 Additional columns are:
 
 note_expression_[Gene1A |Gene2A |Gene1B|Gene2B]: "differential expressed" if zscore is more than threshold or NA "
 
 zscore_[Gene1A |Gene2A |Gene1B|Gene2B] : zscore calculated from using normal expression matrix
 
 ### 4) Project specific filtering to capture recurrent fused genes or genes with specific biological feature from step 3.
 
 Input : Filtered fusion calls annotated with reference gene and/or normal zscore 
 
 Output : Project specific filtered fusion calls optional [summary plot like ](https://github.com/d3b-center/annoFuse/blob/master/vignettes/output.pdf)
 
## Authors
Krutika S. Gaonkar, Komal S. Rathi, Jaclyn N. Taroni, Jo Lynne Rokita

## License

This project is licensed under the MIT License 

## Built With
R version 3.5.1 (2018-07-02) -- "Feather Spray"

## Version 
0.1.0

## Vignette
To browse vignettes use:
`browseVignettes("annoFuse")` 
