# annoFuse v0.91.0 (released 2022-09-26)

## New features
* Add custom type annotation that checks for required columns and maps input columns to output columns using config file
* Fusion standardization function has the ability to accept custom type input that checks for few required columns. Config file has been made optionally and can be use to map input columns with output columns
* Custom type was built on top of existing starfusion and arriba fusion framework. This version will continue to support both of these caller types with same implementation

# annoFuse v0.90.0 (released 2020-09-18)

## New features
* `shinyFuse`, a shiny app for interactive exploration of fusions from _annoFuse_ prioritization
* `reportFuse`, a reproducible R Mardown notebook for fusion analysis using _annoFuse_
* Technical validation vignette using TCGA fusion calls
* Add functionality to visualize fusion breakpoints by fusion or by sample, in addition to viewing breakpoints for a specific fusion for all samples in a cohort
* Add breakpoint plots for 5'-fused gene, 3'-fused gene, in addition to plotting both genes
* Update all functions to _snake_case_
* Add breakpoint location (genic, intergenic, or intragenic) to standardFusioncalls
* Add SpanningDelta to standardFusioncalls
* Add reciprocal status for kinase fusions
* Update OpenPBTA data from `release-v14-20200203` to `release-v16-20200320`
* Add col_type reading function for STAR-Fusion and Arriba
* Add function examples
* Add project-specific fusion filtering example
* Make kinase domain as default for checkDomainStatus 

## Bug fixes

* Remove ConjoinG filter to capture known true positive reciprocal fusion _FIP1L1--PDGFRA_ discovered in OpenPBTA vignette
* Change default SpanningDelta filter from 10 to 100 based on TCGA validation
* Update frameshift domain retention to "Yes" when breakpoint occurs after kinase domain in Gene1A (5' gene)
* Update plot text sizing

## Other notes

* Added a `NEWS.md` file to track changes to the package.
