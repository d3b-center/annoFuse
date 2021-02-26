context("Test fusion driver")

fusionfileArriba <- read_arriba_calls(system.file("extdata", "arriba_example.tsv", package = "annoFuseData"))

sfc <- 
  annoFuse::fusion_standardization(fusion_calls = fusionfileArriba,
                                   caller= "ARRIBA",
                                   tumorID = "BS_W97QQYKQ")
geneListReferenceDataTab <- read.delim(system.file("extdata", "genelistreference.txt", package = "annoFuseData"), stringsAsFactors = FALSE)
fusionReferenceDataTab <- read.delim(system.file("extdata", "fusionreference.txt", package = "annoFuseData"), stringsAsFactors = FALSE)

bioMartDataPfam <- readRDS(system.file("extdata", "pfamDataBioMart.RDS", package = "annoFuseData"))
kinaseid<-unique(bioMartDataPfam$pfam_id[grep("kinase",bioMartDataPfam$NAME)] )
fusion_driver_df <- fusion_driver(sfc,
                                   annotated = FALSE,
                                   geneListReferenceDataTab = geneListReferenceDataTab,
                                   fusionReferenceDataTab = fusionReferenceDataTab,
                                   checkDomainStatus=TRUE,
                                   domainsToCheck=kinaseid)


test_that("Fusion driver annotation for standardized arriba calls", {
  expect_equal(colnames(fusion_driver_df), c("LeftBreakpoint","RightBreakpoint","FusionName","Sample","Caller","Fusion_Type","JunctionReadCount","SpanningFragCount","Confidence","annots","GeneA","Gene1A","Gene2A","GeneB","Gene1B","Gene2B","BreakpointLocation","SpanningDelta","DomainRetainedGene1A","DomainRetainedGene1B","reciprocal_exists","Gene1A_anno","Gene1B_anno","Gene2A_anno","Gene2B_anno","Fusion_anno"))
  expect_equal(nrow(fusion_driver_df), 5)
  expect_equal(fusion_driver_df$DomainRetainedGene1B[1], "Yes")
})
