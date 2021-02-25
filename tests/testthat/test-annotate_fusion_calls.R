context("Test fusion annotation")

fusionfileArriba <- read_arriba_calls(system.file("extdata", "arriba_example.tsv", package = "annoFuse"))

formattedArriba <- 
  annoFuse::fusion_standardization(fusion_calls = fusionfileArriba,
                                   caller= "ARRIBA",
                                   tumorID = "BS_W97QQYKQ")


fusionfileStarfusion <- read_starfusion_calls(system.file("extdata", "starfusion_example.tsv", package = "annoFuse"))

formattedStarFusion <- 
  annoFuse::fusion_standardization(fusion_calls = fusionfileStarfusion,
                                   caller= "STARFUSION",
                                   tumorID = "BS_W97QQYKQ")

standardFusioncalls <- as.data.frame(rbind(formattedStarFusion, formattedArriba))

geneListReferenceDataTab <- read.delim(system.file("extdata", "genelistreference.txt", package = "annoFuse"), stringsAsFactors = FALSE)

fusionReferenceDataTab <- read.delim(
  system.file("extdata", "fusionreference.txt", package = "annoFuse"), stringsAsFactors = FALSE)

filteredFusionAnnotated <- annotate_fusion_calls(
  standardFusioncalls = standardFusioncalls,
  geneListReferenceDataTab = geneListReferenceDataTab,
  fusionReferenceDataTab = fusionReferenceDataTab,
  checkReciprocal = TRUE)  

test_that("Standardizing arriba calls", {
  expect_equal(colnames(filteredFusionAnnotated), c("LeftBreakpoint","RightBreakpoint","FusionName","Sample","Caller","Fusion_Type","JunctionReadCount","SpanningFragCount","Confidence","annots","GeneA","Gene1A","Gene2A","GeneB","Gene1B","Gene2B","BreakpointLocation","SpanningDelta","reciprocal_exists","Gene1A_anno","Gene1B_anno","Gene2A_anno","Gene2B_anno","Fusion_anno" ))
  expect_equal(nrow(filteredFusionAnnotated), 26)
  expect_equal(unique(filteredFusionAnnotated$Sample), "BS_W97QQYKQ")
})