context("Test fusion QC filtering")

fusionfileArriba <- read.delim(system.file("extdata", "arriba_example.tsv", package = "annoFuseData"), check.names = FALSE, stringsAsFactors = FALSE)

formattedArriba <-
  annoFuse::fusion_standardization(
    fusion_calls = fusionfileArriba,
    caller = "ARRIBA",
    tumorID = "BS_W97QQYKQ"
  )


fusionfileStarfusion <- read.delim(system.file("extdata", "starfusion_example.tsv", package = "annoFuseData"), check.names = FALSE, stringsAsFactors = FALSE)

formattedStarFusion <-
  annoFuse::fusion_standardization(
    fusion_calls = fusionfileStarfusion,
    caller = "STARFUSION",
    tumorID = "BS_W97QQYKQ"
  )

standardFusioncalls <- as.data.frame(rbind(formattedStarFusion, formattedArriba))

fusionQCFiltered <- fusion_filtering_QC(
  standardFusioncalls = standardFusioncalls,
  readingFrameFilter = "in-frame|frameshift|other",
  artifactFilter = "GTEx_Recurrent|DGD_PARALOGS|Normal|BodyMap|ConjoinG",
  junctionReadCountFilter = 1,
  spanningFragCountFilter = 100,
  readthroughFilter = TRUE
)


test_that("Fusion filtering for standardized arriba and starfusion calls with default", {
  expect_equal(colnames(fusionQCFiltered), c("LeftBreakpoint", "RightBreakpoint", "FusionName", "Sample", "Caller", "Fusion_Type", "JunctionReadCount", "SpanningFragCount", "Confidence", "annots", "GeneA", "Gene1A", "Gene2A", "GeneB", "Gene1B", "Gene2B", "BreakpointLocation", "SpanningDelta"))
  expect_equal(nrow(fusionQCFiltered), 16)
})


fusionQCFiltered <- fusion_filtering_QC(
  standardFusioncalls = standardFusioncalls,
  readingFrameFilter = "in-frame|frameshift|other",
  artifactFilter = "GTEx_Recurrent|DGD_PARALOGS|Normal|BodyMap|ConjoinG",
  junctionReadCountFilter = 1,
  spanningFragCountFilter = 100,
  readthroughFilter = FALSE
)

test_that("Fusion filtering for standardized arriba and starfusion calls default but keeping readthroughs", {
  expect_equal(colnames(fusionQCFiltered), c("LeftBreakpoint", "RightBreakpoint", "FusionName", "Sample", "Caller", "Fusion_Type", "JunctionReadCount", "SpanningFragCount", "Confidence", "annots", "GeneA", "Gene1A", "Gene2A", "GeneB", "Gene1B", "Gene2B", "BreakpointLocation", "SpanningDelta"))
  expect_equal(nrow(fusionQCFiltered), 26)
})

fusionQCFiltered <- fusion_filtering_QC(
  standardFusioncalls = standardFusioncalls,
  readingFrameFilter = "in-frame|frameshift|other",
  artifactFilter = "GTEx_Recurrent|DGD_PARALOGS|Normal|BodyMap|ConjoinG",
  junctionReadCountFilter = 10,
  spanningFragCountFilter = 100,
  readthroughFilter = FALSE
)

test_that("Fusion filtering for standardized arriba and starfusion calls default but keeping readthroughs and filter junction read support >= 10", {
  expect_equal(colnames(fusionQCFiltered), c("LeftBreakpoint", "RightBreakpoint", "FusionName", "Sample", "Caller", "Fusion_Type", "JunctionReadCount", "SpanningFragCount", "Confidence", "annots", "GeneA", "Gene1A", "Gene2A", "GeneB", "Gene1B", "Gene2B", "BreakpointLocation", "SpanningDelta"))
  expect_equal(nrow(fusionQCFiltered), 5)
})

fusionQCFiltered <- fusion_filtering_QC(
  standardFusioncalls = standardFusioncalls,
  readingFrameFilter = "in-frame",
  artifactFilter = "GTEx_Recurrent|DGD_PARALOGS|Normal|BodyMap|ConjoinG",
  junctionReadCountFilter = 1,
  spanningFragCountFilter = 100,
  readthroughFilter = FALSE
)

test_that("Fusion filtering for standardized arriba and starfusion calls default but keeping readthroughs and inframe only ", {
  expect_equal(colnames(fusionQCFiltered), c("LeftBreakpoint", "RightBreakpoint", "FusionName", "Sample", "Caller", "Fusion_Type", "JunctionReadCount", "SpanningFragCount", "Confidence", "annots", "GeneA", "Gene1A", "Gene2A", "GeneB", "Gene1B", "Gene2B", "BreakpointLocation", "SpanningDelta"))
  expect_equal(nrow(fusionQCFiltered), 7)
})
