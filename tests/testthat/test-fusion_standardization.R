context("Test fusion standardization")

fusionfileArriba <- read.delim(system.file("extdata", "arriba_example.tsv", package = "annoFuseData"), check.names = FALSE)

standardFusioncalls <-
  annoFuse::fusion_standardization(
    fusion_calls = fusionfileArriba,
    caller = "ARRIBA",
    tumorID = "BS_W97QQYKQ"
  )

test_that("Standardizing arriba calls", {
  expect_equal(colnames(standardFusioncalls), c("LeftBreakpoint", "RightBreakpoint", "FusionName", "Sample", "Caller", "Fusion_Type", "JunctionReadCount", "SpanningFragCount", "Confidence", "annots", "GeneA", "Gene1A", "Gene2A", "GeneB", "Gene1B", "Gene2B", "BreakpointLocation", "SpanningDelta"))
  expect_equal(nrow(standardFusioncalls), 25)
  expect_equal(unique(standardFusioncalls$Sample), "BS_W97QQYKQ")
})


fusionfileStarfusion <- read.delim(system.file("extdata", "starfusion_example.tsv", package = "annoFuseData"), check.names = FALSE)

standardFusioncalls <-
  annoFuse::fusion_standardization(
    fusion_calls = fusionfileStarfusion,
    caller = "STARFUSION",
    tumorID = "BS_W97QQYKQ"
  )

test_that("Standardizing starfusion calls (with --examine_coding_effect) calls", {
  expect_equal(colnames(standardFusioncalls), c("LeftBreakpoint", "RightBreakpoint", "FusionName", "Sample", "Caller", "Fusion_Type", "JunctionReadCount", "SpanningFragCount", "Confidence", "annots", "GeneA", "Gene1A", "Gene2A", "GeneB", "Gene1B", "Gene2B", "BreakpointLocation", "SpanningDelta"))
  expect_equal(nrow(standardFusioncalls), 1)
  expect_equal(unique(standardFusioncalls$Sample), "BS_W97QQYKQ")
})


fusionfileStarfusion <- read.delim(system.file("extdata", "starfusion_example.tsv", package = "annoFuseData"), check.names = FALSE)

fusionfileStarfusion <- fusionfileStarfusion[, 1:15]
standardFusioncalls <-
  annoFuse::fusion_standardization(
    fusion_calls = fusionfileStarfusion,
    caller = "STARFUSION",
    tumorID = "BS_W97QQYKQ"
  )


test_that("Standardizing starfusion calls (without --examine_coding_effect) calls", {
  expect_equal(colnames(standardFusioncalls), c("LeftBreakpoint", "RightBreakpoint", "FusionName", "Sample", "Caller", "Fusion_Type", "JunctionReadCount", "SpanningFragCount", "Confidence", "annots", "GeneA", "Gene1A", "Gene2A", "GeneB", "Gene1B", "Gene2B", "BreakpointLocation", "SpanningDelta"))
  expect_equal(nrow(standardFusioncalls), 1)
  expect_equal(unique(standardFusioncalls$Sample), "BS_W97QQYKQ")
})
