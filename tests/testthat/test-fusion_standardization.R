context("Test fusion standardization")

fusionfileArriba <- read.delim(
  system.file("extdata", "arriba_example.tsv", package = "annoFuseData"), check.names = FALSE)

standardFusioncalls <-
  annoFuse::fusion_standardization(
    fusion_calls = fusionfileArriba,
    caller = "ARRIBA",
    tumorID = "BS_W97QQYKQ"
  )

test_that("Standardizing arriba calls", {
  expect_equal(colnames(standardFusioncalls), 
               c("LeftBreakpoint", "RightBreakpoint", "FusionName", 
                 "Sample", "Caller", "Fusion_Type", "JunctionReadCount", 
                 "SpanningFragCount", "Confidence", "annots", 
                 "GeneA", "Gene1A", "Gene2A", "GeneB", "Gene1B", "Gene2B", 
                 "BreakpointLocation", "SpanningDelta"))
  expect_equal(nrow(standardFusioncalls), 25)
  expect_equal(unique(standardFusioncalls$Sample), "BS_W97QQYKQ")
})


fusionfileStarfusion <- read.delim(
  system.file("extdata", "starfusion_example.tsv", package = "annoFuseData"), check.names = FALSE)

standardFusioncalls <-
  annoFuse::fusion_standardization(
    fusion_calls = fusionfileStarfusion,
    caller = "STARFUSION",
    tumorID = "BS_W97QQYKQ"
  )

test_that("Standardizing starfusion calls (with --examine_coding_effect) calls", {
  expect_equal(colnames(standardFusioncalls), 
               c("LeftBreakpoint", "RightBreakpoint", "FusionName", 
                 "Sample", "Caller", "Fusion_Type", "JunctionReadCount", "SpanningFragCount", 
                 "Confidence", "annots", 
                 "GeneA", "Gene1A", "Gene2A", "GeneB", "Gene1B", "Gene2B", 
                 "BreakpointLocation", "SpanningDelta"))
  expect_equal(nrow(standardFusioncalls), 1)
  expect_equal(unique(standardFusioncalls$Sample), "BS_W97QQYKQ")
})


fusionfileStarfusion <- read.delim(
  system.file("extdata", "starfusion_example.tsv", package = "annoFuseData"), check.names = FALSE)
fusionfileStarfusion <- fusionfileStarfusion[, 1:15]
standardFusioncalls <-
  annoFuse::fusion_standardization(
    fusion_calls = fusionfileStarfusion,
    caller = "STARFUSION",
    tumorID = "BS_W97QQYKQ"
  )

test_that("Standardizing starfusion calls (without --examine_coding_effect) calls", {
  expect_equal(colnames(standardFusioncalls), 
               c("LeftBreakpoint", "RightBreakpoint", "FusionName", 
                 "Sample", "Caller", "Fusion_Type", "JunctionReadCount", "SpanningFragCount", 
                 "Confidence", "annots", 
                 "GeneA", "Gene1A", "Gene2A", "GeneB", "Gene1B", "Gene2B", 
                 "BreakpointLocation", "SpanningDelta"))
  expect_equal(nrow(standardFusioncalls), 1)
  expect_equal(unique(standardFusioncalls$Sample), "BS_W97QQYKQ")
})

customtypefile <- read.delim(
  system.file("extdata", "custom_type_example.tsv", package = "annoFuseData"), check.names = FALSE)
configcustomtypepath <- system.file("extdata", "config", package = "annoFuseData")
CustomFusioncalls <-
  annoFuse::fusion_standardization(
    customtypefile,
    "CUSTOM",
    "random",
    input_json_file=configcustomtypepath
  )

test_that("Standardizing Custom type calls with config file", {
  expect_equal(colnames(CustomFusioncalls), 
               c("Sample", "FusionName", "Gene1A", "Gene1B", "Gene2A", "Gene2B", 
                 "Fusion_Type", "annots"))
  expect_equal(nrow(CustomFusioncalls), 283)
})