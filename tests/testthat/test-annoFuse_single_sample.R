context("Test single sample annoFuse filtering")


standardFusioncalls <- annoFuse::annoFuse_single_sample(
  # Example files are provided in extdata, at-least 1 fusionfile is required along
  # with its rsem expression file
  fusionfileArriba = system.file("extdata", "arriba_example.tsv", package = "annoFuseData"),
  fusionfileStarFusion = system.file("extdata", "starfusion_example.tsv", package = "annoFuseData"),
  expressionFile = system.file("extdata", "example.rsem.genes.results.gz", package = "annoFuseData"),
  tumorID = "BS_W97QQYKQ",
  # multiple read flag values for filtering using FusionAnnotator values
  artifactFilter = "GTEx_Recurrent|DGD_PARALOGS|Normal|BodyMap|ConjoinG",
  # keep all in-frame , frameshift and other types of Fusion_Type
  readingFrameFilter = "in-frame|frameshift|other",
  # keep all fusions with atleast 1 junction read support
  junctionReadCountFilter = 1,
  # keep only fusions where spanningFragCount-junctionReadCountFilter less than equal to 10
  spanningFragCountFilter = 10,
  # keep read throughs
  readthroughFilter = FALSE
)


test_that("Standard fusion single sample output from non-empty arriba and starfusion fusion calls", {
  expect_equal(colnames(standardFusioncalls), c("LeftBreakpoint", "RightBreakpoint", "FusionName", "Sample", "Caller", "Fusion_Type", "JunctionReadCount", "SpanningFragCount", "Confidence", "annots", "Gene1A", "Gene2A", "Gene1B", "Gene2B", "BreakpointLocation", "SpanningDelta", "reciprocal_exists", "Gene1A_anno", "Gene1B_anno", "Gene2A_anno", "Gene2B_anno", "Fusion_anno"))
  expect_equal(nrow(standardFusioncalls), 26)
  expect_equal(unique(standardFusioncalls$Sample), "BS_W97QQYKQ")
})

standardFusioncalls <- annoFuse::annoFuse_single_sample(
  # Example files are provided in extdata, at-least 1 fusionfile is required along
  # with its rsem expression file
  fusionfileArriba = system.file("extdata", "arriba_example.tsv", package = "annoFuseData"),
  fusionfileStarFusion = system.file("extdata", "starfusion_example.tsv", package = "annoFuseData"),
  expressionFile = system.file("extdata", "example.rsem.genes.results.gz", package = "annoFuseData"),
  tumorID = "BS_W97QQYKQ",
  # multiple read flag values for filtering using FusionAnnotator values
  artifactFilter = "GTEx_Recurrent|DGD_PARALOGS|Normal|BodyMap|ConjoinG",
  # keep all in-frame , frameshift and other types of Fusion_Type
  readingFrameFilter = "in-frame|frameshift|other",
  # keep all fusions with atleast 100 junction read support (no fusion satisfies this condition so expecting empty file)
  junctionReadCountFilter = 1,
  # keep only fusions where spanningFragCount-junctionReadCountFilter less than equal to 10
  spanningFragCountFilter = 100,
  # keep read throughs
  readthroughFilter = TRUE
)

test_that("Standard fusion single sample run removing read-throughs", {
  expect_equal(colnames(standardFusioncalls), c("LeftBreakpoint", "RightBreakpoint", "FusionName", "Sample", "Caller", "Fusion_Type", "JunctionReadCount", "SpanningFragCount", "Confidence", "annots", "Gene1A", "Gene2A", "Gene1B", "Gene2B", "BreakpointLocation", "SpanningDelta", "reciprocal_exists", "Gene1A_anno", "Gene1B_anno", "Gene2A_anno", "Gene2B_anno", "Fusion_anno"))
  expect_equal(nrow(standardFusioncalls), 16)
})

standardFusioncalls <- annoFuse::annoFuse_single_sample(
  # Example files are provided in extdata, at-least 1 fusionfile is required along
  # with its rsem expression file
  fusionfileArriba = system.file("extdata", "arriba_example.tsv", package = "annoFuseData"),
  fusionfileStarFusion = system.file("extdata", "starfusion_example.tsv", package = "annoFuseData"),
  expressionFile = system.file("extdata", "example.rsem.genes.results.gz", package = "annoFuseData"),
  tumorID = "BS_W97QQYKQ",
  # multiple read flag values for filtering using FusionAnnotator values
  artifactFilter = "GTEx_Recurrent|DGD_PARALOGS|Normal|BodyMap|ConjoinG",
  # keep all in-frame , frameshift and other types of Fusion_Type
  readingFrameFilter = "in-frame|frameshift|other",
  # keep all fusions with atleast 10
  junctionReadCountFilter = 10,
  # keep only fusions where spanningFragCount-junctionReadCountFilter less than equal to 10
  spanningFragCountFilter = 100,
  # keep read throughs
  readthroughFilter = FALSE
)

test_that("Standard fusion single sample run junction read count >= 10", {
  expect_equal(colnames(standardFusioncalls), c("LeftBreakpoint", "RightBreakpoint", "FusionName", "Sample", "Caller", "Fusion_Type", "JunctionReadCount", "SpanningFragCount", "Confidence", "annots", "Gene1A", "Gene2A", "Gene1B", "Gene2B", "BreakpointLocation", "SpanningDelta", "reciprocal_exists", "Gene1A_anno", "Gene1B_anno", "Gene2A_anno", "Gene2B_anno", "Fusion_anno"))
  expect_equal(nrow(standardFusioncalls), 5)
})

standardFusioncalls <- annoFuse::annoFuse_single_sample(
  # Example files are provided in extdata, at-least 1 fusionfile is required along
  # with its rsem expression file
  fusionfileArriba = system.file("extdata", "arriba_example.tsv", package = "annoFuseData"),
  fusionfileStarFusion = system.file("extdata", "starfusion_example.tsv", package = "annoFuseData"),
  expressionFile = system.file("extdata", "example.rsem.genes.results.gz", package = "annoFuseData"),
  tumorID = "BS_W97QQYKQ",
  # multiple read flag values for filtering using FusionAnnotator values
  artifactFilter = "GTEx_Recurrent|DGD_PARALOGS|Normal|BodyMap|ConjoinG",
  # keep all in-frame , frameshift and other types of Fusion_Type
  readingFrameFilter = "in-frame",
  # keep all fusions with atleast 100 junction read support (no fusion satisfies this condition so expecting empty file)
  junctionReadCountFilter = 1,
  # keep only fusions where spanningFragCount-junctionReadCountFilter less than equal to 10
  spanningFragCountFilter = 100,
  # keep read throughs
  readthroughFilter = FALSE
)

test_that("Standard fusion single sample run to only keep in-frame fusions", {
  expect_equal(colnames(standardFusioncalls), c("LeftBreakpoint", "RightBreakpoint", "FusionName", "Sample", "Caller", "Fusion_Type", "JunctionReadCount", "SpanningFragCount", "Confidence", "annots", "Gene1A", "Gene2A", "Gene1B", "Gene2B", "BreakpointLocation", "SpanningDelta", "reciprocal_exists", "Gene1A_anno", "Gene1B_anno", "Gene2A_anno", "Gene2B_anno", "Fusion_anno"))
  expect_equal(nrow(standardFusioncalls), 7)
})



standardFusioncalls <- annoFuse::annoFuse_single_sample(
  # Example files are provided in extdata, at-least 1 fusionfile is required along
  # with its rsem expression file
  fusionfileArriba = system.file("extdata", "arriba_example.tsv", package = "annoFuseData"),
  fusionfileStarFusion = system.file("extdata", "starfusion_example.tsv", package = "annoFuseData"),
  expressionFile = system.file("extdata", "example.rsem.genes.results.gz", package = "annoFuseData"),
  tumorID = "BS_W97QQYKQ",
  # multiple read flag values for filtering using FusionAnnotator values
  artifactFilter = "GTEx_Recurrent|DGD_PARALOGS|Normal|BodyMap|ConjoinG",
  # keep all in-frame , frameshift and other types of Fusion_Type
  readingFrameFilter = "in-frame|frameshift|other",
  # keep all fusions with atleast 100 junction read support (no fusion satisfies this condition so expecting empty file)
  junctionReadCountFilter = 100,
  # keep only fusions where spanningFragCount-junctionReadCountFilter less than equal to 10
  spanningFragCountFilter = 100,
  # keep read throughs
  readthroughFilter = FALSE
)

test_that("Standard fusion single sample run where no fusion passes QC", {
  expect_equal(colnames(standardFusioncalls), c("LeftBreakpoint", "RightBreakpoint", "FusionName", "Sample", "Caller", "Fusion_Type", "JunctionReadCount", "SpanningFragCount", "Confidence", "annots", "Gene1A", "Gene2A", "Gene1B", "Gene2B", "BreakpointLocation", "SpanningDelta", "reciprocal_exists", "Gene1A_anno", "Gene1B_anno", "Gene2A_anno", "Gene2B_anno", "Fusion_anno"))
  expect_equal(nrow(standardFusioncalls), 0)
})
