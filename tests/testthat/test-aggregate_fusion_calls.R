context("Test aggregate fusion calls")

standardFusioncalls <- standardFusioncalls <- annoFuse::annoFuse_single_sample(
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

aggregate_fusion_calls_out <- aggregate_fusion_calls(standardFusioncalls)

test_that("Aggregate fusion calls output from non-empty arriba and starfusion fusion calls", {
  expect_equal(
    colnames(aggregate_fusion_calls_out), 
    c("LeftBreakpoint","RightBreakpoint","FusionName","Sample","Caller","Fusion_Type",
      "JunctionReadCount","SpanningFragCount","Confidence","annots","Gene1A","Gene2A",
      "Gene1B","Gene2B","BreakpointLocation","SpanningDelta","reciprocal_exists",
      "Gene1A_anno","Gene1B_anno","Gene2A_anno","Gene2B_anno" ,"Fusion_anno", 
      "CalledBy","caller_count" )
    )
  expect_equal(aggregate_fusion_calls_out[21,"caller_count"],2)
})

test_that("keep only in-frame and frameshift fusions",{
  afc_noother <- aggregate_fusion_calls(standardFusioncalls, removeother = TRUE)
  dim(afc_noother)
  expect_equal(nrow(afc_noother), 9)
})
