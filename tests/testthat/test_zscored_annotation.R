context("Test zscore comparison fusion calls")

standardFusioncalls <- standardFusioncalls <- annoFuse::annoFuse_single_sample(
  # Example files are provided in extdata, at-least 1 fusionfile is required along 
  # with its rsem expression file
  fusionfileArriba = system.file("extdata", "arriba_example.tsv", package = "annoFuse"),
  fusionfileStarFusion = system.file("extdata", "starfusion_example.tsv", package = "annoFuse"),
  expressionFile = system.file("extdata", "example.rsem.genes.results.gz", package = "annoFuse"),
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

expressionMatrix<-readRDS(system.file("extdata", "expr_collapsed.rds", package = "annoFuse"))
normData<-readRDS(system.file("extdata", "gtex_collapsed.rds", package = "annoFuse"))
zscoredStandardFusioncalls<-zscored_annotation(standardFusioncalls,
                                                zscoreFilter=2,
                                                normData=normData,
                                                expressionMatrix=expressionMatrix)
                                                


test_that("Zscore comparison in standardized fusions", {
  expect_equal(colnames(zscoredStandardFusioncalls), c("FusionName","Sample","zscore_Gene1A","zscore_Gene1B","zscore_Gene2A","zscore_Gene2B","note_expression_Gene1A","note_expression_Gene1B", "note_expression_Gene2A", "note_expression_Gene2B", "LeftBreakpoint","RightBreakpoint" ,"Caller","Fusion_Type","JunctionReadCount","SpanningFragCount","Confidence","annots","Gene1A","Gene2A","Gene1B","Gene2B","BreakpointLocation","SpanningDelta","reciprocal_exists","Gene1A_anno","Gene1B_anno","Gene2A_anno","Gene2B_anno" ,"Fusion_anno"))
})
