context("Testing plot functions - plot_recurrent_genes")

test_that("Plotting recurrent genes", {
  out_annofuse <- system.file("extdata", "PutativeDriverAnnoFuse.tsv", package = "annoFuse")
  sfc <- read.delim(out_annofuse, stringsAsFactors = FALSE)
  # keep only in-frame and fusions where both breakpoints are within genes
  sfc <- as.data.frame(
    sfc[ which(sfc$Fusion_Type == "in-frame" & sfc$BreakpointLocation == "Genic"),])
  
  p <- plot_recurrent_genes(sfc, 
                            groupby = "broad_histology", 
                            countID = "Kids_First_Participant_ID")
  
  expect_is(p, "gg")
})
