context("Test called by N caller filter")

out_annofuse <- system.file("extdata", "PutativeDriverAnnoFuse.tsv", package = "annoFuse")
sfc <- read.delim(out_annofuse)
sfc_called <- called_by_n_callers(sfc, numCaller = 2)

test_that("Called by N caller output for standard fusion file input", {
  expect_equal(colnames(sfc_called), c("Sample","FusionName","Caller","Fusion_Type","CalledBy","caller_count", "note" ))
})
