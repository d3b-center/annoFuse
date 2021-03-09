context("Test groupcount in cohort fusion calls")

out_annofuse <- system.file("extdata", "PutativeDriverAnnoFuse.tsv", package = "annoFuseData")
sfc <- read.delim(out_annofuse)
sfc_groupcount <- groupcount_fusion_calls(sfc, group = "broad_histology", 1)

test_that("Count fusion calls per group in a given cohort", {
  expect_equal(colnames(sfc_groupcount), c("FusionName", "group.ct", "Groups"))
  expect_equal(any(sfc_groupcount$group.ct == 1), FALSE)
})
