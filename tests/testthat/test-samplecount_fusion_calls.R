context("Test counting samples with fusion call in cohort ")

out_annofuse <- system.file("extdata", "PutativeDriverAnnoFuse.tsv", package = "annoFuseData")
sfc <- read.delim(out_annofuse)
sfc_samplecount <- samplecount_fusion_calls(sfc, group = "broad_histology")

test_that("Count samples with fusion calls per group in a given cohort", {
  expect_equal(colnames(sfc_samplecount), c("FusionName", "Sample", "broad_histology", "sample.count", "note"))
  expect_equal(sfc_samplecount[which(sfc_samplecount$broad_histology == "Low-grade astrocytic tumor" & sfc_samplecount$FusionName == "KIAA1549--BRAF"), "sample.count"], 117)
})
