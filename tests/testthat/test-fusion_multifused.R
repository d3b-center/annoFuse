context("Test fusion driver")

fusionfileArriba <- read_arriba_calls(system.file("extdata", "arriba_example.tsv", package = "annoFuse"))

sfc <- 
  annoFuse::fusion_standardization(fusion_calls = fusionfileArriba,
                                   caller= "ARRIBA",
                                   tumorID = "BS_W97QQYKQ")

sfc_multifused <- fusion_multifused(sfc, limitMultiFused = 2)

test_that("Multifused gene annotation for standardized arriba calls", {
  expect_equal(colnames(sfc_multifused), c("Sample","GeneSymbol","Gene.ct","note"))
  expect_equal(nrow(sfc_multifused), 2)
  expect_equal(sfc_multifused$GeneSymbol, c("ANK2","TENM1"))
})