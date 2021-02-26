context("Test domain annotation for Gene1A and Gene1B ")

fusionfileArriba <- read.delim(system.file("extdata", "arriba_example.tsv", package = "annoFuseData"),check.names = FALSE,stringsAsFactors = F)

sfc <- 
  annoFuse::fusion_standardization(fusion_calls = fusionfileArriba,
                                   caller= "ARRIBA",
                                   tumorID = "BS_W97QQYKQ")
geneListReferenceDataTab <- read.delim(system.file("extdata", "genelistreference.txt", package = "annoFuseData"), stringsAsFactors = FALSE)
fusionReferenceDataTab <- read.delim(system.file("extdata", "fusionreference.txt", package = "annoFuseData"), stringsAsFactors = FALSE)

bioMartDataPfam <- readRDS(system.file("extdata", "pfamDataBioMart.RDS", package = "annoFuseData"))
domain_list_df <- get_Pfam_domain(standardFusioncalls = sfc, bioMartDataPfam = bioMartDataPfam)



test_that("Fusion domain annotation for standardized arriba calls", {
  expect_equal(names(domain_list_df), c("Gene1A","Gene1B"))
  expect_equal(nrow(domain_list_df$Gene1A), 42)
  expect_equal(nrow(domain_list_df$Gene1B), 199)
  expect_equal(nrow(domain_list_df$Gene1B[grep("BRAF", domain_list_df$Gene1B$Gene1B),]),4)
})