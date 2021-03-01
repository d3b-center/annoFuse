context("Test domain annotation for Gene1A and Gene1B ")

fusionfileArriba <- read.delim(system.file("extdata", "arriba_example.tsv", package = "annoFuse"),check.names = FALSE,stringsAsFactors = F)

sfc <- 
  annoFuse::fusion_standardization(fusion_calls = fusionfileArriba,
                                   caller= "ARRIBA",
                                   tumorID = "BS_W97QQYKQ")
geneListReferenceDataTab <- read.delim(system.file("extdata", "genelistreference.txt", package="annoFuse"), stringsAsFactors = FALSE)
fusionReferenceDataTab <- read.delim(system.file("extdata", "fusionreference.txt", package="annoFuse"), stringsAsFactors = FALSE)

ensembl <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl",host="ensembl.org")

bioMartDataPfam <- get_biomart_pfam_merge(ensembl = ensembl,pfamDesc_path = "pfamDesc.txt.gz",ucscGenePfam_path = "ucscGenePfam.txt.gz")
domain_list_df <- get_Pfam_domain(standardFusioncalls = sfc, bioMartDataPfam = bioMartDataPfam)



test_that("Fusion domain annotation for standardized arriba calls", {
  expect_equal(names(domain_list_df), c("Gene1A","Gene1B"))
  expect_equal(nrow(domain_list_df$Gene1A), 38)
  expect_equal(nrow(domain_list_df$Gene1B), 262)
  expect_equal(nrow(domain_list_df$Gene1B[grep("BRAF", domain_list_df$Gene1B$Gene1B),]),5)
})