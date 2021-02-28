context("Test biomart pfam data")

ensembl <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl",host="ensembl.org")

download.file("http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/pfamDesc.txt.gz","pfamDesc.txt.gz")
download.file("http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ucscGenePfam.txt.gz","ucscGenePfam.txt.gz")

bioMartDataPfam <- get_biomart_pfam_merge(ensembl = ensembl,pfamDesc_path = "pfamDesc.txt.gz",ucscGenePfam_path = "ucscGenePfam.txt.gz")

test_that("Generates pfam biomart data", {
  expect_equal(colnames(bioMartDataPfam), c("hgnc_symbol","pfam_id" ,"chromosome_name" ,"gene_start"  ,"gene_end" ,"strand", "NAME" ,"DESC","domain_chr","domain_start","domain_end" ))
})
