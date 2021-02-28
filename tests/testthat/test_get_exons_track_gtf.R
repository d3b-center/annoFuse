context("Test exons dataframe generation")

library(EnsDb.Hsapiens.v86)
download.file("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.primary_assembly.annotation.gtf.gz","gencode.v27.primary_assembly.annotation.gtf.gz")

exons <- get_exons_tracks_gtf(gtf_path = file.path("gencode.v27.primary_assembly.annotation.gtf.gz"))

test_that("Generates exons dataframe", {
  expect_equal(colnames(exons), c( "geneID","contig","strand", "start" ,"end","exonNumber", "transcript", "type", "geneName"  ))
})
