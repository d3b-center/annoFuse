context("Testing plot functions - plot_breakpoints")

test_that("Plot of breakpoints work as expected", {
  out_annofuse <- system.file("extdata", "PutativeDriverAnnoFuse.tsv", package = "annoFuseData")
  sfc <- read.delim(out_annofuse)

  exons <- readRDS(system.file("extdata", "exonsToPlot.RDS", package = "annoFuseData"))

  bioMartDataPfam <-
    readRDS(system.file("extdata", "pfamDataBioMart.RDS", package = "annoFuseData"))

  domainDataFrame <- get_Pfam_domain(
    standardFusioncalls = sfc,
    bioMartDataPfam = bioMartDataPfam,
    keepPartialAnno = TRUE
  )
  left <- plot_breakpoints(
    sampleid = "BS_044XZ8ST",
    domainDataFrame = domainDataFrame,
    exons = exons,
    geneposition = "Left",
    fusionname = "ANTXR1--BRAF",
    leftBreakpoint = "2:69193415"
  )

  right <- plot_breakpoints(
    sampleid = "BS_044XZ8ST",
    domainDataFrame = domainDataFrame,
    exons = exons,
    geneposition = "Right",
    fusionname = "ANTXR1--BRAF",
    rightBreakpoint = "7:140787584"
  )
  expect_is(left, "gg")
  expect_is(right, "gg")
})
