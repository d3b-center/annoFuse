context("Testing launch of the shiny app")

test_that("Testing shiny app and info generation elements", {
  out_annofuse <- system.file("extdata", "PutativeDriverAnnoFuse.tsv", package = "annoFuse")
  
  app_obj <- shinyFuse(out_annofuse)
  expect_is(app_obj, "shiny.appobj")
  
  out_annofuse_content <- read.delim(out_annofuse, stringsAsFactors = FALSE)
  expect_error(shinyFuse(out_annofuse_content))
  
  expect_error(shinyFuse("/some/data/not/there"))
  
  ginfo2 <- annoFuse:::doublegeneinfo_2_html("ACTB", "GAPDH")
  expect_is(ginfo2, "html")
  
  ginfo <- annoFuse:::geneinfo_2_html("ACTB")
  expect_is(ginfo, "html")
})
