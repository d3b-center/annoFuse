context("Running report")

test_that("Reporting works as expected, and main errors are checked", {
  out_annofuse <- system.file("extdata", "PutativeDriverAnnoFuse.tsv", package = "annoFuseData")
  report_location <- reportFuse(
    out_annofuse = out_annofuse,
    knitr_show_progress = TRUE,
    force_overwrite = TRUE,
    open_after_creating = FALSE
  )

  expect_is(report_location, "character")

  # triggering for a second time without overwrite...
  expect_error(reportFuse(out_annofuse = out_annofuse))

  expect_error(
    reportFuse(
      out_annofuse = out_annofuse,
      output_file = "this.pdf",
      output_format = "html_document"
    )
  )

  expect_error(
    reportFuse(
      out_annofuse = out_annofuse,
      output_format = "word_document"
    )
  )

  expect_error(
    reportFuse(
      out_annofuse = out_annofuse,
      input_rmd = "some/non/existing/template.Rmd"
    )
  )
})
