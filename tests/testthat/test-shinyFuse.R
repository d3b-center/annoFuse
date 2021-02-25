context("Testing launch of the shiny app")

out_annofuse <- system.file("extdata", "PutativeDriverAnnoFuse.tsv", package = "annoFuse")

app_obj <- shinyFuse(out_annofuse)

expect_is(app_obj, "shiny.appobj")

