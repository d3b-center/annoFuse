

#' report_fuse
#'
#' @details When `report_fuse` is called, a RMarkdown template file could be copied
#' into the output directory, and [rmarkdown::render()] will be called to
#' generate the final report.
#'
#' As a default template, `report_fuse` uses the one delivered together with the
#' `annoFuse` package. Experienced users can take that as a starting point to further
#' edit and customize.
#'
#' @param out_annofuse TODO
#' @param project_id A character string, which can be considered as an identifier
#' for the set/session, and will be e.g. used in the title of the report created
#' via [report_fuse()]
#' @param input_rmd Character string with the path to the RMarkdown (.Rmd) file
#' that will be used as the template for generating the report. Defaults to NULL,
#' which will then use the one provided with the `annoFuse` package.
#' @param output_file Character string, specifying the file name of the output
#' report. The file name extension must be either `.html` or `.pdf`, and consistent
#' with the value of `output_format`.
#' @param output_dir Character, defining the path to the output directory where
#' the report will be generated. Defaults to the temp directory (`tempdir()`).
#' @param output_format The format of the output report. Either `html_document`
#' or `pdf_document`. The file name extension of `output_file` must be consistent
#' with this choice. Can also be left empty and determined accordingly.
#' @param force_overwrite Logical, whether to force overwrite an existing report
#' with the same name in the output directory. Defaults to FALSE.
#' @param knitr_show_progress Logical, whether to display the progress of `knitr`
#' while generating the report. Defaults to FALSE.
#' @param open_after_creating Logical, whether to open the report in the default
#' browser after being generated. Defaults to TRUE.
#' @param ... Other arguments that will be passed to [rmarkdown::render()].
#'
#' @return Generates a fully fledged report in the `output_dir` directory, called
#' `output_file` and returns (invisibly) the name of the generated report.
#' @export
#'
#' @examples
#'
#' out_annofuse <- "PutativeDriverAnnoFuse_test_v14.tsv"
#' \dontrun{
#' report_fuse(out_annofuse = out_annofuse)
#' }
report_fuse <- function(out_annofuse,
                        project_id = "projectID",
                        input_rmd = NULL,
                        output_file = "my_first_annoFuse_report.html",
                        output_dir = tempdir(),
                        output_format = NULL,
                        force_overwrite = FALSE,
                        knitr_show_progress = FALSE,
                        open_after_creating = TRUE,
                        ...) {
  # generates a nice number of outputs, plots, and so on, placed in a report. Boom :)

  # If possible, set output format based on the extension of output_file, if the output format is not provided
  if (is.null(output_format)) {
    if (tools::file_ext(output_file) == "pdf") {
      output_format <- "pdf_document"
    } else {
      output_format <- "html_document"
    }
  }

  # Raise an error if output_format is not one of the allowed
  if (!(output_format %in% c("pdf_document", "html_document"))) {
    stop("The provided output_format is currently not supported. Please ",
      "use either 'html_document' or 'pdf_document'.",
      call. = FALSE
    )
  }

  # Raise an error if the output format and file name extension don't match
  if (output_format != paste0(tools::file_ext(output_file), "_document")) {
    stop("File name extension of output_file (.",
      tools::file_ext(output_file),
      ") doesn't agree with the ",
      "output_format, should be .",
      gsub("_document$", "", output_format),
      call. = FALSE
    )
  }

  # output files
  output_report <- file.path(output_dir, basename(output_file)) # no need of normalizePath?
  output_rmd <- file.path(
    output_dir,
    paste0(tools::file_path_sans_ext(basename(output_file)), ".Rmd")
  )

  # report
  if (file.exists(output_report)) {
    if (!force_overwrite) {
      stop("The file ", output_report,
        " already exists. Please remove or rename the file, provide ",
        "another value of output_file, or set force_overwrite = TRUE.",
        call. = FALSE
      )
    } else {
      warning("The file ", output_report,
        " already exists and will be overwritten, since ",
        "force_overwrite = TRUE.",
        immediate. = TRUE,
        call. = FALSE
      )
    }
  }

  # Rmd template
  if (is.null(input_rmd)) {
    template_rmd <- system.file("extdata",
      "report_template_annoFuse.Rmd",
      package = "annoFuse"
    )
  } else {
    template_rmd <- input_rmd
  }

  if (file.exists(template_rmd)) {
    if (file.exists(output_rmd)) {
      # stop("There is already an .Rmd file called ", output_rmd,
      #      ". Please remove or rename this file, or choose another ",
      #      "output_file name.", call. = FALSE)
    } else {
      # another possible thought: work in a tempdir, that is probably even more elegant
      file.copy(from = template_rmd, to = output_rmd, overwrite = FALSE)
    }
  } else {
    stop("The Rmd template file ", template_rmd, " does not exist.",
      call. = FALSE
    )
  }

  # annofuse_tbl is then used in the Rmd report
  annofuse_tbl <- read.delim(normalizePath(out_annofuse))

  # Process the arguments
  args <- list(...)
  args$input <- output_rmd
  args$output_format <- output_format
  args$output_file <- output_file
  args$quiet <- !knitr_show_progress



  # Render the report
  output_file <- do.call("render", args = args)

  # Remove temporary file
  file.remove(output_rmd)

  # Open up in a browser
  if (open_after_creating) {
    browseURL(output_file)
  }

  invisible(output_file)
}
