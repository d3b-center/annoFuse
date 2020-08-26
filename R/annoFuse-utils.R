


.check_annoFuse_calls <- function(standardFusioncalls) {
  stopifnot(is(standardFusioncalls, "data.frame"))
  
  cols_fusioncalls <- c("LeftBreakpoint", "RightBreakpoint", "FusionName",
                        "Gene1A", "Gene1B")
  # TODO: are these all the columns that are specified?
  
  stopifnot(all(cols_fusioncalls %in% colnames(standardFusioncalls)))
  
  invisible(standardFusioncalls)
}
