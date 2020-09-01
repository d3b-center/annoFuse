


.check_annoFuse_calls <- function(standardFusioncalls) {
  stopifnot(is(standardFusioncalls, "data.frame"))
  
  # with a minimal set of columns to be there...
  cols_fusioncalls <- c("LeftBreakpoint", "RightBreakpoint", "FusionName",
                        "Gene1A", "Gene1B")
  
  stopifnot(all(cols_fusioncalls %in% colnames(standardFusioncalls)))
  
  invisible(standardFusioncalls)
}
