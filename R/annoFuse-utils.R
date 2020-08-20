


.check_annoFuse_calls <- function(standardFusioncalls) {
  stopifnot(is(standardFusioncalls, "data.frame"))
  
  cols_fusioncalls <- c("LeftBreakpoint", "RightBreakpoint", "FusionName",
                        "Fusion_Type", "JunctionReadCount", "SpanningFragCount",
                        "Gene1A", "Gene1B", "Gene1A_anno", "Gene1B_anno",
                        "Gene2A_anno", "Gene2B_anno","Fusion_anno")
  # TODO: are these all the columns that are specified?
  
  stopifnot(all(cols_fusioncalls %in% colnames(standardFusioncalls)))
  
  invisible(standardFusioncalls)
}
