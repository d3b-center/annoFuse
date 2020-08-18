


.check_annoFuse_calls <- function(standardFusioncalls) {
  stopifnot(is(standardFusioncalls, "data.frame"))
  
  cols_fusioncalls <- c("Sample", "LeftBreakpoint", "RightBreakpoint", "FusionName", "Caller",
                        "Fusion_Type", "JunctionReadCount", "SpanningFragCount", "Confidence", 
                        "annots", "Gene1A", "Gene2A", "Gene1B", "Gene2B", 
                        "Gene1A_anno", "Gene1B_anno", "Gene2A_anno", "Gene2B_anno",
                        "Fusion_anno", "CalledBy", "caller.count", "Kids_First_Participant_ID")
  # TODO: are these all the columns that are specified?
  
  stopifnot(all(colnames(standardFusioncalls) %in% cols_fusioncalls))
  
  invisible(standardFusioncalls)
}
