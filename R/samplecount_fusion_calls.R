#' Function to identify fusions called in n samples

#' @param standardFusioncalls A dataframe from star fusion or arriba standardized to run through the filtering steps
#' @param numSample Least number of samples per group that have the fusion. Defaults to 1.
#' @param group column name for grouping variable
#'
#' @export
#'
#' @return Fusions found in atleast n samples
#'
#' @examples
#' out_annofuse <- system.file("extdata", "PutativeDriverAnnoFuse.tsv", package = "annoFuse")
#' sfc <- read.delim(out_annofuse)
#' samplecount_fusion_calls(sfc, group = "Kids_First_Participant_ID")
samplecount_fusion_calls <- function(standardFusioncalls,
                                     numSample = 1,
                                     group) {
  
  standardFusioncalls <- .check_annoFuse_calls(standardFusioncalls)
  stopifnot(is.numeric(numSample))
  stopifnot(is.character(group))

  # Found in at least n samples in each group
  sample.count <- standardFusioncalls %>%
    dplyr::select(.data$FusionName, .data$Sample, !!as.name(group), -.data$Fusion_Type) %>%
    unique() %>%
    group_by(.data$FusionName, !!as.name(group)) %>%
    dplyr::mutate(sample.count = n(), Sample = toString(.data$Sample)) %>%
    dplyr::filter(.data$sample.count > numSample) %>%
    unique() %>%
    mutate(note = paste0("Found in at least ", numSample, " samples in a group")) %>%
    as.data.frame()

  return(sample.count)
}
