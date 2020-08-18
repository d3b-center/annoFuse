#' Function to aggregate Caller and read support per fusions calls
#'
#' @param standardFusioncalls A dataframe from star fusion or arriba standardized to run through the filtering steps
#' @param removeother TRUE to remove Fusion_Type="other" and keep only in-frame and frameshift Default: FALSE
#' @param filterAnnots regex to remove from annots column eg. LOCAL_REARRANGEMENT|LOCAL_INVERSION 
#' ## TODO: should this mentioned above be e.g. the default value?
#'
#' @export
#'
#' @return Standardized fusion calls with aggregated Caller and read support
#'
#' @examples
#' out_annofuse <- system.file("extdata", "PutativeDriverAnnoFuse_test_v14.tsv", package = "annoFuse")
#' sfc <- read.delim(out_annofuse)
#' sfc_aggr <- aggregate_fusion_calls(sfc, removeother = TRUE, 
#'                                    filterAnnots = "LOCAL_REARRANGEMENT|LOCAL_INVERSION")
aggregate_fusion_calls <- function(standardFusioncalls,
                                   removeother = FALSE,
                                   filterAnnots) {
  standardFusioncalls <- .check_annoFuse_calls(standardFusioncalls)
  stopifnot(is.logical(removeother))
  
  if (removeother) {
    # aggregate caller per FusionName, Sample and Fusion_Type
    fusion_caller.summary <- standardFusioncalls %>%
      dplyr::filter(.data$Fusion_Type != "other") %>%
      dplyr::select(.data$Sample, .data$FusionName, .data$Caller, .data$Fusion_Type) %>%
      group_by(.data$FusionName, .data$Sample, .data$Fusion_Type) %>%
      unique() %>%
      dplyr::mutate(CalledBy = toString(.data$Caller), caller.count = n())

    if (!missing(filterAnnots)) {
      # remove fusion within local rearrangement if required for your project
      standardFusioncalls <- standardFusioncalls %>%
        # remove local rearrangement/adjacent genes use "LOCAL_REARRANGEMENT|LOCAL_INVERSION"
        dplyr::filter(!grepl(filterAnnots, .data$annots))
    }

    # to add aggregated caller from fusion_caller.summary
    standardFusioncalls <- standardFusioncalls %>%
      dplyr::filter(.data$Fusion_Type != "other") %>%
      left_join(fusion_caller.summary, by = (c("Sample", "FusionName", "Fusion_Type", "Caller"))) %>%
      unique()

    return(standardFusioncalls)
  }
  else {
    # aggregate caller per FusionName, Sample and Fusion_Type
    fusion_caller.summary <- standardFusioncalls %>%
      dplyr::select(.data$Sample, .data$FusionName, .data$Caller, .data$Fusion_Type) %>%
      group_by(.data$FusionName, .data$Sample, .data$Fusion_Type) %>%
      unique() %>%
      dplyr::mutate(CalledBy = toString(.data$Caller), caller.count = n())

    if (!missing(filterAnnots)) {
      # remove fusion within local rearrangement if required for your project
      standardFusioncalls <- standardFusioncalls %>%
        # remove local rearrangement/adjacent genes use "LOCAL_REARRANGEMENT|LOCAL_INVERSION"
        dplyr::filter(!grepl(filterAnnots, .data$annots))
    }

    # to add aggregated caller from fusion_caller.summary
    standardFusioncalls <- standardFusioncalls %>%
      left_join(fusion_caller.summary, by = (c("Sample", "FusionName", "Fusion_Type", "Caller"))) %>%
      unique()

    return(standardFusioncalls)
  }
}
