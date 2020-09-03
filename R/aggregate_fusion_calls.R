#' Function to aggregate Caller and read support per fusions calls
#'
#' @param standardFusioncalls A dataframe from star fusion or arriba standardized 
#' to run through the filtering steps
#' @param removeother TRUE to remove Fusion_Type="other" and keep only in-frame and 
#' frameshift Default: FALSE
#' @param filterAnnots regex to remove from annots column eg. 
#' LOCAL_REARRANGEMENT|LOCAL_INVERSION 
#' ## TODO: should this mentioned above be e.g. the default value?
#'
#' @export
#'
#' @return Standardized fusion calls with aggregated Caller and read support
#'
#' @examples
#' #' # standardize
#' fusionfileArriba <- read_arriba_calls(
#'   system.file("extdata", "arriba_example.tsv", package = "annoFuse"))
#' fusionfileStarFusion <- read_starfusion_calls(
#'   system.file("extdata", "starfusion_example.tsv", package = "annoFuse"))
#' formattedArriba <- fusion_standardization(fusionfileArriba,
#'                                           caller = "ARRIBA",
#'                                           tumorID = "tumorID")
#' formattedStarFusion <- fusion_standardization(fusionfileStarFusion,
#'                                               caller = "STARFUSION",
#'                                               tumorID = "tumorID")
#' # merge standardized fusion calls
#' standardFusioncalls <- as.data.frame(rbind(formattedStarFusion, formattedArriba))
#' fusionQCFiltered <- fusion_filtering_QC(
#'   standardFusioncalls = standardFusioncalls, 
#'   readingFrameFilter = "in-frame|frameshift|other",
#'   artifactFilter = "GTEx_Recurrent|DGD_PARALOGS|Normal|BodyMap|ConjoinG",
#'   junctionReadCountFilter = 1,
#'   spanningFragCountFilter = 10,
#'   readthroughFilter = TRUE)
#' # aggregate calls 
#' aggregate_fusion_calls(fusionQCFiltered)    
#' 
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
      dplyr::mutate(CalledBy = toString(.data$Caller), caller_count = n())

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
      dplyr::mutate(CalledBy = toString(.data$Caller), caller_count = n())

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
