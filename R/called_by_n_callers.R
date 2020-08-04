#' Function to identify fusions called by at least n callers

#' @param standardFusioncalls A dataframe from star fusion or arriba standardized to run through the filtering steps
#' @param numCaller Least number of callers that have the fusion
#' 
#' @export
#' 
#' @return Fusions called by n callers

called_by_n_callers<-function(standardFusioncalls=standardFusioncalls,numCaller=numCaller){

# aggregate caller per Sample,FusionName and Fusion_Type
fusion_caller.summary <- standardFusioncalls %>%
    dplyr::select(.data$Sample,.data$FusionName,.data$Caller,.data$Fusion_Type) %>%
    group_by(.data$FusionName, .data$Sample,.data$Fusion_Type ) %>%
    unique() %>%
    dplyr::mutate(CalledBy = toString(.data$Caller), caller.count = n())

# Called by at least n callers
fusion_calls.summary <- fusion_caller.summary %>%
  dplyr::filter(.data$caller.count >= numCaller) %>%
  unique() %>%
  mutate(note=paste0("Called by",numCaller, "callers")) %>%
  as.data.frame()

return (fusion_calls.summary)

}
