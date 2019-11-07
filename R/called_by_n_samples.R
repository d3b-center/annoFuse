#' Function to identify fusions called by at least n callers

#' @param standardFusioncalls A dataframe from star fusion or arriba standardized to run through the filtering steps
#' @param numCaller Least number of callers that have the fusion
#' @return Fusions called by n callers

called_by_n_samples<-function(standardFusioncalls=standardFusioncalls,numCaller=numCaller){

# aggregate caller
fusion_caller.summary <- standardFusioncalls %>%
    dplyr::filter(Fusion_Type != "other") %>%
    dplyr::select(Sample,FusionName,Caller) %>%
    group_by(FusionName, Sample ) %>%
    unique() %>%
    dplyr::mutate(CalledBy = toString(Caller), caller.count = n())

# Called by at least n callers
fusion_calls.summary <- fusion_caller.summary %>%
  dplyr::filter(caller.count >= numCaller) %>%
  unique() %>%
  mutate(note=paste0("Called by",numCaller, "callers")) %>%
  as.data.frame()

return (fusion_calls.summary)

}
