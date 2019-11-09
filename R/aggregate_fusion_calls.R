#' Function to aggregate Caller and read support per fusions calls
#'
#' @param standardFusioncalls A dataframe from star fusion or arriba standardized to run through the filtering steps
#' @param removeother TRUE to remove Fusion_Type="other" and keep only in-frame and frameshift
#' @param filterAnnots regex to remove from annots column eg. LOCAL_REARRANGEMENT|LOCAL_INVERSION
#' @return Standardized fusion calls with aggregated Caller and read support

aggregate_fusion_calls<-function(standardFusioncalls=standardFusioncalls,removeother=TRUE,filterAnnots=filterAnnots){

  if (!missing(filterAnnots)){
    # remove fusion within local rearrangement if required for your project
    standardFusioncalls <- standardFusioncalls %>%
      # remove local rearrangement/adjacent genes use "LOCAL_REARRANGEMENT|LOCAL_INVERSION"
      dplyr::filter(!grepl(filterAnnots,annots))
  }


  if( removeother){
  # aggregate caller
  fusion_caller.summary <- standardFusioncalls %>%
    dplyr::filter(Fusion_Type != "other") %>%
    dplyr::select(Sample,FusionName,Caller) %>%
    group_by(FusionName, Sample ) %>%
    unique() %>%
    dplyr::mutate(CalledBy = toString(Caller), caller.count = n())

  # aggregate by read count
  fusion_read.summary <- standardFusioncalls %>%
    dplyr::filter(Fusion_Type != "other") %>%
    dplyr::select(Sample,FusionName,JunctionReadCount,SpanningFragCount) %>%
    group_by(FusionName, Sample ) %>%
    unique() %>%
    dplyr::mutate(JunctionReadCountSum=sum(JunctionReadCount),SpanningFragCountSum=sum(SpanningFragCount))

  fusion_summary<-left_join(fusion_caller.summary,fusion_read.summary,by=c("Sample","FusionName"))

  standardFusioncalls<-standardFusioncalls %>%
    dplyr::filter(Fusion_Type != "other") %>%
    left_join(fusion_summary,by=(c("Sample","FusionName","Caller","JunctionReadCount","SpanningFragCount")))

  return(standardFusioncalls)
  }
  else {
    fusion_caller.summary <- standardFusioncalls %>%
      dplyr::select(Sample,FusionName,Caller) %>%
      group_by(FusionName, Sample ) %>%
      unique() %>%
      dplyr::mutate(CalledBy = toString(Caller), caller.count = n())

    # aggregate by read count
    fusion_read.summary <- standardFusioncalls %>%
      dplyr::select(Sample,FusionName,JunctionReadCount,SpanningFragCount) %>%
      group_by(FusionName, Sample ) %>%
      unique() %>%
      dplyr::mutate(JunctionReadCountSum=sum(JunctionReadCount),SpanningFragCountSum=sum(SpanningFragCount))

    fusion_summary<-left_join(fusion_caller.summary,fusion_read.summary,by=c("Sample","FusionName"))

    standardFusioncalls<-standardFusioncalls %>%
      left_join(fusion_summary,by=(c("Sample","FusionName","Caller","JunctionReadCount","SpanningFragCount")))


    return(standardFusioncalls)

  }

}
