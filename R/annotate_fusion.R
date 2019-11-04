#' Function to annotate fusion calls

#' @param standardFusioncalls A dataframe from star fusion or arriba standardized to run through the filtering steps
#' @param geneListReferenceDataTab A dataframe with column 1 as GeneName 2 source file 3 type; collapse to summarize type
#' @param fusionReferenceDataTab A dataframe with column 1 as FusionName 2 source file 3 type; collapse to summarize type
#' @return Standardized fusion calls annotated with gene list and fusion list provided in reference folder

annotate_fusion_calls<-function(standardFusioncalls=standardFusioncalls,geneListReferenceDataTab=geneListReferenceDataTab,fusionReferenceDataTab=fusionReferenceDataTab){
    annotated_filtered_fusions<-standardFusioncalls %>%
    # annotate Gene1A
    dplyr::left_join(geneListReferenceDataTab,by=c("Gene1A"="Gene_Symbol")) %>% dplyr::rename(Gene1A_anno=type) %>%
    # annotate Gene1B
    dplyr::left_join(geneListReferenceDataTab,by=c("Gene1B"="Gene_Symbol")) %>% dplyr::rename(Gene1B_anno=type) %>%
    # annotate Gene2A
    dplyr::left_join(geneListReferenceDataTab,by=c("Gene2A"="Gene_Symbol")) %>% dplyr::rename(Gene2A_anno=type) %>%
    # annotate Gene2B
    dplyr::left_join(geneListReferenceDataTab,by=c("Gene2B"="Gene_Symbol")) %>% dplyr::rename(Gene2B_anno=type) %>%
    # annotate FusionName
    dplyr::left_join(fusionReferenceDataTab,by=c("FusionName"="FusionName")) %>% dplyr::rename(Fusion_anno=type) %>%
    as.data.frame()
  return(annotated_filtered_fusions)
}


