#' Filter standardized fusion calls for driver fusions
#'
#' If standardized fusion calls are annotated using the geneListReferenceDataTab and fusionReferenceDataTab filters out fusion calls where partner genes are not annotated.
#' If standardized fusion is not annotated it will be annotated with geneListReferenceDataTab and fusionReferenceDataTab provided
#'
#' @param standardFusioncalls A dataframe from star fusion or arriba (more callers to be added)
#' @param annotated Boolean if annotated
#' @param geneListReferenceDataTab A dataframe with column 1 as GeneName 2 source file 3 type; collapse to summarize type
#' @param fusionReferenceDataTab A dataframe with column 1 as FusionName 2 source file 3 type; collapse to summarize type
#'
#' @export
#'
#' @return Putative Driver standardized fusion calls annotated with gene list and fusion list provided in reference folder

fusion_driver <- function(standardFusioncalls = standardFusioncalls, 
                          annotated = TRUE, 
                          geneListReferenceDataTab = geneListReferenceDataTab,
                          fusionReferenceDataTab = fusionReferenceDataTab) {
  # fusion_recurrent5_per_sample <- fusion_multifused(standardFusioncalls,limitMultiFused)

  if (annotated) {
    putative_driver_fusions <- standardFusioncalls %>%
      #  dplyr::filter(!Gene1A %in% fusion_recurrent5_per_sample$GeneSymbol |
      #                  !Gene2A %in% fusion_recurrent5_per_sample$GeneSymbol |
      #                  !Gene1B %in% fusion_recurrent5_per_sample$GeneSymbol |
      #                  !Gene2B %in% fusion_recurrent5_per_sample$GeneSymbol) %>%
      dplyr::filter(!is.na(.data$Gene1A_anno) | !is.na(.data$Gene1B_anno) | !is.na(.data$Gene2A_anno) | !is.na(.data$Gene2B_anno))
  } else {
    standardFusioncalls <- annotate_fusion_calls(standardFusioncalls = standardFusioncalls, geneListReferenceDataTab = geneListReferenceDataTab, fusionReferenceDataTab = fusionReferenceDataTab)
    putative_driver_fusions <- standardFusioncalls %>%
      #    dplyr::filter(!Gene1A %in% fusion_recurrent5_per_sample$GeneSymbol |
      #                    !Gene2A %in% fusion_recurrent5_per_sample$GeneSymbol |
      #                    !Gene1B %in% fusion_recurrent5_per_sample$GeneSymbol |
      #                    !Gene2B %in% fusion_recurrent5_per_sample$GeneSymbol) %>%
      dplyr::filter(!is.na(.data$Gene1A_anno) | !is.na(.data$Gene1B_anno) | !is.na(.data$Gene2A_anno) | !is.na(.data$Gene2B_anno))
  }
  
  return(putative_driver_fusions)
}
