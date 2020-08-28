#' Function to annotate fusion calls

#' @param standardFusioncalls A dataframe from star fusion or arriba standardized to run through the filtering steps
#' @param geneListReferenceDataTab A dataframe with column 1 as GeneName 2 source file 3 type; collapse to summarize type
#' @param fusionReferenceDataTab A dataframe with column 1 as FusionName 2 source file 3 type; collapse to summarize type
#'
#' @export
#'
#' @return Standardized fusion calls annotated with gene list and fusion list provided in reference folder
#'
#' @examples
#' out_annofuse <- system.file("extdata", "PutativeDriverAnnoFuse_test_v16.tsv", package = "annoFuse")
#' sfc <- read.delim(out_annofuse)
#' # TODOTODO: what are some good values for geneListReferenceDataTab and fusionReferenceDataTab?
#' 
annotate_fusion_calls <- function(standardFusioncalls,
                                  geneListReferenceDataTab,
                                  fusionReferenceDataTab) {
  # TODO: I think the check here is not to be done!
  # standardFusioncalls <- .check_annoFuse_calls(standardFusioncalls)
  stopifnot(is(geneListReferenceDataTab, "data.frame"))
  stopifnot(is(fusionReferenceDataTab, "data.frame"))
  
  geneListReferenceDataTab <- geneListReferenceDataTab %>% dplyr::select(-file)
  fusionReferenceDataTab <- fusionReferenceDataTab %>% dplyr::select(-file)
  annotated_filtered_fusions <- standardFusioncalls %>%
    # annotate Gene1A
    dplyr::left_join(geneListReferenceDataTab, by = c("Gene1A" = "Gene_Symbol")) %>%
    dplyr::rename(Gene1A_anno = .data$type) %>%
    # annotate Gene1B
    dplyr::left_join(geneListReferenceDataTab, by = c("Gene1B" = "Gene_Symbol")) %>%
    dplyr::rename(Gene1B_anno = .data$type) %>%
    # annotate Gene2A
    dplyr::left_join(geneListReferenceDataTab, by = c("Gene2A" = "Gene_Symbol")) %>%
    dplyr::rename(Gene2A_anno = .data$type) %>%
    # annotate Gene2B
    dplyr::left_join(geneListReferenceDataTab, by = c("Gene2B" = "Gene_Symbol")) %>%
    dplyr::rename(Gene2B_anno = .data$type) %>%
    # annotate FusionName
    dplyr::left_join(fusionReferenceDataTab, by = c("FusionName" = "FusionName")) %>%
    dplyr::rename(Fusion_anno = .data$type) %>%
    as.data.frame()
  annotated_filtered_fusions <- unique(annotated_filtered_fusions)

  return(annotated_filtered_fusions)
}
