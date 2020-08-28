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
#' # standardize
#' fusionfileArriba = read.delim(system.file("extdata", "arriba_example.tsv", package = "annoFuse"),stringsAsFactors = FALSE)
#' fusionfileStarFusion = read.delim(system.file("extdata", "starfusion_example.tsv", package = "annoFuse"),stringsAsFactors = FALSE)
#' formattedArriba = fusion_standardization(fusionfileArriba,caller="ARRIBA",tumorID = "tumorID")
#' formattedStarFusion = fusion_standardization(fusionfileStarFusion,caller="STARFUSION",tumorID = "tumorID")
#' # merge standardized fusion calls
#' standardFusioncalls <- rbind(formattedStarFusion, formattedArriba) %>% as.data.frame()
#' fusionQCFiltered <- fusion_filtering_QC(standardFusioncalls = standardFusioncalls, 
#'                     readingFrameFilter = "in-frame|frameshift|other",
#'                     artifactFilter = "GTEx_Recurrent|DGD_PARALOGS|Normal|BodyMap|ConjoinG",
#'                     junctionReadCountFilter = 1,
#'                     spanningFragCountFilter = 10,
#'                     readthroughFilter = TRUE)
#' # annotated from gene and fusion refrence list
#' # read in gene and fusion reference tab
#' geneListReferenceDataTab <- read.delim(system.file("extdata", "genelistreference.txt", package = "annoFuse"), stringsAsFactors = FALSE)
#'# column 1 as FusionName 2 source file 3 type; collapse to summarize type
#'fusionReferenceDataTab <- read.delim(system.file("extdata", "fusionreference.txt", package = "annoFuse"), stringsAsFactors = FALSE)
#'filteredFusionAnnotated <- annotate_fusion_calls(standardFusioncalls = fusionQCFiltered, geneListReferenceDataTab = geneListReferenceDataTab, fusionReferenceDataTab = fusionReferenceDataTab)                                   
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
