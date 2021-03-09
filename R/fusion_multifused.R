#' Function to identify fusions called by at least n callers

#' @param standardFusioncalls A dataframe from star fusion or arriba standardized to run through the filtering steps
#' @param limitMultiFused Integer to identify a limit of times a gene can be fused per sample
#'
#' @export
#'
#' @return Fusions where gene partner(s) is multifused per sample
#'
#' @examples
#' out_annofuse <- system.file("extdata", "PutativeDriverAnnoFuse.tsv", package = "annoFuseData")
#' sfc <- read.delim(out_annofuse, stringsAsFactors = FALSE)
#' sfc_multifused <- fusion_multifused(sfc, limitMultiFused = 2)
fusion_multifused <- function(standardFusioncalls,
                              limitMultiFused = 3) {
  standardFusioncalls <- .check_annoFuse_calls(standardFusioncalls)
  stopifnot(is.numeric(limitMultiFused))

  # remove multi-fused genes
  fusion_multifused_per_sample <- standardFusioncalls %>%
    # We want to keep track of the gene symbols for each sample-fusion pair
    dplyr::select(.data$Sample, .data$FusionName, .data$Gene1A, .data$Gene1B, .data$Gene2A, .data$Gene2B) %>%
    # We want a single column that contains the gene symbols
    tidyr::gather(Gene1A, Gene1B, Gene2A, Gene2B,
      key = gene_position, value = GeneSymbol
    ) %>%
    # Remove columns without gene symbols
    dplyr::filter(.data$GeneSymbol != "") %>%
    dplyr::arrange(.data$Sample, .data$FusionName) %>%
    # Retain only distinct rows
    dplyr::distinct() %>%
    group_by(.data$Sample, .data$GeneSymbol) %>%
    dplyr::summarise(Gene.ct = n()) %>%
    dplyr::filter(.data$Gene.ct > limitMultiFused) %>%
    mutate(note = paste0("multfused ", limitMultiFused, " times per sample"))

  return(fusion_multifused_per_sample)
}
