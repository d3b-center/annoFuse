#' Expression filtering with user provided expression Matrix and standard fusion calls

#' @param standardFusioncalls A dataframe from star fusion or arriba standardized to run through the filtering steps
#' @param expressionMatrix Expression matrix for samples used in cohort for fusion calls
#' @param expressionFilter FPKM/TPM threshold for not expressed
#'
#' @export
#'
#' @return Standardized fusion calls annotated with gene list and fusion list provided in reference folder
#'
#' @examples
#' out_annofuse <- system.file("extdata", "PutativeDriverAnnoFuse_test_v14.tsv", package = "annoFuse")
#' sfc <- read.delim(out_annofuse)
#' # TODOTODO: what are some good values for expressionMatrix and expressionFilter?
#' 
expression_filter_fusion <- function(standardFusioncalls,
                                     expressionMatrix,
                                     expressionFilter) {
  # TODO: needs to be restored, otherwise the vignette halts
  # TODO: looks like the file created in the vignette is missing one of the expected columns,
  # or I did put too many mandatory columns!
  # standardFusioncalls <- .check_annoFuse_calls(standardFusioncalls)
  
  # TODO: check if e.g. expressionMatrix is a matrix, and if the filter is a numeric?
  
  fusion_sample_gene_df <- standardFusioncalls %>%
    # We want to keep track of the gene symbols for each sample-fusion pair
    dplyr::select(.data$Sample, .data$FusionName, .data$Gene1A, .data$Gene1B, .data$Gene2A, .data$Gene2B) %>%
    # We want a single column that contains the gene symbols
    tidyr::gather(Gene1A, Gene1B, Gene2A, Gene2B,
      key = gene_position, value = GeneSymbol
    ) %>%
    # Get rid of the Gene1A, Gene1B, Gene2A, Gene2B information
    dplyr::select(-.data$gene_position) %>%
    # Remove columns without gene symbols
    dplyr::filter(.data$GeneSymbol != "") %>%
    # This is for illustrations sake, only
    dplyr::arrange(.data$Sample, .data$FusionName) %>%
    # Retain only distinct rows
    dplyr::distinct()

  expression_long_df <- expressionMatrix %>%
    # Keep the gene symbols and the samples themselves
    dplyr::select(-one_of("gene_id", "EnsembleID")) %>%
    # Get the data into long format
    reshape2::melt(
      variable.name = "Sample",
      value.name = "expression_value"
    )

  # Error handling
  if (!all(fusion_sample_gene_df$Sample %in% expression_long_df$Sample)) {
    warning("Not all samples in expression file. Only returning fusions for samples in expressionMatrix.")
  }
  if (length(unique(setdiff(fusion_sample_gene_df$GeneSymbol, expression_long_df$GeneSymbol))) > 0) {
    warning("Not all genes in expression file. Adding genes unique to fusion in output.")
    # expression_long_df_add<-unique(fusion_sample_gene_df[which(fusion_sample_gene_df$GeneSymbol %in% unique(setdiff(fusion_sample_gene_df$GeneSymbol,expression_long_df$GeneSymbol))),c("Sample","GeneSymbol")])
    # expression_long_df_add$expression_value<--Inf
  }

  expression_long_df <- expression_long_df %>%
    # Remove rows with expression that is too low
    dplyr::filter(.data$expression_value > expressionFilter)


  # expression_long_df<-rbind(expression_long_df,expression_long_df_add)

  expression_filtered_fusions <- fusion_sample_gene_df %>%
    # join the filtered expression values to the data frame keeping track of symbols
    # for each sample-fusion name pair
    dplyr::left_join(expression_long_df, by = c("Sample", "GeneSymbol")) %>%
    # for each sample-fusion name pair, are all genes under the expression threshold?
    # keep track in `all_low_expression` column
    dplyr::group_by(.data$FusionName, .data$Sample) %>%
    dplyr::mutate(all_low_expression = all(is.na(.data$expression_value))) %>%
    # only keep the rows that *don't* have all below threshold
    dplyr::filter(!.data$all_low_expression) %>%
    # we only need the FusionName and Sample to filter the entire fusion data.frame
    dplyr::select(.data$FusionName, .data$Sample) %>%
    # unique FusionName-Sample rows
    dplyr::distinct() %>%
    # use this to filter the QC filtered fusion data frame
    dplyr::inner_join(standardFusioncalls, by = c("FusionName", "Sample")) %>%
    # retain the same columns as merge method
    dplyr::select(c(
      "LeftBreakpoint", "RightBreakpoint", "FusionName", "Sample",
      "Caller", "Fusion_Type", "JunctionReadCount", "SpanningFragCount",
      "Confidence", "annots", "Gene1A", "Gene2A", "Gene1B", "Gene2B"
    ))

  return(expression_filtered_fusions)
}
