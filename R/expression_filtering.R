#' Expression filtering with user provided expression Matrix and standard fusion calls

#' @param standardFusioncalls A dataframe from star fusion or arriba standardized to run through the filtering steps
#' @param expressionMatrix Expression matrix for samples used in cohort for fusion calls
#' @return Standardized fusion calls annotated with gene list and fusion list provided in reference folder

expressionFilterFusion<-function(standardFusioncalls=standardFusioncalls,expressionMatrix=expressionMatrix,expressionFilter=expressionFilter){

  fusion_sample_gene_df <- standardFusioncalls %>%
    # We want to keep track of the gene symbols for each sample-fusion pair
    dplyr::select(Sample, FusionName, Gene1A, Gene1B, Gene2A, Gene2B) %>%
    # We want a single column that contains the gene symbols
    tidyr::gather(Gene1A, Gene1B, Gene2A, Gene2B,
                  key = gene_position, value = GeneSymbol) %>%
    # Get rid of the Gene1A, Gene1B, Gene2A, Gene2B information
    dplyr::select(-gene_position) %>%
    # Remove columns without gene symbols
    dplyr::filter(GeneSymbol != "") %>%
    # This is for illustrations sake, only
    dplyr::arrange(Sample, FusionName) %>%
    # Retain only distinct rows
    dplyr::distinct()


  expression_long_df <- expressionMatrix %>%
    # Keep the gene symbols and the samples themselves
    dplyr::select(-one_of("gene_id","EnsembleID")) %>%
    # Get the data into long format
    reshape2::melt(variable.name = "Sample",
                   value.name = "expression_value") %>%
    # Remove rows with expression that is too low
    dplyr::filter(expression_value > expressionFilter)

  # Error handling
  if (!all(fusion_sample_gene_df$Sample %in% expression_long_df$Sample)) {
    warning("Not all samples in expression file. Only returning fusions for samples in expressionMatrix.")
  }

  expression_filtered_fusions <- fusion_sample_gene_df %>%
    # join the filtered expression values to the data frame keeping track of symbols
    # for each sample-fusion name pair
    dplyr::left_join(expression_long_df, by = c("Sample", "GeneSymbol"))  %>%
    # for each sample-fusion name pair, are all genes under the expression threshold?
    # keep track in `all_low_expression` column
    dplyr::group_by(FusionName, Sample) %>%
    dplyr::mutate(all_low_expression = all(is.na(expression_value))) %>%
    # only keep the rows that *don't* have all below threshold
    dplyr::filter(!all_low_expression) %>%
    # we only need the FusionName and Sample to filter the entire fusion data.frame
    dplyr::select(FusionName, Sample) %>%
    # unique FusionName-Sample rows
    dplyr::distinct() %>%
    # use this to filter the QC filtered fusion data frame
    dplyr::inner_join(QCFiltered, by = c("FusionName", "Sample")) %>%
    # retain the same columns as merge method
    dplyr::select(c('LeftBreakpoint','RightBreakpoint','FusionName' , 'Sample' ,
                    'Caller' ,'Fusion_Type' , 'JunctionReadCount' ,'SpanningFragCount' ,
                    'Confidence' ,'annots','Gene1A','Gene2A','Gene1B','Gene2B'))

  return(expression_filtered_fusion)
}
