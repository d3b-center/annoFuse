#' Annotates standardizes fusion calls calls with zscored expression value from either GTEx/cohort
#'
#' The input should have the following standardized
#' columns to run through this GTEx/cohort normalization function
#' "Sample"  Unique SampleIDs used in your RNAseq dataset
#' "FusionName" GeneA--GeneB ,or if fusion is intergenic then Gene1A/Gene2A--GeneB
#'
#' @param standardFusionCalls Annotates standardizes fusion calls from callers STARfusion| Arriba or QC filtered fusion
#' @param zscoreFilter  Zscore value to use as threshold for annotation of differential expression
#' @param saveZscoredMatrix File to save zscored matrix calculated for the normalized data and expression matrix
#' @param normData normalizing expression dataset to calculate zscore
#' @param expressionMatrix Expression matrix associated with the fusion calls
#'
#' @export
#'
#' @return expression_annotated_fusions is a standardized fusion call set with standard
#'
#' @examples
#' out_annofuse <- system.file("extdata", "PutativeDriverAnnoFuse_test_v16.tsv", package = "annoFuse")
#' sfc <- read.delim(out_annofuse)
#' 
#' # TODOTODO
#' # example with full call
zscored_annotation <- function(standardFusionCalls,
                               zscoreFilter,
                               saveZscoredMatrix,
                               normData,
                               expressionMatrix) {

  standardFusioncalls <- .check_annoFuse_calls(standardFusioncalls)
  
  stopifnot(is.numeric(zscoreFilter))
  if (!missing(saveZscoredMatrix)) 
    stopifnot(is.character(saveZscoredMatrix))
  # TODO: checks on other params as well
  
  # expressionMatrix collapsed at gene level like Gtex max rowMean
  expressionMatrixMatched <- expressionMatrix %>%
    unique() %>%
    # means for each row per each gene_id
    dplyr::mutate(means = rowMeans(dplyr::select(.data$., -.data$GeneSymbol, -.data$gene_id, -.data$EnsembleID))) %>%
    # arrange descending mean
    arrange(desc(.data$means)) %>%
    # to keep only first occurence ie. max rowMean per GeneSymbol
    distinct(.data$GeneSymbol, .keep_all = TRUE) %>%
    ungroup() %>%
    dplyr::select(-.data$means, -.data$gene_id, -.data$EnsembleID) %>%
    dplyr::filter(.data$GeneSymbol %in% normData$GeneSymbol) %>%
    tibble::column_to_rownames("GeneSymbol")
  expressionMatrixMatched <- log2(expressionMatrixMatched + 1)

  # gene matched
  # get log transformed GTEx matrix
  normData <- normData %>%
    tibble::column_to_rownames("GeneSymbol") %>%
    as.matrix()
  normData <- normData[rownames(expressionMatrixMatched), ]
  normData <- log2(normData + 1)

  # normData mean and sd
  normData_means <- rowMeans(normData, na.rm = TRUE)
  normData_sd <- apply(normData, 1, .data$sd, na.rm = TRUE)
  # subtract mean
  expressionMatrixzscored <- sweep(expressionMatrixMatched, 1, normData_means, FUN = "-")
  # divide by SD
  expressionMatrixzscored <- sweep(expressionMatrixzscored, 1, normData_sd, FUN = "/") %>%
    mutate(GeneSymbol = rownames(.data)) %>%
    na.omit(.data)

  # To save GTEx/cohort scored matrix
  if (!missing(saveZscoredMatrix)) {
    saveRDS(expressionMatrixzscored, saveZscoredMatrix)
  }

  # get long format to compare to expression
  expression_long_df <- expressionMatrixzscored %>%
    # Get the data into long format
    reshape2::melt(
      variable.name = "Sample",
      value.name = "zscore_value"
    )

  # fusion calls
  fusion_sample_gene_df <- standardFusionCalls %>%
    # We want to keep track of the gene symbols for each sample-fusion pair
    dplyr::select(.data$Sample, .data$FusionName, .data$Gene1A, .data$Gene1B, .data$Gene2A, .data$Gene2B) %>%
    # We want a single column that contains the gene symbols
    tidyr::gather(.data$Gene1A, .data$Gene1B, .data$Gene2A, .data$Gene2B,
      key = .data$gene_position, value = .data$GeneSymbol
    ) %>%
    # Remove columns without gene symbols
    dplyr::filter(.data$GeneSymbol != "") %>%
    dplyr::arrange(.data$Sample, .data$FusionName) %>%
    # Retain only distinct rows
    dplyr::distinct()

  # check columns for Gene1A,Gene2A,Gene1B and Gene2B
  cols <- c(Gene1A = NA, Gene1B = NA, Gene2A = NA, Gene2B = NA)

  expression_annotated_fusions <- fusion_sample_gene_df %>%
    # join the filtered expression values to the data frame keeping track of symbols
    # for each sample-fusion name pair
    dplyr::left_join(expression_long_df, by = c("Sample", "GeneSymbol")) %>%
    # for each sample-fusion name pair, are all genes under the expression threshold?
    dplyr::group_by(.data$FusionName, .data$Sample) %>%
    dplyr::select(.data$FusionName, .data$Sample, .data$zscore_value, .data$gene_position) %>%
    # cast to keep zscore and gene position
    reshape2::dcast(FusionName + Sample ~ gene_position, value.var = "zscore_value") %>%
    # incase Gene2A/B dont exist like in STARfusion calls
    add_column(!!!cols[setdiff(names(cols), names(.data))]) %>%
    # get annotation from z score
    dplyr::mutate(
      note_expression_Gene1A = ifelse((.data$Gene1A > zscoreFilter | .data$Gene1A < -zscoreFilter), "differentially expressed", "no change"),
      note_expression_Gene1B = ifelse((.data$Gene1B > zscoreFilter | .data$Gene1B < -zscoreFilter), "differentially expressed", "no change"),
      note_expression_Gene2A = ifelse((.data$Gene2A > zscoreFilter | .data$Gene2A < -zscoreFilter), "differentially expressed", "no change"),
      note_expression_Gene2B = ifelse((.data$Gene2B > zscoreFilter | .data$Gene2B < -zscoreFilter), "differentially expressed", "no change")
    ) %>%
    dplyr::rename(
      zscore_Gene1A = .data$Gene1A,
      zscore_Gene1B = .data$Gene1B,
      zscore_Gene2A = .data$Gene2A,
      zscore_Gene2B = .data$Gene2B
    ) %>%
    # unique FusionName-Sample rows
    # use this to filter the QC filtered fusion data frame
    dplyr::inner_join(standardFusionCalls, by = c("FusionName", "Sample")) %>%
    dplyr::distinct()

  return(expression_annotated_fusions)
}
