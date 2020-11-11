#' Annotates standardizes fusion calls calls with zscored expression value from either GTEx/cohort
#'
#' The input should have the following standardized
#' columns to run through this GTEx/cohort normalization function
#' "Sample"  Unique SampleIDs used in your RNAseq dataset
#' "FusionName" GeneA--GeneB ,or if fusion is intergenic then Gene1A/Gene2A--GeneB
#'
#' @param standardFusioncalls Annotates standardizes fusion calls from callers STARfusion| Arriba or QC filtered fusion
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
#' standardFusioncalls <- annoFuse::annoFuse_single_sample(
#'   # Example files are provided in extdata, at-least 1 fusionfile is required along 
#'   # with its rsem expression file
#'   fusionfileArriba = system.file("extdata", "arriba_example.tsv", package = "annoFuse"),
#'   fusionfileStarFusion = system.file("extdata", "starfusion_example.tsv", package = "annoFuse"),
#'   expressionFile = system.file(
#'   "extdata", "example.rsem.genes.results.gz", package = "annoFuse"),
#'   tumorID = "BS_W97QQYKQ",
#'   # multiple read flag values for filtering using FusionAnnotator values
#'   artifactFilter = "GTEx_Recurrent|DGD_PARALOGS|Normal|BodyMap|ConjoinG",
#'   # keep all in-frame , frameshift and other types of Fusion_Type
#'   readingFrameFilter = "in-frame|frameshift|other",
#'   # keep all fusions with atleast 1 junction read support
#'   junctionReadCountFilter = 1,
#'   # keep only fusions where spanningFragCount-junctionReadCountFilter less than equal to 10
#'   spanningFragCountFilter = 10,
#'   # keep read throughs
#'   readthroughFilter = FALSE
#' )
#' expressionMatrix <- readRDS(system.file("extdata", "expr_collapsed.rds", package = "annoFuse"))
#' normData <- readRDS(system.file("extdata", "gtex_collapsed.rds", package = "annoFuse"))
#' zscoredStandardFusioncalls <- zscored_annotation(standardFusioncalls,
#'                                                  zscoreFilter=2,
#'                                                  normData=normData,
#'                                                  expressionMatrix=expressionMatrix)

zscored_annotation <- function(standardFusioncalls,
                               zscoreFilter,
                               saveZscoredMatrix,
                               normData,
                               expressionMatrix) {
  
  standardFusioncalls <- .check_annoFuse_calls(standardFusioncalls)
  
  stopifnot(is.numeric(zscoreFilter))
  if (!missing(saveZscoredMatrix)) 
    stopifnot(is.character(saveZscoredMatrix))
  # TODO: checks on other params as well
  
  # expect unique gene_id expressionMatrix collapsed at gene level 
  expressionMatrixMatched <- expressionMatrix %>%
    dplyr::filter(.data$gene_id %in% normData$gene_id) %>%
    tibble::column_to_rownames("gene_id")
  expressionMatrixMatched <- log2(expressionMatrixMatched + 1)
  
  # gene matched
  # expect unique gene_id normData collapsed at gene level   
  # get log transformed GTEx/cohort matrix
  normData <- normData %>%
    tibble::column_to_rownames("gene_id") %>%
    as.matrix()
  normData <- normData[rownames(expressionMatrixMatched), ]
  normData <- log2(normData + 1)
  
  # normData mean and sd
  normData_means <- rowMeans(normData, na.rm = TRUE)
  normData_sd <- apply(normData, 1, stats::sd, na.rm = TRUE)
  # subtract mean
  expressionMatrixzscored <- sweep(expressionMatrixMatched, 1, normData_means, FUN = "-")
  # divide by SD
  expressionMatrixzscored <- sweep(expressionMatrixzscored, 1, normData_sd, FUN = "/") %>%
    tibble::rownames_to_column(var="gene_id") %>%
    na.omit()
  
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
  fusion_sample_gene_df <- standardFusioncalls %>%
    # We want to keep track of the gene symbols for each sample-fusion pair
    dplyr::select(.data$Sample, .data$FusionName, .data$Gene1A, .data$Gene1B, .data$Gene2A, .data$Gene2B) %>%
    # We want a single column that contains the gene symbols
    tidyr::gather(Gene1A, Gene1B, Gene2A, Gene2B,
                  key = gene_position, value = gene_id
    ) %>%
    # Remove columns without gene symbols
    dplyr::filter(.data$gene_id != "") %>%
    dplyr::arrange(.data$Sample, .data$FusionName) %>%
    # Retain only distinct rows
    dplyr::distinct()
  
  # check columns for Gene1A,Gene2A,Gene1B and Gene2B
  cols <- c(Gene1A = NA, Gene1B = NA, Gene2A = NA, Gene2B = NA)
  
  expression_annotated_fusions <- fusion_sample_gene_df %>%
    # join the filtered expression values to the data frame keeping track of symbols
    # for each sample-fusion name pair
    dplyr::left_join(expression_long_df, by = c("Sample", "gene_id")) %>%
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
    dplyr::inner_join(standardFusioncalls, by = c("FusionName", "Sample")) %>%
    dplyr::distinct()
  
  return(expression_annotated_fusions)
}