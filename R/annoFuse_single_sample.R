#' Single Sample use for annoFuse
#'
#' Performs artifact filter to remove readthroughs,red flags and performs expression
#' filtering with user provided expression Matrix and expression threshold
#' 
#' @param fusionfileArriba A dataframe from arriba fusion caller
#' @param fusionfileStarFusion A dataframe from starfusion caller
#' @param expressionFile Expression matrix for samples used in cohort for fusion calls
#' @param expressionFilter FPKM/TPM threshold for not expressed
#' @param tumorID Sample name to be used
#' @param readingFrameFilter A regex to capture readingframe (eg. in-frame|frameshift|other)
#' @param readthroughFilter Boolean for filtering readthroughs
#' @param artifactFilter A red flag filter from Annotation ; in OpenPBTA annotation 
#' is from FusionAnnotator column "annots"
#' @param junctionReadCountFilter An integer threshold for JunctionReadCount
#' @param spanningFragCountFilter An integer threshold for (SpanningFragCount -
#' JunctionReadCount)
#'
#' @export
#'
#' @return Standardized fusion calls annotated with gene list and fusion list provided in reference folder
#'
#' @examples
#' standardFusioncalls <- annoFuse::annoFuse_single_sample(
#'   # Example files are provided in extdata, at-least 1 fusionfile is required along 
#'   # with its rsem expression file
#'   fusionfileArriba = system.file(
#'     "extdata", "arriba_example.tsv", package = "annoFuseData"),
#'   fusionfileStarFusion = system.file(
#'     "extdata", "starfusion_example.tsv", package = "annoFuseData"),
#'   expressionFile = system.file(
#'     "extdata", "example.rsem.genes.results.gz", package = "annoFuseData"),
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
annoFuse_single_sample <- function(fusionfileArriba,
                                   fusionfileStarFusion,
                                   expressionFile = NULL,
                                   expressionFilter = 1,
                                   tumorID = "tumorID",
                                   readingFrameFilter = "in-frame|frameshift|other",
                                   readthroughFilter = FALSE,
                                   artifactFilter = "GTEx_Recurrent|DGD_PARALOGS|Normal|BodyMap|ConjoinG",
                                   junctionReadCountFilter = 1,
                                   spanningFragCountFilter = 10) {
  
  stopifnot(is.character(fusionfileArriba))
  stopifnot(is.character(fusionfileStarFusion))
  stopifnot(is.numeric(expressionFilter))
  stopifnot(is.character(tumorID))
  stopifnot(is.character(readingFrameFilter))
  stopifnot(is.logical(readthroughFilter))
  stopifnot(is.character(artifactFilter))
  stopifnot(is.numeric(junctionReadCountFilter))
  stopifnot(is.numeric(spanningFragCountFilter))
  
  # read files
  STARFusioninputfile <- read_starfusion_calls(fusionfileStarFusion)
  Arribainputfile <- read_arriba_calls(fusionfileArriba)

  # read in gene and fusion reference tab
  geneListReferenceDataTab <- read.delim(system.file("extdata", "genelistreference.txt", package = "annoFuseData"), stringsAsFactors = FALSE)
  
  # column 1 as FusionName 2 source file 3 type; collapse to summarize type
  fusionReferenceDataTab <- read.delim(system.file("extdata", "fusionreference.txt", package = "annoFuseData"), stringsAsFactors = FALSE)
  

  # if StarFusion and Arriba files empty execution stops
  if (is_empty(STARFusioninputfile$`#FusionName`) & is_empty(Arribainputfile$`#gene1`)) {
    warning("StarFusion and Arriba files empty")
    standardFusioncalls<-data.frame("LeftBreakpoint"=character(),
                                       "RightBreakpoint"=character(),
                                       "FusionName"= character(),
                                       "Sample" = character(),
                                       "Caller" = character(),
                                       "Fusion_Type" = character(),
                                       "JunctionReadCount" = numeric(),
                                       "SpanningFragCount" = numeric(),
                                       "Confidence" = character(),
                                       "annots" = character(),
                                       "Gene1A" = character(),
                                       "Gene2A" = character(),
                                       "Gene1B" = character(),
                                       "Gene2B" = character(),
                                       "BreakpointLocation" = character(),
                                       "SpanningDelta" = numeric(),
                                       "reciprocal_exists" = character(),
                                       "Gene1A_anno" = character(),
                                       "Gene1B_anno" = character(),
                                       "Gene2A_anno" = character(),
                                       "Gene2B_anno" = character(),
                                       "Fusion_anno" = character())
    return(standardFusioncalls)
  }

  # if StarFusion and Arriba files are not empty
  if (!is_empty(STARFusioninputfile$`#FusionName`) & !is_empty(Arribainputfile$`#gene1`)) {

    # standardized fusion calls
    standardizedSTARFusion <- fusion_standardization(fusion_calls = STARFusioninputfile, caller = "STARFUSION",tumorID = tumorID)
    standardizedArriba <- fusion_standardization(fusion_calls = Arribainputfile, caller = "ARRIBA",tumorID = tumorID)

    # merge standardized fusion calls
    standardFusioncalls <- rbind(standardizedSTARFusion, standardizedArriba) %>% as.data.frame()
  }

  # if StarFusion file is empty only run standardization for Arriba calls
  if (is_empty(STARFusioninputfile$`#FusionName`) & !is_empty(Arribainputfile$`#gene1`)) {
    warning(paste("No fusion calls in StarFusion "))


    # standardized fusion calls
    standardizedArriba <- fusion_standardization(fusion_calls = Arribainputfile, caller = "ARRIBA",tumorID = tumorID)

    # standardized fusion calls
    standardFusioncalls <- standardizedArriba %>% as.data.frame()
  }

  # if Arriba file is empty only run standardization for StarFusion calls
  if (!is_empty(STARFusioninputfile$`#FusionName`) & is_empty(Arribainputfile$`#gene1`)) {
    warning(paste("No fusion calls in Arriba "))


    # standardized fusion calls
    standardizedSTARFusion <- fusion_standardization(fusion_calls = STARFusioninputfile, caller = "STARFUSION",tumorID = tumorID)

    # standardized fusion calls
    standardFusioncalls <- standardizedSTARFusion %>% as.data.frame()
  }


  # General fusion QC for read support and red flags
  fusionQCFiltered <- fusion_filtering_QC(standardFusioncalls = standardFusioncalls, readingFrameFilter = readingFrameFilter, artifactFilter = artifactFilter, junctionReadCountFilter = junctionReadCountFilter, spanningFragCountFilter = spanningFragCountFilter, readthroughFilter = readthroughFilter)

  if (!is.null(expressionFile)) {
    expressionMatrix <- read_tsv(expressionFile)

    # split gene id and symbol
    expressionMatrix <- expressionMatrix %>%
      dplyr::mutate(gene_id = str_replace(gene_id, "_PAR_Y_", "_"))

    expressionMatrix <- cbind(expressionMatrix, colsplit(expressionMatrix$gene_id, pattern = "_", names = c("EnsembleID", "GeneSymbol")))


    # collapse to matrix of HUGO symbols x Sample identifiers
    # take max expression per row and use the max value for duplicated gene symbols
    expressionMatrix.collapsed <- expressionMatrix %>%
      arrange(desc(FPKM)) %>% # arrange decreasing by FPKM
      distinct(GeneSymbol, .keep_all = TRUE) %>% # keep the ones with greatest FPKM value. If ties occur, keep the first occurencce
      unique() %>%
      remove_rownames() %>%
      dplyr::select(.data$EnsembleID, .data$GeneSymbol, .data$FPKM, .data$gene_id)

    # rename columns
    colnames(expressionMatrix.collapsed)[3] <- tumorID

    expressionFiltered <- annoFuse::expression_filter_fusion(standardFusioncalls = fusionQCFiltered, expressionFilter = expressionFilter, expressionMatrix = expressionMatrix.collapsed)

    # annotated QC and expression filtered fusion calls
    filteredFusionAnnotated <- annotate_fusion_calls(standardFusioncalls = expressionFiltered, geneListReferenceDataTab = geneListReferenceDataTab, fusionReferenceDataTab = fusionReferenceDataTab)
  } else {
    # annotated QC filtered fusion calls
    filteredFusionAnnotated <- annotate_fusion_calls(standardFusioncalls = fusionQCFiltered, geneListReferenceDataTab = geneListReferenceDataTab, fusionReferenceDataTab = fusionReferenceDataTab)
  }

  # return filtered and annotated dataframe
  return(filteredFusionAnnotated)
}
