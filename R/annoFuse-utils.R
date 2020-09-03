


.check_annoFuse_calls <- function(standardFusioncalls) {
  stopifnot(is(standardFusioncalls, "data.frame"))
  
  # with a minimal set of columns to be there...
  cols_fusioncalls <- c("LeftBreakpoint", "RightBreakpoint", "FusionName",
                        "Gene1A", "Gene1B")
  
  stopifnot(all(cols_fusioncalls %in% colnames(standardFusioncalls)))
  
  invisible(standardFusioncalls)
}


#' Read in fusion calls from  Arriba v1.1.0 
#'
#' @param arriba_calls Please refer to software documenation [fusions.tsv](https://arriba.readthedocs.io/en/latest/output-files/)
#'
#' @return A data.frame object with correct column specifications
#' @export
#'
#' @examples
#' fusionfileArriba <- read_arriba_calls(
#'                         system.file("extdata", "arriba_example.tsv", package = "annoFuse"))
#'                         

read_arriba_calls <- function(arriba_calls){
  # set col types
  col_types = readr::cols(
    `#gene1` = readr::col_character(),
    gene2 = readr::col_character(),
    `strand1(gene/fusion)` = readr::col_character(),
    `strand2(gene/fusion)` = readr::col_character(),
    breakpoint1 = readr::col_character(),
    breakpoint2 = readr::col_character(),
    site1 = readr::col_character(),
    site2 = readr::col_character(),
    type = readr::col_character(),
    direction1 = readr::col_character(),
    direction2 = readr::col_character(),
    split_reads1 = readr::col_integer(),
    split_reads2 = readr::col_integer(),
    discordant_mates = readr::col_integer(),
    coverage1 = readr::col_integer(),
    coverage2 = readr::col_integer(),
    confidence = readr::col_character(),
    closest_genomic_breakpoint1 = readr::col_character(),
    closest_genomic_breakpoint2 = readr::col_character(),
    filters = readr::col_character(),
    fusion_transcript = readr::col_character(),
    reading_frame = readr::col_character(),
    peptide_sequence = readr::col_character(),
    read_identifiers = readr::col_character()
  )
  
  arriba_calls<-read_tsv(arriba_calls, col_types = col_types)
  return(arriba_calls)
}

#' Read in fusion calls from STAR-Fusion 1.5.0
#'
#' @param starfusion_calls Please refer to software documenation [star-fusion.fusion_predictions.tsv] https://github.com/STAR-Fusion/STAR-Fusion/wiki#output-from-star-fusion
#'
#' @return A data.frame object with correct column specifications
#' @export
#'
#' @examples
#' fusionfileStarFusion <- read_starfusion_calls(
#'                            system.file("extdata", "starfusion_example.tsv", package = "annoFuse"))
#'                            
read_starfusion_calls <- function(starfusion_calls) {
  # set col types
  col_types = readr::cols(
    `#FusionName` = readr::col_character(),
    JunctionReadCount = readr::col_integer(),
    SpanningFragCount = readr::col_integer(),
    SpliceType = readr::col_character(),
    LeftGene = readr::col_character(),
    LeftBreakpoint = readr::col_character(),
    RightGene = readr::col_character(),
    RightBreakpoint = readr::col_character(),
    LargeAnchorSupport = readr::col_character(),
    FFPM = readr::col_double(),
    LeftBreakDinuc = readr::col_character(),
    LeftBreakEntropy = readr::col_double(),
    RightBreakDinuc = readr::col_character(),
    RightBreakEntropy = readr::col_double(),
    annots = readr::col_character(),
    CDS_LEFT_ID = readr::col_character(),
    CDS_LEFT_RANGE = readr::col_character(),
    CDS_RIGHT_ID = readr::col_character(),
    CDS_RIGHT_RANGE = readr::col_character(),
    PROT_FUSION_TYPE = readr::col_character(),
    FUSION_MODEL = readr::col_character(),
    FUSION_CDS = readr::col_character(),
    FUSION_TRANSL = readr::col_character(),
    PFAM_LEFT = readr::col_character(),
    PFAM_RIGHT = readr::col_character()
  )
  
  starfusion_calls<- read_tsv(starfusion_calls,col_types = col_types)
  return(starfusion_calls)
}

