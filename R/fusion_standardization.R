#' Standardizes fusion calls
#'
#' Various fusion callers have different formats that make aggregating and filtering data difficult.
#' By standardizing fusion callers output we capture the required columns which we use for downstream
#' analysis
#'
#' @param fusion_calls A dataframe from star fusion or arriba (more callers to be added)
#' @param caller string options STARFUSION/ARRIBA
#' @param tumorID string or character vector of same length as fusion_calls
#'
#' @export
#'
#' @return Standardized fusion calls ready for filtering
#'
#' @examples
#' # read in arriba fusion file
#' fusionfileArriba <- read_arriba_calls(
#'   system.file("extdata", "arriba_example.tsv", package = "annoFuse"))
#' # read in starfusion file
#' fusionfileStarFusion <- read_starfusion_calls(
#'   system.file("extdata", "starfusion_example.tsv", package = "annoFuse"))
#' formattedArriba <- fusion_standardization(fusionfileArriba,
#'                                           caller = "ARRIBA",
#'                                           tumorID = "tumorID")
#' formattedStarFusion <- fusion_standardization(fusionfileStarFusion,
#'                                               caller = "STARFUSION",
#'                                               tumorID = "tumorID")
fusion_standardization <- function(fusion_calls,
                                   caller = c("STARFUSION", "ARRIBA"),
                                   tumorID = "tumorID") {
  stopifnot(is(fusion_calls, "data.frame"))
  stopifnot(is.character(caller))
  stopifnot(is.character(tumorID))

  # caller <- match.arg(caller, choices = c("STARFUSION", "ARRIBA"))
  
  if (caller == "STARFUSION") {
    fusion_calls <- fusion_calls %>%
      # standardize fusion type column name
      dplyr::rename(Fusion_Type = PROT_FUSION_TYPE,
                    FusionName = "#FusionName") %>%
      dplyr::mutate(
        # remove chr notation from breakpoint columns
        LeftBreakpoint = gsub("^chr", "", .data$LeftBreakpoint),
        RightBreakpoint = gsub("^chr", "", .data$RightBreakpoint),
        # remove strand information to match breakpoint locations
        LeftBreakpoint = gsub(":[-|+]$", "", .data$LeftBreakpoint),
        RightBreakpoint = gsub(":[-|+]$", "", .data$RightBreakpoint),
        # STARFusion does not return confidence information
        Confidence = NA,
        # standardize fusion types
        Fusion_Type = dplyr::case_when(
          Fusion_Type == "INFRAME" ~ "in-frame",
          Fusion_Type == "FRAMESHIFT" ~ "frameshift",
          TRUE ~ "other"
        ),
        Sample = tumorID,
        Caller = "STARFUSION"
      )
  }
  else if (caller == "ARRIBA") {
    fusion_calls <- fusion_calls %>%
      # standardizing fusion type annotation
      dplyr::rename(
        Fusion_Type = .data$reading_frame,
        Confidence = .data$confidence,
        # SpanningFragCount is equivalent to discordant_mates in Arriba
        SpanningFragCount = .data$discordant_mates
      ) %>%
      dplyr::mutate(
        LeftBreakpoint = gsub("^chr", "", .data$breakpoint1),
        RightBreakpoint = gsub("^chr", "", .data$breakpoint2),
        # readthrough information from arriba
        annots = dplyr::case_when(
          !is_empty(.data$annots) ~ paste(.data$annots, .data$type, sep = ","),
          is_empty(.data$annots) ~ .data$type),
        # Intergenic gene fusion breakpoints in arriba are annotated as
        # "gene1A,gene2A". As comma is used as a common delimiter in files changing
        # it to "/"
        FusionName = paste0(gsub(",", "/", .data$`#gene1`), "--", gsub(",", "/", .data$gene2)),
        # JunctionReadCount is equivalent to split reads in Arriba. Arriba however
        # provides split_reads1 and split_reads2 to provide information of reads
        # anchoring in gene1 or gene2
        JunctionReadCount = .data$split_reads1 + .data$split_reads2,
        Fusion_Type = dplyr::case_when(
          !Fusion_Type %in% c("out-of-frame", "in-frame") ~ "other",
          Fusion_Type == "out-of-frame" ~ "frameshift",
          TRUE ~ "in-frame"
        ),
        Sample = tumorID,
        Caller = "ARRIBA"
      )
  } else {
    stop(paste(caller, "is not a supported caller string."))
  }

  # Get standard columns for filtering
  
  standard_calls <- fusion_calls %>%
    # select columns for standarda fusion format
    dplyr::select(c("LeftBreakpoint",
                    "RightBreakpoint",
                    "FusionName",
                    "Sample",
                    "Caller",
                    "Fusion_Type",
                    "JunctionReadCount",
                    "SpanningFragCount",
                    "Confidence",
                    "annots")) %>%
    # to obtain geneA and geneB for gene search below
    bind_cols(reshape2::colsplit(fusion_calls$FusionName, pattern = "--", names = c("GeneA", "GeneB"))) %>%
    # Intergenic fusion will have Gene1A,Gene2A,Gene1B,Gene2B
    separate(.data$GeneA, sep = "/", into = c("Gene1A", "Gene2A"), remove = FALSE) %>%
    separate(.data$GeneB, sep = "/", into = c("Gene1B", "Gene2B"), remove = FALSE) %>%
    # remove distance to fusion breakpoint from gene names in intergenic fusion
    mutate(
      Gene1A = gsub("[(].*", "", .data$Gene1A),
      Gene2A = gsub("[(].*", "", .data$Gene2A),
      Gene1B = gsub("[(].*", "", .data$Gene1B),
      Gene2B = gsub("[(].*", "", .data$Gene2B),
      BreakpointLocation = case_when(
        Gene1A==Gene1B ~ "Intragenic",
        grepl("/",FusionName) ~ "Intergenic",
        TRUE ~ "Genic"),
      SpanningDelta = SpanningFragCount - JunctionReadCount
    ) %>%
    as.data.frame()

    return(standard_calls)
}
