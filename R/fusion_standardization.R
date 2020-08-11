#' Standardizes fusion calls
#'
#' Various fusion callers have different formats that make aggregating and filtering data difficult.
#' By standardizing fusion callers output we capture the required columns which we use for downstream
#' analysis
#'
#' @param fusion_calls A dataframe from star fusion or arriba (more callers to be added)
#' @param caller string options STARfusion/arriba
#'
#' @export
#'
#' @return Standardized fusion calls ready for filtering

fusion_standardization <- function(fusion_calls = fusion_calls,
                                   caller = caller) {
  if (caller == "STARFUSION") {
    fusion_calls <- fusion_calls %>%
      # standardize fusion type column name
      dplyr::rename(Fusion_Type = .data$PROT_FUSION_TYPE) %>%
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
        )
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
        annots = paste(.data$annots, .data$type, sep = ","),
        # Intergenic gene fusion breakpoints in arriba are annotated as
        # "gene1A,gene2A". As comma is used as a common delimiter in files changing
        # it to "/"
        FusionName = paste0(gsub(",", "/", .data$gene1), "--", gsub(",", "/", .data$gene2)),
        # JunctionReadCount is equivalent to split reads in Arriba. Arriba however
        # provides split_reads1 and split_reads2 to provide information of reads
        # anchoring in gene1 or gene2
        JunctionReadCount = .data$split_reads1 + .data$split_reads2,
        Fusion_Type = dplyr::case_when(
          !Fusion_Type %in% c("out-of-frame", "in-frame") ~ "other",
          Fusion_Type == "out-of-frame" ~ "frameshift",
          TRUE ~ "in-frame"
        )
      )
  } else {
    stop(paste(caller, "is not a supported caller string."))
  }

  # Get standard columns for filtering
  standard_calls <- unique(fusion_calls[, c("LeftBreakpoint", "RightBreakpoint", "FusionName", "Sample", "Caller", "Fusion_Type", "JunctionReadCount", "SpanningFragCount", "Confidence", "annots")])

  return(standard_calls)
}
