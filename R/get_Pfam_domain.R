#' Function to annotate fusion calls with domain terms

#' @param standardFusioncalls A dataframe from star fusion or arriba standardized to run through the filtering steps
#' @param bioMartDataPfam A dataframe with gene and domain coordinate information with chromosome_name,gene_start,gene_end,domain_chr,domain_start,domain_end,hgnc_symbol
#' @param keepPartialAnno TRUE or FALSE to keep partial status; defaults to FALSE
#'
#' @export
#'
#' @return Standardized fusion calls as list Gene1A and Gene1B annotated with domain terms and chromosome location; retained and not retained,optionally partially retained
#'
#' @examples
#' out_annofuse <- system.file("extdata", "PutativeDriverAnnoFuse.tsv", package = "annoFuseData")
#' sfc <- read.delim(out_annofuse)
#' bioMartDataPfam <- readRDS(system.file("extdata", "pfamDataBioMart.RDS", package = "annoFuseData"))
#' domain_list_df <- get_Pfam_domain(standardFusioncalls = sfc, bioMartDataPfam = bioMartDataPfam)
get_Pfam_domain <- function(standardFusioncalls,
                            bioMartDataPfam,
                            keepPartialAnno = FALSE) {
  standardFusioncalls <- .check_annoFuse_calls(standardFusioncalls)
  stopifnot(is(bioMartDataPfam, "data.frame"))
  stopifnot(is.logical(keepPartialAnno))

  # get loci and chromosome as different columns
  standardFusioncalls$RightBreakpointChr <- gsub(":.*$", "", standardFusioncalls$RightBreakpoint)
  standardFusioncalls$RightBreakpoint <- gsub("^.*:", "", standardFusioncalls$RightBreakpoint)
  standardFusioncalls$LeftBreakpointChr <- gsub(":.*$", "", standardFusioncalls$LeftBreakpoint)
  standardFusioncalls$LeftBreakpoint <- gsub("^.*:", "", standardFusioncalls$LeftBreakpoint)


  # merge fusion and gene + domain information
  standardFusioncallsGene1A <- standardFusioncalls %>%
    left_join(bioMartDataPfam, by = c("Gene1A" = "hgnc_symbol")) %>%
    as.data.frame()
  standardFusioncallsGene1B <- standardFusioncalls %>%
    left_join(bioMartDataPfam, by = c("Gene1B" = "hgnc_symbol")) %>%
    as.data.frame()

  # positive strand gene1A annotation
  standardFusioncallsGene1APosStrand <- standardFusioncallsGene1A[which(standardFusioncallsGene1A$strand == "1"), ]
  standardFusioncallsGene1APosStrand[, "Gene1A_DOMAIN_RETAINED_IN_FUSION"] <- ifelse((
    # domain starts after gene start
    standardFusioncallsGene1APosStrand$domain_start > standardFusioncallsGene1APosStrand$gene_start &
      # domain ends before the breakpoint
      standardFusioncallsGene1APosStrand$domain_end < standardFusioncallsGene1APosStrand$LeftBreakpoint) &
    standardFusioncallsGene1APosStrand$chromosome_name == standardFusioncallsGene1APosStrand$domain_chr,
  # we dont check from Fusin_Type for Gene1A since the frameshifts only affect Gene1B translation
  "Yes", "No"
  )

  # negative strand gene1A annotation
  standardFusioncallsGene1ANegStrand <- standardFusioncallsGene1A[which(standardFusioncallsGene1A$strand == "-1"), ]
  standardFusioncallsGene1ANegStrand[, "Gene1A_DOMAIN_RETAINED_IN_FUSION"] <- ifelse((
    # domain starts after breakpoint
    standardFusioncallsGene1ANegStrand$domain_start > standardFusioncallsGene1ANegStrand$LeftBreakpoint &
      # domain ends before gene end
      standardFusioncallsGene1ANegStrand$domain_end < standardFusioncallsGene1ANegStrand$gene_end) &

    standardFusioncallsGene1ANegStrand$chromosome_name == standardFusioncallsGene1ANegStrand$domain_chr,
  # we dont check from Fusin_Type for Gene1A since the frameshifts only affect Gene1B translation
  "Yes", "No"
  )


  # positive strand gene1B annotation
  standardFusioncallsGene1BPosStrand <- standardFusioncallsGene1B[which(standardFusioncallsGene1B$strand == "1"), ]
  standardFusioncallsGene1BPosStrand[, "Gene1B_DOMAIN_RETAINED_IN_FUSION"] <- ifelse((
    # domain starts after breakpoint
    standardFusioncallsGene1BPosStrand$domain_start > standardFusioncallsGene1BPosStrand$RightBreakpoint &
      # domain ends before gene end
      standardFusioncallsGene1BPosStrand$domain_end < standardFusioncallsGene1BPosStrand$gene_end) &
    standardFusioncallsGene1BPosStrand$chromosome_name == standardFusioncallsGene1BPosStrand$domain_chr &
    # we check from Fusin_Type since the frameshifts affect Gene1B translation
    standardFusioncallsGene1BPosStrand$Fusion_Type == "in-frame", "Yes", "No")

  # negative strand gene1B annotation
  standardFusioncallsGene1BNegStrand <- standardFusioncallsGene1B[which(standardFusioncallsGene1B$strand == "-1"), ]
  standardFusioncallsGene1BNegStrand[, "Gene1B_DOMAIN_RETAINED_IN_FUSION"] <- ifelse((
    # domain end is before breakpoint
    standardFusioncallsGene1BNegStrand$domain_end < standardFusioncallsGene1BNegStrand$RightBreakpoint &
      # domain starts after gene start
      standardFusioncallsGene1BNegStrand$domain_start > standardFusioncallsGene1BNegStrand$gene_start) &
    standardFusioncallsGene1BNegStrand$chromosome_name == standardFusioncallsGene1BNegStrand$domain_chr &
    # we check from Fusin_Type since the frameshifts affect Gene1B translation
    standardFusioncallsGene1BNegStrand$Fusion_Type == "in-frame", "Yes", "No")



  # rbind negative and positive
  standardFusioncallsGene1AMapped <- rbind(standardFusioncallsGene1APosStrand, standardFusioncallsGene1ANegStrand)
  standardFusioncallsGene1BMapped <- rbind(standardFusioncallsGene1BPosStrand, standardFusioncallsGene1BNegStrand)

  # keep partial status of domain location overlap
  if (keepPartialAnno) {
    # domain mapped and retained
    standardFusioncallsGene1AMappedRetained <- standardFusioncallsGene1AMapped[grep("Yes", standardFusioncallsGene1AMapped$Gene1A_DOMAIN_RETAINED_IN_FUSION), ]
    # domain mapped and lost or not domain found
    standardFusioncallsGene1AMappedLost <- standardFusioncallsGene1AMapped[which(standardFusioncallsGene1AMapped$Gene1A_DOMAIN_RETAINED_IN_FUSION == "No" | is.na(standardFusioncallsGene1AMapped$Gene1A_DOMAIN_RETAINED_IN_FUSION)), ]

    # add partial status for in-frame fusions
    standardFusioncallsGene1AMappedLost[, "Gene1A_DOMAIN_RETAINED_IN_FUSION"] <- ifelse((
      # breakpoint is within domains
      standardFusioncallsGene1AMappedLost$domain_start < standardFusioncallsGene1AMappedLost$LeftBreakpoint &
        standardFusioncallsGene1AMappedLost$domain_end > standardFusioncallsGene1AMappedLost$LeftBreakpoint) &
      standardFusioncallsGene1AMappedLost$chromosome_name == standardFusioncallsGene1AMappedLost$domain_chr &
      # and in-frame
      standardFusioncallsGene1AMappedLost$Fusion_Type == "in-frame", "Partial", standardFusioncallsGene1AMappedLost$Gene1A_DOMAIN_RETAINED_IN_FUSION)

    # domain mapped and retained
    standardFusioncallsGene1BMappedRetained <- standardFusioncallsGene1BMapped[grep("Yes", standardFusioncallsGene1BMapped$Gene1B_DOMAIN_RETAINED_IN_FUSION), ]
    # domain mapped and lost or not domain found
    standardFusioncallsGene1BMappedLost <- standardFusioncallsGene1BMapped[which(standardFusioncallsGene1BMapped$Gene1B_DOMAIN_RETAINED_IN_FUSION == "No" | is.na(standardFusioncallsGene1BMapped$Gene1B_DOMAIN_RETAINED_IN_FUSION)), ]
    # add partial status for in-frame fusions
    standardFusioncallsGene1BMappedLost[, "Gene1B_DOMAIN_RETAINED_IN_FUSION"] <- ifelse((
      # breakpoint is within domains
      standardFusioncallsGene1BMappedLost$domain_start < standardFusioncallsGene1BMappedLost$RightBreakpoint &
        standardFusioncallsGene1BMappedLost$domain_end > standardFusioncallsGene1BMappedLost$RightBreakpoint) &
      standardFusioncallsGene1BMappedLost$chromosome_name == standardFusioncallsGene1BMappedLost$domain_chr &
      # and in-frame
      standardFusioncallsGene1BMappedLost$Fusion_Type == "in-frame", "Partial", standardFusioncallsGene1BMappedLost$Gene1B_DOMAIN_RETAINED_IN_FUSION)

    standardFusioncallsGene1AMapped <- rbind(standardFusioncallsGene1AMappedRetained, standardFusioncallsGene1AMappedLost)
    standardFusioncallsGene1BMapped <- rbind(standardFusioncallsGene1BMappedRetained, standardFusioncallsGene1BMappedLost)
  }

  # need to add back genes that were not mapped using bioMart dataframe
  if (nrow(standardFusioncallsGene1A[is.na(standardFusioncallsGene1A$strand), ]) > 0) {
    standardFusioncallsGene1ANotMapped <- standardFusioncallsGene1A[is.na(standardFusioncallsGene1A$strand), ]
    standardFusioncallsGene1ANotMapped[, "Gene1A_DOMAIN_RETAINED_IN_FUSION"] <- NA
    # merge mapped and unmapped df
    standardFusioncallsGene1A <- rbind(standardFusioncallsGene1AMapped, standardFusioncallsGene1ANotMapped)
  } else {
    standardFusioncallsGene1A <- standardFusioncallsGene1AMapped
  }

  if (nrow(standardFusioncallsGene1B[is.na(standardFusioncallsGene1B$strand), ]) > 0) {
    standardFusioncallsGene1BNotMapped <- standardFusioncallsGene1B[is.na(standardFusioncallsGene1B$strand), ]
    standardFusioncallsGene1BNotMapped[, "Gene1B_DOMAIN_RETAINED_IN_FUSION"] <- NA
    # merge mapped and unmapped df
    standardFusioncallsGene1B <- rbind(standardFusioncallsGene1BMapped, standardFusioncallsGene1BNotMapped)
  } else {
    standardFusioncallsGene1B <- standardFusioncallsGene1BMapped
  }




  domain <- list("Gene1A" = standardFusioncallsGene1A, "Gene1B" = standardFusioncallsGene1B)

  return(domain)
}
