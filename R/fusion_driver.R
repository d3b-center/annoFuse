#' Filter standardized fusion calls for driver fusions
#'
#' If standardized fusion calls are annotated using the geneListReferenceDataTab and fusionReferenceDataTab filters out fusion calls where partner genes are not annotated.
#' If standardized fusion is not annotated it will be annotated with geneListReferenceDataTab and fusionReferenceDataTab provided
#'
#' @param standardFusioncalls A dataframe from star fusion or arriba (more callers to be added)
#' @param annotated Boolean if annotated
#' @param geneListReferenceDataTab A dataframe with column 1 as GeneName 2 source file 3 type; collapse to summarize type
#' @param fusionReferenceDataTab A dataframe with column 1 as FusionName 2 source file 3 type; collapse to summarize type
#' @param checkDomainStatus Logical value to check if domain status in fused gene for given domansToCheck, default to FALSE
#' @param domainsToCheck pfamID to check for retention status, the IDs can be found here http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/pfamDesc.txt.gz
#'
#'
#' @export
#'
#' @return Putative Driver standardized fusion calls annotated with gene list and fusion list provided in reference folder.If checkDomainStatus =TRUE and domain retention status for given pfamID is also provided along with the gene location corresponding to the domain retention status
#'
#' @examples
#' out_annofuse <- system.file("extdata", "PutativeDriverAnnoFuse.tsv", package = "annoFuse")
#' sfc <- read.delim(out_annofuse)
#' geneListReferenceDataTab <- read.delim(
#'   system.file("extdata", "genelistreference.txt", package="annoFuse"), stringsAsFactors = FALSE)
#' fusionReferenceDataTab <- read.delim(
#'   system.file("extdata", "fusionreference.txt", package="annoFuse"), stringsAsFactors = FALSE)
#' fusion_driver_df <- fusion_driver(sfc,
#'                                   annotated = TRUE,
#'                                   geneListReferenceDataTab = geneListReferenceDataTab,
#'                                   fusionReferenceDataTab = fusionReferenceDataTab)
#' 
fusion_driver <- function(standardFusioncalls,
                          annotated = TRUE,
                          geneListReferenceDataTab,
                          fusionReferenceDataTab,
                          checkDomainStatus=FALSE,
                          domainsToCheck) {
  standardFusioncalls <- .check_annoFuse_calls(standardFusioncalls)
  stopifnot(is.logical(annotated))
  stopifnot(is(geneListReferenceDataTab, "data.frame"))
  stopifnot(is(fusionReferenceDataTab, "data.frame"))
  
  if(checkDomainStatus){
    stopifnot(is.character(domainsToCheck))
    
    # load bioMart Pfam dataframe
    bioMartDataPfam <- readRDS(system.file("extdata", "pfamDataBioMart.RDS", package = "annoFuse"))
    bioMartDataPfam <- bioMartDataPfam %>% 
      # keep only pfamIDs to check
      dplyr::filter(pfam_id %in% domainsToCheck)
    
    # annotate by domain ;gather partial retention as well
    annDomain <- get_Pfam_domain(standardFusioncalls  = standardFusioncalls,bioMartDataPfam = bioMartDataPfam,keepPartialAnno = TRUE)
    
    # domain annotation
    standardFusionGene1ADomain<-annDomain$Gene1A %>% 
      dplyr::rename("DomainRetainedInGene1A"="Gene1A_DOMAIN_RETAINED_IN_FUSION") %>% 
      # mutated Left and Right Breakpoints since the get_Pfam_domain separtaes the 
      # chromosome and genomic location into LeftBreakpointChr,LeftBreakpoint
      # and RightBreakpointChr,RightBreakpoint
      # but since format for standardFusioncalls is Chr:Breakpont we are updating here
      # for merging in the next step
      dplyr::mutate(LeftBreakpoint=paste(LeftBreakpointChr,LeftBreakpoint,sep=":"),
                    RightBreakpoint=paste(RightBreakpointChr,RightBreakpoint,sep = ":")) %>%
      # select only columns required
      dplyr::select("LeftBreakpoint","RightBreakpoint","FusionName","Gene1A","Sample","DomainRetainedInGene1A") %>%
      unique()
    
    standardFusionGene1BDomain<-annDomain$Gene1B %>% 
      dplyr::rename("DomainRetainedInGene1B"="Gene1B_DOMAIN_RETAINED_IN_FUSION") %>% 
      # mutated Left and Right Breakpoints since the get_Pfam_domain separtaes the 
      # chromosome and genomic location into LeftBreakpointChr,LeftBreakpoint
      # and RightBreakpointChr,RightBreakpoint
      # but since format for standardFusioncalls is Chr:Breakpont we are updating here
      # for merging in the next step
      dplyr::mutate(LeftBreakpoint=paste(LeftBreakpointChr,LeftBreakpoint,sep=":"),
                    RightBreakpoint=paste(RightBreakpointChr,RightBreakpoint,sep = ":")) %>%
      # select only columns required
      dplyr::select("LeftBreakpoint","RightBreakpoint","FusionName","Sample","Gene1B","DomainRetainedInGene1B") %>%
      unique()
    
    standardFusionDomain <- full_join(standardFusionGene1ADomain,standardFusionGene1BDomain )
    
    standardFusioncalls <- standardFusioncalls %>%
      # add status for given pfamIDs
      left_join(standardFusionDomain) %>%
      unique()
  } 

  if (annotated) {
    putative_driver_fusions <- standardFusioncalls %>%
      #  dplyr::filter(!Gene1A %in% fusion_recurrent5_per_sample$GeneSymbol |
      #                  !Gene2A %in% fusion_recurrent5_per_sample$GeneSymbol |
      #                  !Gene1B %in% fusion_recurrent5_per_sample$GeneSymbol |
      #                  !Gene2B %in% fusion_recurrent5_per_sample$GeneSymbol) %>%
      dplyr::filter(!is.na(.data$Gene1A_anno) | !is.na(.data$Gene1B_anno) | !is.na(.data$Gene2A_anno) | !is.na(.data$Gene2B_anno))
  } else {
    standardFusioncalls <- annotate_fusion_calls(standardFusioncalls = standardFusioncalls, geneListReferenceDataTab = geneListReferenceDataTab, fusionReferenceDataTab = fusionReferenceDataTab)
    putative_driver_fusions <- standardFusioncalls %>%
      #    dplyr::filter(!Gene1A %in% fusion_recurrent5_per_sample$GeneSymbol |
      #                    !Gene2A %in% fusion_recurrent5_per_sample$GeneSymbol |
      #                    !Gene1B %in% fusion_recurrent5_per_sample$GeneSymbol |
      #                    !Gene2B %in% fusion_recurrent5_per_sample$GeneSymbol) %>%
      dplyr::filter(!is.na(.data$Gene1A_anno) | !is.na(.data$Gene1B_anno) | !is.na(.data$Gene2A_anno) | !is.na(.data$Gene2B_anno))
  }

  return(putative_driver_fusions)
}
