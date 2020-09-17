#' Filter standardized fusion calls for driver fusions
#'
#' If standardized fusion calls are annotated using the geneListReferenceDataTab and fusionReferenceDataTab filters out fusion calls where partner genes are not annotated.
#' If standardized fusion is not annotated it will be annotated with geneListReferenceDataTab and fusionReferenceDataTab provided.
#' Domain retention status for Gene1A and Gene1B for the given pfamIDs is also annotated; defaults to kinase domain retention status information
#'
#' @param standardFusioncalls A dataframe from star fusion or arriba (more callers to be added)
#' @param filterPutativeDriver filter out fusion calls where partner genes are not annotated from with gene and fusion reference list by annnoFuse::annotate_fusion_calls()
#' @param annotated Logical value to specify if input if annotated by annnoFuse::annotate_fusion_calls()
#' @param geneListReferenceDataTab A dataframe with column 1 as GeneName 2 source file 3 type; collapse to summarize type
#' @param fusionReferenceDataTab A dataframe with column 1 as FusionName 2 source file 3 type; collapse to summarize type
#' @param checkDomainStatus Logical value to check if domain status in fused gene for given domansToCheck, default to FALSE
#' @param domainsToCheck pfamID to check for retention status, the IDs can be found here http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/pfamDesc.txt.gz, defaults to using kinase pfam IDs since the fusions with kinase domain are more relevant for therapy
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
#'   
#' bioMartDataPfam <- readRDS(system.file("extdata", "pfamDataBioMart.RDS", package = "annoFuse"))
#' kinaseid<-unique(bioMartDataPfam$pfam_id[grep("kinase",bioMartDataPfam$NAME)] )
#' fusion_driver_df <- fusion_driver(sfc,
#'                                   annotated = TRUE,
#'                                   geneListReferenceDataTab = geneListReferenceDataTab,
#'                                   fusionReferenceDataTab = fusionReferenceDataTab,
#'                                   checkDomainStatus=FALSE,
#'                                   domainsToCheck=kinaseid)
#' 
fusion_driver <- function(standardFusioncalls,
                          filterPutativeDriver=TRUE,
                          annotated = TRUE,
                          geneListReferenceDataTab,
                          fusionReferenceDataTab,
                          checkDomainStatus=FALSE,
                          domainsToCheck) {
  standardFusioncalls <- .check_annoFuse_calls(standardFusioncalls)
  stopifnot(is.logical(annotated))
  stopifnot(is.logical(checkDomainStatus))
   
  if(checkDomainStatus){
    # load bioMart Pfam dataframe
    bioMartDataPfam <- readRDS(system.file("extdata", "pfamDataBioMart.RDS", package = "annoFuse"))
    
    if(missing(domainsToCheck)){
      message("domainsToCheck was not provided; Using default kinase pfam ids to checkDomainStatus")
      domainsToCheck<-unique(bioMartDataPfam$pfam_id[grep("kinase",bioMartDataPfam$NAME)] )
      kinaseDomainRetained <- TRUE
    }
    
    if(!is.character(domainsToCheck) | !(any(domainsToCheck %in% bioMartDataPfam$pfam_id))){
      message("domainsToCheck was not character types or not found in pfam; Using default kinase pfam ids to checkDomainStatus")
      domainsToCheck<-unique(bioMartDataPfam$pfam_id[grep("kinase",bioMartDataPfam$NAME)] )
      kinaseDomainRetained <- TRUE
    }
    
    if(is.character(domainsToCheck)){
      found<-paste(intersect(domainsToCheck,bioMartDataPfam$pfam_id),collapse = ",")
      warning(paste(found," pfamIDs were found in our pfam list"))
    }
    
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
      dplyr::select("LeftBreakpoint","RightBreakpoint","FusionName","Fusion_Type","Gene1A","Sample","DomainRetainedInGene1A") %>%
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
      dplyr::select("LeftBreakpoint","RightBreakpoint","FusionName","Fusion_Type","Sample","Gene1B","DomainRetainedInGene1B") %>%
      unique()
    
    standardFusionDomain <- full_join(standardFusionGene1ADomain,standardFusionGene1BDomain )
    
    if (exists("kinaseDomainRetained") ){
      # default to always print kinase domain retention
      standardFusionDomain <- standardFusionDomain %>% 
        dplyr::rename("kinaseDomainRetainedGene1A"="DomainRetainedInGene1A",
                      "kinaseDomainRetainedGene1B"="DomainRetainedInGene1B")
    }
    
    standardFusioncalls <- standardFusioncalls %>%
      # add status for given pfamIDs
      left_join(standardFusionDomain) %>%
      unique()
  } 

  if (!annotated) {
    # check reference input
    stopifnot(is(geneListReferenceDataTab, "data.frame"))
    stopifnot(is(fusionReferenceDataTab, "data.frame"))
    
    #annotate
    standardFusioncalls <- annotate_fusion_calls(standardFusioncalls = standardFusioncalls, geneListReferenceDataTab = geneListReferenceDataTab, fusionReferenceDataTab = fusionReferenceDataTab)
  }
  
  if (filterPutativeDriver) {
    standardFusioncalls <- standardFusioncalls %>%
      #  dplyr::filter(!Gene1A %in% fusion_recurrent5_per_sample$GeneSymbol |
      #                  !Gene2A %in% fusion_recurrent5_per_sample$GeneSymbol |
      #                  !Gene1B %in% fusion_recurrent5_per_sample$GeneSymbol |
      #                  !Gene2B %in% fusion_recurrent5_per_sample$GeneSymbol) %>%
      dplyr::filter(!is.na(.data$Gene1A_anno) | !is.na(.data$Gene1B_anno) | !is.na(.data$Gene2A_anno) | !is.na(.data$Gene2B_anno))
  }

  return(standardFusioncalls)
}
