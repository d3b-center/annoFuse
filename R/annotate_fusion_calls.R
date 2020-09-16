#' Function to annotate fusion calls

#' @param standardFusioncalls A dataframe from star fusion or arriba standardized to run through the filtering steps
#' @param geneListReferenceDataTab A dataframe with column 1 as GeneName 2 source file 3 type; collapse to summarize type
#' @param fusionReferenceDataTab A dataframe with column 1 as FusionName 2 source file 3 type; collapse to summarize type
#' @param checkReciprocal Logical value to check if fusion also has reciprocal fusion in Sample, default to FALSE
#' @param checkDomainStatus Logical value to check if domain status in fused gene for given domansToCheck, default to FALSE,
#' @param domainsToCheck pfamID to check for retention status, the IDs can be found here http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/pfamDesc.txt.gz
#'
#' @export
#'
#' @return Standardized fusion calls annotated with gene list and fusion list provided in reference folder. If checkReciprocal ==TRUE reciprocal status of fusion is also provided . If checkDomainStatus =TRUE and domain retention status for given pfamID is also provided along with the gene location corresponding to the domain retention status
#'
#' @examples
#' # standardize
#' fusionfileArriba <- read_arriba_calls(
#'   system.file("extdata", "arriba_example.tsv", package = "annoFuse"))
#' fusionfileStarFusion <- read_starfusion_calls(
#'   system.file("extdata", "starfusion_example.tsv", package = "annoFuse"))
#' formattedArriba <- fusion_standardization(fusionfileArriba,
#'                                           caller="ARRIBA", 
#'                                           tumorID = "tumorID")
#' formattedStarFusion <- fusion_standardization(fusionfileStarFusion, 
#'                                               caller="STARFUSION",
#'                                               tumorID = "tumorID")
#' # merge standardized fusion calls
#' standardFusioncalls <- as.data.frame(rbind(formattedStarFusion, formattedArriba))
#' fusionQCFiltered <- fusion_filtering_QC(standardFusioncalls = standardFusioncalls, 
#'                     readingFrameFilter = "in-frame|frameshift|other",
#'                     artifactFilter = "GTEx_Recurrent|DGD_PARALOGS|Normal|BodyMap|ConjoinG",
#'                     junctionReadCountFilter = 1,
#'                     spanningFragCountFilter = 10,
#'                     readthroughFilter = TRUE)
#' # annotated from gene and fusion refrence list
#' # read in gene and fusion reference tab
#' geneListReferenceDataTab <- read.delim(
#'   system.file("extdata", "genelistreference.txt", package = "annoFuse"), stringsAsFactors = FALSE)
#'# column 1 as FusionName 2 source file 3 type; collapse to summarize type
#'fusionReferenceDataTab <- read.delim(
#'  system.file("extdata", "fusionreference.txt", package = "annoFuse"), stringsAsFactors = FALSE)
#'
#'# to check status of "Pkinase" domain per fused gene
#' kinaseid <-"PF00069"
#'
#'filteredFusionAnnotated <- annotate_fusion_calls(
#'  standardFusioncalls = fusionQCFiltered,
#'  geneListReferenceDataTab = geneListReferenceDataTab, 
#'  fusionReferenceDataTab = fusionReferenceDataTab,
#'  checkReciprocal = TRUE,checkDomainStatus = TRUE,
#'  domainsToCheck=kinaseid)                                   
annotate_fusion_calls <- function(standardFusioncalls,
                                  geneListReferenceDataTab,
                                  fusionReferenceDataTab,
                                  checkReciprocal=FALSE,
                                  checkDomainStatus=FALSE,
                                  domainsToCheck) {
  standardFusioncalls <- .check_annoFuse_calls(standardFusioncalls)

  stopifnot(is(geneListReferenceDataTab, "data.frame"))
  stopifnot(is(fusionReferenceDataTab, "data.frame"))
  stopifnot(is.logical(checkReciprocal))
  stopifnot(is.logical(checkDomainStatus))
  
  if(checkDomainStatus){
    stopifnot(is.character(domainsToCheck))
    
    # load bioMart Pfam dataframe
    bioMartDataPfam <- readRDS(system.file("extdata", "pfamDataBioMart.RDS", package = "annoFuse"))
    bioMartDataPfam <- bioMartDataPfam %>% 
      # keep only pfamIDs to check
      dplyr::filter(pfam_id %in% domainsToCheck)
    
    # annotate by domain ;gather partial retention as well
    annDomain<-get_Pfam_domain(standardFusioncalls  = standardFusioncalls,bioMartDataPfam = bioMartDataPfam,keepPartialAnno = TRUE)
    
    # domain annotation
    standardFusionGene1ADomain<-annDomain$Gene1A %>% 
      mutate(GeneLoc="Gene1A") %>% 
      dplyr::rename("DomainRetained"="Gene1A_DOMAIN_RETAINED_IN_FUSION")
    standardFusionGene1BDomain<-annDomain$Gene1B %>% 
      mutate(GeneLoc="Gene1B") %>% 
      dplyr::rename("DomainRetained"="Gene1B_DOMAIN_RETAINED_IN_FUSION")
    
    standardFusionDomainStatus <- unique(rbind(standardFusionGene1ADomain,standardFusionGene1BDomain)) %>% 
      # mutated Left and Right Breakpoints since the get_Pfam_domain separtaes the 
      # chromosome and genomic location into LeftBreakpointChr,LeftBreakpoint
      # and RightBreakpointChr,RightBreakpoint
      # but since format for standardFusioncalls is Chr:Breakpont we are updating here
      # for merging in the next step
      dplyr::mutate(LeftBreakpoint=paste(LeftBreakpointChr,LeftBreakpoint,sep=":"),
                    RightBreakpoint=paste(RightBreakpointChr,RightBreakpoint,sep = ":")) %>%
      # select only columns required
      dplyr::select("LeftBreakpoint","RightBreakpoint","FusionName","Sample","DomainRetained","GeneLoc") %>%
      unique()
    
    standardFusioncalls <- standardFusioncalls  %>%
      # add status for given pfamIDs
      full_join(standardFusionDomainStatus) %>%
      # select only columns required
      dplyr::select("Sample","LeftBreakpoint","RightBreakpoint",
                    "FusionName","Caller","Fusion_Type","JunctionReadCount",
                    "SpanningFragCount","Confidence","annots","Gene1A",
                    "Gene2A","Gene1B","Gene2B", "BreakpointLocation",
                    "SpanningDelta","DomainRetained","GeneLoc") 
  } 
  
  if (checkReciprocal){
    # check for fusions have reciprocal fusions in the same Sample
    # works only for GeneY -- GeneX ; GeneX -- GeneY matches
    recirpocal_fusion <- function(FusionName,Sample,standardFusioncalls ){
      Gene1A <- strsplit(FusionName,"--")[[1]][1]
      Gene1B <- strsplit(FusionName,"--")[[1]][2]
      reciprocal <- paste0(Gene1B,"--",Gene1A)
      check <- any(grepl(reciprocal,standardFusioncalls$FusionName[standardFusioncalls$Sample==Sample]))
      df <- data.frame("FusionName"=FusionName,"Sample"=Sample,"reciprocal_exists"=check)
    }
    
    # run reciprocal_fusion function to get status of fusion if reciprocal
    is_reciprocal <- apply(standardFusioncalls,1,function(x) recirpocal_fusion(x["FusionName"],x["Sample"],standardFusioncalls))
    
    # convert list to dataframe
    is_reciprocal<-data.frame(Reduce(rbind, is_reciprocal))
    
    # merge to standardized fusion calls
    standardFusioncalls <- standardFusioncalls %>%
      dplyr::left_join(is_reciprocal,by=c("Sample","FusionName"))
    
  }
  
  geneListReferenceDataTab <- geneListReferenceDataTab %>% dplyr::select(-file)
  fusionReferenceDataTab <- fusionReferenceDataTab %>% dplyr::select(-file)
  annotated_filtered_fusions <- standardFusioncalls %>%
    # annotate Gene1A
    dplyr::left_join(geneListReferenceDataTab, by = c("Gene1A" = "Gene_Symbol")) %>%
    dplyr::rename(Gene1A_anno = .data$type) %>%
    # annotate Gene1B
    dplyr::left_join(geneListReferenceDataTab, by = c("Gene1B" = "Gene_Symbol")) %>%
    dplyr::rename(Gene1B_anno = .data$type) %>%
    # annotate Gene2A
    dplyr::left_join(geneListReferenceDataTab, by = c("Gene2A" = "Gene_Symbol")) %>%
    dplyr::rename(Gene2A_anno = .data$type) %>%
    # annotate Gene2B
    dplyr::left_join(geneListReferenceDataTab, by = c("Gene2B" = "Gene_Symbol")) %>%
    dplyr::rename(Gene2B_anno = .data$type) %>%
    # annotate FusionName
    dplyr::left_join(fusionReferenceDataTab, by = c("FusionName" = "FusionName")) %>%
    dplyr::rename(Fusion_anno = .data$type) %>%
    as.data.frame()
  annotated_filtered_fusions <- unique(annotated_filtered_fusions)

  return(annotated_filtered_fusions)
}
