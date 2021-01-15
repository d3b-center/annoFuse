#' Function to annotate fusion calls

#' @param standardFusioncalls A dataframe from star fusion or arriba standardized to run through the filtering steps
#' @param geneListReferenceDataTab A dataframe with column 1 as GeneName 2 source file 3 type; collapse to summarize type
#' @param fusionReferenceDataTab A dataframe with column 1 as FusionName 2 source file 3 type; collapse to summarize type
#' @param checkReciprocal Logical value to check if fusion also has reciprocal fusion in Sample, default to TRUE
#' 
#' @export
#'
#' @return Standardized fusion calls annotated with gene list and fusion list provided in reference folder. If checkReciprocal ==TRUE reciprocal status of fusion is also provided . 
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
#'filteredFusionAnnotated <- annotate_fusion_calls(
#'  standardFusioncalls = fusionQCFiltered,
#'  geneListReferenceDataTab = geneListReferenceDataTab, 
#'  fusionReferenceDataTab = fusionReferenceDataTab,
#'  checkReciprocal = TRUE)                                   
annotate_fusion_calls <- function(standardFusioncalls,
                                  geneListReferenceDataTab,
                                  fusionReferenceDataTab,
                                  checkReciprocal=TRUE) {
  standardFusioncalls <- .check_annoFuse_calls(standardFusioncalls)

  stopifnot(is(geneListReferenceDataTab, "data.frame"))
  stopifnot(is(fusionReferenceDataTab, "data.frame"))
  stopifnot(is.logical(checkReciprocal))
  
  if (checkReciprocal & nrow(standardFusioncalls)>0){
    # check for fusions have reciprocal fusions in the same Sample
    # works only for GeneY -- GeneX ; GeneX -- GeneY matches
    recirpocal_fusion <- function(FusionName,Sample,standardFusioncalls ){
      Gene1A <- strsplit(FusionName,"--")[[1]][1]
      Gene1B <- strsplit(FusionName,"--")[[1]][2]
      reciprocal <- paste0(Gene1B,"--",Gene1A)
      check <- any(standardFusioncalls$FusionName[standardFusioncalls$Sample==Sample] == reciprocal)
      df <- data.frame("FusionName"=FusionName,"Sample"=Sample,"reciprocal_exists"=check)
    }
    
    # run reciprocal_fusion function to get status of fusion if reciprocal
    is_reciprocal <- apply(standardFusioncalls,1,function(x) recirpocal_fusion(x["FusionName"],x["Sample"],standardFusioncalls))
    
    # convert list to dataframe
    is_reciprocal<-data.frame(Reduce(rbind, is_reciprocal))
    
    # merge to standardized fusion calls
    standardFusioncalls <- standardFusioncalls %>%
      dplyr::left_join(is_reciprocal,by=c("Sample","FusionName"))
    
  } else if(checkReciprocal & nrow(standardFusioncalls)==0){
    standardFusioncalls$reciprocal_exists <- logical()
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
