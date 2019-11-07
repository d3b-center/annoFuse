#' Filters standardized fusion calls to remove artifacts and false positives.
#'
#' Events such as polymerase read-throughs, mis-mapping due to gene homology, and fusions occurring in healthy normal
#' tissue require stringent filtering, making it difficult for researchers and clinicians to discern true underlying
#' oncogenic drivers of a tumor and in some cases, appropriate therapy
#
#' @param standardFusioncalls A dataframe from star fusion or arriba standardized to run through the filtering steps
#' @param expressionMatrix A rds with expression data (RSEM/TPM counts) for the same cohort/GTEX
#' @param readingFramFilter A regex to capture readingframe (eg. in-frame|frameshift|other)
#' @param artifactFilter A red flag filter from Annotation ; in OpenPBTA annotation is from FusionAnnotator column: "annots"
#' @param junctionReadCountFilter An integer threshold for JunctionReadCount
#' @param spanningFragCountFilter An integer threshold for (SpanningFragCount - JunctionReadCount)
#' @param expressionFilter An integer threshold for TPM/FPKM filter
#' @return Standardized fusion calls filtered to pass QC and remove calls with insufficient read-support and annotation red-flags

############# artifactFiltering #############


fusion_filtering_QC<-function(standardFusioncalls=standardFusioncalls,readingFrameFilter=readingFrameFilter,artifactFilter=artifactFilter,junctionReadCountFilter=junctionReadCountFilter,spanningFragCountFilter=spanningFragCountFilter,reathroughFilter=TRUE){

  if( reathroughFilter & any(grepl('read.*through|NEIGHBORS',standardFusioncalls$annots,ignore.case = TRUE))){
    # Gather read throughs from standardized fusion calls
    rts <- standardFusioncalls[grep('read.*through|NEIGHBORS',standardFusioncalls$annots,ignore.case = TRUE),"FusionName"]
    # Reverse of read throughs to capture
    rts.rev <- unique(unlist(lapply(strsplit(rts, '--'), FUN = function(x) paste0(x[2],'--',x[1]))))
    # Combine read through and reverse fusion genes
    rts <- unique(c(rts, rts.rev))
    standardFusioncalls <- standardFusioncalls[-which(standardFusioncalls$FusionName %in% rts),]
  }

  if( !missing(readingFrameFilter ) ){
    # Error handling
    readingFrameTypes<-c("in-frame","frameshift","other")
    standardFusioncalls <- standardFusioncalls[grep(readingFrameFilter,standardFusioncalls$Fusion_Type),]
    if(!all(readingFrameTypes %in% standardFusioncalls$Fusion_Type)){
      warning(paste("No fusion calls with readingframe:",readingFrameTypes[-which(readingFrameTypes %in% standardFusioncalls$Fusion_Type)]))
    }
  }

  if( !missing(artifactFilter) & any(grepl(artifactFilter,standardFusioncalls$annots)) ) {
    # Error handling
    artifactFilterTypes<-unlist(strsplit(artifactFilter,"|",fixed=TRUE))
    if(any(unlist(lapply(artifactFilterTypes, function(x) !any(str_detect(as.character(standardFusioncalls$annots),x)))))){
      artifactFilterTypesNotFound<-artifactFilterTypes[unlist(lapply(artifactFilterTypes, function(x) !any(str_detect(as.character(standardFusioncalls$annots),x))))]
      warning(paste("No fusion calls with annotation:",artifactFilterTypesNotFound))
    }
    standardFusioncalls <- standardFusioncalls[-grep(artifactFilter,standardFusioncalls$annots),]
  }

  if( !missing(junctionReadCountFilter ) ){
    standardFusioncalls <- standardFusioncalls[which(standardFusioncalls$JunctionReadCount>=junctionReadCountFilter),]
  }

  if( !missing(spanningFragCountFilter ) ){
    # false positive calls either at the start or end of a transcript will be supported by an uneven majority of support from spanning fragments compared to junction reads
    # to remove these calls we are implementing this condition below
    standardFusioncalls <- standardFusioncalls[which( (standardFusioncalls$SpanningFragCount - standardFusioncalls$JunctionReadCount ) <= spanningFragCountFilter),]
  }

  return(standardFusioncalls)
}


############################################



