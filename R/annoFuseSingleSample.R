#' Single Sample use for annoFuse


#' Performs artifact filter to remove readthroughs,red flags and performs expression filtering with user provided expression Matrix and expression threshold

#' @param fusionfileArriba A dataframe from arriba fusion caller
#' @param fusionfileStarFusion A dataframe from starfusion caller
#' @param expressionFile Expression matrix for samples used in cohort for fusion calls
#' @param expressionFilter FPKM/TPM threshold for not expressed
#' @param tumorID Sample name to be used 
#' @param readingFrameFilter A regex to capture readingframe (eg. in-frame|frameshift|other)
#' @param readthroughFilter Boolean for filtering readthroughs
#' @param artifactFilter A red flag filter from Annotation ; in OpenPBTA annotation is from FusionAnnotator column "annots"
#' @param junctionReadCountFilter An integer threshold for JunctionReadCount
#' @param spanningFragCountFilter An integer threshold for (SpanningFragCount - JunctionReadCount)
#' @return Standardized fusion calls annotated with gene list and fusion list provided in reference folder

annoFuseSingleSample<-function(fusionfileArriba=fusionfileArriba,fusionfileStarFusion=fusionfileStarFusion,expressionFile=NULL,expressionFilter=1,tumorID="tumorID",artifactFilter="GTEx_Recurrent|DGD_PARALOGS|Normal|BodyMap|ConjoinG",readingFrameFilter="in-frame|frameshift|other",junctionReadCountFilter=1,spanningFragCountFilter=10,readthroughFilter=FALSE){
  
  # read files
  STARFusioninputfile<-read_tsv(fusionfileStarFusion)
  Arribainputfile<-read_tsv(fusionfileArriba)
  
  # read in gene and fusion reference tab
  geneListReferenceDataTab<-read.delim(system.file("extdata","genelistreference.txt", package="annoFuse"),stringsAsFactors = FALSE)
  geneListReferenceDataTab<-geneListReferenceDataTab %>% dplyr::group_by(Gene_Symbol) %>% dplyr::mutate(type = toString(type)) %>%
    dplyr::distinct(Gene_Symbol, type,file) %>% as.data.frame()
  
  # column 1 as FusionName 2 source file 3 type; collapse to summarize type
  fusionReferenceDataTab<-read.delim(system.file("extdata","fusionreference.txt", package="annoFuse"),stringsAsFactors = FALSE)
  fusionReferenceDataTab<-fusionReferenceDataTab %>%
    dplyr::distinct(FusionName,type,file) %>% as.data.frame()
  
  # if StarFusion and Arriba files empty execution stops
  if(is_empty(STARFusioninputfile$FusionName) & is_empty(Arribainputfile$"gene1--gene2")){
    stop("StarFusion and Arriba files empty")
  }
  
  # if StarFusion and Arriba files are not empty
  if(!is_empty(STARFusioninputfile$FusionName) & !is_empty(Arribainputfile$"gene1--gene2")){
    colnames(Arribainputfile)[27]<-"annots"
    
    # To have a general column with unique IDs associated with each sample
    STARFusioninputfile$Sample <- tumorID
    Arribainputfile$Sample <- tumorID
    Arribainputfile$Caller <- "Arriba"
    STARFusioninputfile$Caller <- "StarFusion"
    
    # standardized fusion calls
    standardizedSTARFusion<-fusion_standardization(fusion_calls = STARFusioninputfile,caller = "STARFUSION")
    standardizedArriba<-fusion_standardization(fusion_calls = Arribainputfile,caller = "ARRIBA")
    
    #merge standardized fusion calls
    standardFusioncalls<-rbind(standardizedSTARFusion,standardizedArriba) %>% as.data.frame()
  }
  
  # if StarFusion file is empty only run standardization for Arriba calls
  if(is_empty(STARFusioninputfile$FusionName) & !is_empty(Arribainputfile$"gene1--gene2")){
    warning(paste("No fusion calls in StarFusion : ",opt$fusionfileStarFusion))
    
    colnames(Arribainputfile)[27]<-"annots"
    # To have a general column with unique IDs associated with each sample
    Arribainputfile$Sample <- tumorID
    Arribainputfile$Caller <- "Arriba"
    
    # standardized fusion calls
    standardizedArriba<-fusion_standardization(fusion_calls = Arribainputfile,caller = "ARRIBA")
    
    # standardized fusion calls
    standardFusioncalls<-standardizedArriba %>% as.data.frame()
  }
  
  # if Arriba file is empty only run standardization for StarFusion calls
  if(!is_empty(STARFusioninputfile$FusionName) & is_empty(Arribainputfile$"gene1--gene2")){
    warning(paste("No fusion calls in StarFusion : ",opt$fusionfileStarFusion))
    
    # To have a general column with unique IDs associated with each sample
    STARFusioninputfile$Sample <- tumorID
    STARFusioninputfile$Caller <- "StarFusion"
    
    # standardized fusion calls
    standardizedSTARFusion<-fusion_standardization(fusion_calls = STARFusioninputfile,caller = "STARFUSION")
    
    # standardized fusion calls
    standardFusioncalls<-standardizedSTARFusion %>% as.data.frame()
  }
  
  
  # General fusion QC for read support and red flags
  fusionQCFiltered<-fusion_filtering_QC(standardFusioncalls=standardFusioncalls,readingFrameFilter=readingFrameFilter,artifactFilter=artifactFilter,junctionReadCountFilter=junctionReadCountFilter,spanningFragCountFilter=spanningFragCountFilter,readthroughFilter=readthroughFilter)
  
  if (!is.null(expressionFile)){
  expressionMatrix<-read_tsv(expressionFile)
  
  # split gene id and symbol
  expressionMatrix <- expressionMatrix %>% 
    dplyr::mutate(gene_id = str_replace(gene_id, "_PAR_Y_", "_")) 
  
  
  expressionMatrix <- cbind(expressionMatrix, colsplit(expressionMatrix$gene_id, pattern = '_', names = c("EnsembleID","GeneSymbol")))
  
  
  # collapse to matrix of HUGO symbols x Sample identifiers
  # take max expression per row and use the max value for duplicated gene symbols
  expressionMatrix.collapsed <-expressionMatrix  %>% 
    arrange(desc(FPKM)) %>% # arrange decreasing by FPKM
    distinct(GeneSymbol, .keep_all = TRUE) %>% # keep the ones with greatest FPKM value. If ties occur, keep the first occurencce
    unique() %>%
    remove_rownames()  %>%
    dplyr::select(EnsembleID,GeneSymbol,FPKM,gene_id)
  
  # rename columns
  colnames(expressionMatrix.collapsed)[3]<-tumorID
  
  expressionFiltered<-annoFuse::expressionFilterFusion(standardFusioncalls = fusionQCFiltered,expressionFilter = expressionFilter,expressionMatrix = expressionMatrix.collapsed)
  
  # annotated QC and expression filtered fusion calls
  filteredFusionAnnotated<-annotate_fusion_calls(standardFusioncalls=expressionFiltered,geneListReferenceDataTab=geneListReferenceDataTab,fusionReferenceDataTab=fusionReferenceDataTab)
  } else{
    # annotated QC filtered fusion calls
    filteredFusionAnnotated<-annotate_fusion_calls(standardFusioncalls=fusionQCFiltered,geneListReferenceDataTab=geneListReferenceDataTab,fusionReferenceDataTab=fusionReferenceDataTab)
    
  }
  
  
  # return filtered and annotated dataframe
  return(filteredFusionAnnotated)
  
  
}
