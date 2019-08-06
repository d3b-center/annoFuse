#' Filter for putative driver fusions
#'
#' Users can use a driver gene and fusion list or biologically relevant gene
#' sets like tumor suppersor genes, oncogenic genes, tcga fusions
#' @param fusion_calls A dataframe with fusion calls and broad histologies
#' @param lit_fusions table of driver fusiions per histologies
#' @param tsgs dataframe for tumor suppressor gene list
#' @param onco dataframe for oncogenic gene list
#' @param tcga dataframe for tcga fusion
#' @param remove_circularRNA boolean
#' @return filtered putative fusion calls and annotation
#' @examples
#' fusion_calls<-system.file("extdata","filt_read.RData",package="annofusion")
#' load(fusion_calls)
#' ###optional##
#' lit_file<-system.file("extdata", "driver-fusions-v4.txt", package = "annofusion")
#' lit_fusions<-read.delim(lit_file,stringsAsFactors=FALSE)
#' tsgs_file<-system.file("extdata", "tsgs.txt", package = "annofusion")
#' tsgs<-read.delim(tsgs_file,stringsAsFactors=FALSE)
#' onco_file<-system.file("extdata", "oncogenes.txt", package = "annofusion")
#' onco<-read.delim(onco_file,stringsAsFactors=FALSE)
#' tcga_file<-system.file("extdata", "tcga-fusions.txt", package = "annofusion")
#' tcga<-read.delim(tcga_file,stringsAsFactors=FALSE)
#' filt_prioritization(fusion_calls=fusion_calls, tsgs=tsgs,onco=onco,tcga=tcga,remove_circularRNA=TRUE)
#'
#'
#'


filt_prioritization<-function(fusion_calls=fusion_calls,lit_fusions=lit_fusions, tsgs=tsgs,onco=onco,tcga=tcga,remove_circularRNA=TRUE){
  to.add<-""
  if(!missing(lit_fusions)){
  for(i in 1:nrow(lit_fusions)){
    genes.to.search <- lit_fusions[i,2]
    fusions.to.search<-gsub("_","--",lit_fusions[i,3])
    genes.to.search<-unlist(strsplit(genes.to.search, ','))
    print(fusions.to.search)
    genes.to.search <- c(paste0('^',genes.to.search,'-'), paste0('-',genes.to.search,'$'))
    print(genes.to.search)
    if(fusions.to.search == ""){
      print("no fusions to check")
    } else {
      fusions.to.search <- paste0('^',fusions.to.search,'$')
      genes.to.search <- c(genes.to.search, fusions.to.search)
    }
    genes.to.search <- paste0(genes.to.search, collapse = '|')
    hist.to.search <- lit_fusions[1,1]
    getfusions <- fusion_calls[grep(genes.to.search, fusion_calls$FusionName),]
    getfusions <- fusion_calls[which(fusion_calls$broad.histology %in% hist.to.search),]
    getfusions <- unique(getfusions)
    if(nrow(getfusions) == 0){
      print(hist.to.search)
      print(genes.to.search)
    }
    if(i == 1){
      to.add <- getfusions
    } else {
      to.add <- rbind(to.add, getfusions)
    }
  }
  to.add <- unique(to.add)
  colnames(to.add) <- colnames(fusion_calls)
  print("Literature gene fusions")
  to.add$note<-"known oncogenic fusion"
  }

  #tsgs,onco,tcga fusion genes


  genes.to.search <- c(paste0('^',tsgs$GeneSymbol,'-'), paste0('-',tsgs$GeneSymbol,'$'),paste0('^',onco$OncogeneName,'-'), paste0('-',onco$OncogeneName,'$'))
  fusion.to.search<-unlist(strsplit(tcga$TCGA_fusions,","))
  tsgs_onco_fusion<-fusion_calls[unlist(lapply(genes.to.search,function(x) grep(x, fusion_calls$FusionName))),]
  tsgs_onco_fusion<-rbind(tsgs_onco_fusion,fusion_calls[unlist(lapply(fusion.to.search,function(x) grep(x, fusion_calls$FusionName))),])
  tsgs_onco_fusion$note<-"Found in onco/tsgs/tcga list"
  to.add<-rbind(to.add,tsgs_onco_fusion)
  # remove GeneA == GeneB fusion calls
  to.add <- cbind(to.add, colsplit(to.add$FusionName, pattern = '--', names = c("GeneA","GeneB")))
  if(nrow(to.add[-which(to.add$GeneA == to.add$GeneB),])>1 & remove_circularRNA == TRUE){
    to.add<- to.add[-which(to.add$GeneA == to.add$GeneB),]
    }
  to.add$Caller_type<-paste(to.add$Caller,to.add$Fusion_Type,to.add$Confidence,sep="_")

  # Gene fusion should be in-frame
  # Called by at least 2 callers
  fusion_calls.summary <- fusion_calls %>%
    filter(Fusion_Type != "other") %>%
    group_by(FusionName, Sample, broad.histology ,JunctionReadCount,SpanningFragCount) %>%
    unique() %>%
    dplyr::mutate(Caller = toString(Caller), caller.count = n()) %>%
    filter(caller.count >= 2) %>%
    select(-caller.count, -Caller, -Fusion_Type) %>%
    unique() %>%
    as.data.frame()

  print("caller count")
  fusion_calls.summary$note<-"Called by both callers"
  nrow(fusion_calls.summary)

  # or found in at least 2 samples of the same histology
  sample.count <- fusion_calls %>%
    dplyr::filter(Fusion_Type != "other") %>%
    group_by(FusionName, broad.histology,JunctionReadCount,SpanningFragCount) %>%
    unique() %>%
    dplyr::mutate(sample.count = n()) %>%
    filter(sample.count > 1) %>%
    select(-Caller, -sample.count, -Fusion_Type) %>%
    unique() %>%
    as.data.frame()
  length(unique(sample.count$FusionName))
  sample.count$note<-"found in atleast 2 samples in same histology"

  # or GeneB or GeneA gene recurrently fused within a histology (>= 5 genes)
  rec <- cbind(fusion_calls, colsplit(fusion_calls$FusionName, pattern = '--', names = c("GeneA","GeneB")))
  rec2 <- rec %>% group_by(broad.histology) %>%
    select(broad.histology,GeneA,GeneB) %>%
    unique() %>% group_by(broad.histology, GeneA) %>%
    dplyr::summarise(GeneA.ct = n()) %>%
    filter(GeneA.ct >= 5) %>% as.data.frame()
  rec3 <- rec %>% group_by(broad.histology) %>%
    select(broad.histology,GeneA,GeneB) %>%
    unique() %>% group_by(broad.histology, GeneB) %>%
    dplyr::summarise(GeneB.ct = n()) %>%
    filter(GeneB.ct >= 5) %>% as.data.frame()
  rec2 <- merge(rec2, rec, by = c('GeneA','broad.histology'))
  rec3 <- merge(rec3, rec, by = c('GeneB','broad.histology'))
  rec2 <- unique(rec2[,c("Sample","FusionName","broad.histology","JunctionReadCount","SpanningFragCount","Confidence")])
  rec3 <- unique(rec3[,c("Sample","FusionName","broad.histology","JunctionReadCount","SpanningFragCount","Confidence")])
  res <- unique(rbind(rec2, rec3))

  res$note<-"recurrently fused in a histology"

  #check
  colnames(fusion_calls.summary)
  colnames(sample.count)
  colnames(res)


  # merge these
  total <- unique(rbind(fusion_calls.summary, sample.count, res))



  # remove GeneA == GeneB fusion calls
  total <- cbind(total, colsplit(total$FusionName, pattern = '--', names = c("GeneA","GeneB")))
  if(nrow(total[-which(total$GeneA == total$GeneB),])>1 & remove_circularRNA == TRUE){
    total<- total[-which(total$GeneA == total$GeneB),]
  }
  total<-total[,c("FusionName","Sample","broad.histology","note","JunctionReadCount","SpanningFragCount")]


  print("total number of calls")
  head(total)

  # remove fusions that are in > 1 histology
  hist.count <- total %>%
    select(FusionName, broad.histology,JunctionReadCount,SpanningFragCount) %>%
    unique() %>%
    group_by(FusionName) %>%
    dplyr::summarise(hist.count = n()) %>%
    filter(hist.count == 1)
  total <- total[which(total$FusionName %in% hist.count$FusionName),]
  length(unique(total$FusionName))



  total <- merge(total, fusion_calls, by = c('Sample','FusionName','broad.histology','JunctionReadCount','SpanningFragCount'))
  total$Caller_type<-paste(total$Caller,total$Fusion_Type,total$Confidence,sep="_")

  #filter for multi_fused
  rec <- cbind(colsplit(total$FusionName, pattern = '--', names = c("Gene1","Gene2")), total)
  rec <- cbind(colsplit(rec$Gene2, pattern = '/', names = c("Gene2A","Gene2B")), rec)
  rec <- cbind(colsplit(rec$Gene1, pattern = '/', names = c("Gene1A","Gene1B")), rec)
  rec[grep("\\(", rec$Gene1B),"Gene1B"]<-gsub("\\(.*","",rec[grep("\\(", rec$Gene1B),"Gene1B"])
  rec[grep("\\(", rec$Gene2B),"Gene2B"]<-gsub("\\(.*","",rec[grep("\\(", rec$Gene2B),"Gene2B"])
  rec[grep("\\(", rec$Gene1A),"Gene1A"]<-gsub("\\(.*","",rec[grep("\\(", rec$Gene1A),"Gene1A"])
  rec[grep("\\(", rec$Gene2A),"Gene2A"]<-gsub("\\(.*","",rec[grep("\\(", rec$Gene2A),"Gene2A"])


  rec4 <- rec %>% group_by(Sample, Gene1A,Gene1B) %>% dplyr::summarise(Gene1.ct = n())  %>% as.data.frame() %>% cbind("multi_fused"="")
  rec5 <- rec %>% group_by(Sample, Gene2A,Gene2B) %>% dplyr::summarise(Gene2.ct = n())  %>% as.data.frame() %>% cbind("multi_fused"="")

  rec4$multi_fused<-ifelse(length(rec4$Sample[rec4$Gene1.ct>5])>0,"Gene1-multi_fused","")
  rec5$multi_fused<-ifelse(length(rec5$Sample[rec5$Gene1.ct>5])>0,"Gene1-multi_fused","")
  rec4<-merge(rec, rec4, by=c("Sample","Gene1A","Gene1B"))
  rec5<-merge(rec, rec5, by=c("Sample","Gene2A","Gene2B"))



  total_multifused<-rec4[,c("FusionName","Sample","Caller_type","broad.histology","note","multi_fused","JunctionReadCount","SpanningFragCount")]
  total_multifused<-rbind(total_multifused,rec5[,c("FusionName","Sample","Caller_type","broad.histology","note","multi_fused","JunctionReadCount","SpanningFragCount")])

  total_multifused<-unique(total_multifused)

  #merged with multifused genes
  total<-merge(rec,total_multifused, by=c("Sample","FusionName","Caller_type","note","broad.histology","JunctionReadCount","SpanningFragCount"),all=T)
  total<-unique(total[,c("FusionName","Sample","Caller_type","broad.histology","JunctionReadCount","SpanningFragCount")])

  #filtered fusion through literature and multi_fused status
  to.add_not_multi<-total[-which(total$multi_fused == "Gene1-multi_fused" | total$multi_fused == "Gene2-multi_fused"),]
  colnames(to.add)
  colnames(to.add_not_multi)


  final<-rbind(to.add[,c("FusionName","Sample","Caller_type","broad.histology","JunctionReadCount","SpanningFragCount")],to.add_not_multi [,c("FusionName","Sample","Caller_type","broad.histology","JunctionReadCount","SpanningFragCount")])
  final<-unique(final)


  filt_fusion_calls<-list("driver"=final,"prioritized"=total)

  return(filt_fusion_calls)


}
