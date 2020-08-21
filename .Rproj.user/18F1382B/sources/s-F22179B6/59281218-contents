#' Filter fusion calls for read depth, frame information and artifacts
#' Inframe and frameshift fusions with high supporting reads are considered true calls for further annotation
#'
#' @param fusion_calls A dataframe from star fusion or arriba (more callers to be added )
#' @param read_support An integer for spanning read filter
#' @param caller string options STARfusion/arriba
#' @param remove_artifacts boolean
#' @return Inframe and frameshift filtered fusion calls for user provided read depth
#' @examples
#' fusion_file<-system.file("extdata", "test_star.tsv", package = "annofusion")
#' fusion_calls<-read.delim(fusion_file,stringsAsFactors=FALSE)
#' read_support<-10
#' caller<-"STARfusion"
#' filt_read(fusion_calls=fusion_calls,read_support=read_support, caller=caller)
#'
#'


filt_read <- function(fusion_calls=fusion_calls,read_support=read_support, caller=caller,remove_artifacts=TRUE) {
  read_support<-as.numeric(read_support)
  if(caller=="STARfusion"){
  fusion_calls<-fusion_calls[-which(fusion_calls$SpanningFragCount-fusion_calls$JunctionReadCount >read_support|fusion_calls$JunctionReadCount==0|fusion_calls$LargeAnchorSupport == "NO_LDAS"),]
  fusion_calls$LeftBreakpoint <- gsub('^chr','',fusion_calls$LeftBreakpoint)
  fusion_calls$RightBreakpoint <- gsub('^chr','',fusion_calls$RightBreakpoint)
  colnames(fusion_calls)[c(6,8)] <- c('Gene1_pos','Gene2_pos')
  fusion_calls$Fusion_Type<-""
  fusion_calls$Fusion_Type[which(fusion_calls$PROT_FUSION_TYPE=="INFRAME")] <- 'in-frame'
  fusion_calls$Fusion_Type[grep("FRAMESHIFT",fusion_calls$PROT_FUSION_TYPE)] <- 'frameshift'
  fusion_calls$Fusion_Type[-which(fusion_calls$PROT_FUSION_TYPE %in% c("INFRAME","FRAMESHIFT"))] <- 'other'
  fusion_calls$Caller <- 'STARFusion'
  fusion_calls$Sample <-paste(fusion_calls$sample_id,fusion_calls$aliqout_id,sep="_")
  fusion_calls$FusionName<-fusion_calls$X.FusionName
  fusion_calls$Confidence<-"NA"
  if(remove_artifacts & length(grep('readthrough|neighbors|GTEx_Recurrent',fusion_calls$annots))>1 ){
    fusion_calls<-fusion_calls[-grep('readthrough|neighbors|GTEx_Recurrent',fusion_calls$annots),]
  }

  fusion_calls.total <- unique(fusion_calls[,c('FusionName','Sample','Caller','Fusion_Type','JunctionReadCount','SpanningFragCount','Confidence')])
  return(fusion_calls.total)

  }
  if( caller=="arriba"){
    #remove false positives through supporting reads
    fusion_calls<-fusion_calls[-(which(fusion_calls$discordant_mates-(fusion_calls$split_reads1+fusion_calls$split_reads2) >read_support |fusion_calls$split_reads1+fusion_calls$split_reads2==0)),]
    fusion_calls$LeftBreakpoint <- gsub('^chr','',fusion_calls$breakpoint1)
    fusion_calls$RightBreakpoint <- gsub('^chr','',fusion_calls$breakpoint2)
    fusion_calls$Fusion_Type<-fusion_calls$reading_frame
    fusion_calls$Fusion_Type[grep("out-of-frame",fusion_calls$Fusion_Type)]<-"frameshift"
    fusion_calls$Fusion_Type[-which(fusion_calls$Fusion_Type %in% c("in-frame","frameshift"))] <- 'other'

    fusion_calls$Caller <- 'arriba'
    fusion_calls$Sample <-paste(fusion_calls$sample_id,fusion_calls$aliqout_id,sep="_")
    fusion_calls$FusionName <-paste0(gsub(",","/",fusion_calls$X.gene1),"--",gsub(",","/",fusion_calls$gene2))
    fusion_calls$SpanningFragCount<-fusion_calls$discordant_mates
    fusion_calls$JunctionReadCount<-fusion_calls$split_reads1+fusion_calls$split_reads2
    fusion_calls$Confidence<-fusion_calls$confidence
    if(remove_artifacts & length(grep('readthrough|neighbors|GTEx_Recurrent',fusion_calls$annots))>1){
      fusion_calls<-fusion_calls[-grep("read-through|non-canonical_splicing",fusion_calls$type),]
    }
    fusion_calls.total <- unique(fusion_calls[,c('FusionName','Sample','Caller','Fusion_Type','JunctionReadCount','SpanningFragCount','Confidence')])
    return(fusion_calls.total)
  }

}






