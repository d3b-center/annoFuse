#' Function to indetify domains retained and lost per breakpoint in gene
#'
#' @param genesymbol gene symbol associated with the fusion breakpoint
#' @param ensembldb annotation databse from ensembldb
#' @param fusionbk chr:position format breakpoint in fusion
#' @param geneposition RIGHT or LEFT gene
#' @return Standardized fusion calls with aggregated Caller and read support



domain_overlap_fusion<-function(genesymbol=genesymbol,ensembldb=ensembldb,fusionbk=fusionbk,geneposition=geneposition){
# get chromosome and breakpoint
fusion<-list("bk"=gsub("^.*:","",fusionbk),"chr"=gsub(":.*$","",fusionbk))

# transcripts for gene of interest
txs <- getGeneRegionTrackForGviz(edb, filter = ~ genename %in% genesymbol )

# protein domain information
pdoms <- proteins(edb, filter = ~ tx_id %in% txs$transcript &
                    protein_domain_source == "pfam",
                  columns = c("protein_domain_id", "prot_dom_start",
                              "prot_dom_end","gene_name","gene_seq_start","gene_seq_end","seq_strand","seq_length"))

if (unique(pdoms$seq_strand)=="-1"){
  temp<-pdoms$gene_seq_start
  pdoms$gene_seq_start<-pdoms$gene_seq_end
  pdoms$gene_seq_end<-temp
  }

# Iramge of domain start stop
pdoms_rng <- IRanges(start = pdoms$prot_dom_start, end = pdoms$prot_dom_end,
                     names = pdoms$protein_id)

# protein domain from prtein to genome coordinates
pdoms_gnm <- proteinToGenome(pdoms_rng, edb)
pdoms_gnm_grng <- unlist(GRangesList(pdoms_gnm))

# fusion grange
if( geneposition=="LEFT" & pdoms$seq_strand== 1){
  fusion_grng<-GRanges(data.frame("start"=unique(pdoms$gene_seq_start),"end"=fusion$bk,"seqnames"=fusion$chr))
  lost_grng<-GRanges(data.frame("start"=fusion$bk,"end"=unique(pdoms$gene_seq_end),"seqnames"=fusion$chr))
} else if (geneposition=="RIGHT" & pdoms$seq_strand== 1) {
  fusion_grng<-GRanges(data.frame("start"=fusion$bk,"end"=unique(pdoms$gene_seq_end),"seqnames"=fusion$chr))
  lost_grng<-GRanges(data.frame("start"=unique(pdoms$gene_seq_start),"end"=fusion$bk,"seqnames"=fusion$chr))
} else if( geneposition=="LEFT" & pdoms$seq_strand== -1){
  lost_grng<-GRanges(data.frame("end"=fusion$bk,"start"=unique(pdoms$gene_seq_end),"seqnames"=fusion$chr))
  fusion_grng<-GRanges(data.frame("end"=unique(pdoms$gene_seq_start),"start"=fusion$bk,"seqnames"=fusion$chr))
} else if (geneposition=="RIGHT" & pdoms$seq_strand== -1) {
  lost_grng<-GRanges(data.frame("end"=unique(pdoms$gene_seq_start),"start"=fusion$bk,"seqnames"=fusion$chr))
  fusion_grng<-GRanges(data.frame("end"=fusion$bk,"start"=unique(pdoms$gene_seq_end),"seqnames"=fusion$chr))
}

#overlap
domainoverlap<-subsetByOverlaps(pdoms_gnm_grng,fusion_grng,type = "within")
#lost
domainlost<-subsetByOverlaps(pdoms_gnm_grng,lost_grng,type = "within")

domains<-list("retained"=unique(domainoverlap),"lost"=unique(domainlost))
return(domains)
}

