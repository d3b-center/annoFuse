

#' Fusion annotation with tsgs oncogenes kinase tcga fusion
#'
#' @param fusion_calls list of putative drivers and total filtered fusion calls
#' @param clin clinical info
#' @param tsgs tsgs gene list
#' @param onco oncogenic gene list
#' @param tcga tcga fusion list
#' @param kinase kinase gene list
#' @param chimerDB chimerDB
#' @param tf transcription factor gene list
#'
#' @return filtered annotated fusion calls
#'
#' @examples
#' driver_fusions<-system.file("extdata","fusion_exp_filt.Rdata",package="annofusion")
#' prioritized_fusions<-system.file("extdata","fusion_exp_filt_prioritized.Rdata",package="annofusion")
#' load(driver_fusions)
#' load(prioritized_fusions)
#' filt_fusion_exp_genes<-list("driver"=driver,"prioritized"=prioritized)
#' tsgs_file<-system.file("extdata", "tsgs.txt", package = "annofusion")
#' tsgs<-read.delim(tsgs_file,stringsAsFactors=FALSE)
#' onco_file<-system.file("extdata", "oncogenes.txt", package = "annofusion")
#' onco<-read.delim(onco_file,stringsAsFactors=FALSE)
#' tcga_file<-system.file("extdata", "tcga-fusions.txt", package = "annofusion")
#' tcga<-read.delim(tcga_file,stringsAsFactors=FALSE)
#' kinase_file<-system.file("extdata", "Kincat_Hsap.08.02.txt", package = "annofusion")
#' kinase<-read.delim(kinase_file,stringsAsFactors=FALSE)
#' chimerDB_file<-system.file("extdata","ChimerDB3.0_ChimerKB.txt", package = "annofusion")
#' chimerDB<-read.delim(chimerDB_file)
#' tf_file<-system.file("extdata","TRANSFAC_TF.txt", package = "annofusion")
#' tf<-read.delim(tf_file,header=FALSE,stringsAsFactor=FALSE)
#' clin_file<-system.file("extdata","2019-06-07-pbta-histologies.tsv",package="annofusion")
#' clin<-read.delim(clin_file,sep="\t",stringsAsFactors=FALSE)
#' clin$Sample<-paste(clin$sample_id,clin$aliquot_id,sep="_")
#' fusion_anno(fusion_calls=filt_fusion_exp_genes, tsgs=tsgs,onco=onco,tcga=tcga,kinase=kinase,chimerDB=chimerDB,tf=tf,clin=clin)
#'


fusion_anno<-function(fusion_calls=filt_fusion_exp_genes,clin=clin,tsgs=tsgs,onco=onco,tcga=tcga,kinase=kinase,chimerDB=chimerDB,tf=tf){
####### separate low expressing fusions and fusions where gene expression not reported
total_filtered<-filt_fusion_exp_genes$prioritized
total_putativedrivers<-filt_fusion_exp_genes$driver
total_filtered<-total_filtered[,which (colnames(total_filtered) %in% colnames(total_putativedrivers))]
total<-rbind(total_filtered,total_putativedrivers)

colnames(total)

#total <- cbind(total, colsplit(total$FusionName, pattern = '--', names = c("Gene1","Gene2")))
for.table<-total

############### fix gene names and add entrez ids ####################
hs <- org.Hs.eg.db
genes<-c(total_filtered$Gene1A,total_filtered$Gene1B,total_filtered$Gene2A,total_filtered$Gene2B)
my.symbols<-genes
#entrez<-select(hs, keys=my.symbols, columns= c("ENTREZID", "SYMBOL"),keytype = "SYMBOL")
#nomatch<-entrez[is.na(entrez$ENTREZID),]
#match<-read.delim("../../references/2019-02-14-Hugo-Symbols-approved (1).txt")
#previous_symbol<-match[which(match$Previous.symbols %in% nomatch[-which(nomatch$SYMBOL==""),"SYMBOL"]),]
#synonym<-match[which(match$Synonyms %in% nomatch[-which(nomatch$SYMBOL==""),"SYMBOL"]),]








#total$Gene1A<-sub(previous_symbol$Previous.symbols, previous_symbol$Approved.symbol, total$Gene1A)
#total$Gene2A<-sub(previous_symbol$Previous.symbols, previous_symbol$Approved.symbol, total$Gene2A)
#total$Gene2B<-sub(previous_symbol$Previous.symbols, previous_symbol$Approved.symbol, total$Gene2B)
#total$Gene1B<-sub(previous_symbol$Previous.symbols, previous_symbol$Approved.symbol, total$Gene1B)

#total$Gene1A<-sub(synonym$Synonyms, synonym$Approved.symbol, total$Gene1A)
#total$Gene2A<-sub(synonym$Synonyms, synonym$Approved.symbol, total$Gene2A)
#total$Gene2B<-sub(synonym$Synonyms, synonym$Approved.symbol, total$Gene2B)
#total$Gene1B<-sub(synonym$Synonyms, synonym$Approved.symbol, total$Gene1B)


#IGH-@ chnaged to IGH to get Entrez ID
total$Gene1A<-sub("IGH-@", "IGH", total$Gene1A)
total$Gene1B<-sub("IGH-@", "IGH", total$Gene1B)
total$Gene2B<-sub("IGH-@", "IGH", total$Gene2B)
total$Gene2A<-sub("IGH-@", "IGH", total$Gene2A)

for.table$Gene1A_entrez<-""
for.table$Gene1B_entrez<-""
for.table$Gene2A_entrez<-""
for.table$Gene2B_entrez<-""

#
# my.symbols <- for.table$Gene1A
# Gene1A_entrez<-dplyr::select(hs,
#                              keys = my.symbols,
#                              columns = c("ENTREZID", "SYMBOL"),
#                              keytype = "SYMBOL")
#
# Gene1A_entrez[is.na(Gene1A_entrez$ENTREZID),"ENTREZID"]=""
#
#
# my.symbols <- for.table$Gene2A
# Gene2A_entrez<-dplyr::select(hs,
#                              keys = my.symbols,
#                              columns = c("ENTREZID", "SYMBOL"),
#                              keytype = "SYMBOL")
#
#
# Gene2A_entrez[is.na(Gene2A_entrez$ENTREZID),"ENTREZID"]=""
#
#
# if(length(for.table$Gene1B[!is.na(for.table$Gene1B)]) >0 & length(for.table$Gene1B[!is.na(for.table$Gene2B)]) >0){
#
#   my.symbols <- for.table$Gene1B
#   Gene1B_entrez<-dplyr::select(hs,
#                                keys = my.symbols,
#                                columns = c("ENTREZID", "SYMBOL"),
#                                keytype = "SYMBOL")
#
#   Gene1B_entrez[is.na(Gene1B_entrez$ENTREZID),"ENTREZID"]=""
#
#   my.symbols <- for.table$Gene2B
#   Gene2B_entrez<-dplyr::select(hs,
#                                keys = my.symbols,
#                                columns = c("ENTREZID", "SYMBOL"),
#                                keytype = "SYMBOL")
#
#   Gene2B_entrez[is.na(Gene2B_entrez$ENTREZID),"ENTREZID"]=""
# }
#
# ####################################################################







# annotate using kinases
kinase <- kinase[-which(kinase$Entrez_Symbol == ""),]
for.table$Gene1A_Kinase <- ifelse(for.table$Gene1A %in% kinase$Entrez_Symbol,'Yes', 'No')
for.table$Gene1B_Kinase <- ifelse(for.table$Gene1B %in% kinase$Entrez_Symbol,'Yes', 'No')
for.table$Gene2A_Kinase <- ifelse(for.table$Gene2A %in% kinase$Entrez_Symbol,'Yes', 'No')
for.table$Gene2B_Kinase <- ifelse(for.table$Gene2B %in% kinase$Entrez_Symbol,'Yes', 'No')
for.table$Kinase<-paste0(for.table$Gene1A_Kinase,";",for.table$Gene1B_Kinase,";",for.table$Gene2A_Kinase,";",for.table$Gene2B_Kinase)

head(for.table)





print("Kinase")
table(for.table$Kinase)

# annotate using transcription factors
for.table$Gene1A_TF <- ifelse(for.table$Gene1A %in% tf$V1,'Yes', 'No')
for.table$Gene1B_TF <- ifelse(for.table$Gene1B %in% tf$V1,'Yes', 'No')
for.table$Gene2A_TF <- ifelse(for.table$Gene2A %in% tf$V1 ,'Yes', 'No')
for.table$Gene2B_TF <- ifelse(for.table$Gene2B %in% tf$V1,'Yes', 'No')
for.table$TF<-paste0("",for.table$Gene1A_TF,";",for.table$Gene1B_TF,";",for.table$Gene2A_TF,";",for.table$Gene2B_TF)

print("TF")
table(for.table$TF)

# annotate oncogenic
for.table$Gene1A_is_onco <- ifelse(for.table$Gene1A %in% onco$OncogeneName,'Yes', 'No')
for.table$Gene1B_is_onco <- ifelse(for.table$Gene1B %in% onco$OncogeneName,'Yes', 'No')
for.table$Gene2A_is_onco <- ifelse(for.table$Gene2A %in% onco$OncogeneName,'Yes', 'No')
for.table$Gene2B_is_onco <- ifelse(for.table$Gene2B %in% onco$OncogeneName,'Yes', 'No')
for.table$is_onco<-paste0("",for.table$Gene1A_is_onco,";",for.table$Gene1B_is_onco,";",for.table$Gene2A_is_onco,";",for.table$Gene2B_is_onco)
print("onco")
table(for.table$is_onco)

# annotate tsgs
#for.table$is_TSG <- ifelse(for.table$Gene1A %in% tsg$GeneSymbol |for.table$Gene1B %in% tsg$GeneSymbol|for.table$Gene2A %in% tsg$GeneSymbol| for.table$Gene2B %in% tsg$GeneSymbol, 'Yes', 'No')
for.table$Gene1A_is_TSG <- ifelse(for.table$Gene1A %in% tsgs$GeneSymbol,'Yes', 'No')
for.table$Gene1B_is_TSG <- ifelse(for.table$Gene1B %in% tsgs$GeneSymbol,'Yes', 'No')
for.table$Gene2A_is_TSG <- ifelse(for.table$Gene2A %in% tsgs$GeneSymbol,'Yes', 'No')
for.table$Gene2B_is_TSG <- ifelse(for.table$Gene2B %in% tsgs$GeneSymbol,'Yes', 'No')

print("TSG")
for.table$is_TSG<-paste0("",for.table$Gene1A_is_TSG,";",for.table$Gene1B_is_TSG,";",for.table$Gene2A_is_TSG,";",for.table$Gene2B_is_TSG)

table(for.table$is_TSG)

#tcga fusion
tcga<-separate_rows(tcga,TCGA_fusions,sep=",")
tcga$TCGA_fusions<-gsub("_","--",tcga$TCGA_fusions)
for.table$in_TCGA <- ifelse(for.table$FusionName %in% tcga$TCGA_fusions,'Yes', 'No')


#for.table$Gene1A_entrez<-Gene1A_entrez$ENTREZID
#for.table$Gene1B_entrez<-Gene1B_entrez$ENTREZID
#for.table$Gene2A_entrez<-Gene2A_entrez$ENTREZID
#for.table$Gene2B_entrez<-Gene2B_entrez$ENTREZID

if(nrow(unique(for.table[-which(for.table$is_TSG == "No;No;No;No" & for.table$is_onco== "No;No;No;No" & for.table$Kinase =="No;No;No;No" & for.table$TF == "No;No;No;No"),])
) >0){
for.table <- unique(for.table[-which(for.table$is_TSG == "No;No;No;No" & for.table$is_onco== "No;No;No;No" & for.table$Kinase =="No;No;No;No" & for.table$TF == "No;No;No;No"),])
}

filt_anno<-for.table




#clin<-clin[,c("Sample","broad.histology")]


#histology labels
hist.ct <- unique(clin[,c('broad.histology','Sample')])
hist.dt.ct <- unique(clin[,c('broad.histology','Sample')])
hist.ct<-as.data.frame(hist.ct)
hist.dt.ct <- unique(hist.ct %>% dplyr::select(broad.histology) %>% group_by(broad.histology) %>% dplyr::mutate(freq = n()) %>% as.data.frame())
hist.dt.ct$freq <- paste0(hist.dt.ct$broad.histology,' (n=', hist.dt.ct$freq, ')')
colnames(hist.dt.ct)[2] <- 'broad.histology.Label'



#get DriverFusion.txt
# excel tab 1
extab1<-total_filtered

extab1_note<-aggregate(list(extab1$Caller_type), list(extab1$FusionName,extab1$Sample,extab1$broad.histology), paste, collapse=",")
colnames(extab1_note)<-c("FusionName","Sample","broad.histology","Caller_type")

extab1$SpanningFragCount<-as.numeric(extab1$SpanningFragCount)
extab1$JunctionReadCount<-as.numeric(extab1$JunctionReadCount)
extab1_read<-aggregate(list(extab1$SpanningFragCount,extab1$JunctionReadCount), list(extab1$Sample,extab1$broad.histology,extab1$FusionName) ,sum)
colnames(extab1_read)<-c('Sample','broad.histology','FusionName','SpanningFragCount','JunctionReadCount')

extab1_merge<-merge(extab1_read,extab1_note,by=c('Sample','broad.histology','FusionName'))
extab1<-merge(extab1[,-which(colnames(extab1) %in% c("Caller_type","SpanningFragCount","JunctionReadCount"))],extab1_merge,by=c('Sample','broad.histology','FusionName'))

extab1 <- merge(extab1, hist.dt.ct, by = 'broad.histology')
#extab1 <- extab1[,c("Sample","FusionName","broad.histology","Caller_type","Fusion_Type")]
extab1 <- extab1[order(extab1$FusionName, extab1$Sample),]


#Driver only
outextab1 <- merge(extab1, unique(for.table[,c("FusionName","is_onco","is_TSG","TF","Kinase","in_TCGA","Gene1A_entrez","Gene2A_entrez","Gene1B_entrez","Gene2B_entrez")]), by = c('FusionName'))
head(outextab1)
outextab1<-separate(outextab1,is_onco,into=c("Gene1A_is_onco","Gene1B_is_onco","Gene2A_is_onco","Gene2B_is_onco"),sep=";")
outextab1<-separate(outextab1,is_TSG,into=c("Gene1A_is_TSG","Gene1B_is_TSG","Gene2A_is_TSG","Gene2B_is_TSG"),sep=";")
outextab1<-separate(outextab1,Kinase,into=c("Gene1A_Kinase","Gene1B_Kinase","Gene2A_Kinase","Gene2B_Kinase"),sep=";")
outextab1<-separate(outextab1,TF,into=c("Gene1A_TF","Gene1B_TF","Gene2A_TF","Gene2B_TF"),sep=";")

anno_driver<-unique(outextab1)

#add chimerDB
chimerDB<-chimerDB[,-which(colnames(chimerDB) %in% c("H_gene","H_chr","H_position","H_strand","T_gene","T_chr","T_position","T_strand"))]
colnames(chimerDB)<-paste("chimerDB",colnames(chimerDB),sep=":")
chimerDB$`chimerDB:Fusion_pair`<-sub("_","--",chimerDB$`chimerDB:Fusion_pair`)


outextab1<-merge(outextab1,chimerDB,by.x="FusionName",by.y="chimerDB:Fusion_pair",all.x=T)


# excel tab 2
print("fortable")

#extab2 <- merge(filt_anno, hist.dt.ct, by =c('broad.histology'))
#extab2 <- extab2[,c("Sample","FusionName","broad.histology","in_TCGA","is_onco","is_TSG","TF","Kinase")]

extab2<-unique(filt_anno)

print("cols in extab2")
extab2<-extab2[,-which(colnames(extab2) %in% c('SpanningFragCount','JunctionReadCount','Caller_type'))]
extab2<-unique(extab2)



extab2 <- extab2 %>%
  group_by(FusionName, in_TCGA, is_onco,is_TSG, TF, Kinase,Gene1A_entrez,Gene2A_entrez,Gene1B_entrez,Gene2B_entrez) %>%
  dplyr::summarise(Samples = toString(Sample), Samples.Found.In = n()) %>%
  unique() %>% as.data.frame()
#extab2 <- unique(extab2[,c("FusionName","is_onco","is_TSG","TF","Kinase","in_TCGA")])

total_note<-aggregate(list(total$Caller_type), list(total$FusionName,total$Sample,total$broad.histology), paste, collapse=",")
colnames(total_note)<-c("FusionName","Sample","broad.histology","Caller_type")

total$JunctionReadCount<-as.numeric(total$JunctionReadCount)
total$SpanningFragCount<-as.numeric(total$SpanningFragCount)
total_read<-aggregate(list(total$SpanningFragCount,total$JunctionReadCount), list(total$Sample,total$broad.histology,total$FusionName) ,sum)
colnames(total_read)<-c('Sample','broad.histology','FusionName','SpanningFragCount','JunctionReadCount')


total_merge<-merge(total_read,total_note,by=c('Sample','broad.histology','FusionName'))

print("total_merge")
colnames(total_merge)

print("total")
head(total)
which(colnames(total) %in% c("Caller_type","SpanningFragCount","JunctionReadCount"))

total<-merge(total[,-which(colnames(total) %in% c("Caller_type","SpanningFragCount","JunctionReadCount"))],total_merge,by=c('Sample','FusionName','broad.histology'))






extab2 <- merge(total,extab2, by = c("FusionName"))
#extab2[is.na(extab2)] <- "No"
extab2<-unique(extab2)

head(extab2)

#All Filtered
extab2<-separate(extab2,is_onco,into=c("Gene1A_is_onco","Gene1B_is_onco","Gene2A_is_onco","Gene2B_is_onco"),sep=";")
extab2<-separate(extab2,is_TSG,into=c("Gene1A_is_TSG","Gene1B_is_TSG","Gene2A_is_TSG","Gene2B_is_TSG"),sep=";")
extab2<-separate(extab2,Kinase,into=c("Gene1A_Kinase","Gene1B_Kinase","Gene2A_Kinase","Gene2B_Kinase"),sep=";")
extab2<-separate(extab2,TF,into=c("Gene1A_TF","Gene1B_TF","Gene2A_TF","Gene2B_TF"),sep=";")

anno_prioritised<-merge(extab2,chimerDB,by.x="FusionName",by.y="chimerDB:Fusion_pair",all.x=T)

#########
# plot for fusions only
#########
total$Gene1 <- NULL
total$Gene2 <- NULL


# histology detailed labels
hist.ct <- unique(clin[,c('broad.histology','Sample')])
hist.ct <- unique(hist.ct %>% dplyr::select(broad.histology) %>% group_by(broad.histology) %>% dplyr::mutate(freq = n()) %>% as.data.frame())
hist.ct$freq <- paste0(hist.ct$broad.histology,' (n=', hist.ct$freq, ')')
colnames(hist.ct)[2] <- 'broad.histology.Label'
#write.table(hist.ct,"hist_lable_freq.tsv",sep="\t")


#final <- merge(total, clin[,c('Sample','broad.histology')], by = 'Sample')
final<-total
head(final)

head(hist.ct)

final <- merge(final, hist.ct, by= 'broad.histology')
final <- final %>% group_by(Sample,broad.histology) %>% dplyr::summarise(value = n())
final <- final %>% group_by(broad.histology) %>% mutate(median = median(value)) %>% as.data.frame()
to.include <- setdiff(hist.ct$broad.histology, final$broad.histology)
 final <- rbind(final, data.frame(Sample = c(rep(NA,length(to.include))),
                                  broad.histology = to.include,
                                 value = c(rep(0,length(to.include))),
                                 median = c(rep(0,length(to.include)))))
final$broad.histology <- reorder(final$broad.histology, final$median)
#write.table(final, file = paste0(basePD,'/data/processed/FusionPlot_rawdata.txt'), quote = F, sep = "\t", row.names = F)


p <- ggplot(final, aes(x = broad.histology, y = log2(value), color = broad.histology, alpha = 0.5)) +
  geom_boxplot(outlier.shape = 21, fill = 'white') +
  geom_jitter(position=position_jitter(width=.1, height=0), shape = 21) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  theme_bw() +
  guides(alpha = FALSE, fill = FALSE) +
  xlab("Histology") + ylab('Number of Fusions') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  guides(color = F) +
  scale_y_continuous(breaks = seq(0, 26, by = 5))
p
#ggsave(filename = paste0(basePD,'/reports/figures/FusionPlot.pdf'), plot = p, device = 'pdf', height = 6, width = 10)

medians <- unique(final[,c("broad.histology","median")])
#write.table(medians, file = paste0(basePD,'/data/processed/FusionPlot_medians.txt'), quote = F, sep = "\t", row.names = F)


list_anno_fusions<-list("anno_driver"=anno_driver,"anno_prioritised"=anno_prioritised,"plot"=p)
return(list_anno_fusions)

}

