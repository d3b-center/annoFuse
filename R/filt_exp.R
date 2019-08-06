
#' Filter fusion calls with 5' and 3' genes expressed
#
#' Filtering for expresion of 3' and 5' allows us to look for relevant fusions
#'
#' @param fusion_calls data frame of fusion calls with histologies
#' @param expdata tsv file of rsem data
#' @param clin clinical file to annotate
#'
#' @return filtered fusion calls
#'
#' @examples
#' fusion_calls<-system.file("extdata","filt_fusion_calls.RData",package="annofusion")
#' rsem_file<-system.file("extdata","test_rsem.RData",package="annofusion")
#' clin<-system.file("extdata","2019-06-07-pbta-histologies.tsv",package="annofusion")
#' load(fusion_calls)
#' load(rsem_file)
#' clin<-read.delim(clin,sep="\t",stringsAsFactors=FALSE)
#' fusion_calls<-filt_fusion_calls$driver
#' filt_exp(fusion_calls=fusion_calls,expdata=expdata,clin=clin)


filt_exp<-function(fusion_calls=fusion_calls,expdata=expdata,clin=clin){
rna.mat<-expdata
rna.mat<-as.data.frame(rna.mat)
rna.mat$gene_short_name<-unlist(lapply(rownames(rna.mat),function(x) strsplit(as.character(x),"_")[[1]][2]))

head(fusion_calls)

####### separate low expressing fusions and fusions where gene expression not reported
rna.mat$not_expressed <- apply(rna.mat[,2:ncol(rna.mat)], 1, FUN = function(x) all(x < 1))
table(rna.mat$not_expressed)
df <- fusion_calls
df <- cbind(colsplit(df$FusionName, pattern = '--', names = c("Gene1","Gene2")), df)
df <- cbind(colsplit(df$Gene2, pattern = '/', names = c("Gene2A","Gene2B")), df)
df <- cbind(colsplit(df$Gene1, pattern = '/', names = c("Gene1A","Gene1B")), df)

df[grep("\\(", df$Gene1B),"Gene1B"]<-gsub("\\(.*","",df[grep("\\(", df$Gene1B),"Gene1B"])
df[grep("\\(", df$Gene2B),"Gene2B"]<-gsub("\\(.*","",df[grep("\\(", df$Gene2B),"Gene2B"])
df[grep("\\(", df$Gene1A),"Gene1A"]<-gsub("\\(.*","",df[grep("\\(", df$Gene1A),"Gene1A"])
df[grep("\\(", df$Gene2A),"Gene2A"]<-gsub("\\(.*","",df[grep("\\(", df$Gene2A),"Gene2A"])



genes <- unique(c(df$Gene1A, df$Gene2A,df$Gene1B,df$Gene2B))
to.check <- setdiff(genes, rna.mat$gene_short_name) # 90
print("genes not in exp")
to.check

rna.mat <- rna.mat[which(rna.mat$gene_short_name %in% genes),]
rna.mat <- melt(rna.mat)


#rna.mat$sample_id<-gsub("_.*","",rna.mat$variable)
rna.mat$sample_id<-rna.mat$variable

head(rna.mat)


clin<-unique(clin[,c("sample_id","aliquot_id","disease_type","short.histology","broad.histology","composition","Kids.First.Biospecimen.ID","Kids.First.Participant.ID")])
clin$sample<-paste(clin$sample_id,clin$aliquot_id,sep="_")
head(clin)

rna.mat <- merge(rna.mat, clin[,c("sample","broad.histology")], by.x = 'sample_id', by.y = "sample")
rna.mat$broad.histology <- as.character(rna.mat$broad.histology)
nrow(rna.mat)



fusion_calls<-df
nrow(fusion_calls)


# now add filter for  < 1 TPM


#geneA in fused genes and intergenic breakpoints
df$Gene1A_sample <- NA
df$Gene1A_hist_mean <- NA
df$Gene1A_mean_dataset <- NA
df$Gene1A_expr <- NA
df$Gene2A_sample <- NA
df$Gene2A_hist_mean <- NA
df$Gene2A_mean_dataset <- NA
df$Gene2A_expr <- NA
df$Gene1A_not_expressed <- NA
df$Gene2A_not_expressed <- NA




for(i in 1:nrow(df)){
  print(i)
  genea <- df[i,'Gene1A']
  geneb <- df[i,'Gene2A']
  sample <- df[i,'Sample']
  hist <- df[i,'broad.histology']
  hist <- gsub(' [(].*','', hist)
  genea.expr <- unique(rna.mat[which(rna.mat$gene_short_name == genea),'not_expressed'])
  geneb.expr <- unique(rna.mat[which(rna.mat$gene_short_name == geneb),'not_expressed'])
  genea.val <- rna.mat[which(rna.mat$sample_id %in% sample & rna.mat$gene_short_name == genea),'value']
  geneb.val <- rna.mat[which(rna.mat$sample_id %in% sample & rna.mat$gene_short_name == geneb),'value']
  genea.expr <- genea.val< 1
  geneb.expr <- geneb.val< 1
  genea.hist.mean <- mean(rna.mat[which(rna.mat$gene_short_name == genea & rna.mat$broad.histology == hist),'value'])
  geneb.hist.mean <- mean(rna.mat[which(rna.mat$gene_short_name == geneb & rna.mat$broad.histology == hist),'value'])
  genea.mean <- mean(rna.mat[which(rna.mat$gene_short_name == genea),'value'])
  geneb.mean <- mean(rna.mat[which(rna.mat$gene_short_name == geneb),'value'])
  # gene not expressed/ < 1 TPM
  df[i,'Gene1A_not_expressed'] <- ifelse(length(genea.expr) == 0 , NA, genea.expr)
  df[i,'Gene1A_sample'] <- ifelse(is.na(genea.val) || length(genea.val) == 0  , NA, genea.val)
  df[i,'Gene1A_hist_mean'] <- ifelse(is.na(genea.hist.mean) || length(genea.hist.mean) == 0, NA, genea.hist.mean)
  df[i,'Gene1A_mean_dataset'] <- ifelse(is.na(genea.mean) || length(genea.mean) == 0, NA, genea.mean)
  df[i,'Gene1A_expr'] <- ifelse(df[i,'Gene1A_sample'] < df[i,'Gene1A_mean_dataset'],'Decreased',ifelse(df[i,'Gene1A_sample'] == df[i,'Gene1A_mean_dataset'], 'Same','Increased'))
  # gene not expressed/ < 1 TPM
  df[i,'Gene2A_not_expressed'] <- ifelse(length(geneb.expr) == 0 || geneb.val =="NA" , NA, geneb.expr)
  df[i,'Gene2A_sample'] <- ifelse(is.na(geneb.val) || length(geneb.val) == 0 , NA, geneb.val)
  df[i,'Gene2A_hist_mean'] <- ifelse(is.na(geneb.hist.mean) || length(geneb.hist.mean) == 0, NA, geneb.hist.mean)
  df[i,'Gene2A_mean_dataset'] <- ifelse(is.na(geneb.mean) || length(geneb.mean) == 0, NA, geneb.mean)
  df[i,'Gene2A_expr'] <- ifelse(df[i,'Gene2A_sample'] < df[i,'Gene2A_mean_dataset'],'Decreased',ifelse(df[i,'Gene2A_sample'] == df[i,'Gene2A_mean_dataset'], 'Same','Increased'))
}
df[is.na(df)] <- NA
df$Gene1A_not_expressed[is.na(df$Gene1A_not_expressed)] <- "Not Reported"
df$Gene2A_not_expressed[is.na(df$Gene2A_not_expressed)] <- "Not Reported"
#write.table(df, file = paste0(basePD,'/data/processed/Gene_Fusions_exp.txt'), quote = F, sep = "\t", row.names = F)

separate.fusions <- df[(df$Gene1A_not_expressed %in% c(TRUE,"Not Reported") & df$Gene2A_not_expressed %in% c(TRUE,"Not Reported")),]
if(nrow(separate.fusions) > 0){
  print("Fusions to be separated")
  NotExprReported_GeneA=separate.fusions
  df <- df[-which(df$FusionName %in% separate.fusions$FusionName),]
}

fusion_calls_no_intergenic <- df
####complete####





#geneB of intergenic breakpoints
fusion_calls_intergenic<-NA
NotExprReported_GeneB<-NA
if(nrow(fusion_calls[grep ("/",fusion_calls$FusionName),])>0){
df<-fusion_calls[grep ("/",fusion_calls$FusionName),]
df$Gene1B_sample <- NA
df$Gene1B_hist_mean <- NA
df$Gene1B_mean_dataset <- NA
df$Gene1B_expr <- NA
df$Gene2B_sample <- NA
df$Gene2B_hist_mean <- NA
df$Gene2B_mean_dataset <- NA
df$Gene2B_expr <- NA
df$Gene1B_not_expressed <- NA
df$Gene2B_not_expressed <- NA

head(df)

for(i in 1:nrow(df)){
  print(i)
  genea <- df[i,'Gene1B']
  geneb <- df[i,'Gene2B']
  sample <- df[i,'Sample']
  hist <- df[i,'broad.histology']
  hist <- gsub(' [(].*','', hist)
  genea.expr <- unique(rna.mat[which(rna.mat$gene_short_name == genea),'not_expressed'])
  geneb.expr <- unique(rna.mat[which(rna.mat$gene_short_name == geneb),'not_expressed'])
  genea.val <- rna.mat[which(rna.mat$sample_id %in% sample & rna.mat$gene_short_name == genea),'value']
  genea.expr <- genea.val< 1
  geneb.val <- rna.mat[which(rna.mat$sample_id %in% sample & rna.mat$gene_short_name == geneb),'value']
  geneb.expr <- geneb.val< 1
  genea.hist.mean <- mean(rna.mat[which(rna.mat$gene_short_name == genea & rna.mat$broad.histology == hist),'value'])
  geneb.hist.mean <- mean(rna.mat[which(rna.mat$gene_short_name == geneb & rna.mat$broad.histology == hist),'value'])
  genea.mean <- mean(rna.mat[which(rna.mat$gene_short_name == genea),'value'])
  geneb.mean <- mean(rna.mat[which(rna.mat$gene_short_name == geneb),'value'])
  # gene not expressed/ < 1 TPM
  df[i,'Gene1B_not_expressed'] <- ifelse(length(genea.expr) == 0 , NA, genea.expr)
  df[i,'Gene1B_sample'] <- ifelse(is.na(genea.val) || length(genea.val) == 0  , NA, genea.val)
  df[i,'Gene1B_hist_mean'] <- ifelse(is.na(genea.hist.mean) || length(genea.hist.mean) == 0, NA, genea.hist.mean)
  df[i,'Gene1B_mean_dataset'] <- ifelse(is.na(genea.mean) || length(genea.mean) == 0, NA, genea.mean)
  df[i,'Gene1B_expr'] <- ifelse(df[i,'Gene1B_sample'] < df[i,'Gene1B_mean_dataset'],'Decreased',ifelse(df[i,'Gene1B_sample'] == df[i,'Gene1B_mean_dataset'], 'Same','Increased'))
  # gene not expressed/ < 1 TPM
  df[i,'Gene2B_not_expressed'] <- ifelse(length(geneb.expr) == 0, NA, geneb.expr)
  df[i,'Gene2B_sample'] <- ifelse(is.na(geneb.val) || length(geneb.val) == 0 , NA, geneb.val)
  df[i,'Gene2B_hist_mean'] <- ifelse(is.na(geneb.hist.mean) || length(geneb.hist.mean) == 0, NA, geneb.hist.mean)
  df[i,'Gene2B_mean_dataset'] <- ifelse(is.na(geneb.mean) || length(geneb.mean) == 0, NA, geneb.mean)
  df[i,'Gene2B_expr'] <- ifelse(df[i,'Gene2B_sample'] < df[i,'Gene2B_mean_dataset'],'Decreased',ifelse(df[i,'Gene2B_sample'] == df[i,'Gene2B_mean_dataset'], 'Same','Increased'))
}
df[is.na(df)] <- NA
df$Gene1B_not_expressed[is.na(df$Gene1B_not_expressed)] <- "Not Reported"
df$Gene2B_not_expressed[is.na(df$Gene2B_not_expressed)] <- "Not Reported"
#write.table(df, file = paste0(basePD,'/data/processed/Gene_Fusions_exp.txt'), quote = F, sep = "\t", row.names = F)

separate.fusions <- df[(df$Gene1B_not_expressed %in% c(TRUE,"Not Reported") & df$Gene2B_not_expressed %in% c(TRUE,"Not Reported")),]
if(nrow(separate.fusions) > 0){
  print("Fusions to be separated")
  NotExprReported_GeneB=separate.fusions
  df <- df[-which(df$FusionName %in% separate.fusions$FusionName),]
}


fusion_calls_intergenic <-df
fusion_calls<-merge(fusion_calls_no_intergenic,fusion_calls_intergenic[,c("Sample","FusionName","Gene1B_sample","Gene2B_sample","Gene1B_mean_dataset","Gene2B_mean_dataset","Gene1B_hist_mean","Gene2B_hist_mean","Gene1B_expr","Gene2B_expr","Gene1B_not_expressed","Gene2B_not_expressed")],by=c("Sample","FusionName"),all=T)
}else{
  fusion_calls<-fusion_calls_no_intergenic}

fusion_calls<-unique(fusion_calls)

expfilt_fusion_calls<-list("NotexpressedGeneA"=NotExprReported_GeneA,"NotExprReported_GeneB"=NotExprReported_GeneB,fusion_calls,Expressed_genes_fusion=fusion_calls)

return(expfilt_fusion_calls)

}




