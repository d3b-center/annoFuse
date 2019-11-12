#' Function to annotate fusion calls

#' @param standardFusionCalls A dataframe from star fusion or arriba standardized to run through the filtering steps
#' @param outputpdffile Filename to plot image
#' @param groupby column name with grouping variables
#' @return summary pdf


plotSummary<-function(standardFusionCalls=standardFusionCalls,outputpdffile=outputpdffile,groupby=groupby){

  theme_Publication <- function(base_size=15, base_family="Helvetica") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size=base_size, base_family=base_family)
      + theme(plot.title = element_text(face = "bold",
                                        size = 15, hjust = 0.5),
              text = element_text(family="Helvetica"),
              panel.background = element_rect(colour = NA),
              plot.background = element_rect(colour = NA),
              panel.border = element_rect(colour = NA),
              axis.title = element_text(face = "bold",size = 12),
              axis.title.y = element_text(angle=90,vjust =2,size=12),
              axis.title.x = element_text(size=12),
              axis.text.x = element_text(size=12,color="black",face="bold",hjust = 1),
              axis.text.y = element_text(size=12,color="black",face="bold",hjust = 1),
              axis.line = element_line(colour="black",size=1),
              axis.ticks = element_line(),
              panel.grid.major = element_line(colour="#f0f0f0"),
              panel.grid.minor = element_blank(),
              legend.text = element_text(size=12),
              legend.key = element_rect(colour = NA),
              legend.position = "right",
              legend.direction = "vertical",
              legend.key.size= unit(0.5, "cm"),
              legend.spacing  = unit(1, "mm"),
              legend.title = element_text(family="Helvetica",face="italic",size=rel(0.7)),
              plot.margin=unit(c(10,10,1,10),"mm"),
              strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
              strip.text = element_text(face="bold")
      ))

  }

if(missing(groupby)){
  groupby="FusionCalls"
}
  # Inter and intrachromosomal plot
fusion_calls_interchromosomal<-standardFusionCalls %>% dplyr::filter(grepl("INTERCHROMOSOMAL",annots))
fusion_calls_intrachromosomal<-standardFusionCalls %>% dplyr::filter(grepl("INTRACHROMOSOMAL",annots))

fusion_chrom<-rbind(data.frame(cbind(fusion_calls_interchromosomal,"Distance"=rep("Interchromosomal",nrow(fusion_calls_interchromosomal)))),cbind(fusion_calls_intrachromosomal,"Distance"=rep("Intrachromosomal",nrow(fusion_calls_intrachromosomal)))) %>% unique()

p1<-ggplot(fusion_chrom,aes(x=!!as.name(groupby),fill=Distance,alpha=0.75))+geom_bar()+theme_Publication()+theme(legend.position = "top")+xlab(as.name(groupby))+ylab("Count")+guides(alpha=FALSE)+theme(axis.text.x  = element_text(angle = 45))

  # frame information
p2<-ggplot(standardFusionCalls,aes(fill=Fusion_Type,x=Caller,alpha = 0.75))+geom_bar(aes(y=log2(..count..)))+rotate()+xlab("Caller")+ylab("Count (log2)")+theme_Publication()+theme(legend.position = "top")+guides(alpha=FALSE)+theme(axis.text.x  = element_text(angle = 0,hjust = 1,size=12),axis.text.y = element_text(angle=0,vjust =2,size=12))+labs(fill = "Frame")


# keep fusion if atleast 1 is protein-coding
#BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86

genes<-as.data.frame(genes(edb))

# keep if atleast 1 gene is protein coding gene in fusions
fusion_protein_coding_gene_df <- standardFusionCalls %>%
  # We want to keep track of the gene symbols for each sample-fusion pair
  dplyr::select(Sample, FusionName, Gene1A, Gene1B, Gene2A, Gene2B,!!(as.name(groupby))) %>%
  # We want a single column that contains the gene symbols
  tidyr::gather(Gene1A, Gene1B, Gene2A, Gene2B,
                key = gene_position, value = GeneSymbol) %>%
  # Remove columns without gene symbols
  dplyr::filter(GeneSymbol != "") %>%
  dplyr::arrange(Sample, FusionName) %>%
  # Retain only distinct rows
  dplyr::distinct() %>% left_join(genes,by=c("GeneSymbol"="gene_name"))


# bubble plot for genes in biotype fromensemble
fusion_type_gene_df <- fusion_protein_coding_gene_df %>%
  # We want to keep track of the gene symbols for each sample-fusion pair
  dplyr::select(Sample, FusionName, gene_position, gene_biotype,!!(as.name(groupby))) %>%
  unique() %>%
  group_by(!!(as.name(groupby)),gene_biotype,gene_position) %>%
  dplyr::summarise(Type.ct = n())

fusion_type_gene_df<-fusion_type_gene_df[!is.na(fusion_type_gene_df$gene_biotype),]


p3<-ggplot(fusion_type_gene_df,aes(x=gene_biotype,y=gene_position,size=Type.ct,fill=gene_position,color=gene_position))+geom_point(shape=21)+scale_color_brewer(palette="Pastel2")+scale_fill_brewer(palette="Pastel2")+guides( size = FALSE,color=FALSE)+scale_size_continuous(range = c(1,12))+xlab("GeneBiotype")+ylab("GenePosition")+theme_Publication()+theme(axis.text.x  = element_text(angle = 45,hjust = 1,size=12),axis.text.y = element_text(angle=0,vjust =2,size=12))+theme(legend.position="top")+labs(fill="Gene Position")



# Kinase groups
kinaseGeneList<-read_tsv(system.file("extdata", "Kincat_Hsap.08.02.txt", package = "annoFuse"))
kinase_fusion<-fusion_protein_coding_gene_df %>%
  left_join(kinaseGeneList,by=c("GeneSymbol"="Name")) %>%
  dplyr::filter(!is.na(Family)) %>%
  dplyr::select("gene_position","Group","Sample") %>%
  unique() %>%
  group_by(Group,gene_position) %>%
  dplyr::summarise(Type.ct = n()) %>%
  unique()

p4<-ggplot(kinase_fusion,aes(x=Group,y=Type.ct,fill=gene_position,alpha=0.75))+geom_col()+scale_color_brewer(palette="Pastel2")+scale_fill_brewer(palette="Pastel2")+theme_Publication()+ylab("Count")+theme(legend.position = "right")+guides(alpha=FALSE)+theme(axis.text.x  = element_text(angle = 45,hjust = 1,size=12),axis.text.y = element_text(angle=0,vjust =2,size=12),legend.position = "top")+labs(fill = "Gene Position")+guides(fill=FALSE)


# bubble plot for genes in kinase/oncoplot etc
fusion_type_gene_df <- standardFusionCalls %>%
  # We want to keep track of the gene symbols for each sample-fusion pair
  dplyr::select(Sample, FusionName, Gene1A_anno, Gene1B_anno, Gene2A_anno, Gene2B_anno,!!(as.name(groupby))) %>%
  unique() %>%
  # We want a single column that contains the gene symbols
  tidyr::gather(Gene1A_anno, Gene1B_anno, Gene2A_anno, Gene2B_anno,
                key = gene_position, value = Annot) %>%
  mutate(gene_position=gsub("_anno","",gene_position)) %>%
  # Remove columns without gene symbols
  dplyr::filter(Annot != "") %>%
  mutate(Annotation=gsub(" ","",.$Annot)) %>%
  mutate(Annotation = strsplit(as.character(Annotation), ",")) %>%
  unnest(Annotation)  %>%
  group_by(!!as.name(groupby),Annotation,gene_position) %>%
  dplyr::summarise(Annot.ct = n()) %>%
  # To only plot the broader transcription factor annotation
  dplyr::filter(!Annotation %in% c("curatedTF","predictedTF"))

#To-Do rename in original table
fusion_type_gene_df$Annotation[grep("onco",fusion_type_gene_df$Annotation)]<-"Oncogene"
fusion_type_gene_df$Annotation[grep("tsgs",fusion_type_gene_df$Annotation)]<-"TumorSuppressorGene"
fusion_type_gene_df$Annotation[grep("kinase",fusion_type_gene_df$Annotation)]<-"Kinase"

p5<-ggplot(fusion_type_gene_df,aes(x=!!as.name(groupby),y=Annot.ct,size=Annot.ct,fill=Annotation,color=Annotation,alpha=0.75))+geom_point(shape=21)+scale_color_brewer(palette="Pastel1")+scale_fill_brewer(palette="Pastel1")+guides( size = FALSE,alpha=FALSE)+scale_size_continuous(range = c(1,12))+ylab("Count")+xlab(groupby)+theme_Publication()+theme(legend.position = "top")+theme(axis.text.x  = element_text(angle = 45,hjust = 1,size=12),axis.text.y = element_text(angle=0,vjust =2,size=12))+theme(plot.margin=unit(c(0,10,0,15),"mm"))+facet_wrap(~gene_position)

if (!missing(outputpdffile)){
ggarrange(p1,p2,p3,p4,p5,labels = c("A","B","C","D","E"),heights=c(5,5,6),widths=c(2,1) ,nrow=3,ncol=2,font.label = list(size=30)) %>% ggexport(filename = outputpdffile,width = 20,height = 20)
} else {
  summary<-ggarrange(p1,p2,p3,p4,p5,heights=c(5,5,6,6,6),widths=4,ncol=1)
}

return (summary)

}

