#' Function to plot recurrent fused genes

#' @param standardFusioncalls A dataframe from star fusion or arriba standardized to run through the filtering steps
#' @param groupby column name with grouping variables
#' @param plotn top n recurrent fusions to plot
#' @param countID column name to count recurrent fusions SampleID/ParticipantID/tumorID
#' @return recurrent fused genes
#'


plotRecurrentGenes<-function(standardFusioncalls=standardFusioncalls,groupby=groupby,plotn=plotn,countID=countID){
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

# inframe fusions only
  standardFusioncalls<-unique(standardFusioncalls) %>% dplyr::filter(Fusion_Type=="in-frame")
  
# remove geneA==geneB or intergenic
standardFusioncalls<-standardFusioncalls[-which(standardFusioncalls$Gene1A==standardFusioncalls$Gene2A|standardFusioncalls$Gene1A==standardFusioncalls$Gene2B|standardFusioncalls$Gene2A==standardFusioncalls$Gene1B|standardFusioncalls$Gene1A==standardFusioncalls$Gene1B|standardFusioncalls$Gene2B==standardFusioncalls$Gene2A),]

#gene1A recurrent
rec_gene1A<-standardFusioncalls %>%
  dplyr::select("Gene1A",!!as.name(groupby),!!as.name(countID))%>%
  unique() %>%
  group_by(Gene1A,!!as.name(groupby)) %>%
  dplyr::select(-Sample) %>%
  mutate(count=n())%>% unique() %>% as.data.frame()

colnames(rec_gene1A)<-c("Gene",as.name(groupby),"count")

# gene1B recurrent
rec_gene1B<-standardFusioncalls %>%
  dplyr::select("Gene1B",!!as.name(groupby),!!as.name(countID))%>%
  unique() %>%
  group_by(Gene1B,!!as.name(groupby)) %>%
  dplyr::select(-!!as.name(countID)) %>%
  mutate(count=n())%>% unique() %>% as.data.frame()

colnames(rec_gene1B)<-c("Gene",as.name(groupby),"count")
rec_gene<-rbind(rec_gene1A,rec_gene1B)

rec_gene<-head(rec_gene[order(rec_gene$count,decreasing = TRUE),],plotn)
rec_gene$Gene<-factor(rec_gene$Gene,levels = unique(rec_gene$Gene),ordered = TRUE)

# to match with recurrent fusion
palette_2<-palette[which(names(palette) %in% rec_gene[,groupby] )]

rec_genes<-ggplot(rec_gene)+geom_col(aes(x=Gene,y=count,fill=!!as.name(groupby)),alpha=0.75)+ylab('Number of patients') +
  guides(color = F,alpha=FALSE) + xlab(NULL)+
  scale_y_continuous(breaks = seq(0, 150, by =10))+rotate()+theme_Publication()+theme(legend.title = element_blank(),axis.text.y  = element_text(face="italic",angle = 0,hjust = 1,size = 12))+scale_fill_manual(values = palette_2)+scale_x_discrete(limits = rev(levels(rec_gene$Gene)))

return(rec_genes)

}
