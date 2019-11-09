#' Function to plot recurrent fused genes

#' @param standardFusioncalls A dataframe from star fusion or arriba standardized to run through the filtering steps
#' @param groupby column name with grouping variables
#' @param plotn top n recurrent fusions to plot
#' @param countID column name to count recurrent fusions SampleID/ParticipantID/tumorID
#' @return recurrent fused genes
#'


plotRecurrentFusions<-function(standardFusioncalls=standardFusioncalls,groupby=groupby,plotn=plotn,countID=countID){
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

# in-frame fusions only  
  standardFusioncalls<-unique(fusion_calls) %>% dplyr::filter(Fusion_Type=="in-frame")  
# remove geneA==geneB or intergenic
  standardFusioncalls<-standardFusioncalls[-which(standardFusioncalls$Gene1A==standardFusioncalls$Gene2A|standardFusioncalls$Gene1A==standardFusioncalls$Gene2B|standardFusioncalls$Gene2A==standardFusioncalls$Gene1B|standardFusioncalls$Gene1A==standardFusioncalls$Gene1B|standardFusioncalls$Gene2B==standardFusioncalls$Gene2A),]
  
  # gather recurrent per group
  rec_fusions<-standardFusioncalls %>%
    as.data.frame() %>%
    dplyr::select("FusionName",!!as.name(groupby),!!as.name(countID)) %>%
    unique() %>%
    group_by(FusionName,!!as.name(groupby)) %>%
    dplyr::select(-!!as.name(countID)) %>%
    mutate(count=n()) %>%
    unique() %>%
    as.data.frame()
  
  #Top n recurrent fusion per group
  rec_fusions<-head(rec_fusions[order(rec_fusions$count,decreasing = TRUE),],plotn)
  
  #pallet to match recurrent fusion and recurrent fused genes
  n <- length(levels(as.factor(standardFusioncalls[,groupby])))
  palette <- rainbow(n)
  names(palette)<-levels(as.factor(standardFusioncalls[,groupby]))
  colScale <- scale_colour_manual(name = as.name(groupby),values = palette)
  
  palette_1<-palette[which(names(palette) %in% rec_fusions[,groupby] )]
  
  rec_fusions$FusionName<-factor(rec_fusions$FusionName,levels=unique(rec_fusions$FusionName),ordered = TRUE)
  rec_fusions<-ggplot(rec_fusions)+geom_col(aes(x=FusionName,y=count,fill=!!as.name(groupby)),alpha=0.75)+
    guides(alpha = FALSE) + ylab('Number of patients')+ xlab(NULL)+
    guides(color = F,alpha=FALSE)+
    scale_y_continuous(breaks = seq(0, 200, by = 20))+ggpubr::rotate()+scale_fill_manual(name = as.name(groupby),values = palette_1)+theme_Publication()+theme(legend.title = element_blank(),axis.text.y  = element_text(face="italic",angle = 0,hjust = 1,size = 12),legend.position = "bottom")+scale_x_discrete(limits = rev(levels(rec_fusions$FusionName)))

  return(rec_fusions)

}
