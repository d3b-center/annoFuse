#' Function to plot expression as percentile
#'
#' @param expressiongene matrix with 1 row with expression for gene and N columns for samples to be used for analysis
#' @param geneSymbol geneSymbol
#' @param standardFusionCalls standardized fusion calls
#' @return ggplot grob for plotting median and each value per fusion as x axis




plotPercentilePerGene<-function(expressiongene=expressiongene,geneSymbol=geneSymbol,standardFusionCalls=standardFusionCalls){

  expressionMatrixPercentRank<-percent_rank(expressiongene[,-which(colnames(expressiongene) %in% c("gene_id","geneSymbol","EnsembleID"))])


  # long format pbta
  expressionMatrixPercentRank_GOI_long<-expressionMatrixPercentRank %>% as.data.frame(stringsAsFactors=FALSE) %>%
    t() %>%
    as.data.frame(stringsAsFactors=FALSE) %>%
    mutate(gene_symbol=geneSymbol) %>%
    reshape2::melt() %>%
    mutate (variable=as.character(.$variable)) %>%
    left_join(clinical,by=c("variable"="Kids_First_Biospecimen_ID"))


  checkFusion<-paste0("--",geneSymbol,"$|^",geneSymbol,"--")

  # get gene of interest Fusion status
  expressionMatrixPercentRank_GOI_long <- expressionMatrixPercentRank_GOI_long %>%
    left_join(standardFusionCalls[grep(checkFusion,standardFusionCalls$FusionName),c("Sample","FusionName")],by=c("variable"="Sample"))

  # add FusionName Normal to GTEX samples
  expressionMatrixPercentRank_GOI_long$FusionName <-ifelse(is.na(expressionMatrixPercentRank_GOI_long$FusionName) & grepl("^SRR",expressionMatrixPercentRank_GOI_long$variable),"Normal",expressionMatrixPercentRank_GOI_long$FusionName)

  #remove NA PBTA samples
  expressionMatrixPercentRank_GOI_long <- expressionMatrixPercentRank_GOI_long[!is.na(expressionMatrixPercentRank_GOI_long$FusionName),]

  # expressionMatrixPercentRank_GOI_long<-rbind(expressionMatrixPercentRank_GOI_long,expressionMatrixPercentRank_GTEx_long)

  # get median for Fusion=Yes and Fusion=No
  median_per_gene_fusion<-expressionMatrixPercentRank_GOI_long %>%
    group_by(FusionName) %>%
    summarise(value = list(enframe(median(value))),count=n()) %>%
    unnest %>% rename(median=value) %>%
    as.data.frame()

  median_per_gene_fusion<-median_per_gene_fusion[order(median_per_gene_fusion$median),]
  median_per_gene_fusion$FusionName <- factor(median_per_gene_fusion$FusionName, levels = c(as.character(median_per_gene_fusion$FusionName[grep("Normal",median_per_gene_fusion$FusionName)]),as.character(median_per_gene_fusion$FusionName[-grep("Normal",median_per_gene_fusion$FusionName)])) )


  #merge median
  expressionMatrixPercentRank_GOI_long<-expressionMatrixPercentRank_GOI_long %>% left_join(median_per_gene_fusion,by=c("FusionName"="FusionName"))

  p_value<-ggplot(expressionMatrixPercentRank_GOI_long,aes(y=value))+geom_point(aes(y=value,x=FusionName,color=value,size=12))+scale_colour_gradient2(low = "blue", mid = "white",high = "red", midpoint = 0.5)+theme_Publication()+theme(axis.text.x = element_text(angle = 60))+ggtitle(paste0("Percentile expression of ",geneSymbol, " in fusion compared to GTEx normal"))+guides( size = FALSE)

  median_per_gene_fusion$count[grep("Normal",median_per_gene_fusion$FusionName)]<-1

  p_median<-ggplot(median_per_gene_fusion,aes(y=median))+geom_point(aes(y=median,x=FusionName,color=median,size=12))+scale_colour_gradient2(low = "blue", mid = "white",high = "red", midpoint = 0.5)+theme_Publication()+theme(axis.text.x = element_text(angle = 60))+ggtitle(paste0("Percentile expression of ",geneSymbol, " in fusion compared to GTEx normal"))+labs(size="Count (N of normal not shown here)")+guides( size = FALSE)

  p<-list("median"=p_median,"value"=p_value)

  return(p)

}

theme_Publication <- function(base_size=25, base_family="Helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(0.75), hjust = 0.5),
            text = element_text(family="Helvetica"),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(0.75)),
            axis.title.y = element_text(angle=90,vjust =2,size=rel(0.75)),
            axis.title.x = element_text(size=rel(0.75)),
            axis.text = element_text(size=12,color="black",face="bold",angle = 60,hjust = 1),
            axis.line = element_line(colour="black",size=1),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.text = element_text(size=rel(0.7)),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size= unit(0.7, "cm"),
            legend.spacing  = unit(1, "mm"),
            legend.title = element_text(family="Helvetica",face="italic",size=rel(0.7)),
            plot.margin=unit(c(10,10,1,10),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))

}
