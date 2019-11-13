#' theme_publication for plots

#' @param base_size size of font for plot default "15"
#' @param base_family font to be used default "Helvetica"
#' Publication quality images for summary plots

theme_Publication <- function(base_size=15, base_family="Helvetica") {
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