#' Function to plot breakpoints of a gene in a sample/cohort

#' @param domainDataFrame A dataframe from star fusion or arriba standardized to run through the filtering steps and get_Pfam_domain() function to get pfam domain annotation per gene
#' @param exons exon level information from gtf
#' @param geneposition Left/Right position of gene
#'
#' @export
#'
#' @return ggplot of gene and breakpoints
#'
#' @examples
#' # TODOTODO
plot_breakpoints <- function(domainDataFrame = domainDataFrame,
                             exons = exons,
                             geneposition = geneposition) {

  # find unique
  if (geneposition == "Left") {
    geneExons <- exons[which(exons$geneName %in% unique(domainDataFrame$Gene1A)), ]
    # merge domain and exon annotation to plot
    geneExons <- geneExons %>% left_join(domainDataFrame, by = c("geneName" = "Gene1A"))
    # rename transcript column for gene
    geneExons[grep("gene", geneExons$type), "transcript"] <- geneExons[grep("gene", geneExons$type), "geneName"]
    # get gene name
    gene <- unique(geneExons$transcript[grep("gene", geneExons$type)])
    # find unique transcripts
    uniqtranscript <- unique(geneExons$transcript)
    # order transcripts for plots
    uniqtranscript <- c(gene, uniqtranscript[-grep(gene, uniqtranscript)])
    geneExons$transcript <- factor(geneExons$transcript, levels = uniqtranscript[length(uniqtranscript):1])
    # basic plot for exon and gene
    p1 <- ggplot(geneExons, aes(y = geneExons$start, x = geneExons$transcript)) +
      geom_hline(yintercept = as.numeric(domainDataFrame$LeftBreakpoint), linetype = "dashed", color = "violet") +
      geom_linerange(ymin = geneExons$start, ymax = geneExons$end, size = 3) +
      geom_linerange(ymin = geneExons$start, ymax = geneExons$end, size = 6, aes(col = geneExons$DESC, x = gene), color = "grey") +
      coord_flip() +
      theme_publication() +
      xlab("Transcript") +
      ylab("Coordinate") +
      ggtitle(paste0(geneposition, " fused Gene:", gene)) +
      labs(fill = "Domain")
  }
  if (geneposition == "Right") {
    geneExons <- exons[which(exons$geneName %in% unique(domainDataFrame$Gene1B)), ]
    # merge domain and exon annotation to plot
    geneExons <- geneExons %>% left_join(domainDataFrame, by = c("geneName" = "Gene1B"))
    # rename transcript column for gene
    geneExons[grep("gene", geneExons$type), "transcript"] <- geneExons[grep("gene", geneExons$type), "geneName"]
    # get gene name
    gene <- unique(geneExons$transcript[grep("gene", geneExons$type)])
    # find unique transcripts
    uniqtranscript <- unique(geneExons$transcript)
    # order transcripts for plots
    uniqtranscript <- c(gene, uniqtranscript[-grep(gene, uniqtranscript)])
    geneExons$transcript <- factor(geneExons$transcript, levels = uniqtranscript[length(uniqtranscript):1])
    # basic plot for exon and gene
    p1 <- ggplot(geneExons, aes(y = geneExons$end, x = geneExons$transcript)) +
      geom_hline(yintercept = as.numeric(domainDataFrame$RightBreakpoint), linetype = "dashed", color = "violet") +
      geom_linerange(ymin = geneExons$start, ymax = geneExons$end, size = 3) +
      geom_linerange(ymin = geneExons$start, ymax = geneExons$end, size = 6, aes(col = geneExons$DESC, x = gene), color = "grey") +
      coord_flip() +
      theme_publication() +
      xlab("Transcript") +
      ylab("Coordinate") +
      ggtitle(paste0(geneposition, " fused Gene:", gene))
  }

  # add domain level information
  if (any(!is.na(geneExons$DESC))) {
    p1 <- p1 + geom_linerange(ymin = geneExons$domain_start, ymax = geneExons$domain_end, size = 6, aes(col = geneExons$DESC, x = gene)) + labs(color = "Domain")
  }

  # add direction of strand
  if (any(geneExons$strand.x == "+")) {
    p2 <- p1 + geom_segment(aes(x = gene, xend = gene, y = min(p1$data$start), yend = max(p1$data$end), alpha = 0.1), arrow = (arrow(length = unit(0.2, "cm")))) + geom_line(alpha = 0.8) + guides(alpha = FALSE)
  } else {
    p2 <- p1 + geom_segment(aes(x = gene, xend = gene, y = max(p1$data$end), yend = min(p1$data$start), alpha = 0.1), arrow = (arrow(length = unit(0.2, "cm")))) + geom_line(alpha = 0.8) + guides(alpha = FALSE)
  }

  return(p2)
}
