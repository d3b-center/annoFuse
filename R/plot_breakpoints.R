#' Function to plot breakpoints of a gene in a sample/cohort

#' @param domainDataFrame A dataframe from star fusion or arriba standardized to run through the filtering steps and get_Pfam_domain() function to get pfam domain annotation per gene
#' @param exons exon level information from gtf availbale in the package exonsToPlot.RDS
#' @param geneposition Position of gene in the fusion allowed values are Left (Gene1A) or Right (Gene1B)
#' @param sampleid id (if multiple provide a character vector) from Sample column of standardized fusion file format
#' @param fusionname value from FusionName column of standardized fusion file format (Gene1A--Gene1B)
#' @param leftBreakpoint If a specific breakpoint needs to be plotted then provide genomic coordinate as LeftBreakpoint
#' @param rightBreakpoint  If a specific breakpoint needs to be plotted then provide genomic coordinate as RightBreakpoint
#' @param base_size base size for text in plot ; default is 12
#'
#' @export
#'
#' @return ggplot of gene and breakpoints
#'
#' @examples
#' out_annofuse <- system.file("extdata", "PutativeDriverAnnoFuse_test_v14.tsv", package = "annoFuse")
#' sfc <- read.delim(out_annofuse)
#' exons <- readRDS(system.file("extdata", "exonsToPlot.RDS", package = "annoFuse"))
#' bioMartDataPfam <- 
#'   readRDS(system.file("extdata", "pfamDataBioMart.RDS", package = "annoFuse"))
#' domainDataFrame <- get_Pfam_domain(standardFusioncalls = sfc, 
#'                                    bioMartDataPfam = bioMartDataPfam, 
#'                                    keepPartialAnno = TRUE)
#' left <- plot_breakpoints(sampleid = "BS_044XZ8ST", 
#'                          domainDataFrame = domainDataFrame, 
#'                          exons = exons, 
#'                          geneposition = "Left", 
#'                          fusionname = "ANTXR1--BRAF")
#' right <- plot_breakpoints(sampleid = "BS_044XZ8ST", 
#'                           domainDataFrame = domainDataFrame, 
#'                           exons = exons, 
#'                           geneposition = "Right", 
#'                           fusionname = "ANTXR1--BRAF")
#' ggpubr::ggarrange(left, right, align = "h")
plot_breakpoints <- function(domainDataFrame = NULL,
                             exons,
                             geneposition = c("Left", "Right"),
                             sampleid = NULL,
                             fusionname = NULL,
                             leftBreakpoint = NULL,
                             rightBreakpoint = NULL,
                             base_size = 12) {
  
  geneposition <- match.arg(geneposition, c("Left", "Right"))
  stopifnot(is(exons, "data.frame"))
  # stopifnot(is.character(fusionname))
  stopifnot(is.numeric(base_size))
  
  if (is.null(domainDataFrame)) {
    stop("domainDataFrame not provide; please provide output from get_Pfam_domain() ")
  } else {
    if (is.null(fusionname)) {
      warning("FusionName not provided; using first row of domainDataFrame")
    } else {
      # subset the dataframe with domain information to only fusionname
      domainDataFrame$Gene1B <- domainDataFrame$Gene1B[which(domainDataFrame$Gene1B$FusionName == fusionname), ]
      domainDataFrame$Gene1A <- domainDataFrame$Gene1A[which(domainDataFrame$Gene1A$FusionName == fusionname), ]
    }

    if (nrow(domainDataFrame$Gene1B) == 0 | nrow(domainDataFrame$Gene1A) == 0) {
      stop("domainDataFrame is empty after selecting for fusionname ")
    }

    if (!is.null(sampleid)) {
      domainDataFrame$Gene1B <- domainDataFrame$Gene1B[which(domainDataFrame$Gene1B$Sample %in% sampleid), ]
      domainDataFrame$Gene1A <- domainDataFrame$Gene1A[which(domainDataFrame$Gene1A$Sample %in% sampleid), ]
    }

    if (!is.null(leftBreakpoint)) {
      domainDataFrame$Gene1A <- domainDataFrame$Gene1A[which(domainDataFrame$Gene1A$LeftBreakpoint == leftBreakpoint), ]
    }
    if (!is.null(rightBreakpoint)) {
      domainDataFrame$Gene1B <- domainDataFrame$Gene1B[which(domainDataFrame$Gene1B$RightBreakpoint == rightBreakpoint), ]
    }
  }

  if (base_size != 12) {
    base_size <- base_size
  }

  # find unique
  if (geneposition == "Left") {
    geneExons <- exons[which(exons$geneName %in% unique(domainDataFrame$Gene1A$Gene1A)), ]
    # merge domain and exon annotation to plot
    geneExons <- geneExons %>% left_join(domainDataFrame$Gene1A, by = c("geneName" = "Gene1A"))
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
    p1 <- ggplot(geneExons, aes(y = .data$start, x = .data$transcript)) +
      geom_hline(yintercept = as.numeric(domainDataFrame$Gene1A$LeftBreakpoint), linetype = "dashed", color = "violet") +
      geom_linerange(ymin = geneExons$start, ymax = geneExons$end, size = 3) +
      geom_linerange(ymin = geneExons$start, ymax = geneExons$end, size = 6, aes(col = .data$DESC, x = gene), color = "grey") +
      coord_flip() +
      theme_publication(base_size = base_size) +
      xlab("Transcript") +
      ylab("Coordinate") +
      ggtitle(paste0(geneposition, " fused Gene:", gene)) +
      labs(fill = "Domain")
  }

  if (geneposition == "Right") {
    geneExons <- exons[which(exons$geneName %in% unique(domainDataFrame$Gene1B$Gene1B)), ]
    # merge domain and exon annotation to plot
    geneExons <- geneExons %>% left_join(domainDataFrame$Gene1B, by = c("geneName" = "Gene1B"))
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
    p1 <- ggplot(geneExons, aes(y = .data$end, x = .data$transcript)) +
      geom_hline(yintercept = as.numeric(domainDataFrame$Gene1B$RightBreakpoint), linetype = "dashed", color = "violet") +
      geom_linerange(ymin = geneExons$start, ymax = geneExons$end, size = 3) +
      geom_linerange(ymin = geneExons$start, ymax = geneExons$end, size = 6, aes(col = .data$DESC, x = gene), color = "grey") +
      coord_flip() +
      theme_publication(base_size = base_size) +
      xlab("Transcript") +
      ylab("Coordinate") +
      ggtitle(paste0(geneposition, " fused Gene:", gene))
  }

  # add domain level information
  if (!all(is.na(geneExons$DESC))) {
    p1 <- p1 + geom_linerange(ymin = geneExons$domain_start, ymax = geneExons$domain_end, size = 6, aes(col = .data$DESC, x = gene)) + labs(color = "Domain") + theme(legend.position = "bottom") +
      ggplot2::scale_colour_discrete(na.translate = F)
  }

  # add direction of strand
  if (any(geneExons$strand.x == "+")) {
    p2 <- p1 + geom_segment(aes(x = gene, xend = gene, y = min(p1$data$start), yend = max(p1$data$end), alpha = 0.1), arrow = (arrow(length = unit(0.2, "cm")))) + geom_line(alpha = 0.8) + guides(alpha = FALSE)
  } else {
    p2 <- p1 + geom_segment(aes(x = gene, xend = gene, y = max(p1$data$end), yend = min(p1$data$start), alpha = 0.1), arrow = (arrow(length = unit(0.2, "cm")))) + geom_line(alpha = 0.8) + guides(alpha = FALSE)
  }

  return(p2)
}
