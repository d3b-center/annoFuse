#' Function to plot recurrent fused genes

#' @param standardFusioncalls A dataframe from star fusion or arriba standardized to run through the filtering steps
#' @param groupby column name with grouping variables
#' @param plotn top n recurrent fusions to plot
#' @param countID column name to count recurrent fusions SampleID/ParticipantID/tumorID
#' @param palette_rec colors for grouping variables
#'
#' @export
#'
#' @return A ggplot object containing an overview on the recurrent fused genes
#'
#' @examples
#' out_annofuse <- system.file("extdata", "PutativeDriverAnnoFuse_test_v16.tsv", package = "annoFuse")
#' sfc <- read.delim(out_annofuse)
#' plot_recurrent_genes(sfc, groupby = "Fusion_Type", countID = "Sample")
plot_recurrent_genes <- function(standardFusioncalls,
                                 groupby,
                                 plotn = 20,
                                 countID,
                                 palette_rec = NULL) {
  
  standardFusioncalls <- .check_annoFuse_calls(standardFusioncalls)
  stopifnot(is.character(groupby))
  stopifnot(is.numeric(plotn))
  stopifnot(is.character(countID))
  
  stopifnot(all(c(groupby, countID) %in% colnames(standardFusioncalls)))

  # inframe fusions only
  # standardFusioncalls<-unique(standardFusioncalls) %>% dplyr::filter(.data$Fusion_Type=="in-frame")

  # remove geneA==geneB or intergenic
  # standardFusioncalls<-standardFusioncalls[-which(standardFusioncalls$Gene1A==standardFusioncalls$Gene2A|standardFusioncalls$Gene1A==standardFusioncalls$Gene2B|standardFusioncalls$Gene2A==standardFusioncalls$Gene1B|standardFusioncalls$Gene1A==standardFusioncalls$Gene1B|standardFusioncalls$Gene2B==standardFusioncalls$Gene2A),]

  # gene1A recurrent
  rec_gene1A <- standardFusioncalls %>%
    dplyr::select("Gene1A", !!as.name(groupby), !!as.name(countID)) %>%
    unique() %>%
    group_by(.data$Gene1A, !!as.name(groupby)) %>%
    dplyr::select(-!!as.name(countID)) %>%
    mutate(count = n()) %>%
    unique() %>%
    as.data.frame()

  colnames(rec_gene1A) <- c("Gene", as.name(groupby), "count")

  # gene1B recurrent
  rec_gene1B <- standardFusioncalls %>%
    dplyr::select("Gene1B", !!as.name(groupby), !!as.name(countID)) %>%
    unique() %>%
    group_by(.data$Gene1B, !!as.name(groupby)) %>%
    dplyr::select(-!!as.name(countID)) %>%
    mutate(count = n()) %>%
    unique() %>%
    as.data.frame()

  colnames(rec_gene1B) <- c("Gene", as.name(groupby), "count")
  rec_gene <- rbind(rec_gene1A, rec_gene1B)

  rec_gene <- utils::head(rec_gene[order(rec_gene$count, decreasing = TRUE), ], plotn)
  rec_gene$Gene <- factor(rec_gene$Gene, levels = unique(rec_gene$Gene), ordered = TRUE)

  if (!is.null(palette_rec)) {
    # provided palette is rownames=groupby values and color in column "color"
    palette_2 <- palette_rec
  } else {
    # pallet to match recurrent fusion and recurrent fused genes
    n <- length(levels(as.factor(standardFusioncalls[, groupby])))
    palette <- grDevices::rainbow(n)
    names(palette) <- levels(as.factor(standardFusioncalls[, groupby]))
    colScale <- scale_colour_manual(name = as.name(groupby), values = palette)
    # to match with recurrent fusion
    palette_2 <- palette[which(names(palette) %in% rec_gene[, groupby])]
  }

  rec_genes <- ggplot(rec_gene) +
    geom_col(aes(x = .data$Gene, y = .data$count, fill = !!as.name(groupby)), alpha = 0.75) +
    ylab("Number of patients") +
    guides(color = FALSE, alpha = FALSE) +
    xlab(NULL) +
    scale_y_continuous(breaks = seq(0, 150, by = 10)) +
    rotate() +
    theme_publication() +
    theme(legend.title = element_blank(), axis.text.y = element_text(face = "italic", angle = 0, hjust = 1, size = 12)) +
    scale_fill_manual(name = as.name(groupby), values = palette_2) +
    scale_x_discrete(limits = rev(levels(rec_gene$Gene)))

  return(rec_genes)
}
