#' Function to plot recurrent fused genes

#' @param standardFusioncalls A dataframe from star fusion or arriba standardized to run through the filtering steps
#'
#' @param groupby column name with grouping variables
#' @param plotn top n recurrent fusions to plot
#' @param countID column name to count recurrent fusions SampleID/ParticipantID/tumorID
#' @param palette_rec colors for grouping variables
#' @param base_size Numeric, size of font for plot
#'
#' @export
#'
#' @return A ggplot object containing an overview on the recurrent fused genes
#'
#' @examples
#' out_annofuse <- system.file("extdata", "PutativeDriverAnnoFuse.tsv", package = "annoFuse")
#' sfc <- read.delim(out_annofuse, stringsAsFactors = FALSE)
#' # keep only in-frame and fusions where both breakpoints are within genes
#' sfc <- as.data.frame(
#'   sfc[ which(sfc$Fusion_Type == "in-frame" && sfc$BreakpointLocation == "Genic"),])
#' plot_recurrent_genes(sfc, groupby = "broad_histology", countID = "Kids_First_Participant_ID")
plot_recurrent_genes <- function(standardFusioncalls,
                                 groupby,
                                 plotn = 20,
                                 countID,
                                 palette_rec = NULL,
                                 base_size = 20) {
  standardFusioncalls <- .check_annoFuse_calls(standardFusioncalls)
  stopifnot(is.character(groupby))
  stopifnot(is.numeric(plotn))
  stopifnot(is.character(countID))
  stopifnot(is.numeric(base_size))

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

  max_count <- rec_gene %>% 
    group_by(.data$Gene) %>% 
    summarise(max=sum(count)) %>% 
    pull(max) %>% 
    max()

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

  if (base_size != 20) {
    base_size <- base_size
  }

  rec_genes <- ggplot(rec_gene) +
    geom_col(aes(x = .data$Gene, y = .data$count, fill = !!as.name(groupby)), alpha = 0.75) +
    ylab(paste0("Number of ", countID)) +
    guides(color = FALSE, alpha = FALSE) +
    xlab(NULL) +
    scale_y_continuous(limits=c(0, max_count)) +
    rotate() +
    ggtitle("Genes Recurrently Fused") +
    theme_publication(base_size = base_size) +
    theme(
      legend.title = element_blank(), axis.text.y = element_text(face = "italic", angle = 0, hjust = 1),
      legend.position = "bottom"
    ) +
    scale_fill_manual(name = as.name(groupby), values = palette_2) +
    scale_x_discrete(limits = rev(levels(rec_gene$Gene)))

  return(rec_genes)
}
