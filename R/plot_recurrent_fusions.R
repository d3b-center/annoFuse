#' Function to plot recurrent fused genes

#' @param standardFusioncalls A dataframe from star fusion or arriba standardized to run through the filtering steps
#' @param groupby column name with grouping variables
#' @param plotn top n recurrent fusions to plot
#' @param countID column name to count recurrent fusions SampleID/ParticipantID/tumorID
#' @param palette_rec colors for grouping variables
#'
#' @export
#'
#' @return A ggplot object containing an overview on the recurrent fusions
#'
#' @examples
#' out_annofuse <- system.file("extdata", "PutativeDriverAnnoFuse.tsv", package = "annoFuse")
#' sfc <- read.delim(out_annofuse,stringsAsFactors = FALSE)
#' # keep only in-frame and fusions where both breakpoints are within genes
#' sfc <- as.data.frame(sfc[ which(Fusion_Type == "in-frame" && BreakpointLocation == "Genic"),])
#' plot_recurrent_fusions(sfc, groupby = "broad_histology", countID = "Kids_First_Participant_ID")
plot_recurrent_fusions <- function(standardFusioncalls,
                                   groupby,
                                   plotn = 20,
                                   countID,
                                   palette_rec = NULL) {
  
  standardFusioncalls <- .check_annoFuse_calls(standardFusioncalls)
  stopifnot(is.character(groupby))
  stopifnot(is.numeric(plotn))
  stopifnot(is.character(countID))
  
  stopifnot(all(c(groupby, countID) %in% colnames(standardFusioncalls)))
  
  # in-frame fusions only
  #  standardFusioncalls<-unique(standardFusioncalls) %>% dplyr::filter(.data$Fusion_Type=="in-frame")
  # remove geneA==geneB or intergenic
  #  standardFusioncalls<-standardFusioncalls[-which(standardFusioncalls$Gene1A==standardFusioncalls$Gene2A|standardFusioncalls$Gene1A==standardFusioncalls$Gene2B|standardFusioncalls$Gene2A==standardFusioncalls$Gene1B|standardFusioncalls$Gene1A==standardFusioncalls$Gene1B|standardFusioncalls$Gene2B==standardFusioncalls$Gene2A),]

  # gather recurrent per group
  rec_fusions <- standardFusioncalls %>%
    as.data.frame() %>%
    dplyr::select("FusionName", !!as.name(groupby), !!as.name(countID)) %>%
    unique() %>%
    group_by(.data$FusionName, !!as.name(groupby)) %>%
    dplyr::select(-!!as.name(countID)) %>%
    mutate(count = n()) %>%
    unique() %>%
    as.data.frame()

  # Top n recurrent fusion per group
  rec_fusions <- utils::head(rec_fusions[order(rec_fusions$count, decreasing = TRUE), ], plotn)

  if (!is.null(palette_rec)) {
    # provided palette is rownames=groupby values and color in column "color"
    palette_1 <- palette_rec
  } else {
    # palette to match recurrent fusion and recurrent fused genes
    n <- length(levels(as.factor(standardFusioncalls[, groupby])))
    palette <- grDevices::rainbow(n)
    names(palette) <- levels(as.factor(standardFusioncalls[, groupby]))
    colScale <- scale_colour_manual(name = as.name(groupby), values = palette)
    palette_1 <- palette[which(names(palette) %in% rec_fusions[, groupby])]
  }

  rec_fusions$FusionName <- factor(rec_fusions$FusionName, levels = unique(rec_fusions$FusionName), ordered = TRUE)
  rec_fusions_plot <- ggplot(rec_fusions) +
    geom_col(aes(x = FusionName, y = count, fill = !!as.name(groupby)), alpha = 0.75) +
    guides(alpha = FALSE) +
    ylab("Number of patients") +
    xlab(NULL) +
    guides(color = FALSE, alpha = FALSE) +
    scale_y_continuous(breaks = seq(0, 200, by = 20)) +
    ggpubr::rotate() +
    scale_fill_manual(name = as.name(groupby), values = palette_1) +
    theme_publication() +
    theme(legend.title = element_blank(), axis.text.y = element_text(face = "italic", angle = 0, hjust = 1, size = 12), legend.position = "bottom") +
    scale_x_discrete(limits = rev(levels(rec_fusions$FusionName)))

  return(rec_fusions_plot)
}
