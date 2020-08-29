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
#' @return A ggplot object containing an overview on the recurrent fusions
#'
#' @examples
#' out_annofuse <- system.file("extdata", "PutativeDriverAnnoFuse_test_v14.tsv",
#'   package = "annoFuse"
#' )
#' sfc <- read.delim(out_annofuse, stringsAsFactors = FALSE)
#' clinical <- read.delim(
#'   system.file("extdata", "pbta-histologies.tsv", package = "annoFuse")
#' )
#' # Select only in-frame and frameshift
#' library("dplyr")
#' sfc <- sfc %>%
#'   dplyr::filter(Fusion_Type %in% c("in-frame", "frameshift"))
#' sfc <- merge(sfc, clinical[, c("Kids_First_Biospecimen_ID", "broad_histology")],
#'   by.x = "Sample", by.y = "Kids_First_Biospecimen_ID"
#' )
#' # Remove Benign tumor fusions
#' sfc <- sfc[-which(sfc$FusionName %in%
#'   unique(sfc[which(sfc$broad_histology == "Benign tumor"), "FusionName"])), ]
#' # Remove GeneA == GeneB
#' sfc <- sfc[-which(sfc$Gene1A == sfc$Gene1B |
#'   sfc$Gene1A == sfc$Gene2B |
#'   sfc$Gene1B == sfc$Gene2A), ]
#' # Remove intergenic fusions
#' sfc <- sfc[-grep("/", sfc$FusionName), ]
#' plot_recurrent_fusions(sfc,
#'   groupby = "broad_histology",
#'   countID = "Kids_First_Participant_ID"
#' )
plot_recurrent_fusions <- function(standardFusioncalls,
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

  max_count <- rec_fusions %>% 
    group_by(FusionName) %>% 
    summarise(max=sum(count)) %>% 
    pull(max) %>% 
    max()

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

  if (base_size != 20) {
    base_size <- base_size
  }

  rec_fusions$FusionName <- factor(rec_fusions$FusionName, levels = unique(rec_fusions$FusionName), ordered = TRUE)
  rec_fusions_plot <- ggplot(rec_fusions) +
    geom_col(aes(x = FusionName, y = count, fill = !!as.name(groupby)), alpha = 0.75) +
    guides(alpha = FALSE) +
    ylab(paste0("Number of ", countID)) +
    xlab(NULL) +
    guides(color = FALSE, alpha = FALSE) +
    scale_y_continuous(limits=c(0, max_count)) +
    ggpubr::rotate() +
    scale_fill_manual(name = as.name(groupby), values = palette_1) +
    ggtitle("Recurrent Fusions") +
    theme_publication(base_size = base_size) +
    theme(legend.title = element_blank(), axis.text.y = element_text(face = "italic", angle = 0, hjust = 1), legend.position = "bottom") +
    scale_x_discrete(limits = rev(levels(rec_fusions$FusionName)))

  return(rec_fusions_plot)
}
