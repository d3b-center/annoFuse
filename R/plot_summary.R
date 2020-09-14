#' Function to annotate fusion calls

#' @param standardFusioncalls A dataframe from star fusion or arriba standardized to run through the filtering steps
#' @param groupby column name with grouping variables
#' @param outputpdffile Filename to plot image - optional, if not specified, the 
#' plot is simply generated
#' @param base_size size for text in plots, defalt 20 for pdf outputs
#'
#' @export
#'
#' @return A ggplot object containing a summary of the standardized fusion calls
#'
#' @examples
#' out_annofuse <- system.file("extdata", "PutativeDriverAnnoFuse.tsv", package = "annoFuse")
#' sfc <- read.delim(out_annofuse, stringsAsFactors = FALSE)
#' plot_summary(sfc, groupby = "Fusion_Type")
plot_summary <- function(standardFusioncalls,
                         groupby,
                         outputpdffile,
                         base_size=20) {
  standardFusioncalls <- .check_annoFuse_calls(standardFusioncalls)
  if (!missing(outputpdffile))
    stopifnot(is.character(outputpdffile))
  
  if (missing(groupby)) {
    # per sample
    groupby <- "Sample"
  } else {
    stopifnot(is.character(groupby))
  }

  # fusion_calls_interchromosomal<-standardFusioncalls %>% dplyr::filter(grepl("INTERCHROMOSOMAL",annots))
  # fusion_calls_intrachromosomal<-standardFusioncalls %>% dplyr::filter(grepl("INTRACHROMOSOMAL",annots))
  #
  # fusion_chrom<-rbind(data.frame(cbind(fusion_calls_interchromosomal,"Distance"=rep("Interchromosomal",nrow(fusion_calls_interchromosomal)))),cbind(fusion_calls_intrachromosomal,"Distance"=rep("Intrachromosomal",nrow(fusion_calls_intrachromosomal)))) %>% unique()

  # Inter and intrachromosomal plot
  fusion_chrom <- standardFusioncalls %>%
    # get left breakpoints
    group_by(.data$LeftBreakpoint) %>%
    # get left chromosome
    mutate(Leftchr = strsplit(.data$LeftBreakpoint, ":")[[1]][1]) %>%
    # get right breakpoint
    unnest(cols = c()) %>%
    group_by(.data$RightBreakpoint) %>%
    # get right chr
    mutate(Rightchr = strsplit(.data$RightBreakpoint, ":")[[1]][1]) %>%
    unnest(cols = c()) %>%
    mutate(Distance = ifelse(.data$Leftchr == .data$Rightchr, "INTRACHROMOSOMAL", "INTERCHROMOSOMAL")) %>%
    dplyr::select(.data$Distance, .data$LeftBreakpoint, .data$RightBreakpoint, .data$Sample, !!as.name(groupby)) %>%
    unique()

  p1 <- ggplot(fusion_chrom, aes(x = !!as.name(groupby), fill = fusion_chrom$Distance, alpha = 0.75)) +
    geom_bar() +
    theme_publication(base_size = base_size) +
    theme(legend.position = "top") +
    xlab(as.name(groupby)) +
    ylab("Count") +
    guides(alpha = FALSE) +
    theme(axis.text.x = element_text(angle = 45)) +
    labs(fill = "Distance")

  # frame information
  p2 <- ggplot(standardFusioncalls, aes(fill = .data$Fusion_Type, x = .data$Caller, alpha = 0.75)) +
    geom_bar(aes(y = log2(stat(count)))) +
    rotate() +
    xlab("Fusion Caller") +
    ylab("Count (log2)") +
    theme_publication(base_size = base_size) +
    theme(legend.position = "top") +
    guides(alpha = FALSE) +
    theme(axis.text.x = element_text(angle = 0, hjust = 1), axis.text.y = element_text(angle = 0, vjust = 2)) +
    labs(fill = "Frame")


  # keep fusion if atleast 1 is protein-coding
  # BiocManager::install("EnsDb.Hsapiens.v86")
  edb <- EnsDb.Hsapiens.v86

  genes <- as.data.frame(genes(edb))

  # keep if atleast 1 gene is protein coding gene in fusions
  fusion_protein_coding_gene_df <- standardFusioncalls %>%
    # We want to keep track of the gene symbols for each sample-fusion pair
    dplyr::select(.data$Sample, .data$FusionName, .data$Gene1A, .data$Gene1B, .data$Gene2A, .data$Gene2B, !!(as.name(groupby))) %>%
    # We want a single column that contains the gene symbols
    tidyr::gather(Gene1A, Gene1B, Gene2A, Gene2B,
      key = gene_position, value = GeneSymbol
    ) %>%
    # Remove columns without gene symbols
    dplyr::filter(.data$GeneSymbol != "") %>%
    dplyr::arrange(.data$Sample, .data$FusionName) %>%
    # Retain only distinct rows
    dplyr::distinct() %>%
    left_join(genes, by = c("GeneSymbol" = "gene_name"))


  # bubble plot for genes in biotype fromensemble
  fusion_type_gene_df <- fusion_protein_coding_gene_df %>%
    # We want to keep track of the gene symbols for each sample-fusion pair
    dplyr::select(.data$Sample, .data$FusionName, .data$gene_position, .data$gene_biotype, !!(as.name(groupby))) %>%
    unique() %>%
    group_by(!!(as.name(groupby)), .data$gene_biotype, .data$gene_position) %>%
    dplyr::summarise(Type.ct = n()) %>%
    dplyr::filter(gene_position %in% c("Gene1A", "Gene1B"))


  fusion_type_gene_df$gene_position[grep("Gene1A", fusion_type_gene_df$gene_position)] <- "5'-Gene"
  fusion_type_gene_df$gene_position[grep("Gene1B", fusion_type_gene_df$gene_position)] <- "3'-Gene"
  fusion_type_gene_df$gene_position <- factor(fusion_type_gene_df$gene_position, levels = c("5'-Gene", "3'-Gene"), ordered = TRUE)


  p3 <- ggplot(fusion_type_gene_df, aes(x = gene_biotype, y = gene_position, size = Type.ct, fill = gene_position, color = gene_position, alpha = 0.75, width = 2.5, height = 5)) +
    geom_point(shape = 21) +
    guides(size = FALSE, color = FALSE, alpha = FALSE) +
    scale_size_continuous(range = c(1, 12)) +
    xlab("Gene Biotype") +
    ylab("Gene Position") +
    theme_publication(base_size = base_size) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(angle = 0, vjust = 1), legend.position = "left") +
    guides(fill = FALSE)


  # Kinase groups
  bioMartDataPfam <- readRDS(system.file("extdata", "pfamDataBioMart.RDS", package = "annoFuse"))
  annDomain <- get_Pfam_domain(standardFusioncalls = standardFusioncalls, bioMartDataPfam = bioMartDataPfam, keepPartialAnno = TRUE)
  gene1AKinaseDomain <- annDomain$Gene1A %>%
    dplyr::filter(Gene1A_DOMAIN_RETAINED_IN_FUSION == "Yes" & grepl("kinase", NAME)) %>%
    dplyr::select("Gene1A", "NAME", !!(as.name(groupby)), Sample) %>%
    dplyr::rename("GeneSymbol" = "Gene1A", "Domain" = "NAME") %>%
    dplyr::mutate("gene_position" = "Gene1A")
  gene1BKinaseDomain <- annDomain$Gene1B %>%
    dplyr::filter(Gene1B_DOMAIN_RETAINED_IN_FUSION == "Yes" & grepl("kinase", NAME)) %>%
    dplyr::select("Gene1B", "NAME", !!(as.name(groupby)), Sample) %>%
    dplyr::rename("GeneSymbol" = "Gene1B", "Domain" = "NAME") %>%
    dplyr::mutate("gene_position" = "Gene1B")

  fusion_kinase_protein_df <- rbind(gene1AKinaseDomain, gene1BKinaseDomain)

  if (!is_empty(fusion_kinase_protein_df$GeneSymbol)) {
    kinase_fusion <- fusion_kinase_protein_df %>%
      dplyr::select("gene_position", "Domain", !!(as.name(groupby)), GeneSymbol, Sample) %>%
      unique() %>%
      group_by(.data$Domain, .data$gene_position) %>%
      dplyr::summarise(Type.ct = n()) %>%
      unique() %>%
      dplyr::filter(gene_position %in% c("Gene1A", "Gene1B"))

    kinase_fusion$gene_position[grep("Gene1A", kinase_fusion$gene_position)] <- "5'-Gene"
    kinase_fusion$gene_position[grep("Gene1B", kinase_fusion$gene_position)] <- "3'-Gene"
    kinase_fusion$gene_position <- factor(kinase_fusion$gene_position, levels = c("5'-Gene", "3'-Gene"), ordered = TRUE)

    p4 <- ggplot(kinase_fusion, aes(x = Domain, y = Type.ct, fill = gene_position, alpha = 0.95)) +
      geom_col() +
      theme_publication(base_size = base_size) +
      ylab("Count") +
      xlab("Kinase domain retained") +
      guides(alpha = FALSE) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(angle = 0, vjust = 2)) +
      labs(fill = "Gene Position")
  } else {
    p4 <- ggplot() +
      theme_void()
  }

  # bubble plot for genes in kinase/oncoplot etc
  fusion_type_bubblegene_df <- standardFusioncalls %>%
    # We want to keep track of the gene symbols for each sample-fusion pair
    dplyr::select(.data$Sample, .data$FusionName, .data$Gene1A_anno, .data$Gene1B_anno, !!(as.name(groupby))) %>%
    unique() %>%
    # We want a single column that contains the gene symbols
    tidyr::gather(Gene1A_anno, Gene1B_anno,
      key = gene_position, value = Annot
    ) %>%
    mutate(gene_position = gsub("_anno", "", .data$gene_position)) %>%
    # Remove columns without gene symbols
    dplyr::filter(.data$Annot != "") %>%
    mutate(Annotation = gsub(" ", "", .data$Annot)) %>%
    mutate(Annotation = strsplit(as.character(.data$Annotation), ",")) %>%
    unnest(.data$Annotation) %>%
    group_by(!!as.name(groupby), .data$Annotation, .data$gene_position) %>%
    dplyr::summarise(Annot.ct = n()) %>%
    # To only plot the broader transcription factor annotation
    dplyr::filter(!.data$Annotation %in% c("curatedTF", "predictedTF"))

  # To-Do rename in original table
  fusion_type_bubblegene_df$Annotation[grep("onco", fusion_type_bubblegene_df$Annotation)] <- "Oncogene"
  fusion_type_bubblegene_df$Annotation[grep("tsgs", fusion_type_bubblegene_df$Annotation)] <- "TumorSuppressorGene"
  fusion_type_bubblegene_df$Annotation[grep("kinase", fusion_type_bubblegene_df$Annotation)] <- "Kinase"

  fusion_type_bubblegene_df$gene_position[grep("Gene1A", fusion_type_bubblegene_df$gene_position)] <- "5'-Gene"
  fusion_type_bubblegene_df$gene_position[grep("Gene1B", fusion_type_bubblegene_df$gene_position)] <- "3'-Gene"
  fusion_type_bubblegene_df$gene_position <- factor(fusion_type_bubblegene_df$gene_position, levels = c("5'-Gene", "3'-Gene"), ordered = TRUE)

  p5 <- ggplot(fusion_type_bubblegene_df, aes(x = !!as.name(groupby), y = Annot.ct, size = Annot.ct, fill = Annotation, color = Annotation, alpha = 0.75)) +
    geom_point(shape = 21) +
    guides(size = FALSE, alpha = FALSE) +
    scale_size_continuous(range = c(1, 12)) +
    ylab("Count") +
    xlab(as.name(groupby)) +
    theme_publication(base_size = base_size) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(angle = 0, vjust = 2), legend.position = "left") +
    facet_wrap(~gene_position)


  if (!missing(outputpdffile)) {
    ggarrange(p1, p2, p3, p4, p5, 
              labels = c("A", "B", "C", "D", "E"), 
              heights = c(5, 4, 6),
              widths = c(2, 1),
              nrow = 3, ncol = 2,
              font.label = list(size = base_size )) %>% ggexport(filename = outputpdffile, width = 20, height = 20)
  } else {
    p1 <- ggplot(fusion_chrom, aes(x = !!as.name(groupby), fill = fusion_chrom$Distance, alpha = 0.75)) +
      geom_bar() +
      theme_publication(base_size = base_size) +
      theme(legend.position = "top") +
      xlab(as.name(groupby)) +
      ylab("Count") +
      guides(alpha = FALSE) +
      theme(axis.text.x = element_text(angle = 45), plot.margin = unit(c(1, 1, 0, 10), "mm")) +
      labs(fill = "Distance")

    # frame information
    p2 <- ggplot(standardFusioncalls, aes(fill = .data$Fusion_Type, x = .data$Caller, alpha = 0.75)) +
      geom_bar(aes(y = log2(stat(count)))) +
      rotate() +
      xlab("Fusion Caller") +
      ylab("Count (log2)") +
      theme_publication(base_size = base_size) +
      theme(legend.position = "top") +
      guides(alpha = FALSE) +
      theme(axis.text.x = element_text(angle = 0, hjust = 1), axis.text.y = element_text(angle = 0, vjust = 2), plot.margin = unit(c(1, 10, 10, 1), "mm")) +
      labs(fill = "Frame")

    p3 <- ggplot(fusion_type_gene_df, aes(x = gene_biotype, y = gene_position, size = Type.ct, fill = gene_position, color = gene_position, alpha = 0.75, width = 2.5, height = 5)) +
      geom_point(shape = 21) +
      guides(size = FALSE, color = FALSE, alpha = FALSE) +
      scale_size_continuous(range = c(0.1, 3)) +
      xlab("Gene Biotype") +
      ylab("Gene Position") +
      theme_publication(base_size = base_size) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, ), axis.text.y = element_text(angle = 0, vjust = 1), legend.position = "left", plot.margin = unit(c(1, 1, 1, 10), "mm")) +
      guides(fill = FALSE)

    if (!is_empty(fusion_kinase_protein_df$GeneSymbol)) {
      p4 <- ggplot(kinase_fusion, aes(x = Domain, y = Type.ct, fill = gene_position, alpha = 0.95)) +
        geom_col() +
        theme_publication(base_size = base_size) +
        ylab("Count") +
        xlab("Kinase domain retained") +
        guides(alpha = FALSE) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(angle = 0, vjust = 2), plot.margin = unit(c(1, 1, 1, 1), "mm")) +
        labs(fill = "Gene Position")
    } else {
      p4 <- ggplot() +
        theme_void()
    }

    p5 <- ggplot(fusion_type_bubblegene_df, aes(x = !!as.name(groupby), y = Annot.ct, size = Annot.ct, fill = Annotation, color = Annotation, alpha = 0.75, width = 2, height = 3)) +
      geom_point(shape = 21) +
      guides(size = FALSE, alpha = FALSE) +
      scale_size_continuous(range = c(0.1, 3)) +
      ylab("Count") +
      xlab(as.name(groupby)) +
      theme_publication(base_size = base_size) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(angle = 0, vjust = 2), legend.position = "left") +
      theme(plot.margin = unit(c(1, 1, 1, 10), "mm")) +
      facet_wrap(~gene_position)

    p <- ggarrange(p1, p2, p3, p4, p5, 
                   labels = c("A", "B", "C", "D", "E"), 
                   heights = c(5, 3, 5), 
                   widths = c(2, 1), 
                   nrow = 3, ncol = 2, 
                   font.label = list(size = base_size))

    return(p)
  }
}
