#' Given pfam id(s) and an ensembl mart
#' creates a grange for the pfam domains
#' 
#' @param ensembl a S4 class from biomaRt::useMart() call
#' @param pfamDesc_path path to pfam description file
#' @param ucscGenePfam_path path to ucsc gene pfam genomic location
#' @param pfam_id pfam id(s) to filter
#' @param return_pfam_gr logical value to return a pfam location granges
#' 
#'
#' @return granges for genomic locations of given pfam id(s)
#'
#' @examples
#' \dontrun{
#' # hsapiens_gene_ensembl == Human genes (GRCh38.p13) 10th Feb 2021
#' ensembl <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "ensembl.org")
#'
#' download.file("http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/pfamDesc.txt.gz", "pfamDesc.txt.gz")
#' download.file("http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ucscGenePfam.txt.gz", "ucscGenePfam.txt.gz")
#'
#' bioMartDataPfam <- get_biomart_pfam_merge(ensembl = ensembl, pfamDesc_path = "pfamDesc.txt.gz", ucscGenePfam_path = "ucscGenePfam.txt.gz")
#' }
#'
get_biomart_pfam_merge <- function(ensembl, pfamDesc_path, ucscGenePfam_path, pfam_id, return_pfam_gr = NULL) {
  if (!missing(pfam_id)) {
    dataBioMart <- getBM(
      attributes = c(
        "hgnc_symbol", "pfam",
        "chromosome_name", "start_position", "end_position",
        "strand"
      ),
      filters = "pfam",
      values = pfam_id,
      mart = ensembl
    ) %>%
      dplyr::mutate(strand = dplyr::if_else(strand == 1, "+", "-"))
  } else {
    dataBioMart <- getBM(
      attributes = c(
        "hgnc_symbol",
        "chromosome_name", "start_position", "end_position",
        "strand"
      ),
      filters = "chromosome_name",
      values = c(
        "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
        "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
        "21", "22", "X", "Y", "MT"
      ),
      mart = ensembl
    ) %>%
      dplyr::mutate(strand = dplyr::if_else(strand == 1, "+", "-"))
  }

  # Pfam file downloaded from
  # contains pfam IDs and associated names
  pfam <- read_tsv(pfamDesc_path, col_names = FALSE) %>% as.data.frame()
  colnames(pfam) <- c("ID", "NAME", "DESC")

  # UCSC file downloaded from
  # contains pfam names and chromosome positions
  pfam_location <- read.delim(ucscGenePfam_path, header = F)
  colnames(pfam_location) <- c(
    "BIN", "CHROM", "CHROM_START", "CHROM_END", "NAME", "SCORE",
    "STRAND", "THICK_START", "THICK_END", "RESERVED", "BLOCK_COUNT",
    "BLOCK_SIZES", "CHROM_STARTS"
  )

  # merge using pfam names tp get pfam ID
  pfamDescLoc <- pfam %>%
    left_join(pfam_location, by = "NAME") %>%
    dplyr::select("ID", "NAME", "DESC", "CHROM", "CHROM_START", "CHROM_END", "STRAND") %>%
    dplyr::rename(
      "domain_chr" = "CHROM",
      "domain_start" = "CHROM_START",
      "domain_end" = "CHROM_END",
      "domain_strand" = "STRAND"
    ) %>%
    dplyr::filter(!is.na(domain_start) & !is.na(domain_end)) %>%
    dplyr::mutate(domain_chr = gsub("chr", "", domain_chr))


  # granges
  pfam_gr <- makeGRangesFromDataFrame(df = pfamDescLoc, keep.extra.columns = TRUE, seqnames.field = "domain_chr", start.field = "domain_start", end.field = "domain_end", strand.field = "domain_strand")

  gene_gr <- makeGRangesFromDataFrame(df = dataBioMart, keep.extra.columns = TRUE, seqnames.field = "chromosome_name", start.field = "start_position", end.field = "end_position", strand.field = "strand")

  if (!is.null(return_pfam_gr)) {
    # granges
    return(pfam_gr)
  } else {
    # get overlap gene and pfam location
    # dataframe
    pfamGene_gr <- mergeByOverlaps(pfam_gr, gene_gr, type = "within")
    pfamGene <- data.frame(
      "hgnc_symbol" = pfamGene_gr$hgnc_symbol,
      "pfam_id" = pfamGene_gr$ID,
      "chromosome_name" = unlist(lapply(seqnames(pfamGene_gr$gene_gr), as.character)),
      "gene_start" = start(pfamGene_gr$gene_gr),
      "gene_end" = end(pfamGene_gr$gene_gr),
      "strand" = strand(pfamGene_gr$gene_gr),
      "NAME" = pfamGene_gr$NAME,
      "DESC" = pfamGene_gr$DESC,
      "domain_chr" = unlist(lapply(seqnames(pfamGene_gr$pfam_gr), as.character)),
      "domain_start" = start(pfamGene_gr$pfam_gr),
      "domain_end" = end(pfamGene_gr$pfam_gr),
      stringsAsFactors = FALSE
    )


    return(pfamGene)
  }
}
