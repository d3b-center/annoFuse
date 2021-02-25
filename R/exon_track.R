#' Function gather exon locations per gene to plot and annotate breakpoints

#' @param genome human genome version (hg38 and hg19 supported)
#' @param gtf gencode gtf file path,  gencode.v27.primary_assembly.annotation.gtf.gz can be downloaded from "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.primary_assembly.annotation.gtf.gz"
#' 
#' @export
#'
#' @return exon location dataframe
#'
#' @examples
#' library(EnsDb.Hsapiens.v86)
#' out_annofuse <- system.file("extdata", "PutativeDriverAnnoFuse.tsv", package = "annoFuse")
#' sfc <- read.delim(out_annofuse)
#' download.file("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.primary_assembly.annotation.gtf.gz","gencode.v27.primary_assembly.annotation.gtf.gz")
#' # using hg38 gencode gtf to get exon tracks
#' exons <- get_exons_tracks_gtf(gtf_path = file.path("gencode.v27.primary_assembly.annotation.gtf.gz"))
#' 
                                  
get_exons_tracks_gtf <- function(genome="hg38",
                                 gtf_path
                                 ){
  if(genome=="hg38"){
    # to match geneID hg38
    edb <- EnsDb.Hsapiens.v86
  }
  if (genome == "hg19"){
    # to match geneID hg19
    edb <- EnsDb.Hsapiens.v75
  }
ids_to_match <- mapIds(edb, keys=ensembldb::keys(edb, "GENEID"),
                       column = c("GENENAME","GENEID"),
                               keytype = "GENEID") %>%
  as.data.frame() %>%
  rownames_to_column( var = "geneID") %>%
  rename("geneName"=".")


# transcript db 
txdb <- GenomicFeatures::makeTxDbFromGFF(
  file =gtf_path,
  format = "gtf"
)

# exon locations
exons_to_plot <- ensembldb::select(txdb, 
                           ensembldb::keys(txdb, "GENEID"),
                           keytype="GENEID",
                           columns=c("GENEID",
                                     "TXNAME",
                                     "EXONCHROM",
                                     "EXONRANK",
                                     "EXONSTART",
                                     "EXONEND",
                                     "EXONSTRAND")) %>%
  dplyr::rename(geneID=GENEID,
                transcript =TXNAME,
                contig = EXONCHROM,
                start = EXONSTART,
                end = EXONEND,
                strand = EXONSTRAND,
                exonNumber = EXONRANK) %>%
  dplyr::mutate(type="exon")

# gene locations
genes_to_plot <- genes(txdb) 
genes_to_plot <- data.frame(
  contig = seqnames(genes_to_plot),
  type = "gene",
  start = start(genes_to_plot),
  end = end(genes_to_plot),
  strand = strand(genes_to_plot),
  geneID = genes_to_plot$gene_id,
  transcript = "",
  exonNumber = 0
) 

# combined   
tracks_to_plot <- dplyr::bind_rows(exons_to_plot,
                                   genes_to_plot) %>%   
  # to be match with ids_to_match
  dplyr::mutate(geneID = gsub("\\..*", "", geneID)) %>%
  dplyr::left_join(unique(ids_to_match),
                   by=c("geneID")) 

}

