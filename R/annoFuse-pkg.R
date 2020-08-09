#' annoFuse
#'
#' annoFuse is an R package for ...
#'
#' @importFrom dplyr arrange bind_cols case_when count desc distinct filter
#' group_by inner_join left_join mutate n rename select summarise ungroup %>%
#' @importFrom ggplot2 aes coord_flip element_blank element_line element_rect
#' element_text facet_wrap geom_bar geom_col geom_hline geom_line
#' geom_linerange geom_point geom_segment ggplot ggtitle guides labs
#' rel scale_colour_manual scale_fill_manual scale_size_continuous
#' scale_x_discrete scale_y_continuous stat theme theme_void xlab ylab .data
#' @importFrom ggpubr ggarrange ggexport rotate
#' @importFrom ggthemes theme_foundation
#' @importFrom grid arrow unit
#' @importFrom qdapRegex rm_between
#' @importFrom readr col_character cols read_tsv
#' @importFrom reshape2 colsplit dcast melt
#' @importFrom stats na.omit
#' @importFrom stringr str_detect str_replace
#' @importFrom tibble add_column column_to_rownames remove_rownames
#' @importFrom tidyr gather separate unnest one_of
#' @importFrom utils head read.delim browseURL
#' @importFrom EnsDb.Hsapiens.v86 EnsDb.Hsapiens.v86
#' @importFrom purrr is_empty
#' @importFrom rmarkdown render
#'
#' @name annoFuse-pkg
#' @docType package
NULL

globalVariables(c(
  "Gene1A", "Gene1B", "Gene2A", "Gene2B", "GeneSymbol", "Gene_Symbol",
  "gene_position", "type", "FusionName", "gene_id", "FPKM",
  "gene_biotype", "Type.ct", "NAME", "Sample", "Domain",
  "Gene1A_DOMAIN_RETAINED_IN_FUSION", "Gene1B_DOMAIN_RETAINED_IN_FUSION",
  "Gene1A_anno", "Gene1B_anno", "Annot", "Annot.ct", "Annotation"
))
