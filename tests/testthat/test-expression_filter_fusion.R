context("Test expression filtering")

fusionfileArriba <- read_arriba_calls(system.file("extdata", "arriba_example.tsv", package = "annoFuseData"))

formattedArriba <- 
  annoFuse::fusion_standardization(fusion_calls = fusionfileArriba,
                                   caller= "ARRIBA",
                                   tumorID = "BS_W97QQYKQ")


fusionfileStarfusion <- read_arriba_calls(system.file("extdata", "starfusion_example.tsv", package = "annoFuseData"))

formattedStarFusion <- 
  annoFuse::fusion_standardization(fusion_calls = fusionfileStarfusion,
                                   caller= "STARFUSION",
                                   tumorID = "BS_W97QQYKQ")

standardFusioncalls <- as.data.frame(rbind(formattedStarFusion, formattedArriba))

expressionFile <- system.file("extdata", "example.rsem.genes.results.gz", package = "annoFuseData")
expressionMatrix <- read_tsv(expressionFile)
library(reshape2)
# split gene id and symbol
expressionMatrix <- cbind(expressionMatrix, 
   colsplit(expressionMatrix$gene_id, pattern = "_", names = c("EnsembleID", "GeneSymbol")))
# collapse to matrix of HUGO symbols x Sample identifiers
# take max expression per row and use the max value for duplicated gene symbols
expressionMatrix.collapsed <- expressionMatrix %>%
   arrange(desc(FPKM)) %>% # arrange decreasing by FPKM
   distinct(GeneSymbol, .keep_all = TRUE) %>% # keep the ones with greatest FPKM value. 
                                              # If ties occur, keep the first occurencce
   unique() %>%
   remove_rownames() %>%
   dplyr::select(.data$EnsembleID, .data$GeneSymbol, .data$FPKM, .data$gene_id)
 # rename columns

colnames(expressionMatrix.collapsed)[3] <- "BS_W97QQYKQ"
expressionFiltered <- expression_filter_fusion(
   standardFusioncalls = standardFusioncalls, 
   expressionMatrix = expressionMatrix.collapsed,
   expressionFilter = 1)

test_that("Called by N caller output for standard fusion file input", {
  expect_equal(colnames(expressionFiltered), c("LeftBreakpoint","RightBreakpoint","FusionName","Sample","Caller","Fusion_Type","JunctionReadCount","SpanningFragCount","Confidence","annots","Gene1A","Gene2A","Gene1B","Gene2B","BreakpointLocation","SpanningDelta" ))
  expect_equal(nrow(standardFusioncalls), 26)
  expect_equal(unique(standardFusioncalls$Sample), "BS_W97QQYKQ")
})
