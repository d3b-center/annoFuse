#' Standardizes fusion calls
#'
#' Various fusion callers have different formats that make aggregating and filtering data difficult.
#' By standardizing fusion callers output we capture the required columns which we use for downstream
#' analysis
#'
#' @param fusion_calls A dataframe from star fusion or arriba (more callers to be added)
#' @param caller string options STARFUSION/ARRIBA
#' @param tumorID string or character vector of same length as fusion_calls
#' @param input_json_file (optional) json format config file to provide input and output columns headers required for CUSTOM type and not required for other callers
#'
#' @developer Krutika S Gaonkar, Saksham Phul(phuls@chop.edu)
#' @export
#'
#' @return Standardized fusion calls ready for filtering
#'
#' @examples
#' # read in arriba fusion file
#' fusionfileArriba <- read_arriba_calls(
#'   system.file("extdata", "arriba_example.tsv", package = "annoFuseData")
#' )
#' # read in starfusion file
#' fusionfileStarFusion <- read_starfusion_calls(
#'   system.file("extdata", "starfusion_example.tsv", package = "annoFuseData")
#' )
#' formattedArriba <- fusion_standardization(fusionfileArriba,
#'   caller = "ARRIBA",
#'   tumorID = "tumorID"
#' )
#' formattedStarFusion <- fusion_standardization(fusionfileStarFusion,
#'   caller = "STARFUSION",
#'   tumorID = "tumorID"
#' )
#  read in CUSTOM type file
#' formattedCUSTOM <- fusion_standardization(fusionfileStarFusion,
#'   caller = "CUSTOM",
#'   tumorID = "All",
#'   input_json_file = "config"
#' )
#' format of the input_json_file ("Input_header" : "Output_header")
#  {
#  "CUSTOM":{
#	  	"Sample": "Sample_output",	
#		  "FusionName": "FusionName_output",
#		  "Gene1A": "Gene1A_output",
#		  "Gene1B": "Gene1B_output",
#		  "Gene2A": "Gene2A_output",
#		  "Gene2B": "Gene2B_output",
#		  "Fusion_Type":"Fusion_Type_output",
#		  "annots":"annots_output"
#	    } 
#  }

fusion_standardization <- function(fusion_calls,
                                   caller = c("STARFUSION", "ARRIBA","CUSTOM"),
                                   tumorID = "tumorID",
                                  input_json_file = "No file exists")
{
  stopifnot(is(fusion_calls, "data.frame"))
  stopifnot(is.character(caller))
  stopifnot(is.character(tumorID))

  # caller <- match.arg(caller, choices = c("STARFUSION", "ARRIBA"))
 #print(input_json_file)
  if (caller == "STARFUSION") {
    if (!any(grepl("PROT_FUSION_TYPE", colnames(fusion_calls)))) {
      # add NA if --examine_coding_effect was not used while running starfusion
      fusion_calls$PROT_FUSION_TYPE <- NA
    }
    fusion_calls <- fusion_calls %>%
      # standardize fusion type column name
      dplyr::rename(
        Fusion_Type = PROT_FUSION_TYPE,
        FusionName = "#FusionName"
      ) %>%
      dplyr::mutate(
        # remove chr notation from breakpoint columns
        LeftBreakpoint = gsub("^chr", "", .data$LeftBreakpoint),
        RightBreakpoint = gsub("^chr", "", .data$RightBreakpoint),
        # remove strand information to match breakpoint locations
        LeftBreakpoint = gsub(":[-|+]$", "", .data$LeftBreakpoint),
        RightBreakpoint = gsub(":[-|+]$", "", .data$RightBreakpoint),
        # STARFusion does not return confidence information
        Confidence = NA,
        # standardize fusion types
        Fusion_Type = dplyr::case_when(
          Fusion_Type == "INFRAME" ~ "in-frame",
          Fusion_Type == "FRAMESHIFT" ~ "frameshift",
          TRUE ~ "other"
        ),
        Sample = tumorID,
        Caller = "STARFUSION"
      )
       formatted_data=shape_output(fusion_calls)
   return(formatted_data)
  }
  else if (caller == "ARRIBA") {
    if (!any(colnames(fusion_calls) == "annots")) {
      fusion_calls$annots <- ""
    }
    fusion_calls <- fusion_calls %>%
      # standardizing fusion type annotation
      dplyr::rename(
        Fusion_Type = .data$reading_frame,
        Confidence = .data$confidence,
        # SpanningFragCount is equivalent to discordant_mates in Arriba
        SpanningFragCount = .data$discordant_mates
      ) %>%
      dplyr::mutate(
        LeftBreakpoint = gsub("^chr", "", .data$breakpoint1),
        RightBreakpoint = gsub("^chr", "", .data$breakpoint2),
        # readthrough information from arriba
        annots = paste(.data$annots, .data$type, sep = ","),
        # Intergenic gene fusion breakpoints in arriba are annotated as
        # "gene1A,gene2A". As comma is used as a common delimiter in files changing
        # it to "/"
        FusionName = paste0(gsub(",", "/", .data$`#gene1`), "--", gsub(",", "/", .data$gene2)),
        # JunctionReadCount is equivalent to split reads in Arriba. Arriba however
        # provides split_reads1 and split_reads2 to provide information of reads
        # anchoring in gene1 or gene2
        JunctionReadCount = .data$split_reads1 + .data$split_reads2,
        Fusion_Type = dplyr::case_when(
          !Fusion_Type %in% c("out-of-frame", "in-frame") ~ "other",
          Fusion_Type == "out-of-frame" ~ "frameshift",
          TRUE ~ "in-frame"
        ),
        Sample = tumorID,
        Caller = "ARRIBA"
      )
      formatted_data=shape_output(fusion_calls)
   return(formatted_data)
  }
  else if (caller == "CUSTOM") 
  {
        #Condition to check if required columns exists in the input and config file not provided by the user
        if (!file.exists(input_json_file) &&
            any(colnames(fusion_calls) == "Sample") &&
            any(colnames(fusion_calls) == "FusionName") &&
            any(colnames(fusion_calls) == "Gene1A") &&
            any(colnames(fusion_calls) == "Gene1B") &&
            any(colnames(fusion_calls) == "Gene2A") &&
            any(colnames(fusion_calls) == "Gene2B") &&
            any(colnames(fusion_calls) == "Fusion_Type") &&
            any(colnames(fusion_calls) == "annots"))
          {
            print("All required columns exists! ")
            return(fusion_calls)
          }
      else if(file.exists(input_json_file)) # if config file is provided
        { 
          print("Annotating based on config file") 
      
          json_data_frame <- as.data.frame(rjson::fromJSON(file = input_json_file)) #Get data from json file
          json_cols <- colnames(json_data_frame) # extract cols from json data frame
          caller <- strsplit(json_cols, split = "[.]")[[1]][1] #Get the caller from json file
  
          input_columns <- list() #define list to store input columns from json file
          output_columns <- list() #define list to store output columns from json file
          for(i in json_cols){
            input_columns <- append(input_columns,strsplit(i, split = "[.]")[[1]][2])  #extract input columns
          }
          
          if ("Sample" %in% input_columns &&
            "FusionName" %in% input_columns &&
            "Gene1A" %in% input_columns &&
            "Gene1B" %in% input_columns  &&
            "Gene2A" %in% input_columns &&
            "Gene2B" %in% input_columns &&
            "Fusion_Type" %in% input_columns &&
            "annots" %in% input_columns) #check if required column exists in config file else throw an exception
            {
            print("All required columns exists!")   
            for(i in 1:ncol(json_data_frame)) {       # for-loop over columns
              output_columns <- append(output_columns,json_data_frame[ , i])
            }
           # json_dataframe <- do.call(rbind, Map(data.frame, input_name=input_columns, output_name=output_columns))
          
            input_columns_without_None=input_columns 
            output_columns_without_None=output_columns
            output_columns_with_None <- list()
            index_none <-list() 

            for(i in 1:length(input_columns)){ #loop to extract None from input json
                if(input_columns[i] == "None"){
                  output_columns_with_None <-append(output_columns_with_None,output_columns[i])
                  index_none <- append(index_none,i)
                }
              }
            
            if(length(index_none > 0)){ #check if none exist as an input in the config file
                input_columns_without_None <- input_columns[- unlist(index_none)]
                output_columns_without_None <- output_columns[- unlist(index_none)]
            }
          
            standard_calls <- fusion_calls %>% dplyr::select(unlist(input_columns_without_None))
            standard_calls <- gdata::rename.vars(standard_calls, from = unlist(input_columns_without_None), to = unlist(output_columns_without_None))
            for(i in output_columns_with_None)
            {
                if( i != "Caller"){ #set caller for all the rows if json has None::Caller
                standard_calls[i] <- NA
                }
                else{
                    standard_calls[i] <- caller
                }
            }
            standard_calls <- standard_calls[, unlist(output_columns)]
            #print(standard_calls)
            return(standard_calls)
          }
          else{
            stop(paste("Provide all the required columns in the config file. Required columns for Custom caller are: Sample FusionName Gene1A Gene1B Gene2A Gene2B Fusion_Type annots."))
          }
        }  
       else #if user do not meet input expectations 
       {
        stop(paste(caller, " caller requires specific columns or config file do not exists. Required columns are: Sample FusionName Gene1A Gene1B Gene2A Gene2B Fusion_Type annots."))
       } 
    } 
  
 else 
 {
    stop(paste(caller, "is not a supported caller string."))
  }
}  

#function used by STARFUSION and ARRIBA to shape the final output columns as required
shape_output <- function(fusion_calls){
  # Get standard columns for filtering

  standard_calls <- fusion_calls %>%
    # select columns for standard fusion format
    dplyr::select(c(
      "LeftBreakpoint",
      "RightBreakpoint",
      "FusionName",
      "Sample",
      "Caller",
      "Fusion_Type",
      "JunctionReadCount",
      "SpanningFragCount",
      "Confidence",
      "annots"
    )) %>%
    # to obtain geneA and geneB for gene search below
    bind_cols(reshape2::colsplit(fusion_calls$FusionName, pattern = "--", names = c("GeneA", "GeneB"))) %>%
    # Intergenic fusion will have Gene1A,Gene2A,Gene1B,Gene2B
    separate(.data$GeneA, sep = "/", into = c("Gene1A", "Gene2A"), remove = FALSE) %>%
    separate(.data$GeneB, sep = "/", into = c("Gene1B", "Gene2B"), remove = FALSE) %>%
    # remove distance to fusion breakpoint from gene names in intergenic fusion
    mutate(
      Gene1A = gsub("[(].*", "", .data$Gene1A),
      Gene2A = gsub("[(].*", "", .data$Gene2A),
      Gene1B = gsub("[(].*", "", .data$Gene1B),
      Gene2B = gsub("[(].*", "", .data$Gene2B),
      BreakpointLocation = case_when(
        Gene1A == Gene1B & !grepl("/", FusionName) ~ "Intragenic",
        grepl("/", FusionName) ~ "Intergenic",
        TRUE ~ "Genic"
      ),
      SpanningDelta = SpanningFragCount - JunctionReadCount
    ) %>%
    as.data.frame()

  return(standard_calls)  
}
