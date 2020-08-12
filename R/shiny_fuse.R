

#' shiny_fuse
#'
#' TODO
#'
#' @param out_annofuse TODO
#'
#' @return TODO
#' @export
#'
#' @rawNamespace import(shiny, except = c(renderDataTable, dataTableOutput))
#' @import shinydashboard
#' @import rintrojs
#' @import shinythemes
#' @importFrom DT datatable renderDataTable dataTableOutput
#'
#' @examples
#' out_annofuse <- "/Users/fede/Development/annoFuse/PutativeDriverAnnoFuse_test_v14.tsv"
#' if (interactive()) {
#'   shiny_fuse(out_annofuse)
#' }
shiny_fuse <- function(out_annofuse = NULL) {

  # checks on the objects provided
  # ?

  # UI definition -----------------------------------------------------------
  shinyfuse_ui <- shinydashboard::dashboardPage(
    skin = "black",

    # header definition -------------------------------------------------------
    header = shinydashboard::dashboardHeader(
      title = "ShinyFuse",
      titleWidth = 350,
      shinydashboard::dropdownMenu(
        type = "tasks",
        icon = icon("question-circle fa-1g"),
        badgeStatus = NULL,
        headerText = "Documentation",
        shinydashboard::notificationItem(
          text = actionButton(
            "interface_overview", "Overview of the interface",
            icon("hand-o-right")
          ),
          icon = icon(""), # tricking it to not have additional icon
          status = "primary"
        )
      )
    ),

    # sidebar definition ------------------------------------------------------
    sidebar = shinydashboard::dashboardSidebar(
      width = 250,
      collapsed = !is.null(out_annofuse),
      
      uiOutput("choose_annofusedata_file")
      # shinydashboard::menuItem(
      #   text = "Input Settings", icon = icon("cog"),
      #   startExpanded = TRUE,
      #   numericInput(
      #     inputId = "whatevs",
      #     label = "number of genesets",
      #     value = 15, min = 1, max = 50
      #   )
      # )
    ),

    # body definition ---------------------------------------------------------
    body = shinydashboard::dashboardBody(
      rintrojs::introjsUI(),
      # useShinyjs(),
      id = "main-app",
      navbarPage(
        title = div(),
        windowTitle = "shinyfuse",
        footer = "",
        theme = shinytheme("cosmo"),
        selected = "TableExplorer",
        tabPanel(
          title = "TableExplorer", icon = icon("table"),
          fluidPage(
            fluidRow(
              column(
                width = 8,
                DT::dataTableOutput("table_annofuse")
              ),
              column(
                width = 4,
                h4("Some content, for example linked to the selected row"),
                uiOutput("geneinfo_ui"),
                h4("Expanding more on this..."),
                uiOutput("geneplots_ui")
              )
            )
          )
        ),
        tabPanel(
          title = "TableSummary", icon = icon("dna"),
          fluidPage(
            fluidRow(
              column(
                width = 12,
                plotOutput("af_overview")
              )
            ),
            fluidRow(
              column(
                width = 6,
                plotOutput("af_recurrentfusions")
              ),
              column(
                width = 6,
                plotOutput("af_recurrentgenes")
              )
            )
          )
        ),
        navbarMenu(
          title = "About", icon = icon("question-circle"),
          tabPanel(
            title = "apanel1", icon = icon("laptop-code"),
            fluidPage(
              h1("About - panel1")
            )
          ),
          tabPanel(
            title = "apanel2", icon = icon("university"),
            fluidPage(
              h1("About - panel2")
            )
          ),
          tabPanel(
            title = "news", icon = icon("newspaper"),
            fluidPage(
              # includeMarkdown("NEWS.md")
            )
          )
        )
      )
    )
  )


  # Server definition -------------------------------------------------------
  shinyfuse_server <- function(input, output, session) {
    
    # Initializing data storage ------------------------------------------------
    values <- reactiveValues()
    values$annofuse_tbl <- NULL
    values$enhanced_annofuse_tbl <- NULL
    
    # Define data file if annoFuse data is not provided ------------------------
    if (is.null(out_annofuse)) {
      output$choose_annofusedata_file <- renderUI({
        menuItem(
          icon = icon("file-alt"),
          fileInput(inputId = "annofusedatasel",
                    label = "Load sample data file (tab-separated)",
                    accept = c("text/tab-separated-values", "text/plain",
                               ".tsv", ".tab", ".txt"),
                    multiple = FALSE)
        )
      })
    } else {
      annofuse_tbl <- read.delim(out_annofuse)
      values$annofuse_tbl <- annofuse_tbl
      
      enhanced_annofuse_tbl <- annofuse_tbl
      # enhancing the content of the table
      enhanced_annofuse_tbl$Gene1A <- .multilink(enhanced_annofuse_tbl$Gene1A)
      enhanced_annofuse_tbl$Gene1B <- .multilink(enhanced_annofuse_tbl$Gene1B)
      values$enhanced_annofuse_tbl <- enhanced_annofuse_tbl
    }
    
    # Load annoFuse data file --------------------------------------------------
    observeEvent(input$annofusedatasel, {
      values$annofuse_tbl <- read.delim(input$annofusedatasel$datapath)
      values$enhanced_annofuse_tbl <- values$annofuse_tbl
     
      # enhancing the content of the table
      values$enhanced_annofuse_tbl$Gene1A <- .multilink(values$enhanced_annofuse_tbl$Gene1A)
      values$enhanced_annofuse_tbl$Gene1B <- .multilink(values$enhanced_annofuse_tbl$Gene1B)
    })
    
    
    
    
    # TODO? link to the DB where the info was taken from

    output$geneinfo_ui <- renderUI({
      row_id <- input$table_annofuse_rows_selected
      message(row_id)
      gene_for_content <- values$annofuse_tbl[row_id, "Gene1A"]
      gene_for_content_2 <- values$annofuse_tbl[row_id, "Gene1B"]


      doublegeneinfo_2_html(gene_for_content, gene_for_content_2)
      # geneinfo_2_html(gene_for_content)
    })
    
    output$geneplots_ui <- renderUI({
      validate(
        need(
          length(input$table_annofuse_rows_selected) > 0,
          message = "Select a row in the table"
        )
      )
      
      tagList(
        h4("Some general info"),
        tabsetPanel(
          tabPanel("Plot left",
                   plotOutput("geneplots_left")),
          tabPanel("Plot right",
                   plotOutput("geneplots_right"))
          
        )
      )
    })
    
    output$geneplots_right <- renderPlot({
      row_id <- input$table_annofuse_rows_selected
      message(row_id)
      # gene_for_content <- values$annofuse_tbl[row_id, "Gene1A"]
      
      fusion_for_content <- values$annofuse_tbl[row_id, "FusionName"]
      rightfused_for_content <- values$annofuse_tbl[row_id, "Gene1B"]
      
      # currently needs some things to replicate the use case situation:
      # 
      ### TODO: these objects below need to be in the R session - in the final
      ### implementation, this should happen seamlessly, and ideally the app could check 
      ### upon starting that these are available
      
      ### These need to be prepped in advance... (e.g. upon starting the app, or
      # in advance before launching it)
      ###### bioMartDataPfam<-readRDS(system.file("extdata","pfamDataBioMart.RDS", package="annoFuse"))
      ###### standardFusioncalls <- values$annofuse_tbl
      ###### annDomain<-annoFuse::getPfamDomain(
      ######   standardFusioncalls  = standardFusioncalls,
      ######   bioMartDataPfam = bioMartDataPfam,
      ######   # partial overlapping domains are retained == "Partial" with keepPartialAnno=TRUE; 
      ######   # if keepPartialAnno=FALSE then domain retained status == "No"
      ######   keepPartialAnno = TRUE)
      ###### # read in exonsToPlot with exon and gene boundaries from gencode.v27.primary_assembly.annotation.gtf.gz
      ###### exons<-readRDS(system.file("extdata", "exonsToPlot.RDS", package = "annoFuse"))
      

      
            # plot BRAF breakpoint in sample for KIAA1549--BRAF fusion
      breakpoints_info <- annDomain$Gene1B[which(annDomain$Gene1B$FusionName==fusion_for_content & annDomain$Gene1B$Gene1B==rightfused_for_content),] %>% 
        dplyr::filter(!is.na(DESC))
      ## Plot breakpoint
      
      plotBreakpoints(domainDataFrame = breakpoints_info,
                      exons = exons,
                      geneposition = "Right") + 
        theme_Publication(base_size = 12)
    })
    
    output$geneplots_left <- renderPlot({
      row_id <- input$table_annofuse_rows_selected
      message(row_id)
      # gene_for_content <- values$annofuse_tbl[row_id, "Gene1A"]
      
      fusion_for_content <- values$annofuse_tbl[row_id, "FusionName"]
      leftfused_for_content <- values$annofuse_tbl[row_id, "Gene1A"]
      
      # plot BRAF breakpoint in sample for KIAA1549--BRAF fusion
      breakpoints_info <- annDomain$Gene1A[which(annDomain$Gene1A$FusionName==fusion_for_content & annDomain$Gene1A$Gene1A==leftfused_for_content),] %>% dplyr::filter(!is.na(DESC))
      ## Plot breakpoint
      
      plotBreakpoints(domainDataFrame = breakpoints_info,
                      exons = exons,
                      geneposition = "Left") + 
        theme_Publication(base_size = 12)
    })


    output$table_annofuse <- DT::renderDataTable({
      DT::datatable(
        values$enhanced_annofuse_tbl,
        style = "bootstrap",
        rownames = FALSE,
        filter = "top",
        selection = "single",
        escape = FALSE,
        options = list(
          scrollX = TRUE,
          pageLength = 25,
          lengthMenu = c(5, 10, 25, 50, 100, nrow(values$enhanced_annofuse_tbl))
        )
      )
    })
    
    
    
    output$af_overview <- renderPlot({
      withProgress({
        plotSummary(values$annofuse_tbl)
      }, message = "Rendering summary...")
    })
    
    # TODO: spinner for when the plot is loading?
    
    output$af_recurrentfusions <- renderPlot({
      ## TODO: ideally we could have the groupby variable even selectable in shiny
      ## i.e. from an input$ widget
      ## it would even make it more robust when specifying the columns - plus, examples would be cool!
      gby_rf <- "Fusion_Type"
      plotn_rf <- 30
      cid_rf <- "Sample"
      palette_rf <- c("blue","green","orange") # I had to specify this 
      plotRecurrentFusions(values$annofuse_tbl, 
                           groupby = gby_rf, 
                           plotn = plotn_rf,
                           countID = cid_rf, 
                           palette_rec = palette_rf)
    })
    
    output$af_recurrentgenes <- renderPlot({
      ## TODO: ideally we could have the groupby variable even selectable in shiny
      ## i.e. from an input$ widget
      ## it would even make it more robust when specifying the columns - plus, examples would be cool!
      gby_rg <- "Fusion_Type"
      plotn_rg <- 30
      cid_rg <- "Sample"
      palette_rg <- c("blue","green","orange") # I had to specify this 
      plotRecurrentGenes(values$annofuse_tbl, 
                         groupby = gby_rg, 
                         plotn = plotn_rg,
                         countID = cid_rg, 
                         palette_rec = palette_rg)
    })
    
      
    # observeEvent(input$interface_overview, {
    #   tour <- read.delim(system.file("extdata", "interface_overview.txt",
    #                                  package = "GeneTonic"),
    #                      sep = ";", stringsAsFactors = FALSE,
    #                      row.names = NULL, quote = "")
    #   rintrojs::introjs(session, options = list(steps = tour))
    # })
  }

  shinyApp(ui = shinyfuse_ui, server = shinyfuse_server)
}



.actionbutton_biocstyle <- "color: #ffffff; background-color: #0092AC"


## TODO: should be able to handle the case where

## ENSEMBL? needs species info



#' Link to NCBI database
#'
#' @param val Character, the gene symbol
#'
#' @return HTML for an action button
#' @noRd
.link2ncbi <- function(val) {
  sprintf(
    '<a href = "http://www.ncbi.nlm.nih.gov/gene/?term=%s[sym]" target = "_blank" class = "btn btn-primary" style = "%s"><i class="fa fa-database"></i>%s</a>',
    val,
    .actionbutton_biocstyle,
    val
  )
}

#' Link to the GTEx Portal
#'
#' @param val Character, the gene symbol
#'
#' @return HTML for an action button
#' @noRd
.link2gtex <- function(val) {
  sprintf(
    '<a href = "https://www.gtexportal.org/home/gene/%s" target = "_blank" class = "btn btn-primary" style = "%s"><i class="fa fa-dna"></i>%s</a>',
    val,
    .actionbutton_biocstyle,
    val
  )
}

#' Link to the Uniprot Portal
#'
#' @param val Character, the gene symbol
#'
#' @return HTML for an action button
#' @noRd
.link2uniprot <- function(val) {
  sprintf(
    '<a href = "https://www.uniprot.org/uniprot/?query=%s&sort=score" target = "_blank" class = "btn btn-primary" style = "%s"><i class="fa fa-spinner"></i>%s</a>',
    val,
    .actionbutton_biocstyle,
    val
  )
}

#' Link to human protein atlas Portal
#'
#' @param val Character, the gene symbol
#'
#' @return HTML for an action button
#' @noRd
.link2hpa <- function(val) {
  sprintf(
    '<a href = "https://www.proteinatlas.org/search/%s" target = "_blank" class = "btn btn-primary" style = "%s"><i class="fa fa-cubes"></i>%s</a>',
    val,
    .actionbutton_biocstyle,
    val
  )
}

.multilink <- function(val) {
  b1 <- sprintf(
    '<a href = "http://www.ncbi.nlm.nih.gov/gene/?term=%s[sym]" target = "_blank" class = "btn btn-primary" style = "%s"><i class="fa fa-database"></i>%s</a>',
    val,
    .actionbutton_biocstyle,
    val
  )
  b2 <- sprintf(
    '<a href = "https://www.gtexportal.org/home/gene/%s" target = "_blank" class = "btn btn-primary" style = "%s"><i class="fa fa-dna"></i>%s</a>',
    val,
    .actionbutton_biocstyle,
    val
  )
  return(paste(b1, b2))
}


geneinfo_2_html <- function(gene_id) {
  gene_ncbi_button <- .link2ncbi(gene_id)
  gene_gtex_button <- .link2gtex(gene_id)
  gene_uniprot_button <- .link2uniprot(gene_id)
  gene_hpa_button <- .link2hpa(gene_id)

  mycontent <- paste0(
    shiny::tags$b(gene_id), shiny::tags$br(),
    "Link to the NCBI Gene database: ", gene_ncbi_button, shiny::tags$br(),
    "Link to the GTEx Portal: ", gene_gtex_button, shiny::tags$br(),
    "Link to the Uniprot Portal: ", gene_uniprot_button, shiny::tags$br(),
    "Link to the Human Protein Atlas: ", gene_hpa_button, shiny::tags$br()
  )

  return(HTML(mycontent))
}



doublegeneinfo_2_html <- function(gene_id1, gene_id2) {
  gene_ncbi_button_1 <- .link2ncbi(gene_id1)
  gene_gtex_button_1 <- .link2gtex(gene_id1)
  gene_uniprot_button_1 <- .link2uniprot(gene_id1)
  gene_hpa_button_1 <- .link2hpa(gene_id1)

  gene_ncbi_button_2 <- .link2ncbi(gene_id2)
  gene_gtex_button_2 <- .link2gtex(gene_id2)
  gene_uniprot_button_2 <- .link2uniprot(gene_id2)
  gene_hpa_button_2 <- .link2hpa(gene_id2)

  # mycontent <- paste0(
  #   shiny::tags$table(
  #     shiny::tags$tr(
  #       shiny::tags$td(width = "50%", gene_ncbi_button_1),
  #       shiny::tags$td(width = "50%", gene_ncbi_button_2)
  #     ),
  #     shiny::tags$tr(
  #       shiny::tags$td(width = "50%", gene_gtex_button_1),
  #       shiny::tags$td(width = "50%", gene_gtex_button_2)
  #     ),
  #     shiny::tags$tr(
  #       shiny::tags$td(width = "50%", gene_uniprot_button_1),
  #       shiny::tags$td(width = "50%", gene_uniprot_button_2)
  #     ),
  #     shiny::tags$tr(
  #       shiny::tags$td(width = "50%", gene_hpa_button_1),
  #       shiny::tags$td(width = "50%", gene_hpa_button_2)
  #     )
  #   )
  # )

  mycontent <- paste0(
    "<table><tr>",
    '<td width="33%">', "NCBI", "</td>",
    '<td width="33%">', gene_ncbi_button_1, "</td>",
    '<td width="33%">', gene_ncbi_button_2, "</td>",
    "</tr></table>",

    "<table><tr>",
    '<td width="33%">', "GTEx", "</td>",
    '<td width="33%">', gene_gtex_button_1, "</td>",
    '<td width="33%">', gene_gtex_button_2, "</td>",
    "</tr></table>",

    "<table><tr>",
    '<td width="33%">', "UniProt", "</td>",
    '<td width="33%">', gene_uniprot_button_1, "</td>",
    '<td width="33%">', gene_uniprot_button_2, "</td>",
    "</tr></table>",

    "<table><tr>",
    '<td width="33%">', "HPA", "</td>",
    '<td width="33%">', gene_hpa_button_1, "</td>",
    '<td width="33%">', gene_hpa_button_2, "</td>",
    "</tr></table>"
  )

  #     shiny::tags$tr(
  #       shiny::tags$td(width = "50%", gene_ncbi_button_1),
  #       shiny::tags$td(width = "50%", gene_ncbi_button_2)
  #     ),
  #     shiny::tags$tr(
  #       shiny::tags$td(width = "50%", gene_gtex_button_1),
  #       shiny::tags$td(width = "50%", gene_gtex_button_2)
  #     ),
  #     shiny::tags$tr(
  #       shiny::tags$td(width = "50%", gene_uniprot_button_1),
  #       shiny::tags$td(width = "50%", gene_uniprot_button_2)
  #     ),
  #     shiny::tags$tr(
  #       shiny::tags$td(width = "50%", gene_hpa_button_1),
  #       shiny::tags$td(width = "50%", gene_hpa_button_2)
  #     )
  #   )
  # )


  return(HTML(mycontent))
}


## TODO? these functions can be directly used also to simply create a report, whose content is just the table with the nicely enhanced many buttons and co.!
