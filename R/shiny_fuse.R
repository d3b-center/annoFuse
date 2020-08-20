#' shiny_fuse
#'
#' Exploring interactively the results of the annoFuse pipeline
#'
#' The application can also be started without specifying the location of the data,
#' which can be provided (via upload) at runtime.
#'
#' @param out_annofuse The character string specifying the location of the file
#' output by the annoFuse pipeline. This file needs to be structured with the set
#' of columns required for the later exploration steps in the interactive app.
#'
#' @return A Shiny app object
#' @export
#'
#' @rawNamespace import(shiny, except = c(renderDataTable, dataTableOutput))
#' @import shinydashboard
#' @import rintrojs
#' @import shinythemes
#' @importFrom DT datatable renderDataTable dataTableOutput
#'
#' @examples
#' out_annofuse <- system.file("extdata", "PutativeDriverAnnoFuse_test_v14.tsv", package = "annoFuse")
#' if (interactive()) {
#'   shiny_fuse(out_annofuse)
#' }
shiny_fuse <- function(out_annofuse = NULL) {

  # Checks on the objects provided ---------------------------------------------

  if (!is.null(out_annofuse)) {
    if (!is(out_annofuse, "character")) {
      stop("'out_annofuse' has to be a character string")
    }
    if (!file.exists(out_annofuse)) {
      stop("File specified by 'out_annofuse' not found")
    }
  }
  
  oopt <- options(spinner.type = 6, spinner.color = "#0092AC")
  # play nice with other previously chosen options
  on.exit(options(oopt))


  ### TODO: maybe check here that pfam and exons objects are available?
  ### "slight issue": it takes a while to load, so maybe do this in advance? On the server,
  ### it would still need to be done at each session
  ### NOTE: this is not optimal, but it is to give an idea of how it could be ;)
  # if(!exists("bioMartDataPfam")) {
  #   message("Loading pfam data...")
  #   bioMartDataPfam <- readRDS(system.file("extdata","pfamDataBioMart.RDS", package="annoFuse"))
  # }
  # read in exonsToPlot with exon and gene boundaries from gencode.v27.primary_assembly.annotation.gtf.gz
  # if(!exists("exons")) {
  #   message("Loading exons data...")
  #   exons <- readRDS(system.file("extdata", "exonsToPlot.RDS", package = "annoFuse"))
  # }

  # UI definition -----------------------------------------------------------
  shinyfuse_ui <- shinydashboard::dashboardPage(
    skin = "black",

    # header definition -------------------------------------------------------
    header = shinydashboard::dashboardHeader(
      title = "shiny_fuse",
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
      
      actionButton(inputId = "btn_load_demo",
                   label = "Load demo data"),

      uiOutput("choose_annofusedata_file"),
      # shinydashboard::menuItem(
      #   text = "Input Settings", icon = icon("cog"),
      #   startExpanded = TRUE,
      #   numericInput(
      #     inputId = "whatevs",
      #     label = "number of genesets",
      #     value = 15, min = 1, max = 50
      #   )
      # )
      uiOutput("plot_controls")
    ),

    # body definition ---------------------------------------------------------
    body = shinydashboard::dashboardBody(
      rintrojs::introjsUI(),
      ## Define output size and style of error messages
      shiny::tags$head(
        shiny::tags$style(
          shiny::HTML(".shiny-output-error-validation {
                 font-size: 15px;
                 color: forestgreen;
                 text-align: center;
                 }
                 ")
        )
      ),
      id = "main-app",
      navbarPage(
        title = div(),
        windowTitle = "shinyfuse",
        footer = "",
        theme = shinytheme("cosmo"),
        selected = "TableExplorer",

        # ui TableExplorer -----------------------------------------------------
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
                fluidRow(
                  column(
                    width = 4,
                    offset = 4,
                    actionButton("btn_load_exonsdata", "load exons")
                  ),
                  column(
                    width = 4,
                    actionButton("btn_load_pfamdata", "load pfam")
                  )
                ),
                h4("Some content, for example linked to the selected row"),
                uiOutput("geneinfo_ui"),
                uiOutput("geneplots_ui")
              )
            )
          )
        ),

        # ui TableSummary -----------------------------------------------------
        tabPanel(
          title = "TableSummary", icon = icon("dna"),
          fluidPage(
            fluidRow(
              column(
                width = 12,
                box(
                  title = "annoFuse summary",
                  status = "success",
                  width = 12,
                  collapsible = TRUE,
                  collapsed = TRUE,
                  withSpinner(
                    plotOutput("af_overview")
                  ),
                  downloadButton("btn_dl_summary", label = "", 
                                 class = "btn btn-success")
                )
              )
            ),
            fluidRow(
              column(
                width = 6,
                withSpinner(
                  plotOutput("af_recurrentfusions")
                ),
                downloadButton("btn_dl_recufusions", label = "", 
                               class = "btn btn-success")
              ),
              column(
                width = 6,
                withSpinner(
                  plotOutput("af_recurrentgenes")
                ),
                downloadButton("btn_dl_recugenes", label = "",
                               class = "btn btn-success")
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
    values$ann_domain <- NULL

    values$data_exons <- NULL
    values$data_pfam <- NULL
    # currently needs some things to replicate the use case situation:
    #
    ### TODO: these objects below need to be in the R session - in the final
    ### implementation, this should happen seamlessly, and ideally the app could check
    ### upon starting that these are available

    ### These need to be prepped in advance... (e.g. upon starting the app, or
    # in advance before launching it)




    # Define data file if annoFuse data is not provided ------------------------
    if (is.null(out_annofuse)) {
      output$choose_annofusedata_file <- renderUI({
        menuItem(
          icon = icon("file-alt"),
          fileInput(
            inputId = "annofusedatasel",
            label = "Load sample data file (tab-separated)",
            accept = c(
              "text/tab-separated-values", "text/plain",
              ".tsv", ".tab", ".txt"
            ),
            multiple = FALSE
          )
        )
      })
    } else {
      annofuse_tbl <- read.delim(out_annofuse)
      annofuse_tbl <- .check_annoFuse_calls(annofuse_tbl)
      values$annofuse_tbl <- annofuse_tbl

      enhanced_annofuse_tbl <- annofuse_tbl
      # enhancing the content of the table
      enhanced_annofuse_tbl$Gene1A <- .multilink(enhanced_annofuse_tbl$Gene1A)
      enhanced_annofuse_tbl$Gene1B <- .multilink(enhanced_annofuse_tbl$Gene1B)
      values$enhanced_annofuse_tbl <- enhanced_annofuse_tbl
    }

    # Controls for plot panels -------------------------------------------------
    output$plot_controls <- renderUI({
      if (is.null(values$annofuse_tbl)) {
        return(NULL)
      } else {
        tagList(
          selectInput(
            inputId = "af_cols",
            label = "Grouping column",
            choices = c("", colnames(values$annofuse_tbl)),
            selectize = TRUE, multiple = FALSE, selected = "Fusion_Type"
          ),
          numericInput(
            inputId = "af_n_topfusions",
            label = "Nr top fusions",
            value = 20, min = 1, max = 300, step = 1
          ),
          selectInput(
            inputId = "af_countcol",
            label = "Counting column",
            choices = c("", colnames(values$annofuse_tbl)),
            selectize = TRUE, multiple = FALSE, selected = "Sample"
          ),
        )
      }
    })

    # Load annoFuse data file --------------------------------------------------
    observeEvent(input$annofusedatasel, {
      message("Reading in...")
      values$annofuse_tbl <- .check_annoFuse_calls(
        read.delim(input$annofusedatasel$datapath)
      )
      values$enhanced_annofuse_tbl <- values$annofuse_tbl

      # enhancing the content of the table
      values$enhanced_annofuse_tbl$Gene1A <- .multilink(values$enhanced_annofuse_tbl$Gene1A)
      values$enhanced_annofuse_tbl$Gene1B <- .multilink(values$enhanced_annofuse_tbl$Gene1B)

      if (!is.null(values$data_pfam)) {
        message("Creating domain information...")
        values$ann_domain <- annoFuse::get_Pfam_domain(
          standardFusioncalls = values$annofuse_tbl,
          bioMartDataPfam = values$data_pfam, # must be pre-loaded
          # partial overlapping domains are retained == "Partial" with keepPartialAnno=TRUE;
          # if keepPartialAnno=FALSE then domain retained status == "No"
          keepPartialAnno = TRUE
        )
      }
    })
    
    # Load demo data
    observeEvent(input$btn_load_demo, {
      message("Loading demo data...")
      demodata_location <- system.file("extdata", "PutativeDriverAnnoFuse_test_v14.tsv", package = "annoFuse")
      values$annofuse_tbl <- 
        .check_annoFuse_calls(read.delim(demodata_location))
      values$enhanced_annofuse_tbl <- values$annofuse_tbl
      
      # enhancing the content of the table
      values$enhanced_annofuse_tbl$Gene1A <- .multilink(values$enhanced_annofuse_tbl$Gene1A)
      values$enhanced_annofuse_tbl$Gene1B <- .multilink(values$enhanced_annofuse_tbl$Gene1B)
      
      if (!is.null(values$data_pfam)) {
        message("Creating domain information...")
        values$ann_domain <- annoFuse::get_Pfam_domain(
          standardFusioncalls = values$annofuse_tbl,
          bioMartDataPfam = values$data_pfam, # must be pre-loaded
          # partial overlapping domains are retained == "Partial" with keepPartialAnno=TRUE;
          # if keepPartialAnno=FALSE then domain retained status == "No"
          keepPartialAnno = TRUE
        )
      }
    })
    
    

    # Main interactive table for exploration -----------------------------------
    output$table_annofuse <- DT::renderDataTable({
      validate(
        need(
          !is.null(values$annofuse_tbl),
          "Please upload the results of annoFuse to start the exploration"
        )
      )

      DT::datatable(
        values$enhanced_annofuse_tbl,
        style = "bootstrap",
        rownames = FALSE,
        filter = "top",
        selection = "single",
        escape = FALSE,
        extensions = c("Buttons"),
        options = list(
          scrollX = TRUE,
          pageLength = 25,
          lengthMenu = c(5, 10, 25, 50, 100, nrow(values$enhanced_annofuse_tbl)),
          dom = "Bfrtip",
          buttons = list("copy", "print", list(
            extend = "collection",
            buttons = c("csv", "excel", "pdf"),
            text = "Download"
          ))
        )
      )
    })

    # TODO? link to the DB where the info was taken from

    # Gene info and plots from TableExplorer
    output$geneinfo_ui <- renderUI({
      validate(
        need(
          length(input$table_annofuse_rows_selected) > 0,
          "Please select a row to display the genes info"
        )
      )

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
          message = ""
        )
      )

      tagList(
        h4("Some general info"),
        tabsetPanel(
          tabPanel(
            "Plot left",
            withSpinner(
              plotOutput("geneplots_left"),
              type = 6, color = "#0092AC"
            ),
            downloadButton("btn_dl_bpleft", label = "", 
                           class = "btn btn-success")
          ),
          tabPanel(
            "Plot right",
            withSpinner(
              plotOutput("geneplots_right"),
              type = 6, color = "#0092AC"
            ),
            downloadButton("btn_dl_bpright", label = "", 
                           class = "btn btn-success")
          ),
          tabPanel(
            "Plot both",
            withSpinner(
              plotOutput("geneplots_both"),
              type = 6, color = "#0092AC"
            ),
            downloadButton("btn_dl_bpboth", label = "", 
                           class = "btn btn-success")
          )
          
        )
      )
    })

    # Breakpoint plots ---------------------------------------------------------
    output$geneplots_right <- renderPlot({
      validate(
        need(
          (!is.null(values$data_exons) & !is.null(values$data_pfam)),
          "Please load the exons and the pfam information via the buttons above to display the plot"
        )
      )

      row_id <- input$table_annofuse_rows_selected
      message(row_id)
      # gene_for_content <- values$annofuse_tbl[row_id, "Gene1A"]

      fusion_for_content <- values$annofuse_tbl[row_id, "FusionName"]
      rightfused_for_content <- values$annofuse_tbl[row_id, "Gene1B"]

      # plot BRAF breakpoint in sample for KIAA1549--BRAF fusion
      # breakpoints_info <- values$ann_domain$Gene1B[which(values$ann_domain$Gene1B$FusionName == fusion_for_content & values$ann_domain$Gene1B$Gene1B == rightfused_for_content), ] %>%
      #   dplyr::filter(!is.na(.data$DESC))
      ## Plot breakpoint

      p <- plot_breakpoints(
        domainDataFrame = values$ann_domain,
        exons = values$data_exons,
        geneposition = "Right",
        fusionname = fusion_for_content
      ) 
      values$plotobj_breakpoint_right <- p
      print(p)
    })

    output$geneplots_left <- renderPlot({
      validate(
        need(
          (!is.null(values$data_exons) & !is.null(values$data_pfam)),
          "Please load the exons and the pfam information via the buttons above to display the plot"
        )
      )

      row_id <- input$table_annofuse_rows_selected
      message(row_id)
      # gene_for_content <- values$annofuse_tbl[row_id, "Gene1A"]

      fusion_for_content <- values$annofuse_tbl[row_id, "FusionName"]
      leftfused_for_content <- values$annofuse_tbl[row_id, "Gene1A"]

      # plot BRAF breakpoint in sample for KIAA1549--BRAF fusion
      # breakpoints_info <- values$ann_domain$Gene1A[which(values$ann_domain$Gene1A$FusionName == fusion_for_content & values$ann_domain$Gene1A$Gene1A == leftfused_for_content), ] %>% dplyr::filter(!is.na(.data$DESC))
      ## Plot breakpoint

      p <- plot_breakpoints(
        domainDataFrame = values$ann_domain,
        exons = values$data_exons,
        geneposition = "Left",
        fusionname = fusion_for_content
      )
      values$plotobj_breakpoint_left <- p
      print(p)
    })
    
    output$geneplots_both <- renderPlot({
      validate(
        need(
          !is.null(values$plotobj_breakpoint_left) & !is.null(values$plotobj_breakpoint_right),
          "Please load the exons and the pfam information via the buttons above to display the plot"
        )
      )  
      pboth <- ggpubr::ggarrange(
        values$plotobj_breakpoint_left, 
        values$plotobj_breakpoint_right, 
        align = "h")
      values$plotobj_breakpoint_both <- pboth
      print(pboth)      
    })

    # Content for TableSummary panel -------------------------------------------

    output$af_overview <- renderPlot({
      withProgress(
        {
          p <- plot_summary(values$annofuse_tbl)
        },
        message = "Rendering summary..."
      )
      values$plotobj_summary <- p
      print(p)
    })

    # TODO: spinner for when the plot is loading?

    output$af_recurrentfusions <- renderPlot({
      validate(
        need(
          !is.null(values$annofuse_tbl),
          "Please provide the results of annoFuse to display the plot"
        )
      )

      gby_rf <- input$af_cols
      plotn_rf <- input$af_n_topfusions
      cid_rf <- input$af_countcol
      palette_rf <- c("blue", "green", "orange") # I had to specify this
      p <- plot_recurrent_fusions(values$annofuse_tbl,
        groupby = gby_rf,
        plotn = plotn_rf,
        countID = cid_rf,
        palette_rec = palette_rf
      )
      values$plotobj_recufusions <- p
      print(p)
    })

    output$af_recurrentgenes <- renderPlot({
      validate(
        need(
          !is.null(values$annofuse_tbl),
          "Please provide the results of annoFuse to display the plot"
        )
      )

      gby_rg <- input$af_cols
      plotn_rg <- input$af_n_topfusions
      cid_rg <- input$af_countcol
      palette_rg <- c("blue", "green", "orange") # I had to specify this
      p <- plot_recurrent_genes(values$annofuse_tbl,
        groupby = gby_rg,
        plotn = plotn_rg,
        countID = cid_rg,
        palette_rec = palette_rg
      )
      values$plotobj_recugenes <- p
      print(p)
    })

    # Tour trigger -------------------------------------------------------------
    observeEvent(input$interface_overview, {
      tour <- read.delim(system.file("extdata", "annoFuse_overview.txt",
        package = "annoFuse"
      ),
      sep = ";", stringsAsFactors = FALSE,
      row.names = NULL, quote = ""
      )
      rintrojs::introjs(session, options = list(steps = tour))
    })


    # Load data via buttons ----------------------------------------------------
    observeEvent(input$btn_load_pfamdata, {
      if (is.null(values$data_pfam)) {
        message("Loading pfam data...")
        showNotification("Loading pfam data...", type = "message")
        values$data_pfam <- readRDS(system.file("extdata", "pfamDataBioMart.RDS", package = "annoFuse"))
        showNotification("Done loading pfam data!", type = "message")

        if (!is.null(values$annofuse_tbl)) {
          showNotification("Correctly created domain information!", type = "message")
          values$ann_domain <- annoFuse::get_Pfam_domain(
            standardFusioncalls = values$annofuse_tbl,
            bioMartDataPfam = values$data_pfam, # must be pre-loaded
            # partial overlapping domains are retained == "Partial" with keepPartialAnno=TRUE;
            # if keepPartialAnno=FALSE then domain retained status == "No"
            keepPartialAnno = TRUE
          )
        }
      } else {
        showNotification("pfam data already loaded", type = "default")
      }
    })

    observeEvent(input$btn_load_exonsdata, {
      if (is.null(values$data_exons)) {
        message("Loading exons data...")
        showNotification("Loading exons data...", type = "message", duration = 10)
        values$data_exons <- readRDS(system.file("extdata", "exonsToPlot.RDS", package = "annoFuse"))
        showNotification("Done loading exons data!", type = "message")
      } else {
        showNotification("exons data already loaded", type = "default")
      }
    })
    
    
    # Defining behaviors for downloading the plots -----------------------------
    output$btn_dl_bpleft <- downloadHandler(
      filename = "annofuse_bpleft.pdf",
      content = function(file) {
        ggsave(file, plot = values$plotobj_breakpoint_left #, 
               # width = input$export_width,
               # height = input$export_height, units = "cm"
        )
      })
    
    output$btn_dl_bpright <- downloadHandler(
      filename = "annofuse_bpright.pdf",
      content = function(file) {
        ggsave(file, plot = values$plotobj_breakpoint_right #, 
               # width = input$export_width,
               # height = input$export_height, units = "cm"
        )
      })
    
    output$btn_dl_bpboth <- downloadHandler(
      filename = "annofuse_bpboth.pdf",
      content = function(file) {
        ggsave(file, plot = values$plotobj_breakpoint_both #, 
               # width = input$export_width,
               # height = input$export_height, units = "cm"
        )
      })
    
    output$btn_dl_summary <- downloadHandler(
      filename = "annofuse_summary.pdf",
      content = function(file) {
        ggsave(file, plot = values$plotobj_summary #, 
               # width = input$export_width,
               # height = input$export_height, units = "cm"
        )
      })
    
    output$btn_dl_recufusions <- downloadHandler(
      filename = "annofuse_recurrent_fusions.pdf",
      content = function(file) {
        ggsave(file, plot = values$plotobj_recufusions #, 
               # width = input$export_width,
               # height = input$export_height, units = "cm"
        )
      })
    
    output$btn_dl_recugenes <- downloadHandler(
      filename = "annofuse_recurrent_genes.pdf",
      content = function(file) {
        ggsave(file, plot = values$plotobj_recugenes #, 
               # width = input$export_width,
               # height = input$export_height, units = "cm"
        )
      })
    
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
