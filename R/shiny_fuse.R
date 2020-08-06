

#' shiny_fuse
#'
#' TODO
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
#' # TODO
shiny_fuse <- function() {
  
  # checks on the objects provided
  #?
  
  # UI definition -----------------------------------------------------------
  shinyfuse_ui <- shinydashboard::dashboardPage(
    skin = "black",
    
    # header definition -------------------------------------------------------
    header = shinydashboard::dashboardHeader(
      title = "TODOtitle",
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
      shinydashboard::menuItem(
        text = "SomeSettings", icon = icon("cog"),
        startExpanded = TRUE,
        numericInput(inputId = "whatevs",
                     label = "number of genesets",
                     value = 15, min = 1, max = 50)
      )
    ),
    
    # body definition ---------------------------------------------------------
    body = shinydashboard::dashboardBody(
      rintrojs::introjsUI(),
      # useShinyjs(),
      id = "main-app",
      navbarPage(
        title = div(),
        windowTitle = "shinyfuse",
        footer= "", 
        theme = shinytheme("cosmo"),
        selected = "panel1",
        navbarMenu(
          title = "Welcome", icon = icon("home"),
          tabPanel(
            title = "panel1", icon = icon("file-alt"),
            fluidPage(
              h1("welcome - panel1"),
              h2("shinyfuse - version TODO"),
              h3("Motivation"),
              p("test to that"),
              h3("General info on annofuse"),
              p("a plot and tables"),
              h3("Need help?"),
              h3("Behind the curtain"),
              p("linking to source code"),
              h3("Authors"),
              p("yeah well")
            )
          ),
          tabPanel(
            title = "panel2", icon = icon("vials"),
            fluidPage(
              h1("welcome - panel2"),
              h3("The datasets!"),
              DT::dataTableOutput("table_annofuse")
            )
          )
        ),
        navbarMenu(
          title = "ee", icon = icon("dna"),
          tabPanel(
            title = "epanel1", icon = icon("table"),
            fluidPage(
              h1("ee - panel1")
            )
          ),
          tabPanel(
            title = "epanel2", icon = icon("poll-h"),
            fluidPage(
              h1("ee - panel2")
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
    
    annofuse_out <- read.delim("/Users/fede/Development/annoFuse/PutativeDriverAnnoFuse_test_v14.tsv")
    
    output$table_annofuse <- DT::renderDataTable({
      DT::datatable(annofuse_out,
                    style = 'bootstrap', 
                    rownames = FALSE, 
                    filter = 'top',
                    options = list(
                      scrollX = TRUE,
                      pageLength = 10,
                      lengthMenu = c(5, 10, 25, 50, 100, nrow(annofuse_out))
                    )
      )
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
.link2ncbi <- function(val) {
  sprintf('<a href = "http://www.ncbi.nlm.nih.gov/gene/?term=%s[sym]" target = "_blank" class = "btn btn-primary" style = "%s"><i class="fa fa-database"></i>%s</a>',
          val,
          .actionbutton_biocstyle,
          val)
}

#' Link to the GTEx Portal
.link2gtex <- function(val) {
  sprintf('<a href = "https://www.gtexportal.org/home/gene/%s" target = "_blank" class = "btn btn-primary" style = "%s"><i class="fa fa-dna"></i>%s</a>',
          val,
          .actionbutton_biocstyle,
          val)
}

#' Link to the Uniprot Portal
.link2uniprot <- function(val) {
  sprintf('<a href = "https://www.uniprot.org/uniprot/?query=%s&sort=score" target = "_blank" class = "btn btn-primary" style = "%s"><i class="fa fa-spinner"></i>%s</a>',
          val,
          .actionbutton_biocstyle,
          val)
}

#' Link to the human protein atlas Portal
.link2hpa <- function(val) {
  sprintf('<a href = "https://www.proteinatlas.org/search/%s" target = "_blank" class = "btn btn-primary" style = "%s"><i class="fa fa-cubes"></i>%s</a>',
          val,
          .actionbutton_biocstyle,
          val)
}

.multilink <- function(val) {
  b1 <- sprintf('<a href = "http://www.ncbi.nlm.nih.gov/gene/?term=%s[sym]" target = "_blank" class = "btn btn-primary" style = "%s"><i class="fa fa-database"></i>%s</a>',
          val,
          .actionbutton_biocstyle,
          val)
  b2 <- sprintf('<a href = "https://www.gtexportal.org/home/gene/%s" target = "_blank" class = "btn btn-primary" style = "%s"><i class="fa fa-dna"></i>%s</a>',
                val,
                .actionbutton_biocstyle,
                val)
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


## TODO? these functions can be directly used also to simply create a report, whose content is just the table with the nicely enhanced many buttons and co.!