library(shiny)
library(shinydashboard)
library(shinyjs)
library(Seurat)
library(shinydashboardPlus)
library(shinyWidgets)
library(dplyr)

options(shiny.maxRequestSize = 500 * 1024^2)

source('config.R')

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "scRNA-Seq Explorer"),
  dashboardSidebar(
    tags$head(
      tags$style(HTML(
        ".skin-blue .main-header .logo { background-color: #008080; }
         .skin-blue .main-header .navbar { background-color: #008080; }
         .skin-blue .sidebar { background-color: #004d4d; }"
      ))
    ),
    sidebarMenu(id='tab',
                useShinyjs(),
                menuItem("Home", tabName = "home", icon = icon("home")),
                menuItem("Analysis", tabName = "input", icon = icon("flask")),
                conditionalPanel(condition = "input.tab == 'input'",
                                 div(
                                   fileInput("file", "Upload Seurat RDS File", 
                                             multiple=FALSE, 
                                             accept=c('.rds'),
                                             placeholder = "Select Seurat Object"),
                                   actionButton("reset", "Reset", icon = icon("refresh"), 
                                                style = "color: #fff; background-color: #d9534f; width: 87%"),
                                   actionButton("run", "Analyze", icon = icon("play"), 
                                                style = "color: #fff; background-color: #5cb85c; width: 87%")
                                 )
                )
    )
  ), 
  dashboardBody(
    tabItems(
      tabItem(tabName = "input",
              tabsetPanel(id = 'main_tabs',
                          tabPanel("Getting Started",
                                   div(style = "padding: 15px;",
                                       h3("How to Use This App"),
                                       tags$ul(
                                         tags$li("Upload a Seurat object using 'Upload File'."),
                                         tags$li("Press the 'Run' button to analyze the data."),
                                         tags$li("New tabs will appear to explore data visualization."),
                                         tags$li("Press 'Reset' to remove files and clear tabs.")
                                       )
                                   )
                          )
              )
      ),
      tabItem(tabName = "home",
              box(
                title = "Welcome to scRNA-Seq Explorer", 
                status = "primary", 
                solidHeader = TRUE,
                width = 12,
                HTML("<p>This interactive tool allows you to analyze and visualize single-cell RNA sequencing data.</p>
                     <p>Upload a Seurat object RDS file and click 'Analyze' to get started.</p>")
              )
      )
    )
  )         
)

server <- function(input, output, session) {
  options(shiny.maxRequestSize=300*1024^2)
  
  values <- reactiveValues()
  
  # Disable Run by default
  shinyjs::disable("run")
  
  observe({
    if(is.null(input$file) != TRUE) {
      shinyjs::enable("run")
    } else {
      shinyjs::disable("run")
    }
  })
  
  observeEvent(input$run, {
    shinyjs::disable("run")
    
    # Clear tabs before 'Run' is ran another time
    removeTab("main_tabs", "UMAP")
    removeTab("main_tabs", "Gene Expression")
    
    show_modal_spinner(text = "Preparing plots...")
    
    obj <- load_seurat_obj(input$file$datapath)
    if (is.vector(obj)){
      showModal(modalDialog(
        title = "Error with file",
        HTML("<h5>There is an error with the file you uploaded. See below for more details.</h5><br>",
             paste(unlist(obj), collapse = "<br><br>"))
      ))
      shinyjs::enable("run")
      
    } else {
      
      output$umap <- renderPlot({
        if (!is.null(input$metadata_col)) {
          create_metadata_UMAP(obj, input$metadata_col)
        }
      })
      
      output$featurePlot <- renderPlot({
        if (!is.null(input$gene)) {
          create_feature_plot(obj, input$gene)
        }
      })
      
      output$downloadFeaturePlot <- downloadHandler(
        filename = function(){
          paste0(input$gene, '_feature_plot', '.png')
        },
        content = function(file){
          plot <- create_feature_plot(obj, input$gene)
          ggsave(filename=file, width = 10, height = 5, type = "cairo")
        }
      )
      output$download_umap <- downloadHandler(
        filename = function(){
          paste0(input$metadata_col, '_UMAP', '.png')
        },
        content = function(file){
          plot <- create_metadata_UMAP(obj, input$metadata_col)
          ggsave(filename=file, width = 8, height = 5, type = "cairo")
        }
      )
      
      insertTab(
        inputId = "main_tabs",
        tabPanel(
          "UMAP",
          fluidRow(
            column(
              width = 8,
              plotOutput(outputId = 'umap'),
              downloadButton("download_umap", "Download UMAP")
            ),
            column(
              width = 4,
              selectizeInput("metadata_col", 
                             "Metadata Column", 
                             colnames(obj@meta.data)
              )
            )
          ),
          style = "height: 90%; width: 95%; padding-top: 5%;"
        ),
        select = TRUE
      )
      insertTab(
        inputId = "main_tabs",
        tabPanel(
          "Gene Expression",
          fluidRow(
            column(
              width = 8,
              plotOutput(outputId = 'featurePlot'),
              downloadButton("downloadFeaturePlot", "Download Feature Plot")
            ),
            column(
              width = 4,
              selectizeInput("gene", 
                             "Genes", 
                             rownames(obj)
              )
            )
          ),
          style = "height: 90%; width: 95%; padding-top: 5%;"
        )
      )
      
      remove_modal_spinner()
      shinyjs::enable("run")
      
    }
  })
  
  # Clear all sidebar inputs when 'Reset' button is clicked
  observeEvent(input$reset, {
    shinyjs::reset("file")
    removeTab("main_tabs", "UMAP")
    removeTab("main_tabs", "Gene Expression")
    shinyjs::disable("run")
  })
  
}

shinyApp(ui, server)

