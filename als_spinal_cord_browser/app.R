
# ALS Spinal Cord Gene Expression Browser - ASGB
library(shiny)
library(bslib)
load("tpm_metadata.RData")

source("plots.R")

# Define UI for application that draws a histogram
ui <- navbarPage(
    theme = bs_theme(version = 5, bootswatch = "minty"),
    #theme = "layout.css",
    id = "navBarPage",
    windowTitle = "ALS Spinal Cord Gene Expression",
    title = "NYGC ALS Spinal Cord Browser",
    collapsible = FALSE,
    inverse = FALSE,
    #position = "fixed-top"
    # Application title
    #titlePanel("ALS Spinal Cord Gene Expression"),
    tabPanel("Single gene results",
             div(
                 sidebarLayout(
                     sidebarPanel(
                         verticalLayout(fluid=TRUE,
                                        textInput(inputId = "de_gene",value = "CHIT1", label = "Pick a gene"),
                                        actionButton(inputId = "de_log_button", label = "log scale" )
                                        
                         )
                     ),
                     mainPanel(
                         fixedRow(
                             column(5,
                                    #h2(tags$i(textOutput(outputId = "gene_title", inline = TRUE))),
                                    h5(style="text-align: center;",
                                       "ALS vs Control Differential Expression"),
                                    plotOutput("dePlot"),
                                    
                                    DT::dataTableOutput('deTable') 
                             ),
                             column(5,
                                    h5(style="text-align: center;", "Association with ALS duration"),
                                    plotOutput("durPlot"),
                                    
                                    DT::dataTableOutput('durTable')
                             )
                             #downloadButton("downloadPlot", label = "Save plot", class = NULL)
                         )
                     )
                 )
             )
                                        
             ),
   # tabPanel("Gene sets"),
    tabPanel("Custom Plotting",
            #tags$style(type="text/css", href = "layout.css"),
            div(
            sidebarLayout(
                 sidebarPanel(
                     verticalLayout(fluid=TRUE,
                        textInput(inputId = "gene",value = "CHIT1", label = "Pick a gene"),
                        selectInput(inputId = "y_axis", choices = choice_df$long, selected = "Gene", label = "Plot on Y" ),
                        selectInput(inputId = "x_axis", choices = choice_df$long, selected = "Disease", label = "Plot on X" ),
                        selectInput(inputId = "colour", choices = choice_df$long, selected = "Disease", label = "Colour points by" ),
                        checkboxGroupInput(inputId = "tissues", choices = c("Cervical", "Lumbar", "Thoracic"), selected = c("Cervical", "Lumbar"), label = "Tissue"  ),
                        
                        actionButton(inputId = "boxplot_button", label = "Add boxplots" ),
                        actionButton(inputId = "corline_button", label = "Add trend line" ),
                        actionButton(inputId = "stats_button", label = "Add stats" ),
                        actionButton(inputId = "log_button", label = "log2 scale" )
                     )
                 ),
                 
                 # Show a plot of the generated distribution
                 mainPanel(
                     plotOutput("customPlot" ),
                     downloadButton("downloadPlot", label = "Save plot", class = NULL)
                 )
             )
            )
    ),
   tabPanel("About",
            div(
                p("Written by ", tags$a(href="https://jackhump.github.io", "Jack Humphrey.") ),
                p("All data and metadata hosted on ", tags$a(href="https://zenodo.org/record/6385747", "Zenodo.")),
                p("Manuscript published at ", tags$a(href="https://www.nature.com/articles/s41593-022-01205-3", "Nature Neuroscience.", target = "_blank")),
                p("Preprint at ", tags$a(href="https://www.medrxiv.org/content/10.1101/2021.08.31.21262682v1", "medRxiv.", target = "_blank")),
                p("Source code hosted on ", tags$a(href="https://github.com/jackhump/ALS_SpinalCord_QTLs/tree/master/als_spinal_cord_browser", "GitHub.", target = "_blank"),
                  " Please raise any feature suggestions or bugs as an issue.")
            )
   )
)

# Define server logic required to draw plot
server <- function(input, output, session) {
    ## toggle buttons
    
    boxplot <- reactiveVal(FALSE)
    corline <- reactiveVal(FALSE)
    stats <- reactiveVal(FALSE)
    logs <- reactiveVal(TRUE)
    
    
    observeEvent(input$boxplot_button, {
        boxplot(!boxplot())
    })
    
    observeEvent(input$corline_button, {
        corline(!corline())
    })
    
    observeEvent(input$stats_button, {
        stats(!stats())
    })
    
    observeEvent(input$log_button, {
        logs(!logs())
    })
    
    observeEvent(input$de_log_button, {
        logs(!logs())
    })
    
    output$gene_title <- renderText(input$de_gene)
    
    output$customPlot <- renderPlot({
        gene_plot(input$gene, 
                  counts = tpm_df,
                  meta = metadata,
                  x = input$x_axis, 
                  y = input$y_axis, 
                  colourby = input$colour, 
                  tissues = input$tissues, 
                  boxplot = boxplot(), 
                  corline = corline(), 
                  stats = stats(), 
                  log = logs() )  
    })
    
    output$dePlot <- renderPlot({
        gene_plot(input$de_gene, 
                  counts = tpm_df,
                  meta = metadata,
                  x = "disease", 
                  y = "TPM", 
                  colourby = "disease", 
                  tissues = c("Cervical", "Lumbar"), 
                  boxplot = TRUE, 
                  corline = FALSE, 
                  stats = FALSE, 
                  log = logs() )  
    })
    
    output$durPlot <- renderPlot({
        gene_plot(input$de_gene, 
                  counts = tpm_df,
                  meta = metadata,
                  x = "disease_duration", 
                  y = "TPM", 
                  colourby = "disease", 
                  tissues = c("Cervical", "Lumbar"), 
                  boxplot = FALSE, 
                  corline = TRUE, 
                  stats = FALSE, 
                  log = logs() )  
    })
    
    output$deTable <- DT::renderDataTable(#rownames = FALSE,
        gene_table(mygene = input$de_gene,
                   table = de_res, counts = tpm_df ),
        options = list(paging = FALSE, searching = FALSE, 
                       columnDefs = list(list(className = 'dt-center', targets = 0:4)),
                       ordering = FALSE, lengthChange= FALSE, scrollX = FALSE, scrollY = FALSE, info = FALSE)
    )
    
    output$durTable <- DT::renderDataTable(#rownames = FALSE,
        gene_table(mygene = input$de_gene,
                   table = dur_res, counts = tpm_df ),
        options = list(paging = FALSE, searching = FALSE, ordering = FALSE, 
                       columnDefs = list(list(className = 'dt-center', targets = 0:4)),
                       lengthChange= FALSE, scrollX = FALSE, scrollY = FALSE, info = FALSE)
    )
    
    output$downloadPlot <- downloadHandler(
        filename = function() { 'plot.pdf' },
        content = function(file) {
            ggplot2::ggsave(file,
                            plot = gene_plot(input$gene,
                                              counts = tpm_df, 
                                              meta = metadata,
                                              x = input$x_axis, 
                                              y = input$y_axis, 
                                              colourby = input$colour, 
                                              tissues = input$tissues, 
                                              boxplot = boxplot(), 
                                              corline = corline(), 
                                              stats = stats(), 
                                              log = logs() ) , 
                            device = "pdf" )
        }
    )
    
}

# Run the application 
shinyApp(ui = ui, server = server)
