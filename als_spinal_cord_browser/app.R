
# ALS Spinal Cord Gene Expression Browser - ASGB
library(shiny)
library(bslib)
load("tpm_metadata.RData")
metadata$site_of_motor_onset <- gsub(" ", "\n",metadata$site_of_motor_onset )
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
    tabPanel("Plot results",
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
                        actionButton(inputId = "log_button", label = "log scale" )
                     )
                 ),
                 
                 # Show a plot of the generated distribution
                 mainPanel(
                     plotOutput("genePlot" ),
                     downloadButton("downloadPlot", label = "Save plot", class = NULL)
                 )
             )
            )
    ),
    tabPanel("About",
             
             p("written by Jack Humphrey."),
             p(tags$a(href="https://zenodo.org/record/6385747", "All counts, TPMs and metadata")),
             p(tags$a(href="https://www.medrxiv.org/content/10.1101/2021.08.31.21262682v1", "Preprint describing results", target = "_blank")),
             p(tags$a(href="https://github.com/jackhump/als_spinal_cord_browser", "Source code.", target = "_blank") )
             
             )
)

# Define server logic required to draw plot
server <- function(input, output) {
    
    ## toggle buttons
    
    boxplot <- reactiveVal(FALSE)
    corline <- reactiveVal(FALSE)
    stats <- reactiveVal(FALSE)
    logs <- reactiveVal(FALSE)
    
    
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
    
    output$genePlot <- renderPlot({
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