# Code adapted from class code provided in 
# Functional Genomic Technologies 2022 
# by Simon Tomlinson

library(shiny)
library(SummarizedExperiment)
library(org.Hs.eg.db)
library(pheatmap)
library(shinyjs)

#load the example data
load("expression.Rdata")
load("experiment.Rdata")

# Define UI for application that draws a histogram
ui <- fluidPage(

    titlePanel("GSE49448: Heatmap of Top 50 differential gene expression"),

    sidebarLayout(
        sidebarPanel(
            sliderInput("font_row",
                        "Font size row:",
                        min = 6,
                        max = 14,
                        value = 10),
            sliderInput("font_col",
                        "Font size col:",
                        min = 6,
                        max = 14,
                        value = 10),
            
            # Select annotation
            selectInput("select", 
                "Select annotation", 
                choices=c("none",colnames(experiment)), 
                selected = "Group", multiple = T, selectize = TRUE),
            
            checkboxInput("srownames", "Show Row Names", TRUE),
            
            checkboxInput("logtansform", "Log transform values", FALSE),
            
            radioButtons("norm", "Scale by", 
                choices=c("none","row","column"))
            
            ), # end sidebarPanel

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot", height="800")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {

    output$distPlot <- renderPlot(
        # First do the expression that generates the plot
        {
            if(input$logtansform){
                expression <- log2(expression + 1)
                }
       
           if(is.null(input$select)){
                mysel<-NULL
                } else if(input$select[1]=="none"){
                mysel<-NULL
                }else if(length(input$select)==1){
               # Convert df to factor if only has one column
               # Force df type, restore row+col names
                mysel <-as.data.frame(experiment[,input$select[1]])
                rownames(mysel) <-rownames(experiment)
                colnames(mysel) <-input$select[1]
                }else{
                mysel<-experiment[,input$select]
                }
        
        pheatmap(expression,
                fontsize_row = input$font_row,
                fontsize_col = input$font_col,
                show_rownames = input$srownames,
                scale = input$norm,
                annotation_col = mysel)
            }, alt="Heatmap of top 50 expressed genes", execOnResize = F)
    
    observeEvent(input$refresh, {
        session$invalidate
        })
        

}

# Run the application 
shinyApp(ui = ui, server = server)
