#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)

# Define UI for application that draws a histogram
ui <- shinyUI(fluidPage(
    
    h1("GSE49448: Limma DEGs"),
    
    selectInput('inputFile', label='Selected data:', choice='limma_transformed'),
    
    plotOutput('volcanoPlot',
               click='plot_click',
               brush='plot_brush'
               ),
    
    sliderInput('FCCut', label="log2FC cutoff",-8,8, c(0,0), width="400px"),
    sliderInput('Pvalcut', label="Significance", 0,4,c(0,4), width="400px"),

    
    # here the table for the clicked points:
    tableOutput('clickedPoints')

))



# Define server logic required to draw a histogram
server <- function(input, output, session){
    
    #read in the table as a function of the select input
    dataFrame <- reactive({
        filename <- paste0(input$inputFile,".csv")
        read.table(file=filename, header=T, sep=',')
    })
    
    # filter the table by CPM:
    dataFilter <- reactive({
        dataFrame()[dataFrame()$logFC > input$FCCut,]
        dataFrame()[dataFrame()$minus_log10_Pval > input$Pvalcut,]
    })
    
    #plot it normally with ggplot:
    output$volcanoPlot <- renderPlot({ 
        ggplot(dataFilter(),aes(x=logFC, y=minus_log10_Pval, color=sig)) +
            geom_point() +
            coord_cartesian() +
            ylab("-log10(adj.P.Val)") +
            xlab("log2 fold change")
    })
    
    # get the clicked points
    clicked <- reactive({
        # define the x and y variables:
        nearPoints(dataFilter(), input$plot_click, xvar="logFC", yvar="minus_log10_Pval")
    })
    
    
    output$clickedPoints <- renderTable({
        clicked()
    }, rownames=T)
    
    
    # get the brushed area
    brushed <- reactive({
        nearPoints(dataFilter(), input$plot_brush, xvar="logFC", yvar="minus_log10_Pval")
    })
    
    output$brushedPoints <- renderTable({
        brushed()
    }, rownames=T)
}

# Run the application 
shinyApp(ui = ui, server = server)
