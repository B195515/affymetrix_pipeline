# R Shiny Interactive Volcano plots
# By B195515 - 2022
# Functional Genomic Technologies, SEM2, 2022
######################
# Code modified from https://github.com/FredHutch/interactiveVolcano
# and https://2-bitbio.com/2017/12/clickable-volcano-plots-in-shiny.html

library(shiny)
library(ggplot2)

# UI------------------
ui <- shinyUI(fluidPage(
    
    tabsetPanel(
        
        tabPanel(
            "Plot",
            h1("GSE49448: Limma DEGs"),
            
            sidebarLayout(
                sidebarPanel(
                    width=3,
                    
                    sliderInput("pval_thres",
                                "Set significance threshold",
                                min=0, max=5, value=3, step=0.1),
                    
                    sliderInput("fc_thres",
                                "Set FC threshold",
                                min=-8, max=8, value=4, step=0.1)

                    ), # end sidebarPanel
            
                mainPanel(
                    plotOutput('volcanoPlot',
                               width="100%",
                               height="600px",
                               click='plot_click',
                               brush='plot_brush',
                               hover='plot_hover'), #end plotOutput
                    
                    downloadButton('download_plot', 'Download volcano plot as PDF'),
                    br(),
                    br(),
                    # clicked points table
                    tableOutput('clickedPoints')
                    ) # end mainPanel
                ) # end sidebarLayout
            ), # end Plot panel
        
        tabPanel(
            "Data",
            h1("About the data"),
             
            sidebarLayout(
                sidebarPanel(
                    width=3,
                    ("say something here")
                    ), # end sidebarPanel
                 
                mainPanel(
                    dataTableOutput("all_data")
                    ) # end mainPanel
                 
                ) # end sidebarLayout
            ) # end Data panel
        ) # end tabsetPanel
    ) # end fluidPage
) # end UI


# Server--------------
server <- function(input, output, session){
    
    # Read in data
    dataFrame <- reactive({
        read.table("limma_transformed.csv", header=T, sep=',')
        })
    
    # Subset by threshold
    is_de <- reactive ({
        abs(dataFrame()[["logFC"]]) >= input$fc_thres & 
            dataFrame()[["minus_log10_Pval"]] >= input$pval_thres
        })
    
    # Factorize DE points
    de_vec <- reactive({
        factor(is_de(), levels=c("TRUE","FALSE"))
    })
    
    # Make volcano plot
    reactive_plot <- reactive ({
        ggplot(dataFrame(), aes(x=logFC, y=minus_log10_Pval)) +
        geom_point(alpha=.6, aes(color=de_vec())) +
        coord_cartesian() + 
        geom_hline(yintercept=input$pval_thres, 
                   linetype="dashed", col="grey", size=1) +
        geom_vline(xintercept=c(input$fc_thres, -input$fc_thres), 
                   linetype="dashed", col="grey", size=1) + 
        xlab("log2 fold change") +
        ylab("-log10(FDR)") +
        labs(color="Differentially expressed")
        })

    # Show volcano plot    
    output$volcanoPlot <- renderPlot ({
        reactive_plot()
        })
    
    # Retrieve clicked points
    clicked <- reactive({
        # define the x and y variables:
        nearPoints(dataFrame(), 
                   input$plot_click, 
                   xvar="logFC", 
                   yvar="minus_log10_Pval")
        })
    
    output$clickedPoints <- renderTable({
        clicked()
        }, rownames=T)
    
    # Show all data in Data tab
    output$all_data <- renderDataTable(
        dataFrame()
        )
    
    # Setup plot download
    output$download_volcano <- downloadHandler(
        filename = function() {
            paste0("volcano-plot-", Sys.Date(), ".pdf")
        },
        content = function(file) {
            ggsave(file, reactive_volcano(), device="pdf",
                   width=10, height=5, units="in")
        }
    ) # close downloadHandler
} # close server

# Run the application 
shinyApp(ui = ui, server = server)
