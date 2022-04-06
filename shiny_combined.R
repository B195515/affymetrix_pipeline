###################### 
# R Shiny Interactive plots
# By B195515 - 2022
# Functional Genomic Technologies, SEM2, 2022
######################
# Volcano plot code modified from https://github.com/FredHutch/interactiveVolcano
# and https://2-bitbio.com/2017/12/clickable-volcano-plots-in-shiny.html
# Heatmap plot code modified from class code provided 
# in Functional Genomic Technologies 2022 
# by Simon Tomlinson
#####################

library(shiny)
library(ggplot2)
library(DT)
library(pheatmap)

load("differential.Rdata")
load("experiment.Rdata")
load("expression.Rdata")

# UI------------------
ui <- shinyUI(fluidPage(
    
    tabsetPanel(
        
        tabPanel(
            "Volcano plot",
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
            "Heatmap",
            h1("GSE49448: Expression heatmap of Top 50 DEGs"),
            sidebarLayout(
                sidebarPanel(
                    width=3,
                    
                    sliderInput("font_row",
                                "Font size row:",
                                min = 6,max = 14,value = 10),
                    
                    sliderInput("font_col",
                                "Font size col:",
                                min = 6,max = 14,value = 10),
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
                
                mainPanel(
                    plotOutput("distPlot", 
                               height="800",
                               width="80%"),
                    
                    downloadButton('download_hm', 'Download heatmap as PDF'),
                    ) # end mainPanel
                ) # end sidebarlayout
            ), # end panel
        
        tabPanel(
            "Limma data",
            h1("About the data"),
             
            sidebarLayout(
                sidebarPanel(
                    width=3,
                    ("say something here")
                    ), # end sidebarPanel
                 
                mainPanel(
                    dataTableOutput("limma_data")
                    ) # end mainPanel
                 
                ) # end sidebarLayout
            ), # end Data panel

        tabPanel(
            "Expression data",
            h1("About the data"),
            
            sidebarLayout(
                sidebarPanel(
                    width=3,
                    ("say something here")
                    ), # end sidebarPanel
                
                mainPanel(
                    dataTableOutput("exp_data")
                    ) # end mainPanel
                    
                ) # end sidebarLayout
            ), # end Data panel
        
        tabPanel(
            "Experiment",
            h1("Experiment details"),
            
            sidebarLayout(
                sidebarPanel(
                    width=3,
                    ("say something here")
                    ), # end sidebarPanel
                
                mainPanel(
                    dataTableOutput("targets")
                    ) # end mainPanel
                ) # end sidebarLayout
            ) # end tabPanel
        ) # end tabsetPanel
    ) # end fluidPage
) # end UI


# Server--------------
server <- function(input, output, session){
    
    # Subset by threshold
    is_de <- reactive ({
        abs(differential[["logFC"]]) >= input$fc_thres & 
            differential[["minus_log10_Pval"]] >= input$pval_thres
        })
    
    # Factorize DE points
    de_vec <- reactive({
        factor(is_de(), levels=c("TRUE","FALSE"))
    })
    
    # Make volcano plot
    reactive_plot <- reactive ({
        ggplot(differential, aes(x=logFC, y=minus_log10_Pval)) +
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
        nearPoints(differential, 
                   input$plot_click, 
                   xvar="logFC", 
                   yvar="minus_log10_Pval")
        })
    
    output$clickedPoints <- renderTable({
        clicked()
        }, rownames=T)
    
    #------SERVER for heatmap-------------
    
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
    
    
    
    
    
    #------DOWNLOADS----------------------
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

    # Show all data in Data tab
    output$limma_data <- renderDataTable(
        differential
    )
    output$exp_data <- renderDataTable(
        expression
    )
    output$targets <- renderDataTable(
        experiment
    )
    
} # close server

# Run the application 
shinyApp(ui = ui, server = server)
