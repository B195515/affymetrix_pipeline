###################### 
# R Shiny Interactive plots
# By B195515 - 2022
# Functional Genomic Technologies, SEM2, 2022
######################
# Volcano plot code modified from:
# https://github.com/FredHutch/interactiveVolcano
# https://2-bitbio.com/2017/12/clickable-volcano-plots-in-shiny.html
######################
# Heatmap plot code modified from class code provided 
# in FGT 2022 by Simon Tomlinson
######################

library(shiny)
library(ggplot2)
library(DT)
library(pheatmap)
library(Cairo)
library(shinyjs)

# Load data
load("differential.Rdata")
load("experiment.Rdata")
load("expression.Rdata")

# Change file/var names here
GSEnumber <- "GSE49448"
pval_col <- "minus_log10_Pval"
fc_col <- "logFC"

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
                                min=-8, max=8, value=4, step=0.1),
                    
                    downloadButton('download_plot', 'Download plot')

                    ), # end sidebarPanel
            
                mainPanel(
                    plotOutput('volcanoPlot',
                               width="100%",
                               height="600px",
                               click='plot_click',
                               brush='plot_brush',
                               hover='plot_hover'
                       ),
                    
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
                                min=6, max=14, value=10),
                    
                    sliderInput("font_col",
                                "Font size col:",
                                min=6, max=14, value=10),

                    # Select annotation
                    selectInput("select", 
                                "Select annotation", 
                                choices=c("none",colnames(experiment)), 
                                selected="Group", multiple=T, selectize=TRUE),
                    
                    checkboxInput("srownames", "Show Row Names", TRUE),
                    
                    checkboxInput("logtansform", "Log transform values", FALSE),
                    
                    radioButtons("norm", "Scale by", 
                                choices=c("none", "row", "column")),
                    
                    downloadButton('download_hm', 'Download heatmap')
                    ), # end sidebarPanel
                
                mainPanel(
                    plotOutput("distPlot", 
                               height="800",
                               width="80%")
                    ) # end mainPanel
                ) # end sidebarlayout
            ), # end panel
        
        tabPanel(
            "Limma data",
            h1("About the data"),
             
            sidebarLayout(
                sidebarPanel(
                    width=3,

                    h4("say something here"),

                    downloadButton('download_limma', 'Download .csv file')
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

                    h4("say something here"),

                    downloadButton('download_exp', 'Download .csv file')
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

                    h4("say something here"),

                    downloadButton('download_targets', 'Download .csv file')
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

    
    #------SERVER for volcano---------
    # Subset by threshold
    is_de <- reactive ({
        abs(differential[[fc_col]]) >= input$fc_thres & 
            differential[[pval_col]] >= input$pval_thres
        })
    
    # Factorize DE points
    de_vec <- reactive({
        factor(is_de(), levels=c("TRUE", "FALSE"))
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
                   xvar=fc_col, 
                   yvar=pval_col)
        })
    
    output$clickedPoints <- renderTable({
        clicked()
        }, rownames=T)
    
    
    
    #------SERVER for heatmap-------------
    # Store plot in fx for multiple calls
    heatplot <- function(){
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
        }
    
    output$distPlot <- renderPlot({
        heatplot()
        }, alt="Heatmap of top 50 expressed genes", execOnResize = F)
    
    
    
    #------SHOW ALL DATA-------------------
    output$limma_data <- renderDataTable(differential)
    output$exp_data <- renderDataTable(expression)
    output$targets <- renderDataTable(experiment)
    
    
    
    #------DOWNLOADS----------------------
    # output$download_volcano <- downloadHandler(
    #     filename = function(){
    #         paste0("volcano-plot-",GSEnumber,".png")},
    #     content = function(file){
    #         ggsave(file=reactive_plot(), filename=file)}
    #     )
    output$download_plot <- downloadHandler(
        filename <- function() {
            paste("Volcano-plot_",GSEnumber,".pdf", sep="")},
        content <- function(file) {
            pdf(file, width=9, height=9)
            plot(reactive_plot())
            dev.off()},
        contentType = "application/pdf"
        )
    
    
    # TODO: debug
    # output$download_hm <- downloadHandler(
    #     filename <- function() {
    #         paste("HeatmapDEGs_",GSEnumber,".pdf", sep="")},
    #     content <- function(file) {
    #         pdf(file, width=9, height=9)
    #         plot(heatplot())
    #         dev.off()},
    #     contentType = "application/pdf"
    # )
    
    
    output$download_limma <- downloadHandler(
        filename = function(){
            paste0("limmaDE_",GSEnumber,".csv")},
        content = function(file){
            write.csv(differential, file)}
        )
    
    
    output$download_exp <- downloadHandler(
        filename = function(){
            paste0("expressionTop50_",GSEnumber,".csv")},
        content = function(file){
            write.csv(expression, file)}
        )
    
    
    output$download_targets <- downloadHandler(
        filename = function(){
            paste0("targets_",GSEnumber,".csv")},
        content = function(file){
            write.csv(experiment, file)}
        )
        
} # close server

# Run the application 
shinyApp(ui = ui, server = server)
