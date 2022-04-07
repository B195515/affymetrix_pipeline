###################### 
# R Shiny Interactive plots
# By B195515 - 2022
# Functional Genomic Technologies, SEM2, 2022
######################
# Volcano plot code modified from:
# https://github.com/FredHutch/interactiveVolcano
# https://2-bitbio.com/2017/12/clickable-volcano-plots-in-shiny.html
# Hovered points from: https://github.com/JoachimGoedhart/VolcaNoseR
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
gene_col <- "Symbol"
fdr_col <- "adj.P.Val"

# UI------------------
ui <- shinyUI(fluidPage(
    tabsetPanel(
        
        tabPanel(
            "Volcano plot",
            h1("GSE49448: Limma DEGs"),
            
            sidebarLayout(
                sidebarPanel(
                    width=3,
                    h3("Volcano Plot"),
                    "Differentially-expressed genes as output from \"lmfit\"
                    function.",
                    br(),br(),
                    "Hover over a point to show gene symbol, log2 fold change value,
                    and FDR-adjusted p-value.",
                    br(),br(),
                    "Click anywhere on the plot to show a list of genes in the
                    proximity.",
                    
                    h3("Thresholds"),
                    sliderInput("pval_thres",
                                "Set significance threshold",
                                min=0, max=4, value=2.5, step=0.1),
                    textOutput("selected_pval"),
                    br(),
                    sliderInput("fc_thres",
                                "Set FC threshold",
                                min=-8, max=8, value=4, step=0.1),
                    
                    h3("Downloads"),
                    downloadButton('download_plot', 'Download plot')),
                
                mainPanel(
                    plotOutput('volcanoPlot',
                               width="100%",
                               height="600px",
                               click='plot_click',
                               brush='plot_brush',
                               hover=hoverOpts('plot_hover',
                                               delay=10, delayType="debounce")),
                    uiOutput('hover_info'),
                    tableOutput('clickedPoints')) # clicked points table)
                )
            ),
        
        tabPanel(
            "Heatmap",
            h1("GSE49448: Expression heatmap of Top 50 DEGs"),
            
            sidebarLayout(
                sidebarPanel(
                    width=3,
                    
                    h3("Heatmap"),
                    "Heatmap of RMA-normalized expression data for all samples,
                    for the top 50 DEGs as ranked by Limma.",
                    
                    h3("Settings"),
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
                    
                    h3("Downloads"),
                    downloadButton('download_hm', 'Download heatmap')),
                
                mainPanel(
                    plotOutput("distPlot", height="800", width="80%"))
            )
        ),
        
        tabPanel(
            "Limma data",
            h1("Limma differentially expressed genes (ranked)"),
            
            sidebarLayout(
                sidebarPanel(
                    width=3,
                    
                    h3("Data legend"),
                    "Data table shows the gene symbol, log2 fold change, 
                    p-value, FDR-adjusted p-value, B score (B-statistics), 
                    and -log10(adj.p.value).",
                    br(),br(),
                    "Rownames are the probe ids from the 
                    MOE430v2 microarray chips.",
                    
                    h3("External links"),
                    uiOutput("links1"),
                    
                    h3("Downloads"),
                    downloadButton('download_limma', 'Download .csv file')),
                
                mainPanel(
                    dataTableOutput("limma_data"))
            )
        ),
        
        tabPanel(
            "Expression data",
            h1("About the data"),
            
            sidebarLayout(
                sidebarPanel(
                    width=3,
                    h3("Data legend"),
                    "Data table shows the RMA-normalized expression of the
                    top 50 DEGs output from Limma model, for all array samples.",
                    br(),br(),
                    "Rownames show gene symbol and 
                    corresponding chip array probe ID,",
                    h3("Downloads"),
                    downloadButton('download_exp', 'Download .csv file')),
                
                mainPanel(
                    dataTableOutput("exp_data"))
            )
        ),
        
        tabPanel(
            "Experiment",
            h1("Experiment details"),
            sidebarLayout(
                sidebarPanel(
                    width=3,
                    h3("About data",GSEnumber),
                    "Mouse single cardiomyocytes (Ctl) and their derived 
                    single CPCs were captured using microfluidic deviced and 
                    cDNA synthesized and amplified and labelled using 
                    NuGENE kits, and regular Affymetrix hybridization and 
                    wash protocols were used to process the 
                    mouse whole-genome array GeneChip 430 2.0.",
                    h4("Citation"),
                    "PMID: 31231540",
                    h3("Downloads"),
                    downloadButton('download_targets', 'Download .csv file')),
                mainPanel(
                    h3("Transcriptomic Signatures of Mouse Cardiomyocyte 
                       Dedifferentiation In vitro (Affymetrix)"),
                    br(),br(),
                    dataTableOutput("targets"))
            )
        )
    ))
) # end UI


# Server--------------
server <- function(input, output, session){
    #------SERVER for volcano---------
    # Show actual p-value
    output$selected_pval <- renderText({
        paste("adj.P.Val =", format(10**(-input$pval_thres), digits=5))
    })
    
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
        ggplot(differential, aes(x=.data[[fc_col]], y=.data[[pval_col]])
               ) +
            geom_point(alpha=.6, aes(color=de_vec())) +
            coord_cartesian() + 
            geom_hline(yintercept=input$pval_thres, 
                       linetype="dashed", col="grey", size=1) +
            geom_vline(xintercept=c(input$fc_thres, -input$fc_thres), 
                       linetype="dashed", col="grey", size=1) + 
            xlab("log2 fold change") +
            ylab("-log10(FDR)") +
            labs(color="Differentially expressed") +
            theme(legend.position="bottom")
    })
    
    # Show volcano plot    
    output$volcanoPlot <- renderPlot({
        reactive_plot()
    })
    
    # Retrieve clicked points
    clicked <- reactive({
        # define the x and y variables:
        nearPoints(differential, input$plot_click, 
                   xvar=fc_col, yvar=pval_col)
    })
    
    output$clickedPoints <- renderTable({
        clicked()
    }, rownames=T)
    
    # Show hovered point
    hovered <- reactive({
        hover <- input$plot_hover
        point <- nearPoints(differential, hover, 
                            xvar=fc_col, yvar=pval_col,
                            threshold=10, maxpoints=1, addDist=FALSE)
        if (nrow(point)==0) return(NULL)
        
        left_pct <- (hover$x-hover$domain$left) / (hover$domain$right-hover$domain$left)
        top_pct <- (hover$domain$top-hover$y) / (hover$domain$top-hover$domain$bottom)
        
        left_px <- hover$range$left+left_pct * (hover$range$right-hover$range$left)
        top_px <- hover$range$top+top_pct * (hover$range$bottom-hover$range$top)
        
        style <- paste0("position:absolute; padding: 5px; z-index:100; 
                        background-color: rgba(200,200,245,0.65); ",
                        "left:", left_px+20, "px; top:", top_px+32, "px;")
        
        wellPanel(style=style,
                  p(HTML(paste0("<b>Name: </b>", point[gene_col], "<br/>",
                                "<b>FC: </b>", round(point[fc_col],2), "<br/>",
                                "<b>Sig: </b>", round(point[fdr_col],5), "<br/>",
                                NULL))
                    )
                )
    })
    
    output$hover_info <- renderUI({
        hovered()
    })
    
    #------SERVER for heatmap-------------
    # Store plot in fx for multiple calls
    heatplot <- function(){
        if(input$logtansform){
            expression <- log2(expression + 1)}
        
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
        heatplot()}, 
        alt="Heatmap of top 50 expressed genes", 
        execOnResize = F)
    
    #------SHOW ALL DATA-------------------
    bstats <- a("More about B-statistics", 
                href="https://support.bioconductor.org/p/6124/",
                target="_blank")
    chip_ds <- a("Other MOE430v2 datasets", 
                 href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL1261",
                 target="_blank")
    chip_support <- a("Info and chip support", 
                      href="https://www.thermofisher.com/order/catalog/product/900497",
                      target="_blank")
    output$links1 <- renderUI({
        tagList(bstats, br(), br(), 
                chip_ds, br(), br(),
                chip_support)})
    
    output$limma_data <- renderDataTable(differential)
    output$exp_data <- renderDataTable(expression)
    output$targets <- renderDataTable(experiment)
    
    
    #------DOWNLOADS----------------------
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
