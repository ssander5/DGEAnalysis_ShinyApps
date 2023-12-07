library(shiny)
library(variancePartition)
library(sva)

##########################
#         UI       
##########################
ui <- fluidPage(
  titlePanel("Variance Partitioning"),
  
  sidebarLayout(
  
    sidebarPanel(
      fileInput("normCounts", "Normalized Counts (csv)", accept=".csv"),
      fileInput("metaDataTable", "Meta Data (csv)", accept=".csv"),
      actionButton("load", "Load Files"),
    ),
  
    
    mainPanel(
      tabsetPanel(
        tabPanel("Input",
          h1("Welcome"),
          br(),
          p("This app is to help investigate the different factors impacting variance of expression in your data.  We make use of the variancePartition library in R, which in turn uses a mixed model to partition different the variation seen in expression into different sources."),
          p("You will need to upload a normalized counts matrix - which you can get from the edgeR tool.  You will also need a metadata table to define sources of variance.  Please see the example below."),
          em("Example Counts Tables:"),
          br(),
          tableOutput("exampleCountsOne"),
          tableOutput("exampleCountsTwo"),
          br(),
          em("Example MetaData Tables:"),
          tableOutput("exampleMDOne"),
          tableOutput("exampleMDTwo")
        ),
        
        tabPanel("Exploration",
          conditionalPanel(
            condition="output.inputsUploaded",
            tags$h1("Initial Data"),
            p("TBD"),
            plotOutput("correlation"),
            p("TBD"),
            plotOutput("violin"),
            p("TBD"),
            plotOutput("barSplit"),
            p("TBD"),
            selectInput("selectPick", "Select a gene", as.character(rownames(geneExpr))),
            selectInput("selectVect", "Select a source of variation", colnames(info)),
            p("TBD"),
            plotOutput("pickSplit"),
            p("TBD"),
            plotOutput("pickStrat")
          )
        ),
        tabPanel("Correction",
            tableOutput("batchCorrected")
        )

      )
    )
  )
)

##########################
#         SERVER        
##########################

server <- function(input, output) {
  
########################################
#########Example Table Render###########
########################################
    output$exampleCountsOne <- renderTable({
    # create example counts table
    exCountsTable <- data.frame(
      Gene = c("gene-1", "gene-2", "gene-3", "gene-4", "gene-5"),
      SampleOne = c(0, 0, 0, 0, 0),
      SampleTwo = c(10, 20, 30, 40, 50),
      SampleThree = c(111, 222, 333, 444, 555),
      SampleFour = c(1, 2, 3, 4, 5),
      SampleFive = c(0, 0, 0, 0, 0),
      SampleSix = c(1000, 2000, 3000, 4000, 5000),
      SampleSeven = c(11, 12, 13, 14, 15),
      SampleEight = c(0, 0, 0, 0, 0)
    )
  })
  
  # render example gene counts table
  output$exampleCountsTwo <- renderTable({
    # create example counts table
    exCountsTable <- data.frame(
      Gene = c("geneA", "geneB", "geneC"),
      sample_1 = c(0, 0, 0),
      sample_2 = c(10, 20, 30),
      sample_3 = c(111, 222, 333),
      sample_4 = c(1, 2, 3),
      sample_5 = c(3, 3, 3),
      sample_6 = c(1000, 2000, 3000),
      sample_7 = c(11, 12, 13),
      sample_8 = c(1, 1, 1),
      sample_9 = c(123, 12, 1),
      sample_10 = c(3, 32, 321),
      sample_11 = c(33, 333, 33),
      sample_12 = c(2, 2, 2)
    )
  })
  
  # render first example gene counts table
  output$exampleMDOne <- renderTable({
    # create example counts table
    exDesignTable <- data.frame(
      Sample = c("SampleOne", "SampleTwo", "SampleThree", "SampleFour", "SampleFive", "SampleSix"),
      Group = c("cntrl", "cntrl", "cntrl", "treat", "treat", "treat"),
      Batch = c(1,1,1,2,2,2)
    )
  })
  
  # render second example gene counts table
  output$exampleMDTwo <- renderTable({
    # create example counts table
    exDesignTable <- data.frame(
      Individual = c("sample_1", "sample_2", "sample_3", "sample_4", "sample_5", "sample_6", "sample_7", "sample_8", "sample_9", "sample_10", "sample_11", "sample_12"),
      Factors = c("cntrl.high", "cntrl.high", "cntrl.high", "cntrl.low", "cntrl.low", "cntrl.low", "treat.high", "treat.high", "treat.high", "treat.low", "treat.low", "treat.low"),
      Genotype = c("WT","WT","WT","WT","WT","WT","Sc","Sc","Sc","Sc","Sc","Sc"),
      Batch = c(1,1,1,1,2,2,2,3,3,3,4,4)
    )
  })
  
  ########################################
  #########    Input Handling ###########
  ########################################
  
  inputGeneCounts = reactive({
    req(input$normCounts)
    if(is.null(input$normCounts)){return(NULL)}
  
    countInput = read.csv(file = input$normCounts$datapath, row.names=1)
    countNamesSet = rownames(countInput)
  
  })
  
  inputMDTable = reactive({
    req(input$metaDataTable)
    if(is.null(input$metaDataTable)){return(NULL)}
    
    MDInput = read.csv(file = input$metaDataTable$datapath, row.names=1)
    mdNamesSet = rownames(MDInput)
    
  })

  ########################################
  #########     Input Check    ###########
  # Both files uploaded
  ########################################
  
  output$inputsUploaded <- function(){
    # check if the input files are valid
    if(is.null(inputGeneCounts())) {
      return(FALSE)
    }else if(is.null(inputMDTable())) {
      return(FALSE)
    }
    return(TRUE)
  }
  
  outputOptions(output, 'inputsUploaded', suspendWhenHidden=FALSE)
  
  ########################################
  #########  Run Partitioning  ###########
  ########################################
  
  output$correlation = renderPlot({
    data(varPartData)
    form = ~Age +(1| Individual) + (1|Tissue) +(1|Batch)
    C = canCorPairs(form, info)
    plotCorrMatrix(C)
  }) 
  
  output$violin = renderPlot({
    data(varPartData)
    form = ~Age +(1| Individual) + (1|Tissue) +(1|Batch)
    varPart = fitExtractVarPartModel(geneExpr,form,info)
    vp = sortCols(varPart)
    plotVarPart(vp)
  })   
  
  output$barSplit = renderPlot({
    data(varPartData)
    form = ~Age +(1| Individual) + (1|Tissue) +(1|Batch)
    varPart = fitExtractVarPartModel(geneExpr,form,info)
    vp = sortCols(varPart)
    plotPercentBars(vp[1:10,])
  })
  
  output$pickSplit = renderPlot({
    i=input$selectPick
    data(varPartData)
    form = ~Age +(1| Individual) + (1|Tissue) +(1|Batch)
    varPart = fitExtractVarPartModel(geneExpr,form,info)
    vp = sortCols(varPart)
    options(repr.plot.height = 2)
    plotPercentBars(vp[i,], width = 0.1)
  }, height = 400)

  output$pickStrat = renderPlot({
    i=input$selectPick
    data(varPartData)
    form = ~Age +(1| Individual) + (1|Tissue) +(1|Batch)
    varPart = fitExtractVarPartModel(geneExpr,form,info)
    x=info[[input$selectVect]]
    GE = data.frame(Expression = geneExpr[i, ], Vect=x)
    plotStratify(Expression~Vect, GE, main=rownames(geneExpr)[i], text = format(varPart$Individual[i]*100)) + 
      labs(x=input$selectVect)
  })
  
  ########################################
  ##########  Remove Vectors  ############
  ########################################
  output$batchCorrected = renderTable({
    count_matrix <- matrix(rnbinom(400, size=10, prob=0.1), nrow=50, ncol=8)
    batch <- c(rep(1, 4), rep(2, 4))
    adjusted = ComBat_seq(count_matrix, batch=batch, group=NULL)
  })
  
  
  ########################################
  #######  ReRun Partitioning  ###########
  ########################################
  #tbd
  
  ############End of Server###########
} 

##########################
#         Main     
##########################
shinyApp(ui = ui, server = server)
