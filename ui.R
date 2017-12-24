library(shiny)
library(shinydashboard) 
library(shinyBS)
library(GGally)
library(ggplot2)
#library(shinyAce)
library(knitr)
library(rmarkdown)
library(RCurl)
library(shinyjs)
library(DT)

suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(edgeR))
suppressMessages(library(DESeq2))
suppressMessages(library(limma))

options(stringsAsFactors = F) 
#options(encoding = 'UTF-8')
page_Home <- fluidPage(
  includeMarkdown("README.md") 
)
 
page_INPUT <- fluidPage(
  ################################# Data Upload (start) #####################################
  ### three box in a row, for expression matrix / genes annotation /  group information 
  ################################# Data Upload (start) #####################################
  
  fluidRow( 
    box(title = "upload files", 
        width = 12, solidHeader = TRUE, status = "primary",
      
        
      radioButtons("dataInput", "", 
                   list("Upload 3 files" = 2, "Load example data" = 1), selected=2),
      
      conditionalPanel(condition="input.dataInput=='1'",
                       h5("Example dataset:"),
                       radioButtons("sampleData_exprSet", "", 
                                    list("airway(RNA-seq)"=1,'sCLLex(microarray)'=2), 
                                    selected=2)
      ),
      
      conditionalPanel(condition="input.dataInput=='2'",
                       h5("Upload a delimited text file for expression matrix : "),
                       
                       fileInput("upload_exprSet", "", multiple = FALSE),
                       
                       radioButtons("fileSepDF_exprSet", "Delimiter:", 
                                    list("Comma"=1,"Tab"=2,"Semicolon"=3,"Space"=4),
                                    selected=2,inline = T),
                       
                       h5("Upload a delimited text file for  group information : "),
                       fileInput("upload_groupInfo", "", multiple = FALSE),
                       
                       radioButtons("fileSepDF_groupInfo", "Delimiter:", 
                                    list("Comma"=1,"Tab"=2,"Semicolon"=3,"Space"=4),
                                    selected=2,inline = T),
                       
                       h5("Upload a delimited text file for genes annotation : "),
                       fileInput("upload_geneInfo", "", multiple = FALSE),
                       
                       radioButtons("fileSepDF_geneInfo", "Delimiter:", 
                                    list("Comma"=1,"Tab"=2,"Semicolon"=3,"Space"=4),
                                    selected=2,inline = T),
                       
                       HTML('<p>You can upload your data as separated by comma, tab, semicolon or space.</p>'),
                       HTML('<p>Note: First row must be header.</p>')
      ),
      
      br()
    ) ## end for box 
     
    ),  ### end for fluidRow 
  
  ################################# Data Upload (end) #####################################
  
  
  
  ################################# Data tables #####################################
  fluidRow(
    tabBox(
      title="Your input data:",width = 12,
      tabPanel('exprSet',DT::dataTableOutput("exprSet")),
      tabPanel('groupInfo', DT::dataTableOutput("groupInfo")),
      tabPanel('geneInfo',DT::dataTableOutput("geneInfo"))
    )

  ),
  ################################# Data visualization #####################################
  fluidRow(
    tabBox(
      title="Quality control",width = 12,
      tabPanel('boxplot',plotOutput("p1_boxplot")),
      tabPanel('histogram',plotOutput("p1_histogram")),
      tabPanel('density',plotOutput("p1_density")),
       
      tabPanel('gpairs',h4("time consuming,I don't want to draw this figure for you!!!")),
      tabPanel('cluster',plotOutput("p1_cluster")),
      
      tabPanel('PCA',plotOutput("p1_PCA")),
      tabPanel('heatmap',plotOutput("p1_heatmap"))
    )
    
  ),
  
  fluidRow(
    box(
      title = "Box title", width = 6, status = "primary"
    ),
    box(
      status = "warning", width = 6,
      h6("Box content")
    )
  ),
  ### 
  fluidRow(
    column(width = 4,
           box(
             title = "Title 1", width = NULL, solidHeader = TRUE, status = "primary",
             h6("Box content")
           ),
           box(
             width = NULL, background = "black",
             "A box with a solid black background"
           )
    ),
    
    column(width = 4,
           box(
             title = "Title 3", width = NULL, solidHeader = TRUE, status = "warning",
             h6("Box content")
           ),
           box(
             title = "Title 5", width = NULL, background = "light-blue",
             "A box with a solid light-blue background"
           )
    ),
    
    column(width = 4,
           box(
             title = "Title 2", width = NULL, solidHeader = TRUE,
             h6("Box content")
           ),
           box(
             title = "Title 6", width = NULL, background = "maroon",
             "A box with a solid maroon background"
           )
    )
  )
  
)
page_DEG <- fluidPage(
  
)
page_ENRICH <- fluidPage(
  
)
page_GSEA <- fluidPage(
  
) 
page_About <- fluidPage(
  
)

header=dashboardHeader(
  title =p("表达矩阵分析大全！"
           ,style="font-size:90%;font-style:oblique"
  )
)
sidebar = dashboardSidebar(
  conditionalPanel(
    condition = "1",
    sidebarMenu(
      id = "tabs",
      hr(),
      menuItem("简介",tabName = "Home",icon = icon("home")),
      menuItem("step1:上传输入文件",   tabName = "INPUT",icon = icon("flask") 
               ),
      menuItem("step2:自定义差异分析",tabName = "DEG",icon = icon("flask")),
      menuItem("step3:GO/KEGG富集分析",tabName = "ENRICH",icon = icon("flask")),
      menuItem("step4:GSEA",tabName = "GSEA",icon = icon("flask")),
      menuItem("关于我们", tabName = "About", icon = icon("info-circle"))
    ) ## end for sidebarMenu
  ) ## end for conditionalPanel
) ## end for dashboardSidebar

body=dashboardBody(
  tabItems(
    tabItem(tabName = "Home",page_Home),
    tabItem(tabName = "INPUT",page_INPUT), 
    tabItem(tabName = "DEG",page_DEG), 
    tabItem(tabName = "ENRICH",page_ENRICH), 
    tabItem(tabName = "GSEA",page_GSEA), 
    tabItem(tabName = "About",page_About)
  ) 
)

shinyUI(
  dashboardPage(
    header,
    sidebar,
    body,
    title = '表达矩阵分析大全'
  )
)



