library(shiny)
library(shinydashboard) 
suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(edgeR))
suppressMessages(library(DESeq2))
suppressMessages(library(limma))

options(stringsAsFactors = F) 
#options(encoding = 'UTF-8')
page_Home <- fluidPage(
  
)
page_INPUT <- fluidPage(
  ## first row: two boxs 
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
      menuItem("绘图大全简介",tabName = "Home",icon = icon("home")),
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



