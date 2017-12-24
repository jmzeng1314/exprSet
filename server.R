library(shiny)
library(shinydashboard)   
library(RMySQL)
options(stringsAsFactors = F)
createLink <- function(link,char) {
  sprintf('<a href="%s" target="_blank" class="btn btn-primary">%s</a>',link,char)
}

log_cat <- function(info='hello world~',file='log.txt'){
  cat(as.character(Sys.time()),info ,"\n",file=file,append=TRUE)
}

shinyServer(function(input, output, session){
  IP <- reactive({ 
    tmp<-input$getIP 
    as.character(tmp$ip)
  })
  
  ## pre-defined global variables 
  gVs <- reactiveValues(
    userID = NULL,
    projectID    = NULL, 
    exprSet = NULL,
    geneInfo = NULL,
    groupInfo= NULL
  )
  
  dataLoad <- reactive({  ## Data input.
    
    if(input$dataInput==1){  ## Load example data. 
      ##   list("airway"=1,'sCLLex'=2), 
      if(input$sampleData_exprSet == 1){
        
        isolate({ 
          withProgress(message = 'Data upload in progress...',  
                       detail = 'This may take a while...', value = 1,{
            load("www/data/airway.Rdata")
            
          }) 
        }) ## end for isolate 
        
        
        
      } else  if(input$sampleData_exprSet == 2){
        
        isolate({
          
          withProgress(message = 'Data upload in progress...',  
                       detail = 'This may take a while...', value = 1,{
                         load("www/data/sCLLex.Rdata")
                  
                       }) 
        })## end for isolate 
        
      } ## end for elseIf
      
      
      
    }
    else if(input$dataInput==2){  ## Upload data.
      ## First for expression matrix 
      inFile <- input$upload_exprSet
      mySep <- switch(input$fileSepDF_exprSet, '1'=",",'2'="\t",'3'=";", '4'="")
      if (is.null(input$upload_exprSet))  {return(NULL)}
      isolate({
        
        withProgress(message = 'Data upload in progress...', 
                     detail = 'This may take a while...', value = 1,{
          
          # if (file.info(inFile$datapath)$size <= 10485800){
          exprSet <- read.table(inFile$datapath, sep=mySep, header=TRUE, fill=TRUE, na.strings = c("", "NA","."))
          # }
          
          # else print("File is bigger than 10MB and will not be uploaded.")
        })
        
      })
      ## Then for  genes annotation
      inFile <- input$upload_geneInfo
      mySep <- switch(input$fileSepDF_geneInfo, '1'=",",'2'="\t",'3'=";", '4'="")
      if (is.null(input$upload_geneInfo))  {return(NULL)}
      isolate({
        
        withProgress(message = 'Data upload in progress...', 
                     detail = 'This may take a while...', value = 1,{
                       
                       # if (file.info(inFile$datapath)$size <= 10485800){
                       geneInfo <- read.table(inFile$datapath, sep=mySep, header=TRUE, fill=TRUE, na.strings = c("", "NA","."))
                       # }
                       
                       # else print("File is bigger than 10MB and will not be uploaded.")
                     })
        
      })
      ## lastly for group information 
      inFile <- input$upload_groupInfo
      mySep <- switch(input$fileSepDF_groupInfo, '1'=",",'2'="\t",'3'=";", '4'="")
      if (is.null(input$upload_groupInfo))  {return(NULL)}
      isolate({
        
        withProgress(message = 'Data upload in progress...', 
                     detail = 'This may take a while...', value = 1,{
                       
                       # if (file.info(inFile$datapath)$size <= 10485800){
                       groupInfo <- read.table(inFile$datapath, sep=mySep, header=TRUE, fill=TRUE, na.strings = c("", "NA","."))
                       # }
                       colnames(groupInfo)=c('sampleID','group')
                       # else print("File is bigger than 10MB and will not be uploaded.")
                     })
        
      })
    }
     
    exprSet=as.matrix(exprSet)
    return(list(groupInfo=groupInfo,geneInfo=geneInfo,exprSet=exprSet))
    
  })
  
  output$exprSet <- DT::renderDataTable({ 
    tmp=dataLoad() 
    log_cat(paste0('data', Sys.Date(),' ',IP()))
    gVs$groupInfo=tmp$groupInfo
    gVs$exprSet=tmp$exprSet
    gVs$geneInfo=tmp$geneInfo
    exprSet=gVs$exprSet
    DT::datatable(exprSet,
                  extensions = 'FixedColumns',
                  options = list(
                    #dom = 't',
                    scrollX = TRUE,
                    fixedColumns = TRUE
                  )
    )## end for datatable
    
  })  
  
  output$groupInfo <- DT::renderDataTable({
    gVs$groupInfo
  }) 
  output$geneInfo <- DT::renderDataTable({
    gVs$geneInfo
  })  

  
  output$p1_boxplot <- renderPlot({
    exprSet=gVs$exprSet
    if (is.null(exprSet)) return(NULL)
    group_list=gVs$groupInfo$group
    library(reshape2)
    exprSet_L=melt(exprSet)
    colnames(exprSet_L)=c('probe','sample','value')
    exprSet_L$group=rep(group_list,each=nrow(exprSet))
    #head(exprSet_L)
    p=ggplot(exprSet_L,aes(x=sample,y=value,fill=group))+geom_boxplot() 
    p=p+stat_summary(fun.y="mean",geom="point",shape=23,size=3,fill="red")
    p=p+theme_set(theme_set(theme_bw(base_size=20)))
    p=p+theme(text=element_text(face='bold'),axis.text.x=element_text(angle=30,hjust=1),axis.title=element_blank())
    print(p)
    
  })  
  output$p1_histogram <- renderPlot({
    
    exprSet=gVs$exprSet
    if (is.null(exprSet)) return(NULL)
    group_list=gVs$groupInfo$group
    library(reshape2)
    exprSet_L=melt(exprSet)
    colnames(exprSet_L)=c('probe','sample','value')
    exprSet_L$group=rep(group_list,each=nrow(exprSet))
    #head(exprSet_L)
    p=ggplot(exprSet_L,aes(value,fill=group))+geom_histogram(bins = 200)+facet_wrap(~sample, nrow = 4)
    print(p)
    
  })  
  output$p1_density <- renderPlot({
    
    exprSet=gVs$exprSet
    if (is.null(exprSet)) return(NULL)
    group_list=gVs$groupInfo$group
    library(reshape2)
    exprSet_L=melt(exprSet)
    colnames(exprSet_L)=c('probe','sample','value')
    exprSet_L$group=rep(group_list,each=nrow(exprSet))
    #head(exprSet_L)
    p=ggplot(exprSet_L,aes(value,col=group))+geom_density()+facet_wrap(~sample, nrow = 4)
    print(p)
    
  })  
  # output$p1_gpairs <- renderPlot({
  #   exprSet=gVs$exprSet
  #   if (is.null(exprSet)) return(NULL)
  #   library(gpairs)
  #   gpairs(exprSet
  #          #,upper.pars = list(scatter = 'stats') 
  #          #,lower.pars = list(scatter = 'corrgram')
  #   )
  #   
  # })  
  output$p1_cluster <- renderPlot({
    exprSet=gVs$exprSet
    if (is.null(exprSet)) return(NULL)
  
    out.dist=dist(t(exprSet),method='euclidean')
    out.hclust=hclust(out.dist,method='complete')
    plot(out.hclust)
    
  })  
  output$p1_PCA <- renderPlot({
    exprSet=gVs$exprSet
    if (is.null(exprSet)) return(NULL)
    pc <- prcomp(t(exprSet),scale=TRUE)
    pcx=data.frame(pc$x)
    group_list=gVs$groupInfo$group
    pcr=cbind(samples=rownames(pcx),group_list, pcx) 
    p=ggplot(pcr, aes(PC1, PC2))+geom_point(size=5, aes(color=group_list)) +
      geom_text(aes(label=samples),hjust=-0.1, vjust=-0.3)
    print(p)
    
  })  
  output$p1_heatmap <- renderPlot({
    exprSet=gVs$exprSet
    if (is.null(exprSet)) return(NULL)
    
    choose_gene=names(sort(apply(exprSet, 1, mad),decreasing = T)[1:50])
    choose_matrix=exprSet[choose_gene,]
    choose_matrix=scale(choose_matrix)
    heatmap(choose_matrix)
    
  })  
})