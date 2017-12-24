library(shiny)
library(shinydashboard)   
library(RMySQL)
options(stringsAsFactors = F)
createLink <- function(link,char) {
  sprintf('<a href="%s" target="_blank" class="btn btn-primary">%s</a>',link,char)
}

shinyServer(function(input, output, session){
 
  
})