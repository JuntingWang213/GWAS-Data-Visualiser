library(shiny)
library(shinydashboard)
library(plotly)
library(manhattanly)
library(bslib)
library(dplyr)
library(magrittr)
library(forcats)
library(readr)
library(tidyverse)
library(data.table)
library(ggplot2)
library(DT)
library(stringr)

ui <- shinyUI(
dashboardPage(
dashboardHeader(title = "GWAS Data Visualiser", titleWidth = 235),
dashboardSidebar(
width = 235,
sidebarMenu(
# gwas file upload:
fileInput("file1", "Choose csv/txt File", accept = c("text/csv","text/comma-separated-values,text/plain")),
# annotation file upload:
fileInput("file2", "Choose Annotation File", accept = c("text/csv","text/comma-separated-values,text/plain")),
radioButtons('sep', 'Data File Separator Value:', c(Comma=',',Semicolon=';',Tab='\t')),
actionButton(inputId = "DTplot", label = "Merge Data", style='padding:8px; font-size:110%; width:170px'),
# annotation the annotation file:
menuItem("Download Annotation File",     
actionButton("downloadData0", label = a(href="https://www.dropbox.com/scl/fi/amhc1ri1vveh88sei4g5v/annotation.txt?rlkey=onhl168f23w0g5d0ki7zp8gc4&dl=0", 
                                        "Download Annotation File"),
)),   
# column name for manhattan plot:
menuItem("Manhattan Plot",  
uiOutput("manOutput"),
textAreaInput("snp_values", label = "SNP Column Name", value = "",height ='20px',width ='200px'),
textAreaInput("chr_values", label = "Chromosome Column Name", value = "",height ='20px',width ='200px'),
textAreaInput("bp_values", label = "Base-Pair Position Column Name", value = "",height ='20px',width ='200px'),
textAreaInput("p_values", label = "P-Value Column Name", value = "",height ='20px',width ='200px'),
actionButton(inputId = "manhattanplot", label = "Plot",style='padding:8px; font-size:110%; width:170px')),
# 3 filter conditions:               
menuItem("Filtering Conditions",  
# condition 1:
varSelectInput(inputId = "getColumn", label="Filter Condition 1", data = ""),
selectInput(inputId = "level", label="List of choices", choices = NULL),
# condition 2:
varSelectInput(inputId = "getColumn1", label="Filter Condition 2", data = ""),
selectInput(inputId = "level1", label="List of choices", choices = NULL),
# condition 3:
varSelectInput(inputId = "getColumn2", label="Filter Condition 3", data = ""),
selectInput(inputId = "level2", label="List of choices", choices = NULL)),
# filter by range: 
menuItem("Filter by Region",  
textAreaInput("region_values", label = "Condition (chromosome:start bp-end bp/chromosome)", value = "",width ='200px'),
actionButton(inputId = "rangeplot", label = "Plot",style='padding:8px; font-size:110%; width:170px')),

# filter buttons:                 
menuItem("Apply Filtering",  
actionButton(inputId = "filterplot", label = "Filter Manhattan Plot", style='padding:8px; font-size:90%;width:180px'),
actionButton(inputId = "filtertable", label = "Filter data table", style='padding:8px; font-size:90%;width:180px')),
menuItem('Download Filtered File',
actionButton('downloadData', label ='Conditions Filtered File'),
actionButton('downloadRange', label ='Region Filtered File'))
)),
# main panel:                
dashboardBody(
fluidPage(
plotlyOutput(outputId = "mymanhattan"), 
DTOutput('tbl')
)
)
)
)


server <- function(input, output) {
options(shiny.maxRequestSize = 150000*1024^2)


  
output$downloadData <- downloadHandler(filename = function(){paste('data-', Sys.Date(), '.csv', sep='')},
content = function(file) {write.csv(full_join1(), file)})
output$downloadRange <- downloadHandler(filename = function(){paste('data-', Sys.Date(), '.csv', sep='')},
                                       content = function(file) {write.csv(rangefilter(), file)})
  
# read the gwas file： 
rdata <- reactive({
inFile <- input$file1
if(is.null(inFile))
return(NULL) 
data <- fread(inFile$datapath, sep=input$sep)
})
# read the annotation file：   
rdata1 <- reactive({ 
inFile <- input$file2
if(is.null(inFile))
return(NULL) 
fread(inFile$datapath, sep=input$sep)
})
# merge the two files together：
full_join <- reactive({
df <- merge(rdata(), rdata1(), fread, by.x = "RSID", by.y = "Existing_variation", all.x= TRUE)
na.omit(df)
})
#==========================================================================================================
datList0 <- eventReactive(input$manhattanplot, {
full_join() %>% mutate_at(c(input$p_values), as.numeric)

})
datList <- eventReactive(input$manhattanplot, {
datList0() %>% filter(-log10(datList0()[[input$p_values]])>1)
})

#==========================================================================================================
# manhattan plot output：
observeEvent(input$manhattanplot,{
output$mymanhattan <- renderPlotly({
manhattanly(
datList(), 
snp = input$snp_values,
chr = input$chr_values, 
bp = input$bp_values, 
p = input$p_values,
suggestiveline = F, 
genomewideline = F,
col = c("#D2691E","#800080","#6495ED","#9ACD32")
)
})
})

#==========================================================================================================
# data table output：
observeEvent(input$DTplot,{
output$tbl <- DT::renderDataTable(
full_join(),options = list(scrollX = TRUE)
)
})

#==========================================================================================================
# filter parts input:
full_join1 <- reactive({
  filter(full_join(), 
         full_join()[[input$getColumn]] %in% input$level,
         full_join()[[input$getColumn1]] %in% input$level1,
         full_join()[[input$getColumn2]] %in% input$level2
  )
})

# filter the data table:
observeEvent(input$filtertable,{
output$tbl <- renderDT(full_join1(),options = list(scrollX = TRUE))
})

# P column to numeric after filtering:
datList1 <- eventReactive(input$filterplot,{
full_join1() %>% mutate_at(c(input$p_values), as.numeric)
})

#==========================================================================================================
# filter the manhattan plot: 
observeEvent(input$filterplot,{
output$mymanhattan <- renderPlotly({
manhattanly(
datList1(), 
snp = input$snp_values,
chr = input$chr_values, 
bp = input$bp_values, 
p = input$p_values,
suggestiveline = F, 
genomewideline = F,
col=c("#D2691E","#800080","#6495ED","#9ACD32"))
})
})

#==========================================================================================================
# filter the manhattan plot by range: 
observeEvent(input$rangeplot,{
  output$mymanhattan <- renderPlotly({
    manhattanly(
      datList2(), 
      snp = input$snp_values,
      chr = input$chr_values, 
      bp = input$bp_values, 
      p = input$p_values,
      suggestiveline = F, 
      genomewideline = F,
      col=c("#D2691E","#800080","#6495ED","#9ACD32"))
  })
})
  
#==========================================================================================================
# filter condition parts:
observeEvent(input$Filterplot,{
output$tbl <- renderDT(full_join1(),options = list(scrollX = TRUE))
})
#===========================================================================================================================  
# the first filter: 
observeEvent(rdata1(),{
updateVarSelectInput(inputId = "getColumn",data = full_join(),selected = "")
})
  
observeEvent(input$getColumn,{
choice.name <- as.factor(full_join()[[input$getColumn]])
updateSelectInput(inputId = "level", choices =  levels(choice.name),selected = NULL)
})

# the second filter:   
observeEvent(rdata1(),{
updateVarSelectInput(inputId = "getColumn1",data = full_join(),selected = "")
})
  
observeEvent(input$getColumn1, {
choice.name <- as.factor(full_join()[[input$getColumn1]])
updateSelectInput(inputId = "level1", choices =  levels(choice.name),selected = NULL)
})
  

# the third filter:  
observeEvent(rdata1(),{
updateVarSelectInput(inputId = "getColumn2",data = full_join(),selected = "")
})
  
observeEvent(input$getColumn2, {
choice.name <- as.factor(full_join()[[input$getColumn2]])
updateSelectInput(inputId = "level2", choices =  levels(choice.name),selected = NULL)
})

#==========================================================================================================
# the range filter input:
rangefilter <- reactive({
full_join() %>%
      filter(CHR %in% as.character(map(str_split(input$region_values, "[:-]"), 1))) %>% 
      filter(between(POS,as.character(map(str_split(input$region_values, "[chr:-]"), 2)),
                     as.character(map(str_split(input$region_values, "[chr:-]"), 3))
      )
      )
  })
datList2 <- reactive({
rangefilter() %>% mutate_at(c(input$p_values), as.numeric)
})



}
shinyApp(ui = ui, server = server)

