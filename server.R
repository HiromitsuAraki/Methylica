# Load packages ----
library(shiny)
library(shinyBS)
library(shinyjs)
library(shinydashboard)
library(R.utils)
library(data.table)
library(jpeg)
library(png)
library(fields)
library(circlize)
library(ComplexHeatmap)
library(openxlsx)
library(XLConnect)
#library(rtable)
library(rmarkdown)
library(DT)
library(org.Hs.eg.db)
library(ggpubr)
library(sparklyr)
library(flexdashboard)
library(shinycssloaders)
library(highcharter)

shinyInput <- function(FUN, len, id, ...) {
  inputs <- character(len)
  for (i in seq_len(len)) {
    inputs[i] <- as.character(FUN(paste0(id, i), ...))}
  inputs
}

# Source ICA pipeline ----  
#source("ICA_pipeline_shiny_gan2018.R")
source("ICA_pipeline_shiny.R")


###Server
shinyServer(function(input, output, session) {
  
  #acceptable big data size
  options(shiny.maxRequestSize=3000*1024^2) 
  
  
  #data import
  data_methylome <- reactive({
    inFile <- input$data_input
    if (is.null(inFile))
      return(NULL)
      #return(data.frame(fread(inFile$datapath,header=T),row.names=1))
      return(data.frame(fread(inFile$datapath,header=T)))
  })
  
  
  #sample import
  sampleInfo <- reactive({
    inFile <- input$sample_input
    if (is.null(inFile))
      return(NULL)
    return(read.table(inFile$datapath,header=T,row.names=1))
  })
  
  
  #selection of sample property in IC boxplot
  observe({
    updateSelectInput(session, "SampleProperty",
                      label =  "Sample property",
                      choices = colnames(sampleInfo())
    )
  })
  
  #selection of sample property in methylation boxplot of each entry
  observe({
    updateSelectInput(session, "SampleProperty2",
                      label = "Sample property",
                      choices = colnames(sampleInfo())                     
    )
  })
  
  
  output$variables = renderUI({
    colnames(sampleInfo())
  })
  
  
  ica=reactiveValues()
  observeEvent(input$goButton1,{
    Data_methylome <- data_methylome()
    withProgress(message = 'Calculation in progress', {
      if (input$ICnumPCA==FALSE){
        ica$ICAres <- run_ICA(Data_methylome,"fastICA",as.integer(input$ICnum))
      }
      if (input$ICnumPCA==TRUE){
        ica$ICAres <- run_ICA_PC80(Data_methylome,"fastICA")
      }
    })
  })
  
  
  ica_heatmap = reactive({
    ICA_result=ica$ICAres
    GenomicFeature=input$GenomicFeature
    inFile_ICA=sampleInfo()
    make_ICA_heatmap_shiny(ICA_result,inFile_ICA,GenomicFeature)
  })
  
  
  #reset_plot = reactive({
  #  Reset_Plot()
  #})
  
　#output$ReSet <- renderPlot({
  #  reset_plot()
  #})
  
  output$ICA_HC <- renderPlot({
    updateTabItems(session = session, inputId = "ics", selected = "IC")
    ica_heatmap()
  })
  
  
  output$ICboxplot <- renderPlot({
    updateTabItems(session = session, inputId = "ics", selected = "IC") 
    ICA_result=ica$ICAres
    inFile_ICA=sampleInfo()
    #HEight=length(inFile_ICA)
    TargetGroup=as.character(input$SampleProperty)
    ggplot2boxplot_IC_shiny(ICA_result,inFile_ICA,TargetGroup)
  })
  
  
  output$downloadICA_HC <- downloadHandler(
    filename="ICA_Heatmap.pdf",
    content=function(file){
      pdf(file, height= 8, width=12)
      print(ica_heatmap())
      dev.off()
    },
    contentType="image/pdf"
  )
  
  
  output$downloadICAboxplot <- downloadHandler(
    filename="ICA_Boxplot.pdf",
    content=function(file){
      pdf(file, height= 10, width=15)
      print(ggplot2boxplot_IC_shiny(ica$ICAres,sampleInfo(),as.character(input$SampleProperty)))
      dev.off()
    },
    contentType="image/pdf"
  )
  
  
  
  output$downloadICAboxplot <- downloadHandler(
    filename="ICA_Boxplot.pdf",
    content=function(file){
      pdf(file, height= 10, width=15)
      print(ggplot2boxplot_IC_shiny(ica$ICAres,sampleInfo(),as.character(input$SampleProperty)))
      dev.off()
    },
    contentType="image/pdf"
  )
  
  
  loading_heatmap = reactive({
    ICA_result=ica$ICAres
    inFile_ICA=sampleInfo()
    Target=as.integer(input$ic)
    Zscore=as.integer(input$loadings)
    GenomicFeature=input$GenomicFeature
    datamat=data_methylome()
    make_highLF_heatmap_target_shiny(ICA_result,inFile_ICA,datamat,Zscore,Target,GenomicFeature)
  })
  
  
  output$loadingheatmap <- renderPlot({
    #updateTabsetPanel(session = session, inputId = "inTabset", selected = "Loadings") 
    updateTabItems(session = session, inputId = "LOadings") 
    loading_heatmap()
  })
  
  output$downloadDataloadingheatmap <- downloadHandler(
    filename=paste("IC",as.integer(input$ic),"_Loadings",as.integer(input$loadings),".pdf",sep=""),
    content=function(file){
      pdf(file, height= 8, width=12)
      print(loading_heatmap())
      dev.off()
    },
    contentType="image/pdf"
  )
  
  
  region2loadings=reactive({
    ICA_result=ica$ICAres
    Genome=input$checkGenome
    Feature=as.integer(input$GenomicFeature)
    generate_region2loadings(ICA_result,Genome,Feature)
  })
  
  
  highLFlist_fun=reactive({
    ICA_result=ica$ICAres
    Target=as.integer(input$ic)
    Zscore=as.integer(input$loadings)
    Genome=input$checkGenome
    Feature=as.integer(input$GenomicFeature)
    generate_highLF_target_shiny(ICA_result,Target,Zscore,Genome,Feature)
  })
  
  loading_meth = reactive({
    ICA_result=ica$ICAres
    inFile_ICA=sampleInfo()
    Target=as.integer(input$ic)
    Zscore=as.integer(input$loadings)
    datamat=data_methylome()
    make_highLF_meth_target_shiny(ICA_result,inFile_ICA,datamat,Zscore,Target)
  })
  
  my_data<-reactive({
    testdata=highLFlist_fun()
    testdata=testdata[,-5]
    as.data.frame(cbind(View = shinyInput(actionButton, nrow(testdata),'button_', label = "View", onclick = 'Shiny.onInputChange(\"select_button\",  this.id)' ),testdata))
  })
  
  my_data_output<-reactive({
    testdata=highLFlist_fun()
    testdata=testdata[,-4]
    as.data.frame(testdata)
  })
  
  output$region2loadings <- DT::renderDataTable(region2loadings(),selection = 'single',options = list(searching = TRUE,pageLength = 10),server = FALSE, escape = FALSE,rownames= FALSE)
  
  output$highLFlist <- DT::renderDataTable(my_data(),selection = 'single',options = list(searching = TRUE,pageLength = 10),server = FALSE, escape = FALSE,rownames= FALSE)
  
  output$downloadhighLFlist <- downloadHandler(
    filename = function() {"highloadinglist.txt"},
    content = function(file) {
      write.table(my_data_output(), file, row.names = FALSE,sep="\t",quote=F)
    }
  )
  
  
  plotInput <- function(){loading_meth()}
  
  
  SelectedRow <- eventReactive(input$select_button,{
    as.numeric(strsplit(input$select_button, "_")[[1]][2])
  })
  
  # This is needed so that the button is clicked once for modal to show, a bug reported here
  # https://github.com/ebailey78/shinyBS/issues/57
  observeEvent(input$select_button, {
    toggleModal(session, "modalExample", "open")
  })
  
  
  #collist=c("blue","cyan","gray","yellow","red","pink","purple","lightgreen","magenta","green","orange")
  #collist=c("magenta","green","red","blue")
  set.seed(0414) 
  collist=grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)][sample(433,100,replace=F)]


  output$popup <- renderUI({
    print(input$select_button)
    ####bsModal("modalExample", paste("Methylation Plot",input$GenomicFeature, highLFlist_fun()$genesymbol[SelectedRow()], "p=",anova_p(),"(one-way anova)",sep="  "), "", size = "large",
    
    if (input$GenomicFeature==1){genomicFeature="CpG island"; Title=paste("Methylation Plot",genomicFeature,sep="  ")}
    if (input$GenomicFeature==2){genomicFeature="Gene body";  Title=paste("Methylation Plot",highLFlist_fun()$genesymbol[SelectedRow()], genomicFeature,sep="  ")}
    if (input$GenomicFeature==3){genomicFeature="Promoter";   Title=paste("Methylation Plot",highLFlist_fun()$genesymbol[SelectedRow()], genomicFeature,sep="  ")}
    
    TargetGroup=as.character(input$SampleProperty2)
    #dF=data.frame(sampleInfo()[,TargetGroup], as.numeric(plotInput()[SelectedRow(),])*100)
    dF=data.frame(sampleInfo()[,TargetGroup], as.numeric(plotInput()[SelectedRow(),]))
    colnames(dF)=c(TargetGroup,"meth")
    
    gp<-ggboxplot(dF,
                  colnames(dF)[1],colnames(dF)[2],
                  color = colnames(dF)[1],
                  palette=collist[1:NROW(sampleInfo()[,TargetGroup])],
                  add = "jitter",
                  add.params = list(size=2),
                  xlab="",
                  size=1,
                  ylim=c(0,100),ylab="Methylation level (%)") + 
                  font("ylab", size = 20, color = "black")+
                  font("y.text", size = 20, color = "black")+
                  #ylim=c(0,1),ylab="Methylation level") + 
    
      theme(axis.text.x = element_blank())
    
    if (NROW(unique(sampleInfo()[,TargetGroup]))>2){
      gp <- gp + stat_compare_means(method = "anova")
    }
    
    if (NROW(unique(sampleInfo()[,TargetGroup]))==2){
      gp <- gp + stat_compare_means(method = "t.test")
    }
    
    #bsModal("modalExample", paste("Methylation Plot",highLFlist_fun()$genesymbol[SelectedRow()], GenomicFeature,sep="  "), "", size = "large",
    bsModal("modalExample", Title, "", size = "large",
            
            column(12,
                   renderPlot({gp})
            )
    )
    
    
    
  })
})
