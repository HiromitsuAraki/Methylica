librarylist_cran=c("shinyBS","shinyjs","shinydashboard","R.utils","data.table","png",
                   "fields","rmarkdown","DT","ggpubr","sparklyr","flexdashboard","shinycssloaders","highcharter","fastICA","rgl","pals","circlize")

librarylist_bc=c("ComplexHeatmap","org.Hs.eg.db","org.Mm.eg.db")

is.exist_cran <- librarylist_cran %in% rownames(installed.packages())
is.exist_bc <- librarylist_bc   %in% rownames(installed.packages())
if(length(librarylist_cran[!(is.exist_cran)])>0){install.packages(librarylist_cran[!(is.exist_cran)])}
if(length(librarylist_bc[!(is.exist_bc)])>0){
  install.packages("BiocManager")
  BiocManager::install(librarylist_bc[!(is.exist_bc)])
}


# Load packages ----
library(shiny)
library(shinyBS)
library(shinyjs)
library(shinydashboard)
library(R.utils)
library(data.table)
library(png)
library(fields)
library(circlize)
library(ComplexHeatmap)
library(rmarkdown)
library(DT)
library(org.Hs.eg.db)
library(ggpubr)
library(sparklyr)
library(flexdashboard)
library(shinycssloaders)
library(highcharter)


#Source ICA pipeline ----
source("Methylica_pipeline.R")


###UI
shinyUI(dashboardPage(skin="blue",

  ##Header Start
  dashboardHeader(title = "Methylica"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data upload", tabName = "DAtaupload", icon = icon("th")),
      menuItem("Parameters",  tabName = "PArameters", icon = icon("th")),
      menuItem("ICs",         tabName = "ics",        icon = icon("th")),
      menuItem("Loadings",    tabName = "LOadings",   icon = icon("th"))
      #menuItem("Loadings",    tabName = "LOadings",   icon = icon("th")),
      #menuItem("README",      tabName = "REadme",     icon = icon("th"))
    )
  ),
  ##Header End



  ##Main body Start
  dashboardBody(
    tabItems(
 
      # 1st tab content
      tabItem(tabName = "DAtaupload",
          fileInput("data_input",  label="Upload methylome data",    accept=c('text/csv', 'text/comma-separated-values,text/plain', '.txt')),
          fileInput("sample_input",label="Upload sample information",accept=c('text/csv', 'text/comma-separated-values,text/plain', '.txt'))
      ),
      
      
      # 2nd tab content
      tabItem(tabName = "PArameters",

              #Platforms              
              radioButtons("checkPlatform", 
                           label = h3("Select Platform"), 
                           choices = list("Sequence-based" = 1, 
                                          "Infinium array (450K or EPIC)" = 2),                           
                           selected = ""),
              br(),
              
              #Reference version              
              radioButtons("checkGenome", 
                           label = h3("Select Genome"), 
                           choices = list("hg38" = 1, 
                                          "hg19" = 2, 
                                          "mm10" = 3,
                                          "mm9"  = 4),                           
                           selected = ""),
              br(),


              #Genomic feature (select from CGI, Gene body or Promoter)
              radioButtons("GenomicFeature", 
                           label = h3("Select Feature"), 
                           choices = list("CGI"       = 1, 
                                          "Gene body" = 2, 
                                          "1st intron"= 3,
                                          "Promoter"  = 4),                           
                           selected = ""),


              
              br(),
              tags$h3("IC num"),
              numericInput("ICnum",     label = h5("User definition"),value = 1,width='20%'),
              checkboxInput("ICnumPCA", label ="PCA based definition"),
                              
              br(),
              br(),
              br(),
              actionButton("goButton1", "Run")
              
      ),
      
      
      # 3rd tab content (Result of ICA)
      tabItem(tabName = "ics",

              #visualization of ICA clustering
              h3("Heatmap clustering "),
              fluidRow(       
                column(5,

                         #conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                        　#tags$h1("Loading...")),
                                 
                         #imageOutput('ICA_HC',height = "500px",width = "1000px")  %>% withSpinner(type = 4,color="blue"),
                         withSpinner(imageOutput('ICA_HC',height = "500px",width = "1000px"),type = 4,color="blue"),
                         #imageOutput('ICA_HC',height = "500px",width = "1000px"),
                        
                         downloadButton('downloadICA_HC', 'Download')
                )
               ),
               
               br(),br(),br(),
               
              #visualization of A boxplot
              h3("Box plot of each IC"),
              fluidRow(
                column(5,
                        selectInput("SampleProperty", label = h5("Sample property"), choices = 'variables'),
                        withSpinner(imageOutput('ICboxplot',height='600px',width =  '800px'),type = 4,color="blue"),
                        #imageOutput('ICboxplot',height='1000px',width =  '1000px'),
                        downloadButton('downloadICAboxplot', 'Download')
                )
              ),
              br(),br(),
              
              #Table list of loading factor of each regions in all ICs
              h3("IC table"),
              fluidRow(
                column(5,
                      DT::dataTableOutput("region2loadings") 
                )
              )
              
      ),
      
                        
      # 4th tab content (Visuzaliation of methylation level of regions with high loadings)
      tabItem(tabName = "LOadings",

              #visualation of methylation heatmap      
              fluidRow(
                column(1,
                       numericInput("ic", label = h3("IC"), value = 1)),
                
                column(1,
                       numericInput("loadings", label = h3("Loadings"), value = 3))
                
              ),


              #download of methylation heatmap    
              fluidRow(
                column(4,
                   withSpinner(plotOutput('loadingheatmap',height = "500px",width = "800px"), type=4,color="blue"),
                   #plotOutput('loadingheatmap',height = "500px",width = "800px"),
                   downloadButton('downloadDataloadingheatmap', 'Download')
                )
              ),

              br(),br(),
              
              
              #show high loadings table   
              selectInput("SampleProperty2", label = p("Sample property"), choices = 'variables',width='15%'),
              fluidRow(
                column(4,
                  withSpinner(DT::dataTableOutput("highLFlist"),type=4,color="blue"), 
                  uiOutput("popup"),
                  downloadButton('downloadhighLFlist','Download')
                )
              )
       )       
      #),


      # 5th tab content (README)
      #tabItem(tabName = "REadme",
      #        fluidPage(
                ###includeMarkdown("./Untitled.rmd")
              #)
      #)
      
                              
    )
   )
  )
)

  
