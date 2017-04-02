library(shiny)
library(ggplot2)
ALLSAMPLES<-c("babk","bima","eofe","eika","ieki","nusw","vass","pahc","joxm","iisa","hiaf","kuxp","podx","xugn","oapg","pamv","eipl","heja","pelm","hayt","oevr","fafq","sojd","hehd","diku","qolg","eiwy")#"fikt","qaqx"

shinyUI(fluidPage(
  titlePanel(title=h2("MHC Alternative Haplotypes", align="center"),windowTitle = "MHC Alternative Haplotypes"),

  sidebarLayout(
    sidebarPanel(
      
      conditionalPanel(condition="input.tabselected==1",
        textInput("mywindowsize","Insert height of plot (Dynamic)",value = 300),
        selectInput("Haps","Select the Haplotypes to Plot",choices=c("PGF","COX","APD","DBB","MANN","SSTO","QBL","MCF"),multiple = T,selectize = T),
        textInput("Start","Start Coordinate",value = 30100637), textInput("End","End Coordinate",value=30110637),
        radioButtons("redcd",label = "Collapse Alt Haps",choices = c("no","yes"),selected = "no"),
        actionButton("plotme","PLOT"),
        checkboxGroupInput("Sample","Select Samples to plot",as.list(ALLSAMPLES),selected = c("babk","bima","eofe","eika","ieki"))
      ),
      
      conditionalPanel(condition="input.tabselected==2",
                       selectInput("HapsT","Select the Haplotypes",choices=c("All","COX","APD","DBB","MANN","SSTO","QBL","MCF"),multiple = T,selectize = T,selected = "All"),
                       selectInput("State","Condition",choices = c("Macrophages-All","Macrophages-Naive","Macrophages-IFNg stimulated","LCLs"),multiple = F,selectize = T),
                       selectInput("Cats","Select the Category of peaks to show",choices=c("All","Recovery","ToFill","NotPresentInPGF","AltCntg-AltCntg","PGF-AltCntg","Insertion"),multiple = T,selectize = T, selected = "All"),
                       textInput("StartT","Start Coordinate",value = "-"), textInput("EndT","End Coordinate",value = "-"),
                       radioButtons("coordsys","Lifted to PGF Coordinates?",choices = c("yes","no")),
                       textInput("qrygene","Query by Annotated Gene (ENSEMBL ID or Gene Symbol)",value = "All"),
                       radioButtons("collapsetable",label="Collapse Alt Peaks?",choices = c("no","yes"),selected = "no"),
                       actionButton("querytable","Query")
        )
    ),
    
###################################
#conditionalPanel(condition="input.tabselected==1")
###################################
      mainPanel(
        tabsetPanel(
          tabPanel("Tracks", value=1,helpText("White background: Naive (Macrophages), Gray background: IFNg stimulated (Macrophages)"), uiOutput("Peaks.ui"),plotOutput("Genes")),
          tabPanel("Table", value=2,textOutput("dimtable"),tableOutput("TablePeaks")),
          id="tabselected"
        )
    )
  )
))
