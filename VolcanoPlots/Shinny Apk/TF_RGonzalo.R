#El objetivo de esta aplicación es que un usuario pueda subir un archivo resultante de por ejemplo un análisis
#de microarrays, pueda generar un volcano plot, y pueda cambiar los colores de los genes representados según
#los valores de adj-pvalue y Fold change elegidos por él mismo

library(shiny)
library(calibrate)


ui<-shinyUI(fluidPage(
  titlePanel("Configure the parameters of the Volcano Plot"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file","Please load your file (.csv)",accept='csv'),
      sliderInput("Bvalue", "Choose the desired B value", 0, 30, 2,step = 0.5),
      sliderInput("FC", "Choose the positve Fold Change", 0, 5, 2,step = 0.5),
      sliderInput("FCm", "Choose the negative Fold Change", -10, 0, -2,step = 0.5),
      selectInput("color", "Choose the color for the DE genes", c("red","green","blue"))
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Table", tableOutput("tab")),
        tabPanel("Graphic", plotOutput("result"),with=500,height=500, downloadButton('down', 'Download'))
        )
    )
  )
)
)

server<-function(input,output){
  dd<-reactive({
    inFile<-input$file
    if (is.null(inFile)) return(invisible(NULL))
    else
       read.csv2(input$file$datapath,sep = ";",header=TRUE)
    })
  
  output$tab <- renderTable(dd()[c(1:15),c(1:9)])
  
  
  output$result<-renderPlot({
    plot(dd()$logFC, dd()$B, pch=20, main="Volcano plot", xlab="logFC",ylab="B", xlim=c(min(dd()$logFC),max(dd()$logFC))) 
    
    with(subset(dd(), B>input$Bvalue & abs(logFC>input$FC)), points(logFC, B, pch=20, col=input$color))
    with(subset(dd(), B>input$Bvalue & logFC<input$FCm), points(logFC, B, pch=20, col=input$color))
    
    with(subset(dd(), B>input$Bvalue & abs(logFC>input$FC)), textxy(logFC, B, labs=SymbolsA, cex=.7,offset=0.3))
    with(subset(dd(), B>input$Bvalue & logFC<input$FCm), textxy(logFC, B, labs=SymbolsA, cex=.7,offset=0.3))
    
    abline(v=input$FCm,lty=3)
    abline(v=input$FC,lty=3)
  })
   
  
  output$down <- downloadHandler(
    filename = function() "figura.pdf",
    content = function(ff) {
      pdf(ff)
      plot(dd()$logFC, dd()$B, pch=20, main="Volcano plot", xlab="logFC",ylab="B", xlim=c(min(dd()$logFC),max(dd()$logFC))) 
      
      with(subset(dd(), B>input$Bvalue & abs(logFC>input$FC)), points(logFC, B, pch=20, col=input$color))
      with(subset(dd(), B>input$Bvalue & logFC<input$FCm), points(logFC, B, pch=20, col=input$color))
      
      with(subset(dd(), B>input$Bvalue & abs(logFC>input$FC)), textxy(logFC, B, labs=SymbolsA, cex=.7,offset=0.3))
      with(subset(dd(), B>input$Bvalue & logFC<input$FCm), textxy(logFC, B, labs=SymbolsA, cex=.7,offset=0.3))
      
      abline(v=input$FCm,lty=3)
      abline(v=input$FC,lty=3)
      dev.off()
    })
}

shinyApp(ui = ui, server = server)

