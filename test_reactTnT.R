library(TnT)
if (interactive() && require(shiny)) {
  ui <- fluidPage(fluidRow(
    column(width = 2, {
      "A Simple Example Here"
      checkboxInput("Fill", "Fill?", value = T)
    }),
    column(width = 10, {
      TnTOutput("out")
    })
    
  ))
  server <- function (input, output) {
    signal = c(rep(0, 5), table(round(rnorm(100, mean = 10))), rep(0, 5))
    gr <- GRanges("chr12", IRanges(1:length(signal), 1:length(signal)))
    gr$score = signal
   
    
    re.btrack <- reactive({
      
      
      if(input$Fill){
        AreaTrack(gr, value = gr$score, color = "blue")
      }else{
        LineTrack(gr, value = gr$score, color = "blue")
      }
    })
    output$out <- renderTnT({
      TnTGenome(re.btrack())
    })
  }
  shinyApp(ui = ui, server = server)
}