library(TnT)
gen_signal = function(centers, widths = rep(2, length(centers)), heights = 100*widths, size = 100){
  signal = rep(0, size)
  if(length(heights) != length(centers) & length(heights) != 1){
    stop("length of heights must equal lenght of  centers or 1")
  }
  if(length(widths) != length(centers) & length(widths) != 1){
    stop("length of widths must equal lenght of  centers or 1")
  }
  if(length(heights) == 1){
    heights = rep(heights, length(centers))
  }
  if(length(widths) == 1){
    widths = rep(widths, length(centers))
  }

  for(i in 1:length(centers)){
    x = centers[i]
    h = heights[i]
    w = widths[i]
    toadd = table(factor(round(rnorm(h, mean = x, sd = w)), levels = 1:size))
    signal = signal + toadd
  }
  return(as.numeric(signal))
}

peak_pos = c(30, 40, 60)

signal = gen_signal(centers = peak_pos, widths = c(1, 2, 1))
gr <- GRanges("chr12", IRanges(1:length(signal), 1:length(signal)))
gr$value = signal

signal2 = gen_signal(centers = peak_pos[c(1,3)], widths = c(1, 2, 1)[c(1,3)], heights = c(300, 50))
gr2 = GRanges("chr12", IRanges(1:length(signal2), 1:length(signal2)))
gr2$value = signal2

at = AreaTrack(gr, color = "blue", label = "area")
at2 = AreaTrack(gr2, color = "green", label = "area")
lt = LineTrack(gr, color = "blue", label = "line")
lt2 = LineTrack(gr2, color = "green", label = "line")
pt = PinTrack(gr[peak_pos], label = "pin")

tb = TnTBoard(list(at, at2, lt, pt))

library(TnT)
if (interactive() && require(shiny)) {
  ui <- fluidPage(fluidRow(
    column(width = 2, {
      "A Simple Example Here"
      # checkboxInput("Fill", "Fill?", value = T)
    }),
    column(width = 10, {
      TnTOutput("out")
    })
    
  ))
  server <- function (input, output) {
    output$out <- renderTnT({
      tb
    })
  }
  shinyApp(ui = ui, server = server)
}

