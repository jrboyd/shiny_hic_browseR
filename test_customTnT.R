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

new_track <- function (class, data, display = list(), background = NULL, height = NULL, label = NULL) {
  track <- new(Class = class, Data = data, Display = display)
  if (length(label) > 1L)
    label <- paste(label, collapse = " ")
  trackSpec(track, "background") <- background
  trackSpec(track, "height")     <- height
  trackSpec(track, "label")      <- label
  track
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




signal2hic = function(signal, max_dist = 20){
  gr = GRanges("chr12", IRanges(1:length(signal), 1:length(signal)))
  gr1 = NULL
  gr2 = NULL
  for(d in 0:max_dist){
    for(i in 1:(length(gr) - d)){
      s = i
      e = i + d
      g1 = gr[s]
      g2 = gr[e]
      v = min(signal[s], signal[e])
      g1$value = v
      if(is.null(gr1)){
        gr1 = g1
      }else{
        gr1 = c(gr1, g1)
      }
      if(is.null(gr2)){
        gr2 = g2
      }else{
        gr2 = c(gr2, g2)
      }
      
    }
  }
  return(list(ranges1 = gr1, ranges2 = gr2))
}

hic_list = signal2hic(signal, max_dist = 20)
rng_hic = range(hic_list$ranges1$value)
hic_colors = rgb(colorRamp(c("slategray", "blue", "purple", "red"))(hic_list$ranges1$value / max(rng_hic))/255)
hic_list$ranges1$color = hic_colors

rpt = RangePairTrack(hic_list$ranges1, hic_list$ranges2)
tb = TnTBoard(list(at, at2, lt, pt))
tb
tb = TnTBoard(list(at, at2, lt, pt, rpt))
tb
widget <- trackWidget(tb, elementId = NULL)
tntdef <- compileBoard(tntdef)


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

