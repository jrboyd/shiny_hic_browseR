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

set.seed(0)
grs = lapply(c(40,50,60), function(x){
  signal = gen_signal(centers = x)
  gr = GRanges("a", IRanges(1:length(signal), 1:length(signal)), value = signal)
})
#for line tracks (and area), only the last track added animates but with the first track's color
lt1 = LineTrack(grs[[1]], label = "1", color = "red")
lt2 = LineTrack(grs[[2]], label = "2", color = "green")
lt3 = LineTrack(grs[[3]], label = "3", color = "blue")
TnTBoard(list(lt1, lt2, lt3, 
              merge(lt1, lt2), 
              merge(lt2, lt1), 
              merge(lt1, lt2, lt3),
              merge(lt3, lt2, lt1)), coord.range = c(-100, 200))

# #no problem with block tracks
# bt1 = BlockTrack(range(subset(grs[[1]], value > 0)), label = "1", color = "red")
# bt2 = BlockTrack(range(subset(grs[[2]], value > 0)), label = "2", color = "green")
# bt3 = BlockTrack(range(subset(grs[[3]], value > 0)), label = "3", color = "blue")
# TnTBoard(list(bt1, bt2, bt3, 
#               merge(bt1, bt2), 
#               merge(bt2, bt1), 
#               merge(bt1, bt2, bt3),
#               merge(bt3, bt2, bt1)), coord.range = c(-100, 200))