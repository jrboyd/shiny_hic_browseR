require(TnT)

#### Track Constructor      ========
RangePairTrack <- function (range1, range2, label = deparse(substitute(range1)),
                        tooltip = mcols(range1), color = "blue", background = NULL,
                        height = 30) {
  
  data <- RangePairTrackData(range1 = range1, range2 = range2, tooltip = tooltip, color = color)
  new_track("RangePairTrack", background = background, height = height,
            label = label, data = data, display = list())
}

###  TnT Tracks   ##############################################################
#### Track classes          ========
setClass("RangePairTrack", contains = "RangeTrack")#, slots = c(Data = "RangePairTrackData"))

###  Track Data   ##############################################################
#### TrackData classes      ========
setClass("RangePairTrackData", contains = c("RangeTrackData"))

#### TrackData constructors     ========
RangePairTrackData <- function (range1, range2, color = "black", tooltip = mcols(range1), key = seq_along(range1)) {
  range1 <-
    if (is(range1, "IRanges"))
      GRanges(seqnames = "UnKnown", ranges = range1, strand = "*")
  else
    as(range1, "GRanges")
  
  range2 <-
    if (is(range2, "IRanges"))
      GRanges(seqnames = "UnKnown", ranges = range2, strand = "*")
  else
    as(range2, "GRanges")
  
  if(length(range1) != length(range2)){
    stop("length of ranges must be equal")
  }
  
  color <- {
    # TODO: Convert factor to integer or character? And which is better?
    if (!length(range1))
      color <- character(0)
    if (is.numeric(color))
      color <- as.integer(color)
    color
  }
  
  tooltip <-
    if (is.null(tooltip))
      data.frame(matrix( , nrow = length(range1), ncol = 0))
  else 
    as.data.frame(tooltip, optional = TRUE)
  
  mcols(range1) <- NULL
  range1$seqnames2 = seqnames(range2)
  range1$start2 = start(range2)
  range1$end2 = end(range2)
  range1$tooltip <- tooltip
  range1$color <- color
  range1$key <- key
  new("RangePairTrackData", range1)
}

setValidity("RangePairTrackData",
            function (object) {
              if (!is.data.frame(object$tooltip))
                if (is.null(object$tooltip))
                  return("Missing 'tooltip' meta-column in RangePairTrackData")
              else
                return("The 'tooltip' meta-column should be a data frame")
              if (!is.character(object$color) && !is.integer(object$color)) {
                if (is.null(object$color))
                  return("Missing 'color' meta-column in RangePairTrackData")
                else
                  return("The 'color' meta-column should be either character or integer")
              }
              if (is.null(object$seqnames2) | is.null(object$start2) | is.null(object$end2))
                  return("second ranges data is missing")
              k <- object$key
              if (is.null(k))
                return("Missing 'key' meta-column in RangePairTrackData")
              if (length(k) != length(unique(k)))
                return("'key' is not unique in RangePairTrackData")
              TRUE
            }
)

#### TrackData Compilation  ========
setMethod("compileTrackData", signature = "RangePairTrackData",
          function (trackData, full = FALSE) {
            stopifnot(length(unique(seqnames(trackData))) == 1)
            df <- as.data.frame(trackData, optional = TRUE) [
              c("start", "end", colnames(mcols(trackData)))]
            
            if (is(trackData, "TxTrackData"))
              jc.data <- jc(
                tnt.board.track.data.sync = ma(),
                retriever = jc(tnr.range_data_retriever = jc(tnr.cp_tx_color_to_exon = df))
              )
            else
              jc.data <- jc(
                tnt.board.track.data.sync = ma(),
                retriever = jc(tnr.range_data_retriever =
                                 ma(df, if (full) TRUE else FALSE))
              )
            jc.data
          }
)
