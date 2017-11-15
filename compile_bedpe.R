# sapply(dir("TnT/R/", full.names = T), source)
#' 
#' setGeneric("wakeupTrack", function (track) standardGeneric("wakeupTrack"))
#' 
match.class <- function (class, choices) {
  # Match a class to a list of known classes
  #   (i.e. matching the first direct superclass in that list)
  #   It can be used to simulate a method dispatch.
  super.classes <- c(class, getAllSuperClasses(getClass(class)))
  matched <- choices[which.min(match(choices, super.classes))]

  if (!length(matched))
    stop(sprintf("Unmatched class %s", sQuote(class)))
  if (! matched %in% choices) # Sanity check
    stop("<internal> Wrongly matched class")
  matched
}
#' 
setMethod("wakeupTrack", signature = c(track = "RangeTrack"),
          function (track) {
            class <- class(track)
            print("exe")
            # Simulate a method dispatch
            use.class <- match.class(class,
                                     c("BlockTrack", "GeneTrack", "TxTrack", "VlineTrack",
                                       "PinTrack", "LineTrack", "AreaTrack", "RangePairTrack")
            )

            feaname <- switch(use.class,
                              BlockTrack   = "tnt.board.track.feature.block",
                              GeneTrack    = "tnt.board.track.feature.genome.gene",
                              TxTrack      = "tnt.board.track.feature.genome.transcript",
                              VlineTrack   = "tnt.board.track.feature.vline",
                              PinTrack     = "tnt.board.track.feature.pin",
                              LineTrack    = "tnt.board.track.feature.line",
                              AreaTrack    = "tnt.board.track.feature.area",
                              RangePairTrack = "new_bedpe",
                              stop("<internal> Unmatched track class")
            )

            ### TEMP: TO REMOVE IN FUTURE, AND MAKE SURE KEY IS UNIQUE
            if (use.class == "VlineTrack") {
              trackData(track)$key <- start(trackData(track))
            }
            ### TEMP: TO REMOVE IN FUTURE

            track <- .initDisplay(track, feaname = feaname)
            track <- .convertCol(track)
            track
          }
)
#' 
#' compileTrack <- function (tntTrack, wakeup = TRUE) {
#'   # Wake up
#'   if (wakeup)
#'     tntTrack <- wakeupTrack(tntTrack)
#'   
#'   label <- trackSpec(tntTrack, "label")
#'   label <- if (is.null(label)) "" else label
#'   
#'   background <- trackSpec(tntTrack, "background")
#'   background <- if (is.null(background)) col2hex("white")
#'   else                     col2hex(background)
#'   
#'   height <- trackSpec(tntTrack, "height")
#'   height <- if (is.null(height)) 100 else height
#'   
#'   li.spec <- list(
#'     tnt.board.track = ma(),
#'     label  = label,
#'     color  = background,  # rename to "color"
#'     height = height
#'   )
#'   
#'   jc.spec <- asJC(li.spec)
#'   jc.display <- jc(
#'     display = asJC(
#'       tntTrack@Display
#'     )
#'   )
#'   jc.data <- jc(data = compileTrackData(trackData(tntTrack)))
#'   c(jc.spec, jc.display, jc.data)
#' }
#' 
#' 
#' .initDisplay <- function (track, feaname, extra = list()) {
#'   # For RangeTrack
#'   
#'   di.init <- setNames(list(ma()), feaname)
#'   di.color <- list(color = {
#'     ## LineTrack and AreaTrack do not support multiple color values
#'     ## This is a tempoary solution
#'     if (is(track, "LineTrack") || is(track, "AreaTrack")) {
#'       co <- unique(trackData(track)$color)
#'       if (length(co) > 1) {
#'         # TODO: Also add warning in constructor?
#'         warning("LineTrack and AreaTrack do not support multiple color values")
#'         co <- co[1]
#'       }
#'       if (!length(co))
#'         NULL  ## empty track?
#'       else
#'         co
#'     }
#'     else
#'       js('function (d) {return d.color;}')
#'   })
#'   di.index <- list(index = js('function (d) {return d.key;}'))
#'   di.extra <- extra
#'   
#'   # Do not set domain in display any more, but normalize data from domain to 0..1
#'   di.domain <- list(domain = {
#'     if (is(track, "PinTrack")) # TODO: This can be removed when
#'       c(0L, 1L)                    #       default value in upstream is fixed.
#'     else
#'       NULL
#'   })
#'   
#'   di.tooltip <- {
#'     toolti <- tooltip(track)
#'     stopifnot(
#'       is.data.frame(toolti),
#'       #all(sapply(toolti, is.atomic)),
#'       !any(duplicated(names(toolti)))
#'     )
#'     toolti.header <- {
#'       label <- trackSpec(track, "label")
#'       if (is.null(label)) ""
#'       else                label
#'     }
#'     toolti.entries <- colnames(toolti)
#'     list(on = ma("click",
#'                  tooltipCallback(header = toolti.header, entries = toolti.entries)
#'     ))
#'   }
#'   
#'   display <- c(di.init, di.color, di.index, di.extra, di.domain, di.tooltip)
#'   track@Display <- display
#'   track
#' }
#' 
#' .convertCol <- function (track) {
#'   # TODO: if all the colors are identical,
#'   #       we may modify the color callback to reduce the size of file
#'   col <- trackData(track)$color
#'   if (!length(col))
#'     return(track)
#'   col <- col2hex(col)
#'   trackData(track)$color <- col
#'   track
#' }
#' 
#' tooltipCallback <- function (header, entries) {
#'   # Example:
#'   #   tooltipCallback(header = "Tooltip Header",
#'   #                   entries = c("Start", "End", "Description"))
#'   stopifnot(length(header) == 1)
#'   jc(tnr.tooltip_callback = ma(header, entries))
#' }
#' 
#' 
#' #' Save a TnTBoard to an HTML file
#' #' 
#' #' A simple wrapper of \code{\link[htmlwidgets]{saveWidget}}, which saves a
#' #' TnTBoard/TnTGenome object to an HTML file (e.g. for sharing with others).
#' #'
#' #' @param tntdef A TnTBoard/TnTGenome object to save.
#' #' @param file,selfcontained,libdir,background,knitrOptions
#' #'     Passed to \code{\link[htmlwidgets]{saveWidget}}.
#' #'
#' #' @return Return NULL.
#' #' @export
#' #' @examples
#' #' data <- GRanges("chr2", IRanges(c(6,9,42), width = 1),
#' #'                 value = c(0.3, 0.5, 0.9))
#' #' track <- PinTrack(data, label = NULL, background = "green")
#' #' genome <- TnTGenome(list(track))
#' #' destfile <- tempfile(fileext = ".html")
#' #' destfile
#' #' saveTnT(genome, destfile)
#' #' \dontrun{
#' #' utils::browseURL(destfile)
#' #' }
#' saveTnT <- function (tntdef, file, selfcontained = TRUE, libdir = NULL,
#'                      background = "white", knitrOptions = list()) {
#'   stopifnot(is(tntdef, "TnTBoard"))
#'   tntdef <- trackWidget(tntdef)
#'   htmlwidgets::saveWidget(tntdef, file, selfcontained, libdir, background, knitrOptions)
#' }
#' 
#' 
#' 
#' #' Scale Qualitative Values to Color
#' #' 
#' #' A simple util function that scales a factor to color based on the palette function.
#' #'
#' #' @param value A factor or character vector that may have n unique values.
#' #' @param palette.fun The palette function to generate colors.
#' #'     For example, \code{\link[grDevices]{terrain.colors}}.
#' #' @param ... Extra arguments passed to the palette function.
#' #'
#' #' @return
#' #'     A character vector as colors, with the same length of \code{value}. Same values
#' #'     in \code{value} will have the same color.
#' #' @export
#' #' @examples
#' #' mapcol(iris$Species)
#' mapcol <- function (value, palette.fun = grDevices::rainbow, ...) {
#'   uniqueV <- unique(value)
#'   mappedC <- palette.fun(length(uniqueV), ...)
#'   mappedC[match(value, uniqueV)]
#' }
#' 
#' #' Display Labels with Strand
#' #' 
#' #' A simple util function that used internally to generate display labels of GeneTrack and
#' #' TxTrack.
#' #'
#' #' @param labels Character vector, names of each feature.
#' #' @param strands Factor or character vector with the same length of \code{labels}, can
#' #'     be "+", "-" or "*".
#' #' @return
#' #'     A character vector that combines the labels with strand information.
#' #' @export
#' #' @examples
#' #' strandlabel(c("gene1", "gene2", "gene3"), c("+", "-", "*"))
#' strandlabel <- function (labels, strands) {
#'   if (length(labels))
#'     ifelse(strands == "+", paste(labels, ">"),
#'            ifelse(strands == "-", paste("<", labels), labels))
#'   else
#'     ifelse(strands == "+", ">",
#'            ifelse(strands == "-", "<", ""))
#' }
#' 
#' 
#' col2hex <- function (colors) {
#'   mat <- grDevices::col2rgb(colors)
#'   grDevices::rgb(mat[1, ], mat[2, ], mat[3, ], maxColorValue = 255)
#' }
#' 
#' splitdf <- function (df, f) {
#'   # We need to speed up this function
#'   if (requireNamespace("data.table", quietly = TRUE)) {
#'     # relative faster
#'     df[["___FACTOR___"]] <- f
#'     ldf <- split(   # data.table:::split.data.table(
#'       data.table::as.data.table(df),
#'       drop = TRUE, by = "___FACTOR___", keep.by = FALSE, flatten = FALSE
#'     )
#'   }
#'   else {
#'     warning("Install data.table could speed up the spliting process")
#'     ldf <- split(df, f, drop = TRUE)
#'   }
#'   ldf
#' }
#' 
#' 
#' 
#' .consolidateSeqinfo <- function (li.tracks) {
#'   # The function accepts a list of tracks,
#'   # returns a list of tracks with the same Seqinfo
#'   stopifnot(is.list(li.tracks))
#'   
#'   if (!length(li.tracks))
#'     return(li.tracks)
#'   s <- .mergeSeqinfo(li.tracks)
#'   for (i in seq_along(li.tracks))
#'     if (!identical(seqinfo(li.tracks[[i]]), s))
#'       # The merged Seqinfo should be a superset of individual Seqinfo
#'       seqinfo(li.tracks[[i]], new2old = match(seqlevels(s), seqlevels(li.tracks[[i]])),
#'               pruning.mode = "error") <- s
#'   li.tracks
#' }
#' 
#' .mergeSeqinfo <- function (li.seqinfo) {
#'   # It accepts a list of Seqinfo object or a list of tracks, returns the merged Seqinfo.
#'   # It checks whether the objects are identical before applying `merge`,
#'   # so that it is fast in some cases.
#'   if (!length(li.seqinfo))
#'     return(Seqinfo())
#'   for (i in seq_along(li.seqinfo))
#'     if (!is(li.seqinfo[[i]], "Seqinfo"))
#'       li.seqinfo[[i]] <- seqinfo(li.seqinfo[[i]])
#'     
#'     target <- li.seqinfo[[1]]
#'     for (i in seq_along(li.seqinfo)) {
#'       new <- li.seqinfo[[i]]
#'       if (identical(target, new))
#'         target <- new
#'       else
#'         target <- merge(target, new)
#'     }
#'     target
#' }
#' 
#' #### TrackData Compilation  ========
#' 
#' setGeneric("compileTrackData",
#'            function (trackData, ...) standardGeneric("compileTrackData"))
#' 
#' # setMethod("compileTrackData", signature = "NoTrackData",
#' #     function (trackData)
#' #         jc(tnt.board.track.data.empty = na)
#' # )
#' 
#' setMethod("compileTrackData", signature = "RangeTrackData",
#'           function (trackData, full = FALSE) {
#'             stopifnot(length(unique(seqnames(trackData))) == 1)
#'             df <- as.data.frame(trackData, optional = TRUE) [
#'               c("start", "end", colnames(mcols(trackData)))]
#'             
#'             if (is(trackData, "TxTrackData"))
#'               jc.data <- jc(
#'                 tnt.board.track.data.sync = ma(),
#'                 retriever = jc(tnr.range_data_retriever = jc(tnr.cp_tx_color_to_exon = df))
#'               )
#'             else
#'               jc.data <- jc(
#'                 tnt.board.track.data.sync = ma(),
#'                 retriever = jc(tnr.range_data_retriever =
#'                                  ma(df, if (full) TRUE else FALSE))
#'               )
#'             jc.data
#'           }
#' )
#' 
#' setMethod("compileTrackData", signature = "PosTrackData",
#'           function (trackData, full = FALSE) {
#'             stopifnot(length(unique(seqnames(trackData))) == 1)
#'             stopifnot(all(width(trackData) == 1))
#'             
#'             df <- as.data.frame(trackData, optional = TRUE)[c("start", colnames(mcols(trackData)))]
#'             df <- S4Vectors::rename(df, c(start = "pos"))
#'             
#'             jc.data <- jc(
#'               tnt.board.track.data.sync = ma(),
#'               retriever = jc(tnr.pos_data_retriever =
#'                                ma(df, if (full) TRUE else FALSE))
#'             )
#'             jc.data
#'           }
#' )
#' 
#' setMethod("compileTrackData", signature = "PosValTrackData",
#'           function (trackData, full = FALSE) {
#'             stopifnot(length(unique(seqnames(trackData))) == 1)
#'             stopifnot(all(width(trackData) == 1))
#'             validObject(trackData)
#'             
#'             df <- as.data.frame(trackData, optional = TRUE)[c("start", colnames(mcols(trackData)))]
#'             df <- S4Vectors::rename(df, c(start = "pos"))
#'             
#'             domain <- getdomain(trackData)
#'             
#'             jc.data <- jc(
#'               tnt.board.track.data.sync = ma(),
#'               retriever = jc(
#'                 tnr.pos_data_retriever = ma(
#'                   jc(tnr.scale_val = ma(df, domain)), # scale data from domain to [0, 1]
#'                   if (full) TRUE else FALSE  # return full data or not
#'                 )
#'               )
#'             )
#'             jc.data
#'           }
#' )
#' 
#' 
#' #### The following functions are no longer used. -------------------------------
#' 
#' # ul <- function(x)
#' #     unlist(x, recursive = FALSE, use.names = FALSE)
#' # 
#' # `ul<-` <- function (x, value)
#' #     relist(value, x)
#' # 
#' # expose_all <- function (package = "TnT") {
#' #     attachname <- paste0(package, "_all")
#' #     while(attachname %in% search())
#' #         detach(attachname, character.only = TRUE)
#' #     
#' #     pkgns <- loadNamespace(package)
#' #     attach(pkgns, name = attachname)
#' #     invisible(pkgns)
#' # }
#' # 
#' # .removeAsIs <- function (df) {
#' #     # The nested data frame converted from GRanges/DataFrame will have a
#' #     # "AsIs" class, which will cause the data frame can not be shown and
#' #     # can not be converted to JSON correctly.
#' #     for (i in seq_along(df)) {
#' #         element <- df[[i]]
#' #         if (is.data.frame(element)) {
#' #             class(element) <- class(element)[class(element) != "AsIs"]
#' #             df[[i]] <- .removeAsIs(element)
#' #         }
#' #     }
#' #     df
#' # }
#' # 
#' # .JSONFilter <- function (colname) {
#' #     stopifnot(is.character(colname))
#' #     escapeColname <- sapply(colname, function (s) as.character(toJSON(unbox(s))))
#' #     condfilter <- paste(sprintf("[%s]", escapeColname), collapse = "")
#' #     condfilter
#' # }
#' # if (interactive()) .JSONFilter(colname = c("data", "start"))