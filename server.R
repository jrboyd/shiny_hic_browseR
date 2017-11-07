library(shiny)
library(TnT)
library(GenomicRanges)
library(rtracklayer)
# source("class_HiC_matrix_wInsulation.R")
# source("class_integrated_plots.R")
# load(".RData")
# ds_hic = data_source.hic_matrix(my_p_hicsd[[2]])
# ds_ins = data_source.hic_insulation(my_p_hicsd[[2]])
# runx_bws = dir("/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/breast/10A_progression", pattern = "RUNX.+_FE", full.names = T)
# ctcf_bws = dir("/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/breast/10A_progression", pattern = "CTCF.+_FE", full.names = T)
# names(runx_bws) = c("MCF10A", "AT1", "CA1a")
# names(ctcf_bws) = c("MCF10A", "AT1", "CA1a")
# ds_bw = data_source.bigwig(runx_bws[1])
# ds_bw3r = data_source.bigwig(runx_bws)
# ds_bw3c = data_source.bigwig(ctcf_bws)
# 
# load("gencode.v24.gene_dt.save")
# gene_gr = GRanges(dt)
# ds_ref = data_source.ref(gene_gr)

qgr = GRanges("chr6", IRanges(25*10^6, 25.2*10^6))
#bigwig files
bw_files = dir("/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/breast/MDA231_MK_runx/", full.names = T, pattern = "ed_FE.bw")
names(bw_files) = basename(bw_files)
bw_files = as.list(bw_files)
#narrowPeak files
np_files = c("/slipstream/galaxy/uploads/working/qc_framework/output_MK_MDA231_Runx/MDA231_Runx1_pooled/MDA231_Runx1_pooled_peaks.narrowPeak",
             "/slipstream/galaxy/uploads/working/qc_framework/output_MK_MDA231_Runx/MDA231_Runx2_pooled/MDA231_Runx2_pooled_peaks.narrowPeak")
names(np_files) = basename(np_files)
np_files = as.list(np_files)
score2color = function(scores, colors = c("white", "black"), scale_max = max(scores), scale_min = 0){
  scaled = (scores - scale_min) / scale_max
  scaled = ifelse(scaled > 1, 1, scaled)
  scaled = ifelse(scaled < 0, 0, scaled)
  rgb(colorRamp(colors)(scaled / max(scaled))/255)
}
bedgraph2line = function(bdg_gr, win = 50){
  #make bdg_gr amenable to window size
  query_gr = range(bdg_gr)
  start(query_gr) = floor(start(query_gr) / win) * win
  end(query_gr) = ceiling(end(query_gr) / win) * win
  end(query_gr) = end(query_gr) - 1
  
  line_qgr = tile(query_gr, width = win)[[1]]
  mids = mid(ranges(line_qgr)) + 1
  start(line_qgr) = mids
  end(line_qgr) = mids
  olap = findOverlaps(line_qgr, bdg_gr)
  line_qgr = line_qgr[queryHits(olap)]
  line_qgr$score = bdg_gr$score[subjectHits(olap)]
  return(line_qgr)
}
bedgraph2lineTrack = function(bdg_gr, win = 50, label = "Track Name", color = "black", floor_at = 0 ){
  bdg_gr = subset(bdg_gr, score > floor_at)
  line_gr = bedgraph2line(bdg_gr = bdg_gr, win = win)
  LineTrack(line_gr, value = line_gr$score, color = color, label = label)
}

bedgraph2blockTrack = function(bdg_gr, label = "Track Name", color = "black", floor_at = 2){
  bdg_gr$color = score2color(bdg_gr$score, scale_min = 1, colors = c("white", color))
  bdg_gr$tooltip = round(bdg_gr$score, 2)
  end(bdg_gr) = end(bdg_gr) + 1
  bdg_gr = subset(bdg_gr, score > floor_at)
  BlockTrack(bdg_gr, label = label, color = bdg_gr$color, tooltip = data.frame(FE = bdg_gr$tooltip))
}
peak2blockTrack = function(np_f, label = "Track Name", color = "gray"){
  df = read.table(np_f)
  colnames(df) = c("seqnames", "start", "end", "peak_id", "score", "dot", "fe", "p", "p-adj", "rel_summit")
  gr = GRanges(df)
  sum_gr = gr
  start(sum_gr) = start(sum_gr) + sum_gr$rel_summit
  end(sum_gr) = start(sum_gr)
  gr$color = color
  sum_gr$color = "black"
  gr = c(gr, sum_gr)
  gr = sort(gr)
  BlockTrack(gr, label = label, color = gr$color, tooltip = data.frame(ID = gr$peak_id, p = gr$p))
}
# runx1_gr$color = score2color(runx1_gr$score)
# runx1_gr$tooltip = round(runx1_gr$score, 2)
# bt_range = BlockTrack(qgr)
# bt_runx1 = BlockTrack(subset(runx1_gr, score > 2), label = "Runx1")
# bt_runx2 = BlockTrack(subset(runx2_gr, score > 2), label = "Runx2")
# TnTGenome(list(bt_runx1, bt_runx2), view.range = qgr)
# TnTBoard(bt_range)
# lt1 = LineTrack(runx1_gr, value = runx1_gr$score, )
function(input, output, session) {
  
  output$TnTTest = renderTnT({
    qgr = GRanges("chr6", IRanges(25.1*10^6, 25.5*10^6))
    #import using rtracklayer
    runx1_gr = import.bw(bw_files$MDA231_Runx1_pooled_FE.bw, which = qgr)
    runx2_gr = import.bw(bw_files$MDA231_Runx2_pooled_FE.bw, which = qgr)
    #block tracks
    bt_runx1 = bedgraph2blockTrack(runx1_gr, label = "Runx1 FE", color = "green")
    bt_runx2 = bedgraph2blockTrack(runx2_gr, label = "Runx2 FE", color = "purple")
    #line tracks
    lt_runx1 = bedgraph2lineTrack(runx1_gr, label = "Runx1 FE", color = "green", win = 50, floor_at = 0)
    lt_runx2 = bedgraph2lineTrack(runx2_gr, label = "Runx2 FE", color = "purple", win = 50, floor_at = 0)
    #peak tracks 
    np_runx1 = peak2blockTrack(np_files$MDA231_Runx1_pooled_peaks.narrowPeak, label = "Runx1 Peaks", color = "green")
    np_runx2 = peak2blockTrack(np_files$MDA231_Runx2_pooled_peaks.narrowPeak, label = "Runx2 Peaks", color = "purple")
    
    vgr = GRanges("chr6", IRanges(25.2789*10^6-1800, 25.2793*10^6+1800))
    genome(seqinfo(vgr)) = "hg38"
    TnTGenome(tracklist = list(bt_runx1, bt_runx2, lt_runx1, lt_runx2, np_runx1, np_runx2)[c(5,1,3,6,2,4)], 
              zoom.allow = c(1000, 200000),
              coord.range = ranges(qgr),
              view.range = vgr)
  })
  # current_style = reactive({
  #   qgr = GRanges("chr10", IRanges(12*10^6, 14*10^6))
  #   qgr = GRanges("chr10", IRanges(12*10^6, 26*10^6))
  #   bin_size = ds_hic@hic_mat@parameters@bin_size
  #   switch(input$style,
  #          a = list(
  #            plot_from_data_source(ds_hic, qgr),
  #            plot_from_data_source(ds_ins, qgr),
  #            plot_from_data_source(data_source.bigwig(ds_bw3r@bigwig_path[2]), qgr, bin_size = bin_size, alpha = 1, bin_method = "max", bigwig_title = "RUNX1 FE", bigwig_colors = c("red", "green", "blue")),
  #            plot_from_data_source(data_source.bigwig(ds_bw3c@bigwig_path[2]), qgr, bin_size = bin_size, alpha = 1, bin_method = "max", bigwig_title = "CTCF FE", bigwig_colors = c("red", "green", "blue")),
  #            plot_from_data_source(ds_ref, qgr, olap_extension = .1, text_vjust = .5, text_hjust = -.2,  text_x_relative = 1, text_y_relative = 0, text_angle = 0)
  #          ),
  #          b = list(
  #            plot_from_data_source(ds_hic, qgr),
  #            plot_from_data_source(ds_ins, qgr),
  #            plot_from_data_source(ds_bw3r, qgr, bin_size = bin_size, alpha = .3, bin_method = "max", bigwig_title = "RUNX1 FE", bigwig_colors = c("red", "green", "blue")),
  #            plot_from_data_source(ds_bw3c, qgr, bin_size = bin_size, alpha = .3, bin_method = "max", bigwig_title = "CTCF FE", bigwig_colors = c("red", "green", "blue")),
  #            plot_from_data_source(ds_ref, qgr, strand_colors = c("darkgreen", "darkorange"), label_method = "", top_spacer = 0)
  #          ),
  #          c = list(
  #            plot_from_data_source(ds_hic, qgr),
  #            plot_from_data_source(ds_ins, qgr),
  #            plot_from_data_source(ds_bw3r, qgr, bin_size = bin_size, alpha = .3, bin_method = "mean", bigwig_title = "RUNX1 FE", bigwig_colors = c("red", "green", "blue")),
  #            plot_from_data_source(ds_bw3c, qgr, bin_size = bin_size, alpha = .3, bin_method = "mean", bigwig_title = "CTCF FE", bigwig_colors = c("red", "green", "blue")),
  #            plot_from_data_source(ds_ref, qgr, strand_colors = c("darkslateblue", "lightslateblue"), label_method = "text", top_spacer = 1)
  #          ))
  #   
  # })
  
  # output$plot1 <- renderPlot({
  #   
  #   my_plots = current_style()
  #   main = "Testing"
  #   ggplot_list(my_plots = my_plots, "tall", heights = c(3,2,2,2,2))
  # })
  
}