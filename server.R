library(shiny)
library(TnT)
library(GenomicRanges)
library(rtracklayer)
library(Rcpp)
library(data.table)
sourceCpp("straw-R.cpp")
hic_files = dir("/slipstream/home/stein-hic/juicebox_files/hic-pro_resolutions/", pattern = ".hic$", full.names = T)
names(hic_files) = basename(hic_files)
hic_files = as.list(hic_files)

fetch_hic = function(hic_f, chr, s, e, chr2 = NULL, s2 = NULL, e2 = NULL, fill_matrix = F, res = HICPRO_RESOLUTIONS[3]){
  base_cmd = "VC SUB_HIC_FILE SUB_POS_A SUB_POS_B BP SUB_RES"
  pos_a = paste(c(chr, s, e), collapse = ":")
  if(is.null(chr2)) chr2 = chr
  if(is.null(s2)) s2 = s
  if(is.null(e2)) e2 = e
  pos_b = paste(c(chr2, s2, e2), collapse = ":")
  cmd = base_cmd
  cmd = sub("SUB_HIC_FILE", hic_f, cmd)
  cmd = sub("SUB_POS_A", pos_a, cmd)
  cmd = sub("SUB_POS_B", pos_b, cmd)
  cmd = sub("SUB_RES", res, cmd)
  hic_dt <- as.data.table(straw_R(cmd))
  if(fill_matrix){
    hic_dt = rbind(hic_dt, hic_dt[x != y, .(x = y, y = x, counts)])
  }
  return(hic_dt)
}
fetch_hic.gr = function(hic_f, qgr, qgr2 = NULL, fill_matrix = F, res = HICPRO_RESOLUTIONS[3]){
  if(is.null(qgr2)) qgr2 = qgr
  fetch_hic(
    hic_f = hic_f, 
    chr = as.character(seqnames(qgr))[1],
    s = start(qgr)[1],
    e = end(qgr)[1],
    chr2 = as.character(seqnames(qgr2))[1],
    s2 = start(qgr2)[1],
    e2 = end(qgr2)[1],
    fill_matrix = fill_matrix,
    res = res
  )
}


### hic
# qgr = range(subset(gene, symbol == "RUNX1"))
# vgr = qgr
# start(qgr) = start(qgr) - 5*10^6
# end(qgr) = end(qgr) + 5*10^6
# # qgr = GRanges("chr6", IRanges(23*10^6, 28.2*10^6))
# HICPRO_RESOLUTIONS = c(20000, 40000, 150000, 500000, 1000000)
# res = HICPRO_RESOLUTIONS[2]
# dt = fetch_hic.gr(hic_files$MCF10CA1a_pooled_HPr_allValidPairs.hic, qgr, res = res)
# library(ggplot2)
# # ggplot(dt) + geom_tile(aes(x = x, y = y, fill = counts))
# dt[, dist := (y - x) / res]
# dt = dt[dist < 50]
# dt[, xmid := x + res / 2]
# dt[, ymid := y + res / 2]
# dt[, mid := (xmid + ymid) / 2]
# dt[, scaled := counts / max(counts)]
# ggplot(dt) + 
#   geom_point(aes(x = mid, y = dist, col = counts),  size = 1, shape = 16) + 
#   scale_color_gradient(low = "white", high = "red") + theme_minimal()
# hic_colors = rgb(colorRamp(c("white", "red"))(dt$scaled)/255)
# hic_gr = GRanges(seqnames = as.character(seqnames(qgr)), ranges = IRanges(dt$mid, dt$mid), height = dt$dist, colors = hic_colors)
# pt = PinTrack(hic_gr, value = hic_gr$height, color = hic_colors, height = 300) 

### reference
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(EnsDb.Hsapiens.v86)
# txdb =  TxDb.Hsapiens.UCSC.hg38.knownGene
# gtrack <- TnT::GeneTrackFromTxDb(txdb)
gene <- genes(EnsDb.Hsapiens.v86)
gene_fix = GRanges(paste0("chr", seqnames(gene)), IRanges(start(gene), end(gene)))
elementMetadata(gene_fix) = elementMetadata(gene)
gene = gene_fix
ensGeneTrack <- TnT::FeatureTrack(gene, tooltip = as.data.frame(gene),
                                  names = paste(gene$symbol, " (", gene$gene_biotype, ")", sep = ""),
                                  color = TnT::mapcol(gene$gene_biotype, palette.fun = grDevices::rainbow))
# TnTGenome(list(pt, gtrack, ensGeneTrack), view.range = vgr, coord.range = ranges(qgr))
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
bedgraph2lineTrack = function(bdg_gr, win = 50, label = "Track Name", color = "black", floor_at = 0, fill = F ){
  bdg_gr = subset(bdg_gr, score > floor_at)
  line_gr = bedgraph2line(bdg_gr = bdg_gr, win = win)
  line_gr$score = line_gr$score - floor_at
  if(fill){
    AreaTrack(line_gr, value = line_gr$score, color = color, label = label)
  }else{
    LineTrack(line_gr, value = line_gr$score, color = color, label = label)  
  }
  
}
hic2pinTrack = function(hic_file, qgr, label = "Track Name", color = "black", res = 40000, fillmax = NULL){
  qgr = range(qgr)
  HICPRO_RESOLUTIONS = c(20000, 40000, 150000, 500000, 1000000)
  # res = HICPRO_RESOLUTIONS[2]
  dt = fetch_hic.gr(hic_file, qgr, res = res)
  library(ggplot2)
  # ggplot(dt) + geom_tile(aes(x = x, y = y, fill = counts))
  dt[, dist := (y - x) / res]
  dt = dt[dist < 50]
  dt[, xmid := x + res / 2]
  dt[, ymid := y + res / 2]
  dt[, mid := (xmid + ymid) / 2]
  if(!is.null(fillmax)) dt$counts = ifelse(dt$counts > fillmax, fillmax, dt$counts)
  dt[, scaled := counts / max(counts)]
  ggplot(dt) + 
    geom_point(aes(x = mid, y = dist, col = counts),  size = 1, shape = 16) + 
    scale_color_gradient(low = "white", high = "red") + theme_minimal()
  hic_colors = rgb(colorRamp(c("white", color))(dt$scaled)/255)
  hic_gr = GRanges(seqnames = as.character(seqnames(qgr)), ranges = IRanges(dt$mid, dt$mid), height = dt$dist, colors = hic_colors)
  PinTrack(hic_gr, value = hic_gr$height, color = hic_colors, height = 300, label = label) 
}

bedgraph2blockTrack = function(bdg_gr, label = "Track Name", color = "black", floor_at = 2, ceiling_at = 6){
  bdg_gr$color = score2color(bdg_gr$score, scale_min = floor_at, colors = c("white", color), scale_max = ceiling_at)
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
  
  output$TnTTest_bw = renderTnT({
    qgr = GRanges("chr6", IRanges(23.1*10^6, 27.5*10^6))
    #import using rtracklayer
    runx1_gr = import.bw(bw_files$MDA231_Runx1_pooled_FE.bw, which = qgr)
    runx2_gr = import.bw(bw_files$MDA231_Runx2_pooled_FE.bw, which = qgr)
    #block tracks
    bt_runx1 = bedgraph2blockTrack(runx1_gr, label = "Runx1 FE", color = "green")
    bt_runx2 = bedgraph2blockTrack(runx2_gr, label = "Runx2 FE", color = "purple")
    #line tracks
    lt_runx1 = bedgraph2lineTrack(runx1_gr, label = "Runx1 FE", color = "green", win = 50, floor_at = 1, fill = T)
    lt_runx2 = bedgraph2lineTrack(runx2_gr, label = "Runx2 FE", color = "purple", win = 50, floor_at = 1, fill = T)
    #peak tracks 
    np_runx1 = peak2blockTrack(np_files$MDA231_Runx1_pooled_peaks.narrowPeak, label = "Runx1 Peaks", color = "green")
    np_runx2 = peak2blockTrack(np_files$MDA231_Runx2_pooled_peaks.narrowPeak, label = "Runx2 Peaks", color = "purple")
    
    vgr = GRanges("chr6", IRanges(25.972*10^6-1800, 26.228*10^6+1800))
    genome(seqinfo(vgr)) = "hg38"
    TnTGenome(tracklist = list(ensGeneTrack, np_runx1, bt_runx1, lt_runx1,  np_runx2, bt_runx2, lt_runx2), 
              zoom.allow = c(1000, 10^6),
              coord.range = ranges(qgr),
              view.range = vgr)
  })
  
  output$TnTTest_hic = renderTnT({
    qgr = GRanges("chr6", IRanges(24.1*10^6, 27.5*10^6))
    #import using rtracklayer
    runx1_gr = import.bw(bw_files$MDA231_Runx1_pooled_FE.bw, which = qgr)
    runx2_gr = import.bw(bw_files$MDA231_Runx2_pooled_FE.bw, which = qgr)
    #block tracks
    bt_runx1 = bedgraph2blockTrack(runx1_gr, label = "Runx1 FE", color = "green")
    bt_runx2 = bedgraph2blockTrack(runx2_gr, label = "Runx2 FE", color = "purple")
    #line tracks
    lt_runx1 = bedgraph2lineTrack(runx1_gr, label = "Runx1 FE", color = "green", win = 500, floor_at = 1, fill = T)
    lt_runx2 = bedgraph2lineTrack(runx2_gr, label = "Runx2 FE", color = "purple", win = 500, floor_at = 1, fill = T)
    #peak tracks 
    np_runx1 = peak2blockTrack(np_files$MDA231_Runx1_pooled_peaks.narrowPeak, label = "Runx1 Peaks", color = "green")
    np_runx2 = peak2blockTrack(np_files$MDA231_Runx2_pooled_peaks.narrowPeak, label = "Runx2 Peaks", color = "purple")
    
    hic_10a = hic2pinTrack(hic_file = hic_files$MCF10A_pooled_HPr_allValidPairs.hic, qgr, label = "MCF10A hic", fillmax = 40, res = 40000)
    hic_at1 = hic2pinTrack(hic_file = hic_files$MCF10AT1_pooled_HPr_allValidPairs.hic, qgr, label = "AT1 hic", fillmax = 40, res = 40000)
    hic_ca1a = hic2pinTrack(hic_file = hic_files$MCF10CA1a_pooled_HPr_allValidPairs.hic, qgr, label = "CA1a hic", fillmax = 40, res = 40000)
    vgr = GRanges("chr6", IRanges(25*10^6, 26.4*10^6))
    genome(seqinfo(vgr)) = "hg38"
    TnTGenome(tracklist = list(ensGeneTrack, np_runx1, bt_runx1, lt_runx1,  np_runx2, bt_runx2, lt_runx2, hic_10a, hic_at1, hic_ca1a), 
              zoom.allow = c(1000, 10^6),
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