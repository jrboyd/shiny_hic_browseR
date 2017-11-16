
require(magrittr)
bed_path = "~/ShinyApps/shiny_peak_data/beds/"
names(bed_path) = "intersectR"
user_roots = dir("/slipstream/home/", full.names = T) %>% paste0(. , "/ShinyData")
user_roots = subset(user_roots, dir.exists(user_roots))
names(user_roots) = dirname(user_roots) %>% basename()
qcframework_load <<- dir("/slipstream/galaxy/uploads/working/qc_framework", pattern = "^output", full.names = T)
names(qcframework_load) <- basename(qcframework_load)
ucsc_load <<- "/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/"
names(ucsc_load) <- basename(ucsc_load)
hic_load <<- "/slipstream/home/joeboyd/HiC-Pro/shiny_data"
names(hic_load) <- "HiC"
roots_load_set = c(hic_load, bed_path, ucsc_load, user_roots, qcframework_load)

source("module_file_addedit.R")

function(input, output, session) {
  ###file chooser setups
  shinyFileChoose(input, 'FilesDataSource', 
                  roots= roots_load_set, 
                  filetypes=c("bed", "txt", "Peak", "hic", "bigwig", "bw"))
  
  rDataSources = reactiveVal(data.frame(
    file_name = "no data added",
    file_type = "no data added",
    file_path = "no data added", 
    stringsAsFactors = F 
  ))
  
  addDataSource = function(file_name, file_type, file_path){
    ds = rDataSources()
    ds = ds[input$DT_DataSources_rows_all,]
    ds = ds[ds$file_path != "no data added",]
    new_ds = data.frame(
      file_name = file_name,
      file_type = file_type,
      file_path = file_path, 
      stringsAsFactors = F 
    )
    if(any(ds$file_path == new_ds$file_path)){
      ds[ds$file_path == new_ds$file_path, ] = new_ds
    }else{
      ds = rbind(ds, new_ds)  
    }
    
    rDataSources(ds)
  }

    deleteDataSource = function(){
    ds = rDataSources()
    rDataSources(ds[-input$DT_DataSources_rows_selected,])
  }
  
  observeEvent(input$FilesDataSource, {
    print("server find file")
    # file_path = shinyFiles2load(input$FilesLoadSet, roots_load_set)
    file_path = as.character(parseFilePaths(roots_load_set, input$FilesDataSource)$datapath)
    showNotification(paste(basename(file_path), "- what kind of file is this?"))
    # showNotification(format(file.size(file_path, ), units = "MB"))
    showModal(addFileModal(file_path))
  })
  observeEvent(input$BtnCancelFileModal, {
    removeModal()
  })
  observeEvent(input$BtnConfirmFileModal, {
    addDataSource(file_name = input$textFilename, 
                  file_type = input$RadioFileTypes, 
                  file_path = input$textFilepath)
    removeModal()
  })
  
  observeEvent(input$BtnDeleteDataSource, {
    deleteDataSource()
  })
  
  observeEvent(input$BtnEditDataSource, {
    ds = rDataSources()[input$DT_DataSources_rows_selected,]
    showModal(addFileModal(file_path = ds$file_path, file_name = ds$file_name, file_type = ds$file_type))
  })
  
  output$DT_DataSources = DT::renderDataTable({
    # df = data.frame(file_path = "no data added", file_type = "no data added")
    DT::datatable(rDataSources(), rownames = F, selection = "single")
  })
}