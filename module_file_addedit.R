###file chooser setups

module_file_addedit = function(input, output, session, rDF){
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
  
  shinyFileChoose(input, 'FilesDataSource', 
                  roots= roots_load_set, 
                  filetypes=c("bed", "txt", "Peak", "hic", "bigwig", "bw"))
  
  addDataSource = function(file_name, file_type, file_path){
    ds = rDF()
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
    
    rDF(ds)
  }
  
  deleteDataSource = function(){
    ds = rDF()
    rDF(ds[-input$DT_DataSources_rows_selected,])
  }
  
  observeEvent(input$FilesDataSource, {
    print("server find file")
    # file_path = shinyFiles2load(input$FilesLoadSet, roots_load_set)
    file_path = as.character(parseFilePaths(roots_load_set, input$FilesDataSource)$datapath)
    showNotification(paste(basename(file_path), "- what kind of file is this?"))
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
    if(is.null(input$DT_DataSources_rows_selected)){
      showNotification("Select some data to delete.", type = "error")
      return()
    }
    showNotification(paste("delete", as.character(input$DT_DataSources_rows_selected)))
    deleteDataSource()
  })
  
  observeEvent(input$BtnEditDataSource, {
    if(is.null(input$DT_DataSources_rows_selected)){
      showNotification("Select some data to edit", type = "error")
      return()
    }
    showNotification(paste("edit", as.character(input$DT_DataSources_rows_selected)))
    ds = rDF()[input$DT_DataSources_rows_selected,]
    showModal(addFileModal(file_path = ds$file_path, file_name = ds$file_name, file_type = ds$file_type))
  })
  
  output$DT_DataSources = DT::renderDataTable({
    # df = data.frame(file_path = "no data added", file_type = "no data added")
    DT::datatable(rDF(), rownames = F, selection = "single")
  })
}

addFileModal <- function(file_path, file_name = "not_specified", file_type = "not_specified", failed = FALSE) {
  if(file_name == "not_specified"){
    file_name = basename(file_path)
  }
  type_choices = c("bed", "bigwig", "hic sparse matrix",  "juicebox hic")
  get_suffix = function(f){
    ext = strsplit(file_path, split = "\\.")[[1]]
    ext[length(ext)]
  }
  file_suffix = get_suffix(file_path)
  type = switch(file_suffix,
                bed = {
                  "bed"
                },
                txt = {
                  "hic sparse matrix"
                },
                bigwig = {
                  "bigwig"
                },
                bw = {
                  "bigwig"
                }, 
                hic = {
                  "juicebox hic"
                },
                "bed"
  )
  modalDialog(
    if (failed)
      div(tags$b("Hey man, this didn't work! Not sure why.", style = "color: red;")),
    span(paste('Please name sample for file:', file_path)),
    hidden(textInput("textFilepath", label = "shouldn't see", value = file_path)),
    textInput("textFilename", label = "File Name", value = file_name),
    radioButtons(inputId = "RadioFileTypes", label = "File Type", choices = type_choices, selected = type),
    footer = tagList(fluidRow(
      column(width = 8, tags$hr()
      ),
      column(width = 4,
             actionButton("BtnCancelFileModal", "Cancel"),
             actionButton("BtnConfirmFileModal", "Confirm")
      )
    )
    ),
    size = "m",
    title = "File Setup"
  )
}

