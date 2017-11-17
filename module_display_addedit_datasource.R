###file chooser setups

module_display_addedit_datasource = function(input, output, session, rSources, rDisps){
  
  output$SelectDataSource = renderUI({
    
    selectInput("SelDs", label = "Select DataSource", choices = rSources()$file_name)
    
  })
  
  #available types:
  # bigwig
  # hic_arc
  # hic_insulation
  # hic_matrix
  # ref
  display_types = list("bed" = "NYI", "txt" = 1, "Peak" = 2, "hic" = 3, "bigwig" = 1, "bw" = 2)
  
  last_selected = "no data added"
  
  observeEvent(
    {
      input$SelDs
    }, {
      if(is.null(input$SelIDs)){
        return()
      }
      showNotification("maybe update")
      if(input$SelIDs != last_selected){
        return()
      }
      showNotification("do update")
      
      print(getSelectedDataSource())
    })
  
  # getSelectedDisp = function(){
  #   subset(rSources(), file_name == input$SelDs)
  # }
  
  addDisplay = function(plot_id, plot_type, plot_title = "", relative_height = 1, hidden_ds_object){
    ds = rDisps()
    ds = ds[input$DT_Displays_rows_all,]
    ds = ds[ds$plot_id != "no data added",]
    new_ds = data.frame(
      plot_id = plot_id,
      plot_type = "no data added",
      plot_title = plot_type,
      relative_height = relative_height, 
      hidden_ds_object = hidden_ds_object,
      stringsAsFactors = F 
    )
    if(any(ds$plot_id == new_ds$plot_id)){
      ds[ds$plot_id == new_ds$plot_id, ] = new_ds
    }else{
      ds = rbind(ds, new_ds)  
    }
    
    rDisps(ds)
  }
  
  deleteDisplay = function(){
    ds = rDisps()
    rDisps(ds[-input$DT_Displays_rows_selected,])
  }
  
  observeEvent(input$BtnAddNewDisplay, {
    if(is.null(input$DT_Displays_rows_selected)){
      showNotification("")
    }
    
    ds = rSources()[input$DT_Displays_rows_selected,]
    save(ds, file = "ds.save")
    showModal(addDisplayModal(
      plot_id = paste(sample(LETTERS, 1), sample(letters, 1)),
      plot_type = "1", 
      plot_title = "title", 
      relative_height = 1
    ))
  })
  observeEvent(input$BtnCancelDisplayModal, {
    removeModal()
  })
  observeEvent(input$BtnConfirmDisplayModal, {
    addDisplay(plot_id = paste(sample(LETTERS, 1), sample(letters, 1)), 
               file_type = input$RadioDisplayTypes, 
               file_path = input$textFilepath, 
               relative_height = 1)
    removeModal()
  })
  
  observeEvent(input$BtnDeleteDisplay, {
    if(is.null(input$DT_Displays_rows_selected)){
      showNotification("Select some data to delete.", type = "error")
      return()
    }
    showNotification(paste("delete", as.character(input$DT_Displays_rows_selected)))
    deleteDisplay()
  })
  
  observeEvent(input$BtnEditDisplay, {
    if(is.null(input$DT_Displays_rows_selected)){
      showNotification("Select some data to edit", type = "error")
      return()
    }
    showNotification(paste("edit", as.character(input$DT_Displays_rows_selected)))
    disp = rDisps()[input$DT_Displays_rows_selected,]
    showModal(addDisplayModal(plot_id = disp$plot_id, 
                              plot_type = disp$plot_type, 
                              plot_title = disp$plot_title, 
                              relative_height = disp$relative_height))
  })
  
  output$DT_Displays = DT::renderDataTable({
    disp = rDisps()
    disp = disp[,!grepl("hidden", colnames(disp))]
    DT::datatable(disp, 
                  rownames = F, 
                  selection = "single", 
                  extensions = 'Buttons',
                  options = list(
                    dom = 'Bfrtip',
                    buttons = c('csv', 'excel')
                  )
    )
  })
}

addDisplayModal <- function(plot_id, plot_type = "not_specified", plot_title = "", relative_height = 1, failed = FALSE) {
  # plot_id = plot_id,
  # plot_type = "no data added",
  # plot_title = plot_type,
  # relative_height = relative_height, 
  # hidden_ds_object = hidden_ds_object,
  # if(file_name == "not_specified"){
  #   file_name = basename(file_path)
  # }
  # type_choices = c("bed", "bigwig", "hic sparse matrix",  "juicebox hic")
  
  # type = switch(file_suffix,
  #               bed = {
  #                 "bed"
  #               },
  #               txt = {
  #                 "hic sparse matrix"
  #               },
  #               bigwig = {
  #                 "bigwig"
  #               },
  #               bw = {
  #                 "bigwig"
  #               }, 
  #               hic = {
  #                 "juicebox hic"
  #               },
  #               "bed"
  # )
  modalDialog(
    if (failed)
      div(tags$b("Hey man, this didn't work! Not sure why.", style = "color: red;")),
    span(paste('Please name sample for file:', plot_id)),
    hidden(textInput("textPlotID", label = "shouldn't see", value = plot_id)),
    # textInput("textFilename", label = "File Name", value = file_name),
    radioButtons(inputId = "RadioDisplayTypes", label = "Display Type", choices = type_choices, selected = type),
    footer = tagList(fluidRow(
      column(width = 8, tags$hr()
      ),
      column(width = 4,
             actionButton("BtnCancelDisplayModal", "Cancel"),
             actionButton("BtnConfirmDisplayModal", "Confirm")
      )
    )
    ),
    size = "m",
    title = "Display Setup"
  )
}

module_display_addedit_datasource_ui = function(){
  fluidRow(
    column(width = 4, 
           h5("Setup Data Sources"),
           actionButton(inputId = "BtnAddNewDisplay", label = "Add New"),
           actionButton(inputId = "BtnDeleteDisplayFile", "Delete"),
           actionButton(inputId = "BtnEditDisplayFile", "Edit")
    ),
    column(width = 8, 
           DT::dataTableOutput("DT_Displays")
    )
  )
}