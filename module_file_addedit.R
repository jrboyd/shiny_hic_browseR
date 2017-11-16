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