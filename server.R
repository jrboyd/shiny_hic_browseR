

source("module_file_addedit.R")

function(input, output, session) {
  rDataSources = reactiveVal(data.frame(
    file_name = "no data added",
    file_type = "no data added",
    file_path = "no data added", 
    stringsAsFactors = F 
  ))
  
  module_file_addedit(input, output, session, rDataSources)
  
  
}