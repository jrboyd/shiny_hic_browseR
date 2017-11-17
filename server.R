
empty_ds = data.frame(
  file_name = "no data added",
  file_type = "no data added",
  file_path = "no data added", 
  stringsAsFactors = F 
)

default_ds = data.frame(
  file_name = c("MCF10A", "MCF10AT1", "MCF10CA1a"),
  file_type = c("juicebox hic", "juicebox hic", "juicebox hic"),
  file_path = c("/slipstream/home/joeboyd/HiC-Pro/shiny_data/MCF10A_pooled_HPr_allValidPairs.hic", 
                "/slipstream/home/joeboyd/HiC-Pro/shiny_data/MCF10AT1_pooled_HPr_allValidPairs.hic",
                "/slipstream/home/joeboyd/HiC-Pro/shiny_data/MCF10CA1a_pooled_HPr_allValidPairs.hic"),
  stringsAsFactors = F
)

touse_ds = default_ds

function(input, output, session) {
  rDataSources = reactiveVal(touse_ds)
  
  module_file_addedit_datasource(input, output, session, rDataSources)
  
  rDisplays = reactiveVal(data.frame(
    plot_id = "no data added",
    plot_type = "no data added",
    plot_title = "no data added",
    relative_height = "no data added", 
    hidden_ds_object = "no data added",
    stringsAsFactors = F 
  ))
  
  module_display_addedit_datasource(input, output, session, rDataSources, rDisplays)
}