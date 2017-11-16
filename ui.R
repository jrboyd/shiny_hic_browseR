library(shinythemes)
library(shiny)
library(DT)
library(shinyFiles)
library(shinyjs)
shinyUI(
  
  fluidPage(
    useShinyjs(),
    theme = shinytheme(theme = "spacelab"),
    tags$h3("HiC browseR"),
    tags$hr(),
    fluidRow(
      column(width = 4, 
             h5("Setup Data Sources"),
             shinyFilesButton(id = "FilesDataSource", 
                              label = "Find Files on Server", 
                              title = "Find Data Source Files", multiple = T),
             actionButton(inputId = "BtnDeleteDataSource", "Delete"),
             actionButton(inputId = "BtnEditDataSource", "Edit")
      ),
      column(width = 8, 
             DT::dataTableOutput("DT_DataSources")
      )
    ),
    tags$hr()
  )
)
