library(shinythemes)
library(shiny)
library(DT)
library(shinyFiles)
library(shinyjs)

source("module_file_addedit_datasource.R")
source("module_display_addedit_datasource.R")

shinyUI(
  
  fluidPage(
    useShinyjs(),
    theme = shinytheme(theme = "spacelab"),
    tags$h3("HiC browseR"),
    tags$hr(),
    tags$h4("1) Locate DataSources"),
    module_file_addedit_datasource_ui(),
    tags$hr(),
    tags$h4("2) Configure DataSource Plots"),
    uiOutput(outputId = "SelectDataSource"),
    module_display_addedit_datasource_ui()
  )
)
