library(TnT)

shinyUI(
  fluidPage(
    tags$h1('Testing shiny browser'),
    tags$h3("Uses ", tags$a("TnT", href = "http://tnt.marlin.pub/", target="_blank")),
    TnTOutput("TnTTest")
    # sidebarPanel(
    #   # selectInput('style', 'Style', c("a", "b", "c"))#,
    # ),
    # mainPanel(
    #   # plotOutput('plot1')
    # )
  )
)
