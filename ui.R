library(TnT)

shinyUI(
  fluidPage(
    tags$h1('Testing shiny browser'),
    tags$h3("Uses ", tags$a("TnT", href = "http://tnt.marlin.pub/", target="_blank")),
    tags$h3("Takes a bit..."),
    sidebarLayout(
      sidebarPanel(
        h1("bigwigs")
        # checkboxInput("CheckFillLines", label = "Fill Lines", value = F)
      ),
      mainPanel(
        TnTOutput("TnTTest_bw")
      )
    ),
    sidebarLayout(
      sidebarPanel(
        h1("hic")
        # checkboxInput("CheckFillLines", label = "Fill Lines", value = F)
      ),
      mainPanel(
        TnTOutput("TnTTest_hic")
      )
    )

    # sidebarPanel(
    #   # selectInput('style', 'Style', c("a", "b", "c"))#,
    # ),
    # mainPanel(
    #   # plotOutput('plot1')
    # )
  )
)
