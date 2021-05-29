source("R/init.R")

ui <- dashboardPage(
  title = "shinygap",
  dashboardHeader(
     title = tags$a(href='https://github.com/jinghuazhao/', target = '_blank',
                    tags$img(src=paste0("bees.svg"), height = "70%", width = "auto", align = "middle")
     ),
     dropdownMenu(type = "messages",
        tags$li(HTML('<li><a href="https://github.com/jinghuazhao/ShinyApps" target="_blank"><i class="fa fa-code-branch"></i><h4>GitHub</h4></a></li>')),
        tags$li(HTML('<li><a href="mailto:jinghuazhao@hotmail.com" target="_blank"><i class="fa fa-question"></i><h4>email</h4></a></li>'))
     )
  ),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Home", tabName = "landing", icon = icon("home")),
      menuItem("fbDesign", tabName = "fbDesign", icon = icon("upload")),
      menuItem("pbDesign", tabName = "pbDesign", icon = icon("download")),
      menuItem("ccDesign", tabName = "ccDesign", icon = icon("th")),
      menuItem("Report", tabName = "Report", icon = icon("book-open"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "landing",
         div(class = "jumbotron", HTML("<center><h1>Shiny for gap designs</h1></center>")),
         fluidRow(div(class = "col-sm-12",
                      div(class = "box box-primary", style = "padding-right: 5%; padding-left: 5%; font-size:100%", NULL,
                          div(class = "box-body", shiny::includeMarkdown("README.md"))
                      )
                 ),
         withMathJax()
         )
      ),
      tabItem(tabName = "fbDesign",
        h2("Family-based study design"),
        fluidRow(
          h3("Parameters"),
          sliderInput("gamma", "Gamma:", min = 1, max = 50, value = 30),
          sliderInput("p", "p:", min = 0, max = 1, value = 30),
          sliderInput("alpha", "Alpha:", min = 0, max = 1, value = 1e-4),
          sliderInput("beta", "Beta:", min = 0, max = 1, value = 0.2),
          checkboxInput("fb", "Select fb model", FALSE)
        )
      ),
      tabItem(tabName = "pbDesign",
        fluidRow(
          h2("Population-based study design"),
          sliderInput("kp", "Kp:", min = 0, max = 1, value = 0.1),
          sliderInput("gamma", "Gamma:", min = 0, max = 1, value = 4.5),
          sliderInput("p", "p:", min = 0, max = 1, value = 0.15),
          sliderInput("alpha", "Alpha:", min = 0, max = 1, value = 5e-8),
          sliderInput("beta", "Beta:", min = 0, max = 1, value = 0.2),
          checkboxInput("pb", "Select pb model", FALSE)
        )
      ),
      tabItem(tabName = "ccDesign",
        fluidRow(
          h2("Case-cohort study design"),
          sliderInput("n", "n:", min = 0, max = 1, value = 1),
          sliderInput("q", "q:", min = 0, max = 1, value = 0.5),
          sliderInput("pD", "pD:", min = 0, max = 1, value = 0.15),
          sliderInput("p1", "p1:", min = 0, max = 1, value = 0.15),
          sliderInput("alpha", "Alpha:", min = 0, max = 1, value = 5e-8),
          sliderInput("theta", "theta:", min = 0, max = 1, value = 0.2),
          sliderInput("power", "power:", min = 0, max = 1, value = 0.2),
          checkboxInput("cc", "Select cc model", FALSE)
        )
      ),
      tabItem(tabName = "Report",
        helpText("Results of the model."),
        h3(textOutput("fb_caption")),
        plotOutput("fb"),
        h3(textOutput("pb_caption")),
        plotOutput("pb"),
        h3(textOutput("cc_caption")),
        plotOutput("cc"),
        h3("Model summary"),
        verbatimTextOutput("summary"),
        h3("Model fit summary"),
        verbatimTextOutput("fit"),
        radioButtons('reportFormat', 'Report document format:', c('PDF', 'HTML', 'Word'), inline = TRUE),
        downloadButton("report", "Download report")
      )
    )
  )
)
