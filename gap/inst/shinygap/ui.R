source("R/init.R")

ui <- dashboardPage(
  title="shinygap",
  skin="purple",
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
      menuItem("ccDesign", tabName = "ccDesign", icon = icon("th"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "landing",
         div(class = "jumbotron", HTML("<center><h1>Shiny for genetic analysis package (gap) designs</h1></center>")),
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
          sidebarLayout(
               sidebarPanel(
                   uiOutput("fb_var"),
                   sliderInput("fb_n", "n:", min = 1, max = 1000000, value = 500),
                   sliderInput("fb_gamma", "Gamma:", min = 1, max = 100, value = 1.5),
                   sliderInput("fb_p", "p:", min = 0, max = 1, value = 0.2),
                   sliderInput("fb_alpha", "Alpha:", min = 0, max = 1, value = 1e-4),
                   sliderInput("fb_beta", "Beta:", min = 0, max = 1, value = 0.2)
               ),
               mainPanel(
                   h3(verbatimTextOutput("fb_caption")),
                   plotlyOutput("fb"),
                   h3("Header of data:"),
                   tableOutput("fb_preview"),
                   radioButtons('fb_downloadFormat', 'Download file format:', c('bz2', 'gz', 'tsv', 'xz'), inline = TRUE),
                   downloadButton("fb_download", "Download data")
               )
          )
        )
      ),
      tabItem(tabName = "pbDesign",
        h2("Population-based study design"),
        fluidRow(
          h3("Parameters"),
          sidebarLayout(
               sidebarPanel(
                   uiOutput("pb_var"),
                   sliderInput("pb_n", "n:", min = 1, max = 1000000, value = 500),
                   sliderInput("pb_kp", "Kp:", min = 0, max = 1, value = 0.1),
                   sliderInput("pb_gamma", "Gamma:", min = 1, max = 100, value = 4.5),
                   sliderInput("pb_p", "p:", min = 0, max = 1, value = 0.15),
                   sliderInput("pb_alpha", "Alpha:", min = 0, max = 1, value = 5e-8),
                   sliderInput("pb_beta", "Beta:", min = 0, max = 1, value = 0.2)
               ),
               mainPanel(
                   h3(verbatimTextOutput("pb_caption")),
                   plotlyOutput("pb"),
                   h3("Header of data:"),
                   tableOutput("pb_preview"),
                   radioButtons('pb_downloadFormat', 'Download file format:', c('bz2', 'bz', 'tsv', 'xz'), inline = TRUE),
                   downloadButton("pb_download", "Download data")
               )
          )
        )
      ),
      tabItem(tabName = "ccDesign",
        h2("Case-cohort study design"),
        fluidRow(
          h3("Parameters"),
          sidebarLayout(
               sidebarPanel(
                 checkboxInput("cc_power", "sample size/power:", TRUE),
                 uiOutput("cc_var"),
                 sliderInput("cc_n", "n:", min = 0, max = 1, value = 1),
                 sliderInput("cc_q", "q:", min = 0, max = 1, value = 0.5),
                 sliderInput("cc_pD", "pD:", min = 0, max = 1, value = 0.15),
                 sliderInput("cc_p1", "p1:", min = 0, max = 1, value = 0.15),
                 sliderInput("cc_alpha", "Alpha:", min = 0, max = 1, value = 5e-8),
                 sliderInput("cc_hr", "hr:", min = 1, max = 100, value = 1.2)
               ),
               mainPanel(
                   h3(verbatimTextOutput("cc_caption")),
                   plotlyOutput("cc"),
                   h3("Header of data:"),
                   tableOutput("cc_preview"),
                   radioButtons('cc_downloadFormat', 'Download file format:', c('bz2', 'gz', 'tsv', 'xz'), inline = TRUE),
                   downloadButton("cc_download", "Download data")
               )
          ) # siderbarLayout
        ) # fluidRow
      ) # tabitem
    ) # tabItems
  ) # dashbarBody
)
