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
      menuItem("ccDesign", tabName = "ccDesign", icon = icon("th"))
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
          sidebarLayout(
               sidebarPanel(
                   sliderInput("fb_gamma", "Gamma:", min = 1, max = 100, value = 30),
                   sliderInput("fb_p", "p:", min = 0, max = 1, value = 0.5),
                   sliderInput("fb_alpha", "Alpha:", min = 0, max = 1, value = 1e-4),
                   sliderInput("fb_beta", "Beta:", min = 0, max = 1, value = 0.2)
               ),
               mainPanel(
                   h3(verbatimTextOutput("fb_caption")),
                   plotOutput("fb"),
                   radioButtons('fb_reportFormat', 'Report document format:', c('PDF', 'HTML', 'Word'), inline = TRUE),
                   downloadButton("fb_report", "Download report")
               )
          )
        )
      ),
      tabItem(tabName = "pbDesign",
        fluidRow(
          h2("Population-based study design"),
          sidebarLayout(
               sidebarPanel(
                   sliderInput("pb_kp", "Kp:", min = 0, max = 1, value = 0.1),
                   sliderInput("pb_gamma", "Gamma:", min = 1, max = 100, value = 4.5),
                   sliderInput("pb_p", "p:", min = 0, max = 1, value = 0.15),
                   sliderInput("pb_alpha", "Alpha:", min = 0, max = 1, value = 5e-8),
                   sliderInput("pb_beta", "Beta:", min = 0, max = 1, value = 0.2)
               ),
               mainPanel(
                   h3(verbatimTextOutput("pb_caption")),
                   plotOutput("pb"),
                   radioButtons('pb_reportFormat', 'Report document format:', c('PDF', 'HTML', 'Word'), inline = TRUE),
                   downloadButton("pb_report", "Download report")
               )
          )
        )
      ),
      tabItem(tabName = "ccDesign",
        fluidRow(
          h2("Case-cohort study design"),
          sidebarLayout(
               sidebarPanel(
                 sliderInput("cc_n", "n:", min = 0, max = 1, value = 1),
                 sliderInput("cc_q", "q:", min = 0, max = 1, value = 0.5),
                 sliderInput("cc_pD", "pD:", min = 0, max = 1, value = 0.15),
                 sliderInput("cc_p1", "p1:", min = 0, max = 1, value = 0.15),
                 sliderInput("cc_alpha", "Alpha:", min = 0, max = 1, value = 5e-8),
                 sliderInput("cc_theta", "theta:", min = 0, max = 1, value = 0.2),
                 sliderInput("cc_power", "power:", min = 0, max = 1, value = 0.2)
               ),
               mainPanel(
                   h3(verbatimTextOutput("cc_caption")),
                   plotOutput("cc"),
                   radioButtons('cc_reportFormat', 'Report document format:', c('PDF', 'HTML', 'Word'), inline = TRUE),
                   downloadButton("cc_report", "Download report")
               )
          ) # siderbarLayout
        ) # fluidRow
      ) # tabitem
    ) # tabItems
  ) # dashbarBody
)
