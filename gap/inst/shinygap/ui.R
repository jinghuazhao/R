source("R/global.R")

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
        h2("Family-based design"),
        fluidRow(
          h3("Parameters"),
          sidebarLayout(
               sidebarPanel(
                   radioButtons("fb_var", "Variable (x axis of the plot):",
                                c("Genotype relative risk (gamma)"="fb_gamma",
                                  "frequency of disease allele (p)"="fb_p",
                                  "type I error (alpha)"="fb_alpha",
                                  "type II error (beta)"="fb_beta"
                                )
                   ),
                   sliderInput("fb_gamma", "Gamma:", min = 1.2, max = 30, value = 4),
                   sliderInput("fb_p", "p:", min = 0.001, max = 0.8, value = 0.01),
                   sliderInput("fb_alpha", "Alpha:", min = 5e-8, max = 0.05, value = 1e-6),
                   sliderInput("fb_beta", "Beta:", min = 0.01, max = 0.4, value = 0.2)
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
        h2("Population-based design"),
        fluidRow(
          h3("Parameters"),
          sidebarLayout(
               sidebarPanel(
                   radioButtons("pb_var", "Variable (x axis of the plot):",
                                c("Prevalence of disease (K)"="pb_kp",
                                  "Genotype relative risk (gamma)"="pb_gamma",
                                  "Frequency of disease allele (p)"="pb_p",
                                  "type I error (alpha)"="pb_alpha",
                                  "type II error (beta)"="pb_beta"
                                 )
                   ),
                   sliderInput("pb_kp", "Kp:", min = 1e-5, max = 0.4, value = 0.1),
                   sliderInput("pb_gamma", "Gamma:", min = 1, max = 30, value = 4.5),
                   sliderInput("pb_p", "p:", min = 1e-3, max = 0.8, value = 0.1),
                   sliderInput("pb_alpha", "Alpha:", min = 5e-8, max = 0.05, value = 5e-8),
                   sliderInput("pb_beta", "Beta:", min = 0.01, max = 0.4, value = 0.2)
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
        h2("Case-cohort design"),
        fluidRow(
          h3("Parameters"),
          sidebarLayout(
               sidebarPanel(
                 checkboxInput("cc_power", "Power/Sample size:", TRUE),
                 radioButtons("cc_var", "Variable (x axis of the plot):",
                              c("Cohort size (n)"="cc_n",
                                "Fraction for subcohort (q)"="cc_q",
                                "Proportion of failure (pD)"="cc_pD",
                                "Proportion of group 1 (p1)"="cc_p1",
                                "log(HR) for two groups (theta)"="cc_theta",
                                "type I error (alpha)"="cc_alpha",
                                "type II error (beta)"="cc_beta"
                               )
                 ),
                 sliderInput("cc_n", "n:", min = 100, max = 500000, value = 25000),
                 sliderInput("cc_q", "q:", min = 0.01, max = 0.5, value = 0.01),
                 sliderInput("cc_pD", "pD:", min = 1e-5, max = 0.8, value = 0.05),
                 sliderInput("cc_p1", "p1:", min = 1e-5, max = 0.8, value = 0.05),
                 sliderInput("cc_theta", "log(HR):", min = 0.01, max = 2.3, value = 0.3),
                 sliderInput("cc_alpha", "Alpha:", min = 5e-8, max = 0.05, value = 0.05),
                 sliderInput("cc_beta", "Beta:", min = 0.01, max = 0.4, value = 0.2)
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
