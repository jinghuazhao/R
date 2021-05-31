server <- function(input, output) {
# fb curve
  output$fb_var <- renderUI({
     selectInput("fb_var", "Variable (x axis of the plot):",
                 c("Sample size"="N","Genotype relative risk"="gamma","p","type I error"="alpha","type II error"="beta"), selected="N")
  })
  output$fb_caption <- reactive({paste("Figure: family-based design as a function of",input$fb_var)})
  output$fb <- renderPlot({
     g <- input$fb_gamma
     p <- input$fb_p
     alpha <- input$fb_alpha
     beta <- input$fb_beta
     z <- gap::fbsize(g,p)
     plot(z$n3)
  })
  output$fb_report <- downloadHandler(
    filename = function() {
      paste("fb", sep = ".", switch(input$fb_reportFormat, PDF = 'pdf', HTML = 'html', Word = 'docx'))
    },
    content = function(file) {
      src <- normalizePath('report.Rmd')
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      file.copy(src, 'report.Rmd', overwrite = TRUE)
      out <- render('report.Rmd', switch(input$fb_reportFormat, PDF = pdf_document(), HTML = html_document(), Word = word_document()))
      file.rename(out, file)
    }
  )
# pb curve
  output$pb_var <- renderUI({
     selectInput("pb_var", "Variable (x axis of the plot):",
                 c("Sample size"="N","Prevalence of disease"="Kp","Genotype relative risk"="gamma","p","type I error"="alpha","type II error"="beta"),
                 selected="N")
  })
  output$pb_caption <- reactive({paste("Figure: population-based design as a function of", input$pb_var)})
  output$pb <- renderPlot({
     k <- input$pb_kp
     g <- input$pb_gamma
     p <- input$pb_p
     alpha <- input$pb_alpha
     beta <- input$pb_beta
     z <- ceiling(gap::pbsize(k,g,p))
     plot(z, type="b")
  })
  output$pb_report <- downloadHandler(
    filename = function() {
      paste("pb", sep = ".", switch(input$pb_reportFormat, PDF = 'pdf', HTML = 'html', Word = 'docx'))
    },
    content = function(file) {
      src <- normalizePath('report.Rmd')
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      file.copy(src, 'report.Rmd', overwrite = TRUE)
      out <- render('report.Rmd', switch(input$pb_reportFormat, PDF = pdf_document(), HTML = html_document(), Word = word_document()))
      file.rename(out, file)
    }
  )
# cc curve
  output$cc_var <- renderUI({
     selectInput("cc_var", "Variable (x axis of the plot):",
                 c("Cohort size"="n","Fraction for subcohort"="q","Proportion of failure in full cohort"="pD",
                   "proportions of the two groups (p2=1-p1)"="p1","type I error"="alpha","hazard ratio for two groups"="hr",
                   "the power for which sample size is calculated"="power"), selected="n")
  })
  output$cc_caption <- reactive({paste("Figure: case-cohort design as a function of",input$cc_var)})
  output$cc <- renderPlot({
     cc_n <- input$cc_n
     cc_q <- input$cc_q
     cc_pD <- input$cc_pD
     cc_p1 <- input$cc_p1
     cc_alpha <- input$cc_alpha
     cc_hr <- input$cc_hr
     cc_power <- input$cc_power
     if (cc_power) z <- gap::ccsize(cc_n,cc_q,cc_pD,cc_p1,cc_alpha,log(cc_hr),cc_power) else z <- gap::ccsize(cc_n,cc_q,cc_pD,cc_p1,cc_alpha,log(cc_hr))
     plot(z)
  })
  output$cc_report <- downloadHandler(
    filename = function() {
      paste("cc", sep = ".", switch(input$cc_reportFormat, PDF = 'pdf', HTML = 'html', Word = 'docx'))
    },
    content = function(file) {
      src <- normalizePath('report.Rmd')
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      file.copy(src, 'report.Rmd', overwrite = TRUE)
      out <- render('report.Rmd', switch(input$cc_reportFormat, PDF = pdf_document(), HTML = html_document(), Word = word_document()))
      file.rename(out, file)
    }
  )
}
