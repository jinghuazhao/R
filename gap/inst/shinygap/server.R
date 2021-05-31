server <- function(input, output) {
# fb curve
  output$fb_var <- renderUI({
     selectInput("fb_var", "Variable (x axis of the plot):",
                 c("Sample size"="fb_n","Genotype relative risk"="gamma","fb_p","type I error"="fb_alpha","type II error"="fb_beta"),
                 selected="fb_n")
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
                 c("Sample size"="pb_n","Prevalence of disease"="pb_kp",
                   "Genotype relative risk"="pb_gamma","pb_p","type I error"="pb_alpha","type II error"="pb_beta"),
                 selected="pb_gamma")
  })
  output$pb_caption <- reactive({paste("Figure: population-based design as a function of", input$pb_var)})
  output$pb <- renderPlot({
     pb_kp <- input$pb_kp
     pb_gamma <- input$pb_gamma
     pb_p <- input$pb_p
     pb_alpha <- input$pb_alpha
     pb_beta <- input$pb_beta
     data <- data.frame(pb_kp,pb_gamma,pb_p,pb_alpha,pb_beta)
#    if (input$pb_var=="pb_gamma")
     {
       pb_gamma <- seq(1,15,by=1.5)
       x <- pb_gamma
       y <- vector()
       for(z in x) y <- c(y,ceiling(gap::pbsize(pb_kp,z,pb_p)))
       data <- data.frame(x,y)
     }
     plot_ly(data,x=~x) %>%
     add_lines(y=~y,color=I("red"))
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
                 c("Cohort size"="cc_n","Fraction for subcohort"="cc_q","Proportion of failure in full cohort"="cc_pD",
                   "proportions of the two groups (p2=1-p1)"="cc_p1","type I error"="alpha","hazard ratio for two groups"="cc_hr",
                   "the power for which sample size is calculated"="cc_power"),
                 selected="cc_n")
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
