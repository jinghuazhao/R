server <- function(input, output) {
# fb curve
  output$fb_caption <- reactive({print("fb design")})
  output$fb <- renderPlot({plot(1:10)})
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
  output$pb_caption <- reactive({print("pb design")})
  output$pb <- renderPlot({
     k <- input$pb_kp
     g <- input$pb_gamma
     p <- input$pb_p
     alpha <- input$pb_alpha
     beta <- input$pb_beta
     z <- ceiling(pbsize(k,g,p))
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
  output$cc_caption <- reactive({print("cc design")})
  output$cc <- renderPlot({plot(1:10)})
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
