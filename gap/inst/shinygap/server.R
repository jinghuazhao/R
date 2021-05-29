server <- function(input, output) {
# fbDesign
  output$fb <- renderTable({data.frame(a=1,b=2)})
# pbDesign
  output$pb <- renderTable({data.frame(a=1,b=2)})
# ccDesign
  status <- reactive({paste(input$status)})
  time <- reactive({paste(input$time)})
  covariates <- reactive({paste(input$covariates, collapse=" + ")})
  strata <- reactive({paste(input$strata)})
# Report
  # fb curve
  output$fb_caption <- renderText({"fb"})
  # pb curve
  output$pb_caption <- renderText({"pb"})
  # cc curve
  output$cox_caption <- renderText({"cc"})
  output$report <- downloadHandler(
    filename = function() {
      paste(ifelse(input$example,"lung",tools::file_path_sans_ext(input$file)), sep = ".", 
            switch(input$reportFormat, PDF = 'pdf', HTML = 'html', Word = 'docx')
      )
    },
    content = function(file) {
      src <- normalizePath('report.Rmd')
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      file.copy(src, 'report.Rmd', overwrite = TRUE)
      out <- render('report.Rmd', switch(input$reportFormat, PDF = pdf_document(), HTML = html_document(), Word = word_document()))
      file.rename(out, file)
    }
  )
}
