server <- function(input, output) {
# fb curve
  output$fb_var <- renderUI({
     selectInput("fb_var", "Variable (x axis of the plot):",
                 c("Sample size"="fb_n","Genotype relative risk"="gamma","fb_p","type I error"="fb_alpha","type II error"="fb_beta"),
                 selected="fb_gamma")
  })
  output$fb_caption <- reactive({paste("Figure: family-based design as a function of",input$fb_var)})
  fb_data <- reactive({
     fb_gamma <- input$fb_gamma
     fb_p <- input$fb_p
     fb_alpha <- input$fb_alpha
     fb_beta <- input$fb_beta
     x <- seq(1,3,by=0.25)
     z <- gap::fbsize(fb_gamma,fb_p)
     n <- z$n3
     data.frame(x,fb_p,fb_alpha,fb_beta,n,label=paste(paste("Genotype relative risk:", x),paste("ASP+TDT:",n),sep="\n"))
  })
  output$fb_preview <- renderTable({head(fb_data())%>%select(-label)})
  output$fb <- renderPlotly({
     plot_ly(x=fb_data()$x,y=fb_data()$n,type="scatter",mode="markers") %>%
     add_lines(y=fb_data()$z$n3) %>%
     add_markers(text=fb_data()$label) %>%
     layout(scene=list(xaxis=list(title="Genotype relative risk")),yaxis=list(title="Sample size"))
  })
  output$fb_download <- downloadHandler(
    filename = function() {paste("fb", sep=".", switch(input$fb_downloadFormat, bz2="bz2", gz="gz", tsv="tsv", xz="xz"))},
    content = function(file) {vroom::vroom_write(fb_data(), file)}
  )
# pb curve
  output$pb_var <- renderUI({
     selectInput("pb_var", "Variable (x axis of the plot):",
                 c("Sample size"="pb_n","Prevalence of disease"="pb_kp",
                   "Genotype relative risk"="pb_gamma","pb_p","type I error"="pb_alpha","type II error"="pb_beta"),
                 selected="pb_kp")
  })
  output$pb_caption <- reactive({paste("Figure: population-based design as a function of", input$pb_var)})
  pb_data <- reactive({
     pb_kp <- input$pb_kp
     pb_gamma <- input$pb_gamma
     pb_p <- input$pb_p
     pb_alpha <- input$pb_alpha
     pb_beta <- input$pb_beta
#    if (input$pb_var=="pb_kp")
     {
        pb_kp <- seq(0.05,0.4,by=0.05)
        all_params <- data.frame(pb_kp,pb_gamma,pb_p,pb_alpha,pb_beta)
        x <- pb_kp
        y <- ceiling(gap::pbsize(pb_kp,pb_gamma,pb_p,pb_alpha,pb_beta))
     }
     data.frame(x,y,all_params,label=paste(paste("Prevalence of disease:",pb_kp),paste("Sample size:",y),sep="\n"))
  })
  output$pb_preview <- renderTable(head(pb_data()%>%select(-label)))
  output$pb <- renderPlotly({plot_ly(x=pb_data()$x, y=pb_data()$y, type="scatter",mode="markers") %>%
                             add_lines(y=pb_data()$y) %>%
                             add_markers(text=pb_data()$label) %>%
                             layout(scene=list(xaxis=list(title="Prevalence of disease")),yaxis=list(title="Sample size"))
               })
  output$pb_download <- downloadHandler(
    filename = function() {paste("pb", sep=".", switch(input$pb_downloadFormat, bz2="bz2", gz="gz", tsv="tsv", xz="xz"))},
    content = function(file) {vroom::vroom_write(pb_data(), file)}
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
  cc_data <- reactive({
     cc_n <- input$cc_n
     cc_q <- input$cc_q
     cc_pD <- input$cc_pD
     cc_p1 <- input$cc_p1
     cc_alpha <- input$cc_alpha
     cc_hr <- input$cc_hr
     cc_power <- input$cc_power
#    if(input$cc_var=="cc_n")
     {
        cc_n <- seq(1,100000,by=100)
     }
     if (cc_power) z <- gap::ccsize(cc_n,cc_q,cc_pD,cc_p1,cc_alpha,log(cc_hr),cc_power) else z <- gap::ccsize(cc_n,cc_q,cc_pD,cc_p1,cc_alpha,log(cc_hr))
     data.frame(cc_n,cc_q,cc_pD,cc_p1,cc_alpha,cc_hr,cc_power,z,label=paste(paste("Cohort size:",cc_n),paste("Subcohort sample:",z),sep="\n"))
  })
  output$cc_preview <- renderTable({head(cc_data())%>%select(-label)})
  output$cc <- renderPlotly({
     plot_ly(x=cc_data()$cc_n, y=cc_data()$z, type="scatter",mode="markers") %>%
     add_lines(y=cc_data()$z) %>%
     add_markers(text=cc_data()$label) %>%
     layout(scene=list(xaxis=list(title="Cohort size")),yaxis=list(title="Sample size"))
  })
  output$cc_download <- downloadHandler(
    filename = function() {paste("cc", sep=".", switch(input$cc_downloadFormat, bz2="bz2", gz="gz", tsv="tsv", xz="xz"))},
    content = function(file) {vroom::vroom_write(cc_data(), file)}
  )
}
