server <- function(input, output) {
# fb design
  output$fb_var <- renderUI({
     radioButtons("fb_var", "Variable (x axis of the plot):",
                  c("Genotype relative risk"="fb_gamma","Sample size"="fb_n","fb_p","type I error"="fb_alpha","type II error"="fb_beta"))
  })
  output$fb_caption <- reactive({paste("Figure: family-based design as a function of",input$fb_var)})
  fb_data <- reactive({
     fb_gamma <- input$fb_gamma
     fb_p <- input$fb_p
     fb_alpha <- input$fb_alpha
     fb_beta <- input$fb_beta
     fb_gamma_fun <- function()
     {
       x <- seq(1,2,by=0.15)
       z <- gap::fbsize(x,fb_p)
       y <- z$n3
       xlab <- "Genotype relative risk"
       ylab <- "ASP+TDT"
       point.label <- paste(paste(xlab,sep=":",x),paste(ylab,sep=":",y),sep="\n")
       data.frame(x,y,gamma=x,p=fb_p,alpha=fb_alpha,beta=fb_beta,point.label,xlab,ylab)
     }
     fb_p_fun <- function()
     {
       x <- seq(0.05,0.1,by=0.05)
       z <- gap::fbsize(fb_gamma,x)
       y <- z$n3
       xlab <- "Frequency of disease allele"
       ylab <- "Sample size"
       point.label <- paste(paste("p:", x),paste("ASP+TDT:",y),sep="\n")
       data.frame(x,y,gamma=fb_gamma,p=x,alpha=fb_alpha,beta=fb_beta,point.label,xlab,ylab)
     }
     eval(switch(input$fb_var, a=fb_gamma_fun(), b=fb_n_fun(), c=fb_p_fun(), d=fb_alpha_fun(), e=fb_beta_fun(), fb_gamma_fun()))
  })
  output$fb_preview <- renderTable({head(fb_data())%>%select(-point.label,-xlab,-ylab)})
  output$fb <- renderPlotly({with(fb_data(), {
                                  plot_ly(x=x,y=y,type="scatter",mode="markers") %>%
                                  add_lines(x=x,y=y) %>%
                                  add_markers(text=point.label) %>%
                                  layout(xaxis=list(title=xlab[1]),yaxis=list(title=ylab[1]))
                             })
  })
  output$fb_download <- downloadHandler(
    filename = function() {paste("fb", sep=".", switch(input$fb_downloadFormat, bz2="bz2", gz="gz", tsv="tsv", xz="xz"))},
    content = function(file) {vroom_write(fb_data(), file)}
  )
# pb design
  output$pb_var <- renderUI({
     radioButtons("pb_var", "Variable (x axis of the plot):",
                   c("Prevalence of disease"="pb_kp","Sample size"="pb_n",
                     "Genotype relative risk"="pb_gamma","pb_p","type I error"="pb_alpha","type II error"="pb_beta"))
  })
  output$pb_caption <- reactive({paste("Figure: population-based design as a function of", input$pb_var)})
  pb_data <- reactive({
     pb_kp <- input$pb_kp
     pb_gamma <- input$pb_gamma
     pb_p <- input$pb_p
     pb_alpha <- input$pb_alpha
     pb_beta <- input$pb_beta
     if (input$pb_var=="pb_kp")
     {
        pb_kp <- seq(0.05,0.4,by=0.05)
        all_params <- data.frame(pb_kp,pb_gamma,pb_p,pb_alpha,pb_beta)
        x <- pb_kp
        y <- ceiling(gap::pbsize(pb_kp,pb_gamma,pb_p,pb_alpha,pb_beta))
        xlab <- "Prevalance of disease"
        ylab <- "Sample size"
        point.label <- paste(paste(xlab,sep=":",pb_kp),paste(ylab,sep=":",y),sep="\n")
     }
     data.frame(x,y,all_params,point.label,xlab,ylab)
  })
  output$pb_preview <- renderTable(head(pb_data()%>%select(-point.label,-xlab,-ylab)))
  output$pb <- renderPlotly({with(pb_data(), {
                                               plot_ly(x=x, y=y, type="scatter",mode="markers") %>%
                                               add_lines(x=x, y=y) %>%
                                               add_markers(text=point.label) %>%
                                               layout(xaxis=list(title=xlab[1]),yaxis=list(title=ylab[1]))
                             })
               })
  output$pb_download <- downloadHandler(
    filename = function() {paste("pb", sep=".", switch(input$pb_downloadFormat, bz2="bz2", gz="gz", tsv="tsv", xz="xz"))},
    content = function(file) {vroom_write(pb_data(), file)}
  )
# cc design
  output$cc_var <- renderUI({
     radioButtons("cc_var", "Variable (x axis of the plot):",
                   c("Cohort size"="cc_n","Fraction for subcohort"="cc_q","Proportion of failure in full cohort"="cc_pD",
                     "proportions of the two groups (p2=1-p1)"="cc_p1","type I error"="alpha","hazard ratio for two groups"="cc_hr",
                     "the power for which sample size is calculated"="cc_power"))
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
     cc_n_fun <- function()
     {
       x <- seq(1,100000,by=100)
       if (cc_power) z <- gap::ccsize(x,cc_q,cc_pD,cc_p1,cc_alpha,log(cc_hr),cc_power) else z <- gap::ccsize(x,cc_q,cc_pD,cc_p1,cc_alpha,log(cc_hr),cc_power)
       xlab <- "Cohort size"
       ylab <- "Subcohort sample"
       point.label <- paste(paste(xlab, sep=":", x),paste(ylab,sep=":",z),sep="\n")
       data.frame(n=x,q=cc_q,pD=cc_pD,p1=cc_p1,alpha=cc_alpha,hf=cc_hr,power=cc_power,z,point.label,xlab,ylab)
     }
     eval(switch(input$cc_var, a=cc_n_fun(), b=cc_q_fun(), c=cc_pD_fun(), d=cc_p1_fun(), e=cc_alpha_fun, f=cc_beta_fun(), cc_n_fun()))
  })
  output$cc_preview <- renderTable({head(cc_data())%>%select(-point.label,-xlab,-ylab)})
  output$cc <- renderPlotly({with(cc_data(), {
                                               plot_ly(x=n, y=z, type="scatter",mode="markers") %>%
                                               add_lines(x=n, y=z) %>%
                                               add_markers(text=point.label) %>%
                                               layout(xaxis=list(title=xlab[1]),yaxis=list(title=ylab[1]))
                             })
  })
  output$cc_download <- downloadHandler(
    filename = function() {paste("cc", sep=".", switch(input$cc_downloadFormat, bz2="bz2", gz="gz", tsv="tsv", xz="xz"))},
    content = function(file) {vroom_write(cc_data(), file)}
  )
}
