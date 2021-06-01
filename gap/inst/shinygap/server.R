server <- function(input, output) {
# fb design
  output$fb_var <- renderUI({
     radioButtons("fb_var", "Variable (x axis of the plot):",
                  c("Genotype relative risk (gamma)"="fb_gamma",
                    "frequency of disease allele (p)"="fb_p",
                    "type I error (alpha)"="fb_alpha",
                    "type II error (beta)"="fb_beta"
                   )
                 )
  })
  output$fb_caption <- reactive({paste("Figure: family-based design as a function of",input$fb_var)})
  fb_data <- reactive({
     fb_gamma <- input$fb_gamma
     fb_p <- input$fb_p
     fb_alpha <- input$fb_alpha
     fb_beta <- input$fb_beta
     if (input$fb_var=="fb_gamma")
     {
       x <- fb_gamma <- seq(1,5,by=0.15)
       xlab <- "Genotype relative risk"
     }
     else if(input$fb_var=="fb_p")
     {
       x <- fb_p <- seq(0.05,0.4,by=0.05)
       xlab <- "Frequency of disease allele"
     }
     else if(input$fb_var=="fb_alpha")
     {
       x <- fb_alpha <- seq(0.0001,1e-4,by=5e-8)
       xlab <- "type I error"
     }
     else if(input$fb_var=="fb_beta")
     {
       x <- fb_beta <- seq(0.01,0.4,by=0.05)
       xlab <- "type II error"
     }
     y <- with(gap::fbsize(fb_gamma,fb_p,fb_alpha,fb_beta),n3)
     ylab <- "ASP+TDT"
     point.label <- paste(paste(xlab,sep=":",x),paste(ylab,sep=":",y),sep="\n")
     data.frame(x,y,gamma=fb_gamma,p=fb_p,alpha=fb_alpha,beta=fb_beta,point.label,xlab,ylab)
  })
  output$fb_preview <- renderTable({head(fb_data())%>%select(-point.label,-xlab,-ylab)})
  output$fb <- renderPlotly({with(fb_data(), {
                                               plot_ly(x=x,y=y,type="scatter",mode="markers") %>%
                                               add_lines(x=x,y=y) %>%
                                               add_markers(text=point.label) %>%
                                               layout(xaxis=list(title=xlab[1]),yaxis=list(title=ylab[1]))
                                             }
                                 )
                            })
  output$fb_download <- downloadHandler(
    filename = function() {paste("fb", sep=".", switch(input$fb_downloadFormat, bz2="bz2", gz="gz", tsv="tsv", xz="xz"))},
    content = function(file) {vroom_write(fb_data(), file)}
  )
# pb design
  output$pb_var <- renderUI({
     radioButtons("pb_var", "Variable (x axis of the plot):",
                   c("Prevalence of disease (K)"="pb_kp",
                     "Genotype relative risk (gamma)"="pb_gamma",
                     "Frequency of disease allele (p)"="pb_p",
                     "type I error (alpha)"="pb_alpha",
                     "type II error (beta)"="pb_beta"
                    )
                 )
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
        x <- pb_kp <- seq(0.05,0.4,by=0.05)
        xlab <- "Prevalance of disease"
     }
     else if (input$pb_var=="pb_gamma")
     {
        x <- pb_gamma <- seq(1,2,by=0.15)
        xlab <- "Genotype relative risk"
     }
     else if (input$pb_var=="pb_p")
     {
        x <- pb_p <- seq(0.05,0.1,by=0.05)
        xlab <- "frequency of disease allele"
     }
     else if (input$pb_var=="pb_alpha")
     {
        x <- pb_alpha <- seq(0.0001,0.01,by=0.001)
        xlab <- "type I error"
     }
     else if (input$pb_var=="pb_beta")
     {
        x <- pb_beta <- seq(0.01,0.4,by=0.05)
        xlab <- "type II error"
     }
     y <- ceiling(gap::pbsize(pb_kp,pb_gamma,pb_p,pb_alpha,pb_beta))
     ylab <- "Sample size"
     point.label <- paste(paste(xlab,sep=":",x),paste(ylab,sep=":",y),sep="\n")
     data.frame(x,y, kp=pb_kp, gamma=pb_gamma, p=pb_p, alpha=pb_alpha, beta=pb_beta, point.label,xlab,ylab)
  })
  output$pb_preview <- renderTable(head(pb_data()%>%select(-point.label,-xlab,-ylab)))
  output$pb <- renderPlotly({with(pb_data(), {
                                               plot_ly(x=x, y=y, type="scatter",mode="markers") %>%
                                               add_lines(x=x, y=y) %>%
                                               add_markers(text=point.label) %>%
                                               layout(xaxis=list(title=xlab[1]),yaxis=list(title=ylab[1]))
                                             }
                                 )
                            })
  output$pb_download <- downloadHandler(
    filename = function() {paste("pb", sep=".", switch(input$pb_downloadFormat, bz2="bz2", gz="gz", tsv="tsv", xz="xz"))},
    content = function(file) {vroom_write(pb_data(), file)}
  )
# cc design
  output$cc_var <- renderUI({
     radioButtons("cc_var", "Variable (x axis of the plot):",
                   c("Cohort size (n)"="cc_n",
                     "Fraction for subcohort (q)"="cc_q",
                     "Proportion of failure (pD)"="cc_pD",
                     "Proportion of group 1 (p1)"="cc_p1",
                     "type I error (alpha)"="cc_alpha",
                     "log(HR) for two groups (theta)"="cc_theta"
                    )
                 )
  })
  output$cc_caption <- reactive({paste("Figure: case-cohort design as a function of",input$cc_var)})
  cc_data <- reactive({
     cc_n <- input$cc_n
     cc_q <- input$cc_q
     cc_pD <- input$cc_pD
     cc_p1 <- input$cc_p1
     cc_alpha <- input$cc_alpha
     cc_theta <- input$cc_theta
     cc_power <- input$cc_power
     if (input$cc_var=="cc_n")
     {
       x <- cc_n <- seq(1,100000,by=100)
       xlab <- "Cohort size"
     }
     else if(input$cc_var=="cc_q")
     {
       x <- cc_q <- seq(0.01,0.4,by=0.05)
       xlab <- "Sampling fraction"
     }
     else if(input$cc_var=="cc_pD")
     {
       x <- cc_pD <- seq(0.01,0.4,by=0.05)
       xlab <- "Proportion of failure"
     }
     else if(input$cc_var=="cc_p1")
     {
       x <- cc_p1 <- seq(0.01,0.4,by=0.05)
       xlab <- "Proportion of group 1"
     }
     else if(input$cc_var=="cc_theta")
     {
        x <- cc_theta <- seq(1,10,by=1.2)
        xlab <- "log-harzard ratio for two groups"
     }
     else if(input$cc_var=="cc_alpha")
     {
        x <- cc_alpha <- seq(0.0001,1e-4,by=5e-8)
        xlab <- "type I error"
     }
     ylab <- switch(input$cc_power,"Power","Sample size")
     z <- gap::ccsize(cc_n,cc_q,cc_pD,cc_p1,cc_alpha,cc_theta,cc_power)
     point.label <- paste(paste(xlab, sep=":", x),paste(ylab,sep=":",z),sep="\n")
     data.frame(n=cc_n,q=cc_q,pD=cc_pD,p1=cc_p1,alpha=cc_alpha,theta=cc_theta,z,point.label,xlab,ylab)
  })
  output$cc_preview <- renderTable({head(cc_data())%>%select(-point.label,-xlab,-ylab)})
  output$cc <- renderPlotly({with(cc_data(), {
                                               plot_ly(x=n, y=z, type="scatter",mode="markers") %>%
                                               add_lines(x=n, y=z) %>%
                                               add_markers(text=point.label) %>%
                                               layout(xaxis=list(title=xlab[1]),yaxis=list(title=ylab[1]))
                                             }
                                 )
                            })
  output$cc_download <- downloadHandler(
    filename = function() {paste("cc", sep=".", switch(input$cc_downloadFormat, bz2="bz2", gz="gz", tsv="tsv", xz="xz"))},
    content = function(file) {vroom_write(cc_data(), file)}
  )
}
