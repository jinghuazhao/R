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
     if (req(input$fb_var)=="fb_gamma")
     {
       x <- fb_gamma <- seq(1,5,by=0.15)
       xlab <- "Genotype relative risk"
     }
     else if(req(input$fb_var)=="fb_p")
     {
       x <- fb_p <- seq(0.05,0.4,by=0.05)
       xlab <- "Frequency of disease allele"
     }
     else if(req(input$fb_var)=="fb_alpha")
     {
       x <- fb_alpha <- seq(1e-8,1e-4,by=5e-8)
       xlab <- "type I error"
     }
     else if(req(input$fb_var)=="fb_beta")
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
     if (req(input$pb_var)=="pb_kp")
     {
        xlab <- "Prevalance of disease"
        x <- pb_kp <- seq(0.05,0.4,by=0.05)
     }
     else if (req(input$pb_var)=="pb_gamma")
     {
        xlab <- "Genotype relative risk"
        x <- pb_gamma <- seq(1,2,by=0.15)
     }
     else if (req(input$pb_var)=="pb_p")
     {
        xlab <- "frequency of disease allele"
        x <- pb_p <- seq(0.05,0.1,by=0.05)
     }
     else if (req(input$pb_var)=="pb_alpha")
     {
        xlab <- "type I error"
        x <- pb_alpha <- seq(0.0001,0.01,by=0.001)
     }
     else if (req(input$pb_var)=="pb_beta")
     {
        xlab <- "type II error"
        x <- pb_beta <- seq(0.01,0.4,by=0.05)
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
                     "log(HR) for two groups (theta)"="cc_theta",
                     "type I error (alpha)"="cc_alpha",
                     "type II error (beta)"="cc_beta"
                    )
                 )
  })
  output$cc_caption <- reactive({paste("Figure: case-cohort design as a function of",input$cc_var)})
  cc_data <- reactive({
     cc_n <- input$cc_n
     cc_q <- input$cc_q
     cc_pD <- input$cc_pD
     cc_p1 <- input$cc_p1
     cc_theta <- input$cc_theta
     cc_alpha <- input$cc_alpha
     cc_beta <- input$cc_beta
     cc_power <- input$cc_power
     if (req(input$cc_var)=="cc_n")
     {
       x <- cc_n <- seq(1,500000,by=100)
       xlab <- "Cohort size"
     }
     else if(req(input$cc_var)=="cc_q")
     {
       x <- cc_q <- seq(0.01,0.1,by=0.02)
       xlab <- "Sampling fraction"
     }
     else if(req(input$cc_var)=="cc_pD")
     {
       x <- cc_pD <- seq(0.01,0.4,by=0.02)
       xlab <- "Proportion of failure"
     }
     else if(req(input$cc_var)=="cc_p1")
     {
       x <- cc_p1 <- seq(0.01,0.4,by=0.05)
       xlab <- "Proportion of group 1"
     }
     else if(req(input$cc_var)=="cc_theta")
     {
        x <- cc_theta <- seq(0.02,2.3,by=0.1)
        xlab <- "log-harzard ratio for two groups"
     }
     else if(req(input$cc_var)=="cc_alpha")
     {
        x <- cc_alpha <- seq(0.0001,0.4,by=0.05)
        xlab <- "type I error"
     }
     else if(req(input$cc_var)=="cc_beta")
     {
        x <- cc_alpha <- seq(0,0.4,by=0.1)
        xlab <- "type II error"
     }
     power <- req(input$cc_power)
     ylab <- ifelse(power, "Power", "Sample size")
     z <- gap::ccsize(cc_n,cc_q,cc_pD,cc_p1,cc_theta,cc_alpha,cc_beta,power)
     point.label <- paste(paste(xlab, sep=":", x),paste(ylab,sep=":",z),sep="\n")
     data.frame(n=cc_n,q=cc_q,pD=cc_pD,p1=cc_p1,theta=cc_theta,alpha=cc_alpha,beta=cc_beta,z,point.label,xlab,ylab)
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
