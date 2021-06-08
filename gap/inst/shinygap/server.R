server <- function(input, output) {
  options(warn=-1)
  storewarn <- getOption("warn")
# fb design
  output$fb_caption=reactive({paste("Sample size as a function of",gsub("fb_","",input$fb_var))})
  fb_data <- reactive({
     fb_gamma <- req(input$fb_gamma)
     fb_p <- req(input$fb_p)
     fb_alpha <- req(input$fb_alpha)
     fb_beta <- req(input$fb_beta)
     if (req(input$fb_var)=="fb_gamma")
     {
       x <- fb_gamma <- seq(1,fb_gamma,by=0.15)
       xlab <- "Genotype relative risk"
     }
     else if(req(input$fb_var)=="fb_p")
     {
       x <- fb_p <- seq(0.001,fb_p,by=0.01)
       xlab <- "Frequency of disease allele"
     }
     else if(req(input$fb_var)=="fb_alpha")
     {
       x <- fb_alpha
       xlab <- "type I error"
     }
     else if(req(input$fb_var)=="fb_beta")
     {
       x <- fb_beta <- seq(0.01,fb_beta,by=0.05)
       xlab <- "type II error"
     }
     y <- with(gap::fbsize(fb_gamma,fb_p,rep(fb_alpha,3),fb_beta),n3)
     ylab <- "ASP+TDT"
     point.label <- paste(paste(xlab,sep=":",x),paste(ylab,sep=":",y),sep="\n")
     data.frame(x,y,gamma=fb_gamma,p=fb_p,alpha=fb_alpha,beta=fb_beta,point.label,xlab,ylab)
  })
  output$fb_preview <- renderTable({head(fb_data())%>%select(-point.label,-xlab,-ylab)})
  output$fb <- renderPlotly({with(fb_data(), {
                                               plot_ly(x=x,y=y,type="scatter",mode="markers") %>%
                                               add_lines(x=x,y=y) %>%
                                               add_markers(text=point.label) %>%
                                               layout(xaxis=list(title=xlab[1]), yaxis=list(title=ylab[1]))
                                             }
                                 )
                            })
  output$fb_download <- downloadHandler(
    filename = function() {paste("fb", sep=".", switch(input$fb_downloadFormat, bz2="bz2", gz="gz", tsv="tsv", xz="xz"))},
    content = function(file) {vroom_write(fb_data(), file)}
  )
# pb design
  output$pb_caption <- reactive({paste("Sample size as a function of", gsub("pb_","",input$pb_var))})
  pb_data <- reactive({
     pb_kp <- req(input$pb_kp)
     pb_gamma <- req(input$pb_gamma)
     pb_p <- req(input$pb_p)
     pb_alpha <- req(input$pb_alpha)
     pb_beta <- req(input$pb_beta)
     if (req(input$pb_var)=="pb_kp")
     {
        xlab <- "Prevalance of disease"
        x <- pb_kp <- seq(1e-5,pb_kp,by=0.05)
     }
     else if (req(input$pb_var)=="pb_gamma")
     {
        xlab <- "Genotype relative risk"
        x <- pb_gamma <- seq(1,pb_gamma,by=0.15)
     }
     else if (req(input$pb_var)=="pb_p")
     {
        xlab <- "frequency of disease allele"
        x <- pb_p <- seq(1e-3,pb_p,by=0.05)
     }
     else if (req(input$pb_var)=="pb_alpha")
     {
        xlab <- "type I error"
        x <- pb_alpha <- seq(5e-8,pb_alpha,by=1e-4)
     }
     else if (req(input$pb_var)=="pb_beta")
     {
        xlab <- "type II error"
        x <- pb_beta <- seq(0.01,pb_beta,by=0.05)
     }
     y <- ceiling(gap::pbsize(pb_kp,pb_gamma,pb_p,pb_alpha,pb_beta))
     ylab <- "Sample size"
     point.label <- paste(paste(xlab,sep=":",x),paste(ylab,sep=":",y),sep="\n")
     data.frame(x, y, kp=pb_kp, gamma=pb_gamma, p=pb_p, alpha=pb_alpha, beta=pb_beta, point.label,xlab,ylab)
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
  selection <- reactive({input$cc_power})
  output$cc_caption <- reactive({paste(ifelse(selection(),"Power","Sample size"),"as a function of",gsub("cc_","",input$cc_var))})
  cc_data <- reactive({
     cc_n <- req(input$cc_n)
     cc_q <- req(input$cc_q)
     cc_pD <- req(input$cc_pD)
     cc_p1 <- req(input$cc_p1)
     cc_theta <- req(input$cc_theta)
     cc_alpha <- req(input$cc_alpha)
     cc_beta <- req(input$cc_beta)
     if (req(input$cc_var)=="cc_n")
     {
       x <- cc_n <- seq(100,cc_n,by=1000)
       xlab <- "Cohort size"
     }
     else if(req(input$cc_var)=="cc_q")
     {
       x <- cc_q <- seq(0.01,cc_q,by=0.01)
       xlab <- "Sampling fraction"
     }
     else if(req(input$cc_var)=="cc_pD")
     {
       x <- cc_pD <- seq(1e-5,cc_pD,by=0.02)
       xlab <- "Proportion of failure"
     }
     else if(req(input$cc_var)=="cc_p1")
     {
       x <- cc_p1 <- seq(1e-5,cc_p1,by=0.05)
       xlab <- "Proportion of group 1"
     }
     else if(req(input$cc_var)=="cc_theta")
     {
        x <- cc_theta <- seq(0.01,cc_theta,by=0.1)
        xlab <- "log-harzard ratio for two groups"
     }
     else if(req(input$cc_var)=="cc_alpha")
     {
        x <- cc_alpha <- seq(5e-8,cc_alpha,by=0.001)
        xlab <- "type I error"
     }
     else if(req(input$cc_var)=="cc_beta")
     {
        x <- cc_beta <- seq(0.01,cc_beta,by=0.05)
        xlab <- "type II error"
     }
     y <- ifelse(packageVersion("gap")=="1.2.3.1"), gap::ccsize(cc_n,cc_q,cc_pD,cc_p1,cc_alpha,cc_theta,cc_beta,selection())
                                                    gap::ccsize(cc_n,cc_q,cc_pD,cc_p1,cc_theta,cc_alpha,cc_beta,selection())
          )
     ylab <- ifelse(selection(), "Power", "Sample size")
     point.label <- paste(paste(xlab, sep=":", x),paste(ylab,sep=":",y),sep="\n")
     data.frame(x,y,n=cc_n,q=cc_q,pD=cc_pD,p1=cc_p1,theta=cc_theta,alpha=cc_alpha,beta=cc_beta,power=selection(),point.label,xlab,ylab)
  })
  output$cc_preview <- renderTable({filter(cc_data(),y>0)%>%select(-point.label,-xlab,-ylab)})
  output$cc <- renderPlotly({with(cc_data(), {
                                               plot_ly(x=x, y=y, type="scatter",mode="markers") %>%
                                               add_lines(x=x, y=y) %>%
                                               add_markers(text=point.label) %>%
                                               layout(xaxis=list(title=xlab[1]),yaxis=list(title=ylab[1]))
                                             }
                                 )
                            })
  output$cc_download <- downloadHandler(
    filename = function() {paste("cc", sep=".", switch(input$cc_downloadFormat, bz2="bz2", gz="gz", tsv="tsv", xz="xz"))},
    content = function(file) {vroom_write(cc_data(), file)}
  )
  options(warn=storewarn)
}
