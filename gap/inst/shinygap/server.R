server <- function(input, output) {
  storewarn <- getOption("warn")
  options(warn=-1)
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
       x <- fb_alpha <- seq(1e-8,0.05,0.005)
       xlab <- "type I error"
     }
     else if(req(input$fb_var)=="fb_beta")
     {
       x <- fb_beta <- seq(0.01,fb_beta,by=0.05)
       xlab <- "type II error"
     }
     n <- fbsize(fb_gamma,fb_p,fb_alpha,fb_beta)
     n1 <- with(n,n1)
     n2 <- with(n,n2)
     n3 <- with(n,n3)
     ylab <- "ASP+TDT"
     point.label <- paste(paste(xlab,sep=":",x),paste(ylab,sep=":",n3),sep="\n")
     data.frame(x,n1,n2,n3,gamma=fb_gamma,p=fb_p,alpha=fb_alpha,beta=fb_beta,point.label,xlab,ylab)
  })
  output$fb_preview <- renderTable({head(fb_data())%>%select(-point.label,-xlab,-ylab)})
  output$fb <- renderPlotly({with(fb_data(), {
                                               plot_ly(type="scatter",x=~x,y=~n1,name="ASP",mode="markers") %>%
                                               add_trace(y=~n2,name="TDT",mode="marker") %>%
                                               add_trace(y=~n3,name="ASP+TDT",mode="marker",text=point.label) %>%
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
        x <- pb_alpha <- seq(1e-8,pb_alpha,by=1e-4)
     }
     else if (req(input$pb_var)=="pb_beta")
     {
        xlab <- "type II error"
        x <- pb_beta <- seq(0.01,pb_beta,by=0.05)
     }
     y <- ceiling(pbsize(pb_kp,pb_gamma,pb_p,pb_alpha,pb_beta))
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
  cc_selection <- reactive({input$cc_power})
  output$cc_caption <- reactive({paste(ifelse(cc_selection(),"Power","Sample size"),"as a function of",gsub("cc_","",input$cc_var))})
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
       x <- cc_n <- seq(100,cc_n,by=100)
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
     y <- ccsize(cc_n,cc_q,cc_pD,cc_p1,cc_theta,cc_alpha,cc_beta,cc_selection())
     ylab <- ifelse(cc_selection(), "Power", "Sample size")
     point.label <- paste(paste(xlab, sep=":", x),paste(ylab,sep=":",y),sep="\n")
     data.frame(x,y,n=cc_n,q=cc_q,pD=cc_pD,p1=cc_p1,theta=cc_theta,alpha=cc_alpha,beta=cc_beta,power=cc_selection(),point.label,xlab,ylab)
  })
  output$cc_preview <- renderTable({head(filter(cc_data(),y>0))%>%select(-point.label,-xlab,-ylab)})
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
  # tscc design
  tscc_selection <- reactive({input$tscc_model})
  output$tscc_caption <- reactive({paste("Power as a function of",gsub("_",".",gsub("tscc_","",input$tscc_var)))})
  tscc_data <- reactive({
     tscc_GRR <- req(input$tscc_GRR)
     tscc_p1 <- req(input$tscc_p1)
     tscc_n1 <- req(input$tscc_n1)
     tscc_n2 <- req(input$tscc_n2)
     tscc_M <- req(input$tscc_M)
     tscc_alpha_genome <- req(input$tscc_alpha_genome)
     tscc_pi_samples <- req(input$tscc_pi_samples)
     tscc_pi_markers <- req(input$tscc_pi_markers)
     tscc_K <- req(input$tscc_K)
     if (req(input$tscc_var)=="tscc_GRR")
     {
       x <- tscc_GRR <- seq(1,tscc_GRR,by=0.1)
       xlab <- "Genotype relative risk"
     }
     else if(req(input$tscc_var)=="tscc_p1")
     {
       x <- tscc_p1 <- seq(0.01,tscc_p1,by=0.01)
       xlab <- "Estimated risk allele frequency in cases"
     }
     else if(req(input$tscc_var)=="tscc_n1")
     {
       x <- tscc_n1 <- seq(10,tscc_n1,by=100)
       xlab <- "Total number of cases"
     }
     else if(req(input$tscc_var)=="tscc_n2")
     {
       x <- tscc_n2 <- seq(10,tscc_n2,by=100)
       xlab <- "Total number of controls"
     }
     else if(req(input$tscc_var)=="tscc_alpha_genome")
     {
       x <- tscc_alpha_genome <- seq(0.01,tscc_alpha_genome,by=0.01)
       xlab <- "False positive rate at genome level"
     }
     else if(req(input$tscc_var)=="tscc_pi_samples")
     {
       x <- tscc_pi_samples <- seq(1e-5,tscc_pi_samples,by=0.02)
       xlab <- "Sample percentage genotyped at stage 1}"
     }
     else if(req(input$tscc_var)=="tscc_pi_markers")
     {
       x <- tscc_pi_markers <- seq(1e-5,tscc_pi_markers,by=0.02)
       xlab <- "Markers percentage to be selected (also used as the false positive rate at stage 1"
     }
     else if(req(input$tscc_var)=="tscc_K")
     {
       x <- tscc_K <- seq(1e-5,tscc_K,by=0.02)
       xlab <- "The population prevalence"
     }
     z <- tscc(tscc_selection(),tscc_GRR,tscc_p1,tscc_n1,tscc_n2,tscc_M,tscc_alpha_genome,tscc_pi_samples,tscc_pi_markers,tscc_K)
     power1 <- with(z,power[1])
     power2 <- with(z,power[2])
     power3 <- with(z,power[3])
     power4 <- with(z,power[4])
     ylab <- "Power"
     point.label <- paste(paste(xlab,sep=":",x),paste(ylab,sep=":",power4),sep="\n")
     data.frame(x, power1=power1, power2=power2, power3=power3, power4=power4,
                GRR=tscc_GRR, p1=tscc_p1, n1=tscc_n1, n2=tscc_n2, M=tscc_M,
                alpha.genome=tscc_alpha_genome, pi.samples=tscc_pi_samples, pi.markers=tscc_pi_markers,K=tscc_K,
                point.label,xlab,ylab)
  })
  output$tscc_preview <- renderTable(head(tscc_data()%>%select(-point.label,-xlab,-ylab)))
  output$tscc <- renderPlotly({with(tscc_data(), {
                               plot_ly(type="scatter", x=~x, y=~power1, name="no stage", mode="markers") %>%
                               add_trace(y=~power2, name="stage 1", mode="markers") %>%
                               add_trace(y=~power3, name="stage 2", mode="markers") %>%
                               add_trace(y=~power4, name="joint", mode="markers", text=point.label) %>%
                               layout(xaxis=list(title=xlab[1]),yaxis=list(title=ylab[1]))
                               })
                 })
  output$tscc_download <- downloadHandler(
    filename = function() {paste("tscc", sep=".", switch(input$pb_downloadFormat, bz2="bz2", gz="gz", tsv="tsv", xz="xz"))},
    content = function(file) {vroom_write(tscc_data(), file)}
  )
  options(warn=storewarn)
}
