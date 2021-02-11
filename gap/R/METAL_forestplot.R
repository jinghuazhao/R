METAL_forestplot <- function(tbl,all,rsid,package="meta",...)
{
  prot <- MarkerName <- NA
  requireNamespace("dplyr")
  dplyr_rsid <- function(df,rsid)
  {
  # d <- dplyr::nest_join(df,rsid)
    d <- dplyr::left_join(df,rsid)
  # dy <- d["y"]
    m <- within(d, {
  #   rsid <- ifelse(length(lapply(dy,"[[",1)) == 1, unlist(d[["y"]]), unlist(lapply(lapply(dy,"[[",1),"[",1)))
  #   rsid <- unlist(d[["y"]])
      isna <- is.na(rsid)
      rsid[isna] <- MarkerName[isna]
    })
  }
  t <- dplyr_rsid(tbl,rsid)
  a <- dplyr_rsid(all,rsid)
  for(i in 1:nrow(tbl))
  {
     p <- tbl[i,"prot"]
     m <- tbl[i,"MarkerName"]
     A1 <- toupper(tbl[i,"Allele1"])
     A2 <- toupper(tbl[i,"Allele2"])
     print(paste0(i,"-",p,":",m))
     with(subset(all,prot==p & MarkerName==m), {
       print(subset(all,prot==p & MarkerName==m))
       e <- toupper(EFFECT_ALLELE)
       r <- toupper(REFERENCE_ALLELE)
       a1 <- e
       a2 <- r
       c <- rep(1,length(e))
       j <- sapply(a1,'!=',A1)
       a1[j] <- r[j]
       a2[j] <- e[j]
       c[j] <- -1
       print(cbind(A1,A2,EFFECT_ALLELE,REFERENCE_ALLELE,a1,a2,format(BETA,digits=3),format(BETA*c,digits=3)))
       BETA <- BETA * c
       title <- sprintf("%s [%s (%s) (%s/%s) N=%.0f]",p,m,t[i,"rsid"],A1,A2,tbl[i,"N"])
       if (package=="meta")
       {
         requireNamespace("meta")
         mg <- meta::metagen(BETA,SE,sprintf("%s (%.0f)",study,N),title=title)
         meta::forest(mg,colgap.forest.left = "1cm")
         requireNamespace("grid")
         grid::grid.text(title,0.5,0.9)
         with(mg,cat("prot =", p, "MarkerName =", m, "Q =", Q, "df =", df.Q, "p =", pval.Q, "I2 =", I2, "lower.I2 =", lower.I2, "upper.I2 =", upper.I2, "\n"))
       }
       else if(package=="forestplot")
       {
         tabletext <- cbind(c("Study",study,"Summary"),
                              c("Effect",format(BETA,digits=3),format(tbl[i,"Effect"],digits=3)),
                              c("SE",format(SE,digits=3),format(tbl[i,"StdErr"],digits=3)),
                              c("N",N,tbl[i,"N"]))
         print(tabletext)
         requireNamespace("forestplot")
         forestplot::forestplot(tabletext,
                    c(NA,BETA,tbl[i,"Effect"]),
                    c(NA,BETA-1.96*SE,tbl[i,"Effect"]-1.96*tbl[i,"StdErr"]),
                    c(NA,BETA+1.96*SE,tbl[i,"Effect"]+1.96*tbl[i,"StdErr"]),
                    zero=0,
                    is.summary=c(TRUE,rep(FALSE,length(BETA)),TRUE),
                    boxsize=0.75,
                    col=rmeta::meta.colors(box="royalblue",line="darkblue", summary="royalblue"))
         title(title)
       }
       else if(package=="rmeta")
       {
         requireNamespace("rmeta")
         rmeta::metaplot(BETA,SE,N,
                  labels=sprintf("%s (%.3f %.3f %.0f)",study,BETA,SE,N),
                  xlab="Effect distribution",ylab="",xlim=c(-1.5,1.5),
                  summn=tbl[i,"Effect"],sumse=tbl[i,"StdErr"],sumnn=tbl[i,"N"],
                  colors=rmeta::meta.colors(box="red",lines="blue", zero="green", summary="red", text="black"))
         title(title)
       }
     })
  }
}
