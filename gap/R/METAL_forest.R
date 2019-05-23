METAL_forest <- function(tbl,all,rsid)
{
  require(dplyr)
  m <- within(nest_join(tbl,rsid),{rsid <- unlist(lapply(lapply(y,"[[",1),"[",1))})
  isna <- with(m, is.na(rsid))
  t <- within(m, {rsid[isna] <- MarkerName[isna]})
  m <- within(nest_join(all,rsid),{rsid <- unlist(lapply(lapply(y,"[[",1),"[",1))})
  isna <- with(m, is.na(rsid))
  a <- within(m, {rsid[isna] <- MarkerName[isna]})
  for(i in 1:nrow(tbl))
  {
     p <- tbl[i,"prot"]
     m <- tbl[i,"MarkerName"]
     d <- gsub("[?]","",tbl[i,"Direction"])
     s <- unlist(strsplit(d,""))
     f <- as.numeric(paste0(s,1))
     A1 <- toupper(tbl[i,"Allele1"])
     A2 <- toupper(tbl[i,"Allele2"])
     print(paste0(i,"-",p,":",m))
     with(subset(all,prot==p & MarkerName==m), {
       e <- toupper(EFFECT_ALLELE)
       r <- toupper(REFERENCE_ALLELE)
       a1 <- a2 <- vector('character',length(e))
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
       require(meta)
       mg <- metagen(BETA,SE,sprintf("%s (%.0f)",study,N),title=title)
       forest(mg,colgap.forest.left = "1cm")
       require(grid)
       grid.text(title,0.5,0.9)
#      METAL_forestplot(tbl)
     })
  }
}
