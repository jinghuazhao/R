METAL_forestplot <- function(tbl)
{
  tabletext <- cbind(c("Study",study,"Summary"),
                       c("Effect",format(BETA,digits=3),format(tbl[i,"Effect"],digits=3)),
                       c("SE",format(SE,digits=3),format(tbl[i,"StdErr"],digits=3)),
                       c("N",N,tbl[i,"N"]))
  print(tabletext)
  forestplot(tabletext,
             c(NA,BETA,tbl[i,"Effect"]),
             c(NA,BETA-1.96*SE,tbl[i,"Effect"]-1.96*tbl[i,"StdErr"]),
             c(NA,BETA+1.96*SE,tbl[i,"Effect"]+1.96*tbl[i,"StdErr"]),
             zero=0,
             is.summary=c(TRUE,rep(FALSE,length(BETA)),TRUE),
             boxsize=0.75,
             col=meta.colors(box="royalblue",line="darkblue", summary="royalblue"))
  title(title)
  metaplot(BETA,SE,N,
           labels=sprintf("%s (%.3f %.3f %.0f)",study,BETA,SE,N),
           xlab="Effect distribution",ylab="",xlim=c(-1.5,1.5),
           summn=tbl[i,"Effect"],sumse=tbl[i,"StdErr"],sumnn=tbl[i,"N"],
           colors=meta.colors(box="red",lines="blue", zero="green", summary="red", text="black"))
  title(title)
}
