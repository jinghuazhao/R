#' forest plot as R/meta's forest for METAL outputs
#' 
#' @param tbl Meta-anslysis summary statistics.
#' @param all statistics from all contributing studies.
#' @param rsid SNPID-rsid mapping file.
#' @param package "meta" or "metafor" package.
#' @param method an explcit flag for fixed/random effects model.
#' @param split when TRUE, individual prot-MarkerName.pdf will be generated.
#' @param ... Additional arguments to `meta::forest` or `metafor::forest`.
#'
#' @details
#' This functions takes a meta-data from METAL (tbl) and data from contributing studies (all)
#' for forest plot. It also takes a SNPID-rsid mapping (rsid) as contributing studies often
#' involve discrepancies in rsid so it is appropriate to use SNPID, i.e., chr:pos_A1_A2 (A1<=A2).
#'
#' The study-specific and total sample sizes (`N`) can be customised from METAL commands. By default, the input triplets each contain
#' a `MarkerName` variable which is the unique SNP identifier (e.g., chr:pos:a1:a2) and the `tbl` argument has variables
#' `A1` and `A2` as produced by METAL while the `all` argument has `EFFECT_ALLELE` and `REFERENCE_ALLELE` as with a `study` variable
#' indicating study name. Another variable common the `tbl` and `all` is `prot` variable as the function was developed in a protein
#' based meta-analysis. As noted above, the documentation example also has variable `N`.
#' From these all information is in place for generation of a list of forest plots through a batch run.
#'
#' CUSTOMVARIABLE N\cr
#' LABEL N as N\cr
#' WEIGHTLABEL N
#'
#' @export
#' @return
#' It will generate a forest plot specified by pdf for direction-adjusted effect sizes.
#'
#' @references
#' \insertRef{schwarzer07}{gap}
#'
#' \insertRef{willer10}{gap}
#'
#' @examples
#' \dontrun{
#'  data(OPG, package="gap.datasets")
#'  meta::settings.meta(method.tau="DL")
#'  METAL_forestplot(OPGtbl,OPGall,OPGrsid,width=8.75,height=5,digits.TE=2,digits.se=2,
#'                   col.diamond="black",col.inside="black",col.square="black")
#'  METAL_forestplot(OPGtbl,OPGall,OPGrsid,package="metafor",method="FE",xlab="Effect",
#'                   showweights=TRUE)
#' }
#'
#' @author Jing Hua Zhao
#' @keywords hplot distribution

METAL_forestplot <- function(tbl,all,rsid,package="meta",method="REML",split=FALSE,...)
{
  prot <- MarkerName <- NA
  requireNamespace("dplyr")
  requireNamespace("grid")
  requireNamespace("meta")
  requireNamespace("metafor")
  dplyr_rsid <- function(df,rsid)
  {
    d <- dplyr::left_join(df,rsid)
    m <- within(d, {
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
     TITLE <- sprintf("%s [%s (%s) (%s/%s) N=%.0f]",p,m,t[i,"rsid"],A1,A2,tbl[i,"N"])
     sall <- subset(all,prot==p & MarkerName==m)
     with(sall, {
       e <- toupper(EFFECT_ALLELE)
       r <- toupper(REFERENCE_ALLELE)
       a1 <- e
       a2 <- r
       c <- rep(1,length(e))
       j <- sapply(a1,'!=',A1)
       a1[j] <- r[j]
       a2[j] <- e[j]
       c[j] <- -1
       BETA <- BETA * c
       print(cbind(A1,A2,EFFECT_ALLELE,REFERENCE_ALLELE,a1,a2,format(BETA,digits=3),format(BETA*c,digits=3)))
       if (split) pdf(paste0(p,"-",m,".pdf"))
       if (package=="meta")
       {
         meta::settings.meta(method.tau=method)
         mg <- meta::metagen(BETA,SE,sprintf("%s (%.0f)",study,N),title=TITLE,method.tau.ci="")
         meta::forest(mg,colgap.forest.left = "1cm",leftlabs=c("Study","b","SE"),...)
         grid::grid.text(TITLE,0.5,0.9)
         with(mg,cat("prot =", p, "MarkerName =", m, "Q =", Q, "df =", df.Q, "p =", pval.Q,
                     "I2 =", I2, "[", lower.I2, ",", upper.I2, "]\n"))
       } else {
         d <- metafor::escalc(measure="MN",yi=BETA,sei=SE)
         r <- metafor::rma(yi,vi,data=d,method=method,slab=paste0(study," (",N,")"))
         h <- c(TITLE,paste(ifelse(list(...)$showweights,"Weight",""),"Effect (95%CI)",sep="\t"))
         f <- metafor::forest(r, header=h,...)
       # text(f$alim[2],r$k+2,"Weight",cex=1,font=2)
       }
       if (split) dev.off()
     })
  }
}
