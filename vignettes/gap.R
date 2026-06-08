set.seed(0)
options(rmarkdown.html_vignette.check_title = FALSE)
knitr::opts_chunk$set(
  out.extra = 'style="display:block; margin: auto"',
  fig.align = "center",
  fig.path = "./",
  collapse = TRUE,
  comment = "#>",
  dev = "png"
)
if (dev.cur() == 1) {
  grDevices::png(filename = tempfile(fileext = ".png"))
}

p1 <- "
 1  2   3  2  2  7/7  7/10
 2  0   0  1  1  -/-  -/-
 3  0   0  2  2  7/9  3/10
 4  2   3  2  2  7/9  3/7
 5  2   3  2  1  7/7  7/10
 6  2   3  1  1  7/7  7/10
 7  2   3  2  1  7/7  7/10
 8  0   0  1  1  -/-  -/-
 9  8   4  1  1  7/9  3/10
10  0   0  2  1  -/-  -/-
11  2  10  2  1  7/7  7/7
12  2  10  2  2  6/7  7/7
13  0   0  1  1  -/-  -/-
14 13  11  1  1  7/8  7/8
15  0   0  1  1  -/-  -/-
16 15  12  2  1  6/6  7/7
"

p2 <- as.data.frame(scan(file=textConnection(p1),what=list(0,0,0,0,0,"","")))
names(p2) <-c("id","fid","mid","sex","aff","GABRB1","D4S1645")
p3 <- data.frame(pid=10081,p2)

library(gap)
knitr::kable(p3,caption="An example pedigree")
# visible diagram in RStudio / toDOT=TRUE here
gap::pedtodot_verbatim(p3,run=TRUE)
library(DiagrammeR)
gap::pedtodot_verbatim(p3)
DiagrammeR::grViz(readr::read_file("10081.dot"))

data(lukas, package="gap.datasets")
library(kinship2)
ped <- with(lukas,pedigree(id,father,mother,sex))
plot(ped,cex=0.4)

options(width=150)
library(gap)
models <- data.frame(
  gamma = c(4,4,4,4, 2,2,2,2, 1.5,1.5,1.5,1.5),
  p     = c(0.01,0.10,0.50,0.80,
            0.01,0.10,0.50,0.80,
            0.01,0.10,0.50,0.80)
)
res <- t(mapply(function(g, p) {
  z <- fbsize(g, p)
  with(z, c(gamma, p, y, n1, pA, h1, n2, h2, n3, lambdao, lambdas))
}, models$gamma, models$p))
table1 <- as.data.frame(res)
names(table1) <- c("gamma","p","Y","N_asp","P_A","H1",
                   "N_tdt","H2","N_asp_tdt","L_o","L_s")
int_cols <- c("N_asp", "N_tdt", "N_asp_tdt")
dec_cols <- c("Y", "P_A", "H1", "H2", "L_o", "L_s")
table1[int_cols] <- lapply(table1[int_cols], ceiling)
table1[dec_cols] <- lapply(table1[dec_cols], round, 2)
knitr::kable(table1, caption = "Power/Sample size of family-based designs")

g <- 4.5
p <- 0.15
alz <- data.frame(fbsize(g,p))
knitr::kable(alz,caption="Power/Sample size of study on Alzheimer's disease")

library(gap)
kp <- c(0.01, 0.05, 0.10, 0.20)
models <- data.frame(
  gamma = c(4,4,4,4, 2,2,2,2, 1.5,1.5,1.5,1.5),
  p     = c(0.01,0.10,0.50,0.80,
            0.01,0.10,0.50,0.80,
            0.01,0.10,0.50,0.80)
)
res <- t(mapply(function(g, p) {
  ceiling(sapply(kp, function(k) pbsize(k, g, p)))
}, models$gamma, models$p))
table5 <- cbind(models, as.data.frame(res))
names(table5) <- c("gamma", "p", "p1", "p5", "p10", "p20")
knitr::kable(table5, caption = "Sample size of population-based design")

library(gap)
# ARIC study
n <- 15792; pD <- 0.03; p1 <- 0.25; alpha <- 0.05; beta <- 0.2
hr <- c(1.35, 1.40, 1.45); q <- c(1463, 722, 468) / n
aric <- data.frame(
  n, pD, p1, hr, q,
  power = signif(mapply(ccsize, n, q, pD, p1, log(hr),
                        MoreArgs = list(alpha = alpha, beta = beta, power = TRUE)), 3),
  ssize = mapply(ccsize, n, q, pD, p1, log(hr),
                 MoreArgs = list(alpha = alpha, beta = beta, power = FALSE))
)
aric
# EPIC study
n <- 25000; q <- 0.1; alpha <- 5e-8; beta <- 0.2
epic <- subset(
  transform(
    expand.grid(pD = c(0.3, 0.2, 0.1, 0.05),
                p1 = seq(0.1, 0.5, 0.1),
                hr = seq(1.1, 1.4, 0.1)),
    n = n,
    alpha = formatC(alpha, format = "e", digits = 2),
    ssize = mapply(ccsize, n, q, pD, p1, log(hr),
                   MoreArgs = list(alpha = alpha, beta = beta, power = FALSE))
  ),
  !is.na(ssize) & ssize > 0
)
knitr::kable(epic,caption="Sample size of case-cohort designs")

library(gap)
u_obs <- runif(1000)
r <- qqunif(u_obs,pch=21,bg="blue",bty="n")
u_exp <- r$y
hits <- u_exp >= 2.30103
points(r$x[hits],u_exp[hits],pch=21,bg="green")
legend("topleft",sprintf("GC.lambda=%.4f",gc.lambda(u_obs)))

w4 <- w4[order(w4$chr, w4$pos), ]
colors <- c(rep(c("blue","red"),15),"red")
suggestiveline <- -log10(3.60036E-05)
genomewideline <- -log10(1.8E-06)
mhtplot(
  w4,
  control = mht.control(
    colors  = colors,
    gap     = 1000,
    cex     = 0.6,
    cutoffs = c(suggestiveline,genomewideline),
    lab.cex = 0.6,
    xline   = 3,
    yline   = 3
  ),
  pch = 19,
  srt = 0
)

data <- with(mhtdata,cbind(chr,pos,p))
glist <- c("IRS1","SPRY2","FTO","GRIK3","SNED1","HTR1A","MARCH3","WISP3",
           "PPP1R3B","RP1L1","FDFT1","SLC39A14","GFRA1","MC4R")
hdata <- subset(mhtdata,gene%in%glist)[c("chr","pos","p","gene")]
hcolor <- rep("red",length(glist))
ops <- mht.control(axis.tck=0.02,cex=0.4,lab.cex=0.8,xline=4)
hops <- hmht.control(data=hdata,cex=0.8,colors=hcolor,boxed=TRUE)
mhtplot(data,ops,hops,pch=19)

circos.mhtplot(mhtdata, glist)

require(gap.datasets)
library(dplyr)
testdat <- mhtdata[c("chr","pos","p","gene","start","end")] %>%
           rename(log10p=p) %>%
           mutate(chr=paste0("chr",chr),log10p=-log10(log10p))
dat <- mutate(testdat,start=pos,end=pos) %>%
       select(chr,start,end,log10p)
labs <- subset(testdat,gene %in% glist) %>%
        group_by(gene,chr,start,end) %>%
        summarize() %>%
        mutate(cols="blue") %>%
        select(chr,start,end,gene,cols)
labs[2,"cols"] <- "red"
ticks <- 0:2*5
circos.mhtplot2(dat,labs,ticks=ticks,ymax=max(ticks))

mhtdata <- within(mhtdata,{pr=p})
miamiplot(mhtdata,chr="chr",bp="pos",p="p",pr="pr",snp="rsn")
# An alternative implementation
gwas <- select(mhtdata,chr,pos,p) %>%
        mutate(z=qnorm(p/2))
chrmaxpos <- miamiplot2(gwas,gwas,name1="Batch 2",name2="Batch 1",z1="z",z2="z")
labelManhattan(chr=c(2,16),pos=c(226814165,52373776),name=c("AnonymousGene","FTO"),gwas,gwasZLab="z",chrmaxpos=chrmaxpos)

knitr::include_graphics("IL-12B_mhtplot.trunc.png")

asplot(CDKNlocus, CDKNmap, CDKNgenes, best.pval=5.4e-8, sf=c(3,6))
title("CDKN2A/CDKN2B Region")

cnvplot(gap.datasets::cnv)

library(gap)
rs12075 <- data.frame(id=c("CCL2","CCL7","CCL8","CCL11","CCL13","CXCL6","Monocytes"),
                      b=c(0.1694,-0.0899,-0.0973,0.0749,0.189,0.0816,0.0338387),
                      se=c(0.0113,0.013,0.0116,0.0114,0.0114,0.0115,0.00713386))
ESplot(rs12075)

data(OPG,package="gap.datasets")
meta::settings.meta(method.tau="DL")
METAL_forestplot(OPGtbl,OPGall,OPGrsid,width=6.75,height=5,digits.TE=2,digits.se=2,
                 col.diamond="black",col.inside="black",col.square="black")
METAL_forestplot(OPGtbl,OPGall,OPGrsid,package="metafor",method="FE",xlab="Effect",
                 showweights=TRUE)

library(lattice)
z0 <- 1.96
z <- 5
a <- seq(-z, z, length = 10000)
b <- dnorm(a, 0, 1)
xyplot(b ~ a,
       type = "l",
       panel = function(x,y, ...)
               {
                   panel.xyplot(x,y, ...)
                   panel.abline(v = 0, lty = 2)
                   xx <- c(-z, x[x>=-z & x<=-z0], -z0)
                   yy <- c(0, y[x>=-z & x<=-z0], 0)
                   panel.polygon(xx,yy, ..., col='red')
                   xx <- c(z0, x[x>=z0 & x<=z], z)
                   yy <- c(0, y[x>=z0 & x<=z], 0)
                   panel.polygon(xx,yy, ..., col='red')
               },
       xlab="z",
       ylab=expression(1/sqrt(2*pi) * exp(-z^2/2))
)

require(gap)
zlist <- c(5,10,30,40,50,100,500,1000,2000,3000,5000)
zp <- sapply(zlist,function(z) {c(z,pvalue(z),logp(z),log10p(z))})
rownames(zp) <- c("z","P","log(P)","log10(P)")
knitr::kable(t(zp),caption="z, P, log(P) and log10(P)")

oldpar <- par()
par(mfrow=c(1,2))
N <- 2:500
R2 <- 2*(N-2)/(N-1)^2/(N+1)
R2LR <- 2/(N+1)^2
R2t <- 2/(N-1)^2
plot(N,R2,cex=0.6,xaxt="n",xlab="Sample size",ylab=expression(Var(R^2)),col="black",pch=20)
points(N,R2LR,cex=0.6,pch=15,col="red")
points(N,R2t,cex=0.6,pch=17,col="blue")
axis(1,at=c(2,(1:5)*100))
legend(400,0.03,c("Asymptotic","LR","t-statistic"),col=c("black","red","blue"),pch=c(20,15,17))
sLR <- (N-2)*(N+1)/(N-1)^2
st <- (N-2)/(N+1)
plot(N,sLR,cex=0.6,xaxt="n",xlab="Sample size",ylab="Asymptotic/approximation estimator ratio",col="red",pch=20)
points(N,st,cex=0.6,pch=15,col="blue")
abline(h=1,col="black")
axis(1,at=c(2,(1:5)*100))
legend(400,0.2,c("Asymptotic","LR","t-statistic"),col=c("black","red","blue"),pch=c(20,15,17))
par(oldpar)

get_b_se(0.6396966,23991,4.7245)

get_sdy(0.6396966,23991,0.04490488,0.009504684)

txt <- '
 rs188743906  0.6804   0.1104  0.00177 0.01660        NA        NA
   rs2289779 -0.0788   0.0134  0.00104 0.00261 -0.007543 0.0092258
 rs117804300 -0.2281   0.0390 -0.00392 0.00855  0.109372 0.0362219
   rs7033492 -0.0968   0.0147 -0.00585 0.00269  0.022793 0.0119903
  rs10793962  0.2098   0.0212  0.00378 0.00536 -0.014567 0.0138196
    rs635634 -0.2885   0.0153 -0.02040 0.00334  0.077157 0.0117123
    rs176690 -0.0973   0.0142  0.00293 0.00306 -0.000007 0.0107781
 rs147278971 -0.2336   0.0378 -0.01240 0.00792  0.079873 0.0397491
  rs11562629  0.1155   0.0181  0.00960 0.00378 -0.010040 0.0151460
'
v <- c("SNP","b.LIF.R","SE.LIF.R","b.FEV1","SE.FEV1","b.CAD","SE.CAD")
mrdat <- setNames(as.data.frame(scan(text = txt, what = list("",0,0,0,0,0,0), quiet = TRUE)),v)
knitr::kable(mrdat,caption="LIF-R and CAD/FEV1")

res <- mr(mrdat, "LIF.R", c("CAD","FEV1"), other_plots=TRUE)
r <- res$r
rownames(r) <- r[,1]
r <- r[,-1,drop=FALSE]
r <- apply(r, 2, as.numeric)
methods <- c("IVW","EGGER","WM","PWM")
p <- sapply(methods, function(m) 2*pnorm(-abs(r[,paste0("b",m)]/r[,paste0("seb",m)])))
colnames(p) <- paste0("p",methods)
knitr::kable(t(data.frame(round(r,3), format(p,3,scientific=TRUE), check.names=FALSE)),
             align="r", caption="LIFR variant rs635634 and CAD/FEV1"
)

mr_names <- names(mrdat)
LIF.R <- cbind(mrdat[grepl("SNP|LIF.R",mr_names)],trait="LIF.R"); names(LIF.R) <- c("SNP","b","se","trait")
FEV1 <- cbind(mrdat[grepl("SNP|FEV1",mr_names)],trait="FEV1"); names(FEV1) <- c("SNP","b","se","trait")
CAD <- cbind(mrdat[grepl("SNP|CAD",mr_names)],trait="CAD"); names(CAD) <- c("SNP","b","se","trait")
mrdat2 <- within(rbind(LIF.R,FEV1,CAD),{y=b})
library(ggplot2)
p <- ggplot2::ggplot(mrdat2,aes(y = SNP, x = y))+
     ggplot2::theme_bw()+
     ggplot2::geom_point()+
     ggplot2::facet_wrap(~ trait, ncol=3, scales="free_x")+
     ggplot2::geom_segment(aes(x = b-1.96*se, xend = b+1.96*se, yend = SNP))+
     ggplot2::geom_vline(lty=2, ggplot2::aes(xintercept=0), colour = 'red')+
     ggplot2::xlab("Effect size")+
     ggplot2::ylab("")
p

tnfb <- '
              "multiple sclerosis"  0.69058600 0.059270400
    "systemic lupus erythematosus"  0.76687500 0.079000500
          "sclerosing cholangitis"  0.62671500 0.075954700
   "juvenile idiopathic arthritis" -1.17577000 0.160293000
                       "psoriasis"  0.00582586 0.000800016
            "rheumatoid arthritis" -0.00378072 0.000625160
      "inflammatory bowel disease" -0.14334200 0.025272500
          "ankylosing spondylitis" -0.00316852 0.000626225
                  "hypothyroidism" -0.00432054 0.000987324
               "allergic rhinitis"  0.00393075 0.000926002
          "IgA glomerulonephritis" -0.32696600 0.105262000
                   "atopic eczema" -0.00204018 0.000678061
'

tnfb <- as.data.frame(scan(file=textConnection(tnfb),what=list("",0,0))) %>%
        setNames(c("outcome","Effect","StdErr")) %>%
        mutate(outcome=gsub("\\b(^[a-z])","\\U\\1",outcome,perl=TRUE))

mr_forestplot(tnfb, fontsize=12,
             leftcols=c("studlab","effect","seTE","ci"), leftlabs=c("Outcome","b","SE","95% CI"),
             rightcols=c("w.common","w.random"),rightlabs=c("Weight (FE)","Weight (RE)"),
             common=FALSE, random=FALSE, print.I2=FALSE, print.pval.Q=FALSE, print.tau2=FALSE,
             spacing=1.6, digits.TE=2, digits.seTE=2, xlab="Effect size", type.study="square", col.inside="black", col.square="black")

mr_forestplot(tnfb, colgap.forest.left="0.05cm", fontsize=14,
              leftcols="studlab", leftlabs="Outcome", plotwidth="3inch", sm="OR", rightlabs="ci",
              sortvar=tnfb[["Effect"]],
              common=FALSE, random=FALSE, print.I2=FALSE, print.pval.Q=FALSE, print.tau2=FALSE,
              backtransf=TRUE, spacing=1.6,type.study="square",col.inside="black",col.square="black")

mr_forestplot(tnfb,colgap.forest.left="0.05cm", fontsize=14,
              leftcols=c("studlab"), leftlabs=c("Outcome"),
              plotwidth="3inch", sm="OR", sortvar=tnfb[["Effect"]],
              rightcols=c("effect","ci","pval"), rightlabs=c("OR","95%CI","P"),
              digits=3, digits.pval=2, scientific.pval=TRUE,
              common=FALSE, random=FALSE, print.I2=FALSE, print.pval.Q=FALSE, print.tau2=FALSE,
              addrow=TRUE, backtransf=TRUE, spacing=1.6,type.study="square",col.inside="black",col.square="black")

require(gap)
s <- chr_pos_a1_a2(1,c(123,321),letters[1:2],letters[2:1])
s
inv_chr_pos_a1_a2(s)
inv_chr_pos_a1_a2("chr1:123-A_B",seps=c(":","-","_"))

example(ci2ms)

gc.lambda <- function(x, logscale=FALSE, z=FALSE) {
  v <- x[!is.na(x)]
  n <- length(v)
  if (z) {
     obs <- v^2
     exp <- qchisq(log(1:n/n),1,lower.tail=FALSE,log.p=TRUE)
  } else {
    if (!logscale)
    {
      obs <- qchisq(v,1,lower.tail=FALSE)
      exp <- qchisq(1:n/n,1,lower.tail=FALSE)
    } else {
      obs <- qchisq(-log(10)*v,1,lower.tail=FALSE,log.p=TRUE)
      exp <- qchisq(log(1:n/n),1,lower.tail=FALSE,log.p=TRUE)
    }
  }

  lambda <- median(obs)/median(exp)
  return(lambda)
}

# A simplified version is as follows,
# obs <- median(chisq)
# exp <- qchisq(0.5, 1) # 0.4549364
# lambda <- obs/exp
# see also estlambda from GenABEL and qq.chisq from snpStats

# A related function

lambda1000 <- function(lambda, ncases, ncontrols)
  1 + (lambda - 1) * (1 / ncases + 1 / ncontrols)/( 1 / 1000 + 1 / 1000)

invnormal <- function(x) qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))

set.seed(12345)
Ni <- rpois(50, lambda = 4); table(factor(Ni, 0:max(Ni)))
y <- invnormal(Ni)
sd(y)
mean(y)
Ni <- 1:50
y <- invnormal(Ni)
mean(y)
sd(y)

alleles <- c("a","c","G","t")
revStrand(alleles)

library(gap)
search()
lsf.str("package:gap")
data(package="gap")$results
