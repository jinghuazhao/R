#' Score statistics for association of traits with haplotypes
#'
#' @param y Vector of trait values. For  trait.type  =  "binomial",  y  must have values of 1 for event, 0 for no event.
#' @param geno Matrix of alleles, such that each locus has a pair of adjacent columns of alleles, and the order of columns corresponds to the order of loci on a chromosome. If there are K loci, then ncol(geno) = 2*K. Rows represent alleles for each subject.
#' @param trait.type Character string  defining  type  of  trait, with values of "gaussian", "binomial", "poisson", "ordinal".
#' @param offset Vector of offset when trait.type = "poisson".
#' @param x.adj Matrix of non-genetic covariates used to adjust the score statistics. Note that intercept should not be included, as it will be added in this function.
#' @param skip.haplo Skip score statistics for haplotypes with frequencies < skip.haplo.
#' @param locus.label Vector of labels for loci, of length K (see definition of geno matrix).
#' @param miss.val Vector of codes for missing values of alleles.
#' @param n.sim Number of simulations for empirical p-values.  If n.sim=0, no empirical p-values are computed.
#' @param method method of haplotype frequency estimation, "gc" or "hap".
#' @param id an added option which contains the individual IDs.
#' @param handle.miss flag to handle missing genotype data, 0=no, 1=yes.
#' @param mloci maximum number of loci/sites with missing data to be allowed in the analysis.
#' @param sexid flag to indicator sex for data from X chromosome, i=male, 2=female.
#'
#' @details
#' Compute score statistics to evaluate the association of a trait with haplotypes, when linkage phase is unknown and diploid marker 
#' phenotypes are observed among unrelated subjects. For now, only autosomal loci are considered. This package haplo.score
#' which this function is based is greatly acknowledged.
#'
#' @export
#' @return
#' List with the following components:
#' - score.global Global statistic to test association of trait with haplotypes that have frequencies >= skip.haplo.
#' - df Degrees of freedom for score.global.
#' - score.global.p P-value of score.global based on chi-square distribution, with degrees of freedom equal to df.
#' - score.global.p.sim P-value of score.global based on simulations (set equal to NA when n.sim=0).
#' - score.haplo Vector of score statistics for individual haplotypes that have frequencies >= skip.haplo.
#' - score.haplo.p Vector of p-values for score.haplo, based on a chi-square distribution with 1 df.
#' - score.haplo.p.sim Vector of p-values for score.haplo, based on  simulations (set equal to NA when n.sim=0).
#' - score.max.p.sim P-value  of  maximum  score.haplo, based on simulations (set equal to NA when n.sim=0).
#' - haplotype Matrix of hapoltypes  analyzed.  The ith row of haplotype corresponds to the ith item of score.haplo, score.haplo.p, and score.haplo.p.sim.
#' - hap.prob Vector of haplotype probabilies, corresponding to the haplotypes in the matrix haplotype.
#' - locus.label Vector of labels for loci, of length K (same as input argument).
#' - n.sim Number of simulations.
#' - n.val.global Number of valid simulated global statistics.
#' - n.val.haplo Number of valid simulated score statistics (score.haplo) for individual haplotypes.
#'
#' @details This is a version which substitutes haplo.em.
#'
#' @references
#' \insertRef{schaid02}{gap}
#'
#' @examples
#' \dontrun{
#' data(hla)
#' y<-hla[,2]
#' geno<-hla[,3:8]
#' # complete data
#' hap.score(y,geno,locus.label=c("DRB","DQA","DQB"))
#' # incomplete genotype data
#' hap.score(y,geno,locus.label=c("DRB","DQA","DQB"),handle.miss=1,mloci=1)
#' unlink("assign.dat")
#'
#' ### note the differences in p values in the following runs
#' data(aldh2)
#' # to subset the data since hap doesn't handle one allele missing
#' deleted<-c(40,239,256)
#' aldh2[deleted,]
#' aldh2<-aldh2[-deleted,]
#' y<-aldh2[,2]
#' geno<-aldh2[,3:18]
#' # only one missing locus
#' hap.score(y,geno,handle.miss=1,mloci=1,method="hap")
#' # up to seven missing loci and with 10,000 permutations
#' hap.score(y,geno,handle.miss=1,mloci=7,method="hap",n.sim=10000)
#'
#' # hap.score takes considerably longer time and does not handle missing data
#' hap.score(y,geno,n.sim=10000)
#' }
#'
#' @keywords models regression

hap.score <- function(y, geno, trait.type="gaussian",
                      offset = NA, x.adj = NA, skip.haplo=.005,
                      locus.label=NA, miss.val=0, n.sim=0, method="gc", id=NA,
                      handle.miss=0, mloci=NA, sexid=NA)
{
  trait.int <- charmatch(trait.type, c("gaussian", "binomial", "poisson", "ordinal"))
  if(is.na(trait.int)) stop("Invalid trait type")
  if(trait.int == 0)   stop("Ambiguous trait type")
  if(length(y)!=nrow(geno)) stop("Dims of y and geno are not compatible")
  n.loci <- ncol(geno)/2
  if(n.loci != (floor(ncol(geno)/2))) stop("Odd number of cols of geno")
  if(handle.miss==0)
  {
    miss <- apply(is.na(geno),1,any)
    if(!all(is.na(miss.val))) {
       for(mval in miss.val){
          miss <- miss | apply(geno==mval, 1, any)
       }
    }
  }
  else
  {
    if(is.na(mloci)) stop("Maximum number of missing loci (mloci) not specified")
    nmiss <- apply(is.na(geno),1,sum)
    if(!all(is.na(miss.val))) {
       for(mval in miss.val) {
          nmiss <- nmiss + apply(geno==mval, 1, sum)
       }
    }
    if(mloci<0 | mloci >= n.loci) stop("Invalid control for number of missing loci")
    miss <- rep(F, length(y))
    for(i in 1:length(y)) if(nmiss[i] > mloci*2) miss[i] <- T
  }
  adjusted <- T
  if( all(is.na(x.adj))) adjusted <- F
  if(adjusted){
    x.adj <- as.matrix(x.adj)
    if(nrow(x.adj)!=length(y)) stop("Dims of y and x.adj are not compatible")
  }
  miss <- miss | is.na(y) 
  if(adjusted) miss <- miss| apply(is.na(x.adj),1,any)
  if(trait.int==3) {
     if(all(is.na(offset))) stop("Missing offset")
     miss <- miss | is.na(offset)
     offset <- offset[!miss]
  }
  y <- as.numeric(y[!miss])
  geno <- geno[!miss,]
  if(adjusted) x.adj <- x.adj[!miss,,drop=F]
  if(trait.int==2) {
    if(!all(y==1|y==0)) stop("Invalid y values")
    if(all(y==1) | all(y==0)) stop("No variation in y values")
  }
  if(trait.int==4){
     y <- factor(y)
     y.lev <- levels(y)
     y <- as.numeric(y)
     if(max(y) < 3) stop("Less than 3 levels for y values")
  }
  n.subj <- length(y)
  if(all(is.na(id))) id <- 1:n.subj
  method.id<-charmatch(method, c("gc", "hap", "phase"))
  if(is.na(method.id)) stop("Invalid selection of method")
  if(method.id == 0)   stop("Ambiguous method")
  else if(method.id==1) haplo <- gc.em(data=geno, locus.label, converge.eps=0.00001, maxiter=5000, handle.miss=handle.miss, miss.val=miss.val)
  else if(method.id==2) haplo <- hap.em(id, data=geno, locus.label, converge.eps=0.00001, maxiter=5000, miss.val=miss.val)
  if(method.id<3 & !haplo$converge) stop("EM for haplo failed to converge")
  hap1 <- haplo$hap1code
  hap2 <- haplo$hap2code
  indx <- haplo$indx.subj
  post <- haplo$post
  nreps <- as.vector(haplo$nreps)
  uhap<-haplo$uhap
  which.haplo<-haplo$hap.prob>=skip.haplo
  uhap<-uhap[which.haplo]
  x <- outer(hap1,uhap,"==") + outer(hap2,uhap,"==")
  n.x <- ncol(x)
  x.post<-matrix(rep(NA, n.subj * n.x), ncol=n.x)
  for(j in 1:n.x){
     x.post[,j] <- tapply(x[,j]*post, indx, sum)
  }
  if(trait.int <= 3){ 
    if(!adjusted){
       mu <- switch(trait.int, mean(y), mean(y), sum(y)/sum(offset) )
       a  <- switch(trait.int, var(y), 1, 1)
       x.adj <- matrix(rep(1,n.subj),ncol=1)
     }
     if(adjusted){
        reg.out <- glm(y ~ x.adj, family=trait.type)
        x.adj <- cbind(rep(1,n.subj),x.adj)
        mu <- reg.out$fitted.values
        a  <- switch(trait.int,sum(reg.out$residuals^2)/reg.out$df.residual,1, 1)
      }
     v <- switch(trait.int, 1/a, mu*(1-mu), mu )
     tmp <- hap.score.glm(y, mu, a, v, x.adj, nreps, x.post, post, x)
     u.score <- tmp$u.score
     v.score <- tmp$v.score
   }
   if(trait.int ==4) {
      if(adjusted){
         for(p in c("rms")) {
            if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
               if (!requireNamespace(p, quietly = TRUE))
               warning(paste("This function needs package `", p, "' to be fully functional; please install", sep=""))
            }
         }
         reg.out <- rms::lrm(y ~ x.adj)
         K <- max(y)
         n.xadj <- ncol(x.adj)
         alpha <- reg.out$coef[1:(K-1)]
         beta <- reg.out$coeff[K:(K-1 + n.xadj)]

         tmp <- hap.score.podds(y, alpha, beta, x.adj, nreps, x.post, post, x)
       }
      if(!adjusted){
         tbl <- table(y)
         s <- 1- (cumsum(tbl)-tbl)/n.subj
         alpha <-  - log((1-s[-1])/s[-1])
         tmp <- hap.score.podds(y, alpha, beta=NA, x.adj=NA, nreps, x.post, post, x)
       }
      u.score <- tmp$u.score
      v.score <- tmp$v.score
    }
   tmp <- haplo.stats::Ginv(v.score)
   df <- tmp$rank
   g.inv <- tmp$Ginv
   score.global <- u.score%*% g.inv %*%u.score
   score.haplo <- u.score / sqrt(diag(v.score))
   score.max <-  max(score.haplo^2)
   if(n.sim==0){
      score.global.p.sim <- NA
      score.haplo.p.sim <- rep(NA,length(score.haplo))
      score.max.p.sim <- NA
      n.val.global <- NA
      n.val.haplo <- NA
   }
   if(n.sim > 0){
      score.global.rej <- 0
      score.haplo.rej  <- rep(0,length(score.haplo))
      score.max.rej    <- 0
      n.val.global <- 0
      n.val.haplo <- 0
      if(trait.int<=3){
         mu.rand <- mu
         v.rand <- v  
       }
      for(i in 1:n.sim){
         rand.ord <- order(runif(n.subj))
         if(trait.int <=3){ 
           if(adjusted){
              mu.rand <- mu[rand.ord]
              v.rand <- switch(trait.int, v, v[rand.ord], v[rand.ord])
            }
           tmp <- hap.score.glm(y[rand.ord], mu.rand, a, v.rand, 
                                x.adj[rand.ord,], nreps, x.post, post, x)
         }
         if(trait.int ==4){
            if(adjusted){             
               tmp <- hap.score.podds(y[rand.ord], alpha, beta, 
                               x.adj[rand.ord,,drop=F],nreps, x.post, post, x)
            }
            if(!adjusted) {
               tmp <- hap.score.podds(y[rand.ord], alpha, beta=NA, 
                               x.adj=NA,nreps, x.post, post, x)
             }

          }
         u.score <- tmp$u.score
         v.score <- tmp$v.score
         tmp <- haplo.stats::Ginv(v.score)
         g.inv <- tmp$Ginv  
         score.global.sim <- u.score %*% g.inv %*% u.score
         score.haplo.sim  <- (u.score / sqrt(diag(v.score)))^2
         score.max.sim <- max(score.haplo.sim)
         if(!is.na(score.global.sim)) {
            n.val.global <- n.val.global +1
            if(score.global.sim >= score.global) score.global.rej <- score.global.rej +1
          }
         if(!any(is.na(score.haplo.sim))){
            n.val.haplo <- n.val.haplo + 1
            score.haplo.rej <- score.haplo.rej +
                               ifelse(score.haplo.sim >= score.haplo^2, 1, 0)
            if(score.max.sim >= score.max) score.max.rej <- score.max.rej +1
          }
      }
      score.global.p.sim <- score.global.rej /  n.val.global
      score.haplo.p.sim <- score.haplo.rej / n.val.haplo
      score.max.p.sim <- score.max.rej / n.val.haplo
    }
   score.global.p <- 1 - pchisq(score.global,df)
   score.haplo.p <- 1-pchisq(score.haplo^2,1)
   if(all(is.na(locus.label))) {
      locus.label<- paste("loc-",1:n.loci,sep="")
    }
   obj <- (list(score.global=score.global, df=df,score.global.p=score.global.p,
       score.global.p.sim=score.global.p.sim,
       score.haplo=score.haplo,score.haplo.p=score.haplo.p,
       score.haplo.p.sim=score.haplo.p.sim,
       score.max.p.sim=score.max.p.sim,
       haplotype=haplo$haplotype[which.haplo,],
       hap.prob=haplo$hap.prob[which.haplo],
       locus.label=locus.label,
       n.sim=n.sim, n.val.global=n.val.global, n.val.haplo=n.val.haplo))
   class(obj) <- "hap.score"
   return(obj)
}
 
hap.score.glm <- function(y,mu,a,v,x.adj,nreps,x.post, post, x)
{
   u.mtx  <- (y-mu)*x.post / a
   u.score <- apply(u.mtx,2,sum)

   # Var matrix for x.adj covariates
   v.11 <- t(x.adj * v) %*% x.adj

   # Var matrix for covar(x.adj, x.post)
   v.21 <- t(x.post) %*% (x.adj * v)

   # Var matrix for haplo scores
   res <- ( (y - mu)/a ) ^2
   t1 <- rep( (v-res) ,nreps) * post
   v.22 <- t(x*t1) %*% x + t(u.mtx) %*% u.mtx

   # Var matrix for haplo scores, adjusted for x.adj
   v.score <- v.22 - v.21 %*% solve(v.11) %*% t(v.21) 
   return(list(u.score=u.score, v.score=v.score))
}

hap.score.podds <- function(y, alpha, beta=NA, x.adj=NA, nreps, x.post,  post, x)
{
###################################################################
#
# If U=c(u.a, u.e, u.g), where 
#   u.a = score for alpha's
#   u.e = score for unambiguous (x.adj) covariates
#   u.g = score for ambiguous haplotypes
#
# Then the upper triangle of Var(U) can be partitioned as
#
#          | v.aa   v.ae   v.ag |   |           |
#   V(U) = |        v.ee   v.eg | = | v.11 v.12 |
#          |               v.gg |   |      v.gg |
#
# where v.12 is composed of v.aa, v.ae, v.ee
#       v.12 is composed of v.ag, v.eg
#
# and Var(u.g) = v.gg - v.12 * v.12(inv) * t(v.12)
#
# The following computes each of the submatrices as needed
# to determine u.g and Var(u.g)
#
##################################################################
adjusted <- T
if(any(is.na(x.adj))) adjusted <- F
if(adjusted) n.xadj <- ncol(x.adj)
n.x <- ncol(x)
K <- max(y)

# to make suscripting easier, append Inf to front of alpha,
# as place-holder for alpha[1] = Inf
alpha <- c(Inf, alpha)
if(adjusted){
   s   <- ifelse(y==1, 1, 1/(1 + exp(-(alpha[y  ] + x.adj %*% beta ))) )
   s.p <- ifelse(y==K, 0, 1/(1 + exp(-(alpha[y+1] + x.adj %*% beta ))) )
 }
if(!adjusted){
   s   <- ifelse(y==1, 1, 1/(1 + exp(-(alpha[y  ]  ))) )
   s.p <- ifelse(y==K, 0, 1/(1 + exp(-(alpha[y+1]  ))) )
 }
w1 <- (s*(1-s) - s.p*(1-s.p))/(s - s.p)
u.mtx <- w1 * x.post
u.score <- apply(u.mtx,2,sum)

#  compute information matrix for alpha-beta (v.ab) and alpha-alpha (v.aa)
tmp1 <- (s   + s.p^2 - 2*s*s.p)*s.p*(1-s.p)/(s-s.p)^2
tmp2 <- (s.p +   s^2 - 2*s*s.p)*s*(1-s)/(s-s.p)^2
tmp3 <- s.p*(1-s.p)*s*(1-s)/(s-s.p)^2

v.ag <- matrix(rep(0, (K-1)*n.x), ncol=n.x)
if(adjusted) v.ae <- matrix(rep(0, (K-1)*n.xadj), ncol=n.xadj)
v.aa <- matrix(rep(0,(K-1)^2),ncol=(K-1))

n.subj <- length(y)
for(j in 2:K){
   wt <- rep(0,n.subj)
   wt <- ifelse(y==(j-1), (tmp1 - tmp3), wt)
   wt <- ifelse(y==j, (tmp2 - tmp3), wt)
   v.ag[(j-1),] <- apply(wt * x.post, 2,sum)
   if(adjusted) v.ae[(j-1),] <-  apply(wt * x.adj, 2,sum)
   v.aa[(j-1),(j-1)] <- sum(tmp1[y==(j-1)]) + sum(tmp2[y==j])   
   if(j < K) v.aa[(j-1), j] <- -sum(tmp3[y==j])
 }

# fill in lower tri of v.aa to make it symmetric
v.aa <- v.aa + t( (col(v.aa) > row(v.aa))*v.aa )

# Louis' method for v.gg
w2 <- s*(1-s) + s.p*(1-s.p)
t1 <- rep( (w2 - w1^2), nreps) * post
v.gg <- t(x*t1) %*% x + t(u.mtx) %*% u.mtx

if(adjusted){
   v.ee <- t(w2*x.adj) %*% x.adj
   v.eg <- t(w2*x.adj) %*% x.post
   v.11 <- rbind( cbind(v.aa, v.ae), cbind(t(v.ae),v.ee) )
   v.12 <- rbind(v.ag,v.eg)
   v.score <- v.gg - t(v.12) %*% solve(v.11) %*% v.12
 }
if(!adjusted){
   v.score <- v.gg - t(v.ag) %*% solve(v.aa) %*% v.ag
 }

return(list(u.score=u.score, v.score=v.score))
}

#' Plot haplotype frequencies versus haplotype score statistics
#'
#' Method function to plot a class of type hap.score
#'
#' @param x The object returned from hap.score (which has class hap.score).
#' @param ... Optional arguments.
#'
#' @export
#' @return
#' Nothing is returned.
#'
#' This is a plot method function used to plot haplotype frequencies on
#' the x-axis and haplotype-specific scores on the y-axis. Because
#' hap.score is a class, the generic plot function 
#' can be used, which in turn calls this plot.hap.score function.
#'
#' @references
#' Schaid DJ, Rowland CM, Tines DE, Jacobson RM,  Poland  GA (2002)
#' Score tests for association of traits with haplotypes when
#' linkage phase is ambiguous. Amer J Hum Genet 70:425-34
#'
#' @seealso [`hap.score`]
#'
#' @examples
#' \dontrun{
#' save <- hap.score(y, geno, trait.type = "gaussian")
#' 
#' # Example illustrating generic plot function:
#' plot(save)
#'
#' # Example illustrating specific method plot function:
#' plot.hap.score(save)
#' }
#'
#' @keywords hplot

plot.hap.score <- function(x, ...){
   plot(x$hap.prob, x$score.haplo, xlab="Haplotype Frequency", 
        ylab="Haploltype Score Statistic", ...)
   invisible()
}

#' Print a hap.score object
#'
#' Method function to print a class of type hap.score
#'
#' @param x The object returned from hap.score (which has class hap.score).
#' @param ... Optional argunents.
#'
#' @export
#' @return Nothing is returned.
#' 
#' This is a print method function used to print information from
#' hap.score class, with haplotype-specific information given in a
#' table. Because hap.score is a class, the generic print function 
#' can be used, which in turn calls this print.hap.score function.
#'
#' @references
#' Schaid DJ, Rowland CM, Tines DE, Jacobson RM, Poland  GA (2002)
#' Score tests for association of traits with haplotypes when
#' linkage phase is ambiguous. Amer J Hum Genet 70:425-34
#'
#' @seealso [`hap.score`]
#'
#' @examples
#' \dontrun{
#' save <- hap.score(y, geno, trait.type = "gaussian")
#'
#' # Example illustrating generic print function:
#' print(save)
#'
#' # Example illustrating specific method print function:
#' print.hap.score(save)
#' }
#'
#' @keywords print

print.hap.score <- function(x, ...){

# print of global score stats:
   cat("\nGlobal Score Statistics\n\n")
   cat(paste("global-stat = ",round(x$score.global,5),", df = ",x$df,
         ", p-val = ",round(x$score.global.p,5),sep=""))
   if(x$n.sim>0) cat(", sim. p-val = ",x$score.global.p.sim,"\n\n")
   if(x$n.sim>0) cat("max-stat sim. p-val = ",x$score.max.p.sim)
   cat("\n\n")

# create table for haplotype specific stats:
   tbl <- cbind(x$haplotype,round(x$hap.prob,5),round(x$score.haplo,5),
          round(x$score.haplo.p,5))
   if(x$n.sim>0) tbl <- cbind(tbl,x$score.haplo.p.sim)
   ord <- order(x$score.haplo)
   tbl <- tbl[ord,]
   if(x$n.sim == 0) dimnames(tbl) <- list(NULL,c(x$locus.label,"Hap-Freq",
                  "Hap-Score","p-val"))
   if(x$n.sim > 0) dimnames(tbl) <- list(NULL,c(x$locus.label,"Hap-Freq",
                  "Hap-Score","p-val","sim p-val"))
   cat("Haplotype-specific Scores\n\n")
   print(tbl,quote=F)
   cat("\n\n")
   invisible()
}

# 13-9-2003 start to implement
# 14-9-2003 in shape
# 21-9-2003 start extensive checking
# 23-9-2003 rewrite interface to genecounting
# 26-9-2003 done with successful use of by and order
# 17-10-2003 start to implement missing genotype code
# 18-9-2004 to fix S3 class hap.score
