ccsize <- function(n,q,pD,p1,theta,alpha,beta=0.2,power=FALSE,verbose=FALSE)
{
   p2 <- 1 - p1
   if(power)
# equation (5)
   {
     z_alpha <- qnorm(alpha)
     z <- z_alpha + sqrt(n) * theta * sqrt(p1 * p2 / (1 / pD + (1 / q - 1)))
     invisible(pnorm(z))
   }
   else
# equation (6)
   {
     nb <- -999
     z_alpha <- qnorm(alpha, lower.tail=FALSE)
     z_beta <- qnorm(1-beta)
     theta_lon <- (z_alpha + z_beta) / sqrt(p1 * p2 * pD)
     d <- (theta / theta_lon)^2 - (1 - pD) / n
     if (any(d <= 0) & verbose) cat("bad hazard ratio =", exp(theta), "\n")
     else 
     {
       nb <- ceiling(pD / d)
       nb [nb > n] <- -999
       if (any(nb > n) & verbose) cat("bad subcohort size", nb, "\n")
     }
     invisible(nb)
   }
}
