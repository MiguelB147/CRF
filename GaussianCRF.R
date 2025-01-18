library(mvtnorm)
library(copula)

rho <- 0.5
set.seed(123)
sim <- rCopula(100, normalCopula(rho, dim = 2))

u <- sim[,1]
v <- sim[,2]

normalCRF(u,v,rho)

normalCRF <- function (u, v, rho) {
  
  sim <- cbind(u,v)
  
  C00 <- apply(sim, 1, function (x) pmvnorm(lower = c(-Inf, -Inf), upper = x, corr = matrix(c(1,rho,rho,1), byrow=TRUE, ncol=2))) # TODO Ga dit na in de code van de masterproef
  C10 <- pnorm((qnorm(v) - rho*qnorm(u))/sqrt(1-rho^2))
  C01 <- pnorm((qnorm(u) - rho*qnorm(v))/sqrt(1-rho^2))
  C11 <- (1/sqrt(1-rho^2))*(1/dnorm(qnorm(u)))*dnorm((qnorm(u) - rho*qnorm(v))/sqrt(1-rho^2))
  CRF <- C00*C11/(C01*C10)
  
  return(CRF)
}



set.seed(123)
