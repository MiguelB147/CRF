library(mvtnorm)
library(copula)

rho <- 0.5
set.seed(123)
sim <- rCopula(100, normalCopula(rho, dim = 2))

u <- sim[,1]
v <- sim[,2]

C00 <- pmvnorm(upper = c(qnorm(u)), qnorm(v), corr = matrix(rho,1,1,rho), byrow=TRUE)
C10 <- pnorm((qnorm(v) - rho*qnorm(u))/sqrt(1-rho^2))
C01 <- pnorm((qnorm(u) - rho*qnorm(v))/sqrt(1-rho^2))
C11 <- (1/sqrt(1-rho^2))*(1/dnorm(qnorm(u)))*dnorm((qnorm(u) - rho*qnorm(v))/sqrt(1-rho^2))
CRF <- C00*C11/(C01*C10)

set.seed(123)
