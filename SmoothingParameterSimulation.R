library(splines)
library(doParallel)
library(Rcpp)
library(MASS)

setwd("~/Library/CloudStorage/GoogleDrive-miguel-angel.beynaerts@student.uhasselt.be/Mijn Drive/CRF simulations")
# setwd("H:/My Drive/CRF simulations")
source('hunan-functions.R')
sourceCpp('test.cpp')

nsim <- 500
nsets <- 100

degree = 2
df = 7
K <- 1000
unif.ub <- 5 # 5 = 20% censoring, 2.3 = 40% censoring


lambda.grid <- expand.grid(Lambda1 = seq(0.1, 2, 0.1), Lambda2 = seq(0.1, 2, 0.1))



ll <- ll.avg <- c()
betas <- matrix(0, ncol = df^2, nrow = nsim)
lambda.opt <- matrix(0, nrow = nsets, ncol = 2)
dat.list <- list()

set.seed(123)
for (k in 1:nsets) {
  
  dat.list[[k]] <- SimData(K = K, df = df, degree = degree, unif.ub = unif.ub)

  for (i in 1:nrow(lambda.grid)) {
    
    S.lambda <- ginv(lambda.grid[i,1]*Srow(df) + lambda.grid[i,2]*Scol(df))
    betas <- mvrnorm(nsim, mu = rep(0,df^2), Sigma = S.lambda)
    
    for (j in 1:nsim) {
      ll[j] <- as.integer(loglikPenal(coef.vector = betas[j,],
                                      degree = degree,
                                      df = df,
                                      datalist = dat.list[[k]],
                                      lambda = lambda.grid[i,]))
    }
    ll.avg[i] <- mean(ll)
  }
  
  lambda.opt[k,] <- lambda.grid[which.min(ll.avg),]
}







lambda.grid$LogLikelihood <- -ll.avg
ggplot(lambda.grid, aes(Lambda1, Lambda2, fill = LogLikelihood)) + geom_tile() + theme_classic() + xlab(expression(lambda[1])) + ylab(expression(lambda[2]))
