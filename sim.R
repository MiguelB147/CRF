

library(splines)
library(Rcpp)
library(rootSolve)

source('hunan-functions.R')
sourceCpp('test.cpp')

K <- 1000
dim <- 7

plot.grid <- expand.grid(seq(0.25,1.5,by=0.25), seq = seq(0.1,2.2,by = 0.05))
names(plot.grid) <- c("time1","time2")
plot.grid$true <- theta.frank(plot.grid$time1,plot.grid$time2)


B <- 20
CRF <- matrix(NA, nrow = nrow(plot.grid), ncol = B)
Lambda <- matrix(NA, nrow = B, ncol = 2)
for (i in 1:B) {
  
  datalist <- SimData(K = K, unif.ub = NULL)
  
  if (i == 1) {
    beta.start <- rep(1, dim^2)
    lambda.init <- c(10,10)
  } else {
    beta.start <- fit$beta
    lambda.init <- fit$lambda
  }
  
  fit <- EstimatePenal2(datalist = datalist, dim = dim, start = beta.start, type = "ps", lambda.init = lambda.init, scale = FALSE, quantile = FALSE)
  
  CRF[,i] <- WoodTensor.predict(plot.grid$time1, plot.grid$time2, fit)
  Lambda[i,] <- fit$lambda

}