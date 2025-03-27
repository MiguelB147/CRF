

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


B <- 5
CRF <- matrix(NA, nrow = nrow(plot.grid), ncol = B)
Lambda <- matrix(NA, nrow = B, ncol = 2)
for (i in 1:B) {
  
  datalist <- SimData(K = K, unif.ub = NULL)
  
  if (i == 1) {
    beta.start <- rep(1, dim^2)
  } else {
    beta.start <- fit$beta
  }
  
  fit <- EstimatePenal2(datalist = datalist, dim = dim, start = beta.start, type = "bs", scale = FALSE, quantile = FALSE)
  
  CRF[,i] <- WoodTensor.predict(plot.grid$time1, plot.grid$time2, fit)
  Lambda[i,] <- fit$lambda

}

plot.grid$CRF <- rowMeans(CRF)

pdf("sim5ps.pdf")
par(mfrow = c(2,3))
for (i in unique(plot.grid$time1)) {
  
  plottext <- paste0("t[1] == ", i)
  plot(x = plot.grid$time2[plot.grid$time1 == i],
       y = plot.grid$true[plot.grid$time1 == i],
       type = 'l', lwd = 2, col = "grey",
       ylab = "CRF", xlab = expression(t[2]), main = parse(text = plottext),
       ylim=c(0,7))
  lines(x = plot.grid$time2[plot.grid$time1 == i], y = plot.grid$CRF[plot.grid$time1 == i], col = 1, lwd = 2)
}
par(mfrow = c(1,1))
dev.off()
