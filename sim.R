

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
  
  fit <- EstimatePenal2(datalist = datalist, dim = dim)
  
  CRF[,i] <- WoodTensor.predict(plot.grid$time1, plot.grid$time2, fit)
  Lambda[i,] <- fit$lambda

}

for (j in 1:length(Beta)) {
  
  CRF <- apply(Beta[[j]], 2, function(x) WoodTensor.predict2(plot.grid$time1, plot.grid$time2,x))
  
  pdf(paste0("sim10beta",j,".pdf"))
  
  par(mfrow = c(2,3))
  
  for (i in unique(plot.grid$time1)) {
    plottext <- paste0("t[1] == ", i)
    plot(x = plot.grid$time2[plot.grid$time1 == i],
         y = plot.grid$true[plot.grid$time1 == i],
         type = 'l', lwd = 2, col = "grey",
         ylab = "CRF", xlab = expression(t[2]), main = parse(text = plottext),
         ylim=c(0,7))
    for (k in 1:ncol(CRF)){
      lines(x = plot.grid$time2[plot.grid$time1 == i], y = CRF[,k][plot.grid$time1 == i], col = k)
    }

  }

  par(mfrow = c(1,1))
  dev.off()
}


plot.grid$CRF <- rowMeans(CRF)
par(mfrow = c(2,3))

for (i in unique(plot.grid$time1)) {
  plottext <- paste0("t[1] == ", i)
  plot(x = plot.grid$time2[plot.grid$time1 == i],
       y = plot.grid$true[plot.grid$time1 == i],
       type = 'l', lwd = 2, col = "grey",
       ylab = "CRF", xlab = expression(t[2]), main = parse(text = plottext),
       ylim=c(0,7))
  lines(x = plot.grid$time2[plot.grid$time1 == i], y = plot.grid$CRF[plot.grid$time1 == i])
  
}

par(mfrow = c(1,1))



{
  dim = 7
  lambda1 = 5
  lambda2 = 3000
  obj1 <- WoodSpline(datalist$X[,1], dim = dim, type = "bs")
  obj2 <- WoodSpline(datalist$X[,2], dim = dim, type = "bs")
  
  S <- WoodPenalty(obj1, obj2)
  Sl = lambda1*S[[1]] + lambda2*S[[2]]
  
  start = rep(1,dim^2)
  beta <- multiroot(Score2, start = start, rtol = 1e-10, X1 = obj1$X, X2 = obj2$X, Sl = Sl, datalist = datalist)$root
  
  plot.grid <- expand.grid(seq(0.25,1.5,by=0.25), seq = seq(0.1,2.2,by = 0.05))
  names(plot.grid) <- c("time1","time2")
  plot.grid$true <- theta.frank(plot.grid$time1,plot.grid$time2)
  
  X1 <- splines::splineDesign(obj1$knots, plot.grid$time1, ord = 4)
  X2 <- splines::splineDesign(obj2$knots, plot.grid$time2, ord = 4)
  CRF <- exp(row.kronecker(X1, X2) %*% beta)
  plot.grid$CRF <- WoodTensor.predict(plot.grid$time1,plot.grid$time2, fit7ps)
  par(mfrow = c(2,3))
  for (i in unique(plot.grid$time1)) {
    plottext <- paste0("t[1] == ", i)
    plot(x = plot.grid$time2[plot.grid$time1 == i],
         y = plot.grid$true[plot.grid$time1 == i],
         type = 'l', lwd = 2, col = "grey",
         ylab = "CRF", xlab = expression(t[2]),
         main = parse(text = plottext),
         ylim=c(0,7))
    lines(x = plot.grid$time2[plot.grid$time1 == i], y = CRF[plot.grid$time1 == i])
    lines(x = plot.grid$time2[plot.grid$time1 == i], y = plot.grid$CRF[plot.grid$time1 == i], col = "red")
    text(c(6), paste0("lambda1 = ", lambda1))
    text(c(5), paste0("lambda2 = ", lambda2))
    
  }
  par(mfrow = c(1,1))
}

plot_ly(showscale = FALSE) %>%
  add_surface(z = ~xtabs(plot.grid$true ~ plot.grid$time1 + plot.grid$time2)) %>% 
  add_surface(z = ~xtabs(plot.grid$CRF ~ plot.grid$time1 + plot.grid$time2), opacity = 0.8, colorscale = list(c(0,1),c("rgb(255,107,184)","rgb(128,0,64)")))


library(copula)

plot.grid <- expand.grid(seq(0.25,5,by=0.25), seq = seq(0.1,5,by = 0.05))
names(plot.grid) <- c("time1","time2")

plot.grid$CRF <- theta.mix(plot.grid$time1, plot.grid$time2)
plot.grid$CRF2 <- theta.mix2(plot.grid$time1, plot.grid$time2)


par(mfrow = c(2,3))
for (i in unique(plot.grid$time1)) {
  plottext <- paste0("t[1] == ", i)
  plot(x = plot.grid$time2[plot.grid$time1 == i],
       y = plot.grid$CRF2[plot.grid$time1 == i],
       type = 'l', lwd = 2, col = "black",
       ylab = "CRF", xlab = expression(t[2]),
       main = parse(text = plottext))
  lines(x = plot.grid$time2[plot.grid$time1 == i], y = plot.grid$CRF[plot.grid$time1 == i], lwd = 2, col = "grey")
  
}
par(mfrow = c(1,1))
