library(MASS)
library(mgcv)

data("mcycle")

source("motorcycle-functions.R")

gam.fit <- gam(accel ~ s(times, k = 20, bs = "bs"), data = mcycle, optimizer = "efs")
gam.fit$sp
gam.fit$sig2

test <- WoodSpline(t = mcycle$times, dim = 10, type = "ps")

crossprod(test$D)

fit <- EstimatePenal(dim = 20, lambda.init = 1, step.control = F, scale = T, type = "ps", quantile = FALSE)
fit2 <- EstimatePenal(dim = 20, lambda.init = 1, step.control = F, scale = T, type = "bs", quantile = TRUE)
y.ps <- fit$X %*% fit$beta
y.bs <- fit2$X %*% fit2$beta

cbind(fitted(gam.fit), y.ps, y.bs)

pdf("mcycle.pdf")
plot(mcycle$times,mcycle$accel)
lines(sort(mcycle$times), fitted(gam.fit)[order(mcycle$times)])
lines(sort(mcycle$times), y.ps[order(mcycle$times)], col = "blue")
lines(sort(mcycle$times), y.bs[order(mcycle$times)], col = "red")
legend("topleft", legend = c("mgcv - ps", "uniform", "quantile"), lty = rep(1,3), col = c("black", "blue", "red"))
dev.off()


x = mcycle$times
n = length(x)
k = round(quantile(unique(x), seq(0,1,length = 10)),1)
nk = 10
h <- diff(k)

pos <- vector("numeric", length = n)
for (i in 1:n) {
  xi <- x[i]
  j = nk
  while (xi <= k[j]) {
    if (j > 1) {
      j <- j-1
    } else break
  }
  pos[i] <- j
}

am <- (k[pos+1] - x)/h[pos]
ap <- (x - k[pos])/h[pos]
cm <- ( (k[pos+1] - x)^3/h[pos] - h[pos]*(k[pos+1] - x))/6
cp <- ( (x - k[pos])^3/h[pos] - h[pos]*(x - k[pos]))/6

D <- B <- matrix(0, ncol = 8, nrow = nk)
for (i in 1:(nk-2)) {
  D[i,i] <- 1/h[i]
  D[i,i+1] <-  -1/h[i] - 1/h[i+1]
  D[i,i+2] <-  1/h[i+1]
  
  B[i,i] <- (h[i] + h[i+1])/3
  if (i < nk-2) {
    B[i,i+1] <- B[i+1,i] <- h[i+1]/6
  }
}
