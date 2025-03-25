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
