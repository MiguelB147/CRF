library(MASS)
library(mgcv)



data("mcycle")

source("motorcycle-functions.R")

gam.fit <- gam(accel ~ s(times, k = 15, bs = "bs"), data = mcycle, optimizer = "efs")
gam.fit$sp
gam.fit$sig2


fit <- EstimatePenal(dim = 15, lambda.init = 100, step.control = F, type = "ps")
fit2 <- EstimatePenal(dim = 15, lambda.init = 100, step.control = F, type = "bs")
y.ps <- fit$X %*% fit$beta
y.bs <- fit2$X %*% fit2$beta

cbind(fitted(gam.fit), y.ps, y.bs)

plot(mcycle$times,mcycle$accel)
lines(sort(mcycle$times), fitted(gam.fit)[order(mcycle$times)])
lines(sort(mcycle$times), y.ps[order(mcycle$times)], col = "blue")
lines(sort(mcycle$times), y.bs[order(mcycle$times)], col = "red")
