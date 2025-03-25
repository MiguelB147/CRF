
data(mtcars)
head(mtcars)

library(mgcv)
library(splines)

source("cars-functions.R")


fit.gam <- gam(mpg ~ s(hp, k = 20, bs = "bs"), data = mtcars, optimizer = 'efs')
c(fit.gam$coef, fit.gam$sig2)
fit.gam$sp
fit.gam$sig2

test <- EstimatePenal(dim = 10, lambda.init = 1, step.control = F, type = "ps", scale = T, quantile = FALSE)
test2 <- EstimatePenal(dim = 20, lambda.init = 1, step.control = F, type = "bs", scale = T, quantile = FALSE)
y.ps <- test$X %*% test$beta
y.bs <- test2$X %*% test2$beta

cbind(fitted(fit.gam), y.bs, y.ps)

pdf("mtcars.pdf")
plot(mtcars$hp,mtcars$mpg)
lines(sort(mtcars$hp), fitted(fit.gam)[order(mtcars$hp)])
lines(sort(mtcars$hp), y.ps[order(mtcars$hp)], col = "blue")
lines(sort(mtcars$hp), y.bs[order(mtcars$hp)], col = "red")
legend("topright", legend = c("mgcv - ps", "uniform", "quantile"), lty = rep(1,3), col = c("black", "blue", "red"))
dev.off()
