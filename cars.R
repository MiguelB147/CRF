
data(mtcars)
head(mtcars)

library(mgcv)
library(splines)

source("cars-functions.R")


fit.gam <- gam(mpg ~ s(hp, k = 20, bs = "ps"), data = mtcars, optimizer = 'efs')
c(fit.gam$coef, fit.gam$sig2)
fit.gam$sp
fit.gam$sig2

test <- WoodSpline(mtcars$hp, dim = 10, degree = 3)


df <- 10
S <- crossprod(diff(diag(df-1), differences = 2))

test <- EstimatePenal(dim = 30, lambda.init = 1, step.control = F, type = "bs")
test2 <- EstimatePenal(dim = 30, lambda.init = 1, step.control = F, type = "ps")
y.bs <- test$X %*% test$beta
y.ps <- test2$X %*% test2$beta

# nk <- df - 3 + 1 # Number of knots
# xl <- min(hp); xu <- max(hp); xr <- xu - xl
# xl <- xl - 0.001*xr; xu <- xu + 0.001*xr
# k <- seq(xl, xu, length = nk)
# X <- model.matrix(mpg ~ splines::bs(hp, knots = k[-c(1,length(k))], Boundary.knots = k[c(1, length(k))]))
# y.fit <- X %*% test$beta 

cbind(fitted(fit.gam), y.bs, y.ps)

plot(mtcars$hp,mtcars$mpg)
lines(sort(mtcars$hp), fitted(fit.gam)[order(mtcars$hp)])
lines(sort(mtcars$hp), y.ps[order(mtcars$hp)], col = "blue")
lines(sort(mtcars$hp), y.bs[order(mtcars$hp)], col = "red")

