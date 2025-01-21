
data(mtcars)
head(mtcars)

attach(mtcars)
detach(mtcars)

library(mgcv)

source("cars-functions.R")


fit.gam <- gam_model <- gam(mpg ~ s(hp, k = 10, bs = "bs"), optimizer = 'efs')
c(fit.gam$coef, fit.gam$sig2)
fit.gam$sp


df <- 10
S <- crossprod(diff(diag(df-1), differences = 1))

test <- EstimatePenal(S = S, lambda.init = 10)

xl <- min(hp); xu <- max(hp); xr <- xu - xl
xl <- xl - 0.001*xr; xu <- xu + 0.001*xr
X <- model.matrix(mpg ~ 0 + splines::bs(hp, df = df, Boundary.knots = c(xl,xu), intercept = TRUE))
y.fit <- X %*% test$beta

# nk <- df - 3 + 1 # Number of knots
# xl <- min(hp); xu <- max(hp); xr <- xu - xl
# xl <- xl - 0.001*xr; xu <- xu + 0.001*xr
# k <- seq(xl, xu, length = nk)
# X <- model.matrix(mpg ~ splines::bs(hp, knots = k[-c(1,length(k))], Boundary.knots = k[c(1, length(k))]))
# y.fit <- X %*% test$beta 

cbind(fitted(fit.gam), y.fit)

plot(hp,mpg)
lines(sort(hp), fitted(fit.gam)[order(hp)])
lines(sort(hp), y.fit[order(hp)], col = "red")

