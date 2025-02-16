library(MASS)
library(mgcv)

attach(mcycle)

source("motorcycle-functions.R")

gam.fit <- gam(accel ~ s(times, k = 20, bs = "bs"), optimizer = "efs")
gam.fit$sp
gam.fit$sig2

df <- 40

S <- crossprod(diff(diag(df-1), differences = 2))

# Positioning of the knots
xl <- min(times); xu <- max(times); xr <- xu - xl
xl <- xl - 0.001*xr; xu <- xu + 0.001*xr
X <- model.matrix(accel ~ 0 + splines::bs(times, df = df, Boundary.knots = c(xl,xu), intercept = TRUE))

fit <- EstimatePenal(S, lambda.init = 10)
y.fit <- X %*% fit$beta

cbind(fitted(gam.fit), y.fit)

plot(times,accel)
lines(sort(times), fitted(gam.fit)[order(times)])
lines(sort(times), y.fit[order(times)], col = "red")



