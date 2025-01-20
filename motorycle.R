library(MASS)
library(mgcv)

attach(mcycle)

source("motorcycle-functions.R")

gam.fit <- gam(accel ~ s(times, k = 30, bs = "ps", sp = 0))
gam.fit$sp
gam.fit$sig2

df <- 30

S <- crossprod(diff(diag(df), differences = 2))

# Positioning of the knots
nk <- df - 3 + 1# Number of knots
xl <- min(times); xu <- max(times); xr <- xu - xl
xl <- xl - 0.001*xr; xu <- xu + 0.001*xr
k <- seq(xl, xu, length = nk)

X <- model.matrix(accel ~ splines::bs(times, knots = k[-c(1,length(k))], Boundary.knots = k[c(1, length(k))]))

fit <- EstimatePenal(S)
y.fit <- X %*% fit$beta

cbind(fitted(gam.fit), y.fit)

plot(times,accel)
lines(sort(times), fitted(gam.fit)[order(times)])
lines(sort(times), y.fit[order(times)], col = "red")



