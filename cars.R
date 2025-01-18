
data(mtcars)
head(mtcars)

attach(mtcars)

library(mgcv)
library(splines)

source("cars-functions.R")


fit.gam <- gam_model <- gam(mpg ~ s(hp, k = 10), optimizer = 'efs')
c(fit.gam$coef, fit.gam$sig2)
fit.gam$sp

S <- crossprod(diff(diag(10)))
Sl <- fit.gam$sp*S


EstimatePenal(S = S, lambda.init = 10)

test <- nlm(loglikpenal, c(fit.gam$coef, fit.gam$sig2), Sl = Sl)

plot(hp,mpg)
lines(sort(hp), fitted(fit.gam)[order(hp)])


S <- crossprod(diff(diag(10)))
lambda <- seq(1,100, by = 0.5)

mLogLik <- avg <- c()
for (i in 1:length(lambda)) {
  
  S.lambda <- lambda[i]*S
  S.lambda.inv <- ginv(S.lambda)
  for(j in 1:1000) {
    beta <- mvrnorm(1, mu = rep(0,10), Sigma = S.lambda.inv)
    mLogLik[j] <- loglik(c(beta,3)) + t(beta) %*% S.lambda %*% beta
  }
  avg[i] <- mean(mLogLik)
}

lambda[which.min(avg)]

test <- nlm(loglik, p = , hessian = FALSE)

test <- EstimatePenal(S, lambda.init = 20)

coef(fit.gam)

