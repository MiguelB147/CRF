
library(Rcpp)
library(splines)
library(expm)


posthoc <- function(grid, MLE, hessian, S) {
  
  V <- solve(hessian)
  Vinv <- hessian
  
  Dm <- Dw <- Dw.plot <- rep(NA, nrow(grid))
  for (i in 1:nrow(grid)) {
    
    Q <- grid[i,1]*S[[1]] + grid[i,2]*S[[2]]
    A <- solve(Vinv + Q) %*% Vinv
    Vp <- A %*% V %*% t(A)
    MLEp <- A %*% MLE
    
    U <- chol(Vinv)
    N <- length(MLE)
    
    UVpU <- U %*% Vp %*% t(U)
    UdiagVpU <- U %*% diag(Vp) %*% t(U)
    
    Dm[i] <- t(MLEp - MLE) %*% Vinv %*% (MLEp - MLE)
    Dw[i] <- Dm[i] + N + sum(diag(UVpU - 2*sqrtm(UVpU)))
    Dw.plot[i] <- Dm[i] + N + sum(diag(UdiagVpU - 2*sqrtm(UdiagVpU)))
    
    # Dw[i] <- sum((MLEp - MLE)^2) + sum(diag(V + Vp - 2*sqrtm(sqrtm(V) %*% Vp %*% sqrtm(V))))
    # Dw.plot[i] <- sum((MLEp - MLE)^2) + sum(diag(V + diag(Vp) - 2*sqrtm(sqrtm(V) %*% diag(Vp) %*% sqrtm(V))))
  }
  
  return(list(Dm = Dm, Dw = Dw, Dw.plot = Dw.plot))
}


source('hunan-functions.R')
sourceCpp('test.cpp')

degree = 2
df = 6
K <- 1000
unif.ub <- NULL # 5 = 20% censoring, 2.3 = 40% censoring

datalist <- SimData(K = K, df = df, degree = degree, unif.ub = unif.ub)



fit <- nlm(f = wrapper, p = rep(1,df^2), degree = degree, datalist = datalist, hessian = TRUE, steptol = 1e-10)
fit2 <-  optim(par = rep(1,df^2), fn = wrapper, method = "BFGS", degree = degree, datalist = datalist, hessian = TRUE, control = list(reltol = 1e-10))
fit3 <- nloptr(eval_f = wrapper, x0 = rep(1,df^2), Sl = NULL, H = NULL, minusLogLik = TRUE, degree = degree, datalist = datalist, opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-10, xtol_abs = 1e-10, maxeval = 1e+05))

test <- derivatives(fit$estimate, degree, datalist, gradient = TRUE, hessian = TRUE)
test2 <- derivatives(fit2$par, degree, datalist, gradient = TRUE, hessian = TRUE)

library(nloptr)



head(fit2$hessian)
head(test2$hessian)

diag(solve(fit$hessian))
diag(solve(test$hessian))
S <- list(Srow(df), Scol(df))
grid <- expand.grid(seq(0.5,20,0.5), seq(0.5,20,0.5))

test <- posthoc(grid, MLE = fit$estimate, fit$hessian, S)
