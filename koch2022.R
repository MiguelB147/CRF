
library(Rcpp)
library(splines)


posthoc <- function(grid, MLE, hessian, S) {
  
  V <- solve(hessian)
  Vinv <- hessian
  
  Dm <- Dw <- Dw.plot <- rep(NA, nrow(grid))
  for (i in 1:nrow(grid)) {
    
    Q <- grid[i,1]*S[[1]] + grid[i,2]*S[[2]]
    A <- solve(Vinv + Q) %*% Vinv
    Vp <- A %*% V %*% t(A)
    MLEp <- A %*% MLE
    
    Dm[i] <- t(MLEp - MLE) %*% Vinv %*% (MLEp - MLE)
    Dw[i] <- sum((MLEp - MLE)^2) + sum(diag(V + Vp - 2*sqrtm(sqrtm(V) %*% Vp %*% sqrtm(V))))
    Dw.plot[i] <- sum((MLEp - MLE)^2) + sum(diag(V + diag(Vp) - 2*sqrtm(sqrtm(V) %*% diag(Vp) %*% sqrtm(V))))
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
S <- list(Srow(df), Scol(df))
grid <- expand.grid(seq(0.5,20,0.5), seq(0.5,20,0.5))

test <- posthoc(grid, MLE = fit$estimate, fit$hessian, S)
