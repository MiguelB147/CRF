library(splines)
library(Rcpp)
library(progress)

source('hunan-functions.R')
sourceCpp('test.cpp')

degree = 2
df = 8
K <- 1000
unif.ub <- NULL # NULL = no censoring, 5 = 20% censoring, 2.3 = 40% censoring

datalist <- SimData(K = K, df = df, degree = degree, unif.ub = unif.ub)

penaltyfunc <- function(coef.vector, degree, df, Sl, datalist) {
  penalty <- t(coef.vector) %*% Sl %*% coef.vector
  ll <- loglikCpp(coef.vector, degree, df, datalist)
  return(ll + penalty/2)
}

S <- list(Srow(df), Scol(df))
Sl <- 5*(S[[1]] + S[[2]])

fit <- nlm(f = penaltyfunc,
           p = rep(1,df^2),
           degree = degree,
           df = df,
           Sl = Sl,
           datalist = datalist,
           hessian = TRUE)
fit$estimate
fit$hessian

sum(apply(fit$hessian, 2, function(x) all(x == 0)))

try <- seq(-4,4, 0.01)

pb <- progress_bar$new(
  format = "Simulation: [:bar] :percent [Elapsed time: :elapsedfull | Estimated time remaining: :eta]",
  total = length(try)*length(fit$estimate),
  clear = TRUE)

ll <- matrix(nrow = length(try), ncol = length(fit$estimate))
for (j in 1:length(fit$estimate)) {
  for (i in 1:length(try)) {
    pb$tick()
    coef <- fit$estimate
    coef[j] <- try[i]
    ll[i,j] <- penaltyfunc(coef.vector = coef, degree = degree, df = df, Sl = Sl, datalist = datalist)
    pb$tick()
  }
}

# cl <- makeCluster(detectCores()-1)
# registerDoParallel(cl)
# 
# ll <- foreach(j = 1:length(fit$estimate), .combine = 'cbind', .packages = c("splines", "Rcpp2doParallel")) %:%
#   foreach(i = 1:length(try), .combine = 'rbind') %dopar% {
#     
#     coef <- fit$estimate
#     coef[j] <- try[i]
#     loglikCpp(coef.vector = coef, degree = degree, df = df, datalist = datalist)
#     
#   }
# 
# stopCluster(cl)

# pdf(paste0("degree",degree,"df",df,".pdf"), paper = "a4")
pdf("testpenalty.pdf", paper = "a4")
# pdf(paste0("degree",degree,"df",df,".pdf"), paper = "a4")
for (i in 1:ncol(ll)) {
  plot(try, ll[,i], type = "l", lwd = 2, ylab = "-log-likelihood", xlab = paste0("beta",i))
  abline(v = fit$estimate[i], col = "red")
}
dev.off()

betas <- expand.grid(try,try)
ll.surf <- rep(NA, nrow(betas))
coef <- fit$estimate
for (i in 1:nrow(betas)) {
  coef[1] <- betas[1,]
  coef[2] <- betas[2,]
  ll.surf[i] <- loglikCpp(coef.vector = coef, degree = degree, df = df, datalist = datalist)
}

library(plotly)
plot_ly(x = betas[,1], y = betas[,2], z = ll.surf)

