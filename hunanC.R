
library(copula)
library(Rcpp)
library(rootSolve)


source('hunan-functions.R')
sourceCpp('test.cpp')

# Simulation ----

{

# cl <- makeCluster(detectCores()-1)
# registerDoParallel(cl)

nsim <- 200

# pb <- progress_bar$new(
#   format = "Simulation = [:bar] :percent [Elapsed time: :elapsedfull | Estimated time remaining: :eta]",
#   total = nsim,
#   clear = FALSE)
# 
# progress <- function(n){
#   pb$tick()
# } 

degree = 3
df = 8
K <- 1000
unif.ub <- NULL # 5 = 20% censoring, 2.3 = 40% censoring

plot.grid <- expand.grid(seq(0.25,1.5,by=0.25), seq = seq(0.1,2.2,by = 0.05))
names(plot.grid) <- c("time1","time2")
plot.grid$true <- theta.frank(plot.grid$time1,plot.grid$time2, alpha = 0.0023)

# results <- foreach(i=1:nsim,
#                    .packages = c('splines','MASS','rootSolve'),
#                    .errorhandling = "remove",
#                    .combine = 'cbind') %dopar% {

CRF <- matrix(NA, ncol = nsim, nrow = nrow(plot.grid))
for (i in 1:nsim) {

datalist <- SimData(K = K, df = df, degree = degree, unif.ub = unif.ub)
  
# fit <- nlm(f = loglikCpp,
#            p = rep(1,df^2),
#            degree = degree,
#            df = df,
#            datalist = datalist,
#            hessian = TRUE)

lambda <- c(2,2)
Sl <- lambda[1]*Srow(df) + lambda[2]*Scol(df)

# NOTE met penalty minder problemen met gradient=0
fit <- multiroot(Score, start = rep(1, df^2), maxiter = 100, rtol = 1e-10, degree = degree, datalist = datalist, Sl = Sl)


CRF[,i] <- mapply(function(x,y) exp(tensor(x,y, coef.vector = fit$root,
                                       degree = degree, df = df, knots = datalist$knots)),
              plot.grid$time1,
              plot.grid$time2)

}

# stopCluster(cl)


plot.grid$tensor.mean <- apply(CRF, MARGIN = 1, mean)
plot.grid$tensor.median <- apply(CRF, MARGIN = 1, median)
plot.grid$tensor.se <- apply(CRF, MARGIN = 1, sd)
plot.grid$tensor.median.ub <- plot.grid$tensor.median + plot.grid$tensor.se
plot.grid$tensor.median.lb <- plot.grid$tensor.median - plot.grid$tensor.se

}


plot.grid$poly.mean <- apply(results, MARGIN = 1, mean)
plot.grid$poly.median <- apply(results, MARGIN = 1, median)
plot.grid$poly.se <- apply(results, MARGIN = 1, sd)
plot.grid$poly.median.ub <- plot.grid$poly.median + plot.grid$poly.se
plot.grid$poly.median.lb <- plot.grid$poly.median - plot.grid$poly.se




## Pointwise tensor plot ----


# pdf(file = "/Users/miguel/Documents/CRF code/Spline fit/K400/3 internal knots - degree 2/pointwisemedian.pdf")

# pdf('/Users/miguel/Library/CloudStorage/GoogleDrive-miguel-angel.beynaerts@student.uhasselt.be/Mijn Drive/CRF simulations/Run200degree2df5K1000cens20.pdf')
# pdf('H:/My Drive/CRF simulations/Run200cens40.pdf')
par(mfrow = c(2,3))
# title(main = paste0("degree=", degree, " ", "df=", df, "K=", K, " ", nsim, " ", "runs"))
for (i in unique(plot.grid$time1)) {
  
  plottext <- paste0("t[1] == ", i)
  plot(x = plot.grid$time2[plot.grid$time1 == i],
       y = plot.grid$true[plot.grid$time1 == i],
       type = 'l',
       lwd = 2,
       ylab = "CRF", xlab = expression(t[2]), main = parse(text = plottext),
       ylim = c(0,7))

  lines(plot.grid$time2[plot.grid$time1 == i], plot.grid$tensor.median[plot.grid$time1 == i], col = 'grey', lwd = 2)
  points(x = plot.grid$time2[plot.grid$time1 == i], y = plot.grid$tensor.median.lb[plot.grid$time1 == i], pch = 20, cex = 0.5)
  points(x = plot.grid$time2[plot.grid$time1 == i], y = plot.grid$tensor.median.ub[plot.grid$time1 == i], pch = 20, cex = 0.5 )
  legend("topright", legend = c("True", "Estimated"), lty = c(1,1), col = c("black","red"))
  
}
# mtext(paste0("K = ", K, 'degree = ', degree, "df = ", df, " ", nsim, " ", "runs"), side = 3, line = -2, outer = TRUE)
par(mfrow = c(1,1))
dev.off()



## Polynomial plot ----

pdf(file = "/Users/miguel/Documents/CRF code/Polynomial fit/K1000/pointwise500.pdf")
pdf("H:/My Drive/CRF simulations/Polynomial/Run200K1000cens20.pdf")
par(mfrow = c(2,3))
for (i in unique(plot.grid$time1)) {
  
plottext <- paste0("t[1] == ", i)
plot(x = plot.grid$time2[plot.grid$time1 == i],
     y = plot.grid$true[plot.grid$time1 == i],
     type = 'l', lwd = 2,
     ylab = "CRF", xlab = expression(t[2]), main = parse(text = plottext),
     ylim = c(0,7))
lines(plot.grid$time2[plot.grid$time1 == i], plot.grid$poly.median[plot.grid$time1 == i], col = 'red', lwd = 2)
points(x = plot.grid$time2[plot.grid$time1 == i], y = plot.grid$poly.median.lb[plot.grid$time1 == i], pch = 20, cex = 0.5)
points(x = plot.grid$time2[plot.grid$time1 == i], y = plot.grid$poly.median.ub[plot.grid$time1 == i], pch = 20, cex = 0.5 )
legend("topright", legend = c("True", "Estimated"), lty = c(1,1), col = c("black","red"))

}
par(mfrow = c(1,1))
dev.off()



# Fellner - Schall method ----

library(copula)
library(Rcpp)
library(rootSolve)

source('hunan-functions.R')
sourceCpp('test.cpp')

degree = 3
df = 7
K <- 1000
unif.ub <- NULL # 5 = 20% censoring, 2.3 = 40% censoring

# 
# spline1 <- WoodSpline(datalist$X[,1], dim = 10, scale = FALSE, type = "ps")
# spline1$S
# spline2 <- WoodSpline(datalist$X[,2], dim = 10, scale = FALSE, type = "ps")
# spline2$S
# 
# D1 <- spline1$D %x% diag(rep(1, 10))
# S1 <- crossprod(D1)
# 
# D2 <- diag(rep(1, 10)) %x% spline2$D
# S2 <- crossprod(D2)


set.seed(123)
datalist <- SimData(K = K, df = df, degree = degree, unif.ub = unif.ub, alpha = 20)

# S <- list(Srow(df), Scol(df))
# lambda <- c(50,50)
# Sl<- lambda[1]*S[[1]] + lambda[2]*S[[2]]
# S <- Srow(df)
# 
# test <- eigen(Sl)
# prod(test$values[test$values > 0])
# log(prod(test$values[test$values > 0]))

# S <- list(Srow(df, diff = degree-1), Scol(df, diff = degree-1))
# fit <- EstimatePenal(datalist = datalist, degree = degree, S = S, lambda.init = c(10,10), step.control = F)

beta.start <- fit7bsscale$beta
fit7bsquant <- EstimatePenal2(datalist = datalist, dim = 7, type = "ps", lambda.init = c(10,10), scale = FALSE, repara = TRUE, quantile = FALSE)


test <- efsud.fit(fit$beta, degree = degree, datalist = datalist, Sl = 2*S[[1]] + 0*S[[2]])

plot.grid <- expand.grid(seq(0.25,1.5,by=0.25), seq = seq(0.1,2.2,by = 0.05))
names(plot.grid) <- c("time1","time2")

plot.grid$true <- theta.frank(plot.grid$time1,plot.grid$time2)
plot.grid$psscale <- WoodTensor.predict(plot.grid$time1, plot.grid$time2, fit7psscale)
plot.grid$bsscale <- WoodTensor.predict(plot.grid$time1, plot.grid$time2, fit7bsscale)
plot.grid$ps <- WoodTensor.predict(plot.grid$time1, plot.grid$time2, fit7ps)
plot.grid$bs <- WoodTensor.predict(plot.grid$time1, plot.grid$time2, fit7bs)
plot.grid$bsquant <- WoodTensor.predict(plot.grid$time1, plot.grid$time2, fit7bsquant)
# plot.grid$CRFbsscalequant <- WoodTensor.predict(plot.grid$time1, plot.grid$time2, fit7bsscalequant)



CRF <- WoodTensor.predict(plot.grid$time1, plot.grid$time2, fit)
# CRF <- mapply(function(x,y) WoodTensor.predict(x, y, fit),
#               plot.grid$time1,
#               plot.grid$time2)

# CRF <- mapply(function(x,y) exp(tensor(x,y, coef.vector = fit$beta,
#                                        degree = degree, df = df, knots = datalist$knots)),
#               plot.grid$time1,
#               plot.grid$time2)
# CRF <- mapply(function(x,y) exp(tensor(x,y, coef.vector = test$beta,
#                                        degree = degree, df = df, knots = datalist$knots)),
#               plot.grid$time1,
#               plot.grid$time2)




par(mfrow = c(2,3))
for (i in unique(plot.grid$time1)) {
  
  plottext <- paste0("t[1] == ", i)
  plot(x = plot.grid$time2[plot.grid$time1 == i],
       y = plot.grid$true[plot.grid$time1 == i],
       type = 'l', lwd = 2, col = "grey",
       ylab = "CRF", xlab = expression(t[2]), main = parse(text = plottext),
       ylim=c(0,7))
  lines(x = plot.grid$time2[plot.grid$time1 == i], y = plot.grid$psscale[plot.grid$time1 == i], col = 1, lwd = 2)
  lines(x = plot.grid$time2[plot.grid$time1 == i], y = plot.grid$bsscale[plot.grid$time1 == i], col = 2, lwd = 2)
  lines(x = plot.grid$time2[plot.grid$time1 == i], y = plot.grid$ps[plot.grid$time1 == i], col = 3, lwd = 2)
  lines(x = plot.grid$time2[plot.grid$time1 == i], y = plot.grid$bs[plot.grid$time1 == i], col = 4, lwd = 2)
  # lines(x = plot.grid$time2[plot.grid$time1 == i], y = plot.grid$bsscalequant[plot.grid$time1 == i], col = 5, lwd = 2)
  lines(x = plot.grid$time2[plot.grid$time1 == i], y = plot.grid$bsquant[plot.grid$time1 == i], col = 6, lwd = 2)
}
par(mfrow = c(1,1))





MatToVec <- function(Matrix) {
  
  rows <- dim(Matrix)[1]
  columns <- dim(Matrix)[2]
  
  vector <- vector(length = rows*columns)
  k <- 1
  for (j in 1:columns) {
    for (i in 1:rows) {
      vector[k] <- Matrix[i,j]
      k <- k + 1
    }
  }
  
  return(vector)
}


