
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

#######################################

# hunan-functions.R

source('hunan-functions.R')
sourceCpp('test.cpp')

set.seed(123)
datalist <- SimData(K = 1000, w = c(0.2, 0.5, 0.3), alpha = c(3, 3, 2))

obj1 <- WoodSpline(t = datalist$X[,1], dim = 7, degree = 3)
obj2 <- WoodSpline(t = datalist$X[,2], dim = 7, degree = 3)

X1 <- obj1$X
X2 <- obj2$X

beta <- multiroot(Score2, start = rep(1, 7^2), rtol = 1e-10, X1 = X1, X2 = X2, Sl = NULL, datalist = datalist)$root
H <- derivatives2(coef.vector = beta, X1 = X1, X2 = X2, datalist = datalist, gradient = FALSE, hessian = TRUE)$hessian

# hunan-functions2.R

source('hunan-functions2.R')
sourceCpp('test2.cpp')

set.seed(123)
datalist <- SimData(K = 1000, w = c(0.2, 0.5, 0.3), alpha = c(3, 3, 2))

obj1 <- WoodSpline(t = datalist$X[,1], dim = 7, degree = 3)
obj2 <- WoodSpline(t = datalist$X[,2], dim = 7, degree = 3)

X1 <- obj1$X
X2 <- obj2$X

Pen <- WoodPenalty(obj1, obj2)

Sl <- Pen[[1]] + Pen[[2]]

deriv <- deriv_comp(X1 = X1, X2 = X2, datalist = datalist)
beta <- multiroot(Score2, start = rep(1, 7^2),
                  jacfunc = Hessian, jactype = "fullusr", rtol = 1e-10,
                  X1 = X1, X2 = X2, deriv = deriv, datalist = datalist, Sl = Sl)$root
test <- nleqslv(x = rep(1, 7^2), fn = Score2, jac = Hessian, method = "Broyden", global = "hook",
                X1 = X1, X2 = X2, deriv = deriv, datalist = datalist, Sl = Sl)
beta2 <- test$x
H1 <- Hessian(coef.vector = beta, X1 = X1, X2 = X2, deriv = deriv, datalist = datalist, Sl = Sl)
H2 <- Hessian(coef.vector = beta2, X1 = X1, X2 = X2, deriv = deriv, datalist = datalist, Sl = Sl)

H1[1:3,1:3]
H2[1:3,1:3]

#######################################





library(CRFpackage)


set.seed(123)
datalist <- SimData(K = 1000, w = c(0.2, 0.4, 0.4), alpha = c(2, 3, 1.25), margin = "exp")

plot.grid <- expand.grid(seq(0.25,5,by=0.25), seq = seq(0.1,5,by = 0.05))
names(plot.grid) <- c("time1","time2")
plot.grid$true <- theta.mix(plot.grid$time1, plot.grid$time2,  w = c(0.2, 0.4, 0.4), alpha = c(2, 3, 1.25), margin = "exp")

# plot.grid$true <- theta.frank(plot.grid$time1, plot.grid$time2)
# fit.poly <- multiroot(f = poly.fit, start = rep(0,10), datalist = datalist); poly <- exp(polynomial(plot.grid$time1, plot.grid$time2, fit.poly$root))

fit2 <- EstimatePenal2(datalist = datalist, dim = 12, lambda.init = c(5,5), quantile = TRUE, type = "gps", repara = FALSE)
fit.poly <- EstimatePoly(datalist = datalist)

poly <- exp(polynomial(plot.grid$time1, plot.grid$time2, fit.poly))
CRF <- WoodTensor.predict(plot.grid$time1, plot.grid$time2, fit)

pdf("mixtureunifgps.pdf")
par(mfrow = c(2,3))
for (i in unique(plot.grid$time1)) {
  plottext <- paste0("t[1] == ", i)
  plot(x = plot.grid$time2[plot.grid$time1 == i],
       y = plot.grid$true[plot.grid$time1 == i],
       type = 'l', lwd = 2, col = "grey",
       xlim = c(0, max(plot.grid$time2)),
       ylim = c(0,7),
       ylab = "CRF", xlab = expression(t[2]), main = parse(text = plottext))
  
  lines(x = plot.grid$time2[plot.grid$time1 == i], y = CRF[plot.grid$time1 == i], col = "red", lwd = 2)
  for (j in 1:length(fit$knots[["knots2"]])) {
    abline(v = fit$knots[["knots2"]][j], col = "lightgrey")
  }
  
  # lines(x = plot.grid$time2[plot.grid$time1 == i], y = poly[plot.grid$time1 == i], col = "blue", lwd = 2)
}
par(mfrow = c(1,1))
dev.off()

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


