
setwd("~/Library/CloudStorage/GoogleDrive-miguel-angel.beynaerts@student.uhasselt.be/Mijn Drive/CRF simulations")
# setwd("H:/My Drive/CRF simulations")
source('hunan-functions.R')
sourceCpp('test.cpp')

cl <- makeCluster(detectCores()-1)
registerDoSNOW(cl)

nsim <- 200

pb <- progress_bar$new(
  format = "Simulation = [:bar] :percent [Elapsed time: :elapsedfull | Estimated time remaining: :eta]",
  total = nsim,
  clear = FALSE)

progress <- function(n){
  pb$tick()
} 

degree = 2
df = 6
K <- 1000
unif.ub <- 5 # 5 = 20% censoring, 2.3 = 40% censoring

plot.grid <- expand.grid(seq(0.25,1.5,by=0.25), seq = seq(0.1,2.2,by = 0.05))
names(plot.grid) <- c("time1","time2")
plot.grid$true <- theta.frank(plot.grid$time1,plot.grid$time2, alpha = 0.0023)

results <- foreach(i=1:nsim,
                   .packages = c('splines','MASS'),
                   .errorhandling = "remove",
                   .combine = 'cbind',
                   .options.snow = list(progress = progress)) %dopar% {

datalist <- SimData(K = K, df = df, degree = degree, unif.ub = unif.ub)
  
fit <- nlm(f = loglikCpp,
           p = rep(1,df^2),
           degree = degree,
           df = df,
           datalist = datalist,
           hessian = TRUE)

lambda <- c(50,50)
S.lambda <- lambda[1]*Srow(df) + lambda[2]*Scol(df)

# NOTE met penalty minder problemen met gradient=0
fit <- nlm(f = wrapper,
           p = rep(1,df^2),
           degree = degree,
           S.lambda = S.lambda,
           datalist = datalist,
           hessian = TRUE)

# test <- nloptr(x0 = rep(1,df^2),
#                eval_f = loglikCpp,
#                degree = degree,
#                df = df,
#                datalist = dat.list,
#                opts = list("algorithm" = "NLOPT_LN_BOBYQA",
#                            "xtol_rel" = 1.0e-8))

# test2 <- optim(par = rep(1,df^2),
#                fn = loglikCpp,
#                degree = degree,
#                df = df,
#                datalist = dat.list,
#                method = "BFGS",
#                control = list(maxit = 10000))


S <- list(Srow(df), Scol(df))

fit <- EstimatePenalty(datalist = datalist, degree = degree, S = S, lambda.init = c(15,15))
fit <- EstimatePenaltyNoControl(datalist = datalist, degree = degree, S = S, lambda.init = c(15,15))

A.hat <- matrix(fit$estimate, ncol = df, byrow = FALSE)
CRF <- mapply(function(x,y) exp(tensor(x,y, coef.matrix = A.hat,
                                       degree = degree, df = df, knots = cbind(knots1,knots2))),
              plot.grid$time1,
              plot.grid$time2)

}

stopCluster(cl)


plot.grid$tensor.mean <- apply(results, MARGIN = 1, mean)
plot.grid$tensor.median <- apply(results, MARGIN = 1, median)
plot.grid$tensor.se <- apply(results, MARGIN = 1, sd)
plot.grid$tensor.median.ub <- plot.grid$tensor.median + plot.grid$tensor.se
plot.grid$tensor.median.lb <- plot.grid$tensor.median - plot.grid$tensor.se


plot.grid$poly.mean <- apply(results, MARGIN = 1, mean)
plot.grid$poly.median <- apply(results, MARGIN = 1, median)
plot.grid$poly.se <- apply(results, MARGIN = 1, sd)
plot.grid$poly.median.ub <- plot.grid$poly.median + plot.grid$poly.se
plot.grid$poly.median.lb <- plot.grid$poly.median - plot.grid$poly.se



###########################
## Pointwise tensor plot ##
###########################

# pdf(file = "/Users/miguel/Documents/CRF code/Spline fit/K400/3 internal knots - degree 2/pointwisemedian.pdf")

# pdf('/Users/miguel/Library/CloudStorage/GoogleDrive-miguel-angel.beynaerts@student.uhasselt.be/Mijn Drive/CRF simulations/Run200degree2df5K1000cens20.pdf')
pdf('H:/My Drive/CRF simulations/Run200cens40.pdf')
par(mfrow = c(2,3))
# title(main = paste0("degree=", degree, " ", "df=", df, "K=", K, " ", nsim, " ", "runs"))
for (i in unique(plot.grid$time1)) {
  
  plottext <- paste0("t[1] == ", i)
  plot(x = plot.grid$time2[plot.grid$time1 == i],
       y = plot.grid$true[plot.grid$time1 == i],
       type = 'l',
       ylab = "CRF", xlab = expression(t[2]), main = parse(text = plottext),
       ylim = c(0,7))

  lines(plot.grid$time2[plot.grid$time1 == i], plot.grid$tensor.median[plot.grid$time1 == i], col = 'red')
  points(x = plot.grid$time2[plot.grid$time1 == i], y = plot.grid$tensor.median.lb[plot.grid$time1 == i], pch = 20, cex = 0.5)
  points(x = plot.grid$time2[plot.grid$time1 == i], y = plot.grid$tensor.median.ub[plot.grid$time1 == i], pch = 20, cex = 0.5 )
  legend("topright", legend = c("True", "Estimated"), lty = c(1,1), col = c("black","red"))
  
}
# mtext(paste0("K = ", K, 'degree = ', degree, "df = ", df, " ", nsim, " ", "runs"), side = 3, line = -2, outer = TRUE)
par(mfrow = c(1,1))
dev.off()


#####################
## Polynomial plot ##
#####################

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

test <- matrix(1:9, ncol = 3, byrow = F)


