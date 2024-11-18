
library(splines2)
library(doParallel)

theta.frank <- function(x,y,alpha) {
  A <- (alpha-1)*log(alpha)*alpha^(2-exp(-x)-exp(-y))
  B <- (alpha^(1-exp(-x))-alpha)*(alpha^(1-exp(-y))-alpha)
  C <- -1 + exp(-x) + exp(-y) + log(1+ (alpha^(1-exp(-x))-1)*(alpha^(1-exp(-y))-1)/(alpha-1), base = alpha)
  return(A*C/B)
}

tensor <- function(t1, t2, coef.matrix, df, degree, knots) {
  
  if (is.null(knots)) {
    
    spline1 <- bs(t1, df = df, degree = degree)
    B1 <- matrix(as.numeric(spline1), nr = nrow(spline1))
    spline2 <- bs(t2, df = df, degree = degree)
    B2 <- matrix(as.numeric(spline2), nr = nrow(spline2))
    spline12 <- B1 %*% coef.matrix %*% t(B2)
    
  } else {
    
    boundary1 <- knots[c(1,df),1]
    boundary2 <- knots[c(1,df),2]
    int1 <- knots[2:(df-1),1]
    int2 <- knots[2:(df-1),1]
    
    spline1 <- bs(t1, knots = int1, degree = degree, Boundary.knots = boundary1)
    B1 <- matrix(as.numeric(spline1), nr = nrow(spline1))
    spline2 <- bs(t2, knots = int2, degree = degree, Boundary.knots = boundary2)
    B2 <- matrix(as.numeric(spline2), nr = nrow(spline2))
    spline12 <- B1 %*% coef.matrix %*% t(B2)
    
  }
  
  return(spline12)
}

# tensor2 <- function(t1, t2, coef.matrix, int.knots, degree, support) {
#   
#   
#   spline1 <- bSpline(t1, knots = int.knots[,1], degree = degree, Boundary.knots = support[,1])
#   spline2 <- bSpline(t2, knots = int.knots[,2], degree = degree, Boundary.knots = support[,2])
#   spline12 <- spline1 %*% coef.matrix %*% t(spline2)
#   return(spline12)
#   
# }
# 
# tensor3 <- function(t1, t2, coef.matrix, degree, df, knots) {
#   
#   if (is.null(knots)) {
#     spline1 <- bSpline(t1, df = df, degree = degree)
#     spline2 <- bSpline(t2, df = df, degree = degree)
#     spline12 <- spline1 %*% coef.matrix %*% t(spline2)
#   } else {
#     boundary1 <- knots[c(1,df),1]
#     boundary2 <- knots[c(1,df),2]
#     int1 <- knots[2:(df-1),1]
#     int2 <- knots[2:(df-1),1]
#     
#     spline1 <- bSpline(t1, df = df, degree = degree, Boundary.knots = boundary1, knots = int1)
#     spline2 <- bSpline(t2, df = df, degree = degree, Boundary.knots = boundary2, knots = int2)
#     spline12 <- spline1 %*% coef.matrix %*% t(spline2)
#   }
#   
#   return(spline12)
# } 

# bSpline(1, degree = 3, Boundary.knots = c(0,2.3), knots = 1.15)
# bs(1, knots = 1.15, degree = 3, Boundary.knots = c(0,2.3))


polynomial <- function(t1,t2,coef.vec) {
  
  logtheta <- coef.vec[1] + coef.vec[2]*t1 + coef.vec[3]*t2 +
    coef.vec[4]*t1^2 + coef.vec[5]*t2^2 + coef.vec[6]*t1*t2 + 
    coef.vec[7]*(t1^2)*t2 + coef.vec[8]*t1*(t2^2) +
    coef.vec[9]*t1^3 + coef.vec[10]*t2^3
  
  return(logtheta)
    
}

# test1 <- mapply(tensor2, X[,1],X[,2],
#                 MoreArgs = list(int.knots = matrix(c(1.15,1.15), ncol = 2), degree = 3, support = data.frame(c1 = c(0,2.3), c2 = c(0,2.3)), coef.matrix = matrix(rep(1,16), ncol = 4)) )
# 
# test2 <- tensor2(X[,1], X[,2], int.knots = matrix(c(1.15,1.15), ncol = 2), degree = 3, support = data.frame(c1 = c(0,2.3), c2 = c(0,2.3)), coef.matrix = matrix(rep(1,16), ncol = 4))
# test3 <- outer(X[,1],X[,2], function (x,y) mapply(tensor2, x,y,
#                                    MoreArgs = list(int.knots = matrix(c(1.15,1.15), ncol = 2), degree = 3, support = data.frame(c1 = c(0,2.3), c2 = c(0,2.3)), coef.matrix = matrix(rep(1,16), ncol = 4)) ))

### N(t1 = x, t2 = y)
riskset <- function(x, y) {
  N <- sum(1*(X[,1] >= x & X[,2] >= y))
  return(N)
}

#### coef.matrix: (degree+int.knots) x (degree+int.knots)
#### penal: vector c(lambda1,lambda2)
# loglik <- function(coef.vector, degree, int.knots, support) {
#   
#   coef.matrix <- matrix(coef.vector, ncol = degree + nrow(int.knots), byrow = FALSE)
#   
#   # theta <- t(tensor(X, degree=degree, coef.matrix = coef.matrix, int.knots = int.knots))
#   logtheta2 <- tensor(X[,1], X[,2], degree = degree, coef.matrix = coef.matrix, int.knots = int.knots, support = support)
#   logtheta1 <- t(logtheta2)
#   
#   B1 <- c(I5*logtheta1)[N1 > 0]
#   B2 <- c(I6*logtheta2)[N2 > 0]
#   
#   C1 <- c(N1 + I2*(exp(logtheta1)-1))[N1 > 0]
#   C2 <- c(N2 + I4*(exp(logtheta2)-1))[N2 > 0]
#   
#   
#   L1 <- sum(A1*(B1 - log(C1)))
#   L2 <- sum(A2*(B2 - log(C2)))
#   
#   return(-(L1+L2))
# }

loglikdf <- function(coef.vector, degree, df, knots) {
  
  coef.matrix <- matrix(coef.vector, ncol = df, byrow = FALSE)
  
  logtheta2 <- tensor(X[,1], X[,2], degree = degree, coef.matrix = coef.matrix, df = df, knots = knots)
  logtheta1 <- t(logtheta2)
  
  B1 <- c(I5*logtheta1)[N1 > 0]
  B2 <- c(I6*logtheta2)[N2 > 0]
  
  C1 <- c(N1 + I2*(exp(logtheta1)-1))[N1 > 0]
  C2 <- c(N2 + I4*(exp(logtheta2)-1))[N2 > 0]
  
  
  L1 <- sum(A1*(B1 - log(C1)))
  L2 <- sum(A2*(B2 - log(C2)))
  
  return(-(L1+L2))
}

loglik.poly <- function(coef.vector) {
  
  logtheta2 <- outer(X[,1], X[,2], function (x,y) polynomial(x,y, coef.vec = coef.vector))
  logtheta1 <- t(logtheta2)
  
  B1 <- c(I5*logtheta1)[N1 > 0]
  B2 <- c(I6*logtheta2)[N2 > 0]
  
  C1 <- c(N1 + I2*(exp(logtheta1)-1))[N1 > 0]
  C2 <- c(N2 + I4*(exp(logtheta2)-1))[N2 > 0]
  
  
  L1 <- sum(A1*(B1 - log(C1)))
  L2 <- sum(A2*(B2 - log(C2)))
  
  return(-(L1+L2))
}

#########################################
## Simulating correlated survival data ##
#########################################

library(splines)
library(doParallel)
library(Rcpp)
# library(splines2)

# setwd("~/Library/CloudStorage/GoogleDrive-miguel-angel.beynaerts@student.uhasselt.be/Mijn Drive/CRF simulations")
setwd("H:/My Drive/CRF simulations")
source('hunan-functions.R')
sourceCpp('test.cpp')

cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

nsim <- 200

degree = 2
df = 4
K <- 1000
unif.ub <- 2.3 # 5 = 20% censoring, 2.3 = 40% censoring

plot.grid <- expand.grid(seq(0.25,1.5,by=0.25), seq = seq(0.1,2.2,by = 0.05))
names(plot.grid) <- c("time1","time2")
plot.grid$true <- theta.frank(plot.grid$time1,plot.grid$time2, alpha = 0.0023)

results <- foreach(icount(nsim), .options.multicore=list(set.seed=FALSE), .packages = c('splines'), .errorhandling = "remove", .combine = 'cbind', .verbose = TRUE) %dopar% {


# set.seed(123)
u1 <- runif(K, 0, 1)
u2 <- runif(K, 0, 1)

alpha <- 0.0023
a <- alpha^u1 + (alpha - alpha^u1)*u2

# Fan 2000
T1 <- -log(u1)
T2 <- -log(log(a/(a+(1-alpha)*u2),base = alpha))

# Hu 2011

C1 <- runif(K, 0, unif.ub) # 40% censoring
C2 <- runif(K, 0, unif.ub) # 40% censoring

X1 <- pmin(T1,C1)
X2 <- pmin(T2,C2)

knots1 <- seq(min(X1)-1, max(X1)+1, length.out = df)
knots2 <- seq(min(X2)-1, max(X2)+1, length.out = df)

X <- as.matrix(cbind(X1,X2))

delta1 <- 1*(T1 <= C1)
delta2 <- 1*(T2 <= C2)

delta <- as.matrix(cbind(delta1,delta2))


###############################
## Check no ties assumptions ##
###############################

# length(unique(X[,1]))
# length(unique(X[,2]))

#########################################
## Check whether first delta1=delta2=1 ##
#########################################


row_index <- which(delta1 == 1 & delta2 == 1 , arr.ind = TRUE)[1] # First observation with delta1=delta2=1

# Switch rows
if (row_index > 1) {
  
  tmp_row_X <- X[1,]
  tmp_row_delta <- delta[1,]
  
  X[1,] <- X[row_index,]
  delta[1,] <- delta[row_index,]
  
  X[row_index,] <- tmp_row_X
  delta[row_index,] <- tmp_row_delta
  
  rm(tmp_row_delta, tmp_row_X, row_index)
  
} else {X <- X; delta <- delta}


##############################
## Calculating the risk set ##
##############################


# {
#   risk_tmp1 <- risk_tmp2 <- c()
#   N1 <- N2 <- matrix(nrow = length(X1), ncol = length(X1))
#   for (i in 1:length(X1)) {
#     for (j in 1:length(X1)) {
#       risk_tmp1[j] <- riskset(X[i,1],X[j,2])
#       risk_tmp2[j] <- riskset(X[j,1],X[i,2])
#     }
#     N1[,i] <- risk_tmp1
#     N2[,i] <- risk_tmp2
#   }
#   rm(risk_tmp1, risk_tmp2, i, j)
# }

N <- outer(X[,1], X[,2], function(x,y) mapply(riskset,x,y))
N <- risksetC(X[,1], X[,2])
N1 <- c(t(N))
N2 <- c(N)

# {
#   riskset1 <- riskset_tmp <- c()
#   for (i in 1:nrow(X) ){
#     for (j in 1:nrow(X)) {
#       riskset_tmp[j] <- sum(1*(X[,1] >= X[i,1]) & X[,2] >= X[j,2])
#     }
#     riskset1[i] <- sum(riskset_tmp)
#   }
#   rm(riskset_tmp)
# }
# 
# {
#   riskset2 <- riskset_tmp <- c()
#   for (i in 1:nrow(X) ){
#     for (j in 1:nrow(X)) {
#       riskset_tmp[j] <- sum(1*(X[,1] >= X[j,1] & X[,2] >= X[i,2]))
#     }
#     riskset2[i] <- sum(riskset_tmp)
#   }
#   rm(riskset_tmp)
# }

###################################################
## Calculating indicator functions in likelihood ##
###################################################

#### I(X1j >= X1i)
I1 <- sapply(X[,1], function(x) 1*(X[,1] >= x)) # col=1,...,i,...,n row=1,...,j,...,n

#### I(X2j <= X2i)
I2 <- sapply(X[,2], function(x) 1*(X[,2] <= x)) # col=1,...,i,...,n row=1,...,j,...,n
#### I(X2j >= X2i)
# I3 <- sapply(X[,2], function(x) 1*(X[,2] >= x)) # col=1,...,i,...,n row=1,...,j,...,n
I3 <- t(I2)

#### I(X1j <= X1i)
# I4 <- sapply(X[,1], function(x) 1*(X[,1] <= x)) # col=1,...,i,...,n row=1,...,j,...,n
I4 <- t(I1)

#### I(X1j = X1i) NOTE THAT THIS IS DIAG(1,500,500) IF NO TIES
I5 <- sapply(X[,1], function(x) 1*(X[,1] == x)) # col=1,...,i,...,n row=1,...,j,...,n

#### I(X2j = X2i) NOTE THAT THIS IS DIAG(1,500,500) IF NO TIES
I6 <- sapply(X[,2], function(x) 1*(X[,2] == x)) # col=1,...,i,...,n row=1,...,j,...,n

#I1 <- lapply(X1, function(x) 1*(X2 >= x))
#test <- matrix(unlist(I1), ncol = 500, byrow = FALSE)


A1 <- c(I1*outer(delta[,2], delta[,1]))[N1 > 0]
A2 <- c(I3*outer(delta[,1], delta[,2]))[N2 > 0]

delta.prod <- DeltaC(delta[,1],delta[,2])



##################################
## Calculating total likelihood ##
##################################


# fit <- nlm(f = loglikdf,
#            p = rep(1,df^2),
#            degree = degree,
#            df = df,
#            knots = cbind(knots1,knots2),
#            iterlim = 10000,
#            hessian = FALSE)


fit <- nlm(f = loglikC,
           p = rep(1,df^2),
           degree = degree,
           df = df,
           knots = cbind(knots1,knots2),
           iterlim = 10000,
           hessian = FALSE)


A.hat <- matrix(fit$estimate, ncol = df, byrow = FALSE)
CRF <- mapply(function(x,y) exp(tensor(x,y, coef.matrix = A.hat,
                                       degree = degree, df = df, knots = cbind(knots1,knots2))),
              plot.grid$time1,
              plot.grid$time2)

# fit <- nlm(f = loglik.poly,
#            p = rep(1,10),
#            iterlim = 10000,
#            hessian = FALSE)
# 
# CRF <- exp(polynomial(plot.grid$time1, plot.grid$time2, fit$estimate))

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



#####################
## Score functions ##
#####################


beta.deriv13 <- polynomial(X[,1],X[,2],coef.vector = rep(1,10))

U1 <- U3 <- (1/n)*sum(delta[,2]*delta[,1]*beta.deriv13)



