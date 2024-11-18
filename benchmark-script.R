library(Rcpp)

setwd("~/Library/CloudStorage/GoogleDrive-miguel-angel.beynaerts@student.uhasselt.be/Mijn Drive/CRF simulations")
source('hunan-functions.R')
sourceCpp("test.cpp")


K <- 1000
unif.ub <- 2.3

set.seed(1234)
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

# knots1 <- seq(min(X1)-1, max(X1)+1, length.out = df)
# knots2 <- seq(min(X2)-1, max(X2)+1, length.out = df)

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






result <- microbenchmark(
  outer(X[,1], X[,2], function(x,y) mapply(riskset,x,y)),
  risksetC(X[,1],X[,2]),
  times = 100
)

N.r <- outer(X[,1], X[,2], function(x,y) mapply(riskset,x,y)) #N2
N.c <- risksetC(X[,1],X[,2])
identical(N.r,N.c)

# N1 <- c(t(N))
# N2 <- c(N)
