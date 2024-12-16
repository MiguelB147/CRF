library(splines)
library(doParallel)
library(Rcpp)
library(mvtnorm)
library(progress)

source('hunan-functions.R')
sourceCpp('test.cpp')

nsim <- 500

degree = 2
df = 5
K <- 1000
unif.ub <- 5 # 5 = 20% censoring, 2.3 = 40% censoring


lambda.grid <- expand.grid(lambda1 = seq(1, 50, 2), lambda2 = seq(1, 50, 2))

set.seed(123)
data.list <- SimData(K = K, df = df, degree = degree, unif.ub = unif.ub)

S1 <- Srow(df)
S2 <- Scol(df)

ll <- ll.avg <- c()

pb <- progress_bar$new(format = "[:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                       total = nsim,
                       complete = "=",   # Completion bar character
                       incomplete = "-", # Incomplete bar character
                       current = ">",    # Current bar character
                       clear = FALSE)    # If TRUE, clears the bar when finish

cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

SList <- foreach(i=1:nrow(lambda.grid), .packages = "MASS") %dopar% {
  
  S <- lambda.grid[i,1]*S1 + lambda.grid[i,2]*S2
  Sinv <- ginv(S)
  
  list(S,Sinv)
}

beta <- foreach(i=1:nrow(lambda.grid)) %:% foreach(j=1:nsim, .packages = "mvtnorm", .combine = 'cbind') %dopar% {
  
  betas <- c(rmvnorm(1, sigma = SList[[i]][[2]], checkSymmetry = FALSE))
  
}

ll <- ll.avg <- c()
for (i in 1:nrow(lambda.grid)) {
  for (j in 1:nsim) {
    beta <- beta[[i]][,j]
    S <- SList[[i]][[1]]
    
    ll[j] <- as.numeric(loglikPenal(coef.vector = beta,
                                    degree = degree,
                                    df = df,
                                    datalist = data.list,
                                    S = S))
  }
  ll.avg[i] <- mean(ll)
}

lambda.grid[which.min(ll.avg),]

# for (i in 1:nrow(lambda.grid)) {
#   
#   pb$tick()
#   
#   S <- lambda.grid[i,1]*S1 + lambda.grid[i,2]*S2
#   Sinv <- MASS::ginv(S)
# 
#   for (j in 1:nsim) {
#     
#     
#     
#     betas <- c(rmvnorm(1, sigma = Sinv, checkSymmetry = FALSE))
#     ll[j] <- as.numeric(loglikPenal(coef.vector = betas,
#                                     degree = degree,
#                                     df = df,
#                                     datalist = data.list,
#                                     S = S))
#   }
#   ll.avg[i] <- mean(ll)
# }




# result <- foreach(i=1:3, .combine = 'rbind') %:%
#   foreach(j=1:3, .packages = c("test","splines"), .combine = 'cbind', .errorhandling = 'stop') %dopar% {
#     
#     beta <- betas[[i]][,j]
#     S <- SList[[i]][[1]]
#     
#     as.numeric(loglikPenal(coef.vector = beta,
#                            degree = degree,
#                            df = df,
#                            datalist = data.list,
#                            S = S))
#     
#   }


  
# foreach(j=1:3, .packages = c("test","splines"), .combine = 'cbind', .errorhandling = 'stop') %dopar% {
#     
#     beta <- betas[[i]][,j]
#     S <- SList[[i]][[1]]
#     
#     as.numeric(loglikPenal(coef.vector = beta,
#                            degree = degree,
#                            df = df,
#                            datalist = data.list,
#                            S = S))
#     
#   }


stopCluster(cl)

ll.avg <- rowMeans(result)

lambda.grid[which.min(ll.avg),]
  







lambda.grid$LogLikelihood <- -ll.avg
ggplot(lambda.grid, aes(lambda1, lambda2, fill = LogLikelihood)) + geom_tile() + theme_classic() + xlab(expression(lambda[1])) + ylab(expression(lambda[2]))
