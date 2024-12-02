loglik <- function (param) {
  
  beta <- param[1:10]
  sigma <- param[11]
  
  X <- model.matrix(mpg ~ 0 + bs(hp, df = 10))
  pred <- X %*% beta
  
  logL <- -sum(dnorm(mpg, mean = pred, sd = sigma, log = TRUE))
  
  return(logL)
}

loglikpenal <- function (param,S.lambda) {
  
  beta <- param[1:10]
  sigma <- param[11]
  
  X <- model.matrix(mpg ~ 0 + bs(hp, df = 10))
  
  logL <- dnorm(mpg - X %*% beta, mean = 0, sd = sigma, log = TRUE)
  
  return(-sum(logL) + t(beta) %*% S.lambda %*% beta)
}

wrapper <- function(coef.vector, S.lambda=NULL, H = NULL, minusLogLik=TRUE) {
  
  beta <- coef.vector[2:10]
  # Check whether penalty is applied
  if (is.null(S.lambda)) {
    penaltyLik <- penaltyGrad <- penaltyHess <- 0
  } else {
    
    # Calculate penalty terms for log f_lambda(y,beta) Wood (2017) p.1076 
    S.lambda.eigenv <- eigen(S.lambda)$values
    
    penaltyLik <- (t(beta) %*% S.lambda %*% beta)/2
    logS.lambda <- log(prod(S.lambda.eigenv[S.lambda.eigenv > 0]))
    constant <- sum(S.lambda.eigenv == 0)*log(2*pi)/2 # Zie Wood (2016) p.1550
    
    # penaltyGrad <- t(t(coef.vector) %*% S.lambda)
    # penaltyHess <- S.lambda
  }
  
  if (is.null(H)) {
    logdetH <- 0
  } else {
    logdetH <- log(det(H))
    logS.lambda <- logS.lambda/2
  }
  
  
  # Merk op dat C++ code geïmplementeerd is voor -loglik
  sign <- ifelse(isTRUE(minusLogLik), 1, -1)
  
  # log f_lambda(y,beta)
  ll <- loglik(coef.vector) + penaltyLik - logS.lambda + logdetH - constant
  
  return(sign*ll)
  
}

EstimatePenal <- function(S, lambda.init = 1, tol = 0.01, maxiter=50) {
  
  
  lambda.new <- lambda.init # In voorbeelden van Wood (2017) is de initiele lambda = 1
  
  lldiff <- 1e10
  beta <- coef(fit.gam)
  sigma <- 3
  
  iter <- 0
  
  if (iter == 0) {print("Algorithm running...")}
  
  while (lldiff > tol & iter <= maxiter) {
    
    # Update number of iterations
    iter = iter + 1
    
    lambda <- lambda.new
    
    # Some calculations to update lambda later...
    S.lambda <- lambda*S
    S.lambda.inv <- ginv(S.lambda)
    
    # Estimate betas for given lambdas
    beta.fit <- nlm(f = wrapper,
                    p = c(beta,sigma),
                    S.lambda = S.lambda,
                    hessian = TRUE)
    
    # New betas to be used as initial values for possible next iteration
    beta <- beta.fit$estimate[1:10]
    sigma <- beta.fit$estimate[11]
    
    # Make sure that observed hessian is positive definite
    decomp <- eigen(beta.fit$hessian[2:10,2:10])
    A <- diag(abs(decomp$values))
    hessian.obs <- decomp$vectors %*% A %*% t(decomp$vectors)
    
    V <- solve(hessian.obs + S.lambda)
    
    # Update lambdas
    tr1 <- sum(diag(S.lambda.inv %*% S))
    trV1 <- sum(diag(V %*% S))
    denom1 <- t(beta) %*% S %*% beta
    
    lambda.new <- as.numeric((tr1 - trV1)/denom1)
    
    # Create new S.lambda matrix
    S.lambda.new <- lambda.new*S
    
    # Step length of update
    diff <- lambda.new - lambda
    
    # Assess whether update is an increase in the log-likelihood
    # If not, apply step length control
    l1 <- wrapper(
      coef.vector = c(beta,sigma),
      S.lambda = S.lambda.new,
      H = hessian.obs + S.lambda.new,
      minusLogLik = FALSE
    )
    
    l0 <- wrapper(
      coef.vector = c(beta,sigma),
      S.lambda = S.lambda,
      H = hessian.obs + S.lambda,
      minusLogLik = FALSE
    )
    
    # Step length control to guarantee increase in loglik
    k = 0
    
    while (l1 < l0) { 
      
      k = k + 1
      
      delta <- diff/(2^k)
      
      S.lambda.delta <- (lambda + delta)*S
      
      l1 <- wrapper(
        coef.vector = c(beta,sigma),
        S.lambda = S.lambda.delta,
        H = hessian.obs + S.lambda.delta,
        minusLogLik = FALSE
      )
      
      l0 <- wrapper(
        coef.vector = c(beta,sigma),
        S.lambda = S.lambda,
        H = hessian.obs + S.lambda,
        minusLogLik = FALSE
      )
      
    } # end of inner while loop
    
    # Final difference in loglikelihood
    lldiff <- l1 - l0
    
    # If step length control is needed, set updated lambda according to new step length
    if (k > 0) {
      lambda.new <- lambda + delta
    }
    
    # Sanity check: lambda must be positive
    if (sum(lambda.new < 0) > 0) {stop("At least 1 lambda is negative")}
    
    # Calculate loglikelihood for new lambda
    loglik.new <- wrapper(coef.vector = c(beta,sigma),
                              S.lambda = lambda.new*S,
                              minusLogLik=TRUE)
    
    # Print information while running...
    print(paste0("Iteration ", iter,
                 ": k = ", k,
                 " lambda = ", lambda.new,
                 " Likelihood increase = ", lldiff,
                 " Final -loglik = ", loglik.new))
    
  } # end of outer while loop
  
  
  
  if (iter == maxiter) {message <- "Maximum iterations reached"} else {message <- "Convergence reached"}
  
  return(list(
    beta = beta,
    lambda = lambda.new,
    iterations = iter,
    status = message,
    loglik = loglik.new))
}







data(mtcars)
head(mtcars)

attach(mtcars)

library(mgcv)
library(bbmle)

fit.gam <- gam_model <- gam(mpg ~ s(hp, k = 10 ), optimizer = 'efs')

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

test <- mle2(loglik,
             start = c(coef(fit.gam), 3),
             optim = "nlm")

test <- EstimatePenal(S, lambda.init = 20)

coef(fit.gam)

