loglik <- function (param, X) {
  
  beta <- param[1:10]
  sigma <- param[11]
  
  pred <- X %*% beta
  
  logL <- -sum(dnorm(mpg, mean = pred, sd = sigma, log = TRUE))
  
  return(logL)
}

# param = c(beta, sig2)
loglikpenal <- function (param, Sl = NULL, minusloglik = TRUE, REML = FALSE) {
  

  df <- length(param) - 1
  beta <- param[1:df]
  sig2 <- param[df+1]
  
  if (is.null(Sl)) {
    logSl <- 0
    logdetXtX <- 0
  } else {
    ev <- eigen(Sl/sig2)$values
    logSl <- log(prod(ev[ev > 0]))
  }
  
  X <- model.matrix(mpg ~ 0 + bs(hp, df = df))
  
  if (REML) {
    XtX <- crossprod(X)
    logdetXtX <- log(det(XtX/sig2 + Sl/sig2))
  } else logdetXtX <- 0
  
  sign <- ifelse(isTRUE(minusloglik), -1, 1)
  
  ll <- - (norm(mpg - X %*% beta, type = "2")^2 + t(beta) %*% Sl %*% beta)/(2*sig2) + logSl/2 - logdetXtX
  
  # logL <- dnorm(mpg - X %*% beta, mean = 0, sd = sigma, log = TRUE)
  
  return(sign*ll)
}


EstimatePenal <- function(S, lambda.init = 1, tol = 0.001, lambda.max = exp(15)) { 
  
  tiny <- .Machine$double.eps^0.5

  df <- ncol(S)
  
  X <- model.matrix(mpg ~ 0 + bs(hp, df = df, intercept = TRUE))
  n <- nrow(X)
  
  lambda.new <- lambda.init
  
  # Initial values
  init.fit <- gam(mpg ~ s(hp, k = df), optimizer = "efs")
  beta <- c(init.fit$coef, init.fit$sig2)
  

  
  print("Algorithm running...")
  
  score <- c()
  for (iter in 1:200) {
    
    lambda <- lambda.new
    
    # Some calculations to update lambda later...
    Sl <- lambda*S
    Sl.inv <- MASS::ginv(Sl)
    
    # Estimate betas for given lambdas
    beta.fit <- optim(par = beta,
                      fn = loglikpenal,
                      Sl = Sl,
                      hessian = FALSE)
    
    coef <- beta.fit$par[1:df]
    
    XtX <- crossprod(X)
    XtXSl.inv <- solve(XtX + Sl)
    
    
    trSSj <- sum(diag(Sl.inv %*% S))
    trVS <- sum(diag(XtXSl.inv %*% S))
    bSb <- t(coef) %*% S %*% coef
    
    # sig2 <- norm(mpg - X %*% coef, type = "2")^2 / (n - sum(diag(XtXSl.inv %*% XtX)))
    sig2 <- beta.fit$par[df+1]
    
    # New betas to be used as initial values for possible next iteration
    beta <- c(coef,sig2)
    
    # Update lambdas
    update <- pmax(tiny, trSSj - trVS)/pmax(tiny, bSb) #lambda.new <- lambdaUpdate(lambda, S.lambda.inv, S, V, beta)
    update[!is.finite(update)] <- 1e6
    lambda.new <- pmin(update*lambda*sig2, lambda.max) 
    
    # Create new S.lambda matrix
    Sl.new <- lambda*S
    
    # Step length of update
    max.step <- max(abs(lambda.new - lambda))
    
    # Assess whether update is an increase in the log-likelihood
    # If not, apply step length control
    l1 <- loglikpenal(param = beta, Sl = Sl.new, minusloglik = FALSE, REML = TRUE)
    l0 <- loglikpenal(param = beta, Sl = Sl, minusloglik = FALSE, REML = TRUE)
    
    k = 1 # Step length
    
    if (l1 >= l0) { # Improvement
      if(max.step < 1.5) { # Consider step extension
        lambda2 <- pmin(update*lambda*k*7*sig2, exp(12))
        Sl2 <- lambda2*S
        l3 <- loglikpenal(param = beta, Sl = Sl2, minusloglik = FALSE, REML = TRUE)
        
      if (l3 > l1) { # Improvement - accept extension
        lambda.new <- lambda2
      } else lambda.new <- lambda.new # No improvement - Accept old step
    }
      } else { # No improvement
      while (l1 < l0) {
        k <- k/2 ## Contract step
        lambda.new <- pmin(update*lambda*k*sig2, lambda.max)
        Sl.new <- lambda.new*S
        l1 <- loglikpenal(param = beta, Sl = Sl.new, minusloglik = FALSE, REML = TRUE)
      }
    }
    
    # If step length control is needed, update lambda accordingly
    if(k < 1) {
      lambda.new <- lambda3
      max.step <- max(abs(lambda.new - lambda))
      }
    
    #save loglikelihood value
    score[iter] <- l1
    
    # Break procedure if REML change and step size are too small
    if (iter > 3 && max.step < 1 && max(abs(diff(score[(iter-3):iter]))) < .1) break
    # Or break is likelihood does not change
    if (l1 == l0) break
    
    # Print information while running...
    print(paste0("Iteration ", iter,
                 ": k = ", k,
                 " lambda = ", lambda.new,
                 " Score increase = ", score[iter] - score[iter-1],
                 " REML = ", score[iter]))
    
  } # End of for loop
  
  return(list(
    beta = coef,
    sig2 = sig2,
    lambda = lambda.new,
    iterations = iter,
    history = score))
}
