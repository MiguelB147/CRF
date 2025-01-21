
# param = c(beta, sig2)
loglikpenal <- function (param, X, Sl = NULL, H = NULL, minusloglik = TRUE) {
  
  beta <- param
  df <- length(param)
  
  if (is.null(Sl)) {
    logSl <- 0
    logdetH <- 0
    penalty <- 0
  } else {
    ev <- eigen(Sl)$values
    logSl <- log(prod(ev[ev > 0]))
    penalty <- (t(beta[-1]) %*% Sl %*% beta[-1])/2
  }
  
  if (!is.null(H)) {
    logdetH <- log(det(H + Sl))
  } else logdetH <- 0
  
  sign <- ifelse(isTRUE(minusloglik), -1, 1)
  
  ll <- sum(dnorm(mpg, mean = X %*% beta, sd = sqrt(9.54816882), log = TRUE)) - penalty  + logSl/2 - logdetH/2
  
  # logL <- dnorm(mpg - X %*% beta, mean = 0, sd = sigma, log = TRUE)
  
  return(sign*ll)
}


EstimatePenal <- function(S, lambda.init = 1, tol = 0.001, lambda.max = exp(15)) { 
  
  tiny <- .Machine$double.eps^0.5

  df <- ncol(S) + 1
  
  # Positioning of the boundary knots
  xl <- min(hp); xu <- max(hp); xr <- xu - xl
  xl <- xl - 0.001*xr; xu <- xu + 0.001*xr

  X <- model.matrix(mpg ~ 0 + splines::bs(hp, df = df, Boundary.knots = c(xl,xu), intercept = TRUE))
  
  lambda.new <- lambda.init
  
  # Initial values
  init.fit <- gam(mpg ~ s(hp, k = df), optimizer = "efs")
  beta <- init.fit$coef
  

  
  print("Extended Fellner-Schall method:")
  
  score <- c()
  for (iter in 1:200) {
    
    lambda <- lambda.new
    
    # Some calculations to update lambda later...
    Sl <- lambda*S
    Sl.inv <- MASS::ginv(Sl)
    
    # Estimate betas for given lambdas
    beta.fit <- optim(par = beta,
                      fn = loglikpenal,
                      X = X,
                      Sl = Sl,
                      hessian = TRUE)
    
    beta <- beta.fit$par
    
    V <- solve(beta.fit$hessian[-1,-1] + Sl)
    
    trSSj <- sum(diag(Sl.inv %*% S))
    trVS <- sum(diag(V %*% S))
    bSb <- t(beta[-1]) %*% S %*% beta[-1]
    
    
    # Update lambdas
    update <- pmax(tiny, trSSj - trVS)/pmax(tiny, bSb)
    update[!is.finite(update)] <- 1e6
    lambda.new <- pmin(update*lambda, lambda.max) 
    
    # Create new S.lambda matrix
    Sl.new <- lambda.new*S
    
    # Step length of update
    max.step <- max(abs(lambda.new - lambda))
    
    # Assess whether update is an increase in the log-likelihood
    # If not, apply step length control
    l1 <- loglikpenal(param = beta, X = X, Sl = Sl.new, H = beta.fit$hessian[-1,-1], minusloglik = FALSE)
    l0 <- loglikpenal(param = beta, X = X, Sl = Sl, H = beta.fit$hessian[-1,-1], minusloglik = FALSE)
    
    k = 1 # Step length
    
    if (l1 >= l0) { # Improvement
      if(max.step < 1.5) { # Consider step extension
        lambda2 <- pmin(update*lambda*k*2, exp(12))
        Sl2 <- lambda2*S
        l3 <- loglikpenal(param = beta, X = X, Sl = Sl2, H = beta.fit$hessian[-1,-1], minusloglik = FALSE)
        
      if (l3 > l1) { # Improvement - accept extension
        lambda.new <- lambda2
      } # No improvement - Accept old step
    }
      } else { # No improvement
        lk <- l1
      while (lk < l0 && k > 0.001) { # Don't contract too much since the likelihood does not need to increase
        k <- k/2 ## Contract step
        lambda3 <- pmin(update*lambda*k, lambda.max)
        Sl.new <- lambda3*S
        lk <- loglikpenal(param = beta, X = X, Sl = Sl.new, H = beta.fit$hessian[-1,-1], minusloglik = FALSE)
      }
    }
    
    # If step length control is needed, update lambda accordingly
    if (k < 1 & k > 0.001) {
      lambda.new <- lambda3
      l1 <- lk
      max.step <- max(abs(lambda.new - lambda))
      } else k <- 1
    
    #save loglikelihood value
    score[iter] <- l1
    
    # Print information while running...
    print(paste0("Iteration ", iter,
                 ": k = ", k,
                 " lambda = ", lambda.new,
                 " Score increase = ", score[iter] - score[iter-1],
                 " REML = ", score[iter]))
    
    # Break procedure if REML change and step size are too small
    if (iter > 3 && max.step < 1 && max(abs(diff(score[(iter-3):iter]))) < .5) break
    # Or break is likelihood does not change
    if (l1 == l0) break
    

    
  } # End of for loop
  
  return(list(
    beta = beta,
    lambda = lambda.new,
    iterations = iter,
    history = score))
}
