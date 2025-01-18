loglik <- function (param) {
  
  beta <- param[1:10]
  sigma <- param[11]
  
  X <- model.matrix(mpg ~ 0 + bs(hp, df = 10, intercept = TRUE))
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


# EFS gebaseerd op de code van Simon Wood in het mgcv package (zie gam.fit4.r op github)
# gam.control() details in mgcv.r op github
EstimatePenal <- function(datalist, degree, S, lambda.init = c(1,1), tol = 0.001, lambda.max = exp(15)) { 
  
  tiny <- .Machine$double.eps^0.5
  
  S1 <- S[[1]]
  S2 <- S[[2]]
  
  df <- sqrt(ncol(S1))
  
  lambda.new <- lambda.init # In voorbeelden van Wood (2017) is de initiele lambda = 1
  lambda <- 0
  
  lldiff <- 1e10
  beta <- rep(1,df^2)
  
  print("Algorithm running...")
  
  score <- c()
  for (iter in 1:200) {
    
    lambda <- lambda.new
    
    # Some calculations to update lambda later...
    Sl <- lambda[1]*S1 + lambda[2]*S2
    Sl.inv <- MASS::ginv(Sl)
    
    # Estimate betas for given lambdas
    beta.fit <- nlm(f = wrapper,
                    p = beta,
                    degree = degree,
                    lambda = lambda,
                    S = S,
                    datalist = datalist,
                    hessian = FALSE)
    
    # New betas to be used as initial values for possible next iteration
    beta <- beta.fit$estimate
    
    # Make sure that hessian is positive definite
    hessian <- derivatives(coef.vector = beta, degree = degree, datalist = datalist)$hessian
    decomp <- eigen(beta.fit$hessian)
    A <- diag(abs(decomp$values))
    hessian <- decomp$vectors %*% A %*% t(decomp$vectors)
    
    # Calculate V
    V <- solve(hessian + Sl)
    
    # Calculate trSSj, trVS and bSb
    trSSj <- trVS <- bSb <- c()
    for (i in length(S)) {
      trSSj[i] <- sum(diag(Sl.inv %*% S[[i]]))
      trVS[i] <- sum(diag(V %*% S[[i]]))
      bSb[i] <- t(beta) %*% S[[i]] %*% beta
    }
    
    # Update lambdas
    update <- pmax(tiny, trSSj - trVS)/pmax(tiny, bSb) #lambda.new <- lambdaUpdate(lambda, S.lambda.inv, S, V, beta)
    update[!is.finite(update)] <- 1e6
    lambda.new <- pmin(update*lambda, lambda.max) 
    
    # Create new S.lambda matrix
    Sl.new <- lambda.new[1]*S1 + lambda.new[2]*S2
    
    # Step length of update
    max.step <- max(abs(lambda.new - lambda))
    
    # TODO Bij Wood wordt nieuwe lambda en oude Sl gebruikt. Hoe doen we dit in deze context?
    # Assess whether update is an increase in the log-likelihood
    # If not, apply step length control
    l1 <- wrapper(
      coef.vector = beta,
      degree = degree,
      lambda = lambda.new,
      S = S, # TODO Na te gaan of dit oude Sl moet zijn
      H = hessian + S.lambda, # Denk dat dit hessian + S.lambda moet zijn ipv hessian + S.lambda.new
      minusLogLik = FALSE,
      datalist = datalist
    )
    
    l0 <- wrapper(coef.vector = beta ,degree = degree, lambda = lambda, S = S, H = hessian + Sl, minusLogLik = FALSE, datalist = datalist)
    
    k = 1 # Step length
    
    if (l1 >= l0) { # Improvement
      if(max.step < 1.5) { # Consider step extension
        lambda2 <- pmin(update*lambda*k*7, exp(12))
        l3 <- wrapper(coef.vector = beta, degree = degree, lambda = lambda2, S = S,
                      H = hessian + Sl, # Denk dat dit hessian + S.lambda moet zijn ipv hessian + S.lambda.new
                      minusLogLik = FALSE,
                      datalist = datalist
        )
      } if (l3 > l1) { # Improvement - accept extension
        lambda.new <- lambda2
      } else lambda.new <- lambda.new # Accept old step
    } else { # No improvement
      while (l1 < l0) {
        k <- k/2 ## Contract step
        lambda3 <- pmin(update*lambda*k*7, lambda.max)
        l1 <- wrapper(coef.vector = beta, degree = degree, lambda = lambda3, S = S, H = hessian + Sl, minusLogLik = FALSE, datalist = datalist)
      }
    }
    
    # If step length control is needed, update lambda accordingly
    if(k < 1) {
      lambda.new <- lambda3
      max.step <- max(abs(lambda.new - lambda))}
    
    #save loglikelihood value
    score[iter] <- l1
    
    # Sanity check: lambda must be positive
    if (sum(lambda.new < 0) > 0) {stop("At least 1 lambda is negative")}
    
    # Break procedure if REML change and step size are too small
    if (iter > 3 && max.step < 1 && max(abs(diff(score[(iter-3):iter]))) < .1) break
    # Or break is likelihood does not change
    if (l1 == l0) break
    
    # Calculate loglikelihood for new lambda
    loglik.new <- wrapper(coef.vector = beta,
                          degree = degree,
                          datalist = datalist,
                          S.lambda = lambda.new[1]*S1 + lambda.new[2]*S2,
                          minusLogLik=TRUE)
    
    # Print information while running...
    print(paste0("Iteration ", iter,
                 ": k = ", k,
                 " lambda1 = ", lambda.new[1],
                 " lambda2 = ", lambda.new[2],
                 " Score increase = ", score[iter] - score[iter-1],
                 " REML = ", score[iter]))
    
  } # End of for loop
  
  
  return(list(
    beta = beta,
    lambda = lambda.new,
    iterations = iter,
    status = message,
    history = score))
}