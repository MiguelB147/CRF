
WoodSpline <- function(t, dim, degree = 3, type = NULL, quantile = FALSE, scale = TRUE, m2 = degree-1) {
  
  nk <- dim - degree + 1 # Number of "interior" knots (internal + boundary)
  
  xl <- min(t)
  xu <- max(t)
  xr <- xu - xl
  xl <- xl-xr*0.001; xu <- xu+xr*0.001
  dx <- (xu-xl)/(nk-1)
  k <- seq(xl-dx*degree,xu+dx*degree,length=nk+2*degree) # Vector of knots
  if(quantile) {
    k.int <- quantile(t, probs = seq(0, 1, length = nk))[-c(1, nk)]
    k[(degree+2):(length(k)-(degree+1))] <- k.int
  }

  
  X <- splines::splineDesign(k, t, degree+1)
  
  if (is.null(type)) {S <- D1 <- NULL} else if (type == "bs") {
    
    pord <- degree - m2
    k0 <- k[(degree+1):(degree+nk)]
    h <- diff(k0)
    h1 <- rep(h/pord, each = pord)
    k1 <- cumsum(c(k0[1],h1))
    
    D <- splines::splineDesign(k,k1,derivs = m2)
  
    P <- solve(matrix(rep(seq(-1,1,length=pord+1),pord+1)^rep(0:pord,each=pord+1),pord+1,pord+1))
    i1 <- rep(1:(pord+1),pord+1)+rep(1:(pord+1),each=pord+1) ## i + j
    H <- matrix((1+(-1)^(i1-2))/(i1-1),pord+1,pord+1)
    W1 <- t(P)%*%H%*%P
    h <- h/2 ## because we map integration interval to to [-1,1] for maximum stability
    ## Create the non-zero diagonals of the W matrix... 
    ld0 <- rep(sdiag(W1),length(h))*rep(h,each=pord+1)
    i1 <- c(rep(1:pord,length(h)) + rep(0:(length(h)-1) * (pord+1),each=pord),length(ld0))
    ld <- ld0[i1] ## extract elements for leading diagonal
    i0 <- 1:(length(h)-1)*pord+1
    i2 <- 1:(length(h)-1)*(pord+1)
    ld[i0] <- ld[i0] + ld0[i2] ## add on extra parts for overlap
    B <- matrix(0,pord+1,length(ld))
    B[1,] <- ld
    for (k in 1:pord) { ## create the other diagonals...
      diwk <- sdiag(W1,k) ## kth diagonal of W1
      ind <- 1:(length(ld)-k)
      B[k+1,ind] <- (rep(h,each=pord)*rep(c(diwk,rep(0,k-1)),length(h)))[ind]  
    }
    ## ... now B contains the non-zero diagonals of W
    B <- mgcv::bandchol(B) ## the banded cholesky factor.
    ## Pre-Multiply D by the Cholesky factor...
    D1 <- B[1,]*D
    for (k in 1:pord) {
      ind <- 1:(nrow(D)-k)
      D1[ind,] <- D1[ind,] + B[k+1,ind] * D[ind+k,]
    }
    S <- crossprod(D1)
    
    # G <- splines::splineDesign(k,k1,derivs = 2)
    # 
    # P <- H <- matrix(0, ncol = p+1, nrow = p+1)
    # for (i in 1:(p+1)) {
    #   for (j in 1:(p+1)) {
    #     P[i,j] <- (-1 + 2*(i-1)/p)^j
    #     H[i,j] <- (1 + (-1)^(i+j-2))/(i+j-1)
    #   }
    # }
    # Wtilde <- t(solve(P)) %*% H %*% solve(P)
    # 
    # W <- matrix(0, ncol = length(k1), nrow = length(k1))
    # for (q in 1:length(h)) {
    #   for (i in 1:(p+1)) {
    #     for (j in 1:(p+1)) {
    #       W[i+p*q-p,j+p*q-p] <- h[q]*Wtilde[i,j]/2
    #     }
    #   }
    # }
    # 
    # S <- t(G) %*% W %*% G
    # 
    # # Banded cholesky decomposition
    # R <- chol(W) # t(R) %*% R = W
    # D <- R %*% G # t(D) %*% D = S
  } else if (type == "ps") {
    D1 <- diff(diag(dim), differences = m2)
    S <- crossprod(D1)
  }
  
  # if(repara) {
  #   qrX <- qr(X)
  #   R <- qr.R(qrX)
  #   Q <- qr.Q(qrX)
  #   Rinv <- solve(R)
  #   eigenS <- eigen(t(Rinv) %*% S %*% Rinv, symmetric = TRUE)
  #   Sprime <- diag(eigenS$values)
  #   Xprime <- Q %*% eigenS$vectors
  #   
  #   S <- Sprime
  #   X <- Xprime
  # }
  
  if (scale) {
    maXX <- norm(X,type="I")^2
    maS <- norm(S)/maXX
    S <- S/maS
  } else S.scale <- NULL

  
  return(list(X = X, knots = k, S = S, D = D1, S.scale = maS))
}


# param = c(beta, sig2)
loglikpenal <- function (beta, X, Sl = NULL, H = NULL, minusloglik = TRUE) {
  
  
  if (is.null(Sl)) {
    logSl <- 0
    logdetH <- 0
    penalty <- 0
  } else {
    ev <- eigen(Sl)$values
    ev[abs(ev) < 1e-6] <- 0
    logSl <- sum(log(ev[ev > 0]))
    penalty <- (t(beta) %*% Sl %*% beta)/2
  }
  
  if (!is.null(H)) {
    logdetH <- log(det(H))
  } else logdetH <- 0
  
  sign <- ifelse(isTRUE(minusloglik), -1, 1)
  
  ll <- sum(dnorm(mtcars$mpg, mean = X %*% beta, sd = sqrt(9.54816882), log = TRUE)) - penalty  + logSl/2 - logdetH/2
  
  # logL <- dnorm(mpg - X %*% beta, mean = 0, sd = sigma, log = TRUE)
  
  return(sign*ll)
}

efsud.fit <- function(start, X, Sl) {
  beta.fit <- nlm(loglikpenal, p = start, X = X, Sl = Sl, hessian = TRUE)
  # H <- beta.fit$hessian
  H <- hessian(X,Sl)
  beta <- beta.fit$estimate
  fit <-  loglikpenal(beta = beta,
                      X = X,
                      Sl = Sl, H = H,
                      minusloglik = FALSE)
  
  return(list(beta = beta, hessian = H, REML = fit))
}

hessian <- function(X, Sl = NULL) {
  
  dim <- ncol(X)
  
  H <- matrix(ncol = dim, nrow = dim)
  for (i in 1:dim) {
    for (j in i:dim) {
      H[i,j] <- H[j,i] <- (1/9.54816882)*sum(X[,i]*X[,j])
    }
  }
  
  if (!is.null(Sl)) H <- H + Sl
  
  return(H)
}


EstimatePenal <- function(dim = 10, lambda.init = 10, type = "bs", quantile = FALSE, scale = TRUE, tol = 0.001, lambda.max = exp(15), step.control = TRUE) { 
  
  tiny <- .Machine$double.eps^0.5
  

  model <- WoodSpline(t = mtcars$hp, dim = dim, degree = 3, type = type, scale = scale, quantile = quantile)
  S <- model$S
  X <- model$X

  lambda.new <- lambda.init # In voorbeelden van Wood (2017) is de initiele lambda = 1
  
  fit <- efsud.fit(start = rep(1,dim), X = X, Sl = lambda.init*S)
  k <- 1
  score <- rep(0, 200)
  for (iter in 1:200) {
    
    l0 <- fit$REML
    
    lambda <- lambda.new
    
    # Some calculations to update lambda later...
    Sl <- lambda*S
    Sl.inv <- MASS::ginv(Sl)
    
    # decomp <- eigen(hessian)
    # A <- diag(abs(decomp$values))
    # hessian <- decomp$vectors %*% A %*% t(decomp$vectors)
    
    
    # Update ----
    
    # Calculate V
    V <- solve(fit$hessian)
    
    trSSj <- sum(diag(Sl.inv %*% S))
    trVS <- sum(diag(V %*% S))
    bSb <- t(fit$beta) %*% S %*% fit$beta
    
    # Update lambdas
    a <- pmax(tiny, trSSj - trVS)
    update <- a/pmax(tiny, bSb)
    update[a==0 & bSb==0] <- 1
    update[!is.finite(update)] <- 1e6
    lambda.new <- pmin(update*lambda, lambda.max)
    
    # Step length of update
    max.step <- max(abs(lambda.new - lambda))
    
    # Create new S.lambda matrix
    Sl.new <- lambda.new*S
    
    fit <- efsud.fit(start = fit$beta, X = X, Sl = Sl.new)
    l1 <- fit$REML
    
    # Start of step control ----
    if (step.control) {
      
      if (l1 > l0) { # Improvement
        if(max.step < 1) { # Consider step extension
          lambda2 <- pmin(lambda*update^(k*2), exp(12))
          fit2 <- efsud.fit(start = fit$beta, X = X, Sl = lambda2*S)
          l2 <- fit2$REML
          if (l2 > l1) { # Improvement - accept extension
            lambda.new <- lambda2
            l1 <- l2
            fit <- fit2
            k <- k*2
          }
        }
      } else { # No improvement
        lk <- l1
        lambda3 <- lambda.new
        while (lk < l0 && k > 1) { # Don't contract too much since the likelihood does not need to increase k > 0.001
          k <- k/2 ## Contract step
          lambda3 <- pmin(lambda*update^k, lambda.max)
          fit <- efsud.fit(start = fit$beta, X = X, Sl = lambda3*S)
          lk <- fit$REML
          
          # k <- k + 1
          # diff <- diff/2
          # lambda3 <- lambda + diff
          # lk <- wrapper(coef.vector = beta, degree = degree,
          #               # Sl = lambda3[1]*S1 + lambda3[2]*S2,
          #               Sl = lambda3*S,
          #               H = hessian, minusLogLik = FALSE, datalist = datalist)
        }
        lambda.new <- lambda3
        l1 <- lk
        max.step <- max(abs(lambda.new - lambda))
        if (k < 1) k <- 1
      }
    } # end of step length control
    
    # save loglikelihood value
    score[iter] <- l1
    
    # Break procedures ----
    
    # Break procedure if REML change and step size are too small
    if (iter > 3 && max.step < 0.5 && max(abs(diff(score[(iter-3):iter]))) < 1) {print("REML not changing"); break}
    # Or break is likelihood does not change
    if (l1 == l0) {print("Loglik not changing"); break}
    # Stop if loglik is not changing
    # if (iter==1) old.ll <- fit$ll else {
    #   if (abs(old.ll-fit$ll)<eps*abs(fit$ll)) {print("Loglik not changing"); break}  # *100
    #   old.ll <- fit$ll
    # }
    
    # Print information while running...
    print(paste0("Iteration ", iter,
                 ": k = ", k,
                 " lambda = ", round(lambda.new,4),
                 " REML = ", score[iter]))
    
  } # End of for loop
  
  if (iter < 200) print("Converged") else print("Number of iterations is too small")
  
  
  return(list(
    beta = fit$beta,
    lambda = lambda.new,
    iterations = iter,
    history = score[1:iter],
    X = X))
}
