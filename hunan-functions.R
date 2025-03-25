theta.frank <- function(x,y,alpha) {
  A <- (alpha-1)*log(alpha)*alpha^(2-exp(-x)-exp(-y))
  B <- (alpha^(1-exp(-x))-alpha)*(alpha^(1-exp(-y))-alpha)
  C <- -1 + exp(-x) + exp(-y) + log(1+ (alpha^(1-exp(-x))-1)*(alpha^(1-exp(-y))-1)/(alpha-1), base = alpha)
  return(A*C/B)
}

WoodSpline <- function(t, dim, degree = 3, type = NULL, quantile = FALSE, scale = TRUE, m2 = degree-1) {
  
  # Create knot sequence for spline ----
  nk <- dim - degree + 1 # Number of "interior" knots (internal + boundary)
  
  xl <- min(t)
  xu <- max(t)
  xr <- xu - xl
  xl <- xl-xr*0.001; xu <- xu+xr*0.001
  dx <- (xu-xl)/(nk-1)
  k <- seq(xl-dx*degree,xu+dx*degree,length=nk+2*degree) # Vector of knots
  if (quantile) {
    k.int <- quantile(t, probs = seq(0, 1, length = nk))[-c(1, nk)]
    k[(degree+2):(length(k)-(degree+1))] <- k.int
  }
  
  X <- splines::splineDesign(k, t, degree+1)
  
  # Create penalty matrix S = t(D1) %*% D1 if necessary ----
  if (is.null(type)) {S <- D1 <- NULL}
  else if (type == "bs") {
    
    ## Integrated squared derivative penalty ----
    
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
    
    
    # p <- degree - 2
    # k0 <- k[(degree+1):(degree+nk)]
    # h <- diff(k0)
    # h1 <- rep(h/p, each = p)
    # k1 <- cumsum(c(k0[1],h1))
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
    
    # Banded cholesky decomposition
    # R <- chol(W) # t(R) %*% R = W
    # D <- R %*% G # t(D) %*% D = S
    
  } else if (type == "ps") {
    ## Discrete penalty ----
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
  
  # Scaling the penalty matrix S ----
  if (scale) {
    maXX <- norm(X,type="I")^2
    maS <- norm(S)/maXX
    S <- S/maS
    D1 <- D1/sqrt(maS)
  } else maS <- NULL
  
  return(list(X = X, knots = k, S = S, D = D1, S.scale = maS))
}

# See Reiss et al. (2014) 
WoodPenalty <- function(object1, object2) {
  
  df1 <- ncol(object1$X)
  df2 <- ncol(object2$X)
  
  D1 <- object1$D %x% diag(rep(1,df2))
  D2 <- diag(rep(1,df1)) %x% object2$D
  
  S1 <- crossprod(D1)
  S2 <- crossprod(D2)
  
  # if (type == "bs") {
  #   
  #   D1 <- object1$D %x% diag(rep(1,df2))
  #   D2 <- diag(rep(1,df1)) %x% object2$D
  #   
  #   S1 <- crossprod(D1)
  #   S2 <- crossprod(D2)
  #   
  # } else if (type == "ps") {
  # 
  #   
  #   S1 <- crossprod(diff(diag(df1^2), lag = df1, differences = degree-1))
  #   
  #   P <- kronecker(diag(df2), diff(diag(df2), differences = degree-1))
  #   S2 <- crossprod(P)
  #   
  # } else stop("Penalty type is not correctly specified")
  
  return(list(S1 = S1,S2 = S2))
}

WoodTensor <- function(X1, X2, coef.vector) {
  
  coef.matrix <- matrix(coef.vector, ncol = ncol(X1), byrow = FALSE)
  
  spline <- X1 %*% coef.matrix %*% t(X2)
  
  return(spline)
}



tensor <- function(t1, t2, coef.vector, df, degree, knots) {
  
  coef.matrix <- matrix(coef.vector, ncol = df, byrow = FALSE)
  
  if (is.null(knots)) {
    
    # spline1 <- bs(t1, df = df, degree = degree)
    # B1 <- matrix(as.numeric(spline1), nr = nrow(spline1))
    # spline2 <- bs(t2, df = df, degree = degree)
    # B2 <- matrix(as.numeric(spline2), nr = nrow(spline2))
    
    B1 <- model.matrix(~ 0 + bs(t1, df = df, degree = degree))
    B2 <- model.matrix(~ 0 + bs(t2, df = df, degree = degree))

    spline12 <- B1 %*% coef.matrix %*% t(B2)
    
  } else {
    
    boundary1 <- knots[c(1,nrow(knots)),1]
    boundary2 <- knots[c(1,nrow(knots)),2]
    int1 <- knots[2:(nrow(knots)-1),1]
    int2 <- knots[2:(nrow(knots)-1),2]
    
    # spline1 <- bs(t1, knots = int1, degree = degree, Boundary.knots = boundary1)
    # B1 <- matrix(as.numeric(spline1), nr = nrow(spline1))
    # spline2 <- bs(t2, knots = int2, degree = degree, Boundary.knots = boundary2)
    # B2 <- matrix(as.numeric(spline2), nr = nrow(spline2))
    
    B1 <- model.matrix(~ 0 + bs(t1, knots = int1, degree = degree, Boundary.knots = boundary1))
    B2 <- model.matrix(~ 0 + bs(t2, knots = int2, degree = degree, Boundary.knots = boundary2))
    spline12 <- B1 %*% coef.matrix %*% t(B2)
    
  }
  
  return(spline12)
}

# tensor.deriv <- function(t1, t2, df, degree, knots) {
# 
#     boundary1 <- knots[c(1,df),1]
#     boundary2 <- knots[c(1,df),2]
#     int1 <- knots[2:(df-1),1]
#     int2 <- knots[2:(df-1),2]
#     
#     spline1 <- bs(t1, knots = int1, degree = degree, Boundary.knots = boundary1)
#     B1 <- matrix(as.numeric(spline1), nr = nrow(spline1))
#     spline2 <- bs(t2, knots = int2, degree = degree, Boundary.knots = boundary2)
#     B2 <- matrix(as.numeric(spline2), nr = nrow(spline2))
#     spline12 <- B1 %*% t(B2)
#   
#   return(spline12)
# }

# polynomial <- function(t1,t2,coef.vec) {
#   
#   logtheta <- coef.vec[1] + coef.vec[2]*t1 + coef.vec[3]*t2 +
#     coef.vec[4]*t1^2 + coef.vec[5]*t2^2 + coef.vec[6]*t1*t2 + 
#     coef.vec[7]*(t1^2)*t2 + coef.vec[8]*t1*(t2^2) +
#     coef.vec[9]*t1^3 + coef.vec[10]*t2^3
#   
#   return(logtheta)
#   
# }
# 
# riskset <- function(x, y) {
#   N <- sum(1*(X[,1] >= x & X[,2] >= y))
#   return(N)
# }
# 
# loglikdf <- function(coef.vector, degree, df, knots) {
# 
#   logtheta2 <- tensor(X[,1], X[,2], degree = degree, coef.vector = coef.vector, df = df, knots = knots)
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
# 
# loglikdf2 <- function(coef.vector, degree, df, knots) {
#   
#   
#   logtheta2 <- tensor(X[,1],X[,2], degree = degree, coef.vector = coef.vector, df = df, knots = knots)
#   logtheta1 <- t(logtheta2)
#   
#   B1 <- c(I5*logtheta1)[N1 > 0]
#   B2 <- c(I6*logtheta2)[N2 > 0]
#   
#   C1 <- c(N1 + I2*(exp(logtheta1)-1))[N1 > 0]
#   C2 <- c(N2 + I4*(exp(logtheta2)-1))[N2 > 0]
#   
#   L1 <- sum(A1*(B1 - log(C1)))
#   L2 <- sum(A2*(B2 - log(C2)))
#   
#   return(-(L1+L2))
# }

loglikCpp <- function(coef.vector, degree, df, datalist) {
  
  logtheta2 <- tensor(datalist$X[,1], datalist$X[,2], degree = degree, coef.vector = coef.vector, df = df, knots = datalist$knots)
  
  L1 <- logLikC(riskset = t(datalist$riskset),
                logtheta = t(logtheta2),
                delta = t(datalist$delta.prod),
                I1 = datalist$I1, I2 = datalist$I2, I3 = datalist$I5)
  L2 <- 0
  # L2 <- logLikC(riskset = datalist$riskset,
  #               logtheta = logtheta2,
  #               delta = datalist$delta.prod,
  #               I1 = t(datalist$I2), I2 = t(datalist$I1), I3 = datalist$I6)
  
  return(L1+L2)
}


# loglikAsym <- function(coef.vector, degree, df, datalist) {
#   
#   logtheta2 <- tensor(datalist$X[,1], datalist$X[,2], degree = degree, coef.vector = coef.vector, df = df, knots = datalist$knots)
#   
#   L1 <- logLikC(riskset = t(datalist$riskset),
#                  logtheta = t(logtheta2),
#                  delta = t(datalist$delta.prod),
#                  I1 = datalist$I1, I2 = datalist$I2, I3 = datalist$I5)
#   L2 <- 0
#   # L2 <- logLikC2(riskset = datalist$riskset,
#   #                logtheta = logtheta2,
#   #                delta = datalist$delta.prod,
#   #                I1 = t(datalist$I2), I2 = t(datalist$I1), I3 = datalist$I6)
#   
#   return(-L1 - L2)
# }

loglikPenal <- function(coef.vector, degree, df, datalist, lambda = c(0,0), S = NULL) {
  
  logtheta2 <- tensor(datalist$X[,1], datalist$X[,2], degree = degree, coef.vector = coef.vector, df = df, knots = datalist$knots)
  
  L1 <- logLikC(riskset = t(datalist$riskset), # LogLikC2 = asym, LogLikC = full
                logtheta = t(logtheta2),
                delta = t(datalist$delta.prod),
                I1 = datalist$I1, I2 = datalist$I2, I3 = datalist$I5)
  L2 <- 0
  # L2 <- logLikC(riskset = datalist$riskset,
  #               logtheta = logtheta2,
  #               delta = datalist$delta.prod,
  #               I1 = t(datalist$I2), I2 = t(datalist$I1), I3 = datalist$I6)
  
  
  if ( is.null(S)) {Penalty <- 0} else {
    Sl <- lambda[1]*S[[1]] + lambda[2]*S[[2]]
    Penalty <- t(coef.vector) %*% Sl %*% coef.vector
    }
  
  return(L1+L2-Penalty/2)
}


derivatives <- function(coef.vector, degree, datalist, Sl = NULL, gradient = FALSE, hessian = TRUE) {
  
  df <- sqrt(length(coef.vector))
  
  # Tensor product spline
  logtheta2 <- tensor(datalist$X[,1], datalist$X[,2], degree = degree, coef.vector = coef.vector, df = df, knots = datalist$knots)
  
  M <- diag(df^2)
  
  # List of gradient matrices for every spline coefficient
  deriv <- apply(M, 2, tensor,
                 t1 = datalist$X[,1], t2 = datalist$X[,2], degree = degree, df = df, knots = datalist$knots,
                 simplify = FALSE)
  
  if (isTRUE(gradient)) {
    
    gradient <- gradientC(riskset = datalist$riskset,
                          logtheta = logtheta2,
                          df = df,
                          delta = datalist$delta.prod,
                          deriv = deriv,
                          I1 = datalist$I1,
                          I2 = datalist$I2,
                          I3 = datalist$I5,
                          I4 = datalist$I6) # gradientC returns vector of derivatives of -loglik
    
  } else {gradient <- NA}

  if (isTRUE(hessian)) {
    
    hessian <- hessianC(riskset = t(datalist$riskset),
                        logtheta = t(logtheta2),
                        deriv = deriv,
                        df = df,
                        delta = t(datalist$delta.prod),
                        I1 = datalist$I1,
                        I2 = datalist$I2,
                        I3 = datalist$I5) # hessianC returns matrix of second derivatives of -loglik
    
  } else {hessian <- NA}
  
  if (!is.null(Sl)) {
    gradient <- gradient + t(coef.vector) %*% Sl
    hessian <- hessian + Sl
  }
  
  
  return(list(gradient = gradient, hessian = hessian))
  # return(gradient)
}

derivatives2 <- function(coef.vector, X1, X2, datalist, Sl = NULL, gradient = FALSE, hessian = TRUE) {
  
  df <- ncol(X1)
  
  # Tensor product spline
  logtheta2 <- WoodTensor(X1, X2, coef.vector = coef.vector)
  
  M <- diag(df^2)
  
  # List of gradient matrices for every spline coefficient
  deriv <- apply(M, 2, WoodTensor,
                 X1 = X1, X2 = X2,
                 simplify = FALSE)
  
  if (isTRUE(gradient)) {
    
    gradient <- gradientC(riskset = datalist$riskset,
                          logtheta = logtheta2,
                          df = df,
                          delta = datalist$delta.prod,
                          deriv = deriv,
                          I1 = datalist$I1,
                          I2 = datalist$I2,
                          I3 = datalist$I5,
                          I4 = datalist$I6) # gradientC returns vector of derivatives of -loglik
    
  } else {gradient <- NA}
  
  if (isTRUE(hessian)) {
    
    hessian <- hessianC(riskset = t(datalist$riskset),
                        logtheta = t(logtheta2),
                        deriv = deriv,
                        df = df,
                        delta = t(datalist$delta.prod),
                        I1 = datalist$I1,
                        I2 = datalist$I2,
                        I3 = datalist$I5) # hessianC returns matrix of second derivatives of -loglik
    
  } else {hessian <- NA}
  
  if (!is.null(Sl)) {
    gradient <- gradient + t(coef.vector) %*% Sl
    hessian <- hessian + Sl
  }
  
  
  return(list(gradient = gradient, hessian = hessian))
  # return(gradient)
}

Score <- function(coef.vector, degree, datalist, Sl = NULL) {
  
  df <- sqrt(length(coef.vector))
  
  # Tensor product spline
  logtheta2 <- tensor(datalist$X[,1], datalist$X[,2], degree = degree, coef.vector = coef.vector, df = df, knots = datalist$knots)
  
  M <- diag(df^2)
  
  # List of gradient matrices for every spline coefficient
  deriv <- apply(M, 2, tensor,
                 t1 = datalist$X[,1], t2 = datalist$X[,2], degree = degree, df = df, knots = datalist$knots,
                 simplify = FALSE)
  
  
  gradient <- gradientC(riskset = datalist$riskset,
                        logtheta = logtheta2,
                        df = df,
                        delta = datalist$delta.prod,
                        deriv = deriv,
                        I1 = datalist$I1,
                        I2 = datalist$I2,
                        I3 = datalist$I5,
                        I4 = datalist$I6) # gradientC returns vector of derivatives of -loglik
  
  if (!is.null(Sl)) {
    penalty <- t(coef.vector) %*% Sl
  } else penalty <- 0
  
  return(gradient + penalty)
}

Score2 <- function(coef.vector, X1, X2, datalist, Sl = NULL) {
  
  # Tensor product spline
  logtheta2 <- WoodTensor(X1 = X1, X2 = X2, coef.vector = coef.vector)
  
  df <- ncol(X1)
  
  M <- diag(df^2)
  
  # List of gradient matrices for every spline coefficient
  deriv <- apply(M, 2, WoodTensor,
                 X1 = X1, X2 = X2,
                 simplify = FALSE)
  
  
  gradient <- gradientC(riskset = datalist$riskset,
                        logtheta = logtheta2,
                        df = df,
                        delta = datalist$delta.prod,
                        deriv = deriv,
                        I1 = datalist$I1,
                        I2 = datalist$I2,
                        I3 = datalist$I5,
                        I4 = datalist$I6) # gradientC returns vector of derivatives of -loglik
  
  if (!is.null(Sl)) {
    penalty <- t(coef.vector) %*% Sl
  } else penalty <- 0
  
  return(gradient + penalty)
}

# hessian <- function(coef.vector, degree, df, datalist, lambda) {
#   
#   logtheta2 <- tensor(datalist$X[, 1], datalist$X[, 2], degree = degree, coef.vector = coef.vector, df = df, knots = datalist$knots)
#   
#   M <- diag(df^2)
#   S1 <- Srow(df)
#   S2 <- Scol(df)
#   
#   # List of (df) matrices B_k (t1) * B_l (t2) which is the derivative of log(theta) wrt coefficient beta_jk. 
#   deriv <- apply(M, 2, tensor,
#     t1 = datalist$X[,1],
#     t2 = datalist$X[,2],
#     degree = degree,
#     df = df,
#     knots = datalist$knots,
#     simplify = FALSE
#   )
#   
#   index.offdiag <- combn(1:df^2, 2)
#   index.ondiag <- 1:df
#   
#   H.offdiag <- apply(index.offdiag, MARGIN = 2, function (x) {
#     hessianC(riskset = t(datalist$riskset),
#              logtheta = t(logtheta2),
#              deriv1 = deriv[[x[1]]],
#              deriv2 = deriv[[x[2]]],
#              delta = t(datalist$delta.prod),
#              I1 = datalist$I1,
#              I2 = datalist$I2) + 
#       hessianC(riskset = datalist$riskset,
#                logtheta = logtheta2,
#                deriv1 = deriv[[x[1]]],
#                deriv2 = deriv[[x[2]]],
#                delta = datalist$delta.prod,
#                I1 = t(datalist$I2),
#                I2 = t(datalist$I1))
#       
#   })
#   
#   H.ondiag <- apply(array(1:df^2), 1,  function (x) {
#     hessianC(riskset = t(datalist$riskset),
#              logtheta = t(logtheta2),
#              deriv1 = deriv[[x]],
#              deriv2 = deriv[[x]],
#              delta = t(datalist$delta.prod),
#              I1 = datalist$I1,
#              I2 = datalist$I2) + 
#       hessianC(riskset = datalist$riskset,
#                logtheta = logtheta2,
#                deriv1 = deriv[[x]],
#                deriv2 = deriv[[x]],
#                delta = datalist$delta.prod,
#                I1 = t(datalist$I2),
#                I2 = t(datalist$I1))
#     
#   })
#   
#   H <- matrix(0, ncol = df^2, nrow = df^2)
#   diag(H) <- H.ondiag
#   H[index.offdiag[1,], index.offdiag[2,]] <- H[index.offdiag[2,], index.offdiag[1,]] <- H.offdiag
# 
#   return(-H + lambda[1]*S1 + lambda[2]*S2)
# }


# loglik.poly <- function(coef.vector) {
#   
#   logtheta2 <- outer(X[,1], X[,2], function (x,y) polynomial(x,y, coef.vec = coef.vector))
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

SimData <- function (K, df, degree, unif.ub, alpha = 0.0023) {
  
  u1 <- runif(K, 0, 1)
  u2 <- runif(K, 0, 1)
  
  a <- alpha^u1 + (alpha - alpha^u1)*u2
  
  # Fan 2000
  T1 <- -log(u1)
  T2 <- -log(log(a/(a+(1-alpha)*u2),base = alpha))
  
  if (is.null(unif.ub)) {
    # Fan 2000
    X1 <- -log(u1)
    X2 <- -log(log(a/(a+(1-alpha)*u2),base = alpha))
    
    X <- as.matrix(cbind(X1,X2))
    
    delta1 <- delta2 <- rep(1, K)
    
  } else {
    # Hu 2011
    C1 <- runif(K, 0, unif.ub)
    C2 <- runif(K, 0, unif.ub)
    
    X1 <- pmin(T1,C1)
    X2 <- pmin(T2,C2)
    
    X <- as.matrix(cbind(X1,X2))
    
    delta1 <- 1*(T1 <= C1)
    delta2 <- 1*(T2 <= C2)
  }
  
  delta <- as.matrix(cbind(delta1,delta2))
  
  # # NOTE max+1 in geval van unif.ub = 5 geeft gradient=0 voor redelijk veel betas.
  # if (!is.null(unif.ub) && unif.ub < 5) {
  #   qq1 <- quantile(X1[delta1 == 1], probs = seq(0,1,length = df - degree + 2))
  #   knots1 <- c(min(X1)-1, qq1[-c(1,length(qq1))], max(X1)+1)
  #   qq2 <- quantile(X2[delta2 == 1], probs = seq(0,1,length = df - degree + 2))
  #   knots2 <- c(min(X2)-1, qq1[-c(1,length(qq2))], max(X2)+1)
  #   
  #   # knots1 <- seq(min(X1)-1, max(X1)+1, length.out = df - degree + 2)
  #   # knots2 <- seq(min(X2)-1, max(X2)+1, length.out = df - degree + 2)
  # } else {
  #   qq1 <- quantile(X1[delta1 == 1], probs = seq(0,1,length = df - degree + 2))
  #   knots1 <- c(min(X1)-1, qq1[-c(1,length(qq1))], max(X1))
  #   qq2 <- quantile(X2[delta2 == 1], probs = seq(0,1,length = df - degree + 2))
  #   knots2 <- c(min(X2)-1, qq1[-c(1,length(qq2))], max(X2))
  #   
  #   # knots1 <- seq(min(X1)-1, max(X1), length.out = df - degree + 2)
  #   # knots2 <- seq(min(X2)-1, max(X2), length.out = df - degree + 2)
  # }

  
  ## Check whether first delta1=delta2=1
  
  
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
  
  
  
  ## Calculating the risk set 
  
  # N <- outer(X[,1], X[,2], function(x,y) mapply(riskset,x,y))
  # N1 <- c(t(N))
  # N2 <- c(N)
  
  N <- risksetC(X[,1],X[,2])
  
  
  ## Calculating indicator functions in likelihood 
  
  #### I(X1j >= X1i)
  # I1 <- sapply(X[,1], function(x) 1*(X[,1] >= x)) # col=1,...,i,...,n row=1,...,j,...,n
  I1 <- IndGreater(X[,1])
  
  #### I(X2j <= X2i)
  # I2 <- sapply(X[,2], function(x) 1*(X[,2] <= x)) # col=1,...,i,...,n row=1,...,j,...,n
  I2 <- IndLess(X[,2])
  
  #### I(X2j >= X2i)
  # I3 <- sapply(X[,2], function(x) 1*(X[,2] >= x)) # col=1,...,i,...,n row=1,...,j,...,n
  # I3 <- t(I2)
  # 
  # #### I(X1j <= X1i)
  # # I4 <- sapply(X[,1], function(x) 1*(X[,1] <= x)) # col=1,...,i,...,n row=1,...,j,...,n
  # I4 <- t(I1)
  
  #### I(X1j = X1i) NOTE THAT THIS IS DIAG(1,500,500) IF NO TIES
  # I5 <- sapply(X[,1], function(x) 1*(X[,1] == x)) # col=1,...,i,...,n row=1,...,j,...,n
  I5 <- IndEqual(X[,1])
  
  #### I(X2j = X2i) NOTE THAT THIS IS DIAG(1,500,500) IF NO TIES
  # I6 <- sapply(X[,2], function(x) 1*(X[,2] == x)) # col=1,...,i,...,n row=1,...,j,...,n
  I6 <- IndEqual(X[,2])
  
  #I1 <- lapply(X1, function(x) 1*(X2 >= x))
  #test <- matrix(unlist(I1), ncol = 500, byrow = FALSE)
  
  
  # A1 <- c(I1*outer(delta[,2], delta[,1]))[N1 > 0]
  # A2 <- c(I3*outer(delta[,1], delta[,2]))[N2 > 0]
  
  delta.prod = DeltaC(delta[,1], delta[,2])

  return(list(X = cbind(X[,1],X[,2]),
              # knots = cbind(knots1,knots2),
              riskset = N,
              I1 = I1,
              I2 = I2,
              I5 = I5,
              I6 = I6,
              delta.prod = delta.prod))
}


# Srow <- function(df) {
#   P <- matrix(0, ncol = df^2, nrow = df^2)
#   index.col <- df - 2
#   for (i in 1:(index.col * df)) {
#     P[i, i] <- 1
#     P[i, i + df] <- -2
#     P[i, i + 2 * df] <- 1
#   }
# 
#   return(crossprod(P))
# 
# }

# Srow <- function(df, diff = 2) {
#   
#   P <- diff(diag(df^2), lag = df, differences = diff)
#   
#   return(crossprod(P))
#   
# }
# 
# Scol <- function(df, diff = 2) {
# 
#     # Block to fill 
#   P1 <- diff(diag(df), differences = diff)
#   
#   P <- kronecker(diag(df), P1)
#   
#   return(crossprod(P))
#   
# }

wrapper <- function(coef.vector, degree, datalist, Sl = NULL, H = NULL, minusLogLik=TRUE) { # H is hier gewoon de unpenalized hessian
  
  # Check whether penalty is applied
  if (is.null(Sl)) {
    penaltyLik <- logSl <- logdetH <- 0
  } else {

    # Calculate penalty terms for log f_lambda(y,beta) Wood (2017) p.1076 
    Sl.eigenv <- eigen(Sl, only.values = TRUE)$values
    Sl.eigenv[abs(Sl.eigenv) < 1e-6] <- 0
    
    penaltyLik <- t(coef.vector) %*% Sl %*% coef.vector
    logSl <- sum(log(Sl.eigenv[Sl.eigenv > 0]))
    # constant <- sum(Sl.eigenv == 0)*log(2*pi)/2 # Zie Wood (2016) p.1550
    
    # penaltyGrad <- t(t(coef.vector) %*% S.lambda)
    # penaltyHess <- S.lambda
  }
  
  # TODO implement sanity check that !is.null(S.lambda) & !is.null(H) & minusLogLik=FALSE implies laplace likelihood
  if (!is.null(H)) {
    logdetH <- log(det(H + Sl))
  } else logdetH <- 0
  
  logtheta2 <- tensor(datalist$X[,1], datalist$X[,2], degree = degree, coef.vector = coef.vector, df = sqrt(length(coef.vector)), knots = datalist$knots)
  
  
  # log f_lambda(y,beta)
  # ll <- logLikC(riskset = datalist$riskset,
  #               logtheta = logtheta2,
  #               delta = datalist$delta.prod,
  #               I1 = datalist$I1,
  #               I2 = datalist$I2,
  #               I3 = datalist$I5,
  #               I4 = datalist$I6) + penaltyLik - logS.lambda + logdetH - constant
  
  # ll <- logLikC(riskset = datalist$riskset,
  #               logtheta = logtheta2,
  #               delta = datalist$delta.prod,
  #               I1 = datalist$I1,
  #               I2 = datalist$I2,
  #               I3 = datalist$I5) + penaltyLik - logSl/2 + logdetH/2
  
  L1 <- logLikC(riskset = t(datalist$riskset),
                logtheta = t(logtheta2),
                delta = t(datalist$delta.prod),
                I1 = datalist$I1, I2 = datalist$I2, I3 = datalist$I5)
  L2 <- 0
  # L2 <- logLikC(riskset = datalist$riskset,
  #               logtheta = logtheta2,
  #               delta = datalist$delta.prod,
  #               I1 = t(datalist$I2), I2 = t(datalist$I1), I3 = datalist$I6)
  
  ll <- L1 + L2 + penaltyLik/2
  REML <- ll - logSl/2 + logdetH/2

  
  # # List of gradient matrices for every spline coefficient
  # M <- diag(df^2)
  # deriv <- apply(M, 2, tensor,
  #                t1 = datalist$X[,1], t2 = datalist$X[,2], degree = degree, df = df, knots = datalist$knots,
  #                simplify = FALSE)
  
  
  # gradient <- gradientC(riskset = datalist$riskset,
  #                      logtheta = logtheta2,
  #                      deriv = deriv,
  #                      df = df,
  #                      delta = datalist$delta.prod,
  #                      I1 = datalist$I1,
  #                      I2 = datalist$I2,
  #                      I3 = datalist$I5,
  #                      I4 = datalist$I6)
  # 
  # # TODO Controleer implementatie van hessian
  # 
  # hessian <- hessianC(riskset = datalist$riskset,
  #                     logtheta = logtheta2,
  #                     deriv = deriv,
  #                     df = df,
  #                     delta = datalist$delta.prod,
  #                     I1 = datalist$I1,
  #                     I2 = datalist$I2,
  #                     I3 = datalist$I5)
  # 
  # attr(ll, "gradient") <- gradient
  # attr(ll, "hessian") <- hessian
  
  # Merk op dat C++ code geïmplementeerd is voor -loglik
  sign <- ifelse(isTRUE(minusLogLik), 1, -1)
  
  return(list(ll = sign*ll, REML = sign*REML))
  
}

wrapper2 <- function(coef.vector, X1, X2, Sl = NULL, H = NULL, minusLogLik=TRUE) { # H is hier gewoon de unpenalized hessian
  
  # Check whether penalty is applied
  if (is.null(Sl)) {
    penaltyLik <- logSl <- logdetH <- 0
  } else {
    
    # Calculate penalty terms for log f_lambda(y,beta) Wood (2017) p.1076 
    Sl.eigenv <- eigen(Sl, only.values = TRUE)$values
    Sl.eigenv[abs(Sl.eigenv) < 1e-6] <- 0
    
    penaltyLik <- t(coef.vector) %*% Sl %*% coef.vector
    logSl <- sum(log(Sl.eigenv[Sl.eigenv > 0]))
    # constant <- sum(Sl.eigenv == 0)*log(2*pi)/2 # Zie Wood (2016) p.1550
    
    # penaltyGrad <- t(t(coef.vector) %*% S.lambda)
    # penaltyHess <- S.lambda
  }
  
  # TODO implement sanity check that !is.null(S.lambda) & !is.null(H) & minusLogLik=FALSE implies laplace likelihood
  if (!is.null(H)) {
    logdetH <- log(det(H + Sl))
  } else logdetH <- 0
  
  logtheta2 <- WoodTensor(X1 = X1, X2 = X2, coef.vector = coef.vector)
  
  
  # log f_lambda(y,beta)
  # ll <- logLikC(riskset = datalist$riskset,
  #               logtheta = logtheta2,
  #               delta = datalist$delta.prod,
  #               I1 = datalist$I1,
  #               I2 = datalist$I2,
  #               I3 = datalist$I5,
  #               I4 = datalist$I6) + penaltyLik - logS.lambda + logdetH - constant
  
  # ll <- logLikC(riskset = datalist$riskset,
  #               logtheta = logtheta2,
  #               delta = datalist$delta.prod,
  #               I1 = datalist$I1,
  #               I2 = datalist$I2,
  #               I3 = datalist$I5) + penaltyLik - logSl/2 + logdetH/2
  
  L1 <- logLikC(riskset = t(datalist$riskset),
                logtheta = t(logtheta2),
                delta = t(datalist$delta.prod),
                I1 = datalist$I1, I2 = datalist$I2, I3 = datalist$I5)
  L2 <- 0
  # L2 <- logLikC(riskset = datalist$riskset,
  #               logtheta = logtheta2,
  #               delta = datalist$delta.prod,
  #               I1 = t(datalist$I2), I2 = t(datalist$I1), I3 = datalist$I6)
  
  ll <- L1 + L2 + penaltyLik/2
  REML <- ll - logSl/2 + logdetH/2
  
  
  # # List of gradient matrices for every spline coefficient
  # M <- diag(df^2)
  # deriv <- apply(M, 2, tensor,
  #                t1 = datalist$X[,1], t2 = datalist$X[,2], degree = degree, df = df, knots = datalist$knots,
  #                simplify = FALSE)
  
  
  # gradient <- gradientC(riskset = datalist$riskset,
  #                      logtheta = logtheta2,
  #                      deriv = deriv,
  #                      df = df,
  #                      delta = datalist$delta.prod,
  #                      I1 = datalist$I1,
  #                      I2 = datalist$I2,
  #                      I3 = datalist$I5,
  #                      I4 = datalist$I6)
  # 
  # # TODO Controleer implementatie van hessian
  # 
  # hessian <- hessianC(riskset = datalist$riskset,
  #                     logtheta = logtheta2,
  #                     deriv = deriv,
  #                     df = df,
  #                     delta = datalist$delta.prod,
  #                     I1 = datalist$I1,
  #                     I2 = datalist$I2,
  #                     I3 = datalist$I5)
  # 
  # attr(ll, "gradient") <- gradient
  # attr(ll, "hessian") <- hessian
  
  # Merk op dat C++ code geïmplementeerd is voor -loglik
  sign <- ifelse(isTRUE(minusLogLik), 1, -1)
  
  return(list(ll = sign*ll, REML = sign*REML))
  
}

# EFS gebaseerd op de code van Simon Wood in het mgcv package (zie gam.fit4.r op github)
# gam.control() details in mgcv.r op github
EstimatePenal <- function(datalist, degree, S, lambda.init = c(10,10), tol = 0.001, eps = 1e-10, lambda.max = exp(15), step.control = TRUE) { 
  
  print("Extended Fellner-Schall method:")
  
  tiny <- .Machine$double.eps^0.5
  
  S1 <- S[[1]]
  S2 <- S[[2]]

  df <- sqrt(ncol(S1))
  # df <- sqrt(ncol(S))
  
  lambda.new <- lambda.init # In voorbeelden van Wood (2017) is de initiele lambda = 1
  
  fit <- efsud.fit(start = rep(1,df^2), degree = degree, datalist = datalist,
                   # Sl = lambda.init*S
                   Sl = lambda.init[1]*S1 + lambda.init[2]*S2
                   )
  k <- 1
  score <- rep(0, 200)
  for (iter in 1:200) {
    
    l0 <- fit$REML
    
    lambda <- lambda.new
    
    # Some calculations to update lambda later...
    Sl <- lambda[1]*S1 + lambda[2]*S2
    # Sl <- lambda*S
    Sl.inv <- MASS::ginv(Sl)

    # decomp <- eigen(hessian)
    # A <- diag(abs(decomp$values))
    # hessian <- decomp$vectors %*% A %*% t(decomp$vectors)
    
    
    # Update ----
    
    # Calculate V
    V <- solve(fit$hessian + Sl)
    
    # Calculate trSSj, trVS and bSb
    trSSj <- trVS <- bSb <- rep(0, length(S))
    for (i in 1:length(S)) {
      trSSj[i] <- sum(diag(Sl.inv %*% S[[i]]))
      trVS[i] <- sum(diag(V %*% S[[i]]))
      bSb[i] <- t(fit$beta) %*% S[[i]] %*% fit$beta
    }
    
    # trSSj <- sum(diag(Sl.inv %*% S))
    # trVS <- sum(diag(V %*% S))
    # bSb <- t(fit$beta) %*% S %*% fit$beta
    
    # Update lambdas
    a <- pmax(tiny, trSSj - trVS)
    update <- a/pmax(tiny, bSb)
    update[a==0 & bSb==0] <- 1
    update[!is.finite(update)] <- 1e6
    lambda.new <- pmin(update*lambda, lambda.max)
    
    # Step length of update
    max.step <- max(abs(lambda.new - lambda))
    
    # Create new S.lambda matrix
    Sl.new <- lambda.new[1]*S1 + lambda.new[2]*S2
    # Sl.new <- lambda.new*S
    
    fit <- efsud.fit(start = fit$beta, degree = degree, datalist = datalist, Sl = Sl.new)
    l1 <- fit$REML
    
    # Start of step control ----
  if (step.control) {
    
    if (l1 > l0) { # Improvement
      if(max.step < 1) { # Consider step extension
        lambda2 <- pmin(lambda*update^(k*2), exp(12))
        fit2 <- efsud.fit(start = fit$beta, degree = degree, datalist = datalist,
                          # Sl = lambda2*S
                          Sl = lambda2[1]*S1 + lambda2[2]*S2
                          )
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
          fit <- efsud.fit(start = fit$beta, degree = degree, datalist = datalist,
                           # Sl = lambda3*S
                           Sl = lambda3[1]*S1 + lambda3[2]*S2
                           )
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
    if (iter > 3 && max.step < 0.5 && max(abs(diff(score[(iter-3):iter]))) < 0.01) {print("REML not changing"); break}
    # Or break is likelihood does not change
    if (l1 == l0) {print("Loglik not changing"); break}
    # Stop if loglik is not changing
    if (iter==1) old.ll <- fit$ll else {
      if (abs(old.ll-fit$ll)<eps*abs(fit$ll)) {print("Loglik not changing"); break}  # *100
      old.ll <- fit$ll
    }
    
    # Print information while running...
    print(paste0("Iteration ", iter,
                 ": k = ", k,
                 # " lambda = ", round(lambda.new,4),
                 " lambda1 = ", round(lambda.new[1],4),
                 " lambda2 = ", round(lambda.new[2],4),
                 " ll = ", round(fit$ll,4),
                 " REML = ", score[iter]))
    
  } # End of for loop
  
  if (iter < 200) print("Converged") else print("Number of iterations is too small")
  
  
  return(list(
    beta = fit$beta,
    lambda = lambda.new,
    iterations = iter,
    ll = fit$ll,
    history = score[1:iter]))
}

EstimatePenal2 <- function(datalist, dim, degree = 3, lambda.init = c(10,10), type = "bs", quantile = FALSE, scale = TRUE, tol = 0.001, eps = 1e-10, lambda.max = exp(15), step.control = TRUE) { 
  
  print("Extended Fellner-Schall method:")
  
  tiny <- .Machine$double.eps^0.5
  
  obj1 <- WoodSpline(t = datalist$X[,1], dim = dim, degree = 3, scale = scale, quantile = quantile, type = type)
  obj2 <- WoodSpline(t = datalist$X[,2], dim = dim, degree = 3, scale = scale, quantile = quantile, type = type)
  
  X1 <- obj1$X
  X2 <- obj2$X
  
  S <- WoodPenalty(obj1,obj2)
  S1 <- S[[1]]
  S2 <- S[[2]]
  
  lambda.new <- lambda.init # In voorbeelden van Wood (2017) is de initiele lambda = 1
  
  fit <- efsud.fit2(start = rep(1,df^2), X1 = X1, X2 = X2, datalist = datalist,
                   # Sl = lambda.init*S
                   Sl = lambda.init[1]*S1 + lambda.init[2]*S2
  )
  k <- 1
  score <- rep(0, 200)
  for (iter in 1:200) {
    
    l0 <- fit$REML
    
    lambda <- lambda.new
    
    # Some calculations to update lambda later...
    Sl <- lambda[1]*S1 + lambda[2]*S2
    # Sl <- lambda*S
    Sl.inv <- MASS::ginv(Sl)
    
    # decomp <- eigen(hessian)
    # A <- diag(abs(decomp$values))
    # hessian <- decomp$vectors %*% A %*% t(decomp$vectors)
    
    
    # Update ----
    
    # Calculate V
    V <- solve(fit$hessian + Sl)
    
    # Calculate trSSj, trVS and bSb
    trSSj <- trVS <- bSb <- rep(0, length(S))
    for (i in 1:length(S)) {
      trSSj[i] <- sum(diag(Sl.inv %*% S[[i]]))
      trVS[i] <- sum(diag(V %*% S[[i]]))
      bSb[i] <- t(fit$beta) %*% S[[i]] %*% fit$beta
    }
    
    # trSSj <- sum(diag(Sl.inv %*% S))
    # trVS <- sum(diag(V %*% S))
    # bSb <- t(fit$beta) %*% S %*% fit$beta
    
    # Update lambdas
    a <- pmax(tiny, trSSj - trVS)
    update <- a/pmax(tiny, bSb)
    update[a==0 & bSb==0] <- 1
    update[!is.finite(update)] <- 1e6
    lambda.new <- pmin(update*lambda, lambda.max)
    
    # Step length of update
    max.step <- max(abs(lambda.new - lambda))
    
    # Create new S.lambda matrix
    Sl.new <- lambda.new[1]*S1 + lambda.new[2]*S2
    # Sl.new <- lambda.new*S
    
    fit <- efsud.fit2(start = fit$beta, X1 = X1, X2 = X2, datalist = datalist, Sl = Sl.new)
    l1 <- fit$REML
    
    # Start of step control ----
    if (step.control) {
      
      if (l1 > l0) { # Improvement
        if(max.step < 1) { # Consider step extension
          lambda2 <- pmin(lambda*update^(k*2), exp(12))
          fit2 <- efsud.fit2(start = fit$beta, X1 = X1, X2 = X2, datalist = datalist,
                            # Sl = lambda2*S
                            Sl = lambda2[1]*S1 + lambda2[2]*S2
          )
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
          fit <- efsud.fit2(start = fit$beta, X1 = X1, X2 = X2, datalist = datalist,
                           # Sl = lambda3*S
                           Sl = lambda3[1]*S1 + lambda3[2]*S2
          )
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
    if (iter > 3 && max.step < 0.5 && max(abs(diff(score[(iter-3):iter]))) < 0.01) {print("REML not changing"); break}
    # Or break is likelihood does not change
    if (l1 == l0) {print("Loglik not changing"); break}
    # Stop if loglik is not changing
    if (iter==1) old.ll <- fit$ll else {
      if (abs(old.ll-fit$ll)<eps*abs(fit$ll)) {print("Loglik not changing"); break}  # *100
      old.ll <- fit$ll
    }
    
    # Print information while running...
    print(paste0("Iteration ", iter,
                 ": k = ", k,
                 # " lambda = ", round(lambda.new,4),
                 " lambda1 = ", round(lambda.new[1],4),
                 " lambda2 = ", round(lambda.new[2],4),
                 " ll = ", round(fit$ll,4),
                 " REML = ", score[iter]))
    
  } # End of for loop
  
  if (iter < 200) print("Converged") else print("Number of iterations is too small")
  
  return(list(
    beta = fit$beta,
    lambda = lambda.new,
    iterations = iter,
    ll = fit$ll,
    history = score[1:iter],
    knots = list(knots1 = obj1$knots, knots2 = obj2$knots)))
}

efsud.fit <- function(start, degree, datalist, Sl) {
  beta <- multiroot(Score, start = start, rtol = 1e-10, degree = degree, datalist = datalist, Sl = Sl)$root
  H <- derivatives(coef.vector = beta, degree = degree, datalist = datalist, gradient = FALSE, hessian = TRUE)$hessian
  fit <-  wrapper(coef.vector = beta,
                   degree = degree,
                   Sl = Sl, H = H,
                   minusLogLik = FALSE,
                   datalist = datalist)
  
  return(list(beta = beta, hessian = H, REML = fit$REML, ll = fit$ll))
}

efsud.fit2 <- function(start, X1, X2, datalist, Sl) {
  beta <- multiroot(Score2, start = start, rtol = 1e-10, X1 = X1, X2 = X2, Sl = Sl)$root
  H <- derivatives2(coef.vector = beta, X1 = X1, X2 = X2, datalist = datalist, gradient = FALSE, hessian = TRUE)$hessian
  fit <-  wrapper2(coef.vector = beta,
                  X1 = X1, X2 = X2,
                  Sl = Sl, H = H,
                  minusLogLik = FALSE,
                  datalist = datalist)
  
  return(list(beta = beta, hessian = H, REML = fit$REML, ll = fit$ll))
}
