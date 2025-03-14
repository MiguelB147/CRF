theta.frank <- function(x,y,alpha) {
  A <- (alpha-1)*log(alpha)*alpha^(2-exp(-x)-exp(-y))
  B <- (alpha^(1-exp(-x))-alpha)*(alpha^(1-exp(-y))-alpha)
  C <- -1 + exp(-x) + exp(-y) + log(1+ (alpha^(1-exp(-x))-1)*(alpha^(1-exp(-y))-1)/(alpha-1), base = alpha)
  return(A*C/B)
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



# eval_function <- function(coef.vector, degree, df, datalist) {
# 
#   logtheta2 <- tensor(datalist$X[,1], datalist$X[,2], degree = degree, coef.vector = coef.vector, df = df, knots = datalist$knots)
#   logtheta.deriv <- tensor.deriv(datalist$X[,1], datalist$X[,2], degree = degree, df = df, knots = datalist$knots)
#     
#   U1 <- sum(diag(datalist$delta.prod)*diag(logtheta.deriv))
#   U2 <- ScoreFunc(riskset = t(datalist$riskset),
#                   logthetaderiv = t(logtheta.deriv),
#                   logtheta = t(logtheta2),
#                   delta = t(datalist$delta.prod),
#                   I1 = datalist$I1,
#                   I2 = datalist$I2)
#   U4 <- ScoreFunc(riskset = datalist$riskset,
#                   logthetaderiv = logtheta.deriv,
#                   logtheta = logtheta2,
#                   delta = datalist$delta.prod,
#                   I1 = t(datalist$I1),
#                   I2 = t(datalist$I2))
#   
#   return((2*U1 - U2 - U4)/nrow(datalist$X))
#   
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
  
  # set.seed(123)
  

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
  
  # NOTE max+1 in geval van unif.ub = 5 geeft gradient=0 voor redelijk veel betas.
  if (!is.null(unif.ub) && unif.ub < 5) {
    qq1 <- quantile(X1[delta1 == 1], probs = seq(0,1,length = df - degree + 2))
    knots1 <- c(min(X1)-1, qq1[-c(1,length(qq1))], max(X1)+1)
    qq2 <- quantile(X2[delta2 == 1], probs = seq(0,1,length = df - degree + 2))
    knots2 <- c(min(X2)-1, qq1[-c(1,length(qq2))], max(X2)+1)
    
    # knots1 <- seq(min(X1)-1, max(X1)+1, length.out = df - degree + 2)
    # knots2 <- seq(min(X2)-1, max(X2)+1, length.out = df - degree + 2)
  } else {
    qq1 <- quantile(X1[delta1 == 1], probs = seq(0,1,length = df - degree + 2))
    knots1 <- c(min(X1)-1, qq1[-c(1,length(qq1))], max(X1))
    qq2 <- quantile(X2[delta2 == 1], probs = seq(0,1,length = df - degree + 2))
    knots2 <- c(min(X2)-1, qq1[-c(1,length(qq2))], max(X2))
    
    # knots1 <- seq(min(X1)-1, max(X1), length.out = df - degree + 2)
    # knots2 <- seq(min(X2)-1, max(X2), length.out = df - degree + 2)
  }

  
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
              knots = cbind(knots1,knots2),
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

Srow <- function(df) {
  
  P <- diff(diag(df^2), lag = df, differences = 2)
  
  return(crossprod(P))
  
}

Scol <- function(df) {
  P <- matrix(0, ncol = df ^ 2, nrow = df ^ 2)
  
  # Block to fill 
  P1 <- diff(diag(df), differences = 2)
  
  P <- kronecker(diag(df), P1)
  
  return(crossprod(P))
  
}

# Srow <- function(df) {
#   P <- matrix(0, ncol = df^2, nrow = df^2)
#   index.col <- df - 2
#   for (j in 0:(index.col-1)) {
#     P[(1+df*j):(df+j*df),(1+j*df):(df*j+df)] <- diag(1, ncol = df, nrow = df)
#     P[(1+df*j):(df+j*df),(j*df+df+1):(2*df+j*df)] <- diag(-2, ncol = df, nrow = df)
#     P[(1+df*j):(df+j*df),(2*df+j*df+1):(3*df+j*df)] <- diag(1, ncol = df, nrow = df)
#     }
#   
#   return(crossprod(P))
#   
# }

# Internal function used in EstimatePenalty that provides general Fellner-Schall update
lambdaUpdate <- function(lambda.old, S.lambda.inv, S, V, beta) {
  
  S1 <- S[[1]]
  S2 <- S[[2]]
  
  tr1 <- sum(diag(S.lambda.inv %*% S1))
  tr2 <- sum(diag(S.lambda.inv %*% S2))
  trV1 <- sum(diag(V %*% S1))
  trV2 <- sum(diag(V %*% S2))
  denom1 <- t(beta) %*% S1 %*% beta
  denom2 <- t(beta) %*% S2 %*% beta
  
  update1 <- (tr1 - trV1)/denom1
  update2 <- (tr2 - trV2)/denom2
  
  return(c(update1,update2)*lambda.old)
  
}

wrapper <- function(coef.vector, degree, datalist, Sl = NULL, H = NULL, minusLogLik=TRUE) { # H is hier gewoon de unpenalized hessian
  
  # Check whether penalty is applied
  if (is.null(Sl)) {
    penaltyLik <- logSl <- logdetH <- 0
  } else {

    # Calculate penalty terms for log f_lambda(y,beta) Wood (2017) p.1076 
    Sl.eigenv <- eigen(Sl)$values
    
    penaltyLik <- t(coef.vector) %*% Sl %*% coef.vector
    logSl <- log(prod(Sl.eigenv[Sl.eigenv > 0]))
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
  
  ll <- L1 + L2 + penaltyLik/2 - logSl/2 + logdetH/2
  
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
  
  return(sign*ll)
  
}

# FIXME Waarom werkt minusLogLik niet?? Final loglik in estimatepenaltest geeft negatieve waarde
# FIXME penaltyLik = 0 !!!!!! Probleem met S
wrapperTest <- function(coef.vector, degree, datalist, lambda, S=NULL, H = NULL, minusLogLik=TRUE) {
  
  # Check whether penalty is applied
  if (is.null(S)) {
    penaltyLik <- penaltyGrad <- penaltyHess <- 0
  } else {
    
    Sl <- lambda[1]*S[[1]] + lambda[2]*S[[2]]
    
    # Calculate penalty terms for log f_lambda(y,beta) Wood (2017) p.1076 
    Sl.eigenv <- eigen(Sl)$values
    
    # penaltyLik <- (t(coef.vector) %*% Sl %*% coef.vector)/2
    logSl <- log(prod(Sl.eigenv[Sl.eigenv > 0]))
    constant <- sum(Sl.eigenv == 0)*log(2*pi)/2 # Zie Wood (2016) p.1550
    
    # penaltyGrad <- t(t(coef.vector) %*% Sl)
    # penaltyHess <- Sl
  }
  
  # TODO implement sanity check that !is.null(S.lambda) & !is.null(H) & minusLogLik=FALSE implies laplace likelihood
  if (is.null(H)) {
    logdetH <- 0
  } else {
    logdetH <- log(det(H))
    logSl <- logSl/2
  }
  
  df <- sqrt(length(coef.vector))
  
  # Merk op dat C++ code geïmplementeerd is voor -loglik
  sign <- ifelse(isTRUE(minusLogLik), 1, -1)
  
  # log f_lambda(y,beta)
  ll <- loglikPenal(coef.vector, degree, df, datalist, S = Sl) - logSl + logdetH - constant
  
  return(sign*ll)
  
}

# LaplaceLogLik <- function(coef.vector, degree, S.lambda, H, datalist) {
#   
#   logdetH <- log(det(H))
#   
#   S.lambda.eigenv <- eigen(S.lambda)$values
#   logS.lambda <- log(prod(S.lambda.eigenv[S.lambda.eigenv > 0]))
#   
#   ll <- wrapper(
#     coef.vector = coef.vector,
#     degree = degree,
#     S.lambda = S.lambda,
#     datalist = datalist,
#     minusLogLik = FALSE # Geeft loglik ipv -loglik
#   )
#   
#   laplace <- ll - logdetH/2
#   
#   return(laplace)
#   
# }

# TODO Waarom wordt er geen step length control gebruikt in de voorbeelden van Wood?
EstimatePenalty <- function(datalist, degree, S, lambda.init = c(1,1), tol = 0.01, maxiter=50) {
  
  S1 <- S[[1]]
  S2 <- S[[2]]
  
  df <- sqrt(ncol(S1))
  
  lambda.new <- lambda.init # In voorbeelden van Wood (2017) is de initiele lambda = 1
  
  lldiff <- 1e10
  beta <- rep(1,df^2)
  
  iter <- 0
  
  if (iter == 0) {print("Algorithm running...")}
  
  while (lldiff > tol & iter <= maxiter) {
    
    # Update number of iterations
    iter = iter + 1
    
    lambda <- lambda.new
    
    # Some calculations to update lambda later...
    S.lambda <- lambda[1]*S1 + lambda[2]*S2
    S.lambda.inv <- ginv(S.lambda)
    
    # Estimate betas for given lambdas
    beta.fit <- nlm(f = wrapper,
                    p = beta,
                    degree = degree,
                    S.lambda = S.lambda,
                    datalist = datalist,
                    hessian = TRUE)
    
    # New betas to be used as initial values for possible next iteration
    beta <- beta.fit$estimate
    
    # Make sure that observed hessian is positive definite
    decomp <- eigen(beta.fit$hessian)
    A <- diag(abs(decomp$values))
    hessian.obs <- decomp$vectors %*% A %*% t(decomp$vectors)
    
    V <- solve(hessian.obs + S.lambda)
    
    # Update lambdas
    lambda.new <- lambdaUpdate(lambda, S.lambda.inv, S, V, beta)
    
    # Create new S.lambda matrix
    S.lambda.new <- lambda.new[1]*S1 + lambda.new[2]*S2
    
    # lambda <- c((sum(diag(S.lambda.inv %*% S1)) - sum(diag(V %*% S1))) /
    #               (t(beta.new) %*% S1 %*% beta.new),
    #             (sum(diag(S.lambda.inv %*% S2)) - sum(diag(V %*% S2))) /
    #               (t(beta.new) %*% S2 %*% beta.new))
   
    # Step length of update
    diff <- lambda.new - lambda
    
    # Assess whether update is an increase in the log-likelihood
    # If not, apply step length control
    l1 <- wrapper(
      coef.vector = beta,
      degree = degree,
      S.lambda = S.lambda.new,
      H = hessian.obs + S.lambda.new,
      minusLogLik = FALSE,
      datalist = datalist
    )
    
    l0 <- wrapper(
      coef.vector = beta,
      degree = degree,
      S.lambda = S.lambda,
      H = hessian.obs + S.lambda,
      minusLogLik = FALSE,
      datalist = datalist
    )
  
    # Step length control to guarantee increase in loglik
    k = 0

    while (l1 < l0) { 
      
      k = k + 1
      
      delta <- diff/(2^k)
      
      S.lambda.delta <- (lambda + delta)[1]*S1 + (lambda + delta)[2]*S2
      
      l1 <- wrapper(
        coef.vector = beta,
        degree = degree,
        S.lambda = S.lambda.delta,
        H = hessian.obs + S.lambda.delta,
        minusLogLik = FALSE,
        datalist = datalist
      )
      
      l0 <- wrapper(
        coef.vector = beta,
        degree = degree,
        S.lambda = S.lambda,
        H = hessian.obs + S.lambda,
        minusLogLik = FALSE,
        datalist = datalist
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
    
    # Print information while running...
    print(paste0("Iteration ", iter,
                 ": k = ", k,
                 " lambda1 = ", lambda.new[1],
                 " lambda2 = ", lambda.new[2],
                 " Likelihood increase = ", lldiff))
    
  } # end of outer while loop
  
  if (iter == maxiter) {message <- "Maximum iterations reached"} else {message <- "Convergence reached"}
  

  return(list(
    beta = beta,
    lambda = lambda.new,
    iterations = iter,
    status = message))
}

EstimatePenaltyNoControl <- function(datalist, degree, S, lambda.init = c(1,1), tol = 0.01, maxiter=50) {
  
  S1 <- S[[1]]
  S2 <- S[[2]]
  
  df <- sqrt(ncol(S1))
  
  lambda.new <- lambda.init # In voorbeelden van Wood (2017) is de initiele lambda = 1
  
  lldiff <- 1e10
  beta <- rep(1,df^2)
  
  iter <- 0
  
  if (iter == 0) {print("Algorithm running...")}
  
  while (lldiff > tol & iter <= maxiter) {
    
    # Update number of iterations
    iter = iter + 1
    
    lambda <- lambda.new
    
    # Some calculations to update lambda later...
    S.lambda <- lambda[1]*S1 + lambda[2]*S2
    S.lambda.inv <- ginv(S.lambda)
    
    # Estimate betas for given lambdas
    beta.fit <- nlm(f = wrapper,
                    p = beta,
                    degree = degree,
                    S.lambda = S.lambda,
                    datalist = datalist,
                    hessian = TRUE)
    
    # New betas to be used as initial values for possible next iteration
    beta <- beta.fit$estimate
    
    # Make sure that observed hessian is positive definite
    decomp <- eigen(beta.fit$hessian)
    A <- diag(abs(decomp$values))
    hessian.obs <- decomp$vectors %*% A %*% t(decomp$vectors)
    
    V <- solve(hessian.obs + S.lambda)
    
    # Update lambdas
    lambda.new <- lambdaUpdate(lambda, S.lambda.inv, S, V, beta)
    
    # Create new S.lambda matrix
    S.lambda.new <- lambda.new[1]*S1 + lambda.new[2]*S2
    
    # lambda <- c((sum(diag(S.lambda.inv %*% S1)) - sum(diag(V %*% S1))) /
    #               (t(beta.new) %*% S1 %*% beta.new),
    #             (sum(diag(S.lambda.inv %*% S2)) - sum(diag(V %*% S2))) /
    #               (t(beta.new) %*% S2 %*% beta.new))
    
    # Step length of update
    diff <- lambda.new - lambda
    
    # Assess whether update is an increase in the log-likelihood
    # If not, apply step length control
    l1 <- wrapper(
      coef.vector = beta,
      degree = degree,
      S.lambda = S.lambda.new,
      H = hessian.obs + S.lambda.new,
      minusLogLik = FALSE,
      datalist = datalist
    )
    
    l0 <- wrapper(
      coef.vector = beta,
      degree = degree,
      S.lambda = S.lambda,
      H = hessian.obs + S.lambda,
      minusLogLik = FALSE,
      datalist = datalist
    )
    
    # Final difference in loglikelihood
    lldiff <- l1 - l0
    
    if(lldiff < 0) {message("The likelihood decreased for new lambda"); break}

    # Sanity check: lambda must be positive
    if (sum(lambda.new < 0) > 0) {stop("At least 1 lambda is negative")}
    
    # Print information while running...
    print(paste0("Iteration ", iter,
                 " lambda1 = ", lambda.new[1],
                 " lambda2 = ", lambda.new[2],
                 " Likelihood increase = ", lldiff))
    
  } # end of outer while loop
  
  if (iter == maxiter) {message <- "Maximum iterations reached"} else {message <- "Convergence reached"}
  
  
  return(list(
    beta = beta,
    lambda = lambda.new,
    iterations = iter,
    status = message))
}

# EFS gebaseerd op de code van Simon Wood in het mgcv package (zie gam.fit4.r op github)
# gam.control() details in mgcv.r op github
EstimatePenalAsym <- function(datalist, degree, S, lambda.init = c(1,1), tol = 0.001, eps = 1e-7, lambda.max = exp(15), step.control = TRUE) { 
  
  print("Extended Fellner-Schall method:")
  
  tiny <- .Machine$double.eps^0.5
  
  S1 <- S[[1]]
  S2 <- S[[2]]

  df <- sqrt(ncol(S1))
  # df <- sqrt(ncol(S))
  
  lambda.new <- lambda.init # In voorbeelden van Wood (2017) is de initiele lambda = 1
  
  fit <- efsud.fit(start = rep(1,df^2), degree = degree, datalist = datalist, Sl = lambda.init[1]*S1 + lambda.init[2]*S2)

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
    
    # Calculate V
    V <- solve(fit$hessian + Sl)
    
    # Calculate trSSj, trVS and bSb
    trSSj <- trVS <- bSb <- rep(NA, length(S))
    for (i in length(S)) {
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
    
  if (step.control) {
    
    if (l1 > l0) { # Improvement
      if(max.step < 1) { # Consider step extension
        lambda2 <- pmin(lambda*update^(k*2), exp(12))
        fit2 <- efsud.fit(start = fit$beta, degree = degree, datalist = datalist, Sl = lambda2[1]*S1 + lambda2[2]*S2)
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
          fit <- efsud.fit(start = fit$beta, degree = degree, datalist = datalist, Sl = lambda3[1]*S1 + lambda3[2]*S2)
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
    
    # Break procedure if REML change and step size are too small
    if (iter > 3 && max.step < 1 && max(abs(diff(score[(iter-3):iter]))) < 1) break # change lik tol to 0.1
    # Or break is likelihood does not change
    if (l1 == l0) break
    # if (iter==1) l0 <- l1 else {
    #   if (abs(l1-l0)<eps*abs(l1)) break
    #   l0 <- l1
    # }
    
    # Print information while running...
    print(paste0("Iteration ", iter,
                 ": k = ", k,
                 # " lambda = ", round(lambda.new,4),
                 " lambda1 = ", round(lambda.new[1],4),
                 " lambda2 = ", round(lambda.new[2],4),
                 " Score increase = ", round(score[iter] - score[iter-1],2),
                 " REML = ", score[iter]))
    
  } # End of for loop
  
  if (iter < 200) print("Converged") else print("Number of iterations is too small")
  
  
  return(list(
    beta = fit$beta,
    lambda = lambda.new,
    iterations = iter,
    history = score[1:iter]))
}

efsud.fit <- function(start, degree, datalist, Sl) {
  beta <- multiroot(Score, start = start, rtol = 1e-10, degree = degree, datalist = datalist, Sl = Sl)$root
  H <- derivatives(coef.vector = beta, degree = degree, datalist = datalist, gradient = FALSE, hessian = TRUE)$hessian
  REML <-  wrapper(coef.vector = beta,
                   degree = degree,
                   Sl = Sl, H = H,
                   minusLogLik = FALSE,
                   datalist = datalist)
  return(list(beta = beta, hessian = H, REML = REML))
}
