
source('hunan-functions.R')

df = 6
lambda <- c(12.0000,12.0000)
S <- list(Srow(df), Scol(df))

test <- eigen(lambda[1]*S[[1]] + lambda[2]*S[[2]], only.values = TRUE)
# test$values[abs(test$values) < 1e-6] <- 0
test1 <- sum(log(test$values[test$values > 0]))
test1

ldetS <- function(S, lambda) {
  
  Sl <- vector("list", length = length(lambda))
  
  for (i in 1:length(lambda)) {
    Sl[[i]] <- lambda[i]*S[[i]]
    assign(paste0("S",i), S[[i]])
    assign(paste0("Sl",i), lambda[i]*S[[i]])
    
  }
  
  ind <- ifelse(norm(Sl1) > norm(Sl2), 1, 2)

  eig <- eigen(Sl[[ind]])
  # eig$values[abs(eig$values) < 1e-6] <- 0
  r <- eig$values > 0
  D <- diag(eig$values[r])
  Ur <- eig$vector[,r]
  Un <- eig$vector[,!r]
  
  M1 <- D + t(Ur) %*% Sl[[-ind]] %*% Ur
  M2 <- t(Ur) %*% Sl[[-ind]] %*% Un
  M3 <- t(Un) %*% Sl[[-ind]] %*% Ur
  M4 <- t(Un) %*% Sl[[-ind]] %*% Un
  
  Slp <- rbind(cbind(M1,M2), cbind(M3,M4))
  # Slp[abs(Slp) < 1e-6] <- 0
  
  eigp <- eigen(Slp)
  # eigp$values[abs(eigp$values) < 1e-6] <- 0
  
  sum(log(eigp$values[eigp$values > 0]))
    
}
