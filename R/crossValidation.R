computeErrorMatrixR <- function(Zrand,Erand,Yrand,Hrand,Q,lambda,K,cov.missing)
{

  n <- nrow(Zrand)
  p <- ncol(Zrand)

  n.lambda <- length(lambda)
  errors <- matrix(nrow = K, ncol = n.lambda)

  for(j in 1:n.lambda){
    for(k in 1:K){
      testIdx <- ((n / K * (k - 1) + 1):(n / K * k))
      trainIdx <- sort(setdiff(1:n, testIdx), decreasing=TRUE)
      Ytrain <- Yrand[trainIdx]
      Ytest <- Yrand[testIdx]
      Ztrain <- Zrand[trainIdx, ]
      Ztest <- Zrand[testIdx, ]
      Etrain <- Erand[trainIdx, ]
      Etest <- Erand[testIdx, ]
      Htrain <- Hrand[trainIdx, trainIdx]
      Htest <- Hrand[  testIdx,  testIdx]

      n.train <- nrow(Ztrain)

      if (cov.missing) P <- diag(n.train)
      else P <- diag(n.train) - (Etrain %*% solve( t(Etrain) %*% Htrain %*% Etrain ) %*% t(Etrain) %*% Htrain)
      Y.p <- P %*% Ytrain
      Z.p <- P %*% Ztrain

      beta.hat <- Q %*% solve( t(Z.p) %*% Htrain %*% Z.p %*% Q + lambda[j] * diag(p) ) %*% t(Z.p) %*% Htrain %*% Y.p # penalized coefficients
      if (cov.missing) eta.hat <- 0
      else eta.hat <- solve(t(Etrain) %*% Htrain %*% Etrain) %*% t(Etrain) %*% Htrain %*% (Ytrain - Ztrain %*% beta.hat) # unpenalized coefficients

      coefficients <- c(beta.hat, eta.hat)
      if(length(testIdx) == 1)
      {
        Ztest <- t(matrix(Ztest))
        Etest <- t(matrix(Etest))
      }
      Xtest <- cbind(Ztest, Etest)

      yhat <- Xtest %*% coefficients
      errors[k, j] <- t(Ytest - yhat) %*% Htest %*%  (Ytest - yhat)
    }
  }

  return(errors)

}
