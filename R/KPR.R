#' Kernel Penalized Regression
#'
#' Fits a kernel penalized regression model using a design matrix X, response vector Y, sample similarity kernel H, and variable similarity kernel Q.
#' The lambda parameter is found through k-folds cross validation.
#'
#' @param X An n x p data matrix. Should be scaled and centered.
#' @param Y An n x 1 response vector. Should be scaled and centered.
#' @param H An n x n sample similarity kernel. Must be symmetric and positive semidefinite. This defaults to an identity matrix.
#' @param Q A p x p variable similarity kernel. Must be symmetric and postive semidefinite. This defaults to an identity matrix.
#' @param n.lambda The number of lambda values to test through the cross validation search. The values are generated internally.
#' @param lambda A vector of lambda values to test through cross validation. This will override the sequence generated with the n.lambda parameter.
#' @param K The number of folds in the cross validation search.
#' @return
#' \item{betahat.kpr}{A vector of estimated coefficients}
#' \item{lambda.min}{The value of lambda that resulted in the minimun squared error. This value is used to generated the vector of estimated coefficients.}
#' \item{lambda.1se}{The largest value of lambda within one standard error of lambda.min.}
#' @importFrom stats sd
#' @export
KPR <- function(X, Y, H = diag(nrow(X)), Q = diag(ncol(X)),
                n.lambda = 200, lambda, K = 5)
{
  n <- nrow(X)
  p <- ncol(X)
  if(missing(lambda))
    lambda <- exp(seq(from = 0, to = 10, length.out = n.lambda))
  n.lambda <- length(lambda)

  ERR.QH.KPR <- matrix(nrow = K, ncol = n.lambda)
  randidx <- sample(1:n, n)
  Yrand <- Y[randidx]
  Xrand <- X[randidx, ]
  Hrand <- H[randidx, randidx]
  for(j in 1:n.lambda){
    for(k in 1:K){
      Ytrain <- Yrand[-((n / K * (k - 1) + 1):(n / K * k))]
      Ytest <- Yrand[((n / K * (k - 1) + 1):(n / K * k))]
      Xtrain <- Xrand[-((n / K * (k - 1) + 1):(n / K * k)), ]
      Xtest <- Xrand[((n / K * (k - 1) + 1):(n / K * k)), ]
      Htrain <- Hrand[-((n / K * (k - 1) + 1):(n / K * k)), -((n / K * (k - 1) + 1):(n / K * k))]
      Htest <- Hrand[  ((n / K * (k - 1) + 1):(n / K * k)),  ((n / K * (k - 1) + 1):(n / K * k))] #needed?

      beta.QH.KPR <- Q %*% solve( t(Xtrain) %*% Htrain %*% Xtrain %*% Q  + lambda[j] * diag(p)) %*% t(Xtrain) %*% Htrain %*% Ytrain

      y.QH.KPR <- Xtest %*% beta.QH.KPR
      ERR.QH.KPR[k, j] <- t(Ytest - y.QH.KPR) %*% Htest %*%  (Ytest - y.QH.KPR) # sum((Ytest - y.QH.KPR)^2)
    }
  }
  CVidx.QH.KPR <- which.min(colSums(ERR.QH.KPR))

  se.QH.KPR <- sd(ERR.QH.KPR[, CVidx.QH.KPR]) * sqrt(K)
  CVidx.QH.KPR1se <- which(colSums(ERR.QH.KPR) < min(colSums(ERR.QH.KPR)) + se.QH.KPR)[1]

  beta.kpr <- beta.QH.KPR <- Q %*% solve( t(X) %*% H %*% X %*% Q  + lambda[CVidx.QH.KPR] * diag(p)) %*% t(X) %*% H %*% Y
  beta.kpr <- as.vector(beta.kpr)
  names(beta.kpr) <- colnames(X)
  return(list(betahat.kpr = beta.kpr, lambda.1se = lambda[CVidx.QH.KPR1se],
              lambda.min = lambda[CVidx.QH.KPR]))

}
