#' Kernel Penalized Regression
#'
#' Fits a kernel penalized regression model using a design matrix X, response vector Y, sample similarity kernel H, and variable similarity kernel Q.
#'
#' @param X An n x p data matrix. Should be scaled and centered.
#' @param Y An n x 1 response vector. Should be scaled and centered.
#' @param H An n x n sample similarity kernel. Must be symmetric and positive semidefinite. This defaults to an identity matrix.
#' @param Q A p x p variable similarity kernel. Must be symmetric and postive semidefinite. This defaults to an identity matrix.
#' @param n.lambda The number of lambda values to test through the cross validation search. The values are generated internally.
#' @param lambda A vector of lambda values to test through cross validation. This will override the sequence generated with the \code{n.lambda} parameter.
#' @param K The number of folds in the cross validation search.
#' @return
#' \item{betahat}{A matrix of estimated coefficients, where each column corresponds to a different value of lambda.}
#' \item{lambda}{The vector of lambda values used in cross validation.}
#' \item{lambda.min}{The value of lambda that results in the minimum sum of squared errors.}
#' \item{lambda.min.index}{The index of \code{lambda.min} in the \code{lambda} vector. This can be used to retrieve the estimated coefficients corresponding
#' to \code{lambda.min} from the \code{betahat} matrix.}
#' \item{lambda.1se}{The largest value of lambda within one standard error of \code{lambda.min}.}
#' \item{lambda.1se.index}{The index of the \code{lambda.1se} value in the \code{lambda} vector.}
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
  CVidx.QH.KPR1se <- which(colSums(ERR.QH.KPR) < min(colSums(ERR.QH.KPR)) + se.QH.KPR) # this gets all of the indices of the lambda values withing one standard error of the min

  beta.kpr <- sapply(lambda, FUN = function(s) {
      Q %*% solve( t(X) %*% H %*% X %*% Q  + s * diag(p)) %*% t(X) %*% H %*% Y
  })

  if (length(lambda) == 1) beta.kpr <- as.vector(beta.kpr) # if only one lambda was given, just return a vector rather than a matrix

  rownames(beta.kpr) <- colnames(X)
  return(list(betahat = beta.kpr,
              lambda = lambda,
              lambda.min = lambda[CVidx.QH.KPR],
              lambda.min.index = CVidx.QH.KPR,
              lambda.1se = max(lambda[CVidx.QH.KPR1se]),
              lambda.1se.index = which.max(lambda[CVidx.QH.KPR1se]))

}
