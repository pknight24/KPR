#' Kernel Penalized Regression
#'
#' Fits a kernel penalized regression model using a design matrix X, response vector Y, sample similarity kernel H, and variable similarity kernel Q.
#'
#' @param designMatrix An n x p data matrix, consisting of variables that should be penalized by \code{Q}. Should be scaled and centered.
#' @param covariates An n x r data matrix, consisting of variables that should not be penalized. Should be scaled and centered.
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
KPR <- function(designMatrix, covariates, Y, H = diag(nrow(designMatrix)), Q = diag(ncol(designMatrix)),
                n.lambda = 200, lambda, K = 5)
{
  X <- cbind(designMatrix, covariates) # X is now n x (p + r)
  n <- nrow(designMatrix)
  p <- ncol(designMatrix) # number of penalized variables
  r <- ncol(covariates) # number of penalized variables
  if(missing(lambda))
    lambda <- exp(seq(from = 0, to = 10, length.out = n.lambda))
  n.lambda <- length(lambda)

  errors <- matrix(nrow = K, ncol = n.lambda)
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
      Htest <- Hrand[  ((n / K * (k - 1) + 1):(n / K * k)),  ((n / K * (k - 1) + 1):(n / K * k))]

      P <- diag(nrow(Xtrain)) - (Xtrain[,-(1:p)] %*% solve( t(Xtrain[,-(1:p)]) %*% Htrain %*% Xtrain[,-(1:p)] ) %*% t(Xtrain[,-(1:p)]) %*% Htrain)
      Y.p <- P %*% Ytrain
      Z.p <- P %*% Xtrain[,(1:p)] # apply P to the variables that should be penalized

      beta.p <- Q %*% solve( t(Z.p) %*% Htrain %*% Z.p %*% Q + lambda[j] * diag(p) ) %*% t(Z.p) %*% Htrain %*% Y.p # penalized coefficients
      beta.r <- solve(t(Xtrain[,-(1:p)]) %*% Htrain %*% Xtrain[,-(1:p)]) %*% t(Xtrain[,-(1:p)]) %*% Htrain %*% (Ytrain - Xtrain[,1:p] %*% beta.p) # unpenalized coefficients

      betahat <- c(beta.p, beta.r)

      yhat <- Xtest %*% betahat
      errors[k, j] <- t(Ytest - yhat) %*% Htest %*%  (Ytest - yhat)
    }
  }
  CVidx.QH.KPR <- which.min(colSums(errors))

  se.QH.KPR <- sd(errors[, CVidx.QH.KPR]) * sqrt(K)
  CVidx.QH.KPR1se <- which(colSums(errors) < min(colSums(errors)) + se.QH.KPR) # this gets all of the indices of the lambda values withing one standard error of the min

  beta.kpr <- sapply(lambda, FUN = function(s) {
      P <- diag(n) - X[,-(1:p)] %*% solve(t(X[,-(1:p)]) %*% H %*% X[,-(1:p)]) %*% t(X[,-(1:p)]) %*% H
      Y.p <- P %*% Y
      Z.p <- P %*% X[,(1:p)] # apply P to the variables that should be penalized

      beta.p <- Q %*% solve( t(Z.p) %*% H %*% Z.p %*% Q + s * diag(p) ) %*% t(Z.p) %*% H %*% Y.p # penalized coefficients
      beta.r <- solve(t(X[,-(1:p)]) %*% H %*% X[,-(1:p)]) %*% t(X[,-(1:p)]) %*% H %*% (Y - X[,1:p] %*% beta.p)
      c(beta.p, beta.r)
  })

  if (length(lambda) == 1) beta.kpr <- as.vector(beta.kpr) # if only one lambda was given, just return a vector rather than a matrix

  rownames(beta.kpr) <- colnames(X)
  return(list(betahat.penalized = beta.kpr[1:p,],
              betahat.unpenalized = beta.kpr[-(1:p),],
              lambda = lambda,
              lambda.min = lambda[CVidx.QH.KPR],
              lambda.min.index = CVidx.QH.KPR,
              lambda.1se = max(lambda[CVidx.QH.KPR1se]),
              lambda.1se.index = which.max(lambda[CVidx.QH.KPR1se])))

}
