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
#' \item{beta.hat}{A matrix of estimated coefficients for the penalized variables, where each column corresponds to a different value of lambda.}
#' \item{eta.hat}{A matrix of estimated coefficients for the unpenalized variables, where each column corresponds to a different value of lambda.}
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
  cov.missing <- missing(covariates)
  Z <- designMatrix # penalized matrix
  n <- nrow(Z)
  p <- ncol(Z) # number of penalized variables
  if (cov.missing) U <- matrix(0, n) # arbitrarily set U to a zero vector if no covariates are given
  else U <- as.matrix(covariates)

  if(missing(lambda))
    lambda <- exp(seq(from = 0, to = 10, length.out = n.lambda))
  n.lambda <- length(lambda)

  errors <- matrix(nrow = K, ncol = n.lambda)
  randidx <- sample(1:n, n)
  Yrand <- Y[randidx]
  Zrand <- Z[randidx, ]
  Urand <- as.matrix(U[randidx, ])
  Hrand <- H[randidx, randidx]
  for(j in 1:n.lambda){
    for(k in 1:K){
      Ytrain <- Yrand[-((n / K * (k - 1) + 1):(n / K * k))]
      Ytest <- Yrand[((n / K * (k - 1) + 1):(n / K * k))]
      Ztrain <- Zrand[-((n / K * (k - 1) + 1):(n / K * k)), ]
      Ztest <- Zrand[((n / K * (k - 1) + 1):(n / K * k)), ]
      Utrain <- Urand[-((n / K * (k - 1) + 1):(n / K * k)), ]
      Utest <- Urand[((n / K * (k - 1) + 1):(n / K * k)), ]
      Htrain <- Hrand[-((n / K * (k - 1) + 1):(n / K * k)), -((n / K * (k - 1) + 1):(n / K * k))]
      Htest <- Hrand[  ((n / K * (k - 1) + 1):(n / K * k)),  ((n / K * (k - 1) + 1):(n / K * k))]

      n.train <- nrow(Ztrain)

      if (cov.missing) P <- diag(n.train)
      else P <- diag(n.train) - (Utrain %*% solve( t(Utrain) %*% Htrain %*% Utrain ) %*% t(Utrain) %*% Htrain)
      Y.p <- P %*% Ytrain
      Z.p <- P %*% Ztrain

      beta.hat <- Q %*% solve( t(Z.p) %*% Htrain %*% Z.p %*% Q + lambda[j] * diag(p) ) %*% t(Z.p) %*% Htrain %*% Y.p # penalized coefficients
      if (cov.missing) eta.hat <- 0
      else eta.hat <- solve(t(Utrain) %*% Htrain %*% Utrain) %*% t(Utrain) %*% Htrain %*% (Ytrain - Ztrain %*% beta.hat) # unpenalized coefficients

      coefficients <- c(beta.hat, eta.hat)
      Xtest <- cbind(Ztest, Utest)

      yhat <- Xtest %*% coefficients
      errors[k, j] <- t(Ytest - yhat) %*% Htest %*%  (Ytest - yhat)
    }
  }
  lambda.min.index <- which.min(colSums(errors))
  lambda.min <- lambda[lambda.min.index]

  se.QH.KPR <- sd(errors[, lambda.min.index]) * sqrt(K)
  lambda.1se.indices <- which(colSums(errors) < lambda.min + se.QH.KPR) # this gets all of the indices of the lambda values withing one standard error of the min

  lambda.1se <- max(lambda[lambda.1se.indices])
  lambda.1se.index <- which(lambda == lambda.1se)

  estimates <- sapply(lambda, FUN = function(s) {
      if (cov.missing) P <- diag(n)
      else P <- diag(n) - U %*% solve(t(U) %*% H %*% U) %*% t(U) %*% H
      Y.p <- P %*% Y
      Z.p <- P %*% Z # apply P to the variables that should be penalized

      beta.hat <- Q %*% solve( t(Z.p) %*% H %*% Z.p %*% Q + s * diag(p) ) %*% t(Z.p) %*% H %*% Y.p # penalized coefficients
      if (cov.missing) eta.hat <- 0
      else eta.hat <- solve(t(U) %*% H %*% U) %*% t(U) %*% H %*% (Y - Z %*% beta.hat)
      c(beta.hat, eta.hat)
  })
  beta.hat <- estimates[1:p,]
  eta.hat <- estimates[-(1:p),]

  rownames(beta.hat) <- colnames(Z)
  rownames(eta.hat) <- colnames(U)
  if (all(eta.hat == rep(0,length(lambda)))) eta.hat <- NULL



  return(list(beta.hat = beta.hat,
              eta.hat = eta.hat,
              lambda = lambda,
              lambda.min = lambda.min,
              lambda.min.index = lambda.min.index,
              lambda.1se = lambda.1se,
              lambda.1se.index = lambda.1se.index))

}
