#' Kernel Penalized Regression
#'
#' Fits a kernel penalized regression model using a design matrix X, response vector Y, sample similarity kernel H, and variable similarity kernel Q.
#'
#' @param designMatrix An n x p data matrix, consisting of variables that should be penalized by \code{Q}. Should be scaled and centered.
#' @param covariates An n x r data matrix, consisting of variables that should not be penalized. Should be scaled and centered.
#' @param Y An n x 1 response vector. Should be scaled and centered.
#' @param H An n x n sample similarity kernel. Must be symmetric and positive semidefinite. This defaults to an identity matrix.
#' @param Q A p x p variable similarity kernel. Must be symmetric and postive semidefinite. This defaults to an identity matrix.
#' @param REML logical, indicates whether to use REML to find the optimal value of \code{lambda}. If a value of \code{lambda} is given, this is ignored.
#' @param n.lambda The number of lambda values to test through the cross validation search. The values are generated internally.
#' @param lambda A vector of lambda values to test through cross validation. This will override the sequence generated with the \code{n.lambda} parameter.
#' @param K The number of folds in the cross validation search.
#' @param useCpp Indicate whether to use the C++ backend for cross-validation.
#' @param seed Set a seed for random number generation.
#' @return
#' \item{beta.hat}{A matrix of estimated coefficients for the penalized variables, where each column corresponds to a different value of lambda.}
#' \item{eta.hat}{A matrix of estimated coefficients for the unpenalized variables, where each column corresponds to a different value of lambda.}
#' \item{lambda}{The vector of lambda values used in cross validation.}
#' \item{lambda.min}{The value of lambda that results in the minimum sum of squared errors.}
#' \item{lambda.min.index}{The index of \code{lambda.min} in the \code{lambda} vector. This can be used to retrieve the estimated coefficients corresponding
#' to \code{lambda.min} from the \code{betahat} matrix.}
#' \item{lambda.1se}{The largest value of lambda within one standard error of \code{lambda.min}.}
#' \item{lambda.1se.index}{The index of the \code{lambda.1se} value in the \code{lambda} vector.}
#' \item{cv.errors}{The error matrix generated in the cross-validation procedure.}
#' \item{REML}{Was REML used to find \code{lambda}?}
#' @importFrom stats sd
#' @importFrom Rcpp sourceCpp
#' @importFrom nlme lme pdIdent VarCorr
#' @useDynLib KPR, .registration = TRUE
#' @export
KPR <- function(designMatrix, covariates, Y, H = diag(nrow(designMatrix)), Q = diag(ncol(designMatrix)),
                REML = TRUE, n.lambda = 200, lambda, K = 5, useCpp = TRUE, seed)
{
  # eigen.Q <- eigen(Q)
  # Q <- (1/eigen.Q$values[1]) * Q # standardize Q
  # P.Q <- eigen.Q$vectors %*% t(eigen.Q$vectors)
  #
  # XU <- designMatrix %*% eigen.Q$vectors
  # XU.std <- apply(X, 2, function(x) x / norm(x, type="2"))



  cov.missing <- missing(covariates)
  Z <-  designMatrix # XU.std %*% t(eigen.Q$vectors)# penalized matrix
  n <- nrow(Z)
  p <- ncol(Z) # number of penalized variables
  if (cov.missing) E <- matrix(0, n) # arbitrarily set E to a zero vector if no covariates are given
  else E <- as.matrix(covariates)

  if (REML & missing(lambda)) # we don't need to use REML if lambda is given.
  {
    lambda <- remlEstimation(Z, E, Y, H, Q, cov.missing)
    lambda.min <- NULL
    lambda.min.index <- NULL
    lambda.1se <- NULL
    lambda.1se.index <- NULL
    errors <- NULL

  }
  else
  {
    REML <- FALSE
    if(missing(lambda))
      lambda <- exp(seq(from = 0, to = 10, length.out = n.lambda))
    n.lambda <- length(lambda)

    if (n.lambda == 1)
    {
      lambda.min <- NULL
      lambda.min.index <- NULL
      lambda.1se <- NULL
      lambda.1se.index <- NULL
      errors <- NULL
    }
    else
    {
      if (!missing(seed)) set.seed(seed)
      randidx <- sample(1:n, n)
      Yrand <- Y[randidx]
      Zrand <- Z[randidx, ]
      Erand <- as.matrix(E[randidx, ])
      Hrand <- H[randidx, randidx]
      if (useCpp) errors <- computeErrorMatrixCpp(Zrand, Erand, Yrand, Hrand, Q, lambda, K, cov.missing)
      else errors <- computeErrorMatrixR(Zrand, Erand, Yrand, Hrand, Q, lambda, K, cov.missing)
      lambda.min.index <- which.min(colSums(errors))
      lambda.min <- lambda[lambda.min.index]

      se.QH.KPR <- sd(errors[, lambda.min.index]) * sqrt(K)
      lambda.1se.indices <- which(colSums(errors) < sum(errors[,lambda.min.index]) + se.QH.KPR) # this gets all of the indices of the lambda values withing one standard error of the min

      lambda.1se <- max(lambda[lambda.1se.indices])
      lambda.1se.index <- which(lambda == lambda.1se)
    }

  }

  estimates <- sapply(lambda, FUN = function(s) {
      if (cov.missing) P <- diag(n)
      else P <- diag(n) - E %*% solve(t(E) %*% H %*% E) %*% t(E) %*% H
      Y.p <- P %*% Y
      Z.p <- P %*% Z # apply P to the variables that should be penalized

      beta.hat <- Q %*% solve( t(Z.p) %*% H %*% Z.p %*% Q + s * diag(p) ) %*% t(Z.p) %*% H %*% Y.p # penalized coefficients
      if (cov.missing) eta.hat <- 0
      else eta.hat <- solve(t(E) %*% H %*% E) %*% t(E) %*% H %*% (Y - Z %*% beta.hat)
      c(beta.hat, eta.hat)
  })
  beta.hat <- estimates[1:p,]
  eta.hat <- estimates[-(1:p),]

  rownames(beta.hat) <- colnames(Z)
  rownames(eta.hat) <- colnames(E)

  if (all(eta.hat == rep(0,length(lambda)))) eta.hat <- NULL
  if (all(E == matrix(0, n))) E <- NULL


  output <- list(Z = Z,
              E = E,
              Y = Y,
              H = H,
              Q = Q,
              beta.hat = beta.hat,
              eta.hat = eta.hat,
              lambda = lambda,
              lambda.min = lambda.min,
              lambda.min.index = lambda.min.index,
              lambda.1se = lambda.1se,
              lambda.1se.index = lambda.1se.index,
              cv.errors = errors,
              REML = REML)
  class(output) <- "KPR"

  return(output)

}




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

remlEstimation <- function(Z, E, Y, H, Q, cov.missing)
{

  n <- nrow(Z)
  p <- ncol(Z)

  if (cov.missing) P <- diag(n)
  else P <- diag(n) - E %*% solve(t(E) %*% H %*% E) %*% t(E) %*% H
  Y.p <- P %*% Y
  Z.p <- P %*% Z # apply P to the variables that should be penalized

  eigen.Q <- eigen(Q)
  L.Q <- eigen.Q$vectors %*% diag(sqrt(eigen.Q$values))
  eigen.H <- eigen(H)
  L.H <- eigen.H$vectors %*% diag(sqrt(eigen.H$values))

  Z.tilde <- t(L.H) %*% Z.p %*% L.Q
  Y.tilde <- t(L.H) %*% Y.p

  dummyID <- factor(rep(1, nrow(Z)))
  lmm.fit <- suppressWarnings(lme(Y.tilde~1, random=list(dummyID = pdIdent(~-1 + Z.tilde))))
  lambda.reml <- lmm.fit$sigma^2 / as.numeric(VarCorr(lmm.fit)[1,1])^2

  return(lambda.reml)

}
