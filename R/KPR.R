#' Kernel Penalized Regression
#'
#' Fits a kernel penalized regression model using a design matrix X, response vector Y, sample similarity kernel H, and variable similarity kernel Q.
#'
#' @param designMatrix An n x p data matrix, consisting of variables that should be penalized by \code{Q}. Should be scaled and centered.
#' @param covariates An n x r data matrix, consisting of variables that should not be penalized. Should be scaled and centered.
#' @param Y An n x 1 response vector. Should be scaled and centered.
#' @param H An n x n sample similarity kernel. Must be symmetric and positive semidefinite. This defaults to an identity matrix.
#' @param Q A p x p variable similarity kernel. Must be symmetric and postive semidefinite. This defaults to an identity matrix.
#' @param lambda A fixed value of lambda to be used in fitting the KPR model. If this is missing, an appropriate value is chosen using REML estimation.
#' @param scale Logical, indicates whether to scale \code{Q} and the design matrix.
#' @param inference Logical, indicates whether to compute p-values for penalized regression coefficients with the GMD inference.
#' @param ... Additional parameters passed to the GMD inference
#' @return
#' \item{beta.hat}{Estimated coefficients for the penalized variables.}
#' \item{eta.hat}{Estimated coefficients for the unpenalized variables.}
#' \item{theta.hat}{The vector of parameter values found with Maximum Likelihood.}
#' \item{p.values}{P-values for each penalized coefficient, resulting from the GMD inference.}
#' \item{bound}{The stochastic bound used to compute each p-value.}
#' \item{sigmaepsi.hat}{Variance component estimate used in the GMD inference.}
#' @importFrom stats sd pnorm optimize
#' @importFrom nlme lme pdIdent VarCorr
#' @importFrom natural olasso_cv
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom optimx optimx
#' @references Randolph et al. (2018) The Annals of Applied Statistics
#' (\href{https://projecteuclid.org/euclid.aoas/1520564483}{Project Euclid})
#'
#' Wang et al. Technical report.
#' @useDynLib KPR, .registration = TRUE
#' @export
KPR <- function(designMatrix, covariates, Y, H = diag(nrow(designMatrix)), Q = diag(ncol(designMatrix)),
                lambda, scale = FALSE,
                inference = FALSE, REML = FALSE, ...)
{
  if (scale)
  {
    eigen.Q <- eigen(Q)
    Q <- (1/eigen.Q$values[1]) * Q # standardize Q

    XU <- designMatrix %*% eigen.Q$vectors
    XU.std <- apply(XU, 2, function(x) sqrt(length(x)) * x / sqrt(as.vector(t(x) %*% H %*% x) ))

    Z <- XU.std %*% t(eigen.Q$vectors)
    colnames(Z) <- colnames(designMatrix)
  }
  else Z <- designMatrix


  cov.missing <- missing(covariates)
  n <- nrow(Z)
  p <- ncol(Z) # number of penalized variables
  if (cov.missing) E <- matrix(0, n) # arbitrarily set E to a zero vector if no covariates are given
  else E <- as.matrix(covariates)

  if (cov.missing) P <- diag(n)
  else P <- diag(n) - E %*% solve(t(E) %*% H %*% E) %*% t(E) %*% H

  Y.p <- P %*% Y
  Z.p <- P %*% Z # apply P to the variables that should be penalized

  if (REML) # this is just for comparing the new method to the old one; should be removed eventually
  {
      theta.hat <- remlEstimation(Z.p, Y.p, H, Q)
      beta.hat <- Q %*% solve(t(Z.p) %*% H %*% Z.p %*% Q + theta.hat * diag(p)) %*% t(Z.p) %*% H %*% Y
  }
  else
  {
      theta.hat <- fitParameters(Z.p, Y.p, H, Q)
      beta.hat <- computeBetaHat(Z.p, Y.p, H, Q, theta.hat)
  }


  if (cov.missing) eta.hat <- 0
  else eta.hat <- solve(t(E) %*% H %*% E) %*% t(E) %*% H %*% (Y - Z %*% beta.hat)

  rownames(beta.hat) <- colnames(Z)
  rownames(eta.hat) <- colnames(E)

  if (cov.missing) eta.hat <- NULL
  if (cov.missing) E <- NULL

  output <- list(Z = Z,
              E = E,
              Y = Y,
              H = H,
              Q = Q,
              beta.hat = beta.hat,
              eta.hat = eta.hat,
              theta.hat = theta.hat)
  if (inference) {
    infer.out <- GMD.inference(output,...)
    output$p.values <- infer.out$p.values
    output$bound <- infer.out$bound
    output$sigmaepsi.hat <- infer.out$sigmaepsi.hat
  }
  else {
     output$p.values  <- NULL
     output$bound <- NULL
     output$sigmaepsi.hat <- NULL
  }
  class(output) <- "KPR"

  return(output)

}


