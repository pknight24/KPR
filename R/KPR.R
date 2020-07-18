#' Kernel Penalized Regression
#'
#' Fits a kernel penalized regression model using a design matrix X, response vector Y, sample similarity kernel H, and variable similarity kernel Q.
#'
#' @param designMatrix An n x p data matrix, consisting of variables that should be penalized by \code{Q}. Should be scaled and centered.
#' @param covariates An n x r data matrix, consisting of variables that should not be penalized. Should be scaled and centered.
#' @param Y An n x 1 response vector. Should be scaled and centered.
#' @param H An n x n sample similarity kernel. Must be symmetric and positive semidefinite. This defaults to an identity matrix.
#' @param Q A p x p variable similarity kernel. Must be symmetric and postive semidefinite. This defaults to an identity matrix.
#' @param scale Logical, indicates whether to scale \code{Q}, \code{H} and the design matrix to have a spectral norm of 1.
#' @param REML Logical, indicates whether to use REML estimation for finding the parameters. This will only work with a single H and Q matrix.
#' @return
#' \item{beta.hat}{Estimated coefficients for the penalized variables.}
#' \item{eta.hat}{Estimated coefficients for the unpenalized variables.}
#' \item{theta.hat}{The vector of tuning parameter values found with Maximum Likelihood.}
#' @importFrom stats sd pnorm optimize
#' @importFrom nlme lme pdIdent VarCorr
#' @importFrom natural olasso_cv
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom alabama constrOptim.nl
#' @references Randolph et al. (2018) The Annals of Applied Statistics
#' (\href{https://projecteuclid.org/euclid.aoas/1520564483}{Project Euclid})
#' @useDynLib KPR, .registration = TRUE
#' @export
KPR <- function(designMatrix, covariates = NULL, Y, H = diag(nrow(designMatrix)), Q = diag(ncol(designMatrix)),
                scale = FALSE, REML = FALSE)
{

  if (!is.list(Q)) Q <- list(Q) # this handles the case when q = h = 1
  if (!is.list(H)) H <- list(H)

  if (scale)
  {
      for (j in 1:length(Q))
      {
           eigen.Q <- eigen(Q[[j]])
           Q[[j]] <- eigen.Q$vectors %*% (eigen.Q$values/eigen.Q$values[1]) %*% t(eigen.Q$vectors) # standardize each Q
      }
      for (i in 1:length(H))
      {
          eigen.H <- eigen(H[[i]])
           H[[i]] <- eigen.H$vectors %*% (eigen.H$values/eigen.H$values[1]) %*% t(eigen.H$vectors) # standardize each H
      }


      eigen.Z <- eigen(designMatrix)
      Z <- eigen.Z$vectors %*% (eigen.Z$values/eigen.Z$values[1]) %*% t(eigen.Z$vectors) # standardize Z
    }
    else Z <- designMatrix

    cov.missing <- is.null(covariates)
    n <- nrow(Z)
    p <- ncol(Z) # number of penalized variables
    E <- covariates

    if (cov.missing) P <- diag(n)
    else P <- diag(n) - E %*% solve(t(E) %*% H %*% E) %*% t(E) %*% H

    Y.p <- P %*% Y
    Z.p <- P %*% Z # apply P to the variables that should be penalized

    if (REML) # this is just for comparing the new method to the old one; should be removed eventually
    {
        theta.hat <- remlEstimation(Z.p, Y.p, H[[1]], Q[[1]])
        beta.hat <- Q[[1]] %*% solve(t(Z.p) %*% H[[1]] %*% Z.p %*% Q[[1]] + theta.hat * diag(p)) %*% t(Z.p) %*% H[[1]] %*% Y
    }
    else # here is our new method
    {
        theta.hat <- findTuningParameters(Z.p, Y.p, H, Q)
        beta.hat <- computeCoefficientEstimates(Z.p, Y.p, H, Q, theta.hat)
        names(beta.hat) <- colnames(Z)
    }
    if (cov.missing)
    {
        eta.hat <- NULL
        E <- NULL
    }
    else
    {
        h <- length(h)
        H.sum <- abs(theta.hat$sigma[1]) * H[[1]]
        if (h > 1) for (i in 2:h) H.sum <- H.sum + abs(theta.hat$sigma[i])*H[[i]]
        eta.hat <- solve(t(E) %*% H.sum %*% E) %*% t(E) %*% H.sum %*% (Y - Z %*% beta.hat)
        names(eta.hat) <- colnames(E)
    }
    output <- list(Z = Z,
                E = E,
                Y = Y,
                H = H,
                Q = Q,
                beta.hat = beta.hat,
                eta.hat = eta.hat)
    if (REML) output$lambda <- theta.hat
    else
    {
        output$sigma <- theta.hat$sigma
        output$alpha <- theta.hat$alpha
        output$lambda <- theta.hat$c_H / theta.hat$c_Q
    }

    class(output) <- "KPR"

    return(output)

}


