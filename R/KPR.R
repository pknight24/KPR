#' Kernel Penalized Regression
#'
#' Fits a kernel penalized regression model using a design matrix X, response vector Y, sample similarity kernels H_1, H_2, ..., H_h, and variable similarity kernels Q_1, ..., Q_q.
#'
#' @param X An n x p data matrix, consisting of variables that should be penalized by the Q matrices. Should be scaled and centered.
#' @param E An n x r data matrix, consisting of variables that should not be penalized. Should be scaled and centered.
#' @param Y An n x 1 response vector. Should be scaled and centered.
#' @param H A list of n x n sample similarity kernels. If only one matrix is included in the model, it does not need to be wrapped as a list. All matrices must be symmetric positive semidefinite. This defaults to a single identity matrix.
#' @param Q A list of p x p variable similarity kernels. If only one matrix is included in the model, it does not need to be wrapped as a list. All matrices must be symmetric positive semidefinite. This defaults to a single identity matrix.
#' @param scale Logical, indicates whether to scale all the Q's, H's and the design matrix to have a spectral norm of 1.
#' @param REML Logical, indicates whether to use REML estimation for finding the parameters. This will only work with a single H and Q matrix, and is the preferred method in this case.
#' @param Q.inv Logical, indicates whether to penalize Q_composite or Q_composite inverse.
#' @param control.outer A list of parameters used by the outer loop in `constrOptim.nl`. This is only used when `REML = FALSE`.
#' @param control.optim A list of parameters used by the inner loop in `constrOptim.nl`.
#' @return
#' \item{beta.hat}{Estimated coefficients for the penalized variables.}
#' \item{eta.hat}{Estimated coefficients for the unpenalized variables.}
#' \item{lambda}{The optimal lambda parameter estimated with maximum likelihood.}
#' \item{alpha}{The vector of optimal weights corresponding to the Q matrices.}
#' \item{sigma}{The vector of optimal weights corresponding to the H matrices.}
#' @importFrom stats sd pnorm median optim
#' @importFrom nlme lme pdIdent VarCorr
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom alabama constrOptim.nl
#' @references Randolph et al. (2018) The Annals of Applied Statistics
#' (\href{https://projecteuclid.org/euclid.aoas/1520564483}{Project Euclid})
#' @export
KPR <- function(X, E = NULL, Y, H = diag(nrow(X)), Q = diag(ncol(X)),
                scale = FALSE, REML = FALSE, Q.inv = TRUE, control.outer = list(trace=FALSE, NMinit = TRUE, method = "BFGS"), 
                control.optim = list())
{
  
   # this handles the case when q = h = 1
  if (!is.list(Q) & !is.list(H)) REML <- TRUE # we should default to REML in this case
  if (!is.list(Q)) Q <- list(Q)
  if (!is.list(H)) H <- list(H)
  
  
  

  if (scale)
  {
      for (j in 1:length(Q))
      {
           eigen.Q <- eigen(Q[[j]])
           Q[[j]] <- eigen.Q$vectors %*% diag((eigen.Q$values/eigen.Q$values[1])) %*% t(eigen.Q$vectors) # standardize each Q
      }
      for (i in 1:length(H))
      {
          eigen.H <- eigen(H[[i]])
           H[[i]] <- eigen.H$vectors %*% diag((eigen.H$values/eigen.H$values[1])) %*% t(eigen.H$vectors) # standardize each H
      }


      svd.Z <- svd(X)
      Z <- svd.Z$u %*% diag((svd.Z$d/svd.Z$d[1])) %*% t(svd.Z$v) # standardize Z
    }
    else Z <- X

    cov.missing <- is.null(E)
    n <- nrow(Z)
    p <- ncol(Z) # number of penalized variables

    if (cov.missing) P <- diag(n)
    else P <- diag(n) - E %*% solve(t(E) %*% H[[1]] %*% E) %*% t(E) %*% H[[1]] # how do we build this projection matrix with multiple H matrices? discuss w Tim

    Y.p <- P %*% Y
    Z.p <- P %*% Z # apply P to the variables that should be penalized

    if (REML)
    {
        theta.hat <- remlEstimation(Z.p, Y.p, H[[1]], Q[[1]])
        beta.hat <- Q[[1]] %*% solve(t(Z.p) %*% H[[1]] %*% Z.p %*% Q[[1]] + theta.hat * diag(p)) %*% t(Z.p) %*% H[[1]] %*% Y
    }
    else # here is our new method
    {
        theta.hat <- findTuningParameters(Z.p, Y.p, H, Q, Q.inv = Q.inv, control.outer = control.outer, 
                                          control.optim = control.optim)
        beta.hat <- computeCoefficientEstimates(Z.p, Y.p, H, Q, theta.hat, Q.inv = Q.inv)

        names(beta.hat) <- colnames(Z)
    }
    if (cov.missing)
    {
        eta.hat <- NULL
        E <- NULL
    }
    else
    {
       # h <- length(H)
       # H.sum <- abs(theta.hat$sigma[1]) * H[[1]]
       # if (h > 1) for (i in 2:h) H.sum <- H.sum + abs(theta.hat$sigma[i])*H[[i]]
        eta.hat <- solve(t(E) %*% H[[1]] %*% E) %*% t(E) %*% H[[1]] %*% (Y - Z %*% beta.hat)
        names(eta.hat) <- colnames(E)
    }
    output <- list(X = Z,
                E = E,
                Y = Y,
                H = H,
                Q = Q,
                beta.hat = beta.hat,
                eta.hat = eta.hat)
    if (REML)
    {
        output$lambda <- theta.hat
        output$sigma <- 1
        output$alpha <- 1
    }
    else
    {
        output$sigma <- theta.hat$sigma
        output$alpha <- theta.hat$alpha
        output$lambda <- theta.hat$c_H / theta.hat$c_Q
    }

    class(output) <- "KPR"

    return(output)

}


