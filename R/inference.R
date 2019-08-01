#' Inference for Kernel Penalized Regression models
#'
#' Implementation of various inference methods for high dimensional regression. Currently, the available methods are the GMD inference (used by default), the Grace inference method, and the ridge projection method, from the \code{hdi} package. To select the Grace method, set \code{method = "Grace"} and to select the ridge projection method, set \code{method = "hdi"}.
#'
#' Before calling \code{ridge.proj}, we transform the data with a Cholesky decomposition, which allows us to formulate the model as a ridge problem while still including information about H and Q. However, this means that when a non-trivial Q is included in the KPR model, \code{ridge.proj} will return p-values corresponding to the DPCoA estimates, not the true KPR estimates. This may slightly affect interpretability.
#'
#' The value of lambda used in all inference methods is determined in the KPR model fit.
#'
#' @param KPR.output The output of running the \code{KPR} function.
#' @param method A string specifying the inference method to use.
#' @param ... Additional parameters corresponding to specific inference methods.
#' @importFrom natural olasso_cv
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom stats pnorm
#' @importFrom hdi ridge.proj
#' @importFrom Grace grace.test
#' @references Wang et al. Technical report.
#'
#' Dezeure et al. (2014) Statistical Science (\href{https://arxiv.org/abs/1408.4026}{arXiv})
#'
#' Zhao and Shojaie. (2016) Biometrics (\href{https://arxiv.org/abs/1506.08339}{arXiv})
#' @export
inference <- function(KPR.output, method = "GMD", ...)
{
  if (class(KPR.output) != "KPR") stop("Input is not of class 'KPR'")

  if (method == "GMD") GMD.inference(KPR.output, ...)
  else if (method == "hdi") hdi.ridge.inference(KPR.output, ...)
  else if (method == "Grace") grace.inference(KPR.output, ...)
  else stop("Inference method not recognized")
}

GMD.inference <- function(KPR.output, mu = 1, r = 0.05, weight = TRUE, numComponents = 10)
{
  Z <- KPR.output$Z
  E <- KPR.output$E # for now, we will ignore the E matrix
  Y <- KPR.output$Y
  H <- KPR.output$H
  Q <- KPR.output$Q
  if (length(KPR.output$lambda) == 1) lambda <- KPR.output$lambda
  else lambda  <- c(KPR.output$lambda.min, KPR.output$lambda.1se)
  n <- dim(Z)[1]
  p <- dim(Z)[2]
  if (length(KPR.output$lambda) == 1) beta.hat.uncorrected <- KPR.output$beta.hat # before correction
  else beta.hat.uncorrected <- KPR.output$beta.hat[,
                                                   c(KPR.output$lambda.min.index,
                                                     KPR.output$lambda.1se.index)]


  if (is.null(E)) P <- diag(n)
  else P <- diag(n) - E %*% solve(t(E) %*% H %*% E) %*% t(E) %*% H
  Y.p <- P %*% Y
  Z.p <- P %*% Z

  gmd.out <- GMD(X = Z.p, H = H, Q = Q, K = numComponents)

  U <- gmd.out$U
  V <- gmd.out$V
  S <- gmd.out$S

  W.long <- sapply(lambda, function(lam){ # each column corresponds to a value of lambda, each row corresponds to a diagonal value in the W matrix
    sapply(S, function(s) s^2 / (s^2 + lam)^2)
  })


  # bias-correction
  eigen.H = eigen(H)
  values.H = eigen.H$values
  vectors.H = eigen.H$vectors
  L.H = vectors.H%*%diag(sqrt(values.H))
  eigen.Q = eigen(Q)
  values.Q = eigen.Q$values
  vectors.Q = eigen.Q$vectors
  L.Q = vectors.Q%*%diag(sqrt(values.Q))

  Z.tilde = t(L.H)%*%Z.p%*%vectors.Q
  Y.tilde = t(L.H)%*%Y.p

  # using natural lasso method

  olasso.fit = natural::olasso_cv(Z.tilde, Y.tilde, nfold = 3)
  sigmaepsi.hat = olasso.fit$sig_obj

  if(weight == TRUE)
  {
    lambda.opt = glmnet::cv.glmnet(Z.tilde, Y.tilde, alpha = 1, intercept = FALSE, standardize = TRUE, penalty.factor = 1/sqrt(values.Q))$lambda.min
    beta.glasso = glmnet::glmnet(Z.tilde, Y.tilde, alpha = 1, lambda = lambda.opt, intercept = FALSE, standardize = TRUE, penalty.factor = 1/sqrt(values.Q))$beta
  }
  if(weight == FALSE)
  {
    lambda.opt = glmnet::cv.glmnet(Z.tilde, Y.tilde, alpha = 1, intercept = FALSE, standardize = TRUE)$lambda.min
    beta.glasso = glmnet::glmnet(Z.tilde, Y.tilde, alpha = 1, lambda = lambda.opt, intercept = FALSE, standardize = TRUE)$beta
  }

  beta.init = as.numeric(vectors.Q%*%beta.glasso)

  Xi.long <- sapply(1:length(lambda), FUN=function(s){ # each column is a lambda value, each row is a Xi value
    diag(Q%*%V%*%diag(W.long[,s])%*%t(V))
  })

  bias.hat <- sapply(1:length(lambda), FUN=function(s){
    Q%*%V%*%diag(W.long[,s])%*%t(V)%*%beta.init - (1 - mu)*Xi.long[,s]*beta.init - mu*beta.init
  })
  beta.hat.cor <- beta.hat.uncorrected  - bias.hat


  # covariance
  diag.cov.hat.long <- sapply(1:length(lambda),function(s){ # this finds the value of Sigma_jj * sigmaepsi_hat^2, used in the calculation of the p values
    diag(sigmaepsi.hat^2*Q%*%V%*%diag(S^(-2)*W.long[,s]*W.long[,s])%*%t(V)%*%Q)
  })

  bound.hat.long <- sapply(1:length(lambda), function(s){ # each column is a lambda value, each row j is the value of Psi_j (defined in Thm 3.1)
    bound.mat <- (Q%*%V%*%diag(W.long[,s])%*%t(V) - (1- mu)*diag(Xi.long[,s]) - mu*diag(rep(1,p)))%*%L.Q
    apply(bound.mat, 1, function(x){max(abs(x))})*(log(p)/n)^(0.5 - r) # sparsity parameter
  })

  # p-values compuated as given in 3.6
  p.mat <- sapply(1:length(lambda), function(s){
    beta.temp = abs(beta.hat.cor[,s]) - bound.hat.long[,s]
    p.vec = rep(0,p)
    for(i in 1:p){

      if(beta.temp[i] > 0){p.vec[i] = 2*(1 - pnorm(beta.temp[i]/sqrt(diag.cov.hat.long[i,s])))}
      if(beta.temp[i] <= 0){p.vec[i] = 1}

    }
    return(p.vec)
  })

  rownames(p.mat) <- colnames(Z)
  if (length(lambda) == 2) colnames(p.mat) <- c("lambda.min", "lambda.1se")

  return(p.mat)
}


hdi.ridge.inference <- function(KPR.output, ...)
{
    Z <- KPR.output$Z
    Y <- KPR.output$Y
    E <- KPR.output$E
    H <- KPR.output$H
    Q <- KPR.output$Q

    n <- nrow(Z)
    p <- ncol(Z)

    if (is.null(E)) P <- diag(n)
    else P <- diag(n) - E %*% solve(t(E) %*% H %*% E) %*% t(E) %*% H
    Y.p <- P %*% Y
    Z.p <- P %*% Z

    eigen.Q <- eigen(Q)
    L.Q <- eigen.Q$vectors %*% diag(sqrt(eigen.Q$values))
    eigen.H <- eigen(H)
    L.H <- eigen.H$vectors %*% diag(sqrt(eigen.H$values))

    Z.tilde <- t(L.H) %*% Z.p %*% L.Q
    Y.tilde <- t(L.H) %*% Y.p

    if (length(KPR.output$lambda) == 1) lambda <- KPR.output$lambda
    else lambda <- c(KPR.output$lambda.min, KPR.output$lambda.1se)

    pval <- sapply(lambda, function(s) ridge.proj(x = Z.tilde, y = Y.tilde, lambda = s, standardize = FALSE, ...)$pval)
    rownames(pval) <- colnames(Z)
    if (length(lambda) > 1) colnames(pval) <- c("lambda.min", "lambd.1se")
    pval

}

grace.inference <- function(KPR.output, ...)
{
    Z <- KPR.output$Z
    Y <- KPR.output$Y
    E <- KPR.output$E
    H <- KPR.output$H
    Q <- KPR.output$Q

    n <- nrow(Z)
    p <- ncol(Z)

    if (is.null(E)) P <- diag(n)
    else P <- diag(n) - E %*% solve(t(E) %*% H %*% E) %*% t(E) %*% H
    Y.p <- P %*% Y
    Z.p <- P %*% Z

    eigen.H <- eigen(H)
    L.H <- eigen.H$vectors %*% diag(sqrt(eigen.H$values))

    Z.tilde <- t(L.H) %*% Z.p
    Y.tilde <- t(L.H) %*% Y.p

    if (length(KPR.output$lambda) == 1) lambda <- KPR.output$lambda
    else lambda <- c(KPR.output$lambda.min, KPR.output$lambda.1se)

    pval <- sapply(lambda, function(s) grace.test(X = Z.tilde, Y = as.vector(Y.tilde), L = solve(Q), lambda.L = s, C = "cv", ...)$pvalue)
    rownames(pval) <- colnames(Z)
    if (length(lambda) > 1) colnames(pval) <- c("lambda.min", "lambd.1se")
    pval

}
