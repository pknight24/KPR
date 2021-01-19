#' Running the GMD inference on a KPR model
#'
#' @param KPR.output Output from running the \code{KPR} function.
#' @param mu GMD inference parameter
#' @param r GMD inference parameter
#' @param weight Logical, indicates whether to include a penalty factor when computing the beta.glasso vector.
#' @param scale Logical, indicates whether to scale the design matrix with respect to the eigenvalues of the composite Q matrix.
#' @param ... Additional parameters passed to the GMD function.
#' @return An object of classes KPR with the following fields added:
#' \item{p.values}{P-values for each penalized coefficient, resulting from the GMD inference.}
#' \item{bound}{The stochastic bound used to compute each p-value.}
#' @export
inference <- function(KPR.output, mu = 1, r = 0.05, weight = TRUE, scale = FALSE, ...)
{

  Z <- KPR.output$X
  E <- KPR.output$E # for now, we will ignore the E matrix
  Y <- KPR.output$Y
  H <- KPR.output$H
  Q <- KPR.output$Q
  alpha <- KPR.output$alpha
  sigma <- KPR.output$sigma
  lambda <- KPR.output$lambda
  n <- dim(Z)[1]
  p <- dim(Z)[2]
  q <- length(Q)
  h <- length(H)
  beta.KPR <- KPR.output$beta.hat # before correction

  # we first need to form the composite H and Q matrices
  H.sum <- sigma[1] * H[[1]]
  if (h > 1) for (i in 2:h) H.sum <- H.sum + sigma[i]*H[[i]]
  Q.sum <- alpha[1] * Q[[1]]
  if (q > 1) for (j in 2:q) Q.sum <- Q.sum + alpha[j]*Q[[j]]

  H <- H.sum
  Q <- Q.sum



  if (scale)
  {
    eigen.Q <- eigen(Q)
    Q <- (1/eigen.Q$values[1]) * Q # standardize Q

    XU <- Z %*% eigen.Q$vectors
    XU.std <- apply(XU, 2, function(x) sqrt(length(x)) * x / sqrt(as.vector(t(x) %*% H %*% x) ))

    Z <- XU.std %*% t(eigen.Q$vectors)
  }



  if (is.null(E)) P <- diag(n)
  else P <- diag(n) - E %*% solve(t(E) %*% H %*% E) %*% t(E) %*% H
  Y.p <- P %*% Y
  Z.p <- P %*% Z

  gmd.out <- GMD(X = Z.p, H = H, Q = Q, K = sum(svd(Z.p)$d > 10^-10), ...)

  U <- gmd.out$U
  V <- gmd.out$V
  D <- gmd.out$S

  W <- D^2 / (D^2 + lambda)

  ## bias-correction
  eigen.H = eigen(H)
  values.H = eigen.H$values
  vectors.H = eigen.H$vectors
  L.H = vectors.H[,values.H > 0]%*%diag(sqrt(values.H[values.H > 0]))%*%t(vectors.H[, values.H > 0])
  L.H.inv = vectors.H[,values.H > 0]%*%diag(1/sqrt(values.H[values.H > 0]))%*%t(vectors.H[, values.H > 0])

  eigen.Q = eigen(Q)
  values.Q = eigen.Q$values[eigen.Q$values > 0]
  vectors.Q = eigen.Q$vectors[, eigen.Q$values > 0]
  L.Q = vectors.Q%*%diag(sqrt(values.Q))

  Z.tilde = t(L.H)%*%Z.p%*%vectors.Q
  Y.tilde = t(L.H)%*%Y.p

  ### using natural lasso method
  ##olasso.fit = natural::olasso_cv(Z.tilde, Y.tilde, nfold = 3)
  ##sigmaepsi.hat = olasso.fit$sig_obj
  ##for (i in 1:49) sigmaepsi.hat <- c(sigmaepsi.hat,
  ##                                    natural::olasso_cv(Z.tilde, Y.tilde, nfold = 3)$sig_obj)
  ##sigmaepsi.hat <- median(sigmaepsi.hat)

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

  beta.lasso = as.numeric(vectors.Q%*%beta.glasso)
  Xi = diag(Q%*%V%*%diag(W)%*%t(V))
  beta.cor = Q%*%V%*%diag(W)%*%t(V)%*%beta.lasso - (1 - mu)*Xi*beta.lasso - mu*beta.lasso
  beta.hat = beta.KPR  - beta.cor

  # covariance
  cov.hat = Q%*%V%*%diag(D^(-2)*W*W)%*%t(V)%*%Q
  diag.cov.hat = diag(cov.hat)
  bound.mat = (Q%*%V%*%diag(W)%*%t(V) - (1- mu)*diag(Xi) - mu*diag(rep(1,p)))%*%vectors.Q
  bound.hat = apply(bound.mat, 1, function(x){max(abs(x))})*(log(p)/n)^(0.5 - r) # sparsity parameter

  # p-values
  beta.temp = abs(beta.hat) - bound.hat
  p.vec = rep(0,p)
  for(i in 1:p){

    if(beta.temp[i] > 0){p.vec[i] = 2*(1 - pnorm(beta.temp[i]/sqrt(diag.cov.hat[i])))}
    if(beta.temp[i] <= 0){p.vec[i] = 1}

  }

    names(p.vec) <- colnames(Z)

    KPR.output$p.values <- p.vec
    KPR.output$bound <- bound.hat

    return(KPR.output)
}
