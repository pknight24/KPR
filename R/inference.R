#' @importFrom natural olasso_cv
#' @export
inference <- function(KPR.output, method = "GMD", ...)
{
  if (method == "GMD") GMD.inference(KPR.output, ...)
}

GMD.inference <- function(KPR.output, mu, r, weight = TRUE, numComponents = 10)
{
  Z <- KPR.output$Z
  E <- KPR.output$E # for now, we will ignore the E matrix
  Y <- KPR.output$Y
  H <- KPR.output$H
  Q <- KPR.output$Q
  lambda <- KPR.output$lambda
  n <- dim(Z)[1]
  p <- dim(Z)[2]
  beta.hat.uncorrected = KPR.output$beta.hat # before correction


  gmd.out <- GMD(X = Z, H = H, Q = Q, K = numComponents)

  U <- gmd.out$U
  V <- gmd.out$V
  S <- gmd.out$D

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

  Z.tilde = L.H%*%Z%*%vectors.Q
  Y.tilde = L.H%*%Y

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
  return(beta.init) # beta.init is being computed without issues

  # TODO: figure out how to incorporate W.long, so that we can keep track of the estimates at each provided lambda value

  Xi <- diag(Q%*%V%*%diag(W)%*%t(V))
  bias.hat <- Q%*%V%*%diag(W)%*%t(V)%*%beta.init - (1 - mu)*Xi*beta.init - mu*beta.init
  beta.hat <- beta.hat.uncorrected  - bias.hat

  # covariance
  cov.hat = sigmaepsi.hat^2*Q%*%V%*%diag(D^(-2)*W*W)%*%t(V)%*%Q
  diag.cov.hat = diag(cov.hat)
  bound.mat = (Q%*%V%*%diag(W)%*%t(V) - (1- mu)*diag(Xi) - mu*diag(rep(1,p)))%*%L.Q
  bound.hat = apply(bound.mat, 1, function(x){max(abs(x))})*(log(p)/n)^(0.5 - r) # sparsity parameter
  # print(bound.hat)

  # p-value
  beta.temp = abs(beta.hat) - bound.hat
  p.vec = rep(0,p)
  for(i in 1:p){
    # print(beta.temp)

    if(beta.temp[i] > 0){p.vec[i] = 2*(1 - pnorm(beta.temp[i]/sqrt(diag.cov.hat[i])))}
    if(beta.temp[i] <= 0){p.vec[i] = 1}

  }
  names(p.vec) <- colnames(X)
  return(p.vec)

}

GMD <- function(X, H, Q, K)
{
  n = dim(X)[1]
  p = dim(X)[2]
  # output matrix/vec
  U = matrix(0, n, K)
  V = matrix(0, p, K)
  D = rep(0,K)

  X_0 = X
  u_0 = c(1, rep(0,n-1))
  v_0 = c(1, rep(0,p-1))

  for(iter in 1:K){


    error = 1

    while(error > 1e-5){

      temp.uv = get_uv(X_0, H, Q, u_0, v_0)

      u = temp.uv$u
      v = temp.uv$v

      error = norm(u - u_0, "2") + norm(v - v_0, "2")
      #print(error)

      u_0 = u
      v_0 = v

    }

    U[,iter] = u
    V[,iter] = v
    d = t(u)%*%H%*%X_0%*%Q%*%v
    D[iter] =  d

    X_0 = X_0 - u%*%d%*%t(v)

  }

  return(list(U = U, V = V, D = D, H = H, Q = Q))

}

get_uv = function(X, H, Q, u_0, v_0){

  u = X%*%Q%*%v_0/as.numeric(sqrt(t(v_0)%*%t(Q)%*%t(X)%*%H%*%X%*%Q%*%v_0))
  v = t(X)%*%H%*%u/as.numeric(sqrt(t(u)%*%t(H)%*%X%*%Q%*%t(X)%*%H%*%u))

  return(list(u = u, v = v))

}
