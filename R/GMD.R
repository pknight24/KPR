#' Generalized Matrix Decomposition
#'
#' Computes the generalized matrix decomposition of a matrix X with respect to H and Q.
#'
#' @param X An n x p data matrix.
#' @param H An n x n sample-wise similarity kernel.
#' @param Q A p x p variable-wise similarity kernel.
#' @param K The number of GMD components to include in the decomposition.
#' @param fastGMD Logical, indicates whether to use the fastGMD implementation. This requires that H and Q are both positive definite.
#' @return A list of matrices involved in the decomposition.
#' @references Wang et al. Technical report.
#' @export
GMD <- function(X, H, Q, K, fastGMD = TRUE)
{
  if (fastGMD)
  { fgmd <- fastGMD(X, H, Q, K)
    U <- fgmd$U
    S <- fgmd$S
    V <- fgmd$V
    return(list(U = U, S = S, V = V))
  }
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

  tv <- sum(diag(X %*% Q %*% t(X) %*% H))

  return(list(U = U, V = V, S = D, H = H, Q = Q, totalVariation = tv))

}

get_uv = function(X, H, Q, u_0, v_0){

  u = X%*%Q%*%v_0/as.numeric(sqrt(t(v_0)%*%t(Q)%*%t(X)%*%H%*%X%*%Q%*%v_0))
  v = t(X)%*%H%*%u/as.numeric(sqrt(t(u)%*%t(H)%*%X%*%Q%*%t(X)%*%H%*%u))

  return(list(u = u, v = v))

}

fastGMD <- function(X, H, Q, K)
{
  eigen.Q <- eigen(Q)
  L.Q <- eigen.Q$vectors %*% diag(sqrt(eigen.Q$values))
  eigen.H <- eigen(H)
  L.H <- eigen.H$vectors %*% diag(sqrt(eigen.H$values))

  X.tilde <- t(L.H) %*% X %*% L.Q
  svd.X <- svd(X.tilde)
  # U.star <- diag(eigen.H$values^(-1)) %*% t(eigen.H$vectors) %*% svd.X$u[,1:K]
  # V.star <- diag(eigen.Q$values^(-1)) %*% t(eigen.Q$vectors) %*% svd.X$v[,1:K]
  U.star <- solve(L.H) %*% svd.X$u[,1:K]
  V.star <- solve(t(L.Q)) %*% svd.X$v[,1:K]
  d.star <- svd.X$d[1:K]
  return(list(U = U.star, V = V.star, S = d.star))
}
