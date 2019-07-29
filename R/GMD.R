#' Generalized Matrix Decomposition
#'
#' Computes the generalized matrix decomposition of a matrix X with respect to H and Q.
#'
#' @param X An n x p data matrix.
#' @param H An n x n sample-wise similarity kernel.
#' @param Q A p x p variable-wise similarity kernel.
#' @param K The number of GMD components to include in the decomposition.
#' @return A list of matrices involved in the decomposition.
#' @export
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

  tv <- sum(diag(X %*% Q %*% t(X) %*% H))

  return(list(U = U, V = V, S = D, H = H, Q = Q, totalVariation = tv))

}

get_uv = function(X, H, Q, u_0, v_0){

  u = X%*%Q%*%v_0/as.numeric(sqrt(t(v_0)%*%t(Q)%*%t(X)%*%H%*%X%*%Q%*%v_0))
  v = t(X)%*%H%*%u/as.numeric(sqrt(t(u)%*%t(H)%*%X%*%Q%*%t(X)%*%H%*%u))

  return(list(u = u, v = v))

}

