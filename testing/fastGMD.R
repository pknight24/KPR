fastGMD <- function(X, H, Q)
{
  eigen.Q <- eigen(Q)
  L.Q <- eigen.Q$vectors %*% diag(sqrt(eigen.Q$values))
  eigen.H <- eigen(H)
  L.H <- eigen.H$vectors %*% diag(sqrt(eigen.H$values))

  X.tilde <- t(L.H) %*% X %*% L.Q
  svd.X <- svd(X.tilde)
  U.star <- solve(L.H) %*% svd.X$u
  V.star <- solve(L.Q) %*% svd.X$v
  d.star <- svd.X$d
  return(list(U = U.star, V = V.star, D = d.star))
}
