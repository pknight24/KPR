remlEstimation <- function(Z, E, Y, H, Q, cov.missing)
{

  n <- nrow(Z)
  p <- ncol(Z)

  if (cov.missing) P <- diag(n)
  else P <- diag(n) - E %*% solve(t(E) %*% H %*% E) %*% t(E) %*% H
  Y.p <- P %*% Y
  Z.p <- P %*% Z # apply P to the variables that should be penalized

  L.H <- t(chol(H))
  L.Q <- t(chol(Q))
  Z.tilde <- t(L.H) %*% Z.p %*% L.Q
  Y.tilde <- t(L.H) %*% Y.p

  u <- svd(Z.tilde)$u
  s <- svd(Z.tilde)$d

  dummyID <- factor(rep(1, n))
  lmm.fit <- lme(Y.tilde~1, random=list(dummyID = pdIdent(~-1 + u %*% diag(s))))

  lambda.reml <- lmm.fit$sigma^2 / as.numeric(VarCorr(lmm.fit)[1,1])

  return(lambda.reml)
}