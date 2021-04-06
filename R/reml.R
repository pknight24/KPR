remlEstimation <- function(Z.p, Y.p, H, Q)
{

  n <- nrow(Z.p)
  p <- ncol(Z.p)

    L.H <- t(chol(H))
  L.Q <- t(chol(Q))
  Z.tilde <- t(L.H) %*% Z.p %*% L.Q
  Y.tilde <- t(L.H) %*% Y.p

  u <- svd(Z.tilde)$u
  s <- svd(Z.tilde)$d

  us <- u %*% diag(s)

  dummyID <- factor(rep(1, n))
  lmm.fit <- lme(Y.tilde~1, random=list(dummyID = pdIdent(~-1 + us)))

  lambda.reml <- lmm.fit$sigma^2 / as.numeric(VarCorr(lmm.fit)[1,1])

  return(lambda.reml)
}
