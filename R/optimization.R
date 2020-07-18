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

  dummyID <- factor(rep(1, n))
  lmm.fit <- lme(Y.tilde~1, random=list(dummyID = pdIdent(~-1 + u %*% diag(s))))

  lambda.reml <- lmm.fit$sigma^2 / as.numeric(VarCorr(lmm.fit)[1,1])

  return(lambda.reml)
}

buildObjectiveFunction <- function(Z.p, Y.p, H, Q)
{

    q <- length(Q)
    h <- length(H)

    # the order of args in theta is as follows:
    ## c_H, sigma_1, ..., sigma_h, c_Q, alpha_1, ..., alpha_q
    function(theta)
    {
        c_H <- abs(theta[1])
        H.sum <- abs(theta[2]) * H[[1]]
        if (h > 1) for (i in 2:h) H.sum <- H.sum + abs(theta[i+1])*H[[i]]

        c_Q <- abs(theta[h+2])
        Q.sum <- abs(theta[h+3]) * Q[[1]]
        if (q > 1) for (j in 2:q) Q.sum <- Q.sum + abs(theta[h+2+j])*Q[[j]]

        Omega <- c_Q * Z.p %*% Q.sum %*% t(Z.p) + c_H * solve(H.sum)
        (t(Y.p) %*% solve(Omega) %*% Y.p + log(det(Omega)))[1,1]
    }


}

equalityConstraintFn <- function(theta, h)
{
    sigma <- theta[2:(h+1)]
    alpha <- theta[(h+3):(length(theta))]
    z1 <- sum(abs(sigma))
    z2 <- sum(abs(alpha))
    c(z1 - 1, z2 - 1)
}

findTuningParameters <- function(Z.p, Y.p, H, Q)
{
    h <- length(H)
    q <- length(Q)
    fn <- buildObjectiveFunction(Z.p, Y.p, H, Q)
    theta0 <- c(1,
                rep(1, h) / h,
                1,
                rep(1, q) / q)
    eq.fn <- function(x) {equalityConstraintFn(x, h)}

    opt.out <- constrOptim.nl(par = theta0, fn = fn, heq = eq.fn,
                         control.outer = list(trace=FALSE,
                                              NMinit = TRUE))

    list(c_H = abs(opt.out$par[1]),
         sigma = abs(opt.out$par[2:(h+1)]),
         c_Q = abs(opt.out$par[h+2]),
         alpha = abs(opt.out$par[(h+3):(h+q+2)]))

}

# we can probably do better with this
computeCoefficientEstimates <- function(Z.p, Y.p, H, Q, theta.hat)
{
    h <- length(H)
    q <- length(Q)

    H.sum <- theta.hat$sigma[1] * H[[1]]
    if (h > 1) for (i in 2:h) H.sum <- H.sum + theta.hat$sigma[i]*H[[i]]
    Q.sum <- theta.hat$alpha[1] * Q[[1]]
    if (q > 1) for (j in 2:q) Q.sum <- Q.sum + theta.hat$alpha[j]*Q[[j]]

    lambda <- theta.hat$c_H / theta.hat$c_Q

    # solve(t(Z.p) %*% ( (1 / theta.hat$c_H) * H.sum) %*% Z.p + (1 / theta.hat$c_Q) * Q.sum) %*% t(Z.p) %*% ( (1 / theta.hat$c_H) * H.sum) %*% Y.p

    Q.sum %*% solve(t(Z.p) %*% H.sum %*% Z.p %*% Q.sum + lambda * diag(ncol(Z.p))) %*% t(Z.p) %*% H.sum %*% Y.p

}
