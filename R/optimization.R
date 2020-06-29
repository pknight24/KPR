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

# H and Q are lists of matrices
#' @export
buildObjectiveFunction <- function(Z.p, Y.p, H, Q)
{

    Q.inv <- lapply(Q, solve) # invert all the Q matrices before building the function
    q <- length(Q)
    h <- length(H)
    function(theta)
    {

        H.sum <- abs(theta[1]) * H[[1]]
        if (h > 1) for (i in 2:h) H.sum <- H.sum + abs(theta[i])*H[[i]]
        Q.sum <- abs(theta[h+1]) * Q.inv[[1]]
        if (q > 1) for (j in 2:q) Q.sum <- Q.sum + abs(theta[h+j])*Q.inv[[j]]
        Omega <- Z.p %*% solve(Q.sum) %*% t(Z.p) + solve(H.sum)
        (t(Y.p) %*% solve(Omega) %*% Y.p + log(det(Omega)))[1,1]
    }


}

#' @export
fitParameters <- function(Z.p, Y.p, H, Q)
{
    if (!is.list(Q)) Q <- list(Q)
    if (!is.list(H)) H <- list(H)


    fn <- buildObjectiveFunction(Z.p, Y.p, H, Q)
    optim.out <- optimx(par = rep(1, length(Q) + length(H)), fn = fn, method = "Nelder-Mead")
    sapply(1:(length(Q) + length(H)), function(i) optim.out[[paste0("p",i)]])

}

computeBetaHat <- function(Z.p, Y.p, H, Q, theta.hat)
{
    if (!is.list(Q)) Q <- list(Q)
    if (!is.list(H)) H <- list(H)
    q <- length(Q)
    h <- length(H)
    Q.inv <- lapply(Q, solve)

    H.sum <- abs(theta.hat[1]) * H[[1]]
    if (h > 1) for (i in 2:h) H.sum <- H.sum + abs(theta.hat[i])*H[[i]]
    Q.sum <- abs(theta.hat[h+1]) * Q.inv[[1]]
    if (q > 1) for (j in 2:q) Q.sum <- Q.sum + abs(theta.hat[h+j])*Q.inv[[j]]


    solve(t(Z.p) %*% H.sum %*% Z.p + Q.sum) %*% t(Z.p) %*% H.sum %*% Y.p

}
