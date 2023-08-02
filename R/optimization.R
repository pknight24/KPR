buildObjectiveFunction <- function(Z.p, Y.p, H, Q, Q.inv)
{

    q <- length(Q)
    h <- length(H)

    # the order of args in theta is as follows:
    ## c_H, sigma_1, ..., sigma_h, c_Q, alpha_1, ..., alpha_q
    function(theta)
    {
        c_H <- abs(theta[1])
        if (h > 1) {
            H.sum <- abs(theta[2]) * H[[1]]
            for (i in 2:h) H.sum <- H.sum + abs(theta[i+1])*H[[i]]
        }
        else H.sum <- H[[1]]

        k <- h > 1
        c_Q <- abs(theta[h+1+k])
        if (q > 1) {
            Q.sum <- abs(theta[h + 2 + k]) * Q[[1]]
            for (j in 2:q) Q.sum <- Q.sum + abs(theta[h+1+j+k]) * Q[[j]]
        }
        else Q.sum <- Q[[1]]

        if (!Q.inv) Q.sum <- solve(Q.sum)

        Omega <- c_Q * Z.p%*% Q.sum %*% t(Z.p) + c_H * solve(H.sum)
        Omega.inv <- solve(Omega)
        output <-  t(Y.p) %*% Omega.inv %*% Y.p + fastLogDet(Omega)
        return(output)
    }


}

equalityConstraintFn <- function(theta, h, q)
{
    if (h > 1) {
        sigma <- theta[2:(h+1)]
        z1 <- sum(abs(sigma)) - 1
    } else z1 <- NULL
    if (q > 1) {
        alpha <- theta[(h+2 + (h>1)):(length(theta))]
        z2 <- sum(abs(alpha)) - 1
    } else z2 <- NULL
    if (is.null(z1) & !is.null(z2)) return(z2)
    if (!is.null(z1) & is.null(z2)) return(z1)
    c(z1, z2)
}

findTuningParameters <- function(Z.p, Y.p, H, Q, Q.inv,
                                 control.outer,
                                 control.optim)
{
    h <- length(H)
    q <- length(Q)
    fn_raw <- buildObjectiveFunction(Z.p, Y.p, H, Q)
    fn <- function(x) as.numeric(fn_raw(x))
    h.theta0 <- 1
    if (h > 1) h.theta0 <- c(h.theta0, rep(1,h) / h)
    q.theta0 <- 1
    if (q > 1) q.theta0 <- c(q.theta0, rep(1,q) / q)
    theta0 <- c(h.theta0, q.theta0)
    eq.fn <- function(x) {equalityConstraintFn(x, h, q)}
    if (h == 1 & q == 1) eq.fn <- NULL
    if (is.null(eq.fn)) opt.out <- optim(par = theta0, fn = fn, method="BFGS")
    else {
        opt.out <- constrOptim.nl(par = theta0, fn = fn,heq = eq.fn,
                                        control.outer = control.outer,
                                        control.optim = control.optim)
    }
    k <- h > 1
    sigma <- ifelse(k, abs(opt.out$par[(2:(h+1))]), 1)
    if (k) sigma <- abs(opt.out$par[(2:(h+1))])
    else sigma <- 1
    if (q > 1) alpha <- abs(opt.out$par[((h+2+k):(h+q+1+k))])
    else alpha <- 1
    list(c_H = abs(opt.out$par[1]),
         sigma = sigma,
         c_Q = abs(opt.out$par[h+1+k]),
         alpha = alpha)

}


computeCoefficientEstimates <- function(Z.p, Y.p, H, Q, theta.hat, Q.inv)
{
    h <- length(H)
    q <- length(Q)

    H.sum <- theta.hat$sigma[1] * H[[1]]
    if (h > 1) for (i in 2:h) H.sum <- H.sum + theta.hat$sigma[i]*H[[i]]
    Q.sum <- theta.hat$alpha[1] * Q[[1]]
    if (q > 1) for (j in 2:q) Q.sum <- Q.sum + theta.hat$alpha[j]*Q[[j]]

    lambda <- theta.hat$c_H / theta.hat$c_Q

    # solve(t(Z.p) %*% ( (1 / theta.hat$c_H) * H.sum) %*% Z.p + (1 /
    # theta.hat$c_Q) * Q.sum) %*% t(Z.p) %*% ( (1 / theta.hat$c_H) * H.sum) %*%
    # Y.p
    
    if (!Q.inv) Q.sum <- solve(Q.sum)

    Q.sum %*% solve(t(Z.p) %*% H.sum %*% Z.p %*% Q.sum + lambda * diag(ncol(Z.p))) %*% t(Z.p) %*% H.sum %*% Y.p

}


fastLogDet <- function(A) return(2 * sum(log(diag(chol(A)))))
