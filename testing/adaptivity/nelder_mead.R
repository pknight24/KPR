rm(list=ls())
devtools::load_all(".")
library(optimx)
library(mdpeer)

data(flax)
attach(flax)

n <- nrow(X)
p <- ncol(X)

# this is the objective function for a KPR model that adds a ridge term to the Q and H penalties
fn_ridgified <- function(theta)
{
  Sigma <- X %*% solve(abs(theta[1]) * solve(Q) + abs(theta[2]) * diag(p)) %*% t(X) +
    solve(abs(theta[3]) * H + abs(theta[4]) * diag(n))
  (t(Y) %*% solve(Sigma) %*% Y + log(det(Sigma)))[1,1]
}

# this is supposed to be the objective function for the standard KPR model
fn_std <- function(theta)
{
  Sigma <- X %*% solve(abs(theta[1]) * solve(Q)) %*% t(X) + solve(abs(theta[2]) * H)
  (t(Y) %*% solve(Sigma) %*% Y + log(det(Sigma)))[1,1]
}

# this should correspond to the likelihood function derived in the riPEER paper
fn_ripeer <- function(theta)
{
    svd.Q <- svd(Q)
    U <- svd.Q$u
    S <- svd.Q$d

    X.tilde <- X %*% U
    D <- theta[1] * S + theta[2] * diag(p)

    A <- (t(Y) %*% Y)[1,1] - t(Y) %*% X.tilde %*% solve(D + t(X.tilde) %*% X.tilde) %*% t(X.tilde) %*% Y
    (n * log(A) + log(det(D + t(X.tilde) %*% X.tilde)) - log(det(D)))[1,1]

}

# testing ridgified KPR
optimx(par = c(1, 1, 1, 1), fn = fn_ridgified, method="Nelder-Mead")


# testing standard KPR
optim.out <- optimx(par = c(1, 1), fn = fn_std, method="Nelder-Mead")
KPR(designMatrix = X, Y = Y, H = H, Q = Q, inference = FALSE)$lambda
beta.test <- KPR(designMatrix = X, Y = Y, H = H, Q = Q, inference = FALSE,
                 lambda = abs(optim.out$p2 / optim.out$p1))$beta.hat
beta.kpr <- KPR(designMatrix = X, Y = Y, H = H, Q = Q, inference = FALSE)$beta.hat
plot(beta.test, beta.kpr)

                                        # comparing to riPEER
optimx(par=c(1,1), fn = fn_ripeer, method = "Nelder-Mead")
riPEER(Q = Q, y = as.matrix(Y), Z = X, optim.metod = "sbplx")
