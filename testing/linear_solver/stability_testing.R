rm(list = ls())
devtools::load_all(".")
library(ggplot2)
library(dplyr)
data(flax)
attach(flax)

lambda <- seq(from = 0.001, to = 1, length.out = 100)

b_linear <- sapply(lambda, FUN = function(s) KPR(designMatrix = X, Y = Y, H = H, Q = Q, lambda = s, inference=FALSE)$beta.hat)
b_old <- sapply(lambda, FUN = function(s) KPR(designMatrix = X, Y = Y, H = H, Q = Q, lambda = s, inference=FALSE, linear_solve=FALSE)$beta.hat)
beta_variances_linear <- apply(b_linear, 2, var)
beta_variances_old <- apply(b_old, 2, var)
qplot(beta_variances_linear, beta_variances_old) + theme_classic()

eigen.Q <- eigen(Q)
Q.pert <- eigen.Q$vectors %*% diag(eigen.Q$values / 10000) %*% t(eigen.Q$vectors) # shrink the eigenvalues of Q

par(mfrow= c(1, 2))
plot(KPR(designMatrix = X, Y =Y, H = H, Q = Q, inference=FALSE)$beta.hat,
     KPR(designMatrix = X, Y =Y, H = H, Q = Q.pert,inference=FALSE)$beta.hat)
plot(KPR(designMatrix = X, Y =Y, H = H, Q = Q, linear_solve=FALSE,inference=FALSE)$beta.hat,
     KPR(designMatrix = X, Y =Y, H = H, Q = Q.pert, linear_solve = FALSE,inference=FALSE)$beta.hat)

I <- diag(nrow = nrow(Q), ncol=ncol(Q))
b_linear <- sapply(lambda, FUN = function(s) KPR(designMatrix = X, Y = Y, H = H, Q = s * I, inference=FALSE)$beta.hat)
b_old <- sapply(lambda, FUN = function(s) KPR(designMatrix = X, Y = Y, H = H, Q = s * I, inference=FALSE, linear_solve=FALSE)$beta.hat)
beta_variances_linear <- apply(b_linear, 2, var)
beta_variances_old <- apply(b_old, 2, var)
qplot(beta_variances_linear) + theme_classic()

