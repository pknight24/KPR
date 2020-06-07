library(ggplot2)
library(glmnet)

devtools::load_all()
rm(list = ls())

set.seed(42)

data(flax)
X <- flax$X
Q <- flax$Q
H <- flax$H

n <- nrow(X)
p <- ncol(X)

vectors.Q <- eigen(Q)$vectors
values.Q <- eigen(Q)$values

beta.true <- sqrt(values.Q)
beta.true[abs(beta.true) < quantile(abs(beta.true), probs = 0.75)] <- 0

group = sort(rep(c(1,0), n/2))
E <- matrix(group)
eta.true <- 5

Y_ <- X %*% beta.true + rnorm(n, sd=0.05)
Y <- Y_ - mean(Y_)


# Q <- diag(p)
# H <- diag(n)

kpr.out <- KPR(designMatrix = X, Y = Y, Q = Q, H = H, REML = TRUE)

betahat.reml <- kpr.out$beta.hat

# riPEER.out <- mdpeer::riPEER(Q = Q, Z = X, y = Y, compute.boot.CI = TRUE)

# rip.ci <- riPEER.out$boot.CI
# rip.ci$signif <- sapply(1:length(rip.ci$lower), function(s)
#   {
#   !dplyr::between(0, rip.ci$lower[s],rip.ci$upper[s])
# })

true <- which(beta.true > 0)
detected.reml <- which(inference(kpr.out) < 0.05)
# detected.rip <- which(rip.ci$signif)
