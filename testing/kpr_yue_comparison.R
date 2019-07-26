library(KPR)
library(glmnet)
rm(list = ls())
data(flax)
X <- flax$X
Y <- flax$Y
Q <- diag(ncol(X))
H <- diag(nrow(X))

Q <- flax$Q
H <- flax$H

L.H <- chol(H)
L.Q <- chol(Q)

X.KPR <- L.H %*% X %*% L.Q
Y.KPR <- L.H %*% Y

KPR.cv <- cv.glmnet(X.KPR, Y.KPR, family = "gaussian", alpha = 0, standardize = FALSE, intercept = FALSE)

K <- solve(t(X) %*% H %*% X + KPR.cv$lambda.min/2 * solve(Q))
yue.betahat <- c(K %*% t(X) %*% H %*% Y)
print(yue.betahat)

kpr.package.fit <- KPR(designMatrix = X, Y = Y, Q = Q, H = H)
print(kpr.package.fit$beta.hat)

plot(yue.betahat, type="l")
plot(kpr.package.fit$beta.hat[,kpr.package.fit$lambda.min.index], type = "l")

