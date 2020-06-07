library(KPR)

data(flax)

X <- flax$X
Y <- flax$Y
Q <- flax$Q

group <- sort(rep(c(0,1), 18))
w <- rpois(36, 5)
Y.grp <- Y + 2 * group
Y.grp.w <- Y + 2 * group + (-3) * w

U <- cbind(group, w)

kpr.cov.out <- KPR(designMatrix = X, covariates = U, Y = Y.grp.w, Q = Q, n.lambda = 100)
kpr.out <- KPR(designMatrix = X, Y = Y.grp.w, Q = Q, n.lambda = 100)

kpr.cov.out$eta.hat[,1:10]
