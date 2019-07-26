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

beta.true <- vectors.Q[,1]
beta.true[abs(beta.true) < quantile(abs(beta.true), probs = 0.75)] <- 0

group = sort(rep(c(1,0), n/2))
E <- matrix(group)
eta.true <- 5

Y_ <- X %*% beta.true + rnorm(n, sd=0.05)
Y <- Y_ - mean(Y_)


# Q <- diag(p)
# H <- diag(n)

kpr.out.reml <- KPR(designMatrix = X, Y = Y, Q = Q, H = H, REML = TRUE)

kpr.out.cv <- KPR(designMatrix = X, Y = Y, Q = Q, H = H, REML = FALSE)


betahat.reml <- kpr.out.reml$beta.hat
betahat.cv <- kpr.out.cv$beta.hat[,kpr.out.cv$lambda.min.index]

estimates.df <- data.frame(true = beta.true, reml = betahat.reml,
                           cv = betahat.cv,
                           index = 1:length(beta.true))

g <- ggplot(data = estimates.df, mapping = aes(x = index)) +
  geom_point(mapping = aes(y =  true), col="black") +
  geom_point(mapping = aes(y = reml), col = "red") +
  geom_point(mapping = aes(y = cv), col = "blue") +
  ylab("estimates") +
  theme_classic()

g
# riPEER.out <- mdpeer::riPEER(Q = Q, Z = X, y = Y, compute.boot.CI = TRUE)

# mass.out <- MASS::lm.ridge(Y ~ X, lambda = kpr.out.reml$lambda)

# rip.ci <- riPEER.out$boot.CI
# rip.ci$signif <- sapply(1:length(rip.ci$lower), function(s)
#   {
#   !dplyr::between(0, rip.ci$lower[s],rip.ci$upper[s])
# })

true <- which(beta.true > 0)
detected.cv <- which(inference(kpr.out.cv)[,kpr.out.cv$lambda.min.index] < 0.05)
detected.reml <- which(inference(kpr.out.reml) < 0.05)
# detected.rip <- which(rip.ci$signif)
