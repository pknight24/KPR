library(MASS)
library(mdpeer)
library(KPR)
rm(list = ls())
source("~/The-GMD-biplot/GMDR_new.R")

data(flax)

X <- flax$X
Y <- as.vector(scale(flax$Y))
H <- flax$H
Q <- flax$Q

# H <- diag(nrow(X))
# Q <- diag(ncol(X))

n <- nrow(X)
p <- ncol(X)

K <- n # K folds cross validation to evaluate the models
errors <- matrix(rep(0, K * 4), K, 4) # a K x 4 matrix to keep track of the mean squared errors
colnames(errors) <- c("KPR", "MASS::ridge", "riPEER", "GMD")

randidx <- sample(1:n, n)
Yrand <- Y[randidx]
Xrand <- X[randidx, ]
Hrand <- H[randidx, randidx]

for (k in 1:K)
{
  Ytrain <- Yrand[-((n / K * (k - 1) + 1):(n / K * k))]
  Ytest <- Yrand[((n / K * (k - 1) + 1):(n / K * k))]
  Xtrain <- Xrand[-((n / K * (k - 1) + 1):(n / K * k)), ]
  Xtest <- Xrand[((n / K * (k - 1) + 1):(n / K * k)), ]
  Htrain <- Hrand[-((n / K * (k - 1) + 1):(n / K * k)),-((n / K * (k - 1) + 1):(n / K * k))]

  kpr.fit <- KPR(Xtrain, Ytrain, Q = Q, H = Htrain, K = 5)
  mass.ridge.fit <- lm.ridge(Ytrain ~ Xtrain, lambda = kpr.fit$lambda.1se)
  riPEER.fit <- riPEER(Q = Q, y = as.matrix(Ytrain), Z = Xtrain)
  GMD.fit <- GMDR.est(Y = Ytrain, X = Xtrain, H = Htrain, Q = Q)

  kpr.b <- kpr.fit$betahat[, kpr.fit$lambda.1se.index]
  mass.b <- mass.ridge.fit$coef
  ripeer.b <- riPEER.fit$b.est
  names(ripeer.b) <- colnames(X)
  GMD.b <- GMD.fit$beta.opt
  names(GMD.b) <- colnames(X)

  errors[k,1] <- mean((Ytest - Xtest %*% kpr.b)^2)
  errors[k,2] <- mean((Ytest - Xtest %*% mass.b)^2)
  errors[k,3] <- mean((Ytest - Xtest %*% ripeer.b)^2)
  errors[k,4] <- mean((Ytest - Xtest %*% GMD.b)^2)

}

print(errors)
print(colSums(errors))
