library(mdpeer)
library(MASS)
library(hdi)
library(glmnet)
library(KPR)
rm(list = ls())
source("~/The-GMD-biplot/GMDR_new.R")

# set.seed(1234)


data("riboflavin")
X.orig <- scale(riboflavin$x)
Y <- scale(riboflavin$y)

cv.fit <- cv.glmnet(X.orig, Y)
signal.param.indices <- coef(cv.fit,s= "lambda.1se")@i[2:length(coef(cv.fit,s= "lambda.1se")@i)]
all.param.indices <- sample(setdiff(1:ncol(X.orig), signal.param.indices), 100)
X <- X.orig[,all.param.indices]

n <- nrow(X)
p <- ncol(X)

K <- 10 # K folds cross validation to evaluate the models
errors <- matrix(rep(0, K * 4), K, 4) # a K x 4 matrix to keep track of the mean squared errors
colnames(errors) <- c("KPR", "MASS::ridge", "riPEER", "GMD")

randidx <- sample(1:n, n)
Yrand <- Y[randidx]
Xrand <- X[randidx, ]

for (k in 1:K)
{
  Ytrain <- Yrand[-((n / K * (k - 1) + 1):(n / K * k))]
  Ytest <- Yrand[((n / K * (k - 1) + 1):(n / K * k))]
  Xtrain <- Xrand[-((n / K * (k - 1) + 1):(n / K * k)), ]
  Xtest <- Xrand[((n / K * (k - 1) + 1):(n / K * k)), ]

  kpr.fit <- KPR(Xtrain, Ytrain, K = 5)
  mass.ridge.fit <- lm.ridge(Ytrain ~ Xtrain, lambda = kpr.fit$lambda.1se)
  riPEER.fit <- riPEER(Q = diag(ncol(Xtrain)), y = as.matrix(Ytrain), Z = Xtrain)
  GMD.fit <- GMDR.est(Y = Ytrain, X = Xtrain, H = diag(nrow(Xtrain)), Q = diag(ncol(Xtrain)))

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
