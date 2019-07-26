# library(KPR)
devtools::load_all()
rm(list = ls())


n <- 20
p <- 20
signif.params <- sample(x = p, size = 10)
signif.coefs <- rpois(n = 10,lambda = 20)

X <- matrix(rnorm(n * p), n, p)


Y <- X[,signif.params] %*% signif.params + rnorm(n, sd = 0.75)

# kpr.out <- KPR(designMatrix = X, Y = Y, K = 5, lambda = 1)
#
# infer.out <- inference(kpr.out)

K <- 2
microbenchmark::microbenchmark(cpp=KPR(designMatrix = X, Y = Y, useCpp = TRUE, K = K, seed=42),
                               r=KPR(designMatrix = X, Y = Y,  useCpp = FALSE, K = K, seed=42))
