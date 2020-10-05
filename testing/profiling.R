rm(list=ls())
library(profvis)
devtools::load_all()
 data(flax)
 attach(flax)

 n <- nrow(X)
 p <- ncol(X)
 H.eigen <- eigen(H)
 H.eigen$values[H.eigen$values < 0.00001] <- min(H.eigen$values[H.eigen$values > 0.00001])
 H <- H.eigen$vectors %*% diag(H.eigen$values) %*% t(H.eigen$vectors)

profvis({
    par1 <- findTuningParameters(X, Y, H = list(H), Q = list(Q), trick=FALSE)
    par2 <- findTuningParameters(X, Y, H = list(H), Q = list(Q), trick=TRUE)
    rem <- remlEstimation(X, Y, H, Q)
})



## now, what if we change size of the data?
library(microbenchmark)
library(R.matlab)
library(clusterGeneration)
runtimes <- as.data.frame(microbenchmark("asdf" = 1 + 1, times=1))
runtimes$n <- 1
runtimes$p <- 1

set.seed(1234)
H <- readMat("~/Desktop/bcsstk16.mat")$Problem[[2]]
Q <- readMat("~/Desktop/bcsstk17.mat")$Problem[[2]]
H <- genPositiveDefMat(dim = nrow(H))$Sigma
X <- matrix(rnorm(n = nrow(H) * nrow(Q)), nrow = nrow(H), ncol = nrow(Q))
Y <- rnorm(n = nrow(H))
ns <- c(100, 500, 1000, nrow(H))
p <- 50
load_all()
for (n in ns)
{
    print(n)
    H.sub <- H[1:n, 1:n]
    X.sub <- X[1:n, 1:p]
    Y.sub <- Y[1:n]
    runtime <- as.data.frame(
        microbenchmark(asdf = findTuningParameters(X.sub, Y.sub, H = list(H.sub),
                                                   Q = list(diag(p)), trick=TRUE),times = 1))
    runtime$n <- n
    runtime$p <- p
    runtimes <- rbind(runtimes,runtime)
}

