devtools::load_all(".")
library(dplyr)
rm(list = ls())

###### Load data ######
data("yatsunenko")

set.seed(1234)

counts <- yatsunenko$raw.counts
age <- yatsunenko$age
unifrac <- yatsunenko$unifrac
patristic <- yatsunenko$patristic
ec <- yatsunenko$ec
geo <- yatsunenko$geography

rm(yatsunenko)

n <- nrow(counts)
p <- ncol(counts)

counts.clr <- t(apply(log(counts + 1), 1, function(x) x - mean(x)))
Z <- apply(counts.clr, 2, function(x) x - mean(x))
Y <- log(age) - mean(log(age)) # centered log age
E <- model.matrix(~ geo)[,-1]


permuteEigenvalue <- function(mat, idx)
{
  eigen.mat <- eigen(mat)
  permuted.evalues <- eigen.mat$values
  temp <- permuted.evalues[idx] + 
    rnorm(n = 1, mean = permuted.evalues[idx], sd = 5)
  permuted.evalues[idx] <- ifelse(temp <= 0, permuted.evalues[idx], temp)
  eigen.mat$vectors %*% diag(permuted.evalues) %*% t(eigen.mat$vectors)
}

addNoise <- function(mat)
{
  mat.data <- as.vector(mat[lower.tri(mat)])
  noisy.data <- sapply(mat.data, FUN=function(x) rnorm(n = 1, mean = x, sd = abs(x/4)))
  noisy.mat <- matrix(nrow = nrow(mat), ncol = ncol(mat))
  noisy.mat[lower.tri(noisy.mat)] <- noisy.data
  noisy.mat <- t(noisy.mat)
  noisy.mat[lower.tri(noisy.mat)] <- noisy.data
  diag(noisy.mat) <- diag(mat)
  return(noisy.mat)
}

compareModels <- function(orig, perm, ...)
{
  orig.beta <- as.vector(orig$beta.hat)  
  perm.beta <- as.vector(perm$beta.hat)
  relative.error <- round(norm(perm.beta - orig.beta, "2") / norm(orig.beta, "2"), digits=5)
  
  orig.p <- as.vector(orig$p.values)
  perm.p <- as.vector(perm$p.values)
  
  par(mfrow=c(2,1))
  
  plot(orig.beta, col="red", pch=1, ylab = "Coefficients", xlab = "", ...)
  points(perm.beta, col="blue", pch=2)
  legend("topleft", legend = c("original", "permuted"), col = c("red", "blue"), pch = c(1, 2))
  text(x=length(orig.beta) / 2, y = max(orig.beta) * 0.7, labels = paste("Relative error:", relative.error))
  
  plot(orig.p, col="red", pch=1, ylab="pvalues", xlab="")
  points(perm.p, col="blue", pch=2)
  
  
}
