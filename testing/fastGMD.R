fastGMD <- function(X, H, Q)
{
  eigen.Q <- eigen(Q)
  L.Q <- eigen.Q$vectors %*% diag(sqrt(eigen.Q$values))
  eigen.H <- eigen(H)
  L.H <- eigen.H$vectors %*% diag(sqrt(eigen.H$values))

  X.tilde <- t(L.H) %*% X %*% L.Q
  svd.X <- svd(X.tilde)
  U.star <- solve(L.H) %*% svd.X$u
  V.star <- solve(L.Q) %*% svd.X$v
  d.star <- svd.X$d
  return(list(U = U.star, V = V.star, D = d.star))
}

library(KPR)

data("yatsunenko")

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

Q <- generateSimilarityKernel(patristic)
gmd <- GMD(X = Z, H = diag(n), Q = Q, K = length(svd(Z)$d))
fgmd <- fastGMD(Z, diag(n), Q)
