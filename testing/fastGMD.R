devtools::load_all(".")
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
gmd <- GMD(X = Z, H = diag(n), Q = Q, K = length(svd(Z)$d), fastGMD = FALSE)
fgmd <- GMD(X = Z, H = diag(n), Q = Q, K = length(svd(Z)$d), fastGMD = TRUE)
