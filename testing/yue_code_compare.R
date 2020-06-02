rm(list = ls())
devtools::load_all(".")
source("../GMDR_new.R")
data(yatsunenko)

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


# First, I will compare a Ridge model (no H or Q)
yue.fit <- KPR.est(Y, Z, diag(n), diag(p), scale = FALSE)
pkg.fit <- KPR(designMatrix = Z, Y = Y, scale = FALSE, lambda = yue.fit$lambda.opt)
plot(pkg.fit$beta.hat, col = "red")
points(yue.fit$beta.opt, col = "blue")

plot(yue.fit$beta.opt, col="black")
points(yue.fit$betahat2, col = "green")
