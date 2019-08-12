library(KPR)
library(dplyr)
rm(list = ls())

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
E <- model.matrix(~ geo)[,-1]

alpha <- 0.01

Q <- (1 - alpha) * generateSimilarityKernel(patristic) + alpha * diag(p)

# K.ec <- solve(ec %*% t(ec)) %>%
#   (function(x) x / svd(x)$d[1])

model.fit <- KPR(designMatrix = Z, Y = Y, Q = Q)


plot(model.fit$beta.hat, col = ifelse(model.fit$p.values < 0.05,
                                        "blue", "grey"),
     ylab="Effect size", main = "", xlab = "",
     pch = 19)
abline(a = 0, b = 0, col = "red")

