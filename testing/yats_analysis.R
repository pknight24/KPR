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

counts.clr <- t(apply(log(counts + 1), 1, function(x) x - mean(x)))
Z <- apply(counts.clr, 2, function(x) x - mean(x))
Y <- log(age) - mean(log(age)) # centered log age

Q <- generateSimilarityKernel(patristic)


model.fit.Q <- KPR(designMatrix = Z, Y = Y, Q = Q)
model.fit <- KPR(designMatrix = Z, Y = Y)

plot(model.fit.Q$beta.hat, col = ifelse(model.fit.Q$p.values < 0.05,
                                        "blue", "black"))
abline(a = 0, b = 0)
