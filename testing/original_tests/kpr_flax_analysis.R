library(ggplot2)
library(glmnet)
library(hdi)

devtools::load_all()
rm(list = ls())

set.seed(42)

name.table <- read.table("R:/randolph_t/Bactocarb/temp-transfer/FlaxFX_example/genus_names.csv")

taxa.names <- unlist(lapply(as.character(name.table$V2), function(s) {
  unlist(strsplit(s, ";"))[6]
}))

data(flax)
X <- flax$X
Y <- flax$Y
Q <- flax$Q
H <- flax$H

colnames(X) <- taxa.names

n <- nrow(X)
p <- ncol(X)

# Q <- diag(p)
# H <- diag(n)


kpr.out.reml <- KPR(designMatrix = X, Y = Y, Q = Q, H = H, REML = TRUE)

kpr.out.cv <- KPR(designMatrix = X, Y = Y, Q = Q, H = H, REML = FALSE)


betahat.reml <- kpr.out.reml$beta.hat
betahat.cv <- kpr.out.cv$beta.hat[,kpr.out.cv$lambda.min.index]

estimates.df <- data.frame(reml = betahat.reml,
                           cv = betahat.cv,
                           index = 1:length(betahat.reml))

g <- ggplot(data = estimates.df, mapping = aes(x = index)) +
  geom_point(mapping = aes(y = reml), col = "red") +
  geom_point(mapping = aes(y = cv), col = "blue") +
  ylab("estimates") +
  theme_classic()

g

mu <- 1

detected.cv <- which(inference(kpr.out.cv, mu = mu)[,kpr.out.cv$lambda.min.index] < 0.05)
detected.reml <- which(inference(kpr.out.reml, mu = mu) < 0.05)

ridge.proj(x = X, y = Y)$pval.corr


marginal.out <- apply(X, 2, function(x) unlist(anova(lm(Y ~ x))["Pr(>F)"]))[1,]
marginal.out[marginal.out < 0.05]
