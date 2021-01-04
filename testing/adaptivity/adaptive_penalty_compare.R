rm(list=ls())
devtools::load_all(".")

data(flax)
attach(flax)

n <- nrow(X)
p <- ncol(X)
H.eigen <- eigen(H)
H.eigen$values[H.eigen$values < 0.00001] <- min(H.eigen$values[H.eigen$values > 0.00001])
H <- H.eigen$vectors %*% diag(H.eigen$values) %*% t(H.eigen$vectors)

new.est <- KPR(X, Y= Y, REML=FALSE)
new.est <- inference(new.est)
old.est <- KPR(X, Y = Y, REML=TRUE)
old.est <- inference(old.est)
colors <- rep("black", length(new.est$beta.hat))
colors[new.est$p.values < 0.05] <- "red"
plot(new.est$beta.hat, pch=19,
     col = colors)

new.est.Q <- KPR(X, Y = Y, Q = Q, REML=FALSE)
new.est.Q <- inference(new.est.Q)
old.est.Q <- KPR(X, Y = Y, Q = Q, REML=TRUE)
old.est.Q <- inference(old.est.Q)
colors <- rep("black", length(new.est.Q$beta.hat))
colors[new.est.Q$p.values < 0.05] <- "red"
plot(new.est.Q$beta.hat, pch=19,
     col = colors)

new.est.H <- KPR(X, Y = Y, H = H, REML=FALSE)
new.est.H <- inference(new.est.H)
old.est.H <- KPR(X, Y = Y, H = H, REML=TRUE)
old.est.H <- inference(old.est.H)
colors <- rep("black", length(new.est.H$beta.hat))
colors[new.est.H$p.values < 0.05] <- "red"
plot(new.est.H$beta.hat, pch=19,
     col = colors)

new.est.H.Q <-KPR(X, Y = Y, H = H, Q=Q,REML=FALSE, scale = FALSE)
new.est.H.Q <- inference(new.est.H.Q)
old.est.H.Q <- KPR(X, Y = Y, H = H, Q=Q,REML=TRUE, scale = FALSE)
old.est.H.Q <- inference(old.est.H.Q)
colors <- rep("black", length(new.est.H.Q$beta.hat))
colors[new.est.H.Q$p.values < 0.05] <- "red"
plot(new.est.H.Q$beta.hat, pch=19,
     col = colors)


# strange inference results with this one
triH <- KPR(X, Y=Y, H = list(H, X %*% Q %*% t(X), diag(n)), Q = list(diag(p)))
triH <- inference(triH)
colors <- rep("black", length(triH$beta.hat))
colors[triH$p.values < 0.05] <- "red"
plot(triH$beta.hat, pch=19,
     col = colors)
