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
old.est <- KPR(X, Y = Y, REML=TRUE)

new.est.Q <- KPR(X, Y = Y, Q = Q, REML=FALSE)
old.est.Q <- KPR(X, Y = Y, Q = Q, REML=TRUE)

new.est.H <- KPR(X, Y = Y, H = H, REML=FALSE)
old.est.H <- KPR(X, Y = Y, H = H, REML=TRUE)

new.est.H.Q <-KPR(X, Y = Y, H = H, Q=Q,REML=FALSE, scale = FALSE)
old.est.H.Q <- KPR(X, Y = Y, H = H, Q=Q,REML=TRUE, scale = FALSE)


triH <- KPR(X, Y=Y, H = list(H, X %*% Q %*% t(X), diag(n)), Q = list(Q, diag(p)))
