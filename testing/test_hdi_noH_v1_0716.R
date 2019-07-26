args = (commandArgs(TRUE))
for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
}

info.det = as.numeric(a);
s = as.numeric(b)
#mis = as.numeric(d)
#eigen.id = as.numeric(d);


# simulation study: type-I control: based on Yatsunenko data
library(glmnet)
library(MASS)
library(Matrix)

library(scalreg)
library(hdi)
library(natural)
library(Grace)

rm(list = ls())

#source('/Users/taryue/Dropbox/SGMD/code/SGMD_test.R')
#source('/Users/taryue/Dropbox/SGMD/code/SGMD_test.R')
#setwd('/Users/taryue/Dropbox/SGMD/code/new')
source('testing/GMDR_new.R')
source('dScore.R')
#source('gram_schmidt.R')
#source('/Users/taryue/Documents/GMD biplots/code/gram_schmidt.R')

n = 200; p = 300

AR.mat = function(p,rho){
   mat<- matrix(rho, p, p)
   diag(mat) <- 1
   return(mat)
}

autocorr.mat <- function(p, rho) {
  mat <- diag(p)
  return(rho^abs(row(mat)-col(mat)))
}

Q = as.matrix(bdiag(AR.mat(p/2, 0.9), AR.mat(p/2, 0.5)))
#Q = 1/info.det*Q

eigen.Q = eigen(Q)
values.Q = eigen.Q$values
vectors.Q = eigen.Q$vectors
values.Q.new = values.Q/values.Q[1]
Q.new = vectors.Q%*%diag(values.Q.new)%*%t(vectors.Q)
Q.new.inv = vectors.Q%*%diag(1/values.Q.new)%*%t(vectors.Q)

#plot(1:100, vectors.Q[,51], type = 'l')

#Q.sqrt = vectors.Q%*%diag(sqrt(values.Q))%*%t(vectors.Q)
#Q.sqrt.n = vectors.Q%*%diag(1/sqrt(values.Q))%*%t(vectors.Q)
L.Q = vectors.Q%*%diag(sqrt(values.Q.new))
L.Q.inv = diag(1/sqrt(values.Q.new))%*%t(vectors.Q)

#H = autocorr.mat(n, 0.3)
#eigen.H = eigen(H)
#B = eigen.H$vectors
#H.sqrt = B%*%diag(sqrt(eigen.H$values))%*%t(B)
#H.sqrt.n = B%*%diag(1/sqrt(eigen.H$values))%*%t(B)
set.seed(1992)
X.raw = matrix(rnorm(n*p), n, p)
X.tilde = scale(X.raw, center = T, scale = F)
sum.tilde = colSums(X.tilde^2)
X.tilde = X.tilde%*%diag(sqrt(n/sum.tilde))
X = X.tilde%*%t(vectors.Q)

#if(s == 1){
  write.table(X, file = paste0('../simu_data/X_0715_v1_matlab_', s), row.names = F, col.names = F, quote = F)
#}


# X: column mean 0; column squared sum to n

#----------------------------------------------
# power
#(2): simulation data: control type-I error: beta.star = 0
#s = 1

#for(info.det in c(0.25, 0.4, 0.6, 1)){

  #for(s in 1:100){


    error = 0.5*rnorm(n)
    epsilon = error
    #beta.Q = rep(0,p)
    #beta.Q[c(1,2,3,4,5)] = 5

    #info.det = 1.8
    beta.star = vectors.Q[,1]*info.det
    #plot(1:p, beta.star, type = 'l')

    Y.star = X%*%beta.star
    #Y.star.tilde = H.sqrt%*%Y.star
    Y = Y.star + epsilon

    var(Y.star)/var(Y)

    write.table(Y, file = paste0('../simu_data/Y_0715_v1_info_',info.det,'_',s), row.names = F, col.names = F, quote = F)


  #}



#}


#Y.tilde = H.sqrt%*%Y
# Rsquare
 # 0.30



#(3): compute p value for GMDR
gmdr.est = GMDR.est(Y, X, diag(rep(1,n)), Q.new)
gmdr.est$set.opt

#lines(1:p, gmdr.est$GMD.x$V[,1], col = 'red')
lines(1:p, gmdr.est$beta.opt, col = 'red', type = 'l')
p.gmdr.w = GMDR.test(gmdr.est, 1, 0.05, weight = TRUE)
which(p.gmdr.w <= 0.05)

p.gmdr.now = GMDR.test(gmdr.est, 1, 0.05, weight = FALSE)
which(p.gmdr.now <= 0.05)


#(3.5) grace: same estimator, difference inference procedure
sum.x = colSums(X^2)
X.grace = X%*%diag(sqrt(n/sum.x))
test.grace = grace.test(Y = as.numeric(Y), X = X, L = Q.new.inv, lambda.L = seq(0.1,1,0.1))
p.grace = test.grace$pvalue

#fit.grace = grace.est(Y, X.grace, Q)
#lines(1:p, fit.grace$beta.opt, col = 'blue')
#p.grace = grace.test(fit.grace, 0.05)
which(p.grace <= 0.05)


#(3.6) KPR:
X.KPR = X%*%L.Q
KPR.cv = cv.glmnet(X.KPR, Y, family = "gaussian", alpha = 0, standardize = F, intercept = F)
Hat <- solve(t(X) %*% X +KPR.cv$lambda.min/2*Q.new.inv)
betahat <- c(Hat %*% t(X) %*% Y)
#plot(1:100, betahat, type = 'l')

KPR.est = list(GMD.x = gmdr.est$GMD.x, beta.opt = betahat, lambda.opt = KPR.cv$lambda.min/2)
p.kpr.w = KPR.test(KPR.est, 1, 0.05, weight = TRUE)
which(p.kpr.w <= 0.05)

p.kpr.now = KPR.test(KPR.est, 1, 0.05, weight = FALSE)
which(p.kpr.now <= 0.05)


#(4): compute p value for simple ridge in two ways, ignore Q for both estimation and inference
X.rel.scale = X.grace
# using HDI to calculate p value for ridge
outridge = ridge.proj(x = X.rel.scale, y = Y)
p.ridge = outridge$pval
which(p.ridge <= 0.05) # more conservative

#(5): LDPE (zhang and zhang, 2014)
outlasso = lasso.proj(x = X.rel.scale, y = Y)
p.lasso = outlasso$pval
which(p.lasso <= 0.05)

#(6): Decorrelated score
p.dscore <- rep(0,p)
for(j in 1:p){

  outdscore <- dScore(as.numeric(Y), X, coi = j)
  p.dscore[j] <- outdscore$pval
}



write.table(t(c(s,as.numeric(p.gmdr.w))), file = paste0('../result/GMDR_weight_hdi_noH_v1_0716_',info.det,'.txt'),row.names = F, col.names = F, quote = F, append = TRUE)
write.table(t(c(s,as.numeric(p.gmdr.now))), file = paste0('../result/GMDR_noweight_hdi_noH_v1_0716_',info.det,'.txt'),row.names = F, col.names = F, quote = F, append = TRUE)
#write.table(t(c(s,as.numeric(p.gmdr.half))), file = paste0('../result/GMDR_half_hdi_noH_v1_0716_',info.det,'.txt'),row.names = F, col.names = F, quote = F, append = TRUE)

write.table(t(c(s,as.numeric(p.kpr.w))), file = paste0('../result/KPR_weight_hdi_noH_v1_0716_',info.det,'.txt'),row.names = F, col.names = F, quote = F, append = TRUE)
write.table(t(c(s,as.numeric(p.kpr.now))), file = paste0('../result/KPR_noweight_hdi_noH_v1_0716_',info.det,'.txt'),row.names = F, col.names = F, quote = F, append = TRUE)
#write.table(t(c(s,as.numeric(p.KPR.half))), file = paste0('../result/KPR_half_hdi_noH_v1_0716_',info.det,'.txt'),row.names = F, col.names = F, quote = F, append = TRUE)

write.table(t(c(s,as.numeric(p.ridge))), file = paste0('../result/Ridge_hdi_noH_v1_0716_',info.det,'.txt'),row.names = F, col.names = F, quote = F, append = TRUE)
write.table(t(c(s,as.numeric(p.grace))), file = paste0('../result/Grace_hdi_noH_v1_0716_',info.det,'.txt'),row.names = F, col.names = F, quote = F, append = TRUE)
write.table(t(c(s,as.numeric(p.lasso))), file = paste0('../result/LASSO_hdi_noH_v1_0716_',info.det,'.txt'),row.names = F, col.names = F, quote = F, append = TRUE)
write.table(t(c(s,as.numeric(p.dscore))), file = paste0('../result/dscore_hdi_noH_v1_0716_',info.det,'.txt'),row.names = F, col.names = F, quote = F, append = TRUE)
