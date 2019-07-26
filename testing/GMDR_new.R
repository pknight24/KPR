#--------------------------------------------------------------
# GMDR: generalized matrix decomposition regression
# Yue Wang: 20181009
#-------------------------------------------------------------
#library(glmnet)
#library(scalreg, lib.loc = '/pine/scr/t/a/taryue/Rlibs/')
#library(natural)
#library(natural, lib.loc = '/pine/scr/t/a/taryue/Rlibs/')


#--------------------------------------------
# Part I: GMD function

get_uv = function(X, H, Q, u_0, v_0){
  
  u = X%*%Q%*%v_0/as.numeric(sqrt(t(v_0)%*%t(Q)%*%t(X)%*%H%*%X%*%Q%*%v_0))
  v = t(X)%*%H%*%u/as.numeric(sqrt(t(u)%*%t(H)%*%X%*%Q%*%t(X)%*%H%*%u))

  return(list(u = u, v = v))
    
}


GMD =function(X, H, Q, K){
  
  n = dim(X)[1]
  p = dim(X)[2]
  # output matrix/vec
  U = matrix(0, n, K)
  V = matrix(0, p, K)
  D = rep(0,K)
  
  X_0 = X
  u_0 = c(1, rep(0,n-1))
  v_0 = c(1, rep(0,p-1))
  
  for(iter in 1:K){
    
    
    error = 1
 
    while(error > 1e-5){
      
      temp.uv = get_uv(X_0, H, Q, u_0, v_0)
      
      u = temp.uv$u
      v = temp.uv$v
      
      error = norm(u - u_0, "2") + norm(v - v_0, "2")
      #print(error)
      
      u_0 = u
      v_0 = v
   
    }
    
    U[,iter] = u
    V[,iter] = v
    d = t(u)%*%H%*%X_0%*%Q%*%v
    D[iter] =  d
    
    X_0 = X_0 - u%*%d%*%t(v)
  
  }
 
  return(list(U = U, V = V, D = D, H = H, Q = Q))
  
}


#----------------------------------------------------------------
# Part II: GMDR estimation
#----------------------------------------------------------------


GMDR.est = function(Y, X, H, Q){
  
  n = dim(X)[1]
  p = dim(X)[2]
  H.r = qr(t(X)%*%H%*%X%*%Q)$rank
  gmd.x = GMD(X, H, Q, H.r) # the full GMD of X
  U = gmd.x$U
  V = gmd.x$V
  D = gmd.x$D
  perc.var = D^2/sum(D^2)
  
  
  #get rid of the 0's
  #U = U[,D>1e-10]
  #V = V[,D>1e-10]
  #D = D[D>1e-10]
  
  U = U[,perc.var >= 0.01]
  V = V[,perc.var >= 0.01]
  D = D[perc.var >= 0.01]
  
  H.r = length(D)
  
  #weighted least square
  eta = U%*%diag(D)
  gamma = diag(1/D)%*%t(U)%*%H%*%Y


  R.2 = D^2*gamma^2
  R.2.sort = sort(R.2, decreasing = T)
  R.2.order = order(R.2, decreasing = T)
 #R.2.sort = R.2.sort/sum(R.2.sort)

  #sort with respect H.2
  V.sort = V[,R.2.order]
  eta.sort = eta[, R.2.order]
  D.sort = D[R.2.order]
  U.sort = U[,R.2.order]
  
  # consider a candidate set of K's.
  GCV.stat = rep(0, H.r)
  W = rep(0,H.r)
  beta.out = matrix(0, p, H.r)
  
  for(i in 1:H.r){
    
    W.cal = W 
    W.cal[1:i] = 1
    gamma.K = diag(W.cal/D.sort,nrow = H.r, ncol = H.r)%*%t(U.sort)%*%H%*%Y
    beta.out[,i] = Q%*%V.sort%*%gamma.K
    
    # GCV:
    A_lambda = U.sort%*%diag(W.cal)%*%t(U.sort)%*%H
    pred.error = as.numeric((diag(n) - A_lambda)%*%Y)
    GCV.stat[i] = n*t(pred.error)%*%H%*%pred.error/(sum(1 - diag(A_lambda)))^2
  
  }
  
  i.opt = order(GCV.stat)[1]
  
  #return
  return(list(GMD.x = gmd.x, set.opt = R.2.order[1:i.opt], beta.opt = beta.out[,i.opt]))
  
}


#--------------------------------------------------------------------------------
# the following is just for comparison purpose (original PC selection procedure)

GMDR.est.org = function(Y, X, H, Q){
  
  n = dim(X)[1]
  p = dim(X)[2]
  H.r = qr(t(X)%*%H%*%X%*%Q)$rank
  gmd.x = GMD(X, H, Q, H.r) # the full GMD of X
  U = gmd.x$U
  V = gmd.x$V
  D = gmd.x$D
  
  #get rid of the 0's
  U = U[,D>1e-10]
  V = V[,D>1e-10]
  D = D[D>1e-10]
  
  H.r = length(D)
  

  # consider a candidate set of K's.
  GCV.stat = rep(0, H.r)
  W = rep(0,H.r)
  beta.out = matrix(0, p, H.r)
  
  for(i in 1:H.r){
    
    W.cal = W 
    W.cal[1:i] = 1
    gamma.K = diag(W.cal/D,nrow = H.r, ncol = H.r)%*%t(U)%*%H%*%Y
    beta.out[,i] = Q%*%V%*%gamma.K
    
    # GCV:
    A_lambda = U%*%diag(W.cal)%*%t(U)%*%H
    pred.error = as.numeric((diag(n) - A_lambda)%*%Y)
    GCV.stat[i] = n*t(pred.error)%*%H%*%pred.error/(sum(1 - diag(A_lambda)))^2
    
  }
  
  i.opt = order(GCV.stat)[1]
  
  #return
  return(list(GMD.x = gmd.x, set.opt = R.2.order[1:i.opt], beta.opt = beta.out[,i.opt]))
  
}

#-----------------------------------------------------------------------
# Part III: GMDR test
#-----------------------------------------------------------------------

# use general mu.
# mu: test parameter
# r: sparsity parameter
GMDR.test = function(Y, X, GMD.est, mu, r, weight = TRUE){

  U = GMD.est$GMD.x$U
  V = GMD.est$GMD.x$V
  D = GMD.est$GMD.x$D
  H = GMD.est$GMD.x$H
  Q = GMD.est$GMD.x$Q
  n = dim(X)[1]
  p = dim(X)[2]
  beta.GMDR = GMD.est$beta.opt # before correction
  set.GMDR = GMD.est$set.opt
  W = rep(0, length(D))
  W[set.GMDR] = 1

  
  # bias-correction
  eigen.H = eigen(H)
  values.H = eigen.H$values
  vectors.H = eigen.H$vectors
  L.H = vectors.H%*%diag(sqrt(values.H))
  eigen.Q = eigen(Q)
  values.Q = eigen.Q$values
  vectors.Q = eigen.Q$vectors
  #print(vectors.Q[,1])
  L.Q = vectors.Q%*%diag(sqrt(values.Q))

  X.tilde = L.H%*%X%*%vectors.Q
  Y.tilde = L.H%*%Y
  
  # using natural lasso method
    
    olasso.fit = olasso_cv(X.tilde, Y.tilde)
    sigmaepsi.hat = olasso.fit$sig_obj
    
    if(weight == TRUE){
      lambda.opt = cv.glmnet(X.tilde, Y.tilde, alpha = 1, intercept = FALSE, standardize = TRUE, penalty.factor = 1/sqrt(values.Q))$lambda.min
      beta.glasso = glmnet(X.tilde, Y.tilde, alpha = 1, lambda = lambda.opt, intercept = FALSE, standardize = TRUE, penalty.factor = 1/sqrt(values.Q))$beta
       # print(beta.glasso)
      }
   
    if(weight == FALSE){
      lambda.opt = cv.glmnet(X.tilde, Y.tilde, alpha = 1, intercept = FALSE, standardize = TRUE)$lambda.min
      beta.glasso = glmnet(X.tilde, Y.tilde, alpha = 1, lambda = lambda.opt, intercept = FALSE, standardize = TRUE)$beta
      # print(beta.glasso) 
      }
    
    #print(beta.glasso)
    beta.lasso = as.numeric(vectors.Q%*%beta.glasso)
    #lines(1:100, beta.lasso, col = 'orange')
    
    Xi = diag(Q%*%V%*%diag(W)%*%t(V))
    beta.cor = Q%*%V%*%diag(W)%*%t(V)%*%beta.lasso - (1 - mu)*Xi*beta.lasso - mu*beta.lasso
    beta.hat = beta.GMDR  - beta.cor 
    # print(beta.hat)
    
  # covariance
  cov.hat = sigmaepsi.hat^2*Q%*%V%*%diag(D^(-2)*W*W)%*%t(V)%*%Q
  diag.cov.hat = diag(cov.hat)
  bound.mat = (Q%*%V%*%diag(W)%*%t(V) - (1- mu)*diag(Xi) - mu*diag(rep(1,p)))%*%L.Q
  bound.hat = apply(bound.mat, 1, function(x){max(abs(x))})*(log(p)/n)^(0.5 - r) # sparsity parameter
  # print(bound.hat)
  
  # p-value
  beta.temp = abs(beta.hat) - bound.hat
  p.vec = rep(0,p)
  for(i in 1:p){
    # print(beta.temp)
    
    if(beta.temp[i] > 0){p.vec[i] = 2*(1 - pnorm(beta.temp[i]/sqrt(diag.cov.hat[i])))}
    if(beta.temp[i] <= 0){p.vec[i] = 1}
    
  }
  names(p.vec) <- colnames(X)
  return(p.vec)
  
}


# use general mu.
# mu: test parameter
# r: sparsity parameter
KPR.test = function(KPR.est, mu, r, weight = TRUE){
  
  U = KPR.est$GMD.x$U
  V = KPR.est$GMD.x$V
  D = KPR.est$GMD.x$D
  H = KPR.est$GMD.x$H
  Q = KPR.est$GMD.x$Q
  beta.KPR = KPR.est$beta.opt # before correction
  lambda.KPR = KPR.est$lambda.opt
  W = D^2/(D^2 + lambda.KPR)
  
  # bias-correction
  eigen.H = eigen(H)
  values.H = eigen.H$values
  vectors.H = eigen.H$vectors
  L.H = vectors.H%*%diag(sqrt(values.H))
  eigen.Q = eigen(Q)
  values.Q = eigen.Q$values
  vectors.Q = eigen.Q$vectors
  #print(vectors.Q[,1])
  L.Q = vectors.Q%*%diag(sqrt(values.Q))
  
  X.tilde = L.H%*%X%*%vectors.Q
  Y.tilde = L.H%*%Y
  
  # using natural lasso method
  
  olasso.fit = olasso_cv(X.tilde, Y.tilde)
  sigmaepsi.hat = olasso.fit$sig_obj
  
  if(weight == TRUE){
    lambda.opt = cv.glmnet(X.tilde, Y.tilde, alpha = 1, intercept = FALSE, standardize = TRUE, penalty.factor = 1/sqrt(values.Q))$lambda.min
    beta.glasso = glmnet(X.tilde, Y.tilde, alpha = 1, lambda = lambda.opt, intercept = FALSE, standardize = TRUE, penalty.factor = 1/sqrt(values.Q))$beta
  }
  
  if(weight == FALSE){
    lambda.opt = cv.glmnet(X.tilde, Y.tilde, alpha = 1, intercept = FALSE, standardize = TRUE)$lambda.min
    beta.glasso = glmnet(X.tilde, Y.tilde, alpha = 1, lambda = lambda.opt, intercept = FALSE, standardize = TRUE)$beta
  }
  
  #print(beta.glasso)
  beta.lasso = as.numeric(vectors.Q%*%beta.glasso)
  #lines(1:100, beta.lasso, col = 'orange')
  
  Xi = diag(Q%*%V%*%diag(W)%*%t(V))
  beta.cor = Q%*%V%*%diag(W)%*%t(V)%*%beta.lasso - (1 - mu)*Xi*beta.lasso - mu*beta.lasso
  beta.hat = beta.KPR  - beta.cor 
  #print(beta.hat)
  
  # covariance
  cov.hat = sigmaepsi.hat^2*Q%*%V%*%diag(D^(-2)*W*W)%*%t(V)%*%Q
  diag.cov.hat = diag(cov.hat)
  bound.mat = (Q%*%V%*%diag(W)%*%t(V) - (1- mu)*diag(Xi) - mu*diag(rep(1,p)))%*%L.Q
  bound.hat = apply(bound.mat, 1, function(x){max(abs(x))})*(log(p)/n)^(0.5 - r) # sparsity parameter
  #print(bound.hat)
  
  # p-value
  beta.temp = abs(beta.hat) - bound.hat
  p.vec = rep(0,p)
  for(i in 1:p){
    
    if(beta.temp[i] > 0){p.vec[i] = 2*(1 - pnorm(beta.temp[i]/sqrt(diag.cov.hat[i])))}
    if(beta.temp[i] <= 0){p.vec[i] = 1}
    
  }
  names(p.vec) <- colnames(X)
  return(p.vec)
  
}






#-------------------------------------------------
# Grace (Zhao and Shojaie, 2016): cannot handle H.
# for comparison purpose
#-------------------------------------------------
# grace.est = function(Y, X, Q){
#   
#   eigen.Q = eigen(Q)
#   values.Q = eigen.Q$values
#   vectors.Q = eigen.Q$vectors
#   L.Q = vectors.Q%*%diag(sqrt(values.Q))
#   #Q.sqrt = vectors.Q%*%diag(sqrt(values.Q))%*%t(vectors.Q) 
#   Q.inv = vectors.Q%*%diag(1/values.Q)%*%t(vectors.Q)
#   
#   X.tilde = X%*%L.Q
#   #Cov_X = t(X)%*%X
#   # ridge regression
#   
#   cv.ridge = cv.glmnet(X.tilde, Y, alpha = 0, intercept = FALSE, standardize = FALSE)
#   fit.ridge = glmnet(X.tilde,Y, family = "gaussian", alpha = 0, lambda = cv.ridge$lambda.min, intercept = FALSE, standardize = FALSE)
#   beta.grace = L.Q%*%fit.ridge$beta
#   
#   return(list(lambda.opt = cv.ridge$lambda.min, beta.opt = beta.grace, X = X, Y = Y, Q.inv = Q.inv))
# }
# 
# # following Zhao's paper.
# grace.test = function(grace.est, r){
#   
#   X = grace.est$X
#   n = dim(X.tilde)[1]
#   p = dim(X.tilde)[2]
#   Y = grace.est$Y
#   Q.inv = grace.est$Q.inv
#   print(dim(Q.inv))
#   lambda.grace = grace.est$lambda.opt/2
#   beta.grace = grace.est$beta.opt
#   Cov_X = t(X)%*%X
#   
#   # using natural lasso method
#   olasso.fit = olasso_cv(X, Y)
#   sigmaepsi.hat = olasso.fit$sig_obj
#   lambda.opt = cv.glmnet(X, Y, alpha = 1, intercept = FALSE, standardize = FALSE)$lambda.min
#   beta.lasso = glmnet(X, Y, alpha = 1, lambda = lambda.opt, intercept = FALSE, standardize = FALSE)$beta
#   
#   # bias-correction
#   beta.hat = beta.grace + lambda.grace*ginv(Cov_X + lambda.grace*Q.inv)%*%Q.inv%*%beta.lasso
#   cov.hat = sigmaepsi.hat^2*ginv(Cov_X + lambda.grace*Q.inv)%*%Cov_X%*%ginv(Cov_X + lambda.grace*Q.inv)
#   diag.cov.hat = sqrt(diag(cov.hat))
#   
#   # bound:
#   Xi = lambda.grace*ginv(Cov_X + lambda.grace*Q.inv)%*%Q.inv
#   bound.mat = Xi - diag(diag(Xi))
#   
#   bound.hat = apply(bound.mat, 1, function(x){max(abs(x))})*(log(p)/n)^(0.5 - r) # sparsity parameter
#   
#   # p-value
#   beta.temp = abs(beta.hat) - bound.hat
#   p.vec = rep(0,p)
#   for(i in 1:p){
#     
#     if(beta.temp[i] > 0){p.vec[i] = 2*(1 - pnorm(beta.temp[i]/sqrt(diag.cov.hat[i])))}
#     if(beta.temp[i] <= 0){p.vec[i] = 1}
#     
#   }
#   
#   return(p.vec)
# }
# 
# 
# 
