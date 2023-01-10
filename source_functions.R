trace <- function(A) {
  n <- dim(A)[1] # get dimension of matrix
  tr <- 0 # initialize trace value
  
  # Loop over the diagonal elements of the supplied matrix and add the element to tr
  for (k in 1:n) {
    l <- A[k,k]
    tr <- tr + l
  }
  return(tr[[1]])
}

FDAimage.est.ho <- function(B,Q2,K,X,Y,lambda){
  n=nrow(Y)
  nx=ncol(X)
  npix=ncol(Y)
  BQ2=B%*%Q2
  BQ2=as.matrix(BQ2)
  K=as.matrix(K) # THIS IS THE PENALTY P!
  J=ncol(BQ2)
  
  WW=kronecker(crossprod(X),crossprod(BQ2))
  rhs=rowSums(sapply(1:n, function(iter) as.matrix(kronecker(X[iter,],crossprod(BQ2,Y[iter,])))))
  P=as.matrix(crossprod(Q2,K)%*%Q2) # K is actually the P matrix in the paper!!
  
  lambda=as.matrix(lambda)
  nlam=nrow(lambda)
  gamma_all=list(nlam); beta_all=list(nlam); sse_all=c();
  for(il in 1:nlam){
    Lam=diag(rep(lambda[il],nx))
    Dlam=as.matrix(kronecker(Lam,P))
    lhs=WW+Dlam
    theta=solve(lhs,rhs)
    theta.mtx=matrix(theta,J,nx)
    gamma=Q2%*%theta.mtx
    gamma_all[[il]]=gamma
    beta=BQ2%*%theta.mtx
    beta_all[[il]]=beta
    Yhat=tcrossprod(X,beta)
    ssei=apply((Y-Yhat)^2,1,sum)
    sse=sum(ssei)
    sse_all=c(sse_all,sse)
  }
  
  sse.min = which(sse_all == min(sse_all))
  lambdac = lambda[sse.min]
  
  Lam=diag(rep(lambdac,nx))
  Dlam=as.matrix(kronecker(Lam,P))
  lhs=WW+Dlam
  theta=solve(lhs,rhs)
  theta.mtx=matrix(theta,J,nx)
  gamma=Q2%*%theta.mtx
  gamma_all=gamma
  beta=BQ2%*%theta.mtx
  beta_all=beta
  Yhat=tcrossprod(X,beta)
  
  list(beta=beta_all,gamma=gamma_all)
}

quasi_GCV <- function(B, Q2, K, intercep, gamma, theta, lambda, surv, ind.inside){
  
  A = matrix(surv$A, ncol = 1)
  C = matrix(surv$C, ncol = 1)
  U = matrix(surv$U, ncol = 1)
  
  M = surv$M[,ind.inside]

  BM = M %*% B
  Q2BM =  BM %*% Q2
  
  G = cbind(rep(1, nrow(surv)), A, C, U, Q2BM) 
  
  eta = intercep*1 + gamma[1]*(A) + gamma[2]*(C)+ gamma[3]*(U) + Q2BM %*% matrix(theta, ncol = 1)
  
  b.2f = exp(eta)/(1+exp(eta))^2
  
  w = diag(as.vector(b.2f))
  w1.2 = diag(sqrt(as.vector(b.2f)))
  
  y = surv$event
  
  b.f = exp(eta)/(1+exp(eta))
  y.tilde = eta + (y - b.f)/(b.f*(1-b.f))
  
  K=as.matrix(K) 
  Q2.2 = cbind(rep(0, nrow(Q2)),rep(0, nrow(Q2)),rep(0, nrow(Q2)),rep(0, nrow(Q2)), Q2) 
  D2.2 = t(Q2.2)%*%K%*%Q2.2
  
  #-----------------------------------------------------#
  lambda <- as.matrix(lambda)
  nlam <- nrow(lambda)
  gcv_all <- c();
  
  G.b = G
  tG = t(G)
  w.b = w
  
  for(il in 1:nlam){
    Lam <- lambda[il]
    Dlam=as.matrix(kronecker(Lam,D2.2))
    
    mid = solve( t(G) %*% w %*% G + (1/2)*Dlam )
    interm = tG %*% w.b
    H = (G.b %*% mid %*% interm)
    yhat = H %*% y.tilde
    res <- (y.tilde - yhat)
    gcv = (1/n)*sum(diag(w)*res^2)/(trace(diag(1, nrow = nrow(H), ncol = ncol(H)) - H)/n)^2
    
    gcv_all[il] = gcv
  }
  
  j <- which.min(gcv_all)
  lambdac <- lambda[j]
  Dlam=as.matrix(kronecker(lambdac,D2.2))
  
  lhs = t(G) %*% w %*% G + (1/2)*Dlam
  rhs = t(G) %*% w %*% y.tilde
  
  theta.est = solve(lhs, rhs)
  
  return(list(theta.est, lambdac))
}

est.GCV = function(data, max.step, lambda){
  surv = data
  M = surv$M[,ind.inside]
  X = as.matrix(cbind(rep(1, nrow(surv)),surv$A, surv$C, surv$U))
  
  #------- M MODEL ----------#
  fit.M = FDAimage.est.ho(B = B, Q2 = Q2, K = K, X = X, Y = M, lambda = lambda)
  
  intercep_M = fit.M$beta[,1]
  alpha_s = (matrix(fit.M$beta[,2], nrow = 1)) # with intercept
  
  Yhat.i = as.matrix(X %*% t(fit.M$beta))
  R.i = M - Yhat.i
  sig2 = mean(apply(R.i^2, 2, sum)/n)
  
  #alpha_s = beta.A
  
  #------- Y MODEL ----------#
  
  theta.init = rep(0.01, ncol(Q2))
  gamma.init = rep(0.01, (ncol(X)-1)) # one for X and one for A
  intc.init = 0.01
  param.init = c(intc.init, gamma.init, theta.init)
  
  threshold = 1e-4
  step = 1
  
  while(step < max.step){
    
    est  = quasi_GCV(B, Q2, K,intercep=param.init[1], gamma = param.init[2:(2+length(gamma.init)-1)],
                     theta = param.init[-c(1:ncol(X))], lambda = lambda,
                     surv, ind.inside)
    
    param.est = est[[1]]
    param.est = as.vector(param.est)
    
    if( abs(max(abs(param.est - param.init))) <= threshold ){
      break
    }else{
      param.init = param.est
      step = step + 1
    }
  }
  
  theta.est =  matrix(param.est[-c(1:ncol(X))], ncol = 1)
  b.est = t(theta.est) %*% t(Q2)
  beta_s = B %*% t(b.est)
  intercep_Y = param.est[1]
  
  gamma.A = param.est[2]; gamma.C = param.est[3];
  return(list(gamma.A, alpha_s, beta_s, est[[2]], intercep_M, intercep_Y, gamma.C) )
  
}

