NonRobust_vc <- function(g, y, u, e,iterations, kn, degree,sparse,debugging){
  p = dim(g)[2]
  
  x = cbind(1,g)
  n = length(y); 
  
  
  
  ## basis expansion
  
  d=kn+degree+1
  u.k = seq(0, 1, length=kn+2)[-c(1,kn+2)]
  Knots = as.numeric(stats::quantile(u, u.k))
  pi.u = splines::bs(u, knots=Knots, intercept=TRUE, degree=degree)[,1:(d)]
  
  
  xx = as.data.frame(matrix(0, n, (p+1)*d))
  for(j in 1:(p+1)){
    last = j*d; first = last-d+1
    xx[,first:last] = pi.u*x[,j]
  }
  xx = as.matrix(xx)
  
  if(!is.null(e)){
    q = dim(e)[2]
    xxwe = cbind(xx, e)
    lasso.cv = glmnet::cv.glmnet(xxwe,y,alpha=1,nfolds=5)
    lambda.cv = lasso.cv$lambda.min;
    lasso.fit = glmnet::glmnet(xxwe, y, family="gaussian",alpha=1,nlambda=50)
    coeff.array = as.vector(stats::predict(lasso.fit, s=lambda.cv, type="coefficients"))[-1];
    
    hat.m = coeff.array[1:d]      ## coeff for varying intercept
    hat.r = coeff.array[(d+1):((p+1)*d)] ## coeff for varying part
    hat.clin = coeff.array[((p+1)*d+1):dim(xxwe)[2]] ## coeff for clinic covariates
    
    
    xx1=xx[,-(1:d)]
    CLC=cbind(pi.u,e)
    hatAlpha=c(hat.m,hat.clin)
    invSigAlpha0 = diag(10^-3, (d+q))
  }else{
    q = 0
    xxwe = xx
    lasso.cv = glmnet::cv.glmnet(xxwe,y,alpha=1,nfolds=5)
    lambda.cv = lasso.cv$lambda.min;
    lasso.fit = glmnet::glmnet(xxwe, y, family="gaussian",alpha=1,nlambda=50)
    coeff.array = as.vector(stats::predict(lasso.fit, s=lambda.cv, type="coefficients"))[-1];
    
    hat.m = coeff.array[1:d]      ## coeff for varying intercept
    hat.r = coeff.array[(d+1):((p+1)*d)] ## coeff for varying part
    
    
    xx1=xx[,-(1:d)]
    CLC=pi.u
    hatAlpha=hat.m
    invSigAlpha0 = diag(10^-3, (d+q))
  }
  
  invTAUsq.star = rep(0.1, p)
  hat.pi.s = 0.9
  lambda.star = 1
  hat.sigma.sq = 1
  a.star=1
  b.star = 1.5
  alpha = 0.2
  gamma = 0.1
  mu.star=nu.star=1 
  progress = ifelse(debugging, 10^(floor(log10(iterations))-1), 0)
  if(sparse){fit=BGLPointMass(xx1, y, CLC, p, d, iterations, hatAlpha, hat.r, invTAUsq.star, invSigAlpha0, hat.pi.s,
                              lambda.star, hat.sigma.sq, a.star, b.star, alpha, gamma, mu.star, nu.star, progress)}
  else{fit=BGL(xx1, y, CLC, p, d, iterations, hat.r, hatAlpha, invTAUsq.star, invSigAlpha0, lambda.star, 
               hat.sigma.sq, a.star, b.star, alpha, gamma, progress)}
  out=list(fit=fit,iterations=iterations,kn=kn,degree=degree)
  return(out)
}
