Robust <- function(g, y, u, e,quant, iterations, kn, degree,sparse, hyper,debugging){
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
  
  xi1 = (1-2*quant)/(quant*(1-quant))
  xi2 = sqrt(2/(quant*(1-quant)))
  hatTau = 1
  hatV = rep(1,n) 
  hatEtaSq = 1
  hatSg = rep(1, p)
  hatPi = 0.5
  
  sh0_1 = ifelse(is.null(hyper$a0), 1, hyper$a0)
  sh0_0 = ifelse(is.null(hyper$b0), 1, hyper$b0)
  
  a = ifelse(is.null(hyper$c1), 1, hyper$c1)
  b = ifelse(is.null(hyper$c2), 1, hyper$c2)
  
  r = ifelse(is.null(hyper$d2), 1, hyper$d2)
  hatbeta=matrix(hat.r,ncol=p) + 10^-5
  progress = ifelse(debugging, 10^(floor(log10(iterations))-1), 0)
  if(sparse){fit=BRGL_SS(xx1, y, CLC, p, d, iterations, hatAlpha, hatbeta, hatTau, hatV, hatSg, invSigAlpha0, hatPi, hatEtaSq,
                         xi1, xi2, r, a, b, sh0_1, sh0_0, progress)}
  else{fit=BRGL(xx1, y, CLC, p, d, iterations, hatAlpha, hatbeta, hatTau, hatV, hatSg, invSigAlpha0, hatEtaSq,
                xi1, xi2, r, a, b, progress)}
  out = fit
  return(out)
}
