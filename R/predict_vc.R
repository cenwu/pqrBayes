predict_vc=function(obj, g.new, u.new, e.new, y.new, quant,...){
  p = dim(g.new)[2]
  iterations = obj$obj$iterations
  kn = obj$obj$kn
  degree = obj$obj$degree
  x = cbind(1,g.new)
  n = dim(g.new)[1]; 
  
  
  
  ## basis expansion
  
  d=kn+degree+1
  u.k = seq(0, 1, length=kn+2)[-c(1,kn+2)]
  Knots = as.numeric(stats::quantile(u.new, u.k))
  pi.u = splines::bs(u.new, knots=Knots, intercept=TRUE, degree=degree)[,1:(d)]
  
  
  
  # 50% quantile of posterior samples
  
  c1.C=rep(0,d)
  for (i in 1:d) {
    c1.C[i]=stats::quantile(obj$coefficients$GS.alpha[(iterations/2+1):iterations,i],0.5)
  }
  
  c2.C=rep(0,dim(obj$coefficients$GS.beta)[2])
  for (i in 1:dim(obj$coefficients$GS.beta)[2]) {
    c2.C[i]=stats::quantile(obj$coefficients$GS.beta[(iterations/2+1):iterations,i],0.5)
  }
  
  coeffmatrix.C=as.matrix(cbind(c1.C,matrix(c2.C,nrow = d)))
  
  gamma.hat=pi.u%*% coeffmatrix.C
  
  if(!is.null(e.new)){
    q = dim(e.new)[2]
    
    
    beta.hat=rep(0,q)
    for (i in 1:q) {
      beta.hat[i]=stats::quantile(obj$coefficients$GS.alpha[(iterations/2+1):iterations,i+d],0.5)
    }
    
    if(y.new){  
      # prediction error
      y.hat=res=rep(0,n)
      for (i in 1:n) {
        y.hat[i]=t(x[i,])%*%gamma.hat[i,]+t(e.new[i,])%*%beta.hat
        if((y.new[i]-y.hat[i]) >= 0){res[i] = quant*(y.new[i]-y.hat[i])}
        else{res[i] = (quant -1)*(y.new[i]-y.hat[i])}
      }
      
      pred.error=mean(res)
    }else{
      pred.error=NULL
    }
    
    
  }else{
    if(y.new){  
      # prediction error
      y.hat=res=rep(0,n)
      for (i in 1:n) {
        y.hat[i]=t(x[i,])%*%gamma.hat[i,]
        if((y.new[i]-y.hat[i]) >= 0){res[i] = quant*(y.new[i]-y.hat[i])}
        else{res[i] = (quant -1)*(y.new[i]-y.hat[i])}
      }
      
      pred.error=mean(res)
    }else{
      pred.error=NULL
    }
  }
  
  
  
  
  
  pqrBayes_vc= list(error=pred.error, y.pred=y.hat)
  #class(pqrBayes.pred) = "pqrBayes.pred"
  return(pqrBayes_vc)
  #pred
}
