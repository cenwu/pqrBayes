predict_lin=function(obj, g.new, e.new, y.new, quant,...){
  p = dim(g.new)[2]
  x = cbind(1,g.new)
  n = dim(g.new)[1]; 
  coeff.est = c()
  for(j in 1:dim(obj$coefficients$GS.beta)[2])
  {
    t1=as.matrix(obj$coefficients$GS.beta[,j])
    coeff.est[j] = stats::median(t1)
  }
  
  if(!is.null(e.new)){
    q = dim(e.new)[2]
  
    beta.hat=rep(0,q+1)
    for (i in 1:(q+1)) {
      t2=as.matrix(obj$coefficients$GS.alpha[,i])
      beta.hat[i] = stats::median(t2)
    }
    gamma.hat = c(beta.hat[1],coeff.est)
    
      # prediction error
      y.hat=res=rep(0,n)
      for (i in 1:n) {
        y.hat[i]=t(x[i,])%*%gamma.hat+t(e.new[i,])%*%beta.hat[-1]
        if((y.new[i]-y.hat[i]) >= 0){res[i] = quant*(y.new[i]-y.hat[i])}
        else{res[i] = (quant -1)*(y.new[i]-y.hat[i])}
      }
      
      pred.error=mean(res)
    
    
  }else{
    
      t2=as.matrix(obj$coefficients$GS.alpha)
      beta.hat = stats::median(t2)
    
    gamma.hat = c(beta.hat,coeff.est)
    
      # prediction error
      y.hat=res=rep(0,n)
      for (i in 1:n) {
        y.hat[i]=t(x[i,])%*%gamma.hat
        if((y.new[i]-y.hat[i]) >= 0){res[i] = quant*(y.new[i]-y.hat[i])}
        else{res[i] = (quant -1)*(y.new[i]-y.hat[i])}
      }
      
      pred.error=mean(res)
    
  }
  
  
  pqrBayes_lin= list(error=pred.error, y.pred=y.hat)
  #class(pqrBayes.pred) = "pqrBayes.pred"
  return(pqrBayes_lin)
  #pred
}
