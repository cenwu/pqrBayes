predict_bin=function(obj, g.new, e.new, y.new, quant,...){
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
    y.hat=rep(0,n)
    for (i in 1:n) {
      y.hat[i]=t(x[i,])%*%gamma.hat+t(e.new[i,])%*%beta.hat[-1]
    }
    y.hat.new = as.numeric(y.hat>0)
    pred.error <- mean(y.hat.new != y.new)
    
  }else{
    
    t2=as.matrix(obj$coefficients$GS.alpha)
    beta.hat = stats::median(t2)
    
    gamma.hat = c(beta.hat,coeff.est)
    
    # prediction error
    y.hat=rep(0,n)
    for (i in 1:n) {
      y.hat[i]=t(x[i,])%*%gamma.hat
    }
    
    y.hat.new = as.numeric(y.hat>0)
    pred.error <- mean(y.hat.new != y.new)
    
  }
  
  
  pqrBayes_bin= list(error=pred.error, y.pred=y.hat.new)
  #class(pqrBayes.pred) = "pqrBayes.pred"
  return(pqrBayes_bin)
  #pred
}
