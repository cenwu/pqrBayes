#' make predictions from a pqrBayes object
#'
#' make predictions from a pqrBayes object
#'
#' @param object pqrBayes object.
#' @param g.new a matrix of new predictors (e.g. genetic factors) at which predictions are to be made.
#' @param u.new a vector of new environmental factor at which predictions are to be made.
#' @param e.new a vector or matrix of new clinic covariates at which predictions are to be made.
#' @param y.new a vector of the response of new observations. If provided, the prediction error will be computed based on Y.new.
#' @param quant the quantile for the response variable.  The default is 0.5.
#' @param kn the number of interior knots for B-spline.
#' @param degree the degree of B-spline basis.
#' @param ... other predict arguments
#' 
#' @details g.new (u.new) must have the same number of columns as g (u) used for fitting the model. By default, the clinic covariates are NULL unless 
#' provided. The predictions are made based on the posterior estimates of coefficients in the pqrBayes object.
#'
#' If y.new is provided, the prediction error will be computed based on the check loss.
#'
#' @return  an object of class `pqrBayes.pred' is returned, which is a list with components:
#' \item{error}{prediction error. error is NULL if y.new=NULL.}
#' \item{y.pred}{predicted values of the new observations.}
#'
#' @rdname predict.pqrBayes
#' @seealso \code{\link{pqrBayes}}
#'
#' @export
predict.pqrBayes=function(object, g.new, u.new, e.new=NULL, y.new=NULL, quant=0.5, kn=2, degree=2,...){
  p = dim(g.new)[2]
  
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
    c1.C[i]=stats::quantile(object$coefficients$GS.alpha[5001:10000,i],0.5)
  }
  
  c2.C=rep(0,dim(object$coefficients$GS.beta)[2])
  for (i in 1:dim(object$coefficients$GS.beta)[2]) {
    c2.C[i]=stats::quantile(object$coefficients$GS.beta[5001:10000,i],0.5)
  }
  
  coeffmatrix.C=as.matrix(cbind(c1.C,matrix(c2.C,nrow = d)))
  
  gamma.hat=pi.u%*% coeffmatrix.C
  
  if(!is.null(e.new)){
    q = dim(e.new)[2]
    
    
    beta.hat=rep(0,q)
    for (i in 1:q) {
      beta.hat[i]=stats::quantile(object$coefficients$GS.alpha[5001:10000,i+d],0.5)
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
  
  
  
  
  
  pqrBayes.pred = list(error=pred.error, y.pred=y.hat)
  class(pqrBayes.pred) = "pqrBayes.pred"
  return(pqrBayes.pred)
  #pred
}

