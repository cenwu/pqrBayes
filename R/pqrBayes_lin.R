pqrBayes_lin <- function(g, y, e,quant=0.5, iterations=10000, robust=TRUE,sparse=TRUE, hyper=NULL,debugging=FALSE){
  
  if(robust){
    out = Robust_lin(g, y, e, quant, iterations,sparse, hyper,debugging)
  }else{
    out = NonRobust_lin(g, y, e,iterations,sparse,debugging)
  }
  coefficient = list(GS.alpha=out$fit$GS.alpha,GS.beta=out$fit$GS.beta)
  fit = list(obj = out,coefficients = coefficient)
  return(fit)
}
