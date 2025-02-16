pqrBayes_vc <- function(g, y, u, e=NULL,quant=0.5, iterations=10000, kn=2, degree=2, robust=TRUE,sparse=TRUE, hyper=NULL,debugging=FALSE){
  
  if(robust){
    out = Robust_vc(g, y, u, e,quant, iterations, kn, degree,sparse, hyper,debugging)
  }else{
    out = NonRobust_vc(g, y, u, e,iterations, kn, degree,sparse,debugging)
  }
  coefficient = list(GS.alpha=out$fit$GS.alpha, GS.beta=out$fit$GS.beta)
  fit = list(obj = out,coefficients = coefficient)
  return(fit)
}
