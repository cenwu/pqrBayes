pqrBayes_g <- function(g, y, e, d, quant=0.5, iterations=10000, burn.in,robust=TRUE,sparse=TRUE, hyper=NULL,debugging=FALSE){
  if(iterations<1) stop("iterations must be a positive integer.")
  if(is.null(burn.in)){
    BI = floor(iterations)/2
    if(iterations<=BI) stop("iterations must be larger than burn.in.")
  }else if(burn.in>=1){
    BI = as.integer(burn.in)
  }else{
    stop("burn.in must be a positive integer.")
  }
  if(robust){
    out = Robust_g(g, y, e, d, quant, iterations,sparse, hyper,debugging)
  }else{
    out = NonRobust_g(g, y, e, d, iterations,sparse,debugging)
  }
  coefficient = list(GS.alpha=out$fit$GS.alpha[-c(1:BI),],GS.beta=out$fit$GS.beta[-c(1:BI),])
  fit = list(obj = out,coefficients = coefficient)
  return(fit)
}
