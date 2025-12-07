pqrBayes_vc <- function(g, y, e=NULL,quant=0.5, iterations=10000, burn.in, robust=TRUE,prior="SS", hyper=NULL,debugging=FALSE){
  if(iterations<1) stop("iterations must be a positive integer.")
  if(is.null(burn.in)){
    BI = floor(iterations)/2
    if(iterations<=BI) stop("iterations must be larger than burn.in.")
  }else if(burn.in>=1){
    BI = as.integer(burn.in)
  }else{
    stop("burn.in must be a positive integer.")
  }
  n = length(y)
  u = stats::runif(n,0.01,0.99)
  kn=2
  degree=2
  if(robust){
    out = Robust_vc(g, y, u, e,quant, iterations, kn, degree, prior, hyper,debugging)
  }else{
    out = NonRobust_vc(g, y, u, e,iterations, kn, degree, prior,debugging)
  }
  coefficient = list(GS.alpha=out$fit$GS.alpha[-c(1:BI),], GS.beta=out$fit$GS.beta[-c(1:BI),])
  fit = list(obj = out,coefficients = coefficient)
  return(fit)
}
