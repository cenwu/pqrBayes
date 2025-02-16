estimation_lin = function(obj,coefficient){
  coeff.est = c()
  iterations = obj$obj$iterations
  posterior = obj$coefficients
  for(j in 1:dim(posterior$GS.beta)[2])
  {
    t1=as.matrix(posterior$GS.beta[,j])
    t1 = t1[seq(iterations/2+1, iterations,1),]
    coeff.est[j] = stats::median(t1)
  }
  MSE = mean((coeff.est-coefficient)^2)
  error.est = list(MSE = MSE)
  est_lin = list(error=error.est, coeff.est=coeff.est)
  return(est_lin)
}