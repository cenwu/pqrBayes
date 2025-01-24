#' Estimation and estimation accuracy for a pqrBayes object
#'
#' calculate estimation and estimation accuracy for varying coefficients
#'
#' @param posterior posterior samples from the MCMC.
#' @param coefficient the matrix of true varying coefficients evaluated on the grid points.
#' @param u.grid the vector of grid points.
#' @param kn the number of interior knots for B-spline.
#' @param degree the degree of B-spline basis.
#' @param iterations the number of MCMC iterations.
#' @usage estimation.pqrBayes(posterior,coefficient,u.grid,kn,degree,iterations)
#' @return  an object of class `pqrBayes.est' is returned, which is a list with components:
#' \item{error}{integrated mean square errors and total integrated mean square error.}
#' \item{coeff.est}{estimated values of the varying coefficients.}
#' @seealso \code{\link{pqrBayes}}

#' @export
estimation.pqrBayes = function(posterior,coefficient,u.grid,kn=2,degree=2,iterations=10000){
  d=kn+degree+1
  u.star = seq(0, 1, length=kn+2)[-c(1,kn+2)]
  Knots.star = as.numeric(stats::quantile(u.grid, u.star))
  pi.star = splines::bs(u.grid, knots=Knots.star, intercept=TRUE, degree=degree)[,1:(d)]
  # 2.5% quantile of posterior samples
  
  c1_25.C=rep(0,d)
  for (i in 1:d) {
    c1_25.C[i]=stats::quantile(posterior$GS.alpha[(iterations/2+1):iterations,i],0.025)
  }
  
  c2_25.C=rep(0,dim(posterior$GS.beta)[2])
  for (i in 1:dim(posterior$GS.beta)[2]) {
    c2_25.C[i]=stats::quantile(posterior$GS.beta[(iterations/2+1):iterations,i],0.025)
  }
  
  coeffmatrix.C1_25=as.matrix(cbind(c1_25.C,matrix(c2_25.C,nrow = d)))
  
  
  # 97.5% quantile of posterior samples
  
  c1_975.C=rep(0,d)
  for (i in 1:d) {
    c1_975.C[i]=stats::quantile(posterior$GS.alpha[(iterations/2+1):iterations,i],0.975)
  }
  
  c2_975.C=rep(0,dim(posterior$GS.beta)[2])
  for (i in 1:dim(posterior$GS.beta)[2]) {
    c2_975.C[i]=stats::quantile(posterior$GS.beta[(iterations/2+1):iterations,i],0.975)
  }
  
  coeffmatrix.C2_975=as.matrix(cbind(c1_975.C,matrix(c2_975.C,nrow = d)))
  
  # 50% quantile of posterior samples
  
  c1.C=rep(0,d)
  for (i in 1:d) {
    c1.C[i]=stats::quantile(posterior$GS.alpha[(iterations/2+1):iterations,i],0.5)
  }
  
  c2.C=rep(0,dim(posterior$GS.beta)[2])
  for (i in 1:dim(posterior$GS.beta)[2]) {
    c2.C[i]=stats::quantile(posterior$GS.beta[(iterations/2+1):iterations,i],0.5)
  }
  
  coeffmatrix.C=as.matrix(cbind(c1.C,matrix(c2.C,nrow = d)))
  
  gamma.hat=pi.star%*% coeffmatrix.C[,1:ncol(coefficient)]
  gamma.hat_25 = pi.star%*% coeffmatrix.C1_25[,1:ncol(coefficient)]
  gamma.hat_975 = pi.star%*% coeffmatrix.C2_975[,1:ncol(coefficient)]
  IMSE = rep(0,ncol(coefficient))
  for(i in 1:ncol(coefficient)){
    IMSE[i] = mean((gamma.hat[,i]-coefficient[,i])^2)
  }
  TIMSE = sum(IMSE)
  error.est = list(IMSE = IMSE,TIMSE = TIMSE)
  coeff.est = list(gamma.hat = gamma.hat,gamma.hat_25 = gamma.hat_25,gamma.hat_975 = gamma.hat_975)
  pqrBayes.est = list(error=error.est, coeff.est=coeff.est)
  class(pqrBayes.est) = "pqrBayes.est"
  return(pqrBayes.est)
}
  