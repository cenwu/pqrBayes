#' Variable selection for a pqrBayes object
#'
#' Variable selection for a pqrBayes object
#'
#' @param obj pqrBayes object.
#' @param sparse logical flag. 
#' @param iterations the number of MCMC iterations.
#' @param kn the number of interior knots for B-spline.
#' @param degree the degree of B-spline basis.
#'
#' @details For class `Sparse', the median probability model (MPM) (Barbieri and Berger, 2004) is used to identify predictors that are significantly associated
#' with the response variable. For class `NonSparse', variable selection is based on 95\% credible interval.
#' Please check the references for more details about the variable selection.
#'
#' @references
#' Ren, J., Zhou, F., Li, X., Ma, S., Jiang, Y. and Wu, C. (2022). Robust Bayesian variable selection for gene-environment interactions. 
#' {\emph{Biometrics}, (in press)} \url{https://doi.org/10.1111/biom.13670}
#'
#' Barbieri, M.M. and Berger, J.O. (2004). Optimal predictive model selection. {\emph{Ann. Statist}, 32(3):870–897}
#'
#' @rdname select.VC
#' @return an object of class `VCselect' is returned, which includes the indices of the selected genetic factors.
#'
#' @seealso \code{\link{pqrBayes}}
#'
#'
#' @export
VCselect <- function(obj,sparse,iterations=10000,kn=2, degree=2){
  if(sparse){
    method="Sparse"
    idgene=obj$idgene
    id=which(idgene>iterations/4)
    VCselect=list(method=method,id=id)
  }else{
    method="Nonsparse"
    
    d=kn+degree+1
    
    # 2.5% quantile of posterior samples

    c2.C=rep(0,dim(obj$GS.beta)[2])
    for (i in 1:dim(obj$GS.beta)[2]) {
      c2.C[i]=stats::quantile(obj$GS.beta[round(iterations):iterations,i],0.025)
    }

    coeffmatrix.C1=matrix(c2.C,nrow = d)

    # 97.5% quantile of posterior samples

    c2.C=rep(0,dim(obj$GS.beta)[2])
    for (i in 1:dim(obj$GS.beta)[2]) {
      c2.C[i]=stats::quantile(obj$GS.beta[round(iterations):iterations,i],0.975)
    }

    coeffmatrix.C2=matrix(c2.C,nrow = d)
    
    idgene=rep(0,dim(coeffmatrix.C1)[2])
    for (i in 1:dim(coeffmatrix.C1)[2]) {
      for (j in 1:d) {
        if(coeffmatrix.C1[j,i]*coeffmatrix.C2[j,i]>0){idgene[i]=1}
      }
    }
    id=which(idgene==1)
    VCselect=list(method=method,id=id)
  }
  class(VCselect)="VCselect"
}
  
  
 