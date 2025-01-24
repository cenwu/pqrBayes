#' fit a regularized Bayesian quantile varying coefficient model
#' 
#' @keywords models
#' @param g the matrix of predictors (subject to selection) without intercept.
#' @param y the response variable. The current version only supports the continuous response.
#' @param u a vector of effect modifying variable of the quantile varying coefficient model.
#' @param e a matrix of clinical covariates not subject to selection.
#' @param quant the quantile level specified by users. The default value is 0.5.
#' @param iterations the number of MCMC iterations.
#' @param kn the number of interior knots for B-spline.
#' @param degree the degree of B-spline basis.
#' @param robust logical flag. If TRUE, robust methods will be used.
#' @param sparse logical flag. If TRUE, spike-and-slab priors will be used to shrink coefficients of irrelevant covariates to zero exactly.
#' @param hyper a named list of hyperparameters.
#' @param debugging logical flag. If TRUE, progress will be output to the console and extra information will be returned.
#'
#' @details The model described in "\code{\link{data}}" is:
#' \deqn{Y_{i}=\sum_{k=1}^{q} E_{ik} \beta_k +\sum_{j=0}^{p}\gamma_j(V_i)X_{ij} +\epsilon_{i},}
#' where \eqn{\beta_k}'s are the regression coefficients for the clinical covariates and \eqn{\gamma_j}'s are the varying coefficients for the intercept and predictors (e.g. genetic factors).
#'
#' 
#' When \{sparse=TRUE\} (default), spike--and--slab priors are adopted. Otherwise, Laplacian shrinkage will be used.
#
#
#' Users can modify the hyper-parameters by providing a named list of hyper-parameters via the argument `hyper'.
#' The list can have the following named components
#' \describe{
#'   \item{a0, b0: }{ shape parameters of the Beta priors (\eqn{\pi^{a_{0}-1}(1-\pi)^{b_{0}-1}}) on \eqn{\pi_{0}}.}
#'   \item{c1, c2: }{ the shape parameter and the rate parameter of the Gamma prior on \eqn{\nu}.}
#' }
#' 
#' Please check the references for more details about the prior distributions.
#' 
#' @return an object of class "pqrBayes" is returned, which is a list with components:
#' \item{posterior}{posterior samples from the MCMC}
#' \item{coefficients}{a list of posterior estimates of coefficients}
#' 
#' @examples
#' data(data)
#' g=data$g
#' y=data$y
#' u=data$u
#' e=data$e
#'
#' ## default method
#' fit1=pqrBayes(g,y,u,e,quant=0.5)
#' fit1
#'
#' \donttest{
#'
#' ## non-sparse
#' sparse=FALSE
#' fit2=pqrBayes(g,y,u,e,quant=0.5,sparse = sparse)
#' fit2
#' 
#' ## non-robust
#' robust = FALSE
#' fit3=pqrBayes(g,y,u,e,quant=0.5,robust = robust)
#' fit3
#' 
#' }
#' @export
#' @references
#' 
#' Zhou, F., Ren, J., Ma, S. and Wu, C. (2023). The Bayesian regularized quantile varying coefficient model.
#'  {\emph{Computational Statistics & Data Analysis}, 107808} \doi{10.1016/j.csda.2023.107808}
#'  
#' Ren, J., Zhou, F., Li, X., Ma, S., Jiang, Y. and Wu, C. (2023). Robust Bayesian variable selection for gene-environment interactions. 
#' {\emph{Biometrics}, 79(2), 684-694} \doi{10.1111/biom.13670}
#'
#' Ren, J., Zhou, F., Li, X., Chen, Q., Zhang, H., Ma, S., Jiang, Y. and Wu, C. (2020) Semi-parametric Bayesian variable selection for gene-environment interactions.
#' {\emph{Statistics in Medicine}, 39: 617– 638} \doi{10.1002/sim.8434}


pqrBayes <- function(g, y, u, e=NULL,quant=0.5, iterations=10000, kn=2, degree=2, robust=TRUE,sparse=TRUE, hyper=NULL,debugging=FALSE){
  
  if(robust){
    out = Robust(g, y, u, e,quant, iterations, kn, degree,sparse, hyper,debugging)
  }else{
    out = NonRobust(g, y, u, e,iterations, kn, degree,sparse,debugging)
  }
  coefficient = list(GS.alpha=out$GS.alpha, GS.beta=out$GS.beta)
  fit = list(posterior = out,coefficients = coefficient)
  return(fit)
}
