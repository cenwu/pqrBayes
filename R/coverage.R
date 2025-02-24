#' 95\% empirical coverage probability for a pqrBayes object
#'
#'
#'@description
#'Calculate 95\% empirical coverage probabilities for regression coefficients under linear and VC models, respectively.
#'
#' @param object the pqrBayes object.
#' @param coefficient the vector of true regression coefficients under a linear model or the matrix of true varying coefficients evaluated on the grid points under a varying coefficient model.
#' @param u.grid the vector of grid points under a varying coefficient model. When assessing empirical coverage probabilities under a linear model, u.grid = NULL.
#' @param model the model to be fitted. Users can also choose "linear" for a linear model or "VC" for a varying coefficient model.
#' @usage coverage(object,coefficient,u.grid=NULL,model="linear")
#' @return c
#' @seealso \code{\link{pqrBayes}}
#' @examples
#' ## The quantile regression model
#' data(data)
#' data = data$data_linear
#' g=data$g
#' y=data$y
#' e=data$e
#' coeff = data$coeff
#' fit1=pqrBayes(g,y,u=NULL,e,quant=0.5,spline=NULL,model="linear")
#' coverage=coverage(fit1,coeff,model="linear")
#' @export
coverage = function(object,coefficient,u.grid=NULL,model="linear"){
  if(model=="VC"){
    c = coverage_vc(object,coefficient,u.grid)
  }
  else if(model=="linear"){
    c = coverage_lin(object,coefficient,u.grid=NULL)
  }
  else{
    stop("model should be either VC or linear")
  }
  return(c)
}
