#' Variable selection for a pqrBayes object
#'
#' Variable selection for a pqrBayes object
#'
#' @param object a pqrBayes object.
#' @param prior the prior used in the pqrBayes function. Users can choose "SS" for the spike-and-slab prior, "HS" for the horseshoe prior, "HS+" for the horseshoe plus prior, "RHS" for the regularized horseshoe prior and "Laplace" for the Laplace prior. The default value is "SS".
#' @param model the model to be fitted. Users can also choose "linear" for a sparse linear model, "VC" for a varying coefficient model or "group" for group LASSO.
#'
#' @details For class `Sparse', the median probability model (MPM) (Barbieri and Berger, 2004) is used to identify predictors that are significantly associated
#' with the response variable. For class `NonSparse', variable selection is based on 95\% credible interval.
#' Please check the references for more details about the variable selection.
#'
#' @references
#' Ren, J., Zhou, F., Li, X., Ma, S., Jiang, Y. and Wu, C. (2023). Robust Bayesian variable selection for gene-environment interactions. 
#' {\emph{Biometrics}, 79(2), 684-694} \doi{10.1111/biom.13670}
#'
#' Barbieri, M.M. and Berger, J.O. (2004). Optimal predictive model selection. {\emph{Ann. Statist}, 32(3):870–897}
#'
#' @rdname select.pqrBayes
#' @usage pqrBayes.select(object,prior="SS",model="linear")
#' @return an object of class `select' is returned, which includes the indices of the selected predictors (e.g. genetic factors).
#'
#' @seealso \code{\link{pqrBayes}}
#'
#' @examples
#' ## The sparse quantile regression model
#' data(data)
#' data = data$data_linear
#' g=data$g
#' y=data$y
#' e=data$e
#' 
#' fit1=pqrBayes(g,y,e,d = NULL,quant=0.5,model="linear")
#' select=pqrBayes.select(obj = fit1,prior = "SS",model="linear")
#' 
#' ## The quantile varying coefficient model
#' data(data)
#' data = data$data_varying
#' g=data$g
#' y=data$y
#' e=data$e
#' fit1=pqrBayes(g,y,e,d = NULL,quant=0.5,model="VC")
#' select=pqrBayes.select(obj = fit1,prior = "SS",model="VC")
#' select
#'
#' \donttest{
#' ## Non-sparse example with VC model
#' fit2 <- pqrBayes(
#'   g = g, y = y, e = e, d = NULL,
#'   quant = 0.5,
#'   prior= "Laplace",
#'   model = "VC"
#' )
#'
#' select <- pqrBayes.select(obj = fit2, prior = "SS", model = "VC")
#' select
#' }

#'
#' @export
pqrBayes.select <- function(object,prior="SS",model="linear"){
  if(model=="VC"){
    select = VCselect(obj = object,prior)
  }
  else if(model=="linear"){
    select = linselect(obj = object,prior)
  }
  else if(model=="binary"){
    select = linselect(obj = object,prior)
  }
  else if(model=="group"){
    select = linselect(obj = object,prior)
  }
  else{
    stop("model should be either VC, linear, binary or group")
  }
  class(select)="select"
  return(select)
}
  
  
 