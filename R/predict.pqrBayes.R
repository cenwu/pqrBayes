#' Make predictions from a pqrBayes object
#'
#' Make predictions from a pqrBayes object
#'
#' @param object a pqrBayes object.
#' @param g.new a matrix of new predictors (e.g. genetic factors) at which predictions are to be made.
#' @param u.new a vector of new environmental factor at which predictions are to be made. When being applied to the linear model, u.new = NULL.
#' @param e.new a vector or matrix of new clinical covariates at which predictions are to be made.
#' @param y.new a vector of the response of new observations. If provided, the prediction error will be computed based on Y.new.
#' @param quant the quantile level.  The default is 0.5.
#' @param model the model to be fitted. The default is "VC" for a quantile varying coefficient model. Users can also specify "linear" for a linear model.
#' @param ... other predict arguments
#' 
#' @details g.new (u.new) must have the same number of columns as g (u) used for fitting the model. By default, the clinical covariates are NULL unless 
#' provided. The predictions are made based on the posterior estimates of coefficients in the pqrBayes object.
#'
#' If y.new is provided, the prediction error will be computed based on the check loss.
#'
#' @usage predict.pqrBayes(object, g.new, u.new, e.new, y.new, quant, model, ...)
#' @return  an object of class `pqrBayes.pred' is returned, which is a list with components:
#' \item{error}{prediction error. error is NULL if y.new=NULL.}
#' \item{y.pred}{predicted values of the new observations.}
#'
#' @rdname predict.pqrBayes
#' @seealso \code{\link{pqrBayes}}
#'
#' @export
predict.pqrBayes=function(object, g.new, u.new, e.new=NULL, y.new=NULL, quant=0.5,model,...){
  if(model=="VC"){
    pqrBayes.pred = predict_vc(object, g.new, u.new, e.new, y.new, quant,...)
  }else if(model=="linear"){
    pqrBayes.pred = predict_lin(object, g.new, e.new, y.new, quant,...)
  }
  else{
    stop("model should be either VC or linear")
  }
  class(pqrBayes.pred) = "pqrBayes.pred"
  return(pqrBayes.pred)
  #pred
}

