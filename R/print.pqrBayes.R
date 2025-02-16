
#' print a pqrBayes result
#'
#' Print a pqrBayes result
#'
#' @param x pqrBayes result.
#' @param digits significant digits in printout.
#' @param ... other print arguments.
#' @usage \method{print}{pqrBayes}(x, digits = max(3, getOption("digits") - 3), \dots)
#' @return No return value, called for side effects.
#' @seealso \code{\link{pqrBayes}}
#' @export
print.pqrBayes=function(x, digits = max(3, getOption("digits") - 3),...){
  print(x, digits)
}


#' print a pqrBayes.pred object
#'
#' Print a summary of a pqrBayes.pred object
#'
#' @param x pqrBayes.pred object.
#' @param digits significant digits in printout.
#' @param ... other print arguments
#' @usage \method{print}{pqrBayes.pred}(x, digits = max(3, getOption("digits") - 3), \dots)
#' @return No return value, called for side effects.
#' @seealso \code{\link{predict.pqrBayes}}
#' @export
print.pqrBayes.pred=function(x, digits = max(3, getOption("digits") - 3),...){
  cat("\ncheck loss:\n")
  print(x$error, digits)
  cat("\npredicted:\n ")
  print(x$y.pred,digits)
}


#' print a select.pqrBayes object
#'
#' Print a summary of a select.pqrBayes object
#'
#' @param x pqrBayes.select object.
#' @param digits significant digits in printout.
#' @param ... other print arguments
#' @usage \method{print}{pqrBayes.select}(x, digits = max(3, getOption("digits") - 3), \dots)
#' @return No return value, called for side effects.
#' @seealso \code{\link{pqrBayes.select}}
#' @export
print.pqrBayes.select=function(x, digits = max(3, getOption("digits") - 3),...){
  cat("\nMethod:\n")
  print(x$method)
  cat("\n")
  print(x$id)
}
