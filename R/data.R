#' simulated data for demonstrating the features of pqrBayes
#'
#' Simulated gene expression data for demonstrating the features of pqrBayes.
#'
#' @docType data
#' @keywords datasets
#' @name data
#' @format The data object consists of five components: g, y, u, e and coeff. coeff contains the true values of parameters used for generating the response variable \eqn{y}.
#'
#' 
#'
#' @details
#'
#' \strong{The model for generating Y}
#'
#' Use subscript \eqn{i} to denote the \eqn{i}th subject. Let \eqn{(\boldsymbol X_{i}, Y_{i}, U_{i}, \boldsymbol E_{i})}, (\eqn{i=1,\ldots,n}) be
#' independent and identically distributed random vectors. \eqn{Y_{i}} is a continuous response variable representing the
#' disease phenotype. \eqn{\boldsymbol X_{i}=(X_{i0},...,X_{ip})^\top} denotes a \eqn{1+p}--dimensional vector of genetic factors with the first element \eqn{X_{i0}=1}. 
#' The environmental factor \eqn{U_i \in \rm I\!R^1} is a univariate index variable. \eqn{\boldsymbol E_{i}=(E_{i1},...,E_{iq})^\top} is the \eqn{q}-dimensional vector 
#' of clinical covariates. The \eqn{\epsilon_i} follows some heavy-tailed distribution.
#' Considering the following varying coefficient model:
#'
#' \deqn{Y_{i}=\sum_{k=1}^{q} E_{ik} \beta_k +\sum_{j=0}^{p}\gamma_j(U_i)X_{ij} +\epsilon_{i},}
#' where \eqn{\beta_k}'s are the regression coefficients for the clinical covariates and \eqn{\gamma_j(\cdot)}'s are unknown smooth varying-coefficient functions. 
#' The regression coefficients of \eqn{\boldsymbol X} vary with the univariate index variable \eqn{\boldsymbol U=(U_1,...,U_n)^\top} 
#' 
#' The true model that we used to generate Y:
#' \deqn{Y_i=\gamma_0(U_i)+\gamma_1(U_i)X_{i1}+\gamma_2(U_i)X_{i2}+\gamma_3(U_i)X_{i3}+\epsilon_i,}
#' where \eqn{\epsilon_i\sim N(0,1)}, \eqn{\gamma_{0}=1.5\sin(0.2\pi*U_i}, \eqn{\gamma_{1}=2\exp(0.2U_i-1)-1.5 }, \eqn{\gamma_{2}=2-2U_i) } and \eqn{\gamma_3=-4+(U_i-2)^3/6}.
#'
#' @examples
#' data(data)
#' g=data$g
#' dim(g)
#' coeff=data$coeff
#' print(coeff)
#'
#'
#' @seealso \code{\link{pqrBayes}}
NULL
