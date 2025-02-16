#' simulated data under sparse linear and quantile varying coefficient models
#'
#' Simulated data under sparse linear and quantile varying coefficient models
#'
#' @docType data
#' @keywords datasets
#' @name data
#' @format 
#' The data_linear object consists of 4 components: g, y, e and coeff. coeff contains the true values of parameters used for generating the response variable \eqn{y}.
#' The data_varying object consists of five components: g, y, u, e and coeff. coeff contains the true values of parameters used for generating the response variable \eqn{y}.
#' 
#'
#' 
#'
#' @details
#'
#' \strong{Generating Y using a sparse linear (quantile) regression model}
#' 
#' The true data generating model under sparse linear regression:
#' \deqn{Y_i=\beta_0+\beta_{1}X_{i1}+\beta_{2}X_{i2}+\beta_{3}X_{i3}+\epsilon_i,}
#' where \eqn{\epsilon_i\sim N(0,1)}, \eqn{\beta_{0}=0}, \eqn{\beta_{1}=1 }, \eqn{\beta_{2}=1.5} and \eqn{\beta_3=2}.
#' 
#' \strong{Generating Y using a (quantile) varying coefficient model}
#'
#' Data generation under sparse (quantile) VC model:
#' \deqn{Y_i=\gamma_0(v_i)+\gamma_1(v_i)X_{i1}+\gamma_2(v_i)X_{i2}+\gamma_3(v_i)X_{i3}+\epsilon_i,}
#' where \eqn{\epsilon_i\sim N(0,1)}, \eqn{\gamma_{0}(v_i)=1.5\sin(0.2\pi*v_i}), \eqn{\gamma_{1}(v_i)=2\exp(0.2v_i-1)-1.5 }, \eqn{\gamma_{2}(v_i)=2-2v_i } and \eqn{\gamma_3(v_i)=-4+(v_i-2)^3/6}.
#'
#' @examples
#' data(data)
#' data = data$data_linear
#' g=data$g
#' dim(g)
#' y=data$y
#' coeff=data$coeff
#' print(coeff)
#' 
#' data = data$data_varying
#' g=data$g
#' dim(g)
#' coeff=data$coeff
#' print(coeff)
#'
#'
#' @seealso \code{\link{pqrBayes}}
NULL
