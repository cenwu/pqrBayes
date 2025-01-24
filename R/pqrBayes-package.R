#' @useDynLib pqrBayes, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

"_PACKAGE"
#' @keywords overview
#' @name pqrBayes-package
#' @title Regularized Bayesian Quantile Varying Coefficient Model
#' @aliases pqrBayes-package
#' @description In this package, we implement a sparse Bayesian quantile varying coefficient model for non-linear gene-environment interactions. The quantile varying 
#' coefficient functions that can capture the non-linear gene-environment interactions are approximated using B-splines. 
#' Quantile regression is adopted as it's robust to long-tailed distributions in the response/phenotype and provides the capability of describing the relationship 
#' between the response variable and predictors at different quantile leves. The default method, Bayesian regularized quantile varying coefficient model with spike-and-slab priors,  
#' adopts the point-mass spike--and--slab priors to achieve exact sparsity by shrinking the coefficients of unimportant effects to exactly zero and
#' facilitate valid Bayesian inferences on quantile varying coefficients. 
#' In addition to the default method, users can also choose the method without robustness and spike--and--slab priors.
#' 
#' @details The user friendly, integrated interface \strong{pqrBayes()} allows users to flexibly choose the fitting methods by specifying the following parameter:
#' \tabular{rl}{
#' robust: \tab whether to fit the robust sparse quantile varying coefficient \cr\tab models or the non-robust counterpart. \cr\cr
#' sparse: \tab whether to use the spike-and-slab priors to impose sparsity. 
#' }
#'
#' The function pqrBayes() returns a pqrBayes object that contains the posterior estimates of each coefficients.
#'
#' @references
#' Zhou, F., Ren, J., Ma, S. and Wu, C. (2023). The Bayesian regularized quantile varying coefficient model.
#'  {\emph{Computational Statistics & Data Analysis}, 107808} \doi{10.1016/j.csda.2023.107808}
#'  
#' Ren, J., Zhou, F., Li, X., Ma, S., Jiang, Y. and Wu, C. (2023). Robust Bayesian variable selection for gene-environment interactions. 
#' {\emph{Biometrics}, 79(2), 684-694} \doi{10.1111/biom.13670}
#'
#' Wu, C., and Ma, S. (2015). A selective review of robust variable selection with applications in bioinformatics.
#' {\emph{Briefings in Bioinformatics}, 16(5), 873–883} \doi{10.1093/bib/bbu046}
#'
#' Zhou, F., Ren, J., Lu, X., Ma, S. and Wu, C. (2021). Gene–Environment Interaction: a Variable Selection Perspective. 
#' {\emph{Epistasis. Methods in Molecular Biology.} 2212:191–223} \url{https://link.springer.com/protocol/10.1007/978-1-0716-0947-7_13}
#'
#' Ren, J., Zhou, F., Li, X., Chen, Q., Zhang, H., Ma, S., Jiang, Y. and Wu, C. (2020) Semi-parametric Bayesian variable selection for gene-environment interactions.
#' {\emph{Statistics in Medicine}, 39: 617– 638} \doi{10.1002/sim.8434}
#'
#' Ren, J., Zhou, F., Li, X., Wu, C. and Jiang, Y. (2019) spinBayes: Semi-Parametric Gene-Environment Interaction via Bayesian Variable Selection.
#' R package version 0.1.0. \url{https://CRAN.R-project.org/package=spinBayes}
#'
#' Wu, C., Jiang, Y., Ren, J., Cui, Y. and Ma, S. (2018). Dissecting gene-environment interactions: A penalized robust approach accounting for hierarchical structures.
#' {\emph{Statistics in Medicine}, 37:437–456} \doi{10.1002/sim.7518}
#'
#' Wu, C., Shi, X., Cui, Y. and Ma, S. (2015). A penalized robust semiparametric approach for gene-environment interactions.
#' {\emph{Statistics in Medicine}, 34 (30): 4016–4030} \doi{10.1002/sim.6609}
#'
#' Wu, C., Cui, Y., and Ma, S. (2014). Integrative analysis of gene–environment interactions under a multi–response partially linear varying coefficient model.
#' {\emph{Statistics in Medicine}, 33(28), 4988–4998} \doi{10.1002/sim.6287}
#'
#' Wu, C., Zhong, P.S. and Cui, Y. (2018). Additive varying–coefficient model for nonlinear gene–environment interactions.
#' {\emph{Statistical Applications in Genetics and Molecular Biology}, 17(2)} \doi{10.1515/sagmb-2017-0008}
#'
#' Wu, C., Zhong, P.S. and Cui, Y. (2013). High dimensional variable selection for gene-environment interactions.
#' {\emph{Technical Report. Michigan State University.}}
#'
#' @seealso \code{\link{pqrBayes}}
NULL