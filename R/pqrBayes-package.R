#' @useDynLib pqrBayes, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' @docType package
#' @keywords overview
#' @name pqrBayes-package
#' @title Regularized Bayesian Quantile Varying Coefficient Model
#' @aliases pqrBayes-package
#' @description In this package, we implement a Bayesian quantile varying coefficient model for non-linear gene-environment interaction analysis. The varying 
#' coefficient functions capture the possible non-linear gene-environment interactions and they are approximated using linear combinations of B-spline basis. 
#' Quantile regression is adopted as it's robust to long-tailed distributions in the response/phenotype and provides the capability of describing the relationship 
#' between the response variable and predictors at different quantiles of the response variable. The default method (the proposed method) conducts variable 
#' selection by accounting for sparsity. In particular, the spike--and--slab priors are adopted to shrink the coefficients of unimportant effects to exactly zero. 
#' In addition to the default method, users can also choose the method without spike--and--slab priors.
#' 
#' @details The user friendly, integrated interface \strong{pqrBayes()} allows users to flexibly choose the fitting methods by specifying the following parameter:
#' \tabular{rl}{
#' sparse: \tab whether to use the spike-and-slab priors to create sparsity. 
#' }
#'
#' The function pqrBayes() returns a pqrBayes object that contains the posterior estimates of each coefficients.
#'
#' @references
#' Zhou, F., Ren, J., Ma, S. and Wu, C. (2022). The Bayesian regularized quantile varying coefficient model. (submitted)
#' 
#' Ren, J., Zhou, F., Li, X., Ma, S., Jiang, Y. and Wu, C. (2022). Robust Bayesian variable selection for gene-environment interactions. 
#' {\emph{Biometrics}, (in press)} \url{https://doi.org/10.1111/biom.13670}
#'
#' Wu, C., and Ma, S. (2015). A selective review of robust variable selection with applications in bioinformatics.
#' {\emph{Briefings in Bioinformatics}, 16(5), 873–883} \url{https://doi.org/10.1093/bib/bbu046}
#'
#' Zhou, F., Ren, J., Lu, X., Ma, S. and Wu, C. (2020). Gene–Environment Interaction: a Variable Selection Perspective. Epistasis. Methods in Molecular Biology.
#' {\emph{Humana Press} (Accepted)} \url{https://arxiv.org/abs/2003.02930}
#'
#' Ren, J., Zhou, F., Li, X., Chen, Q., Zhang, H., Ma, S., Jiang, Y. and Wu, C. (2020) Semi-parametric Bayesian variable selection for gene-environment interactions.
#' {\emph{Statistics in Medicine}, 39: 617– 638} \url{https://doi.org/10.1002/sim.8434}
#'
#' Ren, J., Zhou, F., Li, X., Wu, C. and Jiang, Y. (2019) spinBayes: Semi-Parametric Gene-Environment Interaction via Bayesian Variable Selection.
#' R package version 0.1.0. \url{https://CRAN.R-project.org/package=spinBayes}
#'
#' Wu, C., Jiang, Y., Ren, J., Cui, Y. and Ma, S. (2018). Dissecting gene-environment interactions: A penalized robust approach accounting for hierarchical structures.
#' {\emph{Statistics in Medicine}, 37:437–456} \url{https://doi.org/10.1002/sim.7518}
#'
#' Wu, C., Shi, X., Cui, Y. and Ma, S. (2015). A penalized robust semiparametric approach for gene-environment interactions.
#' {\emph{Statistics in Medicine}, 34 (30): 4016–4030} \url{https://doi.org/10.1002/sim.6609}
#'
#' Wu, C., Cui, Y., and Ma, S. (2014). Integrative analysis of gene–environment interactions under a multi–response partially linear varying coefficient model.
#' {\emph{Statistics in Medicine}, 33(28), 4988–4998} \url{https://doi.org/10.1002/sim.6287}
#'
#' Wu, C., Zhong, P.S. and Cui, Y. (2018). Additive varying–coefficient model for nonlinear gene–environment interactions.
#' {\emph{Statistical Applications in Genetics and Molecular Biology}, 17(2)} \url{https://doi.org/10.1515/sagmb-2017-0008}
#'
#' Wu, C., Zhong, P.S. and Cui, Y. (2013). High dimensional variable selection for gene-environment interactions.
#' {\emph{Technical Report. Michigan State University.}}
#'
#' @seealso \code{\link{pqrBayes}}
NULL