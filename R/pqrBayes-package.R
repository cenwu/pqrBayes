#' @useDynLib pqrBayes, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

"_PACKAGE"
#' @keywords overview
#' @name pqrBayes-package
#' @title Bayesian penalized quantile regression for the sparse linear model, binary LASSO, group LASSO and varying coefficient models based on spike-and-slab priors and/or the horseshoe family of priors
#' @aliases pqrBayes-package
#' @description In this package, we implement Bayesian penalized quantile regression for sparse linear, binary LASSO, group LASSO and quantile varying coefficient (VC) models. 
#' The point-mass spike-and-slab priors and horseshoe family of (horseshoe, horseshoe+ and regularized horseshoe) priors have been incorporated in the Bayesian hierarchical
#' models to facilitate Bayesian shrinkage estimation, variable selection and statistical inference. For the spike-and-slab priors, the four default methods are Bayesian
#' regularized quantile regression with spike-and-slab priors under the sparse linear (i.e. LASSO), binary LASSO, group LASSO and VC model, correspondingly. In addition to the default methods, 
#' users can also choose models without robustness and/or spike--and--slab priors. Furthermore, under sparse linear models, we have implemented the horseshoe family of (horseshoe, horseshoe+
#' and regularized horseshoe) priors and the three non-robust alternatives. Currently, the horseshoe family of priors is only implemented under the sparse linear model.
#' 
#' @details The user friendly, integrated interface \strong{pqrBayes()} allows users to flexibly choose fitting models by specifying the following parameters:
#' \tabular{rl}{
#' robust: \tab whether to fit a robust sparse quantile regression model (the sparse linear model, binary LASSO, \cr\tab group LASSO or Varying Coefficient models) or their non-robust counterparts. \cr\cr
#' prior: \tab specify which prior to use (the spike-and-slab prior, Laplace prior \cr\tab and the horseshoe family of priors). \cr\cr
#' model: \tab whether to fit a sparse linear model, binary LASSO, group LASSO \cr\tab or a varying coefficient model.
#' }
#'
#' The function pqrBayes() returns a pqrBayes object that stores the posterior estimates of regression coefficients.
#'
#' @references
#' Fan, K., Subedi, S., Yang, G., Lu, X., Ren, J. and Wu, C. (2024). Is Seeing Believing? A Practitioner's Perspective on High-dimensional Statistical Inference in Cancer Genomics Studies. {\emph{Entropy}, 26(9).794} \doi{10.3390/e26090794}  
#'
#' Zhou, F., Ren, J., Ma, S. and Wu, C. (2023). The Bayesian regularized quantile varying coefficient model.
#'  {\emph{Computational Statistics & Data Analysis}, 107808} \doi{10.1016/j.csda.2023.107808}
#'  
#' Ren, J., Zhou, F., Li, X., Ma, S., Jiang, Y. and Wu, C. (2023). Robust Bayesian variable selection for gene-environment interactions. 
#' {\emph{Biometrics}, 79(2), 684-694} \doi{10.1111/biom.13670}
#'
#' Fan, K. and Wu, C. (2025). A New Robust Binary Bayesian LASSO. 
#' {\emph{Stat}, 14 (3), e70078} \doi{10.1002/sta4.70078}
#' 
#' Fan, K., Srijana, S., Dissanayake, V. and Wu, C. (2025). Robust Bayesian high-dimensional variable selection and inference with the horseshoe family of priors 
#' {\emph{arXiv:2507.10975}} \url{https://arxiv.org/abs/2507.10975} 
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