
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pqrBayes

> Bayesian Penalized Quantile Regression

<!-- badges: start -->

[![CRAN](https://www.r-pkg.org/badges/version/pqrBayes)](https://cran.r-project.org/package=pqrBayes)
[![Codecov test
coverage](https://codecov.io/gh/cenwu/pqrBayes/branch/master/graph/badge.svg)](https://app.codecov.io/gh/cenwu/pqrBayes?branch=master)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/pqrBayes)](https://www.r-pkg.org:443/pkg/pqrBayes)
[![R-CMD-check](https://github.com/cenwu/pqrBayes/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/cenwu/pqrBayes/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The quantile varying coefficient model is robust to data heterogeneity, 
    outliers and heavy-tailed distributions in the response variable. In addition, 
    it can flexibly model dynamic patterns of regression coefficients through 
    nonparametric varying coefficient functions. In this package, we have implemented 
    the Gibbs samplers of the penalized Bayesian quantile varying coefficient model with 
    spike-and-slab priors [Zhou et al.(2023)]<doi:10.1016/j.csda.2023.107808> for efficient 
    Bayesian shrinkage estimation, variable selection and statistical inference. In particular,
    valid Bayesian inferences on sparse quantile varying coefficient functions can be validated 
    on finite samples. The Markov Chain Monte Carlo (MCMC) algorithms of the proposed
    and alternative models can be efficiently performed by using the package.   

## How to install

<!-- -->

    install.packages("devtools")
    devtools::install_github("cenwu/pqrBayes")

- Released versions of pqrBayes are available on CRAN
  <!-- [(link)](https://cran.r-project.org/package=pqrBayes) --> , and
  can be installed within R via

<!-- -->

    install.packages("pqrBayes")


## News

### pqrBayes 1.0.3 \[2024-12-21\]

- Fixed the issue of no output in the VCselect() function and added examples.
- Updated the list of pqrBayes output objects.
- Added non-robust sparse Bayesian varying coefficient models and corresponding examples.
- Updated the documentation.

## Methods

This package provides implementation for methods proposed in

  - Zhou, F., Ren, J., Ma, S. and Wu, C. (2023). The Bayesian Regularized Quantile Varying Coefficient Model.  {\emph{Computational Statistics & Data Analysis}, 107808} \doi{10.1016/j.csda.2023.107808}
